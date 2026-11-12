"""Numerical translation of the supplied CSP R construction, with independent GLS checks.

The uploaded R sources are retained unchanged. This module does not claim to execute R.
Both fits use the supplied unconditional 0.5 arm correction and intact-study covariance.
"""
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.linalg import block_diag
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
TOL = 1e-12


def read_data(path=ROOT/'pso-psa-data.csv'):
    raw = pd.read_csv(path)
    dat = raw.loc[raw.indication.eq('psoriasis')].copy()
    assert dat[['study','treatment','r','n']].notna().all().all()
    assert not dat.duplicated(['study','treatment']).any()
    assert (dat.n > 0).all() and (dat.r >= 0).all() and (dat.r <= dat.n).all()
    assert ((dat.n % 1 == 0) & (dat.r % 1 == 0)).all()
    assert set(dat.endpoint) == {'pasi75'} and set(dat.time) == {12}
    ids = {s:i+1 for i,s in enumerate(sorted(dat.study.unique()))}
    return pd.DataFrame({'study':'Study_'+dat.study.astype(str),
                         'id':dat.study.map(ids), 't':dat.treatment,
                         'r':dat.r.astype(float),'n':dat.n.astype(float)})


def ginv(a):
    # MASS::ginv default threshold, rather than numpy's different default.
    return np.linalg.pinv(a, rcond=np.sqrt(np.finfo(float).eps))


def prepare(dat):
    trts = sorted(dat.t.unique()); pairs = list(combinations(trts,2))
    names = [f'{a}:{b}' for a,b in pairs]
    B = np.zeros((len(pairs),len(trts)-1))
    for i,(a,b) in enumerate(pairs):
        if a != trts[0]: B[i,trts.index(a)-1] = -1
        if b != trts[0]: B[i,trts.index(b)-1] = 1
    blocks=[]
    for sid,dd in dat.groupby('id',sort=True):
        dd=dd.sort_values('t'); arms=list(dd.t)
        r=dd.r.to_numpy()+.5; n=dd.n.to_numpy()+1.
        z=np.log(r/(n-r)); s2=1/r+1/(n-r)
        local=list(combinations(range(len(dd)),2))
        M=np.zeros((len(local),len(dd))); X=np.zeros((len(local),len(pairs)))
        local_names=[]; row_names=[]
        for i,(a,b) in enumerate(local):
            M[i,a]=-1; M[i,b]=1
            name=f'{arms[a]}:{arms[b]}';local_names.append(name)
            row_names.append(f'{sid}.{dd.study.iloc[0]}::{name}')
            X[i,names.index(name)]=1
        blocks.append(dict(id=int(sid),study=dd.study.iloc[0],arms=arms,z=z,s2=s2,
                           M=M,X=X,y=M@z,V_within=M@np.diag(s2)@M.T,
                           names=local_names,row_names=row_names))
    return dict(trts=trts,pairs=pairs,pair_names=names,B=B,blocks=blocks,
                y=np.concatenate([b['y'] for b in blocks]),
                X=np.vstack([b['X'] for b in blocks]))


def fit(prepared,tau2=0.):
    q=prepared; V=block_diag(*[b['V_within']+.5*tau2*(b['M']@b['M'].T) for b in q['blocks']])
    Vp=ginv(V); XB=q['X']@q['B']; I=XB.T@Vp@XB
    P=q['B']@ginv(I)@XB.T@Vp
    beta=ginv(I)@XB.T@Vp@q['y']; resid=q['y']-XB@beta
    def lpdet(a):
        ev=np.linalg.eigvalsh(a);return np.log(ev[ev>TOL]).sum()
    objective=.5*(lpdet(V)+lpdet(I)+resid@Vp@resid)
    return dict(**q,P=P,theta=P@q['y'],V=V,tau2=float(tau2),objective=float(objective))


def estimate_tau2(q):
    upper=max(1.,float(np.var(q['y'],ddof=1)))
    for _ in range(21):
        opt=minimize_scalar(lambda t:fit(q,t)['objective'],bounds=(0.,upper),
                            method='bounded',options={'xatol':1e-10})
        if opt.x<.95*upper: break
        upper*=2
    if not opt.success: raise RuntimeError(opt.message)
    # R optimize searches the interior. A boundary comparison is additionally recorded.
    zero=fit(q,0.)['objective']
    best=0. if zero<opt.fun else float(opt.x)
    return best,dict(tau2=best,interior_optimum=float(opt.x),objective=float(opt.fun),
                     objective_at_zero=zero,upper_search_bound=upper,boundary_selected=best==0.)


def independent_gls(q,tau2):
    """Nonredundant study-reference contrasts; ordinary inverse, no pseudoinverse."""
    ys=[]; xs=[]; vs=[]
    for b in q['blocks']:
        k=len(b['arms']); C=np.zeros((k-1,k));C[:,0]=-1;C[:,1:]=np.eye(k-1)
        ys.append(C@b['z']);vs.append(C@np.diag(b['s2']+tau2/2)@C.T)
        Xarm=np.zeros((k,len(q['trts'])-1))
        for i,t in enumerate(b['arms']):
            if t!=q['trts'][0]:Xarm[i,q['trts'].index(t)-1]=1
        xs.append(C@Xarm)
    y=np.concatenate(ys);X=np.vstack(xs);V=block_diag(*vs)
    W=np.linalg.inv(V);I=X.T@W@X;beta=np.linalg.solve(I,X.T@W@y)
    resid=y-X@beta
    return q['B']@beta,float(.5*(np.linalg.slogdet(V)[1]+np.linalg.slogdet(I)[1]+resid@W@resid))


def endpoint_edges(result,target_index,tol=TOL):
    """Port of .canonical_endpoint_edges: direct first, largest residual pair, lexical tie."""
    target=result['pair_names'][target_index];a,b=target.split(':');out=[];offset=0
    for block in result['blocks']:
        k=len(block['names']);p=result['P'][target_index,offset:offset+k].copy();offset+=k
        p[np.abs(p)<=tol]=0
        balance=-block['M'].T@p
        supply=np.maximum(balance,0);demand=np.maximum(-balance,0);arms=block['arms']
        def add(i,j,w,kind):
            if w>tol:
                out.append(dict(target=target,id=block['id'],study=block['study'],
                                fr=arms[i],to=arms[j],w=float(w),type=kind))
                supply[i]-=w;demand[j]-=w
                if supply[i]<=tol:supply[i]=0.
                if demand[j]<=tol:demand[j]=0.
        if a in arms and b in arms:
            i,j=arms.index(a),arms.index(b);add(i,j,min(supply[i],demand[j]),'Dir')
        for _ in range(2*len(arms)+10):
            ii=np.flatnonzero(supply>tol);jj=np.flatnonzero(demand>tol)
            if len(ii)==0 and len(jj)==0:break
            if len(ii)==0 or len(jj)==0:raise RuntimeError('Unmatched study balance')
            candidates=[(min(supply[i],demand[j]),i,j) for i in ii for j in jj]
            maximum=max(v[0] for v in candidates)
            tied=[v for v in candidates if abs(v[0]-maximum)<=tol*max(1.,abs(maximum))]
            w,i,j=min(tied,key=lambda v:(arms[v[1]],arms[v[2]]));add(i,j,w,'Ind')
        else:raise RuntimeError('Endpoint allocation did not terminate')
        reconstructed=np.zeros(len(arms))
        for e in out:
            if e['id']==block['id']:
                reconstructed[arms.index(e['fr'])]+=e['w'];reconstructed[arms.index(e['to'])]-=e['w']
        assert np.max(np.abs(reconstructed-balance))<1e-10
    return out
