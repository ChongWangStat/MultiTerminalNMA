import numpy as np, pandas as pd
from itertools import combinations
from scipy.linalg import block_diag
from scipy.optimize import linprog
from collections import defaultdict

TOL=1e-10
CC=0.5

covid=pd.DataFrame({
'study':['RECOVERY']*4+['SOLIDARITY']*4+['ACTT-1']*2+['CoDEX']*2+['REMAP-CAP']*2,
'id':[1]*4+[2]*4+[3]*2+[4]*2+[5]*2,
't':['A','B','C','D','A','C','D','E','A','E','A','B','A','B'],
'r':[1110,482,353,396,440,104,148,301,75,59,91,85,33,78],
'n':[4321,2104,1596,1561,4088,947,1399,2743,521,541,148,151,58,136]
})


def pinv(a): return np.linalg.pinv(a, rcond=1e-12)

def build(dat, tau2=0.0):
    trts=sorted(dat.t.unique()); ref=trts[0]
    pairs=list(combinations(trts,2)); pair_names=[f'{u}:{v}' for u,v in pairs]
    blocks=[]
    for sid,ds in dat.groupby('id',sort=True):
        ds=ds.sort_values('t')
        arms=list(ds.t); r=ds.r.to_numpy(float); n=ds.n.to_numpy(float)
        r_adj=r+CC; n_adj=n+2*CC
        z=np.log(r_adj/(n_adj-r_adj))
        s2=1/r_adj + 1/(n_adj-r_adj) + tau2/2.0
        local=list(combinations(range(len(arms)),2))
        M=np.zeros((len(local),len(arms))); y=np.zeros(len(local)); names=[]; row_names=[]
        for i,(a,b) in enumerate(local):
            M[i,a]=-1; M[i,b]=1; y[i]=z[b]-z[a]
            nm=f'{arms[a]}:{arms[b]}'; names.append(nm); row_names.append(f'{sid}.{ds.study.iloc[0]}::{nm}')
        V=M@np.diag(s2)@M.T
        X=np.zeros((len(local),len(pair_names)))
        for i,nm in enumerate(names): X[i,pair_names.index(nm)]=1
        blocks.append(dict(study=ds.study.iloc[0],id=int(sid),arms=arms,y=y,V=V,X=X,M=M,names=names,row_names=row_names))
    ty=np.concatenate([b['y'] for b in blocks]); tX=np.vstack([b['X'] for b in blocks]); tV=block_diag(*[b['V'] for b in blocks])
    B=np.zeros((len(pairs),len(trts)-1)); nonref=[x for x in trts if x!=ref]
    for i,(u,v) in enumerate(pairs):
        if u==ref: B[i,nonref.index(v)]=1
        else: B[i,nonref.index(u)]=-1; B[i,nonref.index(v)]=1
    Vp=pinv(tV); XB=tX@B; I=XB.T@Vp@XB
    P=B@pinv(I)@XB.T@Vp
    theta=P@ty
    return dict(trts=trts,pairs=pairs,pair_names=pair_names,blocks=blocks,P=P,theta=theta)

def flow_edges_for_study(pk,names,study,sid):
    edges=[]
    for w,nm in zip(pk,names):
        if abs(w)<=TOL: continue
        u,v=nm.split(':')
        if w>0: fr,to,wt=u,v,w
        else: fr,to,wt=v,u,-w
        edges.append(dict(id=sid,study=study,fr=fr,to=to,w=float(wt)))
    return edges

def choose_source(edges):
    nodes=sorted(set([e['fr'] for e in edges]+[e['to'] for e in edges]))
    cands=[]
    for x in nodes:
        infl=sum(e['w'] for e in edges if e['to']==x and e['w']>TOL)
        out=sum(e['w'] for e in edges if e['fr']==x and e['w']>TOL)
        if infl<=TOL and out>TOL: cands.append((-out,x))
    return sorted(cands)[0][1] if cands else None

def dfs_heaviest(edges,source):
    def rec(cur,visited):
        cand=[i for i,e in enumerate(edges) if e['fr']==cur and e['w']>TOL]
        cand=sorted(cand,key=lambda i:(-edges[i]['w'],edges[i]['to'],edges[i]['study'],i))
        if not cand: return []
        for i in cand:
            nxt=edges[i]['to']
            if nxt in visited: continue
            sub=rec(nxt,visited|{nxt})
            # exact R behavior returns c(ix, integer(0)) at sink, where no outgoing => integer(0)
            if sub is not None: return [i]+sub
        return None
    z=rec(source,{source}); return z if z else None

def collapse_one_study(pk,b):
    work=flow_edges_for_study(pk,b['names'],b['study'],b['id'])
    collapsed=[]
    while work:
        src=choose_source(work)
        if src is None: break
        while True:
            idx=dfs_heaviest(work,src)
            if not idx: break
            seg=[work[i] for i in idx]; wt=min(e['w'] for e in seg)
            collapsed.append(dict(study=b['study'],id=b['id'],start=seg[0]['fr'],end=seg[-1]['to'],w=wt))
            for i in idx: work[i]['w']-=wt
            work=[e for e in work if e['w']>TOL]
            if not work or not any(e['fr']==src for e in work): break
    return collapsed,work

def legacy_endpoint_edges(fit,target_idx):
    tp=fit['pair_names'][target_idx]; src,sink=tp.split(':')
    all_edges=[]; direct=[]; leftovers=[]
    offset=0
    for b in fit['blocks']:
        m=len(b['names']); pk=fit['P'][target_idx,offset:offset+m]; offset+=m
        col,left=collapse_one_study(pk,b)
        leftovers.extend(left)
        # aggregate collapsed endpoint edges within study, matching aggregate_edge_df() in the R source
        dd=defaultdict(float)
        for c in col: dd[(c['study'],c['id'],c['start'],c['end'])]+=c['w']
        for (st,sid,fr,to),wt in dd.items():
            c=dict(study=st,id=sid,start=fr,end=to,w=wt)
            if fr==src and to==sink: direct.append(c)
            else: all_edges.append(dict(study=st,id=sid,fr=fr,to=to,w=wt))
    return tp,src,sink,direct,all_edges,leftovers

def endpoint_edges(fit,target_idx):
    # v10 uses the supplied CSP balance-based endpoint rule in both applications.
    # The v09 edge-path-collapse rule remains above for explicit compatibility checks.
    from psoriasis_csp import endpoint_edges as canonical_endpoint_edges
    tp=fit['pair_names'][target_idx];src,sink=tp.split(':')
    full=canonical_endpoint_edges(fit,target_idx,tol=TOL)
    direct=[dict(study=e['study'],id=e['id'],start=e['fr'],end=e['to'],w=e['w'])
            for e in full if e['type']=='Dir']
    indirect=[{key:e[key] for key in ['study','id','fr','to','w']}
              for e in full if e['type']=='Ind']
    return tp,src,sink,direct,indirect,[]

def enumerate_paths(edges,src,sink,study_label=True):
    adj=defaultdict(list)
    for i,e in enumerate(edges): adj[e['fr']].append(i)
    paths=[]
    def rec(v,verts,idxs):
        if v==sink:
            paths.append(tuple(idxs)); return
        for i in adj.get(v,[]):
            e=edges[i]
            if e['to'] in verts: continue
            rec(e['to'],verts|{e['to']},idxs+[i])
    rec(src,{src},[])
    return paths

def coarsen_edges(edges):
    d=defaultdict(float)
    for e in edges:d[(e['fr'],e['to'])]+=e['w']
    return [dict(study='AGG',id=0,fr=k[0],to=k[1],w=w) for k,w in d.items() if w>TOL]

def audit_edges(edges,src,sink):
    paths=enumerate_paths(edges,src,sink)
    if not paths:
        return dict(n_edges=len(edges),n_paths=0,rank=0,nullity=0,unique=True,max_width=0.0,widths=[])
    A=np.zeros((len(edges),len(paths)))
    f=np.array([e['w'] for e in edges])
    for j,p in enumerate(paths):
        for i in p:A[i,j]=1
    rank=np.linalg.matrix_rank(A,tol=1e-9); nullity=len(paths)-rank
    widths=[]
    for j in range(len(paths)):
        c=np.zeros(len(paths)); c[j]=1
        lo=linprog(c,A_eq=A,b_eq=f,bounds=(0,None),method='highs')
        hi=linprog(-c,A_eq=A,b_eq=f,bounds=(0,None),method='highs')
        if not(lo.success and hi.success):
            widths.append(np.nan)
        else: widths.append((-hi.fun)-lo.fun)
    return dict(n_edges=len(edges),n_paths=len(paths),rank=rank,nullity=nullity,unique=(nullity==0),max_width=float(np.nanmax(widths) if widths else 0),widths=widths,A=A,paths=paths)

def one_run(tau2=0.0):
    fit=build(covid,tau2=tau2); rows=[]; edge_rows=[]
    for ti,tp in enumerate(fit['pair_names']):
        tp,src,sink,direct,edges,left=endpoint_edges(fit,ti)
        a=audit_edges(edges,src,sink); ca=audit_edges(coarsen_edges(edges),src,sink)
        rows.append(dict(target=tp,theta=fit['theta'][ti],direct_weight=sum(x['w'] for x in direct),indirect_weight=sum(e['w'] for e in edges if e['fr']==src)-0.0,
                         study_edges=a['n_edges'],study_paths=a['n_paths'],study_rank=a['rank'],study_nullity=a['nullity'],study_unique=a['unique'],study_max_width=a['max_width'],
                         treatment_edges=ca['n_edges'],treatment_paths=ca['n_paths'],treatment_rank=ca['rank'],treatment_nullity=ca['nullity'],treatment_unique=ca['unique'],treatment_max_width=ca['max_width'],leftover_weight=sum(x['w'] for x in left)))
        for e in edges: edge_rows.append(dict(target=tp,**e,tau2=tau2))
    return pd.DataFrame(rows),pd.DataFrame(edge_rows)

if __name__ == '__main__':
    from pathlib import Path
    import argparse
    parser = argparse.ArgumentParser(description='Recompute COVID fixed and conditional covariance-sensitivity flows.')
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent / 'reproduced')
    parser.add_argument('--tau2', type=float, default=0.002103984)
    args = parser.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    for model,tau2 in [('common',0.0),('random',args.tau2)]:
        main,edges = one_run(tau2)
        assert len(main)==10 and main.leftover_weight.max()<1e-9
        main.to_csv(args.output_dir / f'covid_route_audit_{model}.csv',index=False)
        edges.to_csv(args.output_dir / f'covid_route_edges_{model}.csv',index=False)
        print(model,'study unique',int(main.study_unique.sum()),'treatment unique',int(main.treatment_unique.sum()),flush=True)
