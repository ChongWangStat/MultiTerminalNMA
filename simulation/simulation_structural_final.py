from pathlib import Path
import argparse
_parser = argparse.ArgumentParser(description="Frozen v05 simulation; v08 portable output path.")
_parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "reproduced")
OUTPUT = _parser.parse_args().output_dir.resolve()
OUTPUT.mkdir(parents=True, exist_ok=True)

import numpy as np, pandas as pd, math, itertools, time
from collections import defaultdict, deque

SEED=20261004
rng=np.random.default_rng(SEED)

def connected_treatments(studies,T):
    adj=[set() for _ in range(T)]
    for arms in studies:
        for i in arms:
            adj[i].update(j for j in arms if j!=i)
    seen={0}; q=[0]
    while q:
        u=q.pop()
        for v in adj[u]:
            if v not in seen: seen.add(v); q.append(v)
    return len(seen)==T

def gen_design(T, rep_level, p_multi):
    # base connected backbone + random studies. Repetition intensity increases total studies
    mult={'low':1.2,'moderate':2.0,'high':3.2}[rep_level]
    K=max(T-1, int(math.ceil(mult*T)))
    studies=[]
    # random spanning tree as two-arm backbone
    perm=rng.permutation(T)
    for j in range(1,T):
        parent=perm[rng.integers(0,j)]
        studies.append(tuple(sorted((perm[j],parent))))
    while len(studies)<K:
        m=2
        if rng.random()<p_multi:
            m=3 if T<4 or rng.random()<0.75 else 4
            m=min(m,T)
        if rep_level=='high' and rng.random()<0.45 and studies:
            # deliberately replicate a treatment set/comparison often
            base=studies[rng.integers(len(studies))]
            if len(base)==m or m==2:
                if m==2:
                    if len(base)>=2:
                        arms=tuple(sorted(rng.choice(base,2,replace=False)))
                    else: arms=base
                else: arms=base
            else:
                arms=tuple(sorted(rng.choice(T,m,replace=False)))
        else:
            arms=tuple(sorted(rng.choice(T,m,replace=False)))
        studies.append(arms)
    return studies

def study_conductances(studies,dispersion):
    gs=[]
    for arms in studies:
        if dispersion=='balanced':
            x=np.ones(len(arms))
        else:
            x=np.exp(rng.normal(0,0.8,size=len(arms)))
        gs.append(dict(zip(arms,x)))
    return gs

def treatment_laplacian(studies,gs,T):
    L=np.zeros((T,T),float)
    for arms,g in zip(studies,gs):
        arr=np.array(arms,int); gv=np.array([g[a] for a in arms],float); d=gv.sum()
        K=np.diag(gv)-np.outer(gv,gv)/d
        L[np.ix_(arr,arr)] += K
    return L

def solve_potential(L,a,b):
    T=L.shape[0]; q=np.zeros(T); q[a]=1; q[b]=-1
    # pseudoinverse gauge-invariant
    phi=np.linalg.pinv(L,rcond=1e-12)@q
    return phi

def canonical_edges(studies,gs,phi,target,tol=1e-11):
    a,b=target; edges=[]
    for k,(arms,g) in enumerate(zip(studies,gs)):
        gv=np.array([g[x] for x in arms],float); pv=np.array([phi[x] for x in arms])
        psi=float((gv*pv).sum()/gv.sum())
        bal={x: g[x]*(phi[x]-psi) for x in arms}
        s={x:max(v,0.0) for x,v in bal.items()}; d={x:max(-v,0.0) for x,v in bal.items()}
        def add(u,v,w,typ):
            if w>tol: edges.append((u,v,k,float(w),typ))
        if a in s and b in d:
            w=min(s[a],d[b])
            if w>tol:
                add(a,b,w,'Dir'); s[a]-=w; d[b]-=w
        # largest feasible transfer, deterministic lexicographic tie-break
        for _ in range(2*len(arms)+10):
            src=[u for u,v in s.items() if v>tol]; snk=[v for v,x in d.items() if x>tol]
            if not src and not snk: break
            if not src or not snk: raise RuntimeError('unmatched')
            cands=[]
            for u in src:
                for v in snk:
                    cands.append((min(s[u],d[v]),u,v))
            w,u,v=max(cands,key=lambda z:(z[0],-z[1],-z[2]))
            add(u,v,w,'Ind'); s[u]-=w; d[v]-=w
        else: raise RuntimeError('allocation fail')
    # normalize tiny numerical residual of total source outflow should already 1
    return edges

def positive_indirect_edges(edges,target,tol=1e-10):
    # route analysis is for complete indirect flow: exclude direct target edge(s)
    a,b=target
    return [(u,v,k,w) for u,v,k,w,typ in edges if typ=='Ind' and w>tol]

def graph_counts(edges,a,b,T):
    # edges study-labelled, DAG by phi; use active vertices; count complete paths via DP
    out=defaultdict(list); indeg=defaultdict(int); verts={a,b}
    for idx,(u,v,k,w) in enumerate(edges):
        out[u].append((v,idx)); indeg[v]+=1; indeg.setdefault(u,indeg.get(u,0)); verts|={u,v}
    # reachable from a and to b filter
    reach={a}; stack=[a]
    while stack:
        u=stack.pop()
        for v,_ in out[u]:
            if v not in reach: reach.add(v); stack.append(v)
    rev=defaultdict(list)
    for u in out:
        for v,idx in out[u]: rev[v].append(u)
    tob={b}; stack=[b]
    while stack:
        v=stack.pop()
        for u in rev[v]:
            if u not in tob: tob.add(u); stack.append(u)
    active=reach&tob
    eactive=[e for e in edges if e[0] in active and e[1] in active]
    if a not in active or b not in active: return None
    # topological order on active graph (parallel edges don't alter indegree by unique arcs; okay count indegree edges too)
    indeg2={v:0 for v in active}; out2=defaultdict(list)
    for e in eactive:
        u,v,k,w=e; out2[u].append(e); indeg2[v]+=1
    q=deque([v for v in active if indeg2[v]==0]); topo=[]
    while q:
        u=q.popleft(); topo.append(u)
        for _,v,_,_ in out2[u]:
            indeg2[v]-=1
            if indeg2[v]==0: q.append(v)
    if len(topo)!=len(active): return {'cycle':True}
    pref={v:0 for v in active}; pref[a]=1
    for u in topo:
        for _,v,_,_ in out2[u]: pref[v]+=pref[u]
    suff={v:0 for v in active}; suff[b]=1
    for u in topo[::-1]:
        for _,v,_,_ in out2[u]: suff[u]+=suff[v]
    M=pref[b]; unique=all(min(pref[v],suff[v])==1 for v in active)
    E=len(eactive); V=len(active); dim=int(M-E+V-2) if M else 0
    return {'M':int(M),'E':E,'V':V,'dim':dim,'unique':bool(unique),'topo':topo,'pref':pref,'suff':suff,'edges':eactive}

def coarsen(edges):
    d=defaultdict(float)
    for u,v,k,w in edges: d[(u,v)]+=w
    return [(u,v,-1,w) for (u,v),w in d.items()]

def sample_path(edges,a,b,rng):
    out=defaultdict(list)
    for e in edges: out[e[0]].append(e)
    p=[]; u=a; safety=0
    while u!=b and safety<100:
        es=out[u]
        if not es: return None
        ws=np.array([x[3] for x in es],float); ws/=ws.sum()
        e=es[rng.choice(len(es),p=ws)]; p.append(e); u=e[1]; safety+=1
    return p if u==b else None

def route_bounds(path,edges):
    # path is list edges. vertex throughput from all edges at each internal node
    if not path: return (np.nan,np.nan)
    fs=[e[3] for e in path]
    outflow=defaultdict(float)
    for u,v,k,w in edges: outflow[u]+=w
    internal=[e[1] for e in path[:-1]]
    q=sum(fs)-sum(outflow[v] for v in internal)
    return max(0.0,q),min(fs)

rows=[]
settings=[]
for T in [5,10,20]:
  for rep in ['low','moderate','high']:
   for pm in [0.0,0.25]:
    for disp in ['balanced','unequal']:
     settings.append((T,rep,pm,disp))

N=100
for T,rep,pm,disp in settings:
    for s in range(N):
        studies=gen_design(T,rep,pm); gs=study_conductances(studies,disp)
        L=treatment_laplacian(studies,gs,T)
        # 5 random targets per network, or all if T=5 up to 10 targets
        pairs=list(itertools.combinations(range(T),2))
        nt=min(5,len(pairs)); idx=rng.choice(len(pairs),nt,replace=False)
        for ii in idx:
            a,b=pairs[ii]
            phi=solve_potential(L,a,b)
            all_edges=canonical_edges(studies,gs,phi,(a,b))
            ind=positive_indirect_edges(all_edges,(a,b))
            # if no indirect evidence, provenance has no indirect route issue; mark indirect_mass zero
            indirect_mass=sum(w for u,v,k,w in ind if u==a)
            if indirect_mass<1e-9:
                rows.append(dict(T=T,rep=rep,p_multi=pm,disp=disp,K=len(studies),target=f'{a}:{b}',indirect_mass=0.0,
                                 unique_study=True,dim_study=0,M_study=0,unique_treat=True,dim_treat=0,M_treat=0,
                                 sampled_width=0.0,study_edges=0,treat_edges=0))
                continue
            g=graph_counts(ind,a,b,T)
            cg=graph_counts(coarsen(ind),a,b,T)
            if g is None or 'cycle' in g or cg is None or 'cycle' in cg:
                continue
            # Route widths are evaluated in the dedicated exact series-design simulation.
            widths=[]
            rows.append(dict(T=T,rep=rep,p_multi=pm,disp=disp,K=len(studies),target=f'{a}:{b}',indirect_mass=indirect_mass,
                             unique_study=g['unique'],dim_study=g['dim'],M_study=g['M'],unique_treat=cg['unique'],dim_treat=cg['dim'],M_treat=cg['M'],
                             sampled_width=np.nan,study_edges=g['E'],treat_edges=cg['E']))

DF=pd.DataFrame(rows)
DF.to_csv(OUTPUT / 'structural_final_targets.csv',index=False)
summary=(DF.groupby(['T','rep','p_multi','disp'])
    .agg(n_targets=('target','size'),nonunique_study=('unique_study',lambda x:100*(~x).mean()),
         nonunique_treat=('unique_treat',lambda x:100*(~x).mean()),
         median_dim_nonunique=('dim_study',lambda x: float(np.median(x[np.array(x)>0])) if np.any(np.array(x)>0) else 0.0),
         median_sampled_width=('sampled_width','median'),p90_sampled_width=('sampled_width',lambda x: float(np.nanquantile(x,.9))),
         mean_indirect_mass=('indirect_mass','mean'))
    .reset_index())
summary.to_csv(OUTPUT / 'structural_final_summary.csv',index=False)
print('rows',len(DF))
print(summary.to_string(index=False))
# Aggregate concise across dispersion for main paper
main=(DF.groupby(['T','rep','p_multi'])
    .agg(n_targets=('target','size'),nonunique_study=('unique_study',lambda x:100*(~x).mean()),
         nonunique_treat=('unique_treat',lambda x:100*(~x).mean()),
         median_dim_nonunique=('dim_study',lambda x: float(np.median(x[np.array(x)>0])) if np.any(np.array(x)>0) else 0.0),
         median_sampled_width=('sampled_width','median'),p90_sampled_width=('sampled_width',lambda x: float(np.nanquantile(x,.9))))
    .reset_index())
main.to_csv(OUTPUT / 'structural_final_main.csv',index=False)
print('\nMAIN\n',main.to_string(index=False))
