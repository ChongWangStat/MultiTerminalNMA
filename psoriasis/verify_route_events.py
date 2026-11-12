"""Independent, enumerated-route LP checks; exhaustive small supports + random DAGs."""
from pathlib import Path
from itertools import combinations, product
import argparse, json, time
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from route_events import (Edge, Automaton, visit_groups, ordered_groups,
                          contains_segment, complement, sharp_event_bounds)


def enumerate_routes(edges, source, sink):
    out = {}
    for i,e in enumerate(edges): out.setdefault(e.tail,[]).append(i)
    paths = []
    def walk(v,p):
        if v == sink: paths.append(tuple(p)); return
        for i in out.get(v,[]): walk(edges[i].head,p+[i])
    walk(source,[])
    return paths


def incidence(edges,paths):
    a = np.zeros((len(edges),len(paths)))
    for j,p in enumerate(paths): a[list(p),j] = 1.
    return a


def check(edges,s,t,event,predicate,case,rows):
    paths = enumerate_routes(edges,s,t)
    a = incidence(edges,paths); f = np.array([e.flow for e in edges])
    c = np.array([predicate(p,edges) for p in paths],float)
    lo = linprog(c,A_eq=a,b_eq=f,bounds=(0,None),method='highs')
    hi = linprog(-c,A_eq=a,b_eq=f,bounds=(0,None),method='highs')
    if not lo.success or not hi.success: raise AssertionError('Enumerated LP failure')
    ans = sharp_event_bounds(edges,s,t,event,witnesses=True)
    err = max(abs(lo.fun-ans['lower']),abs(-hi.fun-ans['upper']))
    witness_err = 0.
    for endpoint in ('lower','upper'):
        projection = np.zeros(len(edges)); mass = 0.
        for witness in ans[endpoint+'_witness']:
            p = witness['edges']; weight = witness['mass']
            projection[p] += weight
            if predicate(p,edges): mass += weight
        witness_err = max(witness_err,float(np.max(np.abs(projection-f))),abs(mass-ans[endpoint]))
    cert = max(v for side in ('lower','upper') for v in ans[side+'_certificate'].values())
    if max(err,witness_err,cert) > 2e-8:
        raise AssertionError((case,event.name,err,witness_err,cert))
    rows.append(dict(case=case,event=event.name,edges=len(edges),routes=len(paths),
                     lower=ans['lower'],upper=ans['upper'],max_endpoint_error=err,
                     witness_error=witness_err,certificate_error=cert,
                     product_variables=ans['variables']))


def queries():
    any_a = lambda p,e:any(e[i].label=='a' for i in p)
    both = lambda p,e: all(any(e[i].label==lab for i in p) for lab in ('a','b'))
    ordered = lambda p,e:any(e[p[j]].label=='a' and e[p[k]].label=='b'
                            for j in range(len(p)) for k in range(j+1,len(p)))
    return [(visit_groups([{'a'}]),any_a),
            (visit_groups([{'a'},{'b'}]),both),
            (complement(visit_groups([{'a'},{'b'}])),lambda p,e:not both(p,e)),
            (ordered_groups([{'a'},{'b'}]),ordered),
            (Automaton(2,0,frozenset({1}),lambda q,i,e:q ^ int(e.label=='a'),'odd_a'),
             lambda p,e:sum(e[i].label=='a' for i in p)%2==1)]


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    rows = []; segments = []; graph_count = 0; labeled_count = 0
    candidates = list(combinations(range(4),2))
    for mask in range(1,1<<len(candidates)):
        pairs = [edge for j,edge in enumerate(candidates) if mask & (1<<j)]
        raw = [Edge(u,v,'a',1.) for u,v in pairs]
        paths = enumerate_routes(raw,0,3)
        if not paths or set(i for p in paths for i in p) != set(range(len(raw))): continue
        graph_count += 1
        a = incidence(raw,paths)
        pi = np.arange(1,len(paths)+1,dtype=float); pi /= pi.sum()
        f = a@pi
        for labels in product(('a','b'),repeat=len(raw)):
            labeled_count += 1
            edges = [Edge(u,v,lab,float(w)) for (u,v),lab,w in zip(pairs,labels,f)]
            for event,pred in queries():
                check(edges,0,3,event,pred,f'exhaustive_{mask}_{labeled_count}',rows)
        edges = [Edge(u,v,'a',float(w)) for (u,v),w in zip(pairs,f)]
        segs = sorted({p[j:k] for p in paths for j in range(len(p)) for k in range(j+1,len(p)+1)})
        throughput = {v:sum(e.flow for e in edges if e.tail==v) for v in range(4)}
        for segment in segs:
            pred = lambda p,e,sg=segment:any(tuple(p[j:j+len(sg)])==sg for j in range(len(p)-len(sg)+1))
            check(edges,0,3,contains_segment(segment),pred,f'segment_{mask}',segments)
            formula_l = max(0.,sum(f[i] for i in segment)-sum(throughput[edges[i].head] for i in segment[:-1]))
            formula_u = min(f[i] for i in segment)
            segments[-1]['formula_error'] = max(abs(formula_l-segments[-1]['lower']),abs(formula_u-segments[-1]['upper']))
            assert segments[-1]['formula_error'] < 2e-8
    rng = np.random.default_rng(20261004)
    for rep in range(100):
        n = int(rng.integers(5,9)); pairs = [(i,i+1) for i in range(n-1)]
        pairs += [(i,j) for i in range(n) for j in range(i+1,n) if rng.random()<0.25]
        raw = [Edge(u,v,str(rng.choice(['a','b','c'])),1.) for u,v in pairs]
        paths = enumerate_routes(raw,0,n-1); a = incidence(raw,paths)
        pi = rng.dirichlet(np.ones(len(paths))); pi *= 0.37 if rep%2 else 1.
        f = a@pi; edges = [Edge(e.tail,e.head,e.label,float(w)) for e,w in zip(raw,f)]
        for event,pred in queries(): check(edges,0,n-1,event,pred,f'random_{rep}',rows)
        for accepted in (False,True):
            ev = Automaton(1,0,frozenset({0}) if accepted else frozenset(),lambda q,i,e:0,'all' if accepted else 'none')
            check(edges,0,n-1,ev,lambda p,e,z=accepted:z,f'random_{rep}',rows)
    # Full route enumeration intentionally absent from this scalability check.
    stress = []
    for layers in (10,30,60):
        edges = [Edge(i,i+1,'a' if i==0 and j==0 else 'b' if i==layers-1 and j==0 else 'c',0.5)
                 for i in range(layers) for j in (0,1)]
        start=time.perf_counter()
        result = sharp_event_bounds(edges,0,layers,visit_groups([{'a'},{'b'}]),witnesses=True)
        elapsed=time.perf_counter()-start
        assert abs(result['lower'])<1e-9 and abs(result['upper']-0.5)<1e-9
        stress.append(dict(layers=layers,routes=str(2**layers),edges=len(edges),
                           product_variables=result['variables'],equalities=result['equalities'],
                           lower=result['lower'],upper=result['upper'],seconds=elapsed))
    pd.DataFrame(rows).to_csv(out/'route_event_verification.csv',index=False)
    pd.DataFrame(segments).to_csv(out/'route_segment_verification.csv',index=False)
    pd.DataFrame(stress).to_csv(out/'route_event_scalability.csv',index=False)
    summary = dict(seed=20261004,exhaustive_supports=graph_count,exhaustive_labelings=labeled_count,
                   event_queries=len(rows),segment_queries=len(segments),random_multigraphs=100,
                   max_endpoint_error=max(r['max_endpoint_error'] for r in rows+segments),
                   max_witness_error=max(r['witness_error'] for r in rows+segments),
                   max_certificate_error=max(r['certificate_error'] for r in rows+segments),
                   max_segment_formula_error=max(r['formula_error'] for r in segments),
                   scalability_cases=len(stress))
    (out/'route_event_verification_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--output-dir',type=Path,default=Path('verification'))
    run(parser.parse_args().output_dir)
