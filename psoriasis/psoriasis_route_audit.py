"""All 45 psoriasis targets, sharp route bounds, R comparison, and covariance sensitivity."""
from pathlib import Path
from itertools import combinations
from collections import defaultdict
import argparse, json, hashlib
import numpy as np
import pandas as pd
from scipy.optimize import linprog, minimize_scalar
from psoriasis_csp import ROOT, read_data, prepare, fit, estimate_tau2, independent_gls, endpoint_edges
from covid_route_audit import enumerate_paths, coarsen_edges, build as legacy_build, legacy_endpoint_edges as legacy_edges
from route_events import Edge, validate_flow, visit_groups, sharp_event_bounds
from verify_route_events import check


def route_audit(raw,source,sink,model,target,level,route_rows,check_rows):
    if not raw:
        return dict(edges=0,paths=0,rank=0,nullity=0,unique=True,max_width=0.)
    edges=[Edge(e['fr'],e['to'],e['study'],e['w']) for e in raw]
    order,adj,mass=validate_flow(edges,source,sink)
    paths=enumerate_paths(raw,source,sink); A=np.zeros((len(edges),len(paths)))
    f=np.array([e.flow for e in edges]);F={v:sum(edges[i].flow for i in adj[v]) for v in order}
    for j,p in enumerate(paths):A[list(p),j]=1
    rank=len(edges)-len(order)+2
    assert np.linalg.matrix_rank(A,tol=1e-9)==rank
    maximum=0.
    for j,p in enumerate(paths):
        upper=min(edges[i].flow for i in p)
        lower=max(0.,sum(edges[i].flow for i in p)-sum(F[edges[i].head] for i in p[:-1]))
        assert lower<=upper+1e-9
        lower=min(lower,upper);width=max(0.,upper-lower);maximum=max(maximum,width)
        c=np.zeros(len(paths));c[j]=1
        lo=linprog(c,A_eq=A,b_eq=f,bounds=(0,None),method='highs')
        hi=linprog(-c,A_eq=A,b_eq=f,bounds=(0,None),method='highs')
        if not (lo.success and hi.success):raise RuntimeError((model,target,level,j))
        error=max(abs(lo.fun-lower),abs(-hi.fun-upper))
        assert error<2e-8
        check_rows.append(dict(model=model,target=target,level=level,route_id=j+1,endpoint_error=error))
        route_rows.append(dict(model=model,target=target,level=level,route_id=j+1,
                               treatment_route='->'.join([source]+[edges[i].head for i in p]),
                               study_route=' ||| '.join(str(edges[i].label) for i in p),
                               lower=lower,upper=upper,width=width,indirect_mass=mass))
    return dict(edges=len(edges),paths=len(paths),rank=rank,nullity=len(paths)-rank,
                unique=len(paths)==rank,max_width=maximum)


def edge_map(edges):
    d=defaultdict(float)
    for e in edges:d[(e['id'],e['fr'],e['to'])]+=e['w']
    return d


def run(out,r_reference=None):
    out.mkdir(parents=True,exist_ok=True);dat=read_data();q=prepare(dat)
    python_tau,tau_info=estimate_tau2(q)
    r_present=r_reference is not None and (r_reference/'psoriasis_fit_random.csv').exists()
    tau=float(pd.read_csv(r_reference/'psoriasis_fit_random.csv').tau2.iloc[0]) if r_present else python_tau
    independent_tau=minimize_scalar(lambda t:independent_gls(q,t)[1],bounds=(0,tau_info['upper_search_bound']),
                                     method='bounded',options={'xatol':1e-10}).x
    assert abs(independent_tau-python_tau)<1e-6
    if r_present:assert abs(tau-python_tau)<1e-4
    dat.to_csv(out/'psoriasis_formatted_arm_data.csv',index=False)
    counts=dat.groupby(['id','study']).agg(arms=('t','size'),treatments=('t',lambda x:':'.join(sorted(x))))
    counts['all_pair_contrasts']=counts.arms*(counts.arms-1)//2
    counts['independent_contrasts']=counts.arms-1
    counts.to_csv(out/'psoriasis_study_inventory.csv')
    rows=[];all_edges=[];route_rows=[];check_rows=[];summaries=[];verification={}; fits={}
    for model,t in [('common',0.),('random',tau)]:
        result=fit(q,t);fits[model]=result
        theta_ref,obj_ref=independent_gls(q,t)
        gls_error=float(np.max(np.abs(result['theta']-theta_ref)))
        assert gls_error<1e-10
        old=legacy_build(dat,tau2=t);max_old_error=0.;max_reconstruction=0.;max_balance=0.
        target_rows=[];model_edges=[]
        for ti,tp in enumerate(q['pair_names']):
            s,k=tp.split(':');ee=endpoint_edges(result,ti)
            direct=[e for e in ee if e['type']=='Dir'];indirect=[e for e in ee if e['type']=='Ind']
            dmass=sum(e['w'] for e in direct);imass=sum(e['w'] for e in indirect if e['fr']==s)
            assert abs(dmass+imass-1)<1e-10
            positive=imass>1e-10
            a=route_audit(indirect,s,k,model,tp,'study',route_rows,check_rows)
            ca=route_audit(coarsen_edges(indirect),s,k,model,tp,'treatment',route_rows,check_rows)
            reconstruction=0.
            for e in ee:
                b=next(b for b in q['blocks'] if b['id']==e['id'])
                reconstruction+=e['w']*(b['z'][b['arms'].index(e['to'])]-b['z'][b['arms'].index(e['fr'])])
            max_reconstruction=max(max_reconstruction,abs(reconstruction-result['theta'][ti]))
            balance=defaultdict(float)
            for e in ee:balance[e['fr']]+=e['w'];balance[e['to']]-=e['w']
            max_balance=max(max_balance,max(abs(v-(1 if node==s else -1 if node==k else 0)) for node,v in balance.items()))
            _,_,_,_,previous,_=legacy_edges(old,ti);x=edge_map(indirect);y=edge_map(previous)
            if x or y:max_old_error=max(max_old_error,max(abs(x[z]-y[z]) for z in x.keys()|y.keys()))
            row=dict(model=model,target=tp,tau2=t,theta=result['theta'][ti],direct_weight=dmass,indirect_weight=imass,
                     has_indirect=positive,status='positive_indirect_flow' if positive else 'direct_only_no_indirect_flow')
            row.update({'study_'+key:value for key,value in a.items()})
            row.update({'treatment_'+key:value for key,value in ca.items()})
            target_rows.append(row)
            model_edges.extend([dict(model=model,tau2=t,**e) for e in ee])
        df=pd.DataFrame(target_rows);ef=pd.DataFrame(model_edges)
        df.to_csv(out/f'psoriasis_route_audit_{model}.csv',index=False)
        ef.to_csv(out/f'psoriasis_endpoint_edges_{model}.csv',index=False)
        rows.extend(target_rows);all_edges.extend(model_edges)
        positive=df[df.has_indirect]
        summaries.append(dict(model=model,tau2=t,total_targets=len(df),positive_indirect_targets=len(positive),
                              direct_only_targets=int((~df.has_indirect).sum()),
                              study_unique_positive=int(positive.study_unique.sum()),
                              treatment_unique_positive=int(positive.treatment_unique.sum()),
                              study_unique_all=int(df.study_unique.sum()),treatment_unique_all=int(df.treatment_unique.sum()),
                              study_max_width=float(positive.study_max_width.max()),
                              treatment_max_width=float(positive.treatment_max_width.max()),
                              max_study_width_target=positive.loc[positive.study_max_width.idxmax(),'target'],
                              df_routes=int(df.set_index('target').loc['D:F','study_paths'])))
        verification[model]=dict(nonredundant_gls_max_error=gls_error,reml_objective_offset=result['objective']-obj_ref,
                                 endpoint_reconstruction_error=max_reconstruction,flow_balance_error=max_balance,
                                 earlier_vs_supplied_endpoint_max_error=max_old_error)
        if r_present:
            rfit=pd.read_csv(r_reference/f'psoriasis_fit_{model}.csv').set_index('target').loc[q['pair_names']]
            rp=pd.read_csv(r_reference/f'psoriasis_P_{model}.csv').set_index('target').loc[q['pair_names']]
            names=[nm for b in q['blocks'] for nm in b['row_names']]
            p_error=float(np.max(np.abs(rp[names].to_numpy()-result['P'])))
            theta_error=float(np.max(np.abs(rfit.theta.to_numpy()-result['theta'])))
            redges=pd.read_csv(r_reference/f'psoriasis_endpoint_edges_{model}.csv').rename(columns={'from':'fr','weight':'w'})
            keys=['target','id','fr','to','type']
            rvalues=redges.set_index(keys).w.sort_index();pyvalues=ef.set_index(keys).w.sort_index()
            assert rvalues.index.equals(pyvalues.index)
            e_error=float(np.max(np.abs(rvalues-pyvalues)))
            assert max(p_error,theta_error,e_error)<1e-10
            verification[model]['native_R_comparison']=dict(P_max_error=p_error,theta_max_error=theta_error,edge_max_error=e_error)
        assert max(max_old_error,max_reconstruction,max_balance)<1e-10
        print(model,summaries[-1],flush=True)
    routes=pd.DataFrame(route_rows);routes.to_csv(out/'psoriasis_complete_route_bounds.csv',index=False)
    pd.DataFrame(check_rows).to_csv(out/'psoriasis_route_bound_verification.csv',index=False)
    pd.DataFrame(summaries).to_csv(out/'psoriasis_summary.csv',index=False)
    common=pd.DataFrame(rows).query("model=='common'");random=pd.DataFrame(rows).query("model=='random'")
    assert list(common.loc[~common.has_indirect,'target'])==['A:D','A:G','A:J']
    assert summaries[0]['study_unique_positive']==5 and summaries[0]['treatment_unique_positive']==30
    assert abs(summaries[0]['study_max_width']-.4181)<5e-5
    assert abs(summaries[0]['treatment_max_width']-.1978)<5e-5
    # Illustrative event queries are chosen by the widest common-effect D:F route.
    widest=routes.query("model=='common' and target=='D:F' and level=='study'").sort_values(['width','route_id'],ascending=[False,True]).iloc[0]
    labels=widest.study_route.split(' ||| ')
    queries=[('visits_'+labels[0],[{labels[0]}]),('visits_'+labels[-1],[{labels[-1]}]),
             ('visits_both_'+labels[0]+'_'+labels[-1],[{labels[0]},{labels[-1]}])]
    event_rows=[];event_checks=[];witnesses={}
    for model in ['common','random']:
        result=fits[model];ti=q['pair_names'].index('D:F');raw=[e for e in endpoint_edges(result,ti) if e['type']=='Ind']
        edges=[Edge(e['fr'],e['to'],e['study'],e['w']) for e in raw]
        for name,groups in queries:
            event=visit_groups(groups);pred=lambda p,e,gg=groups:all(any(e[j].label in g for j in p) for g in gg)
            check(edges,'D','F',event,pred,model+'_'+name,event_checks)
            ans=sharp_event_bounds(edges,'D','F',event,witnesses=True)
            event_rows.append(dict(model=model,target='D:F',event=name,lower=ans['lower'],upper=ans['upper']))
            witnesses[model+'_'+name]={'edges':[e.__dict__ for e in edges],**ans}
    pd.DataFrame(event_rows).to_csv(out/'psoriasis_DF_event_bounds.csv',index=False)
    pd.DataFrame(event_checks).to_csv(out/'psoriasis_event_verification.csv',index=False)
    (out/'psoriasis_DF_event_witnesses.json').write_text(json.dumps(witnesses,indent=2)+'\n')
    verification.update(dict(input_sha256=hashlib.sha256((ROOT/'pso-psa-data.csv').read_bytes()).hexdigest(),
                             route_queries_checked=len(check_rows),max_route_endpoint_error=max(r['endpoint_error'] for r in check_rows),
                             event_queries_checked=len(event_checks),event_query_selection='Illustrative, selected using widest common-effect D:F route; no inferential coverage claim.',
                             widest_DF_route=widest.to_dict(),python_reml=tau_info,independent_basis_tau2=float(independent_tau),
                             primary_tau2=tau,native_R_reference_compared=bool(r_present),
                             target_reconciliation='45 pairs; A:D, A:G and A:J have zero indirect flow; 42 positive-indirect targets reproduce all archived summary values.',
                             direct_only_targets=['A:D','A:G','A:J'],
                             support_classification_changes=int(np.sum(common.study_unique.to_numpy()!=random.study_unique.to_numpy())),
                             passed=True))
    (out/'psoriasis_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    print('All route checks passed:',len(check_rows),'max error',verification['max_route_endpoint_error'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,default=ROOT/'reproduced/psoriasis')
    p.add_argument('--r-reference-dir',type=Path,default=ROOT/'psoriasis_results/R')
    args=p.parse_args();run(args.output_dir,args.r_reference_dir)
