"""Record the effect of using the supplied CSP endpoint rule instead of the v09 rule."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from covid_route_audit import covid,build,endpoint_edges,legacy_endpoint_edges
from psoriasis_route_audit import edge_map


def run(out,r_reference=None):
    root=Path(__file__).resolve().parent;out.mkdir(parents=True,exist_ok=True)
    rows=[];r_checks={}
    for model,tau2 in [('common',0.),('random',.002103984)]:
        old=pd.read_csv(root/'legacy_v09'/f'covid_route_audit_{model}.csv').set_index('target')
        new=pd.read_csv(root/f'covid_route_audit_{model}.csv').set_index('target')
        f=build(covid,tau2)
        if r_reference is not None and (r_reference/f'covid_CSP_fit_{model}.csv').exists():
            from psoriasis_csp import endpoint_edges as current_edges
            native_fit=pd.read_csv(r_reference/f'covid_CSP_fit_{model}.csv').set_index('target').loc[f['pair_names']]
            native_P=pd.read_csv(r_reference/f'covid_CSP_P_{model}.csv').set_index('target').loc[f['pair_names']]
            p_error=float(np.max(np.abs(native_P.to_numpy()-f['P'])))
            theta_error=float(np.max(np.abs(native_fit.theta.to_numpy()-f['theta'])))
            ee=pd.DataFrame([e for i in range(len(f['pair_names'])) for e in current_edges(f,i,tol=1e-10)])
            rr=pd.read_csv(r_reference/f'covid_CSP_edges_{model}.csv').rename(columns={'from':'fr','weight':'w'})
            keys=['target','id','fr','to','type'];x=ee.set_index(keys).w.sort_index();y=rr.set_index(keys).w.sort_index()
            assert x.index.equals(y.index)
            edge_error=float(np.max(np.abs(x-y)))
            assert max(p_error,theta_error,edge_error)<1e-10
            r_checks[model]=dict(P_max_error=p_error,theta_max_error=theta_error,edge_max_error=edge_error)
        for i,target in enumerate(f['pair_names']):
            x=edge_map(endpoint_edges(f,i)[4]);y=edge_map(legacy_endpoint_edges(f,i)[4])
            e=max(abs(x[k]-y[k]) for k in x.keys()|y.keys())
            rows.append(dict(model=model,target=target,old_study_paths=int(old.loc[target,'study_paths']),
                             supplied_study_paths=int(new.loc[target,'study_paths']),max_endpoint_flow_change=e,
                             fitted_NMA_change=abs(old.loc[target,'theta']-new.loc[target,'theta']),
                             study_width_change=abs(old.loc[target,'study_max_width']-new.loc[target,'study_max_width']),
                             treatment_width_change=abs(old.loc[target,'treatment_max_width']-new.loc[target,'treatment_max_width']),
                             study_unique_unchanged=bool(old.loc[target,'study_unique']==new.loc[target,'study_unique']),
                             treatment_unique_unchanged=bool(old.loc[target,'treatment_unique']==new.loc[target,'treatment_unique'])))
    df=pd.DataFrame(rows);df.to_csv(out/'covid_CSP_version_comparison.csv',index=False)
    assert df.study_unique_unchanged.all() and df.treatment_unique_unchanged.all()
    assert df[['fitted_NMA_change','study_width_change','treatment_width_change']].max().max()<1e-10
    report=dict(passed=True,scope='Current supplied CSP endpoint rule is now used for both applications; earlier COVID outputs retained in legacy_v09.',
                changed_route_counts=df.loc[df.old_study_paths!=df.supplied_study_paths,['model','target','old_study_paths','supplied_study_paths']].to_dict('records'),
                max_endpoint_change=float(df.max_endpoint_flow_change.max()),
                max_fitted_NMA_change=float(df.fitted_NMA_change.max()),
                max_route_width_change=float(df[['study_width_change','treatment_width_change']].max().max()),
                native_R_comparisons=r_checks)
    (out/'covid_CSP_version_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,default=Path('verification'))
    p.add_argument('--r-reference-dir',type=Path,default=Path(__file__).resolve().parent/'psoriasis_results/R')
    args=p.parse_args();run(args.output_dir,args.r_reference_dir)
