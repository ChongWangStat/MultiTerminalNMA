#!/usr/bin/env python3
"""Refit both recovered R audits in separately supplied, version-checked libraries.
Does not install packages or alter a user library. See R_ENVIRONMENTS_v09.md.
"""
from pathlib import Path
import argparse,os,shutil,subprocess,json
ROOT=Path(__file__).resolve().parent

def run(out):
 rscript=shutil.which('Rscript')
 if not rscript: raise RuntimeError('Rscript is not installed; archived evidence is not a new refit.')
 specifications=[('covid','CSP_NETMETA_R_LIB',{'netmeta':'3.7.0'},'covid_netmeta_compare_v4.R'),('nmaforest','CSP_NMAFOREST_R_LIB',{'netmeta':'3.2.0','igraph':'2.1.4','NMAforest':'0.1.3'},'nmaforest_rjournal_reproduce.R')]
 for name,var,versions,source in specifications:
  folder=(out/name).resolve();folder.mkdir(parents=True,exist_ok=True)
  lib=os.environ.get(var)
  if not lib or not Path(lib).is_dir():raise RuntimeError(f'Set {var} to the prepared, separate R library.')
  lines=[f'.libPaths(c({json.dumps(lib)}, .libPaths()))','stopifnot(as.character(getRversion()) == "4.4.1")']
  lines += [f'stopifnot(as.character(packageVersion("{pkg}")) == "{version}")' for pkg,version in versions.items()]
  lines += [f'source({json.dumps(str(ROOT/source))})']
  if name=='covid':
   lines += ['saveRDS(list(nma=nma,Hlong=Hlong,Hfull=Hfull,phi_L1=cm_l1,phi_L2=cm_pi,netcontrib=nc),"covid_raw_objects.rds")',f'source({json.dumps(str(ROOT/"covid_tau_probe.R"))})']
  lines += ['writeLines(capture.output(sessionInfo()),"sessionInfo.txt")','write.csv(installed.packages()[,c("Package","Version")],"installed_packages.csv",row.names=FALSE)']
  script=folder/'run.R';script.write_text('\n'.join(lines)+'\n')
  with (folder/'execution.log').open('w') as log:
   subprocess.run([rscript,str(script)],cwd=folder,stdout=log,stderr=subprocess.STDOUT,check=True)
 print('Both R environments executed successfully.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,default=ROOT/'reproduced/r_audits');run(p.parse_args().output_dir)
