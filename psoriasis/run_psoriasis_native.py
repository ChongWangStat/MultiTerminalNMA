"""Execute the unchanged supplied R computational source and export its fitted objects."""
from pathlib import Path
import argparse,os,shutil,subprocess
ROOT=Path(__file__).resolve().parent


def run(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=True)
    raw=os.environ.get('CSP_R_EXECUTABLE')
    if raw:
        command=[raw,'--no-save','--no-restore','--quiet','--file='+str(ROOT/'run_psoriasis_csp.R'),'--args',str(out)]
    else:
        rscript=shutil.which('Rscript')
        if not rscript:raise RuntimeError('Install R and the packages listed in R_ENVIRONMENTS_v10.md, or specify CSP_R_EXECUTABLE.')
        command=[rscript,str(ROOT/'run_psoriasis_csp.R'),str(out)]
    with (out/'execution.log').open('w') as log:
        subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    print('Native R common and random psoriasis fits completed; fitted objects and session information exported.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,default=ROOT/'reproduced/psoriasis_R');run(p.parse_args().output_dir)
