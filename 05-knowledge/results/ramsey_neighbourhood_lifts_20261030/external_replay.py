"""Replay the pinned author's exact checkers; does not fetch or modify them."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys

HERE=Path(__file__).resolve().parent
PIN='5770baac746878a3a63b9720fdda97726fe02851'
def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--pdf',type=Path,required=True);p.add_argument('--write',action='store_true');args=p.parse_args()
    src=args.source.resolve()
    got=subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip()
    if got!=PIN:raise RuntimeError('external source revision differs')
    env=dict(os.environ);env.pop('PYTHONOPTIMIZE',None)
    jobs=[('hand','ancillary',['verify_hand_proofs.py']),
          ('flag_chain','ancillary/flag_r33k',['check_r33k5.py','certificates/five/k6/excl_86.json','certificates/five/k7/excl_125.json','certificates/five/k8/excl_174.json'])]
    runs={}
    for name,directory,argv in jobs:
        result=subprocess.run([sys.executable,'-X','utf8','-B',*argv],cwd=src/directory,env=env,capture_output=True,text=True,timeout=120)
        if result.returncode:raise RuntimeError(name+' failed: '+result.stdout+result.stderr)
        if name=='hand' and '17985' not in result.stdout:raise RuntimeError('unexpected hand receipt')
        if name=='flag_chain' and result.stdout.count('VALID')<3:raise RuntimeError('missing exact flag approvals')
        runs[name]={'arguments':argv,'directory':directory,'stdout':result.stdout,'stderr':result.stderr,'exit_code':result.returncode}
    paths=['paper/neighbourhood_counting.tex','ancillary/verify_hand_proofs.py','ancillary/flag_r33k/check_r33k.py','ancillary/flag_r33k/check_r33k5.py']
    paths += ['ancillary/flag_r33k/'+x for x in jobs[1][2][1:]]
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    out={'source_url':'https://github.com/Lachy-Dauth/ramsey-neighbourhood-counting','commit':PIN,'source_sha256':{x:digest(src/x) for x in paths},
         'input_pdf':{'filename':args.pdf.name,'bytes':args.pdf.stat().st_size,'sha256':digest(args.pdf)},'runs':runs,
         'scope':'Author checkers in normal Python; not a Lean build, not an independent reimplementation. R338<=174 replayed; 173 is cited only.'}
    if args.write:(HERE/'external_receipt.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    else:
        if out!=json.loads((HERE/'external_receipt.json').read_text()):raise RuntimeError('source or replay mismatch')
    print(json.dumps({'status':'PASS','pin':PIN,'hand_checks':17985,'flag_exclusions':[86,125,174]}))
if __name__=='__main__':main()
