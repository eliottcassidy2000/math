#!/usr/bin/env python3
"""Replay Session256 new claims in a disposable copy; preserve saved receipts."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,tempfile

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--optimized',action='store_true');args=ap.parse_args()
    root=Path(__file__).resolve().parent
    jobs=[('history/probe_bank.py','history/probe_bank.out','stdout'),
          ('spinor/test_spinor.py','spinor/results.json','stdout'),
          ('spinor/test_sorted.py','spinor/sorted_results.json','stdout'),
          ('gram/check.py','gram/results.json','stdout'),
          ('compiler/engine.py','compiler/results.json','file'),
          ('compiler/morphic_optimizer.py','compiler/morphic_results.json','file'),
          ('compiler/verify_morphic_certificate.py','compiler/morphic_certificate_verified.json','file')]
    report={'scope':'Session256 new-claim engines only','optimized':args.optimized,'engines':[]}
    with tempfile.TemporaryDirectory(prefix='session256-replay-') as tmp:
        work=Path(tmp)/'packet';shutil.copytree(root,work,ignore=shutil.ignore_patterns('__pycache__'))
        for script,receipt,kind in jobs:
            path=work/script
            command=[sys.executable,'-B']+(['-O'] if args.optimized else [])+[str(path)]
            p=subprocess.run(command,cwd=path.parent,text=True,capture_output=True)
            if p.returncode:raise RuntimeError(script+' failed: '+p.stderr)
            actual=p.stdout.encode() if kind=='stdout' else (work/receipt).read_bytes()
            expected=(root/receipt).read_bytes()
            if json.loads(actual)!=json.loads(expected):raise ValueError(script+' receipt mismatch')
            report['engines'].append({'script':script,'receipt':receipt,'status':'PASS',
              'source_sha256':hashlib.sha256((root/script).read_bytes()).hexdigest(),
              'receipt_sha256':hashlib.sha256(expected).hexdigest()})
    report['status']='PASS';print(json.dumps(report,indent=2))
if __name__=='__main__':main()
