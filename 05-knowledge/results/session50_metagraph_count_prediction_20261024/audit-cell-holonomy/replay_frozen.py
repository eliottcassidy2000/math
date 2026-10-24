#!/usr/bin/env python3
"""Replay exact frozen producers in normal/-O; no producer output is an input."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
ROOT=Path(__file__).resolve().parent

def require(value,label):
    if not value: raise RuntimeError(label)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--receipt',required=True); args=parser.parse_args()
    reports={}
    with tempfile.TemporaryDirectory(prefix='cell-holonomy-producer-replay-') as temp:
        directory=Path(temp)
        for name,script in [('holonomy','source/holonomy/check.py.txt'),('predictor','source/predictor/predictor.py.txt')]:
            runs=[]
            for optimized in (False,True):
                work=directory/(name+('-O' if optimized else '')); work.mkdir()
                receipt=work/'receipt.json'
                command=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(ROOT/script)]
                if name=='predictor': command+=['--parents',str(ROOT/'source/predictor/parents-through-seven.json'),'--receipt',str(receipt)]
                result=subprocess.run(command,cwd=work,text=False,capture_output=True)
                require(result.returncode==0,name+' replay failed: '+result.stderr.decode())
                runs.append((result.stdout,receipt.read_bytes()))
            require(runs[0]==runs[1],name+' normal/optimized mismatch')
            stdout,raw=runs[0]; receipt=json.loads(raw)
            if name=='holonomy':
                require(raw==(ROOT/'source/holonomy/receipt.json').read_bytes(),'frozen holonomy receipt mismatch')
                require(stdout==(ROOT/'source/holonomy/normal.out').read_bytes(),'frozen holonomy stdout mismatch')
            else:
                require((receipt['pointed_extension_states'],receipt['upper_refinement_trace'],receipt['SCC_total_error_order_eight'],receipt['predicted_A000568_order_eight'])==(54256,55068,28,6880),'frozen predictor output mismatch')
                require(receipt['child_canonicalization_calls']==receipt['deletion_or_reroot_isomorphism_calls']==0,'frozen predictor boundary')
            reports[name]={'status':'PASS','normal_optimized_identical':True,'receipt_sha256':hashlib.sha256(raw).hexdigest(),'stdout_sha256':hashlib.sha256(stdout).hexdigest()}
    Path(args.receipt).write_text(json.dumps(reports,sort_keys=True,indent=2)+'\n')
    print(json.dumps(reports,sort_keys=True))
if __name__=='__main__': main()
