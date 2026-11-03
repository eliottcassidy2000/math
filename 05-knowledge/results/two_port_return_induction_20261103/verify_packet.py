#!/usr/bin/env python3
"""Replay frozen exact jobs; optionally reuse this session's explicit root stdout."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
JOBS=['rooted_resets','mixed_returns','observer_induction_probe',
      'observer_induction_basis','observer_induction_reset_audit','native_reset_induction_probe']

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reuse-root-normal',type=Path)
    args=ap.parse_args();results=[]
    for name in JOBS:
        expected=(HERE/(name+'.json')).read_bytes()
        for opt in [False,True]:
            if name=='rooted_resets' and not opt and args.reuse_root_normal:
                actual=args.reuse_root_normal.read_bytes()
            else:
                command=[sys.executable,'-I','-B']+(['-O'] if opt else [])+[str(HERE/(name+'.py'))]
                actual=subprocess.check_output(command,cwd=HERE)
            if actual!=expected:raise RuntimeError(f'{name}: frozen mismatch, optimized={opt}')
        results.append({'script':name+'.py','receipt':name+'.json',
                        'script_sha256':hashlib.sha256((HERE/(name+'.py')).read_bytes()).hexdigest(),
                        'receipt_sha256':hashlib.sha256(expected).hexdigest(),
                        'normal_exact':True,'optimized_exact':True})
        print('PASS '+name,flush=True)
    (HERE/'verification.json').write_text(json.dumps({'status':'PASS',
        'scope':'Six frozen jobs reproduced in normal and optimized Python; no expanded scan universe',
        'jobs':results},indent=2,sort_keys=True)+'\n')

if __name__=='__main__':main()
