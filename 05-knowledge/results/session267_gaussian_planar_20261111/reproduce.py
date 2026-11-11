#!/usr/bin/env python3
"""Replay only this packet's new finite controls, including Python -O."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
SCRIPTS=[
    'gaussian/gaussian_bridge_exact.py',
    'planar/check.py',
    'planar/check_extension.py',
    'topology/permutation_recursion_controls.py',
    'degree/check.py',
    'degree/state_check.py',
]

def main():
    rows=[]
    for script in SCRIPTS:
        results=[]
        for flags in ([],['-O']):
            result=subprocess.run([sys.executable,*flags,script],cwd=ROOT,
                                  capture_output=True,text=True)
            if result.returncode:
                raise RuntimeError(f'{script} {flags}: {result.stderr}')
            results.append(result.stdout)
        if results[0]!=results[1]:
            raise RuntimeError(f'Optimized output differs for {script}')
        rows.append({'script':script,'normal_exit':0,'optimized_exit':0,
                     'outputs_identical':True,
                     'stdout_sha256':hashlib.sha256(results[0].encode()).hexdigest(),
                     'stdout':results[0]})
    import sympy
    print(json.dumps({'status':'PASS','python':sys.version.split()[0],
                      'sympy':sympy.__version__,'programs':rows,
                      'scope':'New finite controls only; no inherited proof rebuild or repository-wide test.'},
                     indent=2))

if __name__=='__main__':main()
