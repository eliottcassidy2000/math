#!/usr/bin/env python3
"""Replay all portable producer and independent exact checks, including -O."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
CASES=[
 ('graft_images/verify.py','graft_images/CHECK.json'),
 ('image_audit/independent_check.py','image_audit/independent.json'),
 ('multiroot/verify.py','multiroot/CHECK.json'),
 ('construction_audit/independent_multiroot_audit.py','construction_audit/normal.json'),
 ('boundary_escape/verify.py','boundary_escape/CHECK.json'),
 ('escape_audit/independent_check.py','escape_audit/INDEPENDENT.json'),
 ('menu_boundary/verify.py','menu_boundary/CHECK.json'),
 ('menu_audit/independent_check.py','menu_audit/independent.json'),
 ('adaptive_arithmetic/verify.py','adaptive_arithmetic/CHECK.json'),
 ('arithmetic_audit/independent_check.py','arithmetic_audit/INDEPENDENT.json'),
]

def main():
    records=[]
    for script, receipt in CASES:
        path=ROOT/script
        expected=(ROOT/receipt).read_bytes()
        for optimized in (False,True):
            cmd=[sys.executable,'-I','-B']+(['-O'] if optimized else [])+[path.name]
            out=subprocess.run(cmd,cwd=path.parent,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=600)
            if out.returncode:
                raise RuntimeError(f'{script}, optimized={optimized}: exit {out.returncode}: {out.stderr.decode()}')
            if out.stdout!=expected:
                raise RuntimeError(f'{script}, optimized={optimized}: receipt mismatch')
            records.append({'script':script,'optimized':optimized,'receipt':receipt,'sha256':hashlib.sha256(out.stdout).hexdigest(),'status':'PASS'})
    print(json.dumps({'status':'PASS','runs':len(records),'checks':records,'scope':'Exact replay of stated finite universes; all-order claims rest on the accompanying audited proofs.'},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
