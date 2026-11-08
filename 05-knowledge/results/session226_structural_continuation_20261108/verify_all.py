#!/usr/bin/env python3
"""Replay every portable mathematical checker in both Python modes."""
import hashlib,json,subprocess,sys
from pathlib import Path
BASE=Path(__file__).resolve().parent
JOBS=[
 ('strong_grafts/verify.py','strong_grafts/CHECK.json'),
 ('root_descent/verify.py','root_descent/CHECK.json'),
 ('graft/verify.py','graft/CHECK.json'),
 ('graft/inherited_verify.py','graft_audit/inherited.json'),
 ('modular/check.py','modular/CHECK-OUTPUT.txt'),
 ('structural_audit/independent.py','structural_audit/independent.json'),
 ('graft_audit/independent_check.py','graft_audit/independent.json'),
 ('modular_audit/independent.py','modular_audit/independent-normal.txt'),
]
rows=[]
for source,receipt in JOBS:
    want=(BASE/receipt).read_bytes()
    for optimized in (False,True):
        command=[sys.executable,'-I','-B']+(['-O'] if optimized else [])+[str(BASE/source)]
        result=subprocess.run(command,cwd=(BASE/source).parent,capture_output=True,check=False)
        if result.returncode or result.stdout!=want:
            raise RuntimeError((source,optimized,result.returncode,result.stderr.decode(),result.stdout.decode()))
        rows.append(dict(source=source,optimized=optimized,stdout_sha256=hashlib.sha256(want).hexdigest()))
print(json.dumps(dict(status='PASS',runs=len(rows),checks=rows),indent=2,sort_keys=True))
