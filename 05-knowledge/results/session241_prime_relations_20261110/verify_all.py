#!/usr/bin/env python3
"""Replay every session241 producer and independent audit in normal/-O mode."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
JOBS=[
 ('constructions/verify.py','constructions/exact.json'),
 ('coverage/verify.py','coverage/verification.json'),
 ('coverage_audit/verify.py','coverage_audit/verification.json'),
 ('repunit_chart/verify.py','repunit_chart/exact.json'),
 ('embedding_audit/verify.py','embedding_audit/verification.json'),
 ('provenance/verify_smooth_lengths.py',None),
]
report=[]
for script,output in JOBS:
    normal=subprocess.run([sys.executable,'-B',str(ROOT/script)],cwd=ROOT,capture_output=True,text=True,check=True)
    normal_bytes=(ROOT/output).read_bytes() if output else normal.stdout.encode()
    optimized=subprocess.run([sys.executable,'-B','-O',str(ROOT/script)],cwd=ROOT,capture_output=True,text=True,check=True)
    optimized_bytes=(ROOT/output).read_bytes() if output else optimized.stdout.encode()
    if normal_bytes!=optimized_bytes or normal.stdout!=optimized.stdout:
        raise RuntimeError('normal/optimized mismatch: '+script)
    report.append({'script':script,'status':'PASS','modes':['normal','optimized'],'sha256':hashlib.sha256(normal_bytes).hexdigest()})
result={'status':'PASS','processes':2*len(JOBS),'checks':report,'scope':'Six independently scoped producer/audit programs; not unrestricted tournament coverage or a global repository test.'}
text=json.dumps(result,indent=2,sort_keys=True)+'\n'
(ROOT/'verification.json').write_text(text)
print(text,end='')
