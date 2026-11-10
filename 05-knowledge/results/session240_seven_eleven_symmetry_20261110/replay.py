#!/usr/bin/env python3
"""Replay Session240's independent finite certificates without network access."""
from pathlib import Path
import subprocess,sys,json,hashlib,tempfile
base=Path(__file__).resolve().parent
opt='--optimized' in sys.argv
jobs=[('gaps/probe.py','gaps/probe.json'),('paley/verify.py','paley/output.jsonl'),('paley/audit.py','paley/audit.out'),('provenance/verify_small.py','provenance/receipt.json'),('mixed/verify.py','mixed/output.json'),('audit_gap/audit.py','audit_gap/output.json'),('audit_design/verify_independent.py','audit_design/output.jsonl')]
def normalized(s):
    # Timing is diagnostic, never part of an exact certificate.
    try:return json.dumps(json.loads(s),sort_keys=True)
    except json.JSONDecodeError:pass
    rows=[]
    for line in s.splitlines():
        try:
            d=json.loads(line)
            if isinstance(d,dict):d.pop('elapsed_seconds',None)
            rows.append(json.dumps(d,sort_keys=True))
        except json.JSONDecodeError:rows.append(line)
    return '\n'.join(rows)
results=[]
for script,saved in jobs:
    cmd=[sys.executable]+(['-O'] if opt else [])+[str(base/script)]
    p=subprocess.run(cmd,cwd=base,text=True,capture_output=True)
    if p.returncode:raise RuntimeError(script+' failed: '+p.stderr)
    expected=(base/saved).read_text()
    if normalized(p.stdout)!=normalized(expected):raise RuntimeError(script+' output changed')
    results.append({'script':script,'status':'PASS','saved_output':saved})
with tempfile.TemporaryDirectory(prefix='s240-audit-') as tmp:
    binary=str(Path(tmp)/'audit19')
    subprocess.run(['g++','-std=c++17','-O2',str(base/'audit_gap/audit19.cpp'),'-o',binary],check=True,capture_output=True)
    p=subprocess.run([binary],check=True,text=True,capture_output=True)
    expected=(base/'audit_gap/output19.txt').read_text()
    if p.stdout!=expected:raise RuntimeError('19-vertex DP output changed')
    results.append({'script':'audit_gap/audit19.cpp','status':'PASS','saved_output':'audit_gap/output19.txt'})
print(json.dumps({'status':'PASS','python_optimized':opt,'checks':results},indent=2,sort_keys=True))
