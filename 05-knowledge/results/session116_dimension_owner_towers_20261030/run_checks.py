#!/usr/bin/env python3
"""Replay the exact finite checks normally and optimized; compare output bytes."""
from pathlib import Path
import subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parent
SCRIPTS=[
 'lrc/check_second_band.py','audit_lrc/check_independent.py',
 'tower/check_tower.py','tower/check_spectrum.py','tower/check_integer_shells.py',
 'audit_tower/check_independent.py','audit_tower/check_integer_controls.py',
 'arithmetic/check_arithmetic.py','audit_arithmetic/check_independent.py',
 'collatz/check_phase_completion.py','audit_collatz/check_independent.py',
 'check_synthesis.py','audit_sparse_merge/check.py']
rows=[]
for name in SCRIPTS:
 p=ROOT/name;outputs=[]
 for flags in ([],['-O']):
  r=subprocess.run([sys.executable,*flags,p.name],cwd=p.parent,capture_output=True,check=True)
  outputs.append(r.stdout)
 if outputs[0]!=outputs[1]:raise RuntimeError('Normal/optimized mismatch: '+name)
 rows.append({'script':name,'normal_optimized':'identical','stdout_sha256':hashlib.sha256(outputs[0]).hexdigest()})
 print('PASS '+name,file=sys.stderr,flush=True)
print(json.dumps({'status':'PASS','scripts':len(rows),'runs':2*len(rows),'checks':rows},indent=2))
