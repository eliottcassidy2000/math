#!/usr/bin/env python3
"""Replay every published checker normally and optimized, comparing exact output."""
from pathlib import Path
import subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parent
SCRIPTS=[
 'quadratic/check_odd_eighths.py','audit_phases/check_boolean_independent.py',
 'unitary/check_word_unitary.py','audit_phases/check_unitary_independent.py',
 'switching/check_eighth_groups.py','audit_groups/check_group_weights.py',
 'lrc/check_harmonic.py','audit_lrc/independent_harmonic.py',
 'collatz/check_one_reset.py','collatz/check_constant.py','collatz/check_general_reset.py',
 'audit_collatz/check_classes_independent.py','audit_collatz/check_general_independent.py',
 'check_synthesis.py','audit_synthesis/check_cross_maps.py',
 'audit_incoming/check_independent.py']
rows=[]
for name in SCRIPTS:
 p=ROOT/name;outputs=[]
 for flags in ([],['-O']):
  run=subprocess.run([sys.executable,*flags,p.name],cwd=p.parent,capture_output=True,check=True)
  outputs.append(run.stdout)
 if outputs[0]!=outputs[1]:raise RuntimeError('Output differs under optimization: '+name)
 rows.append({'script':name,'normal_optimized':'identical','stdout_sha256':hashlib.sha256(outputs[0]).hexdigest()})
 print('PASS '+name,file=sys.stderr,flush=True)
print(json.dumps({'status':'PASS','scripts':len(rows),'runs':2*len(rows),'checks':rows},indent=2))
