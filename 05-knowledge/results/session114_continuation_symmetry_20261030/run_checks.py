#!/usr/bin/env python3
"""Run the published checks normally and optimized; compare deterministic outputs."""
from pathlib import Path
import subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parent
SCRIPTS=[
 'aliquot/check_aliquot.py','audit/check_aliquot_quotient.py',
 'jacobian/check_one_coordinate.py','charts/check_parent_multiplicity.py',
 'audit/check_symmetry_independent.py',
 'lrc/check_resonant_tail.py','lrc/check_general.py',
 'lrc_audit/independent_core.py','lrc_audit/independent_general.py',
 'collatz/check_periodic_rates.py','collatz_audit/check_periodic_independent.py',
 'collatz_audit/check_literal_means.py']
rows=[]
for rel in SCRIPTS:
 p=ROOT/rel;results=[]
 for flags in ([],['-O']):
  r=subprocess.run([sys.executable,*flags,p.name],cwd=p.parent,capture_output=True,check=True)
  results.append(r.stdout)
 if results[0]!=results[1]:raise RuntimeError('Normal/optimized mismatch: '+rel)
 rows.append({'script':rel,'normal_optimized':'identical','stdout_sha256':hashlib.sha256(results[0]).hexdigest()})
 print('PASS '+rel,file=sys.stderr,flush=True)
print(json.dumps({'status':'PASS','scripts':len(rows),'runs':2*len(rows),'checks':rows},indent=2))
