#!/usr/bin/env python3
"""Replay every mathematical checker with active optimization-safe checks."""
from pathlib import Path
import subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parent
CHECKS={'structural/check_geometric_slots.py':'structural/checks.json','poles/verify_pole_transport.py':'poles/verification.json','readout/check_pole_readout.py':'readout/receipt.json','transfer/check_carry_clock_poles.py':'transfer/checks.json','transfer/check_intercept_moments.py':'transfer/intercept_checks.json','audit/check_independent.py':'audit/independent_checks.json','audit/check_intercept_independent.py':'audit/intercept_independent.json'}
rows=[]
for script,receipt in CHECKS.items():
    expected=(ROOT/receipt).read_bytes().rstrip();outs=[]
    for flags in ([],['-O']):
        p=subprocess.run([sys.executable,*flags,str(ROOT/script)],cwd=(ROOT/script).parent,capture_output=True,check=True)
        outs.append(p.stdout)
    if outs[0]!=outs[1]:raise RuntimeError('normal/-O mismatch '+script)
    if outs[0].rstrip()!=expected:raise RuntimeError('stored receipt mismatch '+script)
    rows.append({'script':script,'normal_optimized_identical':True,'stdout_sha256':hashlib.sha256(outs[0]).hexdigest()})
print(json.dumps({'status':'PASS','count':len(rows),'scripts':rows},indent=2))
