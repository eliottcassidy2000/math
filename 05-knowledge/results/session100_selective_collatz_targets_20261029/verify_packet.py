#!/usr/bin/env python3
"""Replay all exact checkers with optimization-safe receipt comparison."""
from pathlib import Path
import sys,subprocess,json,hashlib
ROOT=Path(__file__).resolve().parent
CHECKS={'bijections/check_valuation_transfer.py':'bijections/checks.json','controllers/check_return_controller_trap.py':'controllers/checks.json','drift/check_drift.py':'drift/checks.json','atlas/check_affine_atlas.py':'atlas/checks.json','audit/check_independent.py':'audit/checks.json'}
rows=[]
for script,receipt in CHECKS.items():
    expected=(ROOT/receipt).read_bytes().rstrip();outputs=[]
    for flags in ([],['-O']):
        p=subprocess.run([sys.executable,*flags,str(ROOT/script)],cwd=(ROOT/script).parent,capture_output=True,check=True)
        outputs.append(p.stdout)
    if outputs[0]!=outputs[1]:raise RuntimeError('normal/-O mismatch '+script)
    if json.loads(outputs[0])!=json.loads(expected):raise RuntimeError('JSON receipt mismatch '+script)
    rows.append({'script':script,'normal_optimized_identical':True,'stdout_sha256':hashlib.sha256(outputs[0]).hexdigest()})
print(json.dumps({'status':'PASS','count':len(rows),'scripts':rows},indent=2))
