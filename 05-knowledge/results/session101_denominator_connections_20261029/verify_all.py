#!/usr/bin/env python3
"""Run only the explicit mathematical packet checks, normally and optimized."""
import hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CASES=[('tournaments/verify_tournament_controls.py','tournaments/verification.json'),('spectra/check_spectra.py','spectra/checks.json'),('continuations/check_continuations.py','continuations/checks.json'),('symmetry/check_swap_observer.py','symmetry/receipt.json'),('squarefree/check_squarefree.py','squarefree/checks.json'),('lattices/check_lattices.py','lattices/checks.json'),('audit/independent_check.py','audit/results.json'),('audit/squarefree_check.py','audit/squarefree_results.json'),('audit/lattice_check.py','audit/lattice_results.json')]
rows=[]
for script,receipt in CASES:
 a=subprocess.check_output([sys.executable,'-B',str(ROOT/script)],cwd=ROOT)
 b=subprocess.check_output([sys.executable,'-B','-O',str(ROOT/script)],cwd=ROOT)
 if a!=b:raise RuntimeError('optimized mismatch '+script)
 if json.loads(a)!=json.loads((ROOT/receipt).read_text()):raise RuntimeError('stored receipt mismatch '+script)
 rows.append({'script':script,'receipt':receipt,'stdout_sha256':hashlib.sha256(a).hexdigest(),'normal_optimized_identical':True,'stored_receipt_matches':True})
print(json.dumps({'status':'VERIFIED','checks':rows,'scope':'Nine mathematical scripts. No unrestricted Collatz conclusion or CI pass.'},indent=2,sort_keys=True))
