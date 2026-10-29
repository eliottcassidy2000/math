#!/usr/bin/env python3
"""Replay the explicit mathematical test universe, without mutating receipts."""
import hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CASES=[
('audit/check_syndrome_independent.py','audit/syndrome-independent.json'),
('audit/check_native_quotient_independent.py','audit/native-quotient-independent.json'),
('grammar/check.py','grammar/results.json'),
('switching/check_signed_cuts.py','switching/results.json'),
('pfaffian/check_pfaffian.py','pfaffian/results.json'),
('metagraph/check_sign_grammar.py','metagraph/checks.json'),
('metagraph/check_extremal_syndrome.py','metagraph/syndrome_checks.json'),
('metagraph/check_cayley_quotient.py','metagraph/cayley_checks.json'),
('audit/check_signed_cuts_independent.py','audit/signed-cuts.json'),
('audit/check_switch_groups_independent.py','audit/switch-groups.json'),
('audit/check_historical_hostiles_independent.py','audit/historical-hostiles.json'),
('audit/check_two_move_boundary_independent.py','audit/two-move.json'),
('audit/check_pfaffian_independent.py','audit/pfaffian-independent.json'),
('audit/check_projective_group_independent.py','audit/projective-independent.json')]
rows=[]
for script,receipt in CASES:
 saved=(ROOT/receipt).read_bytes()
 a=subprocess.check_output([sys.executable,'-B',str(ROOT/script)],cwd=ROOT)
 b=subprocess.check_output([sys.executable,'-B','-O',str(ROOT/script)],cwd=ROOT)
 if a!=b:raise RuntimeError('optimized mismatch '+script)
 if (ROOT/receipt).read_bytes()!=saved:raise RuntimeError('receipt mutated '+script)
 if json.loads(a)!=json.loads(saved):raise RuntimeError('stored receipt mismatch '+script)
 rows.append({'script':script,'receipt':receipt,'stdout_sha256':hashlib.sha256(a).hexdigest(),'normal_optimized_identical':True,'stored_receipt_matches':True})
print(json.dumps({'status':'VERIFIED','checks':rows,'scope':'Fourteen mathematical scripts, with explicit finite universes and separate all-order proofs. No full prime atlas or CI-pass claim.'},indent=2,sort_keys=True))
