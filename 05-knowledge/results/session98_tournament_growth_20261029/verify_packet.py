#!/usr/bin/env python3
"""Replay all mathematical checks, with optimization-safe stdout equality."""
from pathlib import Path
import subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parent
SCRIPTS=['recurrence/verify.py','composition/inherited_recreate.py','composition/verify_surplus_jets.py','composition/observable_quotient.py','composition/search_reachable.py','partials/check_partial_exposure.py','defects/check_defects.py','defects/check_recursive_grammar.py','limits/verify_limits.py','limits/verify_all_poles.py','limits/verify_palette_poles.py','audit/independent_audit.py','audit/check_observation_independent.py']
RECEIPTS={
'recurrence/verify.py':'recurrence/receipt.json',
'composition/verify_surplus_jets.py':'composition/receipt.json',
'composition/observable_quotient.py':'composition/quotient_receipt.json',
'composition/search_reachable.py':'composition/reachable_receipt.json',
'partials/check_partial_exposure.py':'partials/partial_exposure.json',
'defects/check_defects.py':'defects/results.json',
'defects/check_recursive_grammar.py':'defects/recursion.json',
'limits/verify_limits.py':'limits/verification.json',
'limits/verify_all_poles.py':'limits/all_poles_verification.json',
'limits/verify_palette_poles.py':'limits/palette_poles_verification.json',
'audit/independent_audit.py':'audit/normal.json',
'audit/check_observation_independent.py':'audit/observation.normal.json'}
rows=[]
for rel in SCRIPTS:
    path=ROOT/rel
    expected=(ROOT/RECEIPTS[rel]).read_bytes().rstrip() if rel in RECEIPTS else None
    outs=[]
    for flags in ([],['-O']):
        p=subprocess.run([sys.executable,*flags,str(path)],cwd=path.parent,capture_output=True,check=True)
        outs.append(p.stdout)
    if outs[0]!=outs[1]:raise RuntimeError('normal/-O mismatch: '+rel)
    if not outs[0].strip():raise RuntimeError('empty checker output: '+rel)
    if rel in RECEIPTS and outs[0].rstrip()!=expected:
        raise RuntimeError('recorded receipt mismatch: '+rel)
    # Equality covers all stdout bytes, not just a final PASS marker.
    rows.append({'script':rel,'normal_optimized_identical':True,'stdout_sha256':hashlib.sha256(outs[0]).hexdigest()})
print(json.dumps({'status':'PASS','scripts':rows,'count':len(rows)},indent=2))
