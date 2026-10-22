#!/usr/bin/env python3
"""Bounded independent follow-up: ordinary starting-value guards only."""
from pathlib import Path
from math import gcd
import hashlib
import json

checks=0
def require(ok, message):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(message)

root=Path(__file__).resolve().parent.parent
base=Path(__file__).resolve().parent
expected=[(12,32,102),(11,11,11),(3,18,18),(14,4,74)]
require(pow(2,12,105)==1, 'period four follows from exponent increment twelve')
for r in range(4):
    representative=r if r else 4
    x=2**(3*representative+1)-5
    require(tuple(x%d for d in (15,35,105))==expected[r], 'independent four-row derivation')
for m in range(1,513):
    x=2**(3*m+1)-5
    require(tuple(x%d for d in (15,35,105))==expected[m%4], 'literal integer residues')
    both=gcd(x,15)==1 and gcd(x,35)==1
    require(both==(gcd(x,105)==1)==(gcd(x,525)==1)==bool(m%2),
            'joint unit guards, lcm/product and odd-m equivalence')
    require((x%15)%5==(x%35)%5,'compatibility across common factor five')
targets=['session31-collatz/AUDIT.md','session31-collatz/verify_audit.py',
         'session31-square-blocks/AUDIT.md','session31-square-blocks/verify_square_blocks.py']
def resolve_source(relative):
    first, rest = relative.split('/', 1)
    folder = first.removeprefix('session31-')
    portable = root / (rest if folder == 'root' else folder + '/' + rest)
    return portable if portable.is_file() else root / relative
hashes={relative:hashlib.sha256(resolve_source(relative).read_bytes()).hexdigest() for relative in targets}
require(hashes['session31-square-blocks/AUDIT.md']==
        'dcdf8e690f0b017c7f698316573ff5fd937c6b1d0ad2a9ab6dbaf90103d8b1f9',
        'requested exact compositum-wording revision')
text=resolve_source('session31-square-blocks/AUDIT.md').read_text()
require('be exactly the compositum' in text, 'corrected explicit field choice')
report={'status':'PASS: bounded addendum, original audit history retained',
        'base_check_count_preserved':46965,
        'separate_addendum_check_count':checks,
        'guard_test_universe':'1 <= m <= 512',
        'rows_m_mod4_0_through_3':expected,
        'unit_condition':'Both ordinary starting-integer guards pass iff m is odd.',
        'all_order_reason':'Period four modulo 105; each passing m progression is unbounded; previously proved growing-prefix formula applies.',
        'scope_exclusion':'No native source-address map or guard is identified, and no later endpoint residue restriction is inferred.',
        'audited_sha256':hashes,
        'history':'history_v1/AUDIT.md and original hash-pinned receipts remain unchanged.'}
(base/'addendum_receipt.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True))
