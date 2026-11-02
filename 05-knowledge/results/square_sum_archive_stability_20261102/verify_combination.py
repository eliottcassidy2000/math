#!/usr/bin/env python3
"""Exact check that the incoming five-cycle payment uses unspent terms."""
import json
from fractions import Fraction as F
from verify_reduction import symbolic,need

rows=[]
for row in symbolic():
    remainder={k:F(v) for k,v in row['nonnegative_remainder'].items()}
    for term,payment in [('2,5',1),('3,5',3)]:
        need(remainder.get(term,F(0))>=payment,'cycle payment already spent')
        remainder[term]-=payment
    need(all(v>=0 for v in remainder.values()),'negative combined ledger')
    rows.append({'full_cover_multiplier':row['full_cover_multiplier'],
                 'remaining_coefficients':{k:str(v) for k,v in remainder.items() if v}})
print(json.dumps({'status':'passed','scope':'Exact all-matrix cycle allocation; analytic five-vertex identity is cited in INCOMING.md.',
                  'combined_C5_coefficient':'1','retained_C6_coefficient':'1/2',
                  'signed_ledgers':rows},indent=2,sort_keys=True))
