#!/usr/bin/env python3
"""Extend the incoming 78 core bank through the joint two-tail dilation lemma."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
from joint_phase import pair_bank,need

records=[]
for core in combinations(range(1,14),11):
    bank=pair_bank(tuple(2*x for x in core),F(1,14))
    need(bank['all_pairs_accepted'],'doubled core rejected a pair')
    records.append(bank)
result=dict(incoming_revision='373661611e',core_universe='all eleven-subsets of[1,13]',
    theorem_dependency='DILATIONS.md D1-D3 plus incoming base-bank L4',
    families=records,core_count=len(records),
    finite_queries=sum(b['finite_queries'] for b in records),
    maximum_both_large_cutoff=max(b['both_large_cutoff'] for b in records),
    all_dilations_all_positive_pairs=True,arbitrary_primitive_cores=False)
Path('dilation78.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='families'},indent=2))
