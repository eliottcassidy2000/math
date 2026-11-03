#!/usr/bin/env python3
"""Compare independent rational certificates with the inherited interval replay."""
from fractions import Fraction
from pathlib import Path
import json,re
p=Path(__file__).resolve().parent
match=re.search(r'^C = \[([^,]+), ([^\]]+)\]',(p/'inherited/certify_mrb.replay.out').read_text(),re.M)
if match is None:raise RuntimeError('Missing inherited C interval')
old_lo,old_hi=map(Fraction,match.groups())
rows=[]
for mode,name in [('normal','exact_precision_results.json'),('optimized','exact_precision_results_optimized.json')]:
    data=json.loads((p/'audit'/name).read_text())
    for row in data['certificates']:
        lo,hi=Fraction(row['C_lower']),Fraction(row['C_upper'])
        if not lo<old_lo<old_hi<hi:
            raise RuntimeError('Independent C enclosures disagree at p='+str(row['p']))
        rows.append({'mode':mode,'p':row['p'],'inherited_interval_strictly_contained':True})
result={'status':'PASS','comparison':'Exact rational comparisons against parsed decimal endpoints','rows':rows,
        'scope':'Consistency check of two elementary evaluation routes sharing the proved Euler remainder mechanism; not independence of the analytic theorem.'}
(p/'certificate_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS: every executed rational-precision interval strictly contains the independently replayed inherited interval.')
