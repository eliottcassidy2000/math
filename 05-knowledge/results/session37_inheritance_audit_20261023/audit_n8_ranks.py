#!/usr/bin/env python3
"""Contrast common-slope and individual-class affine failure with true fiber recovery."""
import csv,json
from collections import defaultdict
from pathlib import Path
from audit_polynomial_identities import rank
base=Path(__file__).resolve().parent
rows=[{k:int(v) for k,v in r.items()} for r in csv.DictReader((base/'n8_fibers.tsv').open(),delimiter='\t')]
groups=defaultdict(list)
for r in rows:
    char=(-r['c3'],-r['c4'],-r['c5'],r['D33']-r['c6'],r['D34']-r['c7'],r['E8']-r['c8'])
    groups[char].append(r)
pool=[];bad=[];ranktypes=defaultdict(int)
for char,rs in groups.items():
    b=rs[0];d=[]
    for r in rs[1:]:
        # On a fixed char fiber, unsigned e_m deltas = 2*cycle-count deltas.
        d.append([2*(r['c6']-b['c6']),2*(r['c7']-b['c7']),2*(r['c8']-b['c8']),r['H']-b['H']])
    rx,ry=rank([r[:-1] for r in d]),rank(d)
    ranktypes[(rx,ry)]+=1;pool+=d
    if ry>rx:bad.append({'char_nontrivial_signed_coefficients':char,'rows':rs,'ranks':[rx,ry]})
assert len(rows)==3217
out={'n8_joint_polynomial_fibers':len(rows),'characteristic_polynomial_fibers':len(groups),
     'pooled_rank_e_unsigned':rank([r[:-1] for r in pool]),'pooled_rank_e_unsigned_H':rank(pool),
     'classes_with_individual_affine_failure':len(bad),
     'class_rank_histogram':{str(k):v for k,v in sorted(ranktypes.items())},
     'smallest_individual_affine_obstruction':min(bad,key=lambda x:len(x['rows'])) if bad else None}
(base/'audit_n8_ranks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
