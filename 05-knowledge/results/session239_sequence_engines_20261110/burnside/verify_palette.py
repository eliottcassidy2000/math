#!/usr/bin/env python3
"""Independent partition and exhaustive orbit verification of palette DP."""
from itertools import combinations,permutations,product
from pathlib import Path
from hashlib import sha256
import json,time
from palette_engines import palette_series,direct_palette
from engines import integer_response,require

def brute(n,q,t):
    """Enumerate EVERY labelled palette assignment and quotient by S_n.
    Reversal fixes 0..t-1 and swaps adjacent states t,t+1; t+2,t+3; ... .
    """
    edges=list(combinations(range(n),2));index={e:i for i,e in enumerate(edges)}
    transforms=[]
    for p in permutations(range(n)):
        transforms.append([(index[tuple(sorted((p[a],p[b])))],p[a]>p[b]) for a,b in edges])
    def reverse(x):return x if x<t else t+((x-t)^1)
    representatives=set()
    for colors in product(range(q),repeat=len(edges)):
        orbit=[]
        for transform in transforms:
            image=[0]*len(edges)
            for i,(j,flip) in enumerate(transform):image[j]=reverse(colors[i]) if flip else colors[i]
            orbit.append(tuple(image))
        representatives.add(min(orbit))
    return len(representatives)

cases={'tournaments':(2,0),'graphs':(2,2),'digraphs':(4,2),'oriented_graphs':(3,1),'empty_palette':(0,0),'one_color':(1,1),'two_reversal_pairs':(4,0)}
report={}
for name,(q,t) in cases.items():
    start=time.perf_counter();a,stats=palette_series(12,q,t);duration=time.perf_counter()-start
    b=direct_palette(12,q,t);require(a==b,'direct partition equality '+name)
    brute_values=[brute(n,q,t) for n in range(5)];require(a[:5]==brute_values,'exhaustive orbit equality '+name)
    report[name]={'q':q,'t':t,'prefix_0_to_12':a,'brute_0_to_4':brute_values,'dp_seconds':duration,'stats':stats}
require(palette_series(60,2,0)[0]==integer_response(60)[0],'tournament specialization through60')
# Boundary rejection is deliberate: no reversal involution on these palettes.
for q,t in [(2,1),(1,2),(-1,0),(3,-1)]:
 try:palette_series(3,q,t)
 except AssertionError:pass
 else:raise AssertionError('invalid palette accepted')
# Numerical controls are historical well-known values, not newly web-verified.
controls={'graphs':[1,1,2,4,11,34,156,1044,12346], 'digraphs':[1,1,3,16,218,9608], 'oriented_graphs':[1,1,2,7,42,582]}
for name,values in controls.items():require(report[name]['prefix_0_to_12'][:len(values)]==values,'historical sequence control '+name)
print(json.dumps(report,indent=2));print('PASS: 7 palettes through12, exhaustive pair-color orbits through4, tournament specialization through60, invalid-palette controls')
Path(__file__).with_name('palette.json').write_text(json.dumps(report,indent=2)+'\n')
