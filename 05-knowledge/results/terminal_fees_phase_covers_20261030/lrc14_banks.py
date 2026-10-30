#!/usr/bin/env python3
"""Fixed twelve-speed cores, every positive thirteenth speed, radius 1/14."""
from phase_cover import safe_intervals,finite_rejections,need,F,period,attachment,dist
from itertools import combinations
from random import Random
from pathlib import Path
import json

cores=list(combinations(range(1,14),12))
cores += [tuple(26*j for j in range(1,13)),
          (2,3,5,7,11,13,17,19,23,29,31,37),tuple(range(1,24,2)),
          (1,2,3,5,8,13,21,34,55,89,97,101),
          tuple(15*j for j in (1,2,4,5,7,9,11,12,14,16,19,23))]
rng=Random(20261030)
cores += [tuple(sorted(rng.sample(range(1,81),12))) for _ in range(10)]
rows=[]
for core in cores:
    intervals=safe_intervals(core,F(1,14));cutoff,rejected=finite_rejections(intervals,F(1,14))
    need(not rejected, ('actual uncovered input',core,rejected))
    for v in range(1,cutoff):
        need(bool(safe_intervals((*core,v),F(1,14))), 'independent literal new comb')
        t=attachment(intervals,F(1,14),v)
        need(t is not None and min(dist(s*t) for s in (*core,v))>=F(1,14), 'native target witness')
    rows.append(dict(core=core,components=len(intervals),phase_period=period(intervals),
                     all_speeds_from=cutoff,finite_checked=cutoff-1,rejected=[]))
out=dict(status='VERIFIED fixed-core all-positive-speed coverage with proved finite tail',
         target_radius='1/14',core_count=len(rows),finite_queries=sum(r['finite_checked'] for r in rows),
         scope='These 28 twelve-speed cores only; arbitrary-input LRC14 remains OPEN',cores=rows)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='cores'},indent=2))
