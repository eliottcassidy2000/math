#!/usr/bin/env python3
"""Cross-check actual producer against independent enumerated geometry."""
from pathlib import Path
import sys,json
from fractions import Fraction as F
from random import Random
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lrc'))
from witness import choose
# Import the independent geometric oracle; its self-check corpus also runs.
from check_independent import cell_measure,events,bad,require
rng=Random(123321)
accepted=rejected=0
for _ in range(500):
    k=rng.randrange(1,7);vs=[rng.randrange(1,13) for _ in range(k)]
    ps=[F(rng.randrange(-8,9),rng.randrange(1,10)) for _ in vs]
    a=F(rng.randrange(-8,9),rng.randrange(1,9));b=a+F(rng.randrange(1,12),rng.randrange(1,9))
    r=F(1,rng.randrange(3,17))
    ms=[cell_measure(a,b,v,p,r) for v,p in zip(vs,ps)]
    domain=all(m<b-a for m in ms) and sum(ms)<=b-a
    try:t,stats=choose((a,b),vs,r,ps)
    except ValueError:
        require(not domain,'producer rejected admissible input');rejected+=1
    else:
        require(domain,'producer accepted outside promised domain')
        require(a<=t<=b and not any(bad(t,v,p,r) for v,p in zip(vs,ps)),'producer unsafe witness')
        require(stats['bisections']<=stats['depth_bound'],'producer exceeded bound');accepted+=1
print(json.dumps({'producer_crosscheck':'PASS','accepted':accepted,'rejected':rejected},sort_keys=True))
