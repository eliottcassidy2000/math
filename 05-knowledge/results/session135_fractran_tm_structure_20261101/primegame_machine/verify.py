#!/usr/bin/env python3
"""Finite-exact checks of the basin theorem; all gates survive python -O."""
from pathlib import Path
from itertools import product
from collections import Counter
from math import prod
import hashlib,json,random
import primegame_basin as pg

HERE=Path(__file__).resolve().parent
checks=0

def check(c,m):
    global checks
    checks+=1
    if not c: raise RuntimeError(m)

# Independently printed list and literal division; no vector rule import.
RACK=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),
      (77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
def raw(n):
    for i,(a,b) in enumerate(RACK):
        if n%b==0: return n//b*a,'ABCDEFGHIJKLMN'[i]
    raise RuntimeError('total rack')

# A witness at C(n,n-1) is literally on Conway's standard orbit; no
# untested prime-detection shortcut or candidate integer threshold is used.
goals={5**n*7**(n-1)*13 for n in range(2,100)}
cache={n:'prime' for n in goals}
max_raw_steps=0

def raw_basin(n):
    global max_raw_steps
    path=[]; positions={}; x=n
    while x not in cache and x not in positions:
        positions[x]=len(path); path.append(x)
        x,_=raw(x)
        # Finite test safety, not the theorem's total decision procedure.
        check(len(path)<=2000000,'raw control trace exceeded declared test guard')
    if x in cache:
        label=cache[x]
    else:
        period=len(path)-positions[x]
        check(period in (2,6),'unexpected literal cycle')
        label='cycle'+str(period)
    for y in path: cache[y]=label
    max_raw_steps=max(max_raw_steps,len(path))
    return label

pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
for p in pins['inputs']:
    check(hashlib.sha256((HERE/p['local']).read_bytes()).hexdigest()==p['sha256'],'input source pin')

# Exact all-rule algebra. Columns are dT,dQ,dk.
expected=((0,0,0),(0,0,0),(0,-1,0),(0,0,0),(0,0,0),(0,0,0),
          (0,0,0),(0,0,0),(0,-1,-1),(0,0,0),(0,0,0),(0,1,0),
          (0,-1,0),(1,0,1))
deltas=[]
for i,(u,v) in enumerate(RACK):
    a,_=pg.valuations(u); b,_=pg.valuations(v)
    delta=tuple(x-y for x,y in zip(pg.measures(a),pg.measures(b)))
    check(delta==expected[i],'measure increment '+str(i))
    deltas.append(dict(rule=pg.LETTERS[i],numerator=u,denominator=v,dT=delta[0],dQ=delta[1],dk=delta[2]))

onehot=Counter(); onehot_n=0
for q in pg.CONTROLS:
    for a,b,c,d in product(range(6),repeat=4):
        x=[a,b,c,d,0,0,0,0,0,0]
        if q: x[pg.PRIMES.index(q)]=1
        n=pg.integer(x)
        predicted='prime' if pg.onehot_escape(x) else 'cycle2'
        actual=raw_basin(n)
        check(predicted==actual,('onehot',x,predicted,actual))
        onehot[predicted]+=1; onehot_n+=1

box=Counter(); max_classify_steps=0; box_count=0
for x in product(range(3),repeat=10):
    n=pg.integer(x)
    y,letter=pg.step(x); z,raw_letter=raw(n)
    check(pg.integer(y)==z and letter==raw_letter,('step',x))
    T,Q,k=pg.measures(x); t,q,h=pg.measures(y)
    if k>0:
        check(letter in 'ABCDEFGHIJK' and T==t and Q>=q and k>=h,('budget',x))
    result=pg.classify_vector(x)
    check(result['basin']==raw_basin(n),('box basin',x,result))
    check(result['transitions']<=result['bound'],'decision bound')
    max_classify_steps=max(max_classify_steps,result['transitions'])
    box[result['basin']]+=1; box_count+=1

ints=Counter()
for n in range(1,30001):
    result=pg.classify_integer(n)
    # A spectator-scaled trajectory is compared after only fixed-prime extraction.
    x,s=pg.valuations(n)
    check(pg.integer(x,s)==n,'valuation reconstruction')
    check(result['basin']==raw_basin(n//s),('integer basin',n,result))
    ints[result['basin']]+=1

rng=random.Random(135)
random_count=300
for _ in range(random_count):
    x=tuple(rng.randrange(9) for _ in range(10))
    result=pg.classify_vector(x)
    check(result['basin']==raw_basin(pg.integer(x)),('higher box sample',x,result))
    max_classify_steps=max(max_classify_steps,result['transitions'])

# Guard loss hostiles: k cannot be discarded; positive-control premise matters.
hostiles={}
for n in (1,2,3,5,11,13,17,23,29,31*2,31*55,4301,4301*31,121,143,91,247):
    hostiles[str(n)]=pg.classify_integer(n)
check(hostiles['4301']['basin']=='cycle6','three-control cycle hostile')
check(hostiles['2']['basin']=='prime' and hostiles['11']['basin']=='cycle2','total != divergent')
check(hostiles['23']['basin']=='prime','control23 always escape')
# The same data counters (0,0,0,0) can have different one-token outcomes.
check(hostiles['17']['basin']=='cycle2' and hostiles['23']['basin']=='prime','phase loss hostile')
check(pg.measures(pg.step(pg.valuations(2)[0])[0])[1]>pg.measures(pg.valuations(2)[0])[1],
      'Q is not a global decreasing potential')
check(pg.onehot_escape(pg.valuations(2*11)[0]) is False,'inactive JK trap')

report=dict(status='FINITE-EXACT PASS; all-input claims proved in THEORY.md',
            checks=checks,onehot=dict(count=onehot_n,census=dict(sorted(onehot.items())),counter_range=[0,5]),
            capped_box=dict(count=box_count,census=dict(sorted(box.items())),exponent_range=[0,2]),
            integers=dict(range=[1,30000],census=dict(sorted(ints.items()))),
            random=dict(count=random_count,seed=135,exponent_range=[0,8]),
            max_decider_transitions=max_classify_steps,max_uncached_raw_control_path=max_raw_steps,
            raw_control_cache_size=len(cache),measure_changes=deltas,hostiles=hostiles)
print(json.dumps(report,sort_keys=True,indent=2))
