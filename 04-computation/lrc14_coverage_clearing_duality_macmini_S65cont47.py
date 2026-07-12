#!/usr/bin/env python3
"""cont.47 NEW ANGLE: the coverage-clearing duality. Test whether 7-sector coverage
deficiency (density side, my three-gap) predicts mod-q unit-clearing (liveness side,
klein's crux). If spread=bad-7-sector-coverer => bad-mod-q-coverer => clears, then the
three-gap coverage advantage is the master key for the WHOLE endgame."""
from fractions import Fraction as F
from math import gcd
import random
def measS7(v):  # 7-sector coverage deficiency = P(some sector empty), density side
    E=list(v)
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(abs(e)):
            for s in range(8): pts.add(F(m,abs(e))+F(s,7*abs(e)))
    pts=sorted(x for x in pts if 0<=x<=1); defc=F(0)
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        if len(hit)<7: defc+=b-a
    return float(defc)
def clearing_count(v,q):  # liveness side (klein): # multipliers p in [1,q-1] with some runner lonely
    cnt=0
    for p in range(1,q):
        # lonely at p/q iff some sector empty in the DISCRETE sense: min_i ||v_i p/q|| ... 
        # klein: clears iff a unit residue is MISSED by {v_i p mod q}. Use: exists p, all ||v_i p/q||>=1/14
        ok=True
        for e in v:
            r=(e*p)%q
            d=min(r,q-r)
            if d < q/14: ok=False; break
        if ok: cnt+=1
    return cnt
# residual DC-like families: covering, primitive, spread
def make_dc(seed):
    random.seed(seed)
    # a spread divisor-complete-ish family: include small divisors + spread large
    base=[2,3,4,6,12]; extra=sorted(random.sample(range(13,60),8))
    v=sorted(set(base+extra))[:13]
    while len(v)<13: v.append(max(v)+random.randint(1,5)); v=sorted(set(v))
    g=v[0]
    for x in v[1:]: g=gcd(g,x)
    return [x//g for x in v][:13]
print("coverage-clearing duality: 7-sector deficiency (density) vs mod-q clearing count (liveness)")
print(f"{'family type':22s} {'measS7(defic)':>13s} {'clearing q<=31':>15s}  (bad coverer=>clears?)")
def show(nm,v):
    d=measS7(v)
    tot=sum(clearing_count(v,q) for q in range(15,32))
    print(f"  {nm:20s} {d:13.4f} {tot:>15d}")
    return d,tot
# AP (good coverer -- low deficiency, should clear LITTLE)
show("AP {1..13}",list(range(1,14)))
# spread DC (bad coverer -- high deficiency, should clear MUCH)
pts=[]
for s in range(8):
    v=make_dc(s); d,t=show("spread-DC seed%d"%s,v); pts.append((d,t))
# correlation
import statistics
ds=[p[0] for p in pts]; ts=[p[1] for p in pts]
apd=measS7(list(range(1,14))); apt=sum(clearing_count(list(range(1,14)),q) for q in range(15,32))
print(f"\n  AP: deficiency {apd:.4f}, clearing {apt} (GOOD coverer = FEW clears = the WALL)")
print(f"  spread-DC: deficiency mean {statistics.mean(ds):.4f}, clearing mean {statistics.mean(ts):.1f} (BAD coverer = MANY clears = EASY)")
if len(ds)>1:
    cor=statistics.correlation(ds,ts) if hasattr(statistics,'correlation') else 0
    print(f"  correlation(deficiency, clearing) = {cor:+.3f}  => {'DUALITY CONFIRMED (bad coverer clears)' if cor>0.3 else 'weak'}")
print("\n=> if measS7 (density, three-gap) predicts clearing (liveness), ONE property closes both.")
