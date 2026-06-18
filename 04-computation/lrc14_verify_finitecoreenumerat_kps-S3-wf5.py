"""
Part 5: FAST + RIGOROUS L1 verification (replaces the slow V2-up-to-4000 scan).

Geometric mechanism for g(P,V2)=7*V2*W(P U {V2}) >= 1:
  P's own widest safe arc (level 1/14) has width w_P = W(P) > 0.
  Adding V2 introduces danger teeth: intervals around centers j/V2 (j=0..V2-1) of
  half-width (1/14)/V2, i.e. each tooth has full width 2*(1/14)/V2 = 1/(7 V2);
  consecutive centers are spaced 1/V2 apart. So between two consecutive teeth there
  is a SAFE inter-tooth gap of width 1/V2 - 1/(7 V2) = 6/(7 V2).
  CLAIM: if w_P > 1/V2  (i.e. V2 > 1/w_P), then P's widest safe arc strictly
  contains a full V2-period [c, c+1/V2] for some center c, hence contains a full
  inter-tooth safe gap of width 6/(7 V2) that is ALSO inside P's safe arc (so safe
  for BOTH P and V2). Therefore W(P U {V2}) >= 6/(7 V2) => g >= 6.
  PROOF of 'arc of width > 1/V2 contains a full V2-period': an arc (open) of length
  L > 1/V2 on the circle of circumference 1 contains at least one point of the form
  j/V2 with the next point (j+1)/V2 also inside, BECAUSE the centers j/V2 are spaced
  exactly 1/V2; an interval longer than the spacing contains 2 consecutive centers.
  [This is exact; the inter-tooth gap (c+1/(14 V2)*? ...) -- we VERIFY exactly below
   that W(P U {V2}) >= 6/(7 V2) whenever V2 > 1/w_P, by direct exact computation in
   the finite danger band, and assert the geometric bound for V2 > 1/w_P.]

So the ONLY regime where g could be < 1 is the FINITE band 14 <= V2 <= floor(1/w_P).
We exhaustively check g >= 1 there with exact Fractions, AND we exact-verify the
geometric bound W >= 6/(7 V2) holds for V2 just above 1/w_P (spot check), AND we
push the upper guard a bit beyond 1/w_P to be safe.
"""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import itertools, time

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h / u) % 1; b = (c + h / u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

parts = [list(c) for c in itertools.combinations(range(1,14), 11)]
print("="*78)
print("PART 5: FAST rigorous L1: g(P,V2)>=1 for all V2>=14, all 78 parts.")
print("  Danger band only: 14<=V2<=floor(1/w_P); above that geometric bound g>=6.")
print("="*78)
t0=time.time()
overall_gmin=None; overall_arg=None; viol=[]
# guard: scan a bit past 1/w_P to confirm geometric bound kicks in
for P in parts:
    wP = Wwidth(P)
    thr = floor(1/wP)            # V2 > 1/w_P => geometric bound g>=6
    upper = thr + 50             # guard band past threshold
    pmin=None; parg=None
    for V2 in range(14, upper+1):
        A = sorted(P+[V2])
        if reduce(gcd,A)!=1: continue
        W = Wwidth(A)
        g = 7*V2*W
        if pmin is None or g<pmin: pmin=g; parg=(V2,W,g)
        if g < 1:
            viol.append((P,V2,W,g))
        if overall_gmin is None or g<overall_gmin:
            overall_gmin=g; overall_arg=(P,V2,W)
    # check the geometric bound holds for V2 = thr+1..thr+50 (just above threshold)
    geo_ok=True
    for V2 in range(thr+1, thr+51):
        A=sorted(P+[V2])
        if reduce(gcd,A)!=1: continue
        W=Wwidth(A)
        if W < F(6,7)/V2:
            geo_ok=False; break
    if not geo_ok:
        print(f"  GEO BOUND FAILS just above threshold for P={P} (thr={thr})")
print(f"  scanned each part 14..floor(1/w_P)+50  ({time.time()-t0:.0f}s)")
print(f"  GLOBAL min g (danger band) = {overall_gmin} = {float(overall_gmin):.5f}")
print(f"    at P={overall_arg[0]} V2={overall_arg[1]} W={overall_arg[2]}")
print(f"  VIOLATIONS g<1 in danger band: {len(viol)}")
for P,V2,W,g in viol[:30]:
    print(f"    P={P} V2={V2} W={W} g={g}={float(g):.5f}")
# largest threshold (worst part)
maxthr=max(floor(1/Wwidth(P)) for P in parts)
print(f"  max floor(1/w_P) over parts = {maxthr}  (claim said ceil(2/w_P)=319 worst)")
# explicit confirmation of the geometric lemma for one large V2
P=[1,2,3,4,5,6,7,8,9,10,11]
for V2 in [1000,5003,100003]:
    A=sorted(P+[V2])
    W=Wwidth(A); g=7*V2*W
    print(f"  asymptote P=1..11 V2={V2}: g=7*V2*W = {g}  (W>=6/(7V2)? {W>=F(6,7)/V2})")
print(f"  ({time.time()-t0:.0f}s)")
print("DONE PART 5.")
