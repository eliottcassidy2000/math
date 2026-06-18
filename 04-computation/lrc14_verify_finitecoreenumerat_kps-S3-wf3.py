"""
Part 3: ATTACK L1 (drop-max scaling) and the V2>=51 -> 7*V2*W(A)>=1 claim.

L1 says: for A = P U {V2}, the product g(P,V2) = 7*V2*W(A) -> 6 as V2->inf,
and >= 1 for all V2 >= 51. The prompt admits L1 is verified only up to V2~319/400,
NOT proved. We attack:
  (a) Is g(P,V2) >= 1 for ALL V2 >= 14 (not just 51)? Find global min over V2.
  (b) Scan V2 to LARGE values (up to several thousand) for every 11-subset part P,
      hunt for any dip g < 1. If found for some V2 with Vmax>V2 you'd need a real S
      with Vmax slightly above; that's a candidate counterexample direction.
  (c) Verify the geometric mechanism: when W(P) (the part's own widest arc) spans
      >= 2 consecutive V2-tooth-centers (gap 1/V2), a full inter-tooth safe gap of
      width 6/(7 V2) survives so W(A) >= 6/(7 V2), giving g >= 6. The claim's
      'worst ceil(2/w_P)=319'. Verify w_P=W(P) per part and the 2/w_P threshold.

NOTE: W(A) computed exactly. For the >=63 branch the relevant question is the
minimum of 7*Vmax*W(P U {V2}) over admissible (V2 < Vmax, Vmax>=63). Since
Vmax > V2 and the product uses Vmax not V2, having 7*V2*W >= 1 plus Vmax>V2 gives
7*Vmax*W > 7*V2*W >= 1 ONLY IF W is computed with the SAME A. But A=P U {V2} does
NOT contain Vmax. So 7*Vmax*W(P U {V2}) = (Vmax/V2)*(7*V2*W(P U {V2})) >
7*V2*W(P U{V2}) >= 1. GOOD: the factor Vmax/V2 > 1 helps. So L1's '7*V2*W>=1'
DOES imply '7*Vmax*W(A)>1' when Vmax>V2. The load-bearing claim is purely
'7*V2*W(P U {V2}) >= 1 for all V2 in the >=63-branch's V2-range (14..inf)'.
We must verify g(P,V2) >= 1 for ALL V2 >= 14 (since with Vmax>=63 we still might
have V2 as small as 14). For V2<=62 part 2A already showed 7*63*W>1 i.e. with the
factor. But for V2 in 14..62 and Vmax exactly 63, we used W_min. For V2>=63 the
set has V2>=63 so BOTH large speeds >=63... that's still k=2, Vmax>V2>=63. Need g>=1.
So overall: prove g(P,V2)=7*V2*W(P U {V2}) >= 1 for ALL V2>=14, ALL 78 parts P.
(Then 7*Vmax*W(A) = (Vmax/V2)*g >= g >= 1, strict since Vmax>V2.)
We hunt for any (P,V2) with g < 1.
"""
from fractions import Fraction as F
from math import gcd
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
print("PART 3: g(P,V2)=7*V2*W(P U {V2}); need g>=1 for all V2>=14, all 78 parts.")
print("="*78)
t0=time.time()
# first the part's own widest arc w_P = W(P)
wP = {tuple(P): Wwidth(P) for P in parts}
print("  w_P = W(P) for each 11-subset part (should all be > 1/7 since M(P)>=1/12):")
mn_wP = min(wP.values())
print(f"    min w_P over 78 parts = {mn_wP} = {float(mn_wP):.5f};  2/min_wP = {float(2/mn_wP):.2f}")
# global min of g over V2 in 14..VMAX_SCAN
VMAX_SCAN = 4000
glob_g_min = None; glob_arg=None
viol=[]
per_part_min={}
for P in parts:
    pmin=None; parg=None
    for V2 in range(14, VMAX_SCAN+1):
        A = sorted(P+[V2])
        if reduce(gcd,A)!=1: continue
        W = Wwidth(A)
        g = 7*V2*W
        if pmin is None or g<pmin: pmin=g; parg=(P,V2,W)
        if g < 1:
            viol.append((P,V2,W,g))
        if glob_g_min is None or g<glob_g_min: glob_g_min=g; glob_arg=(P,V2,W)
    per_part_min[tuple(P)]=(pmin,parg)
print(f"  scanned V2 in 14..{VMAX_SCAN} for all 78 parts ({time.time()-t0:.0f}s)")
print(f"  GLOBAL min g = {glob_g_min} = {float(glob_g_min):.5f}")
print(f"    at P={glob_arg[0]} V2={glob_arg[1]} W={glob_arg[2]}")
print(f"  VIOLATIONS g<1 : {len(viol)}")
for P,V2,W,g in viol[:30]:
    print(f"    P={P} V2={V2} W={W} g={g}={float(g):.5f}")
# where does the min occur (small V2 or large)?
small_V2_mins = sorted(((float(v[0]), tuple(v[1][0]), v[1][1]) for v in per_part_min.values()))[:5]
print("  5 smallest per-part g-minima:")
for g,P,V2 in small_V2_mins:
    print(f"    g={g:.5f}  P={list(P)} V2={V2}")
print(f"  ({time.time()-t0:.0f}s)")
print("DONE PART 3.")
