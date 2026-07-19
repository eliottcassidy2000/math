# opus-2026-07-19-S397 -- HYP-7830: SETTLING (1/14, 3/41).
#
# M = D/s lies strictly inside (1/14, 3/41) iff
#       1/14 < D/s < 3/41   <=>   41D/3 < s < 14D.
# That pins the admissible (D,s) pairs EXACTLY, and the small-D cases die on
# arithmetic alone:
#     D=1: need 13.67 < s < 14  -> no integer.
#     D=2: need 27.33 < s < 28  -> no integer.
#     D=3: need 41    < s < 42  -> no integer.
# So the interval FORCES D >= 4 -- strengthening boxeph-S123's D >= 3 (which was
# for the wider Farey interval) by one full step, on exact arithmetic.
from fractions import Fraction as F
from math import gcd, ceil, floor
from itertools import combinations
import random
LO=F(1,14); HI=F(3,41)
print("(1) THE EXACT CANDIDATE LIST: integers s with 41D/3 < s < 14D")
cands=[]
for D in range(1,26):
    lo=F(41*D,3); hi=14*D
    ss=[s for s in range(int(lo)+1, hi) if F(D,s)>LO and F(D,s)<HI]
    if ss: cands.extend((D,s) for s in ss)
    tag = ", ".join(f"s={s} -> {F(D,s)}" for s in ss) if ss else "NONE (no integer in range)"
    if D<=10: print(f"    D={D:2d}: ({float(lo):.2f}, {hi})  {tag}")
print(f"    ...")
print(f"    => D = 1, 2, 3 are EXCLUDED BY ARITHMETIC ALONE (empty integer range).")
print(f"       The interval forces D >= 4, one step beyond boxeph's D >= 3.")

print()
print("(2) THE STRUCTURE OF A REALISING FAMILY")
print("    At the maximiser t* = p/s (p coprime to s), M = D/s means every speed v")
print("    satisfies  min(vp mod s, s - vp mod s) >= D,  with equality for the")
print("    active pair.  Since p is invertible we may take p = 1 after relabelling,")
print("    so the speeds' residues mod s must lie in the BAND [D, s-D], one hitting D.")
for D,s in cands[:6]:
    band=s-2*D+1
    print(f"    D={D:2d}, s={s:3d}: band [{D},{s-D}] has {band} of {s} residues "
          f"({100*band/s:.0f}%) -- residues are NOT the obstruction")

print()
print("(3) TARGETED CONSTRUCTION: build families with residues in the band and")
print("    an active pair summing to s, then test whether M is EXACTLY D/s")
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def Mfull(V):
    Dn=set()
    for a,b in combinations(V,2):
        Dn.add(a+b)
        if a!=b: Dn.add(abs(a-b))
    for v in V: Dn.add(2*v)
    best=F(0)
    for q in sorted(Dn):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best=g
    return best
random.seed(397)
for D,s in cands[:4]:
    target=F(D,s)
    band=[u for u in range(D, s-D+1)]
    hits=0; tried=0; best_seen=F(1)
    for _ in range(400):
        # active pair: v_i with v_i mod s = D-ish, v_j = s - v_i
        vi=random.choice([u for u in band if u<=s//2]) if band else None
        if vi is None or vi<=0: continue
        vj=s-vi
        if vj<=0 or vj==vi: continue
        rest=set()
        while len(rest)<11:
            u=random.choice(band)
            k=random.randint(0,3)
            v=u+k*s
            if v>0 and v not in (vi,vj): rest.add(v)
        V=sorted({vi,vj}|rest)
        if len(V)!=13: continue
        tried+=1
        m=Mfull(V)
        if m<best_seen: best_seen=m
        if m==target: hits+=1
    print(f"    D={D:2d}, s={s:3d}  target {target} ({float(target):.6f}): "
          f"{tried} built, {hits} hit the target; smallest M seen {best_seen} "
          f"({float(best_seen):.6f})")
print()
print("    (a hit would REFUTE emptiness; smallest-M-seen shows how close the")
print("     construction gets before some other point beats the intended maximiser)")
