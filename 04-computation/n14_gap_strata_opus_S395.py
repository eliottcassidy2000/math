# opus-2026-07-19-S395 -- HYP-7810: TRANSFERRING boxeph-S123's DETERMINANT
# STRATIFICATION FROM THE n=12 GAP TO n=14.
#
# boxeph-S123 stratified the n=12 uniqueness gap (1/13, 2/25) by the reduced
# numerator p of M = D/s, noting p = D/gcd(D,s) <= D so that stratifying by p IS
# stratifying by determinant (crediting my THM-1210).  Their exclusion needs only
# two lines:
#   p = 1: M = 1/q >= 1/13 (LRC(13)) forces q <= 13, so M = 1/13 or M >= 1/12.
#   p = 2: M = 2/q reduced makes q ODD; 2/q >= 1/13 gives q <= 26 hence q <= 25,
#          so M >= 2/25.
# Inputs: a proved lower bound on M, plus PARITY of a reduced denominator.
#
# THE TRANSFER.  At n = 14 the Farey neighbour of 1/14 is 2/27 (1*27 - 2*14 = -1),
# so the stability gap is (1/14, 2/27).  CONDITIONAL ON LRC(14) the same two lines
# should exclude p <= 2, leaving a discrete p >= 3 ladder whose depth-minimal
# member is the mediant 3/41 -- and 41 = 3*14 - 1, matching boxeph's 38 = 3*13 - 1.
from fractions import Fraction as F
from itertools import combinations
import random
LO=F(1,14); HI=F(2,27)
print(f"(1) THE n=14 STABILITY GAP: ({LO}, {HI}) = ({float(LO):.6f}, {float(HI):.6f})")
print(f"    Farey neighbours?  1*27 - 2*14 = {1*27-2*14}   (|.|=1 confirms adjacency)")
print(f"    mediant = (1+2)/(14+27) = {F(3,41)} = {float(F(3,41)):.6f}"
      f"   and 41 = 3*14 - 1 = {3*14-1}")

print()
print("(2) THE NUMERATOR <= 2 EXCLUSION, transferred (conditional on LRC(14))")
print("    p = 1: M = 1/q >= 1/14 forces q <= 14, so M = 1/14 (boundary) or M >= 1/13.")
print(f"           is 1/13 above the gap?  1/13 = {float(F(1,13)):.6f} > 2/27 = {float(HI):.6f} -> "
      f"{'YES' if F(1,13)>HI else 'NO'}")
print("    p = 2: M = 2/q reduced makes q ODD; 2/q >= 1/14 gives q <= 28, q odd so q <= 27,")
print(f"           hence M >= 2/27 = the upper edge -> never strictly inside.")
print("    => conditional on LRC(14), a family in the open gap needs determinant D >= 3.")

print()
print("(3) THE p >= 3 LADDER inside the gap  (p/q with 27p/2 < q < 14p)")
for p in range(1,8):
    vals=[F(p,q) for q in range(1,14*p+1) if F(p,q)>LO and F(p,q)<HI and F(p,q).numerator==p]
    tag = "EXCLUDED (proved above)" if p<=2 else (f"{len(vals)} value(s)")
    print(f"    p={p}: {tag}   {[str(v) for v in vals][:5]}")

print()
print("(4) IS ANYTHING IN THE GAP AT ALL?  random + targeted search over 13-families")
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def M(V):
    D=set()
    for a,b in combinations(V,2):
        D.add(a+b)
        if a!=b: D.add(abs(a-b))
    for v in V: D.add(2*v)
    best=F(0)
    for q in sorted(D):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best=g
    return best
random.seed(395)
hits=[]; n=0
for _ in range(60):
    V=sorted(random.sample(range(1,40),13)); n+=1
    m=M(V)
    if LO<m<HI: hits.append((m,V))
# targeted: perturbations of the two extremals
base=[list(range(1,14)),[1,2,3,4,5,6,7,8,9,10,11,13,24]]
for B in base:
    for i in range(13):
        for r in range(1,60):
            V=sorted(set([x for j,x in enumerate(B) if j!=i]+[r]))
            if len(V)!=13: continue
            n+=1
            m=M(V)
            if LO<m<HI: hits.append((m,V))
print(f"    {n} families tested (random + all single perturbations of both extremals)")
print(f"    families with M strictly inside the gap: {len(hits)}")
for m,V in hits[:4]: print(f"      M={m} at {V}")
print()
print("    (empty => the gap looks genuinely vacant at n=14, matching the n=12 picture)")
