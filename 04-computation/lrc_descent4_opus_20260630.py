"""
Refinement: is the COVERING-MIN off-cusp (lower than any cusp covering set)? Scan off-cusp covering sets
near the construction {1..12,182}. Compare cusp-covering-min (2/23=0.087) vs off-cusp construction (14/183).
"""
import math
from fractions import Fraction
from itertools import combinations
def M_exact(S,Qmax=400):
    best=Fraction(0)
    for q in range(2,Qmax+1):
        bb=0
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            m=q
            for s in S:
                r=(s*a)%q; d=r if r<=q-r else q-r
                if d<m:m=d
                if m<=bb:break
            if m>bb:bb=m
        v=Fraction(bb,q)
        if v>best:best=v
    return best
def is_cov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
def cusp(S): return set(x%7 for x in S if x%2==1)==set(range(7))
# off-cusp covering construction family: {1..12, K} with K a covering completion for q=13,14 etc.
print("OFF-CUSP covering sets of form {1..12, K} (K provides missing mult of 13/14):")
best=Fraction(1); bestS=None
for K in [182,182*2,84,168,156,364,14*13,12*13,154,170,196]:
    S=list(range(1,13))+[K]
    if not is_cov(S): continue
    c=cusp(S); M=M_exact(S)
    tag="CUSP" if c else "offcusp"
    print(f"   K={K:>4} ({tag:>7}): M={str(M):>8} = {float(M):.5f}")
    if not c and M<best: best=M; bestS=S
print(f"\n   off-cusp covering-min in this family = {best} = {float(best):.5f} at {bestS}")
print(f"   compare: cusp-covering-min (tested) = 2/23 = {2/23:.5f}; off-cusp 14/183 = {14/183:.5f}")
print()
print("THE REFINEMENT (corrects earlier 'cusp=hard' emphasis):")
print(f"   * AP (cusp, NON-covering): M=1/14={1/14:.5f}  -- GLOBAL extremal, measure 0, comb-witness.")
print(f"   * off-cusp covering-min:   M=14/183={14/183:.5f} -- the BINDING covering constraint (n/Phi6(n)).")
print(f"   * cusp covering sets:      M>=2/23={2/23:.5f}   -- EASIER (higher M).")
print("   => the CONJECTURE (covering-min >= 1/n) binds OFF-CUSP, not at the cusp. The cusp (measure 0) is")
print("      the global extremal via the comb (AP), but among COVERING sets the tightest M is OFF-cusp.")
print("   => cusp = measure phenomenon (existence/comb); covering-min = M phenomenon (the n/Phi6(n) floor).")
print("      DIFFERENT extremes. The descent product (measure) and the floor M are genuinely distinct.")
