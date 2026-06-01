from fractions import Fraction as F
from itertools import combinations
from math import gcd
from lrc import measure_intersection, lonely_measure

# HYPOTHESIS to test: if ALL pairwise correlation ratios r_ij <= 1 (no over-correlated
# pair), is mu > 0 guaranteed? Search for counterexample.
def all_ratios_le1(s,n,tol=F(1)):
    for i,j in combinations(range(len(s)),2):
        r=measure_intersection([s[i],s[j]],n)*F(n*n,4)
        if r>tol: return False,r
    return True,None

print("Search m=4,n=5: among sets with ALL pairwise ratios <=1, is mu always >0?")
n=5; cnt=0; cnt0=0; minmu=None
for s in combinations(range(1,20),4):
    ok,_=all_ratios_le1(s,n)
    if ok:
        cnt+=1
        mu=lonely_measure(list(s),n)
        if mu==0: cnt0+=1; print("   COUNTEREXAMPLE mu=0 with all r<=1:",s)
        if minmu is None or mu<minmu: minmu=mu; arg=s
print(f"  sets with all r<=1: {cnt}; of those mu=0: {cnt0}; min mu={minmu} at {arg if cnt else None}")

print("\nSame for m=5,n=6, speeds in 1..18:")
n=6; cnt=0; cnt0=0; minmu=None; arg=None
for s in combinations(range(1,19),5):
    ok,_=all_ratios_le1(s,n)
    if ok:
        cnt+=1
        mu=lonely_measure(list(s),n)
        if mu==0: cnt0+=1; print("   COUNTEREXAMPLE mu=0 with all r<=1:",s)
        if minmu is None or mu<minmu: minmu=mu; arg=s
print(f"  sets with all r<=1: {cnt}; of those mu=0: {cnt0}; min mu={minmu} at {arg}")

# Also test STRICT condition r<1 strictly (no resonance at all).
print("\nWith all ratios STRICTLY <1 (m=4,n=5):")
n=5; cnt=0; cnt0=0
for s in combinations(range(1,20),4):
    bad=False
    for i,j in combinations(range(4),2):
        if measure_intersection([s[i],s[j]],n)*F(n*n,4)>=1: bad=True;break
    if not bad:
        cnt+=1
        if lonely_measure(list(s),n)==0: cnt0+=1; print("  ce:",s)
print(f"  sets all r<1: {cnt}, of those mu=0: {cnt0}")
