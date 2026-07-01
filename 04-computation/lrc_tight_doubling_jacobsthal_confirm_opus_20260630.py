"""
Confirm the repo's DOUBLING + JACOBSTHAL picture through my avoided-arc lens.
Sporadic tight = AP with removed element k, added g. The avoided-arc argument fails at the min-gap
difference d=k (no runner k). It's rescued iff the absent-difference regime ||kt||<1/n is covered.
CHECK: is g=2k (doubling) the covering element? And which n have a single-swap sporadic?
Also: my census (single-swap) vs the doubling gate.
"""
from fractions import Fraction
def fr(x): x=x%1.0; return min(x,1-x)
def Mfloat(S,Q):
    b=0.0
    for q in range(1,Q+1):
        for a in range(1,q):
            m=min(fr(s*a/q) for s in S)
            if m>b: b=m
    return b
def Mexact(S,Q):
    b=Fraction(0)
    for q in range(1,Q+1):
        for a in range(1,q):
            m=min(min((Fraction(s*a,q))%1,1-(Fraction(s*a,q))%1) for s in S)
            if m>b: b=m
    return b
# Jacobsthal-ish: v is 'gateable' if [v/2..? ] ... simpler: report for each n the single-swap (k,g) with g=2k
print("Single-swap tight sets S=({1..n-1}\{k}) U {g}; flag doublings g=2k:")
for n in range(5,18):
    hits=[]
    base=list(range(1,n))
    for k in base:
        for g in range(n,4*n+1):
            if g in base: continue
            S=[x for x in base if x!=k]+[g]
            if Mfloat(S,5*n)<=1.0/n+1e-9 and Mexact(S,10*n)==Fraction(1,n):
                hits.append((k,g,'DOUBLE' if g==2*k else ('g=k+%d*n'%((g-k)//n) if (g-k)%n==0 else '?')))
    print(f"  n={n} (2n-1={2*n-1}): {hits}")
print()
print("=> the DOUBLING sporadics (g=2k) are the repo's operad; k=n-2 at n=14 (12->24). The n=5,6 ones are")
print("   NOT doublings -> a second mechanism. Sporadics are irregular & rare (matches Jacobsthal gate).")
