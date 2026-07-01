"""
Enumerate the TIGHT LOCUS (primitive (n-1)-sets S with M(S)=1/n) for small n; classify.
 - filter with floats (Qmax=5n), confirm candidates with exact Fractions (Qmax=big).
 - classify: difference-closed (=AP)? or 'AP with element k replaced by multiple mk'?
"""
from fractions import Fraction
from itertools import combinations
def fr(x): x=x%1.0; return min(x,1-x)
def Mfloat(S,Qmax):
    best=0.0
    for q in range(1,Qmax+1):
        for a in range(1,q):
            m=min(fr(s*a/q) for s in S)
            if m>best: best=m
    return best
def Mexact(S,Qmax):
    from fractions import Fraction as F
    best=F(0)
    for q in range(1,Qmax+1):
        for a in range(1,q):
            m=min(min((F(s*a,q))%1,1-(F(s*a,q))%1) for s in S)
            if m>best: best=m
    return best
def diff_closed(S):
    Ss=set(S); return all(abs(x-y) in Ss for x in S for y in S if x!=y)
def ap_replace_form(S,n):
    # is S = ({1..n-1} minus one k) plus one big element g, with g a multiple of some absent difference?
    base=set(range(1,n)); Ss=set(S)
    removed=sorted(base-Ss); added=sorted(Ss-base)
    return removed,added
for n in [5,6,7]:
    N=3*n; tights=[]
    for S in combinations(range(2,N), n-2):  # contains 1, plus n-2 more from [2,N)
        S=(1,)+S
        if Mfloat(S,5*n) < 1.0/n + 1e-9:  # candidate M<=1/n (tight, since M>=1/n by LRC)
            Me=Mexact(list(S), 8*n)
            if Me==Fraction(1,n): tights.append(S)
    print(f"n={n}: {len(tights)} primitive tight sets (M=1/n), max elt <= {N-1}:")
    for S in tights:
        dc=diff_closed(S); rem,add=ap_replace_form(S,n)
        tag="AP (diff-closed)" if dc else f"AP minus {rem} plus {add}"
        print(f"    {S}  [{tag}]")
    print()
