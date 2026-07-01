from fractions import Fraction
from itertools import combinations
def frac(x): x=x%1; return min(x,1-x)
def Mval(S,Qmax):
    best=Fraction(0)
    for q in range(1,Qmax+1):
        for a in range(1,q):
            m=min(frac(Fraction(s*a,q)) for s in S)
            if m>best: best=m
    return best
def diff_closed(S):
    Sset=set(S)
    return all((abs(x-y) in Sset) for x in S for y in S if x!=y)
# (1) primitive difference-closed sets of size n-1 -> only the AP (no M needed, fast)
print("(1) primitive (contains 1) difference-closed sets of size n-1 in [1..2n]:")
for n in [5,6,7,8]:
    found=[S for S in combinations(range(1,2*n), n-1) if 1 in S and diff_closed(S)]
    print(f"   n={n}: {found}  (AP {tuple(range(1,n))} only? {found==[tuple(range(1,n))]})")
print("   => PROVED: the only primitive difference-closed (n-1)-set is the AP {1..n-1}.")
print("      (1 in S; a2-1 in S,<a2 => a2=2; a_k-1 in S,<=k-1 => a_k<=k, >k-1 => a_k=k.)")
print()
# (2) dilated APs: difference-closed + tight
print("(2) dilated APs d*{1..n-1} (n=7): difference-closed AND M=1/n?")
n=7
for d in [1,2,3]:
    S=[d*k for k in range(1,n)]; M=Mval(S,4*n)
    print(f"   d={d}: diff-closed={diff_closed(S)}, M={M}=1/n? {M==Fraction(1,n)}")
print()
# (3) Goddyn-Wong {1..11,13,24} (n=14): NOT diff-closed, yet tight -> the OTHER tight component
print("(3) Goddyn-Wong GW={1..11,13,24} (n=14): tight but NOT difference-closed:")
GW=list(range(1,12))+[13,24]
print(f"   |GW|={len(GW)} (n=14 needs 13 speeds); difference-closed? {diff_closed(GW)}")
missing=[abs(x-y) for x in GW for y in GW if x!=y and abs(x-y) not in set(GW)]
print(f"   differences NOT in GW (e.g. 12=13-1): sample {sorted(set(missing))[:8]}")
M=Mval(GW, 200)
print(f"   M(GW) with Qmax=200 = {M} = {float(M):.5f} (1/14={1/14:.5f}); tight={M<=Fraction(1,14)}")
print()
print("=> The avoided-arc-edge argument proves M(S)<=1/n for DIFFERENCE-CLOSED S = the dilated APs (the unique")
print("   primitive one is {1..n-1}). This is ONE provable component of the tight locus. Goddyn-Wong is tight")
print("   but NOT difference-closed (12=13-1 absent) -- its min-gap difference 12 is not a runner, so the")
print("   argument needs a DIFFERENT covering for the '12t' gap. That is the residual for OPEN-Q-108.")
