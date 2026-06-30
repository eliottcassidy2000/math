"""THE EVEN BLOCK: 2*{1..n-1} is covering (even n) with M=M({1..n-1})=1/n. The covering-min = 1/n (TIGHT),
NOT n/Phi_6. The construction was a red herring. Verify at n=8,10,12,14; show odd n differs."""
import math
from fractions import Fraction
def M_exact(S,Qmax):
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
def is_cov(S,n): return len(set(S))==n-1 and all(any(x%q==0 for x in S) for q in range(2,n+1))
print("THE EVEN BLOCK 2*{1..n-1} vs the construction n/Phi_6:")
for n in [7,8,9,10,12,14]:
    Phi6=n*n-n+1; constr=Fraction(n,Phi6); Q=2*Phi6
    eb=[2*k for k in range(1,n)]           # the even block = 2*{1..n-1}
    ap=list(range(1,n))                     # the AP {1..n-1} (non-covering, extremal M=1/n)
    cov_eb=is_cov(eb,n)
    M_eb=M_exact(eb,Q) if cov_eb else None
    M_ap=M_exact(ap,Q)
    par="EVEN" if n%2==0 else "ODD"
    print(f"  n={n:>2} ({par}): AP {{1..{n-1}}} M={M_ap}={float(M_ap):.5f}=1/n ; even-block 2*AP covering={cov_eb}"
          + (f", M={M_eb}={float(M_eb):.5f}" if cov_eb else " (NOT covering: no mult of odd n)"))
    print(f"            construction n/Phi_6={constr}={float(constr):.5f} ; covering-min <= {'1/n (even block!)' if cov_eb else '2/(2n-1) or higher'}")
print()
print("=> EVEN n: covering-min = 1/n EXACTLY (even block = AP doubled, COVERING). The conjecture is TIGHT.")
print("   ODD n: even block fails (no mult of odd n); covering-min > 1/n (n=7: 2/13).")
print("   n=14 is EVEN => covering-min = 1/14, NOT 14/183. The construction (n/Phi_6) is a RED HERRING.")
print("   The 'razor-thin margin n/Phi_6 - 1/n ~ 1/n^2' was about the WRONG (non-extremal) family.")
