"""
Integrate mac-mini HYP-3701: n/Phi_6 is NOT the universal covering-min (small n use 2/(2n-1) via drop-2);
beaters don't scale. Verify small-n beaters, the n=14 drop-2 family (should be >>14/183), and the transition.
"""
import math
from fractions import Fraction
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
def is_cov(S,n): return all(any(x%qq==0 for x in S) for qq in range(2,n+1))
print("mac-mini's small-n drop-2 beaters (covering-min = 2/(2n-1)) vs construction n/Phi_6(n):")
beaters={4:[1,3,4],5:[1,3,4,5],6:[1,3,4,5,18]}
for n,S in beaters.items():
    Phi6=n*n-n+1; cons=list(range(1,n-1))+[(n-1)*n]
    Mb=M_exact(S) if is_cov(S,n) else None
    Mc=M_exact(cons) if is_cov(cons,n) else None
    print(f"  n={n}: drop-2 {S} M={Mb}={float(Mb):.4f} (2/(2n-1)={2}/{2*n-1}={2/(2*n-1):.4f}); construction n/Phi_6={n}/{Phi6}={n/Phi6:.4f} -> covering-min={'drop-2' if Mb<Mc else 'construction'}")
print()
print("the formula comparison: n/Phi_6 vs 2/(2n-1): n(2n-1)-2Phi_6 = n-2, so 2/(2n-1) < n/Phi_6 for ALL n>2")
print("  => IF drop-2 achieved 2/(2n-1) it would ALWAYS win. The content (mac-mini): it DOESN'T SCALE.")
print()
print("n=14 drop-2 family (skip element 2, add tuned large 14m) -- does it beat 14/183=0.07650?")
for big in [14,28,42,56,70,84,98,112,140,168,196]:
    S=sorted(set([1,3,4,5,6,7,8,9,10,11,12,13,big]))
    if len(S)!=13 or not is_cov(S,14): continue
    M=M_exact(S)
    print(f"  drop-2 + {big:>3}: M={str(M):>7}={float(M):.5f}  {'BEATS 14/183!' if M<Fraction(14,183) else '> 14/183 (construction stands)'}")
print()
print(f"  2/(2n-1) at n=14 = 2/27 = {2/27:.5f} (< 14/183={14/183:.5f}) -- the would-be beater IF it scaled")
print(f"  => drop-2 does NOT reach 2/27 at n=14; construction 14/183 stands. covering-min is n-DEPENDENT.")
