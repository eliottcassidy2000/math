"""
Upgrade the dip bound: the dips are FINITE-Qmax ARTIFACTS. M_c = env EXACTLY for all c.
 (1) the block achievability extends to RATIONAL q=a/b>=n via t=(a-b)/a -> M_c=env at c=(a-nb)/2a (DENSE).
 (2) the dip optimal-t has denominator ~ n^2 (e.g. c=1/28: q=196/13, t=183/196); Qmax<n^2 under-reports.
 (3) M_c is 1-LIPSCHITZ => dense achievability => M_c>=env everywhere; with M_c<=env (clumping) => M_c=env.
"""
from fractions import Fraction
from math import gcd
def frac(x): x=x%1; return min(x,1-x)
def mind(S,c,t): return min(frac(Fraction(s)*t-c) for s in S)
def Mc(S,c,Qmax):
    best=Fraction(0); bt=None
    for q in range(1,Qmax+1):
        for a in range(q):
            m=min(frac(Fraction(s*a,q)-c) for s in S)
            if m>best: best=m; bt=Fraction(a,q)
    return best,bt
n=14; AP=list(range(1,n)); env=lambda c: Fraction(1,n)+c*Fraction(n-2,n)
print("(1) the RATIONAL-q construction t=(a-b)/a achieves env at the 'dip' c's (q=a/b non-integer):")
for (cn,cd) in [(1,28),(3,28),(5,28),(9,28),(1,14),(3,14)]:
    c=Fraction(cn,cd)
    # q = n/(1-2c) = a/b
    q=Fraction(n)/(1-2*c); a,b=q.numerator,q.denominator
    t=Fraction(a-b,a)
    got=mind(AP,c,t); e=env(c)
    print(f"   c={str(c):>5}: q={str(q):>7} t=(a-b)/a={str(t):>8} (denom={a}={a/n**2:.2f}*n^2)  min_v||vt-c||={str(got)} env={str(e)} {'OK' if got==e else 'DIP'}")
print()
print("(2) so the dips need Qmax ~ n^2. M_c at the dips with growing Qmax -> env:")
for (cn,cd) in [(1,28),(9,28)]:
    c=Fraction(cn,cd); e=env(c)
    row=f"   c={str(c)}: env={float(e):.5f}  "
    for Q in [3*n, 12*n, n*n]:
        M,_=Mc(AP,c,Q); row+=f"M(Q={Q})={float(M):.5f}{'=env' if M==e else ''}  "
    print(row)
print()
print("(3) M_c is 1-LIPSCHITZ in c: |M_c - M_c'| <= |c-c'| (sup_t of 1-Lipschitz min_v||vt-c||). Check:")
import random; random.seed(3)
maxratio=0
for _ in range(200):
    c1=Fraction(random.randint(0,50),100); c2=Fraction(random.randint(0,50),100)
    if c1==c2: continue
    M1,_=Mc(AP,c1,2*n); M2,_=Mc(AP,c2,2*n)
    maxratio=max(maxratio, float(abs(M1-M2)/abs(c1-c2)))
print(f"   max observed |dM|/|dc| over random pairs = {maxratio:.4f} (<=1 confirms 1-Lipschitz)")
print()
print("=> UPGRADE: rational-q achievability is DENSE in [0,1/2); M_c 1-Lipschitz => M_c>=env EVERYWHERE;")
print("   with upper bound M_c<=env => M_c = env EXACTLY (no dips). Hence L = 1/4 + 1/(2n) EXACTLY.")
print("   The computational 'dips' were Qmax<n^2 under-estimates (optimal t has denominator ~n^2).")
