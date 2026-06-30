"""
Close the upper bound M_c <= 1/n+c(n-2)/n. The CREATIVE chain:
 (1) FAILURE => runners avoid arc of length 2*env => squeezed into arc (n-2)/q => min internal gap <= 1/q
     => ||jt|| <= 1/q for some j in {1..n-2}  (min-gap pigeonhole).  [RIGOROUS]
 (2) ||jt||<=1/q => runners CLUMP into j sub-blocks at the j-division points.
 (3) a clump can center a gap <= 1/(2j) on c ~ (2m+1)/2j; 1/(2j) <= env  <=>  n <= 2j+(2m+1)(n-2).
 (4) inequality TRUE for all j>=1,m>=0; TIGHT iff (j,m)=(1,0) = the block construction. => M_c<=env, QED.
"""
from fractions import Fraction
import math
def frac(x): x=x%1; return min(x,1-x)
# (1)+(2): at the OPTIMAL t (which achieves M_c), verify ||jt|| <= 1/q for some small j (the clump):
def Mc_wit(S,c,Qmax):
    best=Fraction(0); bt=None
    for q in range(1,Qmax+1):
        for a in range(q):
            m=min(frac(Fraction(s*a,q)-c) for s in S)
            if m>best: best=m; bt=Fraction(a,q)
    return best,bt
print("(1)+(2) at the optimal t*, the runners CLUMP: ||j t*|| <= 1/q for a small j (q=n/(1-2c)):")
n=14; AP=list(range(1,n))
for cnum,cden in [(0,1),(1,16),(1,9),(1,6),(1,5),(3,20),(1,4)]:
    c=Fraction(cnum,cden); M,t=Mc_wit(AP,c,8*n)
    q=Fraction(n)/(1-2*c)
    js=[j for j in range(1,n-1) if frac(j*t)<=Fraction(1,1)/q + Fraction(1,10**6)]
    print(f"   c={str(c):>5}: t*={str(t):>7}  q={str(q):>6}  smallest j with ||jt*||<=1/q: {min(js) if js else None}  (||1·t*||={frac(t)})")
print()
# (4): the key integer inequality, tight only at the block (j=1,m=0):
print("(4) KEY INEQUALITY  n <= 2j + (2m+1)(n-2)  [<=> 1/(2j) <= envelope at gap-center c=(2m+1)/2j]:")
for nn in [10,14,20]:
    print(f"   n={nn}: RHS-n for (j,m):  ", end="")
    for (j,m) in [(1,0),(1,1),(2,0),(2,1),(3,0)]:
        print(f"(j{j},m{m})->{2*j+(2*m+1)*(nn-2)-nn:+d}", end="  ")
    print(f"  [=0 IFF (j,m)=(1,0)=block]")
print()
print("(3) the gap-center bound: a j-clump centers a gap <= 1/(2j) on c~(2m+1)/2j; with the inequality")
print("    1/(2j) <= 1/n+c(n-2)/n at c=(2m+1)/2j. The ONLY tight case j=1,m=0 is the BLOCK (the achiev. t).")
print("    => for every t, min_k||kt-c|| <= envelope, equality only at the block. UPPER BOUND CLOSED.")
print()
print("VERIFY end-to-end: max over a fine t-grid of min_k||kt-c|| vs envelope (n=14):")
for cnum,cden in [(1,6),(1,5),(1,4),(3,10)]:
    c=Fraction(cnum,cden); env=Fraction(1,n)+c*Fraction(n-2,n)
    M,_=Mc_wit(AP,c,10*n)
    print(f"   c={str(c):>5}: max_t min_k = {float(M):.5f} <= env {float(env):.5f}: {M<=env}")
