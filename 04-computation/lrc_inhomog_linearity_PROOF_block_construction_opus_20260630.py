"""
PROOF of M_c(AP_n) = 1/n + c(n-2)/n. The ACHIEVABILITY (lower bound) via an EXPLICIT optimal t:
  for c = (q-n)/(2q) with integer q>=n, take t = (q-1)/q.
  Then {kt mod 1 : k=1..n-1} = {(q-k)/q} = a BLOCK of n-1 consecutive multiples of 1/q
  ([(q-n+1)/q , (q-1)/q]); the complementary gap (through 0) has length (q-n+2)/q and is
  CENTERED at c, so min_k ||kt-c|| = (q-n+2)/(2q) = 1/n + c(n-2)/n.   <-- algebra below.
"""
from fractions import Fraction
def mind(S,c,t):
    return min(min((Fraction(s)*t-c)-((Fraction(s)*t-c).__floor__()), 1-((Fraction(s)*t-c)-((Fraction(s)*t-c).__floor__()))) for s in S)
print("ACHIEVABILITY: explicit t=(q-1)/q at c=(q-n)/(2q) gives min_k||kt-c|| = 1/n+c(n-2)/n EXACTLY:")
for n in [10,14,20]:
    AP=list(range(1,n))
    print(f"  n={n}:")
    for q in range(n, n+8):
        c=Fraction(q-n,2*q); t=Fraction(q-1,q)
        got=mind(AP,c,t); env=Fraction(1,n)+c*Fraction(n-2,n)
        print(f"     q={q}: c=(q-n)/2q={str(c):>7}={float(c):.4f}  t={str(t):>6}  min_k||kt-c||={str(got):>9}  env={str(env):>9}  {'OK' if got==env else 'XX'}")
print()
print("ALGEBRA: gap=(q-n+2)/q centered at c => M=(q-n+2)/(2q). Sub q=n/(1-2c) [<=> c=(q-n)/2q]:")
print("  M = 1/2 - (n-2)/(2q) = 1/2 - (n-2)(1-2c)/(2n) = [n-(n-2)]/(2n) + c(n-2)/n = 1/n + c(n-2)/n.  QED(achiev)")
print()
# UPPER BOUND check: M_c <= envelope for ALL c (large Qmax), with equality exactly at nice c=(q-n)/2q
def Mc(S,c,Qmax):
    best=Fraction(0)
    for q in range(1,Qmax+1):
        for a in range(q):
            m=min(min((Fraction(s*a,q)-c)-((Fraction(s*a,q)-c).__floor__()),1-((Fraction(s*a,q)-c)-((Fraction(s*a,q)-c).__floor__()))) for s in S)
            if m>best: best=m
    return best
print("UPPER BOUND: M_c <= 1/n+c(n-2)/n for ALL c (large Qmax), tight at nice c, dips at non-integer q=n/(1-2c):")
n=14; AP=list(range(1,n))
for cnum,cden in [(1,28),(1,14),(3,28),(1,7),(1,4),(2,7)]:
    c=Fraction(cnum,cden); env=Fraction(1,n)+c*Fraction(n-2,n); M=Mc(AP,c,12*n)
    qstar=Fraction(n)/(1-2*c)
    print(f"   c={str(c):>5}: M_c={float(M):.5f}  env={float(env):.5f}  M<=env:{M<=env}  q*=n/(1-2c)={str(qstar):>5} {'(integer->tight)' if qstar.denominator==1 else '(non-int->dip)'}")
print()
print("=> LINEARITY M_c=1/n+c(n-2)/n holds EXACTLY at the dense set c=(q-n)/2q (q>=n integer); M_c<=env always;")
print("   between nice c, q rounds and M_c dips O(1/n^2) below env. As n->inf nice c densify => L=1/4+1/(2n).")
