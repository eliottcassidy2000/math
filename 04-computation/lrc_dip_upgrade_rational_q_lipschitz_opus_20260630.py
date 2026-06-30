"""
RIGOROUS rational-q achievability (the clean min computation):
 q=a/b>=n (gcd(a,b)=1), t=(a-b)/a, c=(q-n)/2q=(a-nb)/2a.  Then
   ||vt - c|| = ||(a+wb)/2a||  with w=2v-n in {-(n-2),..,n-2} step 2  [since vt-c = -(vb/a + c) and algebra]
 a+wb in [a-(n-2)b, a+(n-2)b]; nearest multiple of 2a is 0 or 2a; dist = a-|w|b, minimized at |w|=n-2:
   min_v ||vt-c|| = (a-(n-2)b)/2a = (q-n+2)/2q = ENV.   Test across n and rational q.
"""
from fractions import Fraction
from math import gcd
def frac(x): x=x%1; return min(x,1-x)
def mind(S,c,t): return min(frac(Fraction(s)*t-c) for s in S)
ok_all=True
for n in [8,10,12,14,16,20]:
    AP=list(range(1,n))
    for (a,b) in [(n,1),(2*n+1,2),(3*n+2,3),(n*n,n-1),(5*n+1,5),(7*n+3,7)]:
        if gcd(a,b)!=1 or Fraction(a,b)<n: continue
        q=Fraction(a,b); c=(q-n)/(2*q); t=Fraction(a-b,a)
        env=(q-n+2)/(2*q)
        got=mind(AP,c,t)
        ok = (got==env) and (0<=c<=Fraction(1,2))
        ok_all = ok_all and ok
        if not ok: print(f"  FAIL n={n} q={q}: got {got} env {env}")
print(f"rational-q construction min_v||vt-c||=env for ALL tested (n=8..20, many q): {ok_all}")
print()
print("PROOF (elementary): with w=2v-n, min_v||vt-c|| = min_{|w|<=n-2} (a-|w|b)/2a = (a-(n-2)b)/2a = env.")
print("The achieving c=(a-nb)/2a range DENSELY over [0,1/2) as a/b ranges over rationals >= n. QED achievability.")
print()
print("=> THE DIP UPGRADE (rigorous):")
print("   (i)  achievability M_c=env at the DENSE set c=(q-n)/2q, q in Q, q>=n  [proven above];")
print("   (ii) M_c is 1-LIPSCHITZ in c (sup_t of 1-Lipschitz min_v||vt-c||);")
print("   (iii) (i)+(ii) => M_c >= env EVERYWHERE (no dips); with upper bound M_c<=env => M_c = env EXACTLY.")
print("   (iv)  hence L = int_0^1 M_c dc = 1/4 + 1/(2n) EXACTLY (the O(1/n^2) error term is REMOVED).")
print("   The computational 'dips' were Qmax<n^2 under-estimates: the optimal t at c=odd/2n has denominator n^2.")
