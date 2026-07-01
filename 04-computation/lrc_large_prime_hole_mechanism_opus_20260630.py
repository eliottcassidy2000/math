"""
Make 'large primes forced' rigorous. q prime in (n/2, n-1], remove q from AP: where is the hole (t*=lonely
time, M>1/n)? does adding 2q patch? WHY not (for prime q)? Find the specific escape time t*.
"""
from fractions import Fraction
def frac(x): x=x%1; return min(x,1-x)
def M_wit(S,Qmax):
    best=Fraction(0); bt=None
    for q in range(1,Qmax+1):
        for a in range(1,q):
            m=min(min((Fraction(s*a,q))%1,1-(Fraction(s*a,q))%1) for s in S)
            if m>best: best=m; bt=Fraction(a,q)
    return best,bt
def isprime(m): return m>1 and all(m%d for d in range(2,int(m**.5)+1))
print("q prime in (n/2,n-1]: AP\{q} lonely time t* and M; +2q patch?; the escape structure:")
for n in [10,12,14,18]:
    primes=[q for q in range(n//2+1,n) if isprime(q)]
    for q in primes:
        AP=list(range(1,n))
        Sm=[v for v in AP if v!=q]
        M,t=M_wit(Sm, 6*n)
        Sp=[v for v in AP if v!=q]+[2*q]
        Mp,tp=M_wit(Sp, 6*n)
        # at t*, which AP-runner is closest to 0 (the one that WOULD bind if present)?
        near=sorted(AP, key=lambda v: frac(v*t))[:3]
        print(f"  n={n} q={q}(prime): AP\q -> M={M}={float(M):.4f}>1/n={1/n:.4f}, t*={t}; +2q -> M={Mp}={float(Mp):.4f} (tight:{Mp==Fraction(1,n)})")
        print(f"       at t*={t}: 3 nearest AP-runners to 0: {[(v, str(frac(v*t))) for v in near]}  (q={q} would be at {frac(q*t)})")
