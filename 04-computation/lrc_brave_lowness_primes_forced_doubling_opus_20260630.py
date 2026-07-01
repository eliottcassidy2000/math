"""
BRAVE attack on the lowness: tight => S contains all but O(polylog) of {1..n-1}.
Mechanism: tight=>covering (mult of every q, q-witness). Swapping q->multiple must patch H_q, JACOBSTHAL-gated.
KEY: the gate is CLOSED for primes (g(prime)=2, no 2 consecutive nonunits) => LARGE PRIMES FORCED.
Test: (1) census swappable k are composite w/ Jacobsthal run g(k)>=n-k+2; (2) large primes never swapped;
(3) count swappable k; (4) predict doubling at n=20 (18->36).
"""
from fractions import Fraction
from math import gcd
def frac(x): x=x%1.0; return min(x,1-x)
def Mf(S,Q):
    b=0.0
    for q in range(1,Q+1):
        for a in range(1,q):
            m=min(frac(s*a/q) for s in S)
            if m>b: b=m
    return b
def Me(S,Q):
    b=Fraction(0)
    for q in range(1,Q+1):
        for a in range(1,q):
            m=min(min((Fraction(s*a,q))%1,1-(Fraction(s*a,q))%1) for s in S)
            if m>b: b=m
    return b
def isprime(m): return m>1 and all(m%d for d in range(2,int(m**.5)+1))
def jacobsthal(m):  # max run of consecutive nonunits mod m
    nu=[1 if gcd(i,m)>1 else 0 for i in range(m)]
    best=run=0
    for x in nu+nu:
        run=run+1 if x else 0; best=max(best,run)
    return best
print("Single-swap tight k->g per n; is k composite? prime? Jacobsthal run vs n-k+2:")
for n in range(6,17):
    base=list(range(1,n)); hits=[]
    for k in base:
        for g in range(n,4*n+1):
            if g in base: continue
            S=[v for v in base if v!=k]+[g]
            if Mf(S,5*n)<=1.0/n+1e-9 and Me(S,10*n)==Fraction(1,n):
                hits.append((k,g))
    tagged=[(k,g,"prime" if isprime(k) else "composite", f"J={jacobsthal(k)} vs n-k+2={n-k+2}") for k,g in hits]
    print(f"  n={n}: {tagged if tagged else 'AP only'}")
print()
print("(2) are any LARGE PRIMES q in (n/2,n-1] ever the swapped-OUT k? (should be NONE):")
print("    census swapped-k above: check none is a prime > n/2.")
print()
print("(3) count of swappable k <= ? and (4) n=20 doubling 18->36:")
n=20; base=list(range(1,n)); S=[v for v in base if v!=18]+[36]
print(f"    n=20, swap 18->36: M={float(Mf(S,6*n)):.5f} vs 1/20={1/20:.5f}; tight={Me(S,10*n)==Fraction(1,20)}")
print(f"    (n=20: n-2=18 composite, 6|18, J(18)={jacobsthal(18)}>=n-k+2={20-18+2}=4? {jacobsthal(18)>=4}); predicts doubling since n=20==2 mod 6")
