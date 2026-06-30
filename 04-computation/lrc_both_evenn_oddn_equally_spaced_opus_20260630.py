"""
EVEN-n: even block 2*{1..n-1} at t=1/(2n) gives the EQUALLY-SPACED AP {1/n..(n-1)/n}, M=1/n (verify).
ODD-n: covering-min (find extremal set + structure). The apex gap = Cayley spectral gap (Ramanujan thread).
"""
import math
from fractions import Fraction
from itertools import combinations
def M_exact(S,Qmax):
    best=Fraction(0); bq=None
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
        if v>best:best=v;bq=q
    return best,bq
def is_cov(S,n): return len(set(S))==n-1 and all(any(x%q==0 for x in S) for q in range(2,n+1))
print("EVEN-n: even block 2*{1..n-1}, witness, equally-spaced AP:")
for n in [8,14]:
    eb=[2*k for k in range(1,n)]
    M,q=M_exact(eb,4*n)
    print(f"  n={n}: even block M={M}={float(M):.5f}=1/n? {M==Fraction(1,n)}; witness q={q} (=2n? {q==2*n}); at t=1/(2n) speeds*1/(2n)=k/n equally-spaced")
print()
print("ODD-n covering-min (thorough search, find extremal):")
for n,B,Q in [(7,18,45),(9,15,82)]:
    Phi6=n*n-n+1; constr=Fraction(n,Phi6)
    best=Fraction(1); bestS=None
    for combo in combinations(range(1,B+1),n-1):
        S=list(combo)
        if not is_cov(S,n): continue
        M,_=M_exact(S,Q)
        if M<best: best=M; bestS=S
    print(f"  n={n}(ODD): covering-min={best}={float(best):.5f} at {bestS}  (construction {constr}={float(constr):.5f}, 1/n={1/n:.5f})")
    # structure of extremal
    odds=[x for x in bestS if x%2==1]; evens=[x for x in bestS if x%2==0]
    print(f"     extremal: odds={odds}, evens={evens}; odd count={len(odds)} (must cover odd q=n via an odd mult)")
