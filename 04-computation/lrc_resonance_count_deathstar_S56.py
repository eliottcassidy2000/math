from math import gcd
from fractions import Fraction as F
def totient(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def live_set(fam,n,q): return [p for p in range(1,q) if all(q<=n*((v*p)%q)<=(n-1)*q for v in fam)]
def M_exact(fam,Qmax):
    best=F(0)
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best
n=14
AP=list(range(1,14))

print("=== A. Divisor-covering is NECESSARY for tightness (upper-bound engine) ===")
# AP with 7 removed, replaced by 15 (still 13 elts, NO multiple of 7)
noseven=[1,2,3,4,5,6,8,9,10,11,12,13,15]
print(f"  family with no multiple of 7: min mult7? {any(v%7==0 for v in noseven)}")
print(f"  M = {M_exact(noseven,40)} (>1/14 => NOT tight, witnessed at t=1/7): {M_exact(noseven,40)>F(1,14)}")

print("\n=== B. Base resonance q=n: liveCount(14)=phi(14)=6 for tight families & coprime dilates ===")
for name,fam in [("AP",AP),("GW",[1,2,3,4,5,6,7,8,9,10,11,13,24]),("3*AP",[3*v for v in AP]),("5*AP",[5*v for v in AP])]:
    ls=live_set(fam,n,14); units=[o for o in range(1,14) if gcd(o,14)==1]
    print(f"  {name:5s}: liveCount(14)={len(ls)}  live={ls} == units {units}? {ls==units}")

print("\n=== C. The dilation ghost: 3*AP live@q=42 and reduced denominators ===")
fam=[3*v for v in AP]; ls=live_set(fam,n,42)
print(f"  3*AP liveCount(42)={len(ls)}; witnesses & reduced denom of p/42:")
denoms={}
for p in ls:
    d=42//gcd(p,42); denoms[d]=denoms.get(d,0)+1
print(f"  reduced-denominator histogram of the 18 live times: {denoms}  (12 at denom 42 = the ghosts, 6 at denom 14)")
print(f"  primitive? gcd(3*AP)={__import__('functools').reduce(gcd,fam)}")

print("\n=== D. Primitive tight families (AP,GW): liveCount(nm)=6 at all m (l(nm'')=0, m''>=2) ===")
for name,fam in [("AP",AP),("GW",[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
    ok=all(len(live_set(fam,n,14*m))==6 for m in range(1,31))
    print(f"  {name}: liveCount(14m)=6 for all m<=30? {ok}")
