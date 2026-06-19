import sys
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def covered(S, r):
    """True iff danger arcs {t: ||v t|| <= r} over all v in S cover [0,1).
       i.e. M(S) <= r  (no lonely point with all ||vt|| > r). r a Fraction."""
    # arcs around m/v of half-width r/v in t: [ (m-r)/v, (m+r)/v ], for m=0..v (cover [0,1])
    arcs=[]
    for v in S:
        # m/v +- r/v ; m ranges so that interval meets [0,1]
        # left endpoints (m-r)/v <=1  and (m+r)/v>=0
        mlo = 0
        mhi = v
        for m in range(mlo, mhi+1):
            lo = Fraction(m, v) - Fraction(r, v)
            hi = Fraction(m, v) + Fraction(r, v)
            if hi < 0 or lo > 1: continue
            arcs.append((max(lo,Fraction(0)), min(hi,Fraction(1))))
    arcs.sort()
    # check cover of [0,1]
    cur = Fraction(0)
    for lo,hi in arcs:
        if lo > cur:   # gap
            return False
        if hi > cur: cur = hi
        if cur >= 1: return True
    return cur >= 1

def witness_ge(S, r):
    """True iff exists t with min_v ||vt|| >= r (so M(S) >= r). Use covering complement:
       M(S) >= r iff NOT covered(S, r') for r' just below r. Equivalent: lonely set at
       radius slightly less than r nonempty. We test: not covered(S, r - tiny). Since values
       are k-rational, test covered at r and detect equality via the boundary."""
    # M(S) = inf{ r : covered(S,r) }. So M(S)==v* iff covered(S,v*) and not covered(S, v*-eps).
    # eps: next smaller value. Use a strictly smaller rational: v* - 1/(BIG).
    return None

def M_equals(S, target):
    """Return True iff M(S) == target (Fraction)."""
    S=tuple(sorted(set(S)))
    # covered at target?
    if not covered(S, target):
        return False, ">target"   # M > target
    # not covered just below target?  use target - 1/(2*L) where L = lcm of denominators bound
    # safe eps smaller than any gap between distinct max-min values: use 1/(2*(maxS)^2 *(k+2))
    denom_bound = (max(S)**2) * (len(S)+2) * 4
    below = target - Fraction(1, denom_bound)
    if below <= 0: below = target/2
    if covered(S, below):
        return False, "<target"   # M < target (covered even below)
    return True, "=="

print("=== PRIMORIAL conjecture: level a first appears via F(k,a)={1..k-2,k,a(k-1)} at k-1=primorial ===")
prim=[]
# primorials: 2, 6, 30, 210, 2310
primorials={2:2, 3:6, 5:30, 7:210, 11:2310}
tests=[(3,7),(4,31),(5,211),(3,13),(3,19),(4,31),(6,2311)]
print("\nDirect level checks F(k,a) =?= a/(a(k+1)-1):")
for a,k in [(2,5),(3,7),(3,13),(3,19),(4,31),(5,211),(6,2311)]:
    S=sorted(set(list(range(1,k-1))+[k,a*(k-1)]))
    if len(S)!=k:
        print(f"  a={a} k={k}: size mismatch"); continue
    target=Fraction(a, a*(k+1)-1)
    eq,why=M_equals(S, target)
    floor=Fraction(1,k+1)
    gk2=float(target-floor)*k*k
    fac=sorted(set(p for p in [2,3,5,7,11,13] if (k-1)%p==0))
    print(f"  a={a:2d} k={k:4d} (k-1={k-1}={fac}+...): M==a/(a(k+1)-1)={target}? {eq} [{why}]  g*k^2={gk2:.4f} ~1/a={1/a:.4f}")
