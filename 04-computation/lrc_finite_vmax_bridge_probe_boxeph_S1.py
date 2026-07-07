"""
FINITE-Vmax BRIDGE probe (THM-527.G.3, the 'o(Vmax) arc bound' monad-explorer
HYP-4787 flags as Part-A's remaining need).  boxeph-2026-07-07-S1, HYP-4760.

The density->reach reduction is EXACT in the Vmax->inf limit: good-period
fraction rho_K -> rho* = mu_{1/7}(co-offset config).  For a FIXED instance we
need rho_K > 0 (>=1 good period => M(S)>=1/14).  Question: how large must Vmax be
(V_0) for rho_K > 0, i.e. how fast does rho_K -> the 1/7-density floor?

Since mu_{1/7} is DILATION-INVARIANT, any compact binding family (e.g. the moat
{1..13}) can be dilated to make Vmax huge with mu_{1/7}, M both unchanged.  So the
real question: for the co-offset shape E (0 in E), rho_K(Vmax) vs the floor
mu_{1/7}(E).  If rho_K tracks the floor already at small Vmax, the bridge is
easy (the 0.44 floor dominates the O(#arcs/Vmax) error); V_0 is small.
"""
from fractions import Fraction as F
import math

def dist_int(y):
    f = y - math.floor(y)
    return min(f, 1.0 - f)

def actual_good_fraction(S, Vmax, npts=500):
    """EXACT-ish fraction of Vmax-ruler periods I_j that contain a fully-safe tau."""
    good = 0
    for j in range(Vmax):
        lo = (14*j+1)/(14*Vmax); hi = (14*j+13)/(14*Vmax)
        mx = 0.0
        for a in range(npts):
            t = lo + (hi-lo)*(a+0.5)/npts
            m = 1.0
            for v in S:
                d = dist_int(v*t)
                if d < m:
                    m = d
                    if m < mx: break
            if m > mx: mx = m
        if mx >= 1/14 - 1e-9:
            good += 1
    return good / Vmax

# co-offset shapes E (0 in E). floor = mu_{1/7}(E) computed exactly elsewhere.
shapes = {
    "APcoff {0..12}": (list(range(13)), 0.4425),   # binding: floor 477/1078
    "perf7 {0,2..6,8}": ([0,2,3,4,5,6,8], 1.0),
    "spread {0,2..12,17,28}": ([0,2,3,4,5,6,7,8,9,10,11,12,17,28][:13], None),
}
print("Finite-Vmax bridge: good-period fraction rho_K(Vmax) vs the 1/7-density floor")
print("(instance S = {Vmax} u {Vmax - e : e in E, e>0})")
print("="*70)
for name,(E, floor) in shapes.items():
    E = sorted(set(E))
    print(f"\n{name}  (floor mu_1/7 ~ {floor})")
    print(f"  {'Vmax':>6} | {'#good':>6} | {'rho_K':>7} | {'M(S)>=1/14?':>11}")
    V0 = None
    for Vmax in [14, 15, 20, 30, 50, 100, 300, 1000]:
        S = [Vmax] + [Vmax - e for e in E if e > 0]
        # skip if duplicate speeds
        if len(set(S)) != len([e for e in E if e>0]) + 1:
            continue
        rho = actual_good_fraction(S, Vmax, npts=400 if Vmax<=100 else 200)
        ok = "YES" if rho > 0 else "no"
        if rho > 0 and V0 is None: V0 = Vmax
        print(f"  {Vmax:>6} | {int(round(rho*Vmax)):>6} | {rho:>7.4f} | {ok:>11}")
    print(f"  => smallest Vmax with a good period (V_0) <= {V0}")

print("\n" + "="*70)
print("INTERPRETATION")
print("="*70)
print("""If rho_K(Vmax) tracks the 1/7-density floor already for small Vmax (V_0
small), then the finite-Vmax bridge is EASY given the floor: mu_{1/7}(E) >= 0.32
(union-bound ledger, monad-explorer HYP-4787) dominates the O(#arcs/Vmax) finite
error for Vmax >= V_0, and Vmax < V_0 is a tiny finite check.  This makes the
mu_{1/7} TAIL route's remaining Part-A gap benign -- in contrast to the
reverse-Markov E[maxgap] MEAN route whose T*=0.191 bar leaves only ~0.006 margin.
Conclusion: the comfortable object is mu_{1/7} (the sharp-1/7 good-density),
NOT the mean E[maxgap].""")
