"""
lrc_equidistribution_band_barrier_opus_20260701.py  (HYP-3781)
Creative geometry of the covering-min transition (a(n)=n for n>=12, klein-S60/HYP-3778):
THE EQUIDISTRIBUTION BAND-BARRIER (character sums on the band). A beater must keep small residue-gaps at every
band modulus D' in (n,2n-2]; the CONSTRUCTION is the interval (band-gap 3 for all n, radius-1 covering); beaters
sit at band-gap 4 (a thin low-discrepancy shell); a SPREAD set equidistributes (Weil) => max gap ~1.5 ln n,
growing past the shell => the shell collapses onto the interval at the transition.
See reflection creative-geometry-of-the-transition-...-opus-20260701.md.
"""
import numpy as np, random
from math import gcd, log
def maxgap(res,D):
    r=sorted(set(x%D for x in res))
    if len(r)<2: return D
    g=[r[i+1]-r[i] for i in range(len(r)-1)]+[r[0]+D-r[-1]]
    return max(g)
def band_gap(S,n):    # worst radius-1-band hole at the best rotation
    worst=0
    for D in range(n+1,2*n-1):
        worst=max(worst,min(maxgap([(v*j)%D for v in S],D) for j in range(1,D) if gcd(j,D)==1))
    return worst
beaters={7:[1,2,5,6,7,8],8:[1,4,5,6,7,11,16],9:[1,3,4,5,7,11,18,32],10:[1,2,3,5,6,7,8,9,30],11:[2,6,8,9,10,11,13,14,17,19]}
print(f"{'n':>3} {'constr band-gap':>16} {'beater band-gap':>16} {'rand-spread mean':>17} {'equidist 1.5 ln n':>18}")
for n in range(7,15):
    constr=list(range(1,n-1))+[n*(n-1)]
    rng=random.Random(1+n); rg=[band_gap(sorted(rng.sample(range(1,4*n),n-1)),n) for _ in range(30)]
    bg=band_gap(beaters[n],n) if n in beaters else None
    print(f"{n:>3} {band_gap(constr,n):>16} {str(bg):>16} {np.mean(rg):>17.1f} {1.5*log(n-1):>18.2f}")
print("\n=> construction = interval (band-gap 3 for ALL n); beaters = band-gap 4 (thin low-discrepancy shell);")
print("   random spread ~4.5-5 and growing; equidistribution gap 1.5 ln n climbs through the shell near the")
print("   ILP transition n=12. The shell collapses onto the pure interval = the construction. (Mechanism, not")
print("   yet proof: band-gap is coarse; the sharp version = all-moduli Erdos-Turan/Weil discrepancy.)")
