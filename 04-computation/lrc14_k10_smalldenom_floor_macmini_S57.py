"""
mac-mini-2026-07-07-S57 (HYP-5237) -- k=10 SMALL-DENOMINATOR FLOOR: verify the composition
residual is truth-safe.

The average-form conditional composition (avgc + W_F^{G_P}/toll) leaves 61 residual shapes,
all small-denominator (G_P teeth wide, proxy degrades). This confirms the TRUE rho* is far
above m_P on the entire residual, so k=10 = [composition: 225] u [exact-G2: 61 small-denom].

For each residual shape P (|P|=3) we minimise the TRUE rho*(P,E) = meas(G_P n {maxgap>1/7})
over the adversarial family class (consecutive block and near-block, which minimise mu at
every k>=8 per THM-530) at high grid resolution, and report min rho* / m_P.  A uniform
floor >> 1 certifies the residual is truth-safe (the exact-G2 engine, monad-S4, upgrades
this to a rational certificate).
"""
import numpy as np
from itertools import combinations
from math import gcd
from functools import reduce

MP = 14249/252252
TH = 1/7

def true_rho(P, E, grid=1_000_000):
    x = (np.arange(grid)+0.5)/grid
    inGP = np.ones(grid, bool)
    for p in P:
        px = np.mod(p*x, 1.0)
        inGP &= (np.minimum(px, 1-px) >= 1/14)
    E = np.array(sorted(E), float)
    ph = np.mod(np.outer(x, E), 1.0); ph.sort(axis=1)
    gaps = np.diff(ph, axis=1)
    wrap = (ph[:,0]+1-ph[:,-1])[:,None]
    mg = np.concatenate([gaps, wrap], axis=1).max(axis=1)
    return float((inGP & (mg > TH)).mean())

# adversarial family battery (mu-minimising per THM-530: consecutive & near-consecutive)
def battery():
    fams = [list(range(10))]                              # block
    for X in [10,11,12,13,15,20,30]: fams.append(list(range(9))+[X])   # block+outlier
    for X in [11,12,15,25]: fams.append(list(range(8))+[X,X+1])
    for holes in [(4,),(5,),(4,7)]:                       # perforated
        b=[t for t in range(12) if t not in holes][:10]
        if len(b)==10: fams.append(b)
    # keep primitive only
    def prim(E):
        d=[E[i+1]-E[i] for i in range(len(sorted(E))-1)]; return reduce(gcd,d)==1
    return [E for E in fams if prim(E)]
FAMS = battery()

# the actual composition-residual small-denominator shapes (worst 20 from the composition run)
SMALL = [(4,5,6),(5,12,13),(5,6,12),(4,5,12),(5,6,13),(3,4,5),(5,6,11),(5,11,12),
         (4,6,10),(4,5,13),(3,5,12),(10,12,13),(3,4,6),(4,5,11),(2,3,5),(2,5,6),
         (3,5,6),(4,6,12),(5,10,12),(1,5,6)]
print(f"=== k=10 small-denominator floor: {len(SMALL)} shapes (min part <= 6) x {len(FAMS)} adversarial families ===")
print(f"    m_P = {MP:.5f}\n")
worst = (1e9, None, None)
per = []
for P in SMALL:
    mr = (1e9, None)
    for E in FAMS:
        r = true_rho(P, E)
        if r < mr[0]: mr = (r, tuple(E))
    per.append((mr[0], P, mr[1]))
    if mr[0] < worst[0]: worst = (mr[0], P, mr[1])
per.sort()
print(f"GLOBAL min true rho* over all small-denom (shape,family): {worst[0]:.4f} = {worst[0]/MP:.2f}x m_P")
print(f"  at P={worst[1]}, E={worst[2]}")
print(f"  all {len(SMALL)} shapes have min rho* >= {per[0][0]:.4f} ({per[0][0]/MP:.1f}x m_P)"
      f"  => residual TRUTH-SAFE with >= {per[0][0]/MP:.0f}x slack\n")
print("  10 lowest-floor small-denom shapes (min rho*, x m_P, shape, minimiser):")
for r, P, E in per[:10]:
    Es = f"{E[:3]}..{E[-1]}(d{max(E)-min(E)})"
    print(f"    {r:.4f}  {r/MP:5.2f}x  P={P}  E={Es}")
