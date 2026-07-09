"""
mac-mini-2026-07-08-S58 -- THM-527-A finite-Vmax glue: the SUFFICIENT CONDITION and its
stress test on ADVERSARIAL / structured clusters.

REDUCTION (this session): a good PERIOD exists (=> gap>1/7 witness => M(S)>=1/14) as soon as
  #{good grid pts} >= Vmax*rho* - #arcs(G*) > 0,   i.e.   #arcs(G*) < Vmax * rho*(E).
Here G* = {y: maxgap{frac(e_i y)}>1/7}, rho* = meas(G*) >= m_P (density floor, THM-661),
#arcs is Vmax-INDEPENDENT (bounded-arc-count lemma). Since spread <= Vmax, the BINDING case is
Vmax = spread, giving the clean test:  rho*(E) > #arcs(E)/spread(E)  (large spread), OR
#arcs(E) < Vmax*rho* with #arcs an absolute constant (small spread, finite check).

Two-regime closure:
 (small spread <= S0): #arcs <= A0 absolute => good period for Vmax > A0/m_P (finite check).
 (large spread > S0):  #arcs <= beta*spread with beta < rho*_min => Vmax(rho*-beta)>0 always.
This tests beta = max #arcs/spread and rho*_min on the WORST (structured) shapes.
"""
import numpy as np
from math import gcd
from functools import reduce

def arc_count(E, GRID):
    x = (np.arange(GRID)+0.5)/GRID
    Ea = np.array(sorted(E), float)
    ph = np.mod(np.outer(x, Ea), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0]+1-ph[:, -1])[:, None]], axis=1)
    good = (g.max(axis=1) > 1/7).astype(int)
    tr = int(np.sum((good - np.roll(good, 1)) == 1))
    return tr, good.mean()
def primitive(E):
    E = sorted(E); return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1

print("STRESS TEST: beta=#arcs/spread and rho* on structured clusters (worst case Vmax=spread)")
print("PASS iff rho* > #arcs/spread (large spread) -- i.e. the binding Vmax=spread has a good period.\n")

def shapes(k, s):
    """structured high-additive-energy / high-arc-count k-clusters of spread s."""
    out = {}
    out['block+outlier'] = list(range(k-1)) + [s]
    out['2-block'] = list(range(k//2)) + [s-(k-k//2)+1+i for i in range(k-k//2)]
    if s >= 2*(k-1):
        d = s//(k-1); out[f'AP(d={d})'] = [d*i for i in range(k-1)] + [s]  # near-AP, max energy
    step = max(1, s//(k+1))
    perf = sorted(set(list(range(0, s+1, step))[:k-1]) | {0, s})
    out['perforated-block'] = perf
    return {n: E for n, E in out.items() if len(set(E)) == k and primitive(E)}

worst_beta = {11: 0, 13: 0}
for k in (11, 13):
    print(f"k={k}:")
    for s in [40, 100, 400, 1600, 6400]:
        GRID = min(40_000_000, max(4_000_000, 120*s))
        for name, E in shapes(k, s).items():
            n, meas = arc_count(E, GRID)
            beta = n/s
            worst_beta[k] = max(worst_beta[k], beta if s >= 100 else 0)
            ok = meas > n/s          # rho* > #arcs/spread  (Vmax=spread binding)
            print(f"   s={s:5d} {name:16s}: #arcs={n:5d} beta={beta:.4f} rho*={meas:.4f}  "
                  f"{'PASS' if ok or n < 40 else 'CHECK'}  (rho*-beta={meas-beta:+.4f})")
    print(f"   => worst beta (s>=100) = {worst_beta[k]:.4f}  vs rho*_min(large s) ~ 0.98  "
          f"=> margin {0.98-worst_beta[k]:+.3f}\n")
print("If worst beta << rho*_min for all structured shapes, the large-spread half CLOSES:")
print("#arcs <= beta*spread <= beta*Vmax < rho**Vmax => #{good grid} > 0 => good period exists.")
