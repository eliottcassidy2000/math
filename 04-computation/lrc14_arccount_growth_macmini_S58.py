"""
mac-mini-2026-07-08-S58 -- growth law of #arcs of the good set vs cluster spread.
Finite-Vmax glue: a good PERIOD exists once #arcs < Vmax*rho* (rho* >= m_P). #arcs is
Vmax-independent. If #arcs = o(spread) (sublinear), then since spread <= Vmax, #arcs << Vmax
for all large Vmax => good period exists (small finite check on the rest). DECISIVE.
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(5271)

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

print("GROWTH LAW: #arcs vs spread (m_P=0.056487; good period exists if #arcs < Vmax*rho*)\n")
for k in (11, 13):
    print(f"k={k}:  {'spread':>8} {'#arcs min/med/MAX':>22} {'rho* min':>10} {'MAX#arcs/spread':>16}")
    for spread in [50, 100, 300, 1000, 3000, 10000, 30000]:
        GRID = min(30_000_000, max(3_000_000, 120*spread))
        arcs = []; rhos = []; tries = 0
        while len(arcs) < 25 and tries < 300:
            tries += 1
            mid = sorted(random.sample(range(1, spread), k-2))
            E = [0] + mid + [spread]
            if len(set(E)) != k or not primitive(E): continue
            n, meas = arc_count(E, GRID)
            arcs.append(n); rhos.append(meas)
        arcs.sort()
        amin, amed, amax = arcs[0], arcs[len(arcs)//2], arcs[-1]
        print(f"   {spread:>8} {f'{amin} / {amed} / {amax}':>22} {min(rhos):>10.4f} {amax/spread:>16.5f}")
    print()

print("INTERPRETATION: if MAX#arcs/spread -> 0, then #arcs = o(spread) <= o(Vmax), so")
print("#arcs < Vmax*m_P for all Vmax > (small); the large-spread finite-Vmax case CLOSES.")
