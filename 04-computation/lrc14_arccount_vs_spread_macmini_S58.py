"""
mac-mini-2026-07-08-S58 -- THE decisive computation for THM-527-A large-spread half.

Finite-Vmax glue: #{grid j/Vmax in good set} >= Vmax*rho* - #arcs (discrepancy/variation
bound), rho* >= m_P (density floor). So a good period EXISTS as soon as  #arcs < Vmax*m_P.
#arcs = #arcs of the good set G* = {y: maxgap{frac(e_i y)} > 1/7}, and by the bounded-arc-count
lemma it is Vmax-INDEPENDENT (depends only on cluster-internal differences).

QUESTION: does #arcs grow with the cluster SPREAD, or stay ~ k?
 - If #arcs = O(k) uniformly  => the condition Vmax > #arcs/m_P ~ O(k)/m_P ~ 230 is a SMALL
   finite check, and the large-spread case DISSOLVES (whole finite-Vmax glue closes).
 - If #arcs ~ spread => large-spread needs the high-rho* / decorrelation argument.
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(527)

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

print("#arcs of {y: maxgap{frac(e_i y)}>1/7} vs cluster SPREAD (grid ~ 60*spread, resolves m/(diff))\n")
for k in (11, 13):
    print(f"k={k}:  (m_P=0.056487; good period exists if #arcs < Vmax*m_P)")
    print(f"  {'spread':>7} {'#arcs (min/median/max over 40 random k-subsets)':>50} {'rho* min':>10}")
    for spread in [10, 20, 40, 80, 160, 320]:
        GRID = max(2_000_000, 200*spread)
        arcs = []; rhos = []
        tries = 0
        while len(arcs) < 40 and tries < 400:
            tries += 1
            mid = sorted(random.sample(range(1, spread), k-2))
            E = [0] + mid + [spread]
            if len(set(E)) != k or not primitive(E): continue
            n, meas = arc_count(E, GRID)
            arcs.append(n); rhos.append(meas)
        arcs.sort()
        amin, amed, amax = arcs[0], arcs[len(arcs)//2], arcs[-1]
        print(f"  {spread:>7} {f'{amin} / {amed} / {amax}':>50} {min(rhos):>10.4f}")
    print()

# also the STRUCTURED worst cases: long-AP + outlier (max additive energy) at various spreads
print("structured (block+far-outlier) #arcs vs spread:")
for k in (11, 13):
    print(f"  k={k}:")
    for D in [k, 2*k, 5*k, 20*k, 100*k]:
        E = list(range(k-1)) + [D]
        GRID = max(2_000_000, 200*D)
        n, meas = arc_count(E, GRID)
        print(f"    block+outlier spread={D:5d}: #arcs={n:3d}  rho*={meas:.4f}")
