"""
mac-mini-2026-07-08-S59 -- confirm the WEYL DECORRELATION for THM-527-A large-spread:
E_grid[W] -> (6/7)^k as spread grows, and the grid good-period sum stays > 0.

Fourier identity: E_grid[W] = (6/7)^k + sum_{n!=0: Vmax | n.e} What(n), the SAME resonance sum
as the density-floor tail decorrelation (THM-518 / LEM-009 machinery), now with 'Vmax | n.e'
in place of 'n.e = 0'. As the cluster spreads, resonances thin => E_grid[W] -> (6/7)^k > 0.
Good period exists <=> E_grid[W] > 6/(7Vmax).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
random.seed(591)

def W_at(E, j, Vmax):
    ph = sorted({(e*j) % Vmax for e in E}); m = len(ph)
    if m == 1: return F(6,7)
    W = F(0)
    for i in range(m):
        g = F(ph[i+1]-ph[i], Vmax) if i < m-1 else F(ph[0]+Vmax-ph[-1], Vmax)
        if g > F(1,7): W += g - F(1,7)
    return W
def egrid_W(E, Vmax):
    return sum((W_at(E,j,Vmax) for j in range(Vmax)), F(0)) / Vmax
def has_gap(E, j, Vmax):
    ph = sorted({(e*j) % Vmax for e in E}); m = len(ph)
    if m == 1: return True
    mg = max((ph[(i+1) % m]-ph[i]) % Vmax for i in range(m-1)); mg = max(mg, ph[0]+Vmax-ph[-1])
    return mg*7 > Vmax
def prim(E): E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1

print("DECORRELATION: E_grid[W] -> (6/7)^k as spread grows (random primitive, Vmax=spread+7)\n")
for k in (11, 13):
    iid = float(F(6,7)**k)
    print(f"k={k}: iid (6/7)^k = {iid:.5f}")
    for spread in [20, 40, 80, 160, 320, 640]:
        devs = []; egs = []; mingood = 10**9
        tries = 0
        while len(egs) < 12 and tries < 200:
            tries += 1
            mid = sorted(random.sample(range(1, spread), k-2)); E = [0]+mid+[spread]
            if len(set(E)) != k or not prim(E): continue
            Vmax = spread + 7
            eg = egrid_W(E, Vmax); egs.append(float(eg)); devs.append(float(eg)-iid)
            g = sum(1 for j in range(1,Vmax) if has_gap(E,j,Vmax))   # good-period count (excludes j=0)
            mingood = min(mingood, g)
        adev = sum(abs(d) for d in devs)/len(devs)
        print(f"   spread={spread:4d}: mean|E_grid[W]-iid|={adev:.5f}  meanE_grid[W]={sum(egs)/len(egs):.5f}  "
              f"min #good periods={mingood}")
    print()
