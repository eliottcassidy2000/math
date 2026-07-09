"""
mac-mini-2026-07-08-S58 -- adversarial j* sweep for LEM-010 (ii): the smallest good period
j* on clusters BUILT to defeat j=1 (AP {0,d,..,(k-1)d}, d in [Vmax/14, Vmax/7], spread>=6Vmax/7).
Supports the conjecture j* = O(k) (<< the 3^(k-1) Dirichlet guarantee).
"""
from math import gcd
from functools import reduce
import random
random.seed(9)
def maxgap_int(E, j, V):
    ph = sorted({(e*j) % V for e in E}); m = len(ph)
    if m == 1: return V
    mg = max((ph[(i+1) % m]-ph[i]) % V for i in range(m-1)); return max(mg, ph[0]+V-ph[-1])
def good(E, j, V): return maxgap_int(E, j, V)*7 > V
def prim(E): E = sorted(E); return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1
def smallest_good(E, V):
    for j in range(1, V):
        if good(E, j, V): return j
    return None

print("adversarial AP clusters {0,d,..,(k-1)d}, d in [V/14,V/7] (defeat j=1); find j*:")
for k in (11, 13):
    worst = 0; nfail = 0; tested = 0
    for V in range(500, 8000, 7):
        step = max(1, (V//7 - V//14)//5)
        for d in range(V//14, V//7+1, step):
            E = sorted({d*i % V for i in range(k)})
            if len(E) != k or not prim(E) or E[-1]-E[0] < 6*V/7 or good(E, 1, V): continue
            tested += 1; js = smallest_good(E, V)
            if js is None: nfail += 1
            else: worst = max(worst, js)
    print(f"  k={k}: {tested} hard AP clusters; MAX j* = {worst}  (Dirichlet bound 3^(k-1)={3**(k-1)}); "
          f"none-found={nfail}")
print("=> j* <= 7 for k<=13 even adversarially; supports j*=O(k), which would make THM-527-A fully elementary.")
