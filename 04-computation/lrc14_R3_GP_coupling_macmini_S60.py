"""
mac-mini-2026-07-08-S60 -- R3: the G_P / cluster coupling in the finite-Vmax existence route.

A good PERIOD (=> M(S)>=1/14) is a grid j with BOTH:
  Good_E(j): maxgap{frac(e_i j/Vmax)} > 1/7   (cluster co-offsets),  AND
  G_P(j):    ||p j/Vmax|| >= 1/14  for all p in P   (observer block).
LEM-010 handles Good_E (cluster) -- but its good period can be small j, where small observer
speeds p give ||p/Vmax|| < 1/14, FAILING G_P. So the coupling is NOT automatic. This tests
whether G_P ∩ Good_E has a grid point for admissible (P,E), incl. the hard spread-dense cluster,
and reports WHICH j (small? spread-out?).
"""
from math import gcd
from functools import reduce
from fractions import Fraction as F
import random
random.seed(603)

def cluster_good(E, j, V):
    ph = sorted({(e*j) % V for e in E}); m = len(ph)
    if m == 1: return True
    mg = max((ph[(i+1) % m]-ph[i]) % V for i in range(m-1)); mg = max(mg, ph[0]+V-ph[-1])
    return mg*7 > V
def in_GP(P, j, V):
    # ||p j / V|| >= 1/14  <=>  min(r, V-r) >= V/14 where r = p*j mod V
    for p in P:
        r = (p*j) % V
        if 14*min(r, V-r) < V: return False
    return True
def prim(E):
    E = sorted(E); return len(E) >= 2 and reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1

print("R3: does a grid good period exist in G_P ∩ Good_E?  (admissible P,E; hard spread-dense E)\n")
for k in (8, 11, 13):
    psz = 13 - k
    if psz < 1: continue
    nfail = 0; tested = 0; jsmall = 0; examples = []
    for _ in range(3000):
        V = random.randint(200, 4000)
        # observer block P = small distinct speeds in [1,13]
        P = sorted(random.sample(range(1, 14), psz))
        # cluster co-offsets E: spread-dense (hard: spread>=6V/7) so cluster j=1 fails
        lo = 6*V//7 + 1
        if lo >= V: continue
        s = random.randint(lo, V-1)
        mid = sorted(random.sample(range(1, s), k-2)); E = [0]+mid+[s]
        if len(set(E)) != k or not prim(E): continue
        tested += 1
        # find a grid good period in the intersection
        found = None
        for j in range(1, V):
            if cluster_good(E, j, V) and in_GP(P, j, V):
                found = j; break
        if found is None: nfail += 1; examples.append((P, tuple(E[:5]), V))
        elif found <= k: jsmall += 1
    print(f"k={k} (|P|={psz}): {tested} admissible spread-dense (P,E); "
          f"NO good period in G_P∩Good_E: {nfail}; smallest-good-j <= k in {jsmall} cases")
    if examples: print(f"    FAILURE example: P={examples[0][0]} E~{examples[0][1]} V={examples[0][2]}")
print("\n=> if nfail=0, the coupling HOLDS (a grid good period in the intersection always exists);")
print("   the good j need not be small (small j fails G_P for small p) -- it's a spread-out j.")
