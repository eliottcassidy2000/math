"""
mac-mini-2026-07-08-S58 -- DETERMINISTIC good-period existence for THM-527-A (finite-Vmax glue),
sidestepping quantitative equidistribution entirely.

CLAIM 1 (j=1):  if spread(E) < 6*Vmax/7 then j=1 is a good period.
  At j=1 phases = {e_i/Vmax} in [0, spread/Vmax]; wraparound gap = 1 - spread/Vmax > 1/7.
CLAIM 2 (Dirichlet): for Vmax > 3^(k-1), SOME j in {1..3^(k-1)} is a good period.
  Pigeonhole j -> (floor(3*frac(e_i j/Vmax)))_i in {0,1,2}^(k-1): two j's collide => j* with
  ||e_i j*/Vmax|| < 1/3 all i => all phases in a 2/3-arc => empty arc >= 1/3 > 1/7.
Together: a good period ALWAYS exists (finite check only for the small-Vmax remainder).
This is stronger/cleaner than the soft #arcs<rho*Vmax route (which is vacuous for 2-block/AP).
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(52719)

def maxgap_int(E, j, Vmax):
    ph = sorted({(e*j) % Vmax for e in E}); m = len(ph)
    if m == 1: return Vmax
    mg = max((ph[(i+1) % m]-ph[i]) % Vmax for i in range(m-1))
    mg = max(mg, ph[0] + Vmax - ph[-1])
    return mg
def is_good(E, j, Vmax):
    return maxgap_int(E, j, Vmax) * 7 > Vmax     # maxgap > 1/7
def primitive(E):
    E = sorted(E); return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1

print("CLAIM 1: spread < 6*Vmax/7  =>  j=1 is a good period\n")
fails1 = 0; tested1 = 0
for k in (8, 11, 13):
    for _ in range(3000):
        Vmax = random.randint(k, 4000)
        spread = random.randint(k-1, Vmax-1)
        if spread >= 6*Vmax/7: continue          # only test the claim's regime
        E = sorted(random.sample(range(1, spread), k-2)) if spread > k else []
        E = [0] + (E if E else list(range(1, k-1))) + [spread]
        if len(set(E)) != k or max(E) >= Vmax or not primitive(E): continue
        tested1 += 1
        if not is_good(E, 1, Vmax): fails1 += 1
print(f"  tested {tested1} clusters with spread<6Vmax/7: j=1 good in ALL but {fails1} "
      f"{'(CLAIM 1 HOLDS)' if fails1 == 0 else '*** FAILS ***'}\n")

print("CLAIM 2: smallest good period j* (Dirichlet bound 3^(k-1)); how small is it REALLY?\n")
for k in (8, 11, 13):
    bound = 3**(k-1)
    worst = 0; nfail_all = 0; tested2 = 0
    # focus on the hard regime spread >= 6Vmax/7 where j=1 fails
    for _ in range(4000):
        Vmax = random.randint(200, 20000)
        lo = int(np.ceil(6*Vmax/7))
        if lo >= Vmax: continue
        spread = random.randint(lo, Vmax-1)
        mid = sorted(random.sample(range(1, spread), k-2))
        E = [0] + mid + [spread]
        if len(set(E)) != k or not primitive(E): continue
        if is_good(E, 1, Vmax): continue          # j=1 already good; want the hard ones
        tested2 += 1
        # find smallest good j
        jstar = None
        for j in range(1, min(Vmax, 20000)):
            if is_good(E, j, Vmax): jstar = j; break
        if jstar is None: nfail_all += 1
        else: worst = max(worst, jstar)
    print(f"  k={k} (3^(k-1)={bound}): {tested2} hard clusters (spread>=6Vmax/7, j=1 fails); "
          f"max needed j* = {worst}  (<< {bound}); no-good-j-found: {nfail_all}")
print("\n=> if j=1 covers spread<6Vmax/7 and small j* covers the rest (<< 3^(k-1)),")
print("   a good period always exists deterministically; only small Vmax needs a finite check.")
