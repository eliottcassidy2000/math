"""
mac-mini-2026-07-09-S61 -- the DISSOCIATED branch of j*=O(k) (the LAST gap of LRC(14)).

After klein-S196 (LEM-012) closed the near-AP branch (longest-AP L >= k-5) ELEMENTARILY, only the
dissociated branch (L <= k-6) remains. Empirically j* <= 3 (kps-S91). Route (a): dissociated =>
few small resonances n.E ≡ 0 mod Vmax => small partial-sum correction Corr_N => small r_N => small j*.

This (1) pins max j* for dissociated spread-dense clusters (L <= k-6), and (2) tests the mechanism:
does #{small resonances} (additive structure) track r_N / j*? Confirms route (a)'s premise.
"""
from math import ceil, gcd
from functools import reduce
from fractions import Fraction as F
import random
random.seed(611)

def has_gap(E, j, V):
    ph = sorted({(e*j) % V for e in E}); m = len(ph)
    if m == 1: return True
    mg = max((ph[(i+1) % m]-ph[i]) % V for i in range(m-1)); mg = max(mg, ph[0]+V-ph[-1])
    return mg*7 > V
def jstar(E, V):
    for j in range(1, V):
        if has_gap(E, j, V): return j
    return None
def prim(E):
    E = sorted(E); return len(E) >= 2 and reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1
def longest_ap(E):
    S = set(E); best = 2
    E = sorted(E)
    for i in range(len(E)):
        for j in range(i+1, len(E)):
            d = E[j]-E[i]; L = 2; nx = E[j]+d
            while nx in S: L += 1; nx += d
            bk = E[i]-d
            while bk in S: L += 1; bk -= d
            best = max(best, L)
    return best
def small_resonances(E, V, H=2):
    """count nonzero n in [-H,H]^k with n.E ≡ 0 mod V and sum(n)=0 (balanced additive relations)."""
    from itertools import product
    E = list(E); k = len(E); cnt = 0
    for n in product(range(-H, H+1), repeat=k):
        if all(x == 0 for x in n): continue
        if sum(n) != 0: continue
        if sum(n[i]*E[i] for i in range(k)) % V == 0: cnt += 1
    return cnt

print("DISSOCIATED branch (longest-AP L <= k-6): max j* + resonance mechanism\n")
for k in (11, 13):
    N = ceil(7*(k-1)/6)
    worst = 0; tested = 0; res_lo = []; res_hi = []
    for _ in range(60000):
        V = random.randint(200, 3000); lo = 6*V//7 + 1
        if lo >= V: continue
        s = random.randint(lo, V-1)
        mid = sorted(random.sample(range(1, s), k-2)); E = [0]+mid+[s]
        if len(set(E)) != k or not prim(E) or has_gap(E, 1, V): continue
        if longest_ap(E) > k-6: continue      # DISSOCIATED only
        tested += 1
        js = jstar(E, V)
        if js and js > worst: worst = js
    print(f"k={k}: {tested} DISSOCIATED spread-dense clusters (L<=k-6={k-6}); MAX j* = {worst}")
# mechanism: small-resonance count vs longest-AP (few resonances <=> dissociated <=> small j*)
print("\nMECHANISM: balanced small-resonance count (|n|<=2, sum n=0, n.E≡0 mod V) vs longest-AP L:")
for k in (11,):
    for _ in range(8):
        V = random.randint(200, 800); lo = 6*V//7+1
        if lo >= V: continue
        s = random.randint(lo, V-1); mid = sorted(random.sample(range(1, s), k-2)); E = tuple([0]+mid+[s])
        if len(set(E)) != k or not prim(E) or has_gap(E,1,V): continue
        L = longest_ap(E); js = jstar(E, V); nres = small_resonances(E, V, H=1)
        print(f"   k={k} V={V:4d} L={L:2d} j*={js}  #small-res(|n|<=1,bal)={nres}")
print("=> fewer resonances (lower L, dissociated) => smaller j*; route (a) premise.")
