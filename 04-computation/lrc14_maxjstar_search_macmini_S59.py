"""
mac-mini-2026-07-08-S59 -- HARD adversarial search for the max smallest-good-period j*.
If j* is uniformly bounded (j* <= f(k), small), THM-527-A closes: a good period always exists
at tiny j, so no finite check beyond Vmax <= f(k) (already done exact for Vmax<=1001, kps-S30).
Targets the spread-dense / j=1-fails regime (spread >= 6Vmax/7, all internal gaps <= Vmax/7).
"""
from math import gcd
from functools import reduce
import random
random.seed(5927)

def has_gap(E, j, V):
    ph = sorted({(e*j) % V for e in E}); m = len(ph)
    if m == 1: return True
    mg = max((ph[(i+1) % m]-ph[i]) % V for i in range(m-1)); mg = max(mg, ph[0]+V-ph[-1])
    return mg*7 > V
def jstar(E, V):
    for j in range(1, V):
        if has_gap(E, j, V): return j
    return None
def prim(E): E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1

print("max j* over adversarial spread-dense clusters (j=1 fails), by k:\n")
for k in (8, 11, 13):
    worst = 0; worstE = None; worstV = None; nfail = 0; tested = 0
    # 1) AP-like dense clusters at many scales
    for V in range(100, 30000, 3):
        for _ in range(3):
            # dense cluster: near-AP with small jitter, spread >= 6V/7
            d = random.randint(V//14, V//7)
            base = [ (d*i + random.randint(0,2)) % V for i in range(k)]
            E = sorted(set(base))
            if len(E) != k: continue
            E = [e - min(E) for e in E]
            if max(E) < 6*V//7 or max(E) >= V or not prim(E): continue
            if has_gap(E, 1, V): continue    # want j=1 to fail
            tested += 1
            js = jstar(E, V)
            if js is None: nfail += 1
            elif js > worst: worst, worstE, worstV = js, tuple(E), V
    # 2) fully random dense clusters
    for _ in range(30000):
        V = random.randint(200, 20000); lo = 6*V//7 + 1
        if lo >= V: continue
        s = random.randint(lo, V-1)
        mid = sorted(random.sample(range(1, s), k-2)); E = [0]+mid+[s]
        if len(set(E)) != k or not prim(E) or has_gap(E,1,V): continue
        tested += 1; js = jstar(E, V)
        if js is None: nfail += 1
        elif js > worst: worst, worstE, worstV = js, tuple(E), V
    print(f"k={k}: {tested} adversarial spread-dense clusters; MAX j* = {worst}  "
          f"(no-good-j: {nfail}); worst at Vmax={worstV}")
    if worstE: print(f"       worst-j* cluster (first 6 elts): {worstE[:6]}... Vmax={worstV}")
print("\n=> if max j* stays small (<< Vmax) for all k, a good period exists at tiny j always")
print("   => THM-527-A closes: finite check only Vmax <= max j* (trivially done, kps-S30 Vmax<=1001).")
