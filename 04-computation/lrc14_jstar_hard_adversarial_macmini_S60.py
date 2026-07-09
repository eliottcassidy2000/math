"""
mac-mini-2026-07-08-S60 -- HARD adversarial search for large j* (smallest good period).
Goal: decide whether j* = O(k) is robust or breakable. Constructs clusters DESIGNED to stay
"bad" (phases 1/7-dense) for many consecutive j: perturbed APs, interleaved APs, geometric,
mod-small-prime structured, and Vmax with many factors. If max j* stays ~ k, the conjecture is
robust; if it blows up, j*=O(k) is false and R1 needs a different route.
"""
from math import gcd
from functools import reduce
import random
random.seed(60)

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
def valid(E, V):
    return len(set(E)) == len(E) and max(E)-min(E) >= 6*V//7 and prim([e-min(E) for e in E])

def report(k, cands):
    worst = (0, None, None)
    for E, V in cands:
        E = sorted(set(e % V for e in E))
        if len(E) != k: continue
        E = [e-min(E) for e in E]
        if not valid(E, V) or has_gap(E, 1, V): continue
        js = jstar(E, V)
        if js and js > worst[0]: worst = (js, tuple(E[:6]), V)
    return worst

print("HARD adversarial max j* search (target: break j*=O(k))\n")
for k in (11, 13):
    allc = []
    # 1) perturbed APs at every scale + jitter
    for V in range(60, 4000):
        for d in range(V//14, V//7+1, max(1, (V//7-V//14)//4)):
            for _ in range(2):
                E = [(d*i + random.randint(-1, 1)) for i in range(k)]
                allc.append((E, V))
    # 2) interleaved two-APs
    for V in range(200, 4000, 7):
        d1 = random.randint(V//15, V//8); d2 = random.randint(V//15, V//8)
        h = k//2
        E = [d1*i for i in range(h)] + [d2*i + d1//2 for i in range(k-h)]
        allc.append((E, V))
    # 3) geometric / structured mod small prime
    for V in range(200, 4000, 11):
        g = random.randint(2, 6)
        E = [(g**i) % V for i in range(k)]
        allc.append((E, V))
    # 4) Vmax highly composite (many factors) with AP
    for V in [720, 840, 1260, 2520, 5040, 27720, 55440, 720720 % 100000]:
        for d in range(max(1, V//14), V//7+1, max(1, (V//7-V//14)//8)):
            E = [d*i for i in range(k)]
            allc.append((E, V))
    w = report(k, allc)
    print(f"k={k}: scanned ~{len(allc)} adversarial clusters; MAX j* = {w[0]}  (bound ~k={k}); "
          f"at Vmax={w[2]}, cluster~{w[1]}")
print("\n=> if max j* stays ~k under this hard search, j*=O(k) is robust.")
