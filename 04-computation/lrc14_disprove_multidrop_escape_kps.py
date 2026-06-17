"""
DISPROVE-side probe for LRC(14) / inf_S L(S).

Attack the precise gap the PROVE side left open. The decoupling argument only
handles ONE element -> infinity (floor (6/7)meas(G_C) >= 1/143). The escape, if
any, must use MULTIPLE simultaneously-growing, arithmetically-coordinated speeds.

Leads:
  - drop-6 single family minimizes at the LARGE w=69 (L=19/10626 ~ 0.00179).
  - Can a k-drop family (k>=2) push L below 1/1260 (=0.000794), toward L->0?
"""
from fractions import Fraction as Fr
from math import gcd, lcm
from functools import reduce
from itertools import combinations

def danger_arcs(v):
    w = Fr(1, 14*v); A = []
    for k in range(v+1):
        c = Fr(k, v); lo = c - w; hi = c + w
        if lo < 0:
            A += [(Fr(0), hi), (1+lo, Fr(1))]
        elif hi > 1:
            A += [(lo, Fr(1)), (Fr(0), hi-1)]
        else:
            A.append((lo, hi))
    return A

_cache = {}
def darcs(v):
    r = _cache.get(v)
    if r is None:
        r = danger_arcs(v); _cache[v] = r
    return r

def L_exact(S):
    A = []
    for v in S:
        A.extend(darcs(v))
    A.sort(key=lambda t: (t[0], t[1]))
    tot = Fr(0); cl = ch = None
    for a, b in A:
        if b <= a: continue
        if ch is None:
            cl, ch = a, b
        elif a <= ch:
            if b > ch: ch = b
        else:
            tot += ch - cl; cl, ch = a, b
    if ch is not None:
        tot += ch - cl
    return 1 - tot

def gcd_list(S):
    return reduce(gcd, S)

THRESH = Fr(1, 1260)
AP = list(range(1, 14))

print("=== (A) Single-drop minimizers, extended window Wmax=600 ===")
Wmax = 600
single_best = {}
for e in range(1, 14):
    core = [x for x in AP if x != e]
    best = None; bestw = None; bestlcm = None
    for w in range(1, Wmax+1):
        if w in core: continue
        S = sorted(core + [w])
        L = L_exact(S)
        if L == 0:  # skip tight (trivial re-tightening or sporadic)
            continue
        if best is None or L < best:
            best = L; bestw = w; bestlcm = lcm(*S)
    single_best[e] = (best, bestw)
    print(f"  drop e={e:2d}: min L = {best} = {float(best):.6f} at w={bestw}, lcm={bestlcm}, below 1/1260? {best < THRESH}")

print()
print(f"  Global single-drop min: {min(v[0] for v in single_best.values())}")
print(f"  THRESH = 1/1260 = {float(THRESH):.6f}")
