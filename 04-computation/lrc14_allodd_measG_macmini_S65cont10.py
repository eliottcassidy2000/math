from itertools import combinations
from math import gcd
from functools import reduce
from fractions import Fraction as F
import random
random.seed(92)
def measG(S, theta=F(1,7)):
    """meas{x in [0,1) : all ||o x|| >= theta} exactly: complement = union of arcs."""
    arcs = []
    for o in S:
        for k in range(o + 1):
            lo, hi = F(k, o) - theta / o, F(k, o) + theta / o
            arcs.append((max(lo, F(0)), min(hi, F(1))))
    arcs.sort()
    bad = F(0); cur_lo, cur_hi = arcs[0]
    for lo, hi in arcs[1:]:
        if lo <= cur_hi: cur_hi = max(cur_hi, hi)
        else: bad += cur_hi - cur_lo; cur_lo, cur_hi = lo, hi
    bad += cur_hi - cur_lo
    return 1 - bad
# exhaustive small caps
for cap in (21, 31, 41):
    worst = None; n = 0
    odds = list(range(1, cap + 1, 2))
    for S in combinations(odds, 6):
        if reduce(gcd, S) != 1: continue
        n += 1
        m = measG(list(S))
        if worst is None or m < worst[0]: worst = (m, S)
    print(f"cap {cap}: {n} sets; min meas(G_1/7) = {worst[0]} = {float(worst[0]):.5f} at {worst[1]}  (vs 1/7 = 0.14286)")
# adversarial larger
worst = None
for _ in range(4000):
    S = sorted(random.sample(range(1, 200, 2), 6))
    if reduce(gcd, S) != 1: continue
    m = measG(S)
    if worst is None or m < worst[0]: worst = (m, S)
print(f"adversarial cap 199 (4k samples): min meas = {float(worst[0]):.5f} at {worst[1]}")
