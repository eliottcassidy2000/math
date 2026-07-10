from itertools import combinations
from math import gcd
from functools import reduce
# exact M for 6-runner sets via pair-sum events (THM-668, n=7 => threshold 1/7)
def exact_M(S):
    bestn, bestd = 0, 1
    qs = {S[i] + S[j] for i in range(len(S)) for j in range(i, len(S))}
    for q in qs:
        for p in range(1, q):
            m = min(min(v * p % q, q - v * p % q) for v in S)
            if m * bestd > bestn * q:
                bestn, bestd = m, q
    return bestn, bestd
# all-odd primitive 6-sets, cap sweep: min M and whether 1/7 is attained
for cap in (25, 35, 45):
    worst = None
    n = 0
    odds = list(range(1, cap + 1, 2))
    for S in combinations(odds, 6):
        if reduce(gcd, S) != 1: continue
        n += 1
        bn, bd = exact_M(list(S))
        # compare M = bn/bd with 1/7: cushion = M - 1/7
        if worst is None or bn * worst[1] * 7 < worst[0] * bd * 7 - 0:
            if worst is None or bn * worst[1] < worst[0] * bd:
                worst = (bn, bd, S)
    bn, bd, S = worst
    from fractions import Fraction as F
    M = F(bn, bd)
    print(f"cap {cap}: {n} all-odd primitive 6-sets; min M = {M} = {float(M):.5f} at {S}; "
          f"cushion over 1/7: {M - F(1,7)} ({float(M - F(1,7)):.5f})")
