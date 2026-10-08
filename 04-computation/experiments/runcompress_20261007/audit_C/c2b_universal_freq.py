#!/usr/bin/env python3
"""audit_C item 2 (supplement): ladder-type frequencies of merges from the universal state, independent replication of the
session's sampling design (K in [5,40], J in [3,12], r in {1,2}, t = 3^-(K-1) mod 2^(j-1) + 2^(j-1) * (200 random bits)),
own code; horizon 3000 post-run Terras steps. Compare with THM-4603: child i=1,2,3: 870, 50, 18; source i=1..4: 672, 33, 7, 1
(1651 merges in 3000 sources)."""
import random
from collections import Counter
rnd = random.Random(99991)
def v2(x): return (x & -x).bit_length() - 1
def T(x): return (3 * x + 1) >> 1 if x & 1 else x >> 1
kinds = Counter(); merges = 0; N = 3000
for trial in range(N):
    K = rnd.randint(5, 40); J = rnd.randint(3, 12); r = rnd.choice((1, 2)); j = 2 * J + r
    M = 1 << (j - 1)
    t0 = pow(3, -(K - 1), M)
    while True:
        t = t0 + M * rnd.getrandbits(200)
        if t & 1 and v2(3 ** (K - 1) * t - 1) == j - 1:
            break
    x = 2 * 3 ** (K - 1) * t - 1; y = (x + 1) // 27 - 1
    u, v, k = x, y, 3
    for s in range(2 * J):
        k += (u & 1) - (v & 1); u, v = T(u), T(v)
    assert k == 3 and u - 1 == 27 * (v - 1)
    zu = zv = None
    for s in range(3000):
        if u & 1: zu = u
        if v & 1: zv = v
        k += (u & 1) - (v & 1); u, v = T(u), T(v)
        if u == v:
            if k == 0:
                merges += 1
                cu, cv = v2(3 * zu + 1), v2(3 * zv + 1)
                i = abs(cv - cu) // 2
                kinds[('child' if cv > cu else 'source', i)] += 1
            break
        if u == 1 or v == 1:
            break
print(f"{merges} merges in {N} sources ({merges/N:.3f}; session 1651/3000 = 0.550)")
print("by (kind, i):", dict(sorted(kinds.items())))
tot = sum(kinds.values())
print(f"i = 1 share {sum(c for (kd, i), c in kinds.items() if i == 1)/tot:.3f} (session {(870+672)/1651:.3f}); "
      f"source share {sum(c for (kd, i), c in kinds.items() if kd == 'source')/tot:.3f} (session {(672+33+7+1)/1651:.3f})")
