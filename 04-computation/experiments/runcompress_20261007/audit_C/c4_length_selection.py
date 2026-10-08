#!/usr/bin/env python3
"""audit_C item 4: is the length-tercile enrichment of Mersenne orphans a selection effect of the partner criterion?

Fact: K' is a partner of K iff o(K') = o(K) and sigma_T(K) - K = sigma_T(K') - K'. So the post-run length c = sigma_T - K
is a CLASS INVARIANT: every non-orphan inherits c from a smaller exponent, and its normalised length c/K is deflated by
K_0/K (K_0 = the class start), whereas an orphan (= a class start) carries a fresh c.
Null model: keep the class membership of every odd K in [3, 12800] exactly as observed, but replace each class's offset by
the post-run length of an independent random odd integer R of the same size as 3^(K_0) - 1 (R uniform in [3^K_0/2, 3^K_0]).
Report the tercile orphan rates under the null (several replicas) next to the observed ones, and the orphans' observed
offsets as z-scores against the random-orbit law (mean and sd estimated from random R).
"""
import random, math, statistics
from collections import defaultdict
exec(open('c4_general_sources.py').read().split('N = int(sys.argv[1])')[0])   # ot()

data = {}
for line in open('c4_mersenne_12800_audit.txt'):
    K, o, t = map(int, line.split())
    data[K] = (o, t)
cls = {}            # K -> class key (o, offset)
start = {}          # class key -> first K
for K in range(2, 12801):
    o, t = data[K]
    key = (o, t - K)
    if key not in start:
        start[key] = K
    cls[K] = key
oddK = [K for K in range(3, 12801, 2)]
orphan = {K: start[cls[K]] == K for K in oddK}
assert sum(orphan.values()) == 135


def terc(lo, hi, offset_of):
    rows = [(offset_of(cls[K]) / K, orphan[K]) for K in oddK if lo <= K <= hi]
    xs = sorted(x for x, _ in rows)
    q1, q2 = xs[len(xs) // 3], xs[2 * len(xs) // 3]
    out = []
    for a, b in ((-1e9, q1), (q1, q2), (q2, 1e18)):
        sel = [o for x, o in rows if a <= x < b]
        out.append(100 * sum(sel) / len(sel))
    return out


obs = lambda key: key[1]
for lo, hi in ((1000, 6400), (1000, 12800)):
    print(f"observed terciles [{lo},{hi}]: " + ", ".join(f"{x:.2f}%" for x in terc(lo, hi, obs)))

rnd = random.Random(99)
# random-orbit law of the post-run length c for size 3^K0: c = sigma_T(R), R ~ 3^K0
needed = sorted(set(start[key] for key in set(cls[K] for K in oddK)))
for rep in range(5):
    newoff = {}
    for key, K0 in start.items():
        if K0 < 3:
            newoff[key] = key[1]
            continue
        hi = 3 ** K0
        R = rnd.randrange(hi // 2, hi) | 1
        newoff[key] = ot(R)[1]
    nl = lambda key: newoff[key]
    print(f"null replica {rep}: [1000,6400] " + ", ".join(f"{x:.2f}%" for x in terc(1000, 6400, nl)) +
          " | [1000,12800] " + ", ".join(f"{x:.2f}%" for x in terc(1000, 12800, nl)))

# z-scores of the orphans' actual offsets against the random-orbit law (mean/sd estimated per size)
zs = []
for K0 in [K for K in oddK if orphan[K] and K >= 1000]:
    hi = 3 ** K0
    samp = [ot(rnd.randrange(hi // 2, hi) | 1)[1] for _ in range(40)]
    mu, sd = statistics.mean(samp), statistics.stdev(samp)
    c = data[K0][1] - K0
    zs.append((c - mu) / sd)
print(f"orphans K >= 1000 (n = {len(zs)}): z-score of the actual post-run length vs random orbits of the same size: "
      f"mean {statistics.mean(zs):+.2f}, median {statistics.median(zs):+.2f}, fraction z > 0: {sum(z > 0 for z in zs)/len(zs):.2f}")
# the same z-score for ALL odd K (non-orphans use their inherited offset, compared to random orbits of size 3^K)
zall = []
for K in rnd.sample([K for K in oddK if K >= 1000 and not orphan[K]], 150):
    hi = 3 ** K
    samp = [ot(rnd.randrange(hi // 2, hi) | 1)[1] for _ in range(40)]
    mu, sd = statistics.mean(samp), statistics.stdev(samp)
    zall.append(((data[K][1] - K) - mu) / sd)
print(f"non-orphans (150 sampled, K >= 1000): z-score mean {statistics.mean(zall):+.2f}, median {statistics.median(zall):+.2f}")
