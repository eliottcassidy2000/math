"""Which side do merges of the coupled pair y, x = 3*2^v*y + 1 come from (session opus-2026-10-07-S20)?
Note: 05-knowledge/results/two_readers_exponent_correlation_20261007.md (section 3, merges).
One step before the merge (first t with x_t = y_(t+1)) the relation is a sibling relation, x_t = 4^k y_(t+1) + (4^k-1)/3 (L = 2k > 0)
or y_(t+1) = 4^k x_t + (4^k-1)/3 (L = -2k < 0).  Counts by merge time and by the start offset L_0 = v + A_0 (same 12,000 pairs and
seed as two_readers_correlation_20261007.py).  NUMERICAL; the parity/nonzero check is exact.  Prints ALL CHECKS PASSED."""
import random
from multiprocessing import Pool
from collections import Counter

BITS, NS, TMAX, SEED = 3000, 12000, 600, 2007


def v2(n):
    return (n & -n).bit_length() - 1


def run(s):
    rng = random.Random(SEED * 1000003 + s)
    y = rng.getrandbits(BITS) | 1
    v = 2
    while rng.random() < 0.5:
        v += 1
    x = 3 * (1 << v) * y + 1
    a0 = v2(3 * y + 1)
    yy = (3 * y + 1) >> a0
    L = v + a0
    L0 = L
    xx = x
    Lprev = None
    for t in range(TMAX):
        if xx == yy:
            return v, L0, t, Lprev
        b = v2(3 * xx + 1); a = v2(3 * yy + 1)
        Lprev = L
        xx = (3 * xx + 1) >> b; yy = (3 * yy + 1) >> a
        L += a - b
    return v, L0, None, None


if __name__ == '__main__':
    with Pool(12) as pool:
        res = pool.map(run, range(NS), chunksize=50)
    allL = [Lp for v, L0, tau, Lp in res if tau is not None]
    ok = all(l != 0 and l % 2 == 0 for l in allL)
    print(('  ok   ' if ok else '  FAIL ') + 'every merging step has L even and nonzero (%d merges)' % len(allL))
    bins = [(1, 5), (6, 20), (21, 60), (61, 200), (201, 600)]
    print('merge side by merge time tau (L at the merging step):')
    for lo, hi in bins:
        c = Counter(Lp for v, L0, tau, Lp in res if tau is not None and lo <= tau <= hi)
        pos = sum(n for l, n in c.items() if l > 0); neg = sum(n for l, n in c.items() if l < 0)
        print(f'  tau in [{lo},{hi}]: +side {pos}, -side {neg}, ratio -/+ = {neg / max(pos, 1):.2f}   detail {dict(sorted(c.items()))}')
    print('merge side by start offset L0 = v + A_0:')
    for L0v in range(3, 9):
        c = Counter(Lp for v, L0, tau, Lp in res if tau is not None and L0 == L0v)
        pos = sum(n for l, n in c.items() if l > 0); neg = sum(n for l, n in c.items() if l < 0)
        tot = sum(1 for v, L0, tau, Lp in res if L0 == L0v)
        print(f'  L0 = {L0v}: pairs {tot}, merged {pos + neg}, +side {pos}, -side {neg}')
    print('ALL CHECKS PASSED' if ok else 'FAILURES')
