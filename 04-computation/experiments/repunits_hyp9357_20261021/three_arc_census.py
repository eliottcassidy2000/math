"""Three reversed arcs (session opus-2026-10-21-S51; THM-4748 (6), note section 2).
For the transitive order on m vertices with three arcs reversed (spans >= 2), every configuration with m <= 13
(exact Hamiltonian-path DP): which repunits 2^k - 1 occur and at which least order, and the odd parts of H + 1 among
the values = 7 (mod 8) (with two reversed arcs they lie on the lines c*2^j - 1 with c in {3, 5} plus 15, 23;
THM-4748 (2)).  Prints ALL CHECKS PASSED."""
import itertools, time
from collections import Counter
import numpy as np
import numba as nb

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


@nb.njit(cache=True)
def ham(n, out):
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            m = out[v] & ~S
            while m:
                w = 0
                while not (m >> w) & 1:
                    w += 1
                m &= m - 1
                dp[S | (1 << w), w] += c
    t = 0
    for v in range(n):
        t += dp[N - 1, v]
    return t


@nb.njit(cache=True)
def run_all(m, arcs, combos, out_h):
    base = np.zeros(m, dtype=np.int64)
    for i in range(m):
        for j in range(i + 1, m):
            base[i] |= 1 << j
    for t in range(combos.shape[0]):
        out = base.copy()
        for k in range(3):
            j = arcs[combos[t, k], 0]; i = arcs[combos[t, k], 1]
            out[i] &= ~(1 << j)
            out[j] |= 1 << i
        out_h[t] = ham(m, out)


def oddpart(x):
    while x % 2 == 0:
        x //= 2
    return x


mers = {2 ** k - 1: k for k in range(2, 60)}
first = {}
maxc = {}
share = {}
for m in range(5, 14):
    arcs = np.array([(j, i) for i in range(m) for j in range(i + 2, m)], dtype=np.int64)
    combos = np.array(list(itertools.combinations(range(len(arcs)), 3)), dtype=np.int64)
    out_h = np.zeros(len(combos), dtype=np.int64)
    run_all(m, arcs, combos, out_h)
    vals = set(int(x) for x in out_h)
    for h in sorted(vals):
        if h in mers and h not in first:
            first[h] = m
    v7 = [h for h in vals if h % 8 == 7]
    maxc[m] = max(oddpart(h + 1) for h in v7)
    share[m] = len(v7) / len(vals)
    print(f'     m = {m:2d}: {len(combos):6d} configurations, {len(vals):5d} distinct H, {len(v7):4d} of them = 7 mod 8 '
          f'(share {share[m]:.3f}); largest odd part of H + 1 among them {maxc[m]}; repunits so far {dict(sorted(first.items()))} [{time.time() - T0:.0f}s]', flush=True)
check(all(first.get(2 ** k - 1) for k in range(4, 12)) and 2 ** 15 - 1 in first,
      'three reversed arcs reach 2^k - 1 for every 4 <= k <= 11 and for k = 15 (m <= 13): the bound of THM-4748 (3) is attained')
check(all(maxc[m] < maxc[m + 1] for m in range(7, 13)), 'with three reversed arcs the odd parts of H + 1 (H = 7 mod 8) grow without visible bound: no finite set of lines')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
