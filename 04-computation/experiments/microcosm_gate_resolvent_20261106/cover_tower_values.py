# mac-mini-2026-11-06-microcosm. HYP-9490 engine: exact values H(Gamma^r T) via the cover-polynomial value recursion (Session210 operator) up to total order N ~ 1201.
import itertools, math, sys
from fractions import Fraction
from collections import Counter
def forward_counts(n, arc):
    cnt = Counter()
    for p in itertools.permutations(range(n)):
        cnt[sum(1 for i in range(n - 1) if arc[p[i]][p[i + 1]])] += 1
    return cnt
def P_value(cnt, n, q):
    """P_T(q) = sum_pi binom(q + f(pi), n) evaluated at integer q (generalized binomial, exact)."""
    tot = Fraction(0)
    for f, c in cnt.items():
        x = q + f; num = 1
        for i in range(n): num *= (x - i)
        tot += Fraction(num, math.factorial(n)) * c
    assert tot.denominator == 1
    return int(tot)
def tower_H(cnt, n, r):
    lo, hi = 1 - r, 1 + r
    V = {j: P_value(cnt, n, j) for j in range(lo, hi + 1)}
    for i in range(r):
        lo += 1; hi -= 1
        V = {j: j * (V[j + 1] + j * V[j] - V[j - 1]) for j in range(lo, hi + 1)}
    return V[1]
cores = {
    'K1': (1, [[False]]),
    'TT2': (2, [[False, True], [False, False]]),
    'TT3': (3, [[i < j for j in range(3)] for i in range(3)]),
    'C3': (3, [[(j - i) % 3 == 1 for j in range(3)] for i in range(3)]),
    'TT5': (5, [[i < j for j in range(5)] for i in range(5)]),
    'R5': (5, [[(j - i) % 5 in (1, 2) for j in range(5)] for i in range(5)]),
    'QR7': (7, [[(j - i) % 7 in (1, 2, 4) for j in range(7)] for i in range(7)]),
    'TT7': (7, [[i < j for j in range(7)] for i in range(7)]),
}
cnts = {k: forward_counts(n, a) for k, (n, a) in cores.items()}
# check against brute-force tower values
print('K1 tower:', [tower_H(cnts['K1'], 1, r) for r in range(8)])
print('TT3 tower:', [tower_H(cnts['TT3'], 3, r) for r in range(6)])
for N in (101, 201, 401, 801, 1201):
    ref = tower_H(cnts['K1'], 1, (N - 1) // 2)
    row = []
    for name in ('C3', 'TT3', 'R5', 'TT5', 'QR7', 'TT7'):
        n = cores[name][0]; h = tower_H(cnts[name], n, (N - n) // 2)
        dev = Fraction(h, ref) - 1
        row.append('%s %+.3e (N^2 dev %+.3f)' % (name, float(dev), float(dev) * N * N))
    print('N = %4d: ' % N + ' | '.join(row), flush=True)
