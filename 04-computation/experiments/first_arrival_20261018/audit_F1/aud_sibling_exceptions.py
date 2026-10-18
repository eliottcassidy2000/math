#!/usr/bin/env python3
"""audit F1: exact description of the exceptions to the sibling case analysis (A)/(B) in the proof of THM-4690(4).
With y = 3x+1 (x = 1 mod 8, so y = 4 mod 24), Y = 6N+1, level L >= 1 of x:
  (A) x ends at level L, sibling NOT (same level, not ending)  <=>  y = ceil(Y / (2*4^(L-1))) - 1
  (B) x not ending at level L, sibling NOT (level L-1, ending) <=>  y = ceil(Y / 4^(L-1)) - 1
(derived in REPORT.md). Hence at most one exception per level and case. Brute-force check over all N <= NMAX."""
import sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
def lev(x, N):
    L, b = 0, x
    while b < 2 * N: L += 1; b = 4 * b + 1
    return L, b <= 4 * N
def ceil_div(a, b): return -(-a // b)
tot = bad = mx = 0
for N in range(1, NMAX + 1):
    Y = 6 * N + 1; per = {}
    for x in range(1, 4 * N + 1, 8):
        if 2 * x + 1 > 4 * N: break
        (L, e), (L2, e2) = lev(x, N), lev(2 * x + 1, N)
        if L == 0: continue
        y = 3 * x + 1
        if e:  exc = not (L2 == L and not e2);     pred = (y == ceil_div(Y, 2 * 4 ** (L - 1)) - 1)
        else:  exc = not (L2 == L - 1 and e2);     pred = (y == ceil_div(Y, 4 ** (L - 1)) - 1)
        if exc != pred: bad += 1
        if exc: tot += 1; per[(L, e)] = per.get((L, e), 0) + 1
    if per: mx = max(mx, max(per.values()))
print('N <= %d: %d exceptions in total, max per (level, case) = %d, mismatches with the exact criterion: %d' % (NMAX, tot, mx, bad))
sys.exit(1 if bad else 0)
