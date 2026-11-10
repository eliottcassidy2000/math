#!/usr/bin/env python3
"""Digit-degree law on a vertex cube.  H(R + v_X) = H(R) + sum_{i in X} a_i - sum_{{i,j} in X} mu_ij is an integer
polynomial of degree 2 in the 0/1 coordinates of X.  Claim: for k >= 1 the k-th binary digit of H(R + v_X), as a
function of X, is an F_2-polynomial of degree at most 2^k.  This script computes the exact algebraic degree of digits
1, 2, 3 (Moebius transform over F_2) for every tournament class R of order m, with H(R + v_X) counted DIRECTLY
(subset DP on the (m+1)-vertex tournament), and reports the distribution of degrees."""
import subprocess, sys
from collections import Counter

def ham(n, out):
    inn = [sum(1 << w for w in range(n) if out[w] >> v & 1) for v in range(n)]
    f = [[0] * n for _ in range(1 << n)]
    for S in range(1, 1 << n):
        for v in range(n):
            if not S >> v & 1: continue
            R = S & ~(1 << v)
            if R == 0: f[S][v] = 1; continue
            a = 0; mm = R & inn[v]
            while mm: w = (mm & -mm).bit_length() - 1; mm &= mm - 1; a += f[R][w]
            f[S][v] = a
    return sum(f[(1 << n) - 1])

def degree(vals, m):
    """algebraic degree over F_2 of the Boolean function X -> vals[X]"""
    a = list(vals)
    for i in range(m):
        for X in range(1 << m):
            if X >> i & 1: a[X] ^= a[X ^ (1 << i)]
    return max((bin(X).count('1') for X in range(1 << m) if a[X]), default=-1)

m = int(sys.argv[1])
res = subprocess.run(['/opt/homebrew/bin/gentourng', '-q', str(m)], capture_output=True, text=True).stdout.split()
deg = {1: Counter(), 2: Counter(), 3: Counter()}; bad = 0
for line in res:
    out = [0] * m; k = 0
    for i in range(m):
        for j in range(i + 1, m):
            if line[k] == '1': out[i] |= 1 << j
            else: out[j] |= 1 << i
            k += 1
    vals = []
    for X in range(1 << m):
        o2 = [out[u] | (0 if X >> u & 1 else 1 << m) for u in range(m)] + [X]
        vals.append(ham(m + 1, o2))
    for kk in (1, 2, 3):
        d = degree([(v >> kk) & 1 for v in vals], m); deg[kk][d] += 1
        if d > 2 ** kk: bad += 1
print("order m = %d, %d classes: degree of digit 1: %s ; digit 2: %s ; digit 3: %s ; violations of degree <= 2^k: %d"
      % (m, len(res), dict(sorted(deg[1].items())), dict(sorted(deg[2].items())), dict(sorted(deg[3].items())), bad))
