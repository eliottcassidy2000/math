#!/usr/bin/env python3
"""Deletion-pair form of the second moment of the mod-4 bias (README section 4.1, remark after the tables).

For R on m labelled vertices let b_R = E_X (-1)^{q_R(X)}, q_R(X) = (H(R+v_X) - H(R))/2 mod 2 (Theorem C).
Since T - v' = R + v_X and T - v = R + v'_X' are conditionally independent given R = T - v - v',

    E_R[b_R^2] = E_T (-1)^{(H(T-v) - H(T-v'))/2},     T uniform on m+2 labelled vertices, v != v' fixed.

Route 1 (right side): every labelled tournament on n = m+2 vertices, H(T-v), H(T-v') by subset DP.
Route 2 (left side):  every labelled R on m vertices, every out-set X, H(R+v_X) by subset DP.
Route 3 (orders m = 2..7): from the rank/bias table out/betarank_labelled_2_7.out (b_R^2 = 2^-rank or 0).
Exact rational arithmetic.  Usage: python3 deletion_pair_corr.py [max n for routes 1-2, default 6]
"""
import sys, re, os
from fractions import Fraction
from itertools import combinations

def ham(n, out):
    """number of Hamiltonian paths; out[v] = bitmask of out-neighbours"""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, full + 1):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if c:
                nx = out[v] & ~S
                while nx:
                    w = (nx & -nx).bit_length() - 1; nx &= nx - 1
                    dp[S | (1 << w)][w] += c
    return sum(dp[full])

def labelled(n):
    prs = list(combinations(range(n), 2))
    for code in range(1 << len(prs)):
        out = [0] * n
        for k, (i, j) in enumerate(prs):
            if code >> k & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield out

def delete(out, v):
    """tournament with vertex v removed, relabelled 0..n-2"""
    n = len(out); keep = [u for u in range(n) if u != v]
    pos = {u: k for k, u in enumerate(keep)}
    return [sum(1 << pos[w] for w in keep if out[u] >> w & 1) for u in keep]

def route1(n):
    tot = 0; cnt = 0
    for out in labelled(n):
        a = ham(n - 1, delete(out, 0)); b = ham(n - 1, delete(out, 1))
        tot += -1 if ((a - b) // 2) % 2 else 1; cnt += 1
    return Fraction(tot, cnt)

def route2(m):
    tot = Fraction(0); cnt = 0
    for out in labelled(m):
        h = ham(m, out); s = 0
        for X in range(1 << m):
            o2 = [out[u] | (0 if X >> u & 1 else 1 << m) for u in range(m)] + [X]
            s += -1 if ((ham(m + 1, o2) - h) // 2) % 2 else 1
        tot += Fraction(s, 1 << m) ** 2; cnt += 1
    return tot / cnt

def route3(path):
    res = {}; m = None
    for line in open(path):
        g = re.match(r'order (\d+): classes (\d+)', line)
        if g: m = int(g.group(1)); res[m] = [Fraction(0), int(g.group(2))]; continue
        g = re.match(r'\s+rank (\d+): (\d+) classes with zero bias, (\d+) with bias', line)
        if g: res[m][0] += Fraction(int(g.group(3)), 1 << int(g.group(1)))
    return {m: s / c for m, (s, c) in res.items()}

nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
here = os.path.dirname(os.path.abspath(__file__))
r3 = route3(os.path.join(here, 'out', 'betarank_labelled_2_7.out'))
print("m   E_R[b_R^2] (table)      ratio to 2^-m   route 1 (two deletions, n = m+2)   route 2 (vertex cubes)")
ok = True
for m in sorted(r3):
    v = r3[m]; a = b = None
    if m + 2 <= nmax: a = route1(m + 2)
    if m + 1 <= nmax: b = route2(m)
    ok &= (a is None or a == v) and (b is None or b == v)
    print("%d   %-22s  %.4f          %-32s  %s" % (m, v, float(v * 2 ** m), a if a is not None else '-', b if b is not None else '-'))
print("ALL ROUTES AGREE" if ok else "MISMATCH")
