# Test of the "self-blow-up principle": for a regular tournament T on n vertices with skew matrix S = A - A^T,
#   kappa(T) = H(T) 2^(n-1)/n!   versus   rho(T) = det(nI - S)/|det(inI - S)| = prod_k (n^2 + s_k^2)/(n^2 - s_k^2)
# (Wormald / Bhattacharjee blow-up gain).  All circulant classes, odd n <= 25.  mac-mini-2026-10-22-sierp.
import re, glob, math
from math import gcd, factorial
from collections import defaultdict
rows = defaultdict(dict)
for f in sorted(glob.glob('../tournament_pseudoprimes_20261021/circ_*.out')):
    for line in open(f):
        m = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
        if m: rows[int(m.group(1))][tuple(sorted(int(x) for x in m.group(2).split(',')))] = int(m.group(4))
for line in open('../lattice_bridges_20261021/circ25_H.txt'):
    m = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
    if m: rows[25][tuple(int(x) for x in m.group(1).split(','))] = int(m.group(2))
def rho(n, S):
    r = 1.0
    for j in range(1, (n - 1) // 2 + 1):
        s = 2 * sum(math.sin(2 * math.pi * j * x / n) for x in S)
        r *= (n * n + s * s) / (n * n - s * s)
    return r
def spearman(a, b):
    def rk(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
        for k, i in enumerate(o): r[i] = k
        return r
    ra, rb = rk(a), rk(b); n = len(a); ma = (n - 1) / 2
    num = sum((x - ma) * (y - ma) for x, y in zip(ra, rb)); den = sum((x - ma) ** 2 for x in ra)
    return num / den
for n in sorted(rows):
    if n < 7: continue
    ks, rs = [], []
    for S, H in rows[n].items():
        ks.append(H * 2 ** (n - 1) / factorial(n)); rs.append(rho(n, S))
    ratios = [k / r for k, r in zip(ks, rs)]
    imax_k = max(range(len(ks)), key=lambda i: ks[i]); imax_r = max(range(len(rs)), key=lambda i: rs[i])
    print(f"n={n:2d} ({len(ks):3d} classes): kappa/rho in [{min(ratios):.4f}, {max(ratios):.4f}], spearman(kappa, rho) = {spearman(ks, rs):+.3f}, argmax agree: {imax_k == imax_r}")
