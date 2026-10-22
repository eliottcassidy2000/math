# Audit of HYP-9364 evidence bullet 2: kappa/rho band over circulant classes and rank correlation.
# rho for Cay(Z_n, D): skew eigenvalues i*sigma_j, sigma_j = 2 sum_{s in D} sin(2 pi j s / n), j = 1..(n-1)/2.
import re, glob
from fractions import Fraction
from math import factorial
from mpmath import mp, mpf, sin, pi
mp.dps = 30
base = '/Users/e/Documents/GitHub/math-wt-fanopsp/04-computation/experiments/'
data = {}
for f in glob.glob(base + 'tournament_pseudoprimes_20261021/circ_*.out'):
    for line in open(f):
        m = re.match(r'n=\s*(\d+) S=\{([\d,]+)\} mult=\d+ H=(\d+)', line)
        if m:
            n = int(m.group(1)); D = tuple(int(x) for x in m.group(2).split(',')); H = int(m.group(3))
            data.setdefault(n, {})[D] = H
for line in open(base + 'lattice_bridges_20261021/circ25_H.txt'):
    m = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
    if m:
        D = tuple(int(x) for x in m.group(1).split(',')); data.setdefault(25, {})[D] = int(m.group(2))
def rho(n, D):
    r = mpf(1)
    for j in range(1, (n - 1) // 2 + 1):
        s = 2 * sum(sin(2 * pi * j * x / n) for x in D)
        r *= (n * n + s * s) / (n * n - s * s)
    return r
def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v); i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]: j += 1
        for k in range(i, j + 1): r[order[k]] = (i + j) / 2
        i = j + 1
    return r
def spearman(a, b):
    ra, rb = ranks(a), ranks(b); n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
    return num / den
for n in sorted(data):
    if n < 7: continue
    Ds = list(data[n]); H = [data[n][D] for D in Ds]
    kap = [mpf(Fraction(h * 2 ** (n - 1), factorial(n)).numerator) / Fraction(h * 2 ** (n - 1), factorial(n)).denominator for h in H]
    rh = [rho(n, D) for D in Ds]
    q = [k / r for k, r in zip(kap, rh)]
    Rn = tuple(range(1, (n - 1) // 2 + 1))
    iR = Ds.index(Rn) if Rn in Ds else None
    argmax_rho = max(range(len(Ds)), key=lambda i: rh[i]); argmax_H = max(range(len(Ds)), key=lambda i: H[i])
    print(f"n={n:2d} classes={len(Ds):3d} kappa/rho in [{float(min(q)):.4f}, {float(max(q)):.4f}]  n(1-min)={float(n*(1-min(q))):.3f} n(1-max)={float(n*(1-max(q))):.3f}  spearman(kappa,rho)={spearman([float(x) for x in kap],[float(x) for x in rh]):+.3f}  argmax rho is R_n: {argmax_rho == iR}  argmax H is R_n: {argmax_H == iR}")
