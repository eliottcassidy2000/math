#!/usr/bin/env python3
"""Which affine-linear relations over F_2 hold for the pair form beta_R (and its linear part ell_R) of EVERY
tournament R of order m?  Lemma 4.5 gives  sum_{i<j} beta_R(i,j) = 1  and  sum_i ell_i = 1.
This script takes tournaments of order m (all labelled ones for m <= 5, random labelled ones otherwise; seeded),
computes (ell_R, beta_R) from the subset DP, and computes the dimension of the affine span of
  (a) the vectors beta_R in F_2^C(m,2),      (b) the vectors (ell_R, beta_R) in F_2^(m + C(m,2)).
An affine span of dimension C(m,2) - 1 in (a) is a CERTIFICATE that (4.5) is the only affine relation on beta_R
at that order (a lower bound on a span needs no completeness of the sample).
Usage: python3 affine_span.py m [samples]"""
import sys, random
from itertools import combinations

def forms(n, out):
    full = (1 << n) - 1
    inn = [sum(1 << w for w in range(n) if out[w] >> v & 1) for v in range(n)]
    f = [[0] * n for _ in range(1 << n)]; g = [[0] * n for _ in range(1 << n)]
    for S in range(1, full + 1):
        for v in range(n):
            if not S >> v & 1: continue
            R = S & ~(1 << v)
            if R == 0: f[S][v] = g[S][v] = 1; continue
            a = 0; mm = R & inn[v]
            while mm: w = (mm & -mm).bit_length() - 1; mm &= mm - 1; a += f[R][w]
            b = 0; mm = R & out[v]
            while mm: w = (mm & -mm).bit_length() - 1; mm &= mm - 1; b += g[R][w]
            f[S][v] = a; g[S][v] = b
    M = [[0] * n for _ in range(n)]          # M[j][i]: P1 ends at j, P2 starts at i
    for j in range(n):
        for i in range(n):
            if i == j: continue
            rest = full & ~(1 << i) & ~(1 << j); tot = 0; S0 = rest
            while True:
                S = S0 | (1 << j); tot += f[S][j] * g[full ^ S][i]
                if S0 == 0: break
                S0 = (S0 - 1) & rest
            M[j][i] = tot
    ell = []
    for i in range(n):
        a = g[full][i] - f[full][i] + sum(M[j][i] for j in range(n) if j != i)
        assert a % 2 == 0; ell.append((a >> 1) & 1)
    beta = []
    for i, j in combinations(range(n), 2):
        mu = M[i][j] + M[j][i]; assert mu % 2 == 0; beta.append((mu >> 1) & 1)
    return ell, beta

class Span:
    """affine span over F_2 of bit-vectors given as ints: keeps a reduced basis of differences"""
    def __init__(self): self.base = None; self.piv = {}
    def add(self, v):
        if self.base is None: self.base = v; return
        d = v ^ self.base
        while d:
            h = d.bit_length() - 1
            if h in self.piv: d ^= self.piv[h]
            else: self.piv[h] = d; return
    def dim(self): return len(self.piv)

m = int(sys.argv[1]); samples = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
prs = list(combinations(range(m), 2)); npr = len(prs)
rng = random.Random(20261110 + m)
def tournaments():
    if m <= 5:
        for code in range(1 << npr): yield code
    else:
        for _ in range(samples): yield rng.getrandbits(npr)
A = Span(); B = Span(); cnt = 0; viol = 0
for code in tournaments():
    out = [0] * m
    for k, (i, j) in enumerate(prs):
        if code >> k & 1: out[i] |= 1 << j
        else: out[j] |= 1 << i
    ell, beta = forms(m, out); cnt += 1
    if sum(beta) % 2 != 1 or sum(ell) % 2 != 1: viol += 1
    vb = sum(b << k for k, b in enumerate(beta)); vl = sum(b << k for k, b in enumerate(ell))
    A.add(vb); B.add(vb | (vl << npr))
    if cnt % 200 == 0 and A.dim() == npr - 1 and B.dim() == npr + m - 2 and m > 5: break
print("m = %d: %d tournaments (%s); violations of (4.5): %d; affine span of beta_R: dimension %d of at most %d (= C(m,2) - 1);"
      " of (ell_R, beta_R): dimension %d of at most %d (= m + C(m,2) - 2)"
      % (m, cnt, "all labelled" if m <= 5 else "seeded random labelled", viol, A.dim(), npr - 1, B.dim(), npr + m - 2))
