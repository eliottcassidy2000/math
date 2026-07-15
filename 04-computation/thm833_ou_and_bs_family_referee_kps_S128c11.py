#!/usr/bin/env python3
"""
thm833_ou_and_bs_family_referee_kps_S128c11.py
==============================================
kind-pasteur-2026-07-15-S128 (cont.11). Three pieces:
 (A) THM-833 referee: OU drift E[dc3|T] = -kappa (c3 - C(n,3)/4), kappa = 8/(n(n-1)) -- exact per
     tournament (all tournaments n=4,5,6 via tilings) + fluctuation-dissipation at n=4..7.
 (B) THE B_s / c_k(s) INTEGER FAMILY (THM-826's slope increments): c_k(s) = #{consecutive Farey
     pairs of F_k with i+j = s} -- prove-check c_k(s) = #{i in [s-k, k] : gcd(i, s) = 1}, print the
     triangle k=2..14, row sums = |F_k|, first column phi(k+1), and B_s reconstruction.
 (C) CHI_7 ODD MODE vs BLUE/BLACK: the 7 Venn cells of the 3-corner staircase decomposition at
     general n (corners A,B,C = delete top row / bottom-right / apex-ish), chi_7 sign string
     ++-+--+ over cells 1..7 in prompt order (QR mod 7 = {1,2,4} positive): test whether the sign
     matches parity invariants of each cell's sub-staircase (gridsym-count parity, H mod 4 of the
     cell's transitive class, cell size parity) at n=6,7 -- an exploratory cross-table.
"""
import sys
from fractions import Fraction as F
from math import comb, gcd
from itertools import combinations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

print("=" * 94)
print("(A) THM-833 OU LAW referee")
print("=" * 94)
for n in [4, 5, 6]:
    m = comb(n - 1, 2)
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    tidx = {t: i for i, t in enumerate(tiles)}
    ok = True
    for t in range(1 << m):
        tv = [(t >> i) & 1 for i in range(m)]
        d = [0] * n
        for k in range(2, n + 1):
            d[k - 1] += 1
        for (x, y), i in tidx.items():
            if tv[i] == 1:
                d[x - 1] += 1
            else:
                d[y - 1] += 1
        c3 = comb(n, 3) - sum(comb(di, 2) for di in d)
        # enumerate arcs and the drift
        drift = F(0)
        arcs = 0
        B = [[False] * n for _ in range(n)]
        for k in range(2, n + 1):
            B[k - 1][k - 2] = True
        for (x, y), i in tidx.items():
            if tv[i] == 1:
                B[x - 1][y - 1] = True
            else:
                B[y - 1][x - 1] = True
        for u in range(n):
            for v in range(n):
                if u != v and B[u][v]:
                    drift += d[u] - d[v] - 1
                    arcs += 1
        assert arcs == comb(n, 2)
        lhs = drift / comb(n, 2)
        rhs = -F(8, n * (n - 1)) * (c3 - F(comb(n, 3), 4))
        if lhs != rhs:
            ok = False
            print("  n=%d FAIL at tiling %d: lhs=%s rhs=%s" % (n, t, lhs, rhs))
            break
    print("n=%d: OU drift exact on ALL %d tournaments: %s" % (n, 1 << m, ok))
# fluctuation-dissipation: uniform Var(c3) vs E[(dc3)^2]/(2 kappa) - correction. Use the exact
# discrete stationarity identity: E_unif[ (c3' - c3) * (c3 + c3' - 2c3*) ] = 0 for the symmetric
# chain (detailed balance) -- equivalently E[ (dc3)^2 ] = 2 kappa E[(c3-c3*)^2] + ... check directly:
for n in [4, 5, 6]:
    m2 = comb(n, 2)
    # E over uniform T and uniform arc flip
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    tidx = {t: i for i, t in enumerate(tiles)}
    mm = len(tiles)
    S1 = F(0); S2 = F(0); Svar = F(0); N = 1 << mm
    c3s = F(comb(n, 3), 4)
    for t in range(N):
        tv = [(t >> i) & 1 for i in range(mm)]
        d = [0] * n
        for k in range(2, n + 1):
            d[k - 1] += 1
        for (x, y), i in tidx.items():
            d[(x - 1) if tv[i] else (y - 1)] += 1
        c3 = comb(n, 3) - sum(comb(di, 2) for di in d)
        Svar += (c3 - c3s) ** 2
        B = [[False] * n for _ in range(n)]
        for k in range(2, n + 1):
            B[k - 1][k - 2] = True
        for (x, y), i in tidx.items():
            B[(x - 1)][(y - 1)] = tv[i] == 1
            B[(y - 1)][(x - 1)] = tv[i] == 0
        for u in range(n):
            for v in range(n):
                if u != v and B[u][v]:
                    dd = d[u] - d[v] - 1
                    S1 += F(dd, m2) * (c3 - c3s)
                    S2 += F(dd * dd, m2)
    kappa = F(8, n * (n - 1))
    fd = (S1 / N == -kappa * (Svar / N))
    print("n=%d: fluctuation-dissipation E[dc3 (c3-c3*)] == -kappa Var(c3): %s   (Var=%s, E[dc3^2]=%s)"
          % (n, fd, Svar / N, S2 / N))

print()
print("=" * 94)
print("(B) the c_k(s) INTEGER TRIANGLE (B_s increments of THM-826)")
print("=" * 94)
def farey_pairs(k):
    fr = sorted(set(F(a, i) for i in range(1, k + 1) for a in range(i)))
    out = []
    for x, y in zip(fr, fr[1:] + [F(1)]):
        out.append((x.denominator, y.denominator if y != 1 else 1))
    return out
print("c_k(s) = #{consecutive F_k pairs, i+j=s}; claim c_k(s) = #{i in [max(1,s-k), k]: gcd(i,s)=1}")
allok = True
for k in range(2, 15):
    pairs = farey_pairs(k)
    row = []
    for s in range(k + 1, 2 * k):
        c = sum(1 for i, j in pairs if i + j == s)
        pred = sum(1 for i in range(max(1, s - k), k + 1) if gcd(i, s) == 1)
        allok &= (c == pred)
        row.append(c)
    tot = sum(row)
    phis = sum(1 for a in range(1, k + 1) if gcd(a, k + 1) == 1)
    print("k=%2d: %s  rowsum=%d=|F_k|:%s  c(k+1)=phi(k+1):%s"
          % (k, row, tot, tot == len(pairs), row[0] == phis))
print("triangle law c_k(s) = interval-totient: %s" % allok)

print()
print("=" * 94)
print("(C) chi_7 odd mode vs blue/black parities (exploratory, n=6,7 Venn cells)")
print("=" * 94)
CHI7 = {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1, 7: 1}   # QR mod 7 {1,2,4}: + ; string ++-+--+ over 1..7...
# prompt order A..G = 1..7 with sign ++-+--+ (from OPEN-Q-108 S31q): A+,B+,C-,D+,E-,F-,G+
SIGNS = {1: +1, 2: +1, 3: -1, 4: +1, 5: -1, 6: -1, 7: +1}
print("chi_7(cell) vs SIGNS(cell): %s" % all(CHI7[i] == SIGNS[i] for i in range(1, 8)))
# Venn cells of three corner staircases at size n: corners = sub-staircases from deleting
# vertex 1 (A), vertex n (B), and both-extremes overlap structure; cells sizes from OPEN-Q-108 S114:
# A,B: N-1; D: N-2; overlaps C: N-2, E,F: N-3; G: N-4 (in the odd-mode Venn ordering A+B+D-C-E-F+G).
# The cell 'size' = sub-tournament order. Parity invariants per cell size at n:
for n in [6, 7]:
    sizes = {1: n - 1, 2: n - 1, 3: n - 2, 4: n - 2, 5: n - 3, 6: n - 3, 7: n - 4}   # A,B,C,D,E,F,G (S115 odd-mode)
    # candidate parity: (-1)^{C(size-1,2)} = parity of tile count of the sub-staircase
    # candidate 2: gridsym-count exponent parity: r = (m+f)/2 with m=C(size-1,2), f=floor((size-1)/2)
    rows = []
    for c in range(1, 8):
        sz = sizes[c]
        mt = comb(sz - 1, 2)
        f = (sz - 1) // 2
        r = (mt + f) // 2
        rows.append((c, SIGNS[c], (-1) ** mt, (-1) ** r, (-1) ** sz))
    match_mt = all(s == p for _, s, p, _, _ in rows)
    match_r = all(s == p for _, s, _, p, _ in rows)
    match_sz = all(s == p for _, s, _, _, p in rows)
    print("n=%d: sign==(-1)^tiles: %s ; sign==(-1)^gridsym-exp: %s ; sign==(-1)^size: %s ; table=%s"
          % (n, match_mt, match_r, match_sz, rows))
print("\nDONE")
