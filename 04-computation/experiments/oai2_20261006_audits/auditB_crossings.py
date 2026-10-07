#!/usr/bin/env python3
"""Audit B, item 3: book-drawing crossing counts, by brute force over 4-sets (independent code).

Spine = Z_N in cyclic order; an edge {i, j} drawn in a page crosses another edge of the same page iff
their endpoints interleave. For a 4-set a < b < c < d the only interleaved pairing is {a,c},{b,d}.

(A) DDS drawing of K_n: page(i,j) = 0 iff (i + j mod n) in a contiguous block of floor(n/2) or ceil(n/2)
    consecutive sums; count crossings, compare with Z(n). Also the minimum over all cyclic offsets/sizes.
(B) Parity-bipartite drawing of K_{m,m} on Z_{2m}: edges = odd-sum pairs; the m odd classes
    s = 2r+1 (r in Z_m) split contiguously: r in [0, floor(m/2)) on page 0, the rest on page 1.
    Count crossings, compare with Z(m,m) = d_m^2 and with the telescoping formula; check the
    per-class-pair crossing number G(kappa) = 2 kappa (m - kappa) - m and the same-page pair count
    m - 2 kappa directly.
(C) DDS drawing of K_{2m} restricted: number of crossings whose two edges are both even-odd.
"""
import itertools, sys
from math import comb


def Zn(n):
    return (n // 2) * ((n - 1) // 2) * ((n - 2) // 2) * ((n - 3) // 2) // 4


def d(r):
    return (r // 2) * ((r - 1) // 2)


def crossings(N, page, edge_ok):
    """page(i,j) -> 0/1; edge_ok(i,j) -> bool. Count interleaved same-page edge pairs."""
    c = 0
    P = {}
    for i in range(N):
        for j in range(i + 1, N):
            if edge_ok(i, j):
                P[(i, j)] = page(i, j)
    for a, b, cc, dd in itertools.combinations(range(N), 4):
        e1, e2 = (a, cc), (b, dd)
        if e1 in P and e2 in P and P[e1] == P[e2]:
            c += 1
    return c


def dds_page(n, start, size):
    blk = {(start + k) % n for k in range(size)}
    return lambda i, j: 0 if (i + j) % n in blk else 1


ok = True
# (A) DDS for K_n
bad = []
for n in range(4, 61):
    c = crossings(n, dds_page(n, 0, n // 2), lambda i, j: True)
    c2 = crossings(n, dds_page(n, 0, (n + 1) // 2), lambda i, j: True)
    if c != Zn(n) or c2 != Zn(n):
        bad.append((n, c, c2, Zn(n)))
print("(A) DDS contiguous half-split of K_n has exactly Z(n) crossings for 4 <= n <= 60:", not bad, bad[:5])
ok &= not bad

# (B) parity-bipartite K_{m,m}
badB = []
for m in range(2, 41):
    N = 2 * m
    p0 = m // 2
    def page(i, j, m=m, p0=p0):
        s = (i + j) % (2 * m)
        r = (s - 1) // 2
        return 0 if r < p0 else 1
    c = crossings(N, page, lambda i, j: (i + j) % 2 == 1)
    K = (m - 1) // 2
    tele = sum((2 * k * (m - k) - m) * (m - 2 * k) for k in range(1, K + 1))
    if not (c == d(m) ** 2 == tele == K * K * (m - K - 1) ** 2):
        badB.append((m, c, d(m) ** 2, tele))
print("(B) parity-bipartite contiguous split has exactly d_m^2 = Z(m,m) crossings, = telescoping sum, 2 <= m <= 40:",
      not badB, badB[:5])
ok &= not badB

# (B') per-class-pair crossing counts G(kappa) and same-page pair counts
badG = []
for m in range(2, 26):
    N = 2 * m
    # crossings between class r1 and class r2 (odd sums 2r+1)
    X = {}
    for a, b, cc, dd in itertools.combinations(range(N), 4):
        if (a + cc) % 2 == 1 and (b + dd) % 2 == 1:
            r1 = (((a + cc) % N) - 1) // 2
            r2 = (((b + dd) % N) - 1) // 2
            key = tuple(sorted((r1, r2)))
            X[key] = X.get(key, 0) + 1
    for r1 in range(m):
        for r2 in range(r1, m):
            kap = min(r2 - r1, m - (r2 - r1))
            want = 0 if kap == 0 else 2 * kap * (m - kap) - m
            if X.get((r1, r2), 0) != want:
                badG.append((m, r1, r2, X.get((r1, r2), 0), want))
    p0 = m // 2
    for kap in range(1, m // 2 + 1):
        same = sum(1 for r1 in range(m) for r2 in range(r1 + 1, m)
                   if min(r2 - r1, m - r2 + r1) == kap and ((r1 < p0) == (r2 < p0)))
        if same != m - 2 * kap:
            badG.append(('same-page', m, kap, same, m - 2 * kap))
print("(B') class pair at cyclic distance kappa has G(kappa) = 2kappa(m-kappa)-m crossings (0 within a class); "
      "same-page pairs at distance kappa = m - 2kappa; 2 <= m <= 25:", not badG, badG[:5])
ok &= not badG

# (C) DDS of K_{2m}: crossings with both edges even-odd
for n in (8, 14):
    m = n // 2
    P = {}
    page = dds_page(n, 0, n // 2)
    for i in range(n):
        for j in range(i + 1, n):
            P[(i, j)] = page(i, j)
    tot = bb = 0
    for a, b, cc, dd in itertools.combinations(range(n), 4):
        if P[(a, cc)] == P[(b, dd)]:
            tot += 1
            if (a + cc) % 2 == 1 and (b + dd) % 2 == 1:
                bb += 1
    print(f"(C) n = {n}: DDS crossings {tot} (Z(n) = {Zn(n)}), bipartite-bipartite {bb} (d_m^2 = {d(m)**2})")
    ok &= (tot == Zn(n) and bb == d(m) ** 2)

# (D) the minimum over ALL 2-colourings of the m odd classes (class-colouring minimum), small m
for m in range(3, 9):
    N = 2 * m
    X = {}
    for a, b, cc, dd in itertools.combinations(range(N), 4):
        if (a + cc) % 2 == 1 and (b + dd) % 2 == 1:
            r1 = (((a + cc) % N) - 1) // 2
            r2 = (((b + dd) % N) - 1) // 2
            X[(r1, r2)] = X.get((r1, r2), 0) + 1
    best = min(sum(v for (r1, r2), v in X.items() if ((col >> r1) & 1) == ((col >> r2) & 1)) for col in range(1 << m))
    print(f"(D) m = {m}: min over all 2^{m} class colourings = {best}, Z(m,m) = {d(m)**2}")
    ok &= best == d(m) ** 2
print("ALL OK" if ok else "SOME FAILED")
