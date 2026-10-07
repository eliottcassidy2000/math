#!/usr/bin/env python3
"""Book drawings of K_n and K_{m,m} by endpoint-sum classes, against the Harary-Hill / Zarankiewicz theorems claimed in
openai/math #165 (mac-mini, 2026-10-06).

Spine = a circle with vertices 0..N-1 in order; two chords on the same page cross iff their endpoints alternate.
DDS construction (Damiani-D'Antona-Salemi 1994; Blazek-Koman 1964; Shahrokhi et al.; de Klerk-Pasechnik-Salazar 2012
Section 5.1): edge ij lies in the matching M_s, s = i + j mod n; page 1 = M_0..M_(m-1), page 2 = the rest, m = floor(n/2).
This is the repo's THM-913 "parallel-class book drawing" (prior art recorded).

  1. The DDS 2-page drawing of K_n has exactly Z(n) crossings, 3 <= n <= 120 (proved for all n in DPS 2012 and in
     openai/math #165-K_n, Prop. "matching two-page drawing").
  2. Parity-bipartite version (THM-922): on Z_(2m) with parts = parities, K_{m,m} = the odd sum classes; the contiguous
     split of the odd classes has exactly d_m^2 = Z(m,m) crossings for every m. Proof: a crossing is a 4-set whose two
     alternating diagonals are both bipartite and on one page; with cyclic gaps g1..g4 this forces k = g1 + g3 = 2*kappa,
     there are G(kappa) = 2 kappa (m - kappa) - m gap choices and m - 2 min(kappa, m - kappa) same-page residues, so
     c = sum_(kappa=1)^K (2 kappa (m-kappa) - m)(m - 2 kappa) with K = floor((m-1)/2), and the sum telescopes:
     sum_(kappa=1)^K (2 kappa (m - kappa) - m)(m - 2 kappa) = K^2 (m - K - 1)^2 = d_m^2.
     The drawing is the DDS drawing of K_(2m) restricted to the even-odd edges.
  3. Consequences: with Abrego et al. 2012 (nu_2(K_n) = Z(n)), the class-coloring minimum for K_n is Z(n) for every n
     (unconditional; THM-922 (III) had n <= 14). CONDITIONAL on #165 (cr(K_{m,n}) = Z(m,n)): nu_2(K_{m,n}) = Z(m,n)
     for all m, n (the 2-page Zarankiewicz conjecture of de Klerk-Pasechnik-Salazar 2014), and THM-922's class-coloring
     minimum is Z(m,m) for every m.
  4. The rank inequality behind #165 (sum_(i,j) dim(U_i n V_j) <= floor(m^2/4) N when U_i n V_i = 0) at N = 1 is
     |A||B| <= floor(m^2/4) for disjoint A, B in [m]; random tests at small N (NUMERICAL sanity only).
Run: python3 oai2_20261006_crossing_books.py   (about 10 s)
"""
import itertools, random
import numpy as np
import sympy as sp

OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg, flush=True)
    OK &= bool(cond)


def Z(n):
    return (n // 2) * ((n - 1) // 2) * ((n - 2) // 2) * ((n - 3) // 2) // 4


def d(r):
    return (r // 2) * ((r - 1) // 2)


def crossings(pages):
    tot = 0
    for edges in pages:
        if len(edges) < 2:
            continue
        E = np.array(sorted((min(a, b), max(a, b)) for a, b in edges))
        a, b = E[:, 0][:, None], E[:, 1][:, None]
        c, e = E[:, 0][None, :], E[:, 1][None, :]
        tot += int((((a < c) & (c < b) & (b < e)) | ((c < a) & (a < e) & (e < b))).sum()) // 2
    return tot


def dds_kn(n):
    m = n // 2
    p1, p2 = [], []
    for i, j in itertools.combinations(range(n), 2):
        (p1 if (i + j) % n < m else p2).append((i, j))
    return p1, p2


def parity_kmm(m):
    N = 2 * m
    first = set(list(range(1, N, 2))[: m // 2])
    E = [(i, j) for i in range(0, N, 2) for j in range(1, N, 2)]
    return [e for e in E if (e[0] + e[1]) % N in first], [e for e in E if (e[0] + e[1]) % N not in first]


print("1. DDS two-page drawing of K_n")
check(all(crossings(dds_kn(n)) == Z(n) for n in range(3, 121)), "c(DDS K_n) = Z(n) = (1/4) floor(n/2) floor((n-1)/2) floor((n-2)/2) floor((n-3)/2), 3 <= n <= 120")

print("2. parity-bipartite endpoint-sum drawing of K_{m,m}")
check(all(crossings(parity_kmm(m)) == d(m) ** 2 for m in range(1, 61)), "c = d_m^2 = Z(m,m) for 1 <= m <= 60")
k, mm, K = sp.symbols("k m K", integer=True)
S = sp.factor(sp.summation((2 * k * (mm - k) - mm) * (mm - 2 * k), (k, 1, K)))
check(sp.simplify(S - K ** 2 * (mm - K - 1) ** 2) == 0, f"sum_(kappa=1)^K (2 kappa (m-kappa) - m)(m - 2 kappa) = {S}")
s = sp.Symbol("s", integer=True, positive=True)
check(sp.simplify(S.subs({mm: 2 * s + 1, K: s}) - (s * s) ** 2) == 0 and sp.simplify(S.subs({mm: 2 * s, K: s - 1}) - (s * (s - 1)) ** 2) == 0,
      "at K = floor((m-1)/2): K^2 (m-K-1)^2 = d_m^2 for m odd and even")
good = True
for m in range(2, 31):
    N = 2 * m
    p1, p2 = dds_kn(N)
    bip = [[e for e in p if (e[0] + e[1]) % 2 == 1] for p in (p1, p2)]
    good &= crossings(bip) == d(m) ** 2 and crossings((p1, p2)) == Z(N)
check(good, "the parity drawing is the DDS drawing of K_(2m) restricted to the even-odd edges: d_m^2 of its Z(2m) crossings, 2 <= m <= 30")
print("   e.g. Z(2m) vs d_m^2:", [(2 * m, Z(2 * m), d(m) ** 2) for m in range(3, 9)])

print("3. class-coloring minima (exhaustive, small cases; THM-922 re-check)")


def classcolor_min_kn(n):
    best = None
    for mask in range(1 << (n - 1)):
        cls = [(mask >> s) & 1 for s in range(n - 1)] + [0]
        p = ([], [])
        for i, j in itertools.combinations(range(n), 2):
            p[cls[(i + j) % n]].append((i, j))
        c = crossings(p)
        best = c if best is None else min(best, c)
    return best


check(all(classcolor_min_kn(n) == Z(n) for n in range(5, 12)), "min over all 2-page class colorings of K_n = Z(n), 5 <= n <= 11 (with Abrego et al. 2012: every n)")

print("4. the rank inequality of #165 (sanity)")
rng = random.Random(1)
good = True
for trial in range(300):
    m = rng.randint(2, 5)
    N = rng.randint(1, 4)
    # random subspaces of Q^N given by random integer bases of random dimension; keep pairs with U_i n V_i = 0
    def rand_sub():
        r = rng.randint(0, N)
        return sp.Matrix(N, r, lambda i, j: rng.randint(-1, 1)) if r else sp.zeros(N, 0)
    Us, Vs = [], []
    while len(Us) < m:
        U, V = rand_sub(), rand_sub()
        if U.shape[1] and V.shape[1] and sp.Matrix.hstack(U, V).rank() < U.rank() + V.rank():
            continue
        Us.append(U); Vs.append(V)
    tot = 0
    for U in Us:
        for V in Vs:
            if U.shape[1] == 0 or V.shape[1] == 0:
                continue
            tot += U.rank() + V.rank() - sp.Matrix.hstack(U, V).rank()
    good &= tot <= (m * m // 4) * N
check(good, "300 random instances (m <= 5, N <= 4): sum_(i,j) dim(U_i n V_j) <= floor(m^2/4) N whenever U_i n V_i = 0")
print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED")
