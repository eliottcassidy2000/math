#!/usr/bin/env python3
"""Faithful embeddings of EI17 in the finite unit-distance graphs G_p = Cay(F_(p^2), mu_(p+1)) (S39; twisted model).

Same direction structure as ei17_finite_field_realizations_20261017.py (gauge U5 = 1; free U0, U1, U2, U3, U4, U6;
dependent U7..U12), now with every direction in the norm-one torus mu_(p+1) of F_(p^2) = F_p(sqrt n) and positions in
F_(p^2): a faithful embedding of EI17 in G_p up to the affine automorphisms fixing the gauge edge.  If p = 54 (mod 55)
the primes of L = Q(zeta55) above p have residue field F_(p^2) on which complex conjugation acts as Frobenius
(THM-4628), so a realization of EI17 with all directions in L reduces, at a prime of good reduction, to such an embedding.
The quadratic step: for c in F_(p^2), the units x with N(c - x) = 1 are x = (N(c) +- s sqrt n)/(2 conj(c)) when
D = N(c)(N(c) - 4) = n s^2 is a nonsquare in F_p, the single x = N(c)/(2 conj(c)) when D = 0, none when D is a nonzero
square (THM-4663, codegree lemma); c = 0 makes every unit a solution, and each such degenerate branch identifies two
vertices (U1 + U2 = 1: v5 = v11; U1 + U2 = 2: v7 = v10; U1 + U2 + U3 + U6 = 1: v2 = v11; U4 = U2 - 1: v1 = v6), so it is
skipped.
Usage: python ei17_twisted_unit_graph_embeddings_20261017.py p1,p2,...
Session opus-2026-10-17-S39.
"""
import sys
import numpy as np
from numba import njit

VERT = {1: (-1, 0, -1, 0, -1, 0, -1), 2: (0, 0, -1, -1, -1, 0, -1), 3: (0, 0, -1, -1, 0, 0, -1), 4: (-1, 1, -1, 0, -1, 0, -1),
        5: (0, 0, -1, -1, 0, 0, 0), 6: (0, 0, -1, -1, -1, 1, -1), 7: (-1, 0, 0, 0, -1, 0, -1), 8: (-1, 0, -1, 1, -1, 0, -1),
        9: (-1, 0, 0, 0, -1, 0, 0), 10: (1, 0, -1, -1, -1, 0, -1), 11: (-1, 0, 0, 0, 0, 0, 0), 12: (1, 0, -1, -1, 0, 0, -1),
        13: (-1, 1, -1, 1, -1, 0, -1), 14: (0, 0, 0, -1, 0, 0, 0), 15: (-1, 0, 0, 1, -1, 0, -1), 16: (1, 0, -1, -1, -1, 1, -1),
        17: (0, 0, 0, 0, 0, 0, 0)}
VM = np.array([VERT[v] for v in range(1, 18)], dtype=np.int64)
PE = [(1, 4), (1, 7), (1, 8), (2, 3), (2, 6), (2, 10), (3, 5), (3, 7), (3, 12), (4, 6), (4, 13), (4, 14), (5, 8),
      (5, 9), (5, 14), (6, 11), (6, 16), (7, 9), (7, 15), (8, 13), (8, 15), (9, 10), (9, 11), (10, 12), (10, 16),
      (11, 12), (11, 17), (13, 17), (14, 15), (14, 17), (16, 17)]
EDGES = np.array([(i - 1, j - 1) for i, j in PE], dtype=np.int64)


@njit(cache=True)
def pw(a, e, p):
    r = 1; a %= p
    while e > 0:
        if e & 1:
            r = r * a % p
        a = a * a % p; e >>= 1
    return r


@njit(cache=True)
def mul(a0, a1, b0, b1, n, p):
    return (a0 * b0 + n * (a1 * b1 % p)) % p, (a0 * b1 + a1 * b0) % p


@njit(cache=True)
def norm(a0, a1, n, p):
    return (a0 * a0 - n * (a1 * a1 % p)) % p


@njit(cache=True)
def units_on(c0, c1, n, p, sq, outr, outi):
    """units x with N(c - x) = 1; returns the count (p + 1 for the degenerate c = 0, then outr/outi unused)."""
    c0 %= p; c1 %= p
    if c0 == 0 and c1 == 0:
        return p + 1
    N = norm(c0, c1, n, p)
    if N == 0:
        return 0
    D = N * ((N - 4) % p) % p
    # 1/(2 conj c) = c / (2 N)
    i2N = pw(2 * N % p, p - 2, p)
    k0, k1 = c0 * i2N % p, c1 * i2N % p          # = 1/(2 conj(c))
    if D == 0:
        outr[0], outi[0] = N * k0 % p, N * k1 % p
        return 1
    if sq[D] >= 0:
        return 0
    s = sq[D * pw(n, p - 2, p) % p]
    cnt = 0
    for sg in (1, -1):
        r0, r1 = N, (sg * s) % p                    # N + sg s sqrt n
        x0, x1 = mul(r0, r1, k0, k1, n, p)
        outr[cnt], outi[cnt] = x0, x1
        cnt += 1
    return cnt


@njit(cache=True)
def count(p, n, sq, mu0, mu1, VM, EDGES):
    m = mu0.shape[0]
    faithful = 0; total = 0; edge_fail = 0; degenerate = 0
    r3 = np.zeros(2, dtype=np.int64); i3_ = np.zeros(2, dtype=np.int64)
    r6 = np.zeros(2, dtype=np.int64); i6_ = np.zeros(2, dtype=np.int64)
    r4 = np.zeros(2, dtype=np.int64); i4_ = np.zeros(2, dtype=np.int64)
    r0 = np.zeros(2, dtype=np.int64); i0_ = np.zeros(2, dtype=np.int64)
    zr = np.zeros(17, dtype=np.int64); zi = np.zeros(17, dtype=np.int64)
    U3r = np.zeros(p + 1, dtype=np.int64); U3i = np.zeros(p + 1, dtype=np.int64)
    for a in range(m):
        U1r, U1i = mu0[a], mu1[a]
        for b in range(m):
            U2r, U2i = mu0[b], mu1[b]
            n3 = units_on(U1r + U2r - 1, U1i + U2i, n, p, sq, r3, i3_)
            n6 = units_on(2 - U1r - U2r, -U1i - U2i, n, p, sq, r6, i6_)
            if n3 == p + 1 or n6 == p + 1:
                degenerate += 1        # c = 0 identifies two vertices (v5 = v11, resp. v7 = v10): not faithful
                continue
            for s3 in range(n3):
                for s6 in range(n6):
                    Ur = (1, 0, U1r, U2r, r3[s3], 0, r6[s6]); Ui = (0, 0, U1i, U2i, i3_[s3], 0, i6_[s6])
                    t0 = (2 * U2r - r3[s3] - r6[s6] - 1) % p; t1 = (2 * U2i - i3_[s3] - i6_[s6]) % p
                    if norm(t0, t1, n, p) != 1:
                        continue
                    n4 = units_on(U1r + U2r + r3[s3] + r6[s6] - 1, U1i + U2i + i3_[s3] + i6_[s6], n, p, sq, r4, i4_)
                    if n4 == p + 1:
                        degenerate += 1
                        continue
                    for s4 in range(n4):
                        n0 = units_on(1 - U2r + r4[s4], -U2i + i4_[s4], n, p, sq, r0, i0_)
                        if n0 == p + 1:
                            degenerate += 1
                            continue
                        for s0 in range(n0):
                            t0 = (-r0[s0] + U1r - U2r + r3[s3] + r6[s6] + 1) % p
                            t1 = (-i0_[s0] + U1i - U2i + i3_[s3] + i6_[s6]) % p
                            if norm(t0, t1, n, p) != 1:
                                continue
                            total += 1
                            vr = (1, r0[s0], U1r, U2r, r3[s3], r4[s4], r6[s6])
                            vi = (0, i0_[s0], U1i, U2i, i3_[s3], i4_[s4], i6_[s6])
                            for v in range(17):
                                s_r = 0; s_i = 0
                                for t in range(7):
                                    s_r += VM[v, t] * vr[t]; s_i += VM[v, t] * vi[t]
                                zr[v] = s_r % p; zi[v] = s_i % p
                            ok = True
                            for v in range(17):
                                for w in range(v + 1, 17):
                                    if zr[v] == zr[w] and zi[v] == zi[w]:
                                        ok = False
                            if ok:
                                for e in range(EDGES.shape[0]):
                                    x, y = EDGES[e, 0], EDGES[e, 1]
                                    if norm(zr[x] - zr[y], zi[x] - zi[y], n, p) != 1:
                                        edge_fail += 1
                                faithful += 1
    return faithful, total, degenerate, edge_fail


if __name__ == '__main__':
    plist = [int(x) for x in sys.argv[1].split(',')]
    print('  p  p mod 55  faithful embeddings  all solutions  degenerate branches skipped')
    for p in plist:
        n = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
        sq = -np.ones(p, dtype=np.int64)
        for x in range(p):
            if sq[x * x % p] < 0:
                sq[x * x % p] = x
        mu = [(a, b) for a in range(p) for b in range(p) if (a * a - n * b * b) % p == 1]
        assert len(mu) == p + 1
        mu0 = np.array([u[0] for u in mu], dtype=np.int64); mu1 = np.array([u[1] for u in mu], dtype=np.int64)
        f, tot, deg, ef = count(p, n, sq, mu0, mu1, VM, EDGES)
        assert ef == 0
        print('%5d  %3d     %6d            %6d        %6d' % (p, p % 55, f, tot, deg), flush=True)
