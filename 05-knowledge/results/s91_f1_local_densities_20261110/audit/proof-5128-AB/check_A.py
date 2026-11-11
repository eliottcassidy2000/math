#!/usr/bin/env python3
"""Audit proof-5128-AB: exact check of THM-5128 (A) (class law of a bilinear gluing).

For primes p, ranks r, random nondegenerate bilinear forms B(u, v) = u^T M v on (Z/p)^r and random probability
distributions D1, D2 with rational weights, everything below is computed in exact rational arithmetic:

  P(c)      = P(B(u, v) = c), by enumeration of all pairs (u, v);
  Rhat(c)  := P(c) - 1/p - (pi1 + pi2 - p^-r)([c = 0] - 1/p)                      (definition through (1));
  (a) sum_c Rhat(c) = 0;
  (b) Rhat(c) = sum_{u != 0} delta1(u) (P_v(B(u, v) = c) - 1/p)                   (the proof's R, after Fourier inversion);
  (c) Rhat(c) = sum_{u, v} (D1 - U)(u) (D2 - U)(v) [B(u, v) = c]                  (bilinear in the two deviations);
  (d) Rhat(c)^2 <= ((p-1)/p)^2 p^r |D1 - U|^2 |D2 - U|^2                          (bound (2), squared: exact);
  (e) Rhat = 0 if D1 or D2 is uniform; both uniform: P(c) = 1/p + p^-r([c = 0] - 1/p);
  (f) point mass at u0 != 0: P(c) = D2(affine hyperplane).
In floating point (complex exponentials):
  (g) the Fourier formulas of the proof: Phi(t) = pi1 + pi2 - p^-r + R(t), Rhat(c) = (1/p) sum_{t != 0} e(-tc/p) R(t),
      Parseval sum_{xi != 0} |Dhat(xi)|^2 = p^r |D - U|^2, |R(t)| <= p^(r/2) |D1 - U| |D2 - U|;
  (h) |D - U|_2 <= E and the second form of (2); point mass: |P(c) - 1/p| <= E_2.
Also: the largest value of |Rhat(c)| / (|D1 - U|_2 |D2 - U|_2) met, against the stated constant, and the exact
operator norm of the kernel on the sum-zero subspace (the sharp constant of an inequality of the shape (2)).
A degenerate form is tested to see that nondegeneracy is needed.
"""
import itertools
import math
import random
import sys
from fractions import Fraction

import numpy as np


def vectors(p, r):
    return list(itertools.product(range(p), repeat=r))


def det_mod(M, p):
    M = [row[:] for row in M]
    n = len(M)
    det = 1
    for i in range(n):
        piv = None
        for j in range(i, n):
            if M[j][i] % p:
                piv = j
                break
        if piv is None:
            return 0
        if piv != i:
            M[i], M[piv] = M[piv], M[i]
            det = -det
        det = det * M[i][i] % p
        inv = pow(M[i][i], -1, p)
        for j in range(i + 1, n):
            f = M[j][i] * inv % p
            for k in range(i, n):
                M[j][k] = (M[j][k] - f * M[i][k]) % p
    return det % p


def random_invertible(p, r, rng):
    while True:
        M = [[rng.randrange(p) for _ in range(r)] for _ in range(r)]
        if det_mod(M, p):
            return M


def bil(M, u, v, p):
    r = len(u)
    return sum(u[i] * M[i][j] * v[j] for i in range(r) for j in range(r)) % p


def random_distribution(V, rng, kind):
    n = len(V)
    if kind == "generic":
        w = [rng.randrange(0, 50) for _ in range(n)]
    elif kind == "sparse":
        w = [0] * n
        for _ in range(rng.randrange(1, 5)):
            w[rng.randrange(n)] += rng.randrange(1, 9)
    elif kind == "near_uniform":
        w = [1000 + rng.randrange(-30, 31) for _ in range(n)]
    elif kind == "heavy_zero":
        w = [rng.randrange(0, 5) for _ in range(n)]
        w[0] += rng.randrange(20, 200)
    elif kind == "no_zero":
        w = [rng.randrange(1, 20) for _ in range(n)]
        w[0] = 0
    elif kind == "uniform":
        w = [1] * n
    elif kind == "point0":
        w = [0] * n
        w[0] = 1
    elif kind == "point":
        w = [0] * n
        w[rng.randrange(1, n)] = 1
    else:
        raise ValueError(kind)
    if sum(w) == 0:
        w[rng.randrange(n)] = 1
    s = sum(w)
    return [Fraction(a, s) for a in w]


KINDS = ["generic", "sparse", "near_uniform", "heavy_zero", "no_zero", "uniform", "point0", "point"]


def check_instance(p, r, M, D1, D2, V, stats):
    n = len(V)
    idx0 = 0  # V[0] is the zero vector
    assert all(a == 0 for a in V[0])
    U = Fraction(1, n)
    pi1, pi2 = D1[idx0], D2[idx0]
    # table of B(u, v)
    Btab = [[bil(M, u, v, p) for v in V] for u in V]
    # P_v(B(u, v) = c) for each u
    cond = [[Fraction(0)] * p for _ in range(n)]
    for i in range(n):
        row = Btab[i]
        ci = cond[i]
        for j in range(n):
            if D2[j]:
                ci[row[j]] += D2[j]
    P = [sum(D1[i] * cond[i][c] for i in range(n)) for c in range(p)]
    assert sum(P) == 1
    main = [Fraction(1, p) + (pi1 + pi2 - U) * ((1 if c == 0 else 0) - Fraction(1, p)) for c in range(p)]
    Rhat = [P[c] - main[c] for c in range(p)]
    # (a)
    assert sum(Rhat) == 0, "sum Rhat"
    # (b)
    d1 = [D1[i] - U for i in range(n)]
    d2 = [D2[i] - U for i in range(n)]
    for c in range(p):
        alt = sum(d1[i] * (cond[i][c] - Fraction(1, p)) for i in range(1, n))
        assert alt == Rhat[c], "formula (b)"
    # (c)
    for c in range(p):
        alt = Fraction(0)
        for i in range(n):
            if d1[i] == 0:
                continue
            row = Btab[i]
            alt += d1[i] * sum(d2[j] for j in range(n) if row[j] == c)
        assert alt == Rhat[c], "formula (c)"
    # (d)
    n1 = sum(a * a for a in d1)
    n2 = sum(a * a for a in d2)
    bound_sq = Fraction(p - 1, p) ** 2 * (p ** r) * n1 * n2
    for c in range(p):
        assert Rhat[c] ** 2 <= bound_sq, "bound (2)"
    if n1 and n2:
        ratio = max(abs(float(x)) for x in Rhat) / math.sqrt(float(n1 * n2))
        stats["max_ratio"] = max(stats.get("max_ratio", 0.0), ratio)
        stats["const"] = (p - 1) / p * p ** (r / 2)
    # (e)
    if all(a == U for a in D1) or all(a == U for a in D2):
        assert all(x == 0 for x in Rhat), "uniform => Rhat = 0"
    if all(a == U for a in D1) and all(a == U for a in D2):
        for c in range(p):
            assert P[c] == Fraction(1, p) + U * ((1 if c == 0 else 0) - Fraction(1, p))
    # floating point Fourier formulas
    Vn = np.array(V, dtype=np.int64)
    Mn = np.array(M, dtype=np.int64)
    D1f = np.array([float(a) for a in D1])
    D2f = np.array([float(a) for a in D2])
    dots = (Vn @ Vn.T) % p                                # dots[xi, v] = xi . v
    hat2 = np.exp(2j * np.pi * dots / p) @ D2f            # D2hat(xi) = sum_v D2(v) e(xi.v/p)
    hat1 = np.exp(2j * np.pi * dots / p) @ D1f
    index = {v: i for i, v in enumerate(V)}
    Bstar = (Vn @ Mn) % p                                 # row u: (M^T u)^T = u^T M
    n1f, n2f = float(n1), float(n2)
    # Parseval
    assert abs(np.sum(np.abs(hat2[1:]) ** 2) - p ** r * n2f) < 1e-9
    assert abs(np.sum(np.abs(hat1[1:]) ** 2) - p ** r * n1f) < 1e-9
    Rt = np.zeros(p, dtype=np.complex128)
    for t in range(1, p):
        tb = (t * Bstar) % p
        ids = [index[tuple(int(a) for a in row)] for row in tb]
        Phi = sum(D1f[i] * hat2[ids[i]] for i in range(n))
        Rt[t] = sum((D1f[i] - 1.0 / n) * hat2[ids[i]] for i in range(1, n))
        assert abs(Phi - (float(pi1 + pi2 - U) + Rt[t])) < 1e-9, "Phi(t)"
        assert abs(Rt[t]) <= p ** (r / 2) * math.sqrt(n1f * n2f) + 1e-9, "|R(t)|"
    for c in range(p):
        val = sum(np.exp(-2j * np.pi * t * c / p) * Rt[t] for t in range(1, p)) / p
        assert abs(val - float(Rhat[c])) < 1e-9, "Rhat Fourier"
    E1 = float(np.max(np.abs(hat1[1:])))
    E2 = float(np.max(np.abs(hat2[1:])))
    assert n1f <= E1 ** 2 + 1e-12 and n2f <= E2 ** 2 + 1e-12, "|D - U|_2 <= E"
    for c in range(p):
        assert abs(float(Rhat[c])) <= (p - 1) / p * p ** (r / 2) * E1 * E2 + 1e-9, "second form of (2)"
    # (f) point mass at u0 != 0
    supp = [i for i in range(n) if D1[i]]
    if len(supp) == 1 and supp[0] != 0:
        i = supp[0]
        for c in range(p):
            hyper = sum(D2[j] for j in range(n) if Btab[i][j] == c)
            assert hyper == P[c]
            assert abs(float(P[c]) - 1.0 / p) <= E2 + 1e-12, "point mass, 1/p up to E_2"
            assert abs(float(P[c]) - 1.0 / p) <= (p - 1) / p * E2 + 1e-12
    stats["instances"] = stats.get("instances", 0) + 1


def sharp_constant(p, r, M, V):
    """Largest singular value of N_c(u, v) = [B(u, v) = c] on the sum-zero subspaces, for each c."""
    n = len(V)
    Q = np.eye(n) - np.ones((n, n)) / n
    out = []
    Bt = np.array([[bil(M, u, v, p) for v in V] for u in V])
    for c in range(p):
        N = (Bt == c).astype(np.float64)
        s = np.linalg.svd(Q @ N @ Q, compute_uv=False)
        out.append(float(s[0]))
    return out


def main():
    rng = random.Random(20261111)
    cases = [(2, 1, 60), (2, 2, 60), (2, 3, 40), (3, 1, 120), (3, 2, 120), (3, 3, 30), (5, 1, 120), (5, 2, 80),
             (7, 1, 120), (7, 2, 40), (11, 1, 40), (5, 3, 4)]
    print("exact check of THM-5128 (A): identity (1), sum Rhat = 0, bound (2), uniform and point-mass cases")
    print(f"{'p':>3} {'r':>2} {'instances':>9} {'max |Rhat|/(|d1||d2|)':>22} {'stated const':>13} "
          f"{'sharp const (c=0 | c!=0)':>26}")
    for p, r, trials in cases:
        V = vectors(p, r)
        stats = {}
        for trial in range(trials):
            if p > 2 and r == 2 and trial % 3 == 0:
                s = rng.randrange(0, 6)
                M = [[1, 0], [0, pow(2, s, p)]]            # the F1 form a1 a2 + 2^s b1 b2
            else:
                M = random_invertible(p, r, rng)
            k1 = KINDS[trial % len(KINDS)]
            k2 = KINDS[(trial // len(KINDS) + trial) % len(KINDS)] if trial % 5 else rng.choice(KINDS)
            D1 = random_distribution(V, rng, k1)
            D2 = random_distribution(V, rng, k2)
            check_instance(p, r, M, D1, D2, V, stats)
        M = random_invertible(p, r, rng)
        sc = sharp_constant(p, r, M, V)
        nz = max(sc[1:]) if p > 1 else 0.0
        print(f"{p:3d} {r:2d} {stats['instances']:9d} {stats.get('max_ratio', 0.0):22.6f} "
              f"{stats.get('const', (p - 1) / p * p ** (r / 2)):13.6f}        {sc[0]:9.6f} | {nz:9.6f}")
    # nondegeneracy is needed
    p, r = 3, 2
    V = vectors(p, r)
    M = [[1, 0], [0, 0]]
    Uv = [Fraction(1, len(V))] * len(V)
    P0 = sum(Uv[i] * Uv[j] for i, u in enumerate(V) for j, v in enumerate(V) if bil(M, u, v, p) == 0)
    law0 = Fraction(1, p) + Fraction(1, p ** r) * (1 - Fraction(1, p))
    print(f"degenerate form a1 a2 on (Z/3)^2, both uniform: P(B = 0) = {P0}, formula (1) would give {law0}: "
          f"{'differs (nondegeneracy needed)' if P0 != law0 else 'equal'}")
    assert P0 != law0
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
