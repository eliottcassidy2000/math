#!/usr/bin/env python3
"""S37 audit 2, script b1b: the limitation of THM-4658 statement 7 is essential.

The excess fixed-point equation of statement 4(b),
    beta(A) = -2^(-rho) beta(A/2) - (3/2)^rho [beta((3A+1)/2) - Haar(A) beta(2 + 3Z_3)],
relates, on a level-q cylinder A, two level-q values of beta to one level-(q+1) value.  Upward construction: given the
level-q masses, the level-q equations fix the level-(q+1) masses of the classes = 2 (mod 3) (each is (3a+1)/2 for exactly one
a mod 3^q); their sums over siblings reproduce the parent masses because of the level-(q-1) equations; the children of
classes = 0, 1 (mod 3) are free.  So for EVERY rho and every mean-zero level-1 marginal there are excess eigen-distributions
(eigenvalue 1) of unbounded depth.  The script builds such a distribution to depth 7 with random free splits and a
prescribed level-1 marginal, then checks every cylinder equation of levels 0..6 directly in the measure form (images of
cylinders by brute force), for several rho including rho = 1/2 and a zero of D.  It also shows the depth-7 masses are not
uniformly split inside level-Q classes for any Q < 7 (so they are not finite-depth densities, as statement 7 requires).
"""
import cmath
import numpy as np

rng = np.random.default_rng(3)
DEPTH = 7


def build(rho, beta1):
    a2 = 2.0 ** (-rho)
    b = 1.5 ** rho
    levels = [np.array([sum(beta1)], dtype=complex), np.array(beta1, dtype=complex)]   # level 0: beta(Z_3) = 0
    for q in range(1, DEPTH):
        m, M = 3 ** q, 3 ** (q + 1)
        bq = levels[q]
        nxt = np.full(M, np.nan + 0j)
        i2m = pow(2, -1, m); i2M = pow(2, -1, M)
        beta_2mod3 = levels[1][2]
        for a in range(m):
            c = ((3 * a + 1) * i2M) % M
            # equation at level q for A = a + 3^q Z_3:  bq[a] + a2 bq[a/2] + b (nxt[c] - 3^-q beta(2+3Z)) = 0
            nxt[c] = -(bq[a] + a2 * bq[(a * i2m) % m]) / b + beta_2mod3 / m
        for p in range(m):
            kids = [p + m * k for k in range(3)]
            if p % 3 == 2:
                assert all(not np.isnan(nxt[k].real) for k in kids)
            else:
                w = rng.normal(size=3) + 1j * rng.normal(size=3)
                w = w - w.mean()                       # random split with the right total
                for j, k in enumerate(kids):
                    nxt[k] = bq[p] / 3 + w[j] * abs(bq[p] + 1)
        levels.append(nxt)
    return levels


def check(rho, levels):
    a2 = 2.0 ** (-rho)
    b = 1.5 ** rho
    worst_eq = 0.0
    worst_marg = 0.0
    for q in range(0, DEPTH):
        m = 3 ** q
        fine = 3 ** (q + 2)
        bq = levels[q]
        for a in range(m):
            reps = [a + m * k for k in range(fine // m)]
            half = sorted(set((x * pow(2, -1, fine)) % m for x in reps))
            odd = sorted(set(((3 * x + 1) * pow(2, -1, 3 * fine)) % (3 * m) for x in reps))
            assert len(half) == 1 and len(odd) == 1
            lhs = bq[a]
            rhs = -a2 * bq[half[0]] - b * (levels[q + 1][odd[0]] - levels[1][2] / m)
            worst_eq = max(worst_eq, abs(lhs - rhs))
        if q + 1 <= DEPTH:
            kids = levels[q + 1].reshape(3, m).sum(axis=0) if q + 1 <= DEPTH else None
            worst_marg = max(worst_marg, np.abs(kids - bq).max())
    return worst_eq, worst_marg


beta1 = np.array([0.3 + 0.1j, -0.7 + 0.2j, 0.4 - 0.3j])      # mean zero, nonzero
for rho in (0.5, 0.75 + 5.0j, 1.2 - 2.0j, 0.99820922076525681 + 240.24783014069124j, 0.3 + 0.0j):
    L = build(rho, beta1)
    eq, mg = check(rho, L)
    # departure from uniform splitting inside level-Q classes (finite-depth densities have zero departure)
    dep = []
    for Q in range(1, DEPTH):
        m = 3 ** Q
        bD = L[DEPTH].reshape(3 ** (DEPTH - Q), m)          # rows: k, columns: class mod 3^Q
        dep.append(float(np.abs(bD - bD.mean(axis=0)).max()))
    print(f'rho = {rho}: max |beta(A) - (P~ beta)(A)| over all cylinders of levels 0..{DEPTH - 1} = {eq:.1e}; marginal '
          f'consistency {mg:.1e}; level-1 marginal {np.round(L[1], 3)}; max departure from uniform splitting inside level-Q '
          f'classes, Q = 1..{DEPTH - 1}: ' + ', '.join(f'{d:.2f}' for d in dep))
