#!/usr/bin/env python3
"""Eigenvalue-1 points of the excess transfer operator of THM-4658 4(b) on step densities of level Q (argument principle).

On densities f constant on residue classes modulo 3^Q (beta = f dHaar), statement 4(b) acts by
    (P f)(a) = -2^(-s) f(a/2) - (1/3)(3/2)^s f((3a+1)/2) + (3/2)^s <f, 1_(2+3Z)>,   a in Z/3^Q,
which preserves the mean-zero subspace (dimension 3^Q - 1); the step densities of level Q form an invariant subspace,
so an eigenvalue-1 point rho of this matrix is an exponent at which a power-law excess main term beta(A) N^rho with a
level-Q step density is consistent with the exact recursion.  The script counts the zeros of
det(P_Q(s) - I) in rectangles [SIG0, 1.5] x [T0, T1] by the argument principle (winding number of the determinant
along the boundary) for Q = 1..4, and checks the exact spectrum of THM-4658 statement 7 (eigenvalues -2^(-s) zeta,
zeta^(2 3^k) = 1, k < Q) at sample points.  NUMERICAL control of a PROVED statement; the winding counts are integers
stable under refinement of the boundary sampling, which is checked.
Session opus-2026-10-15-S37.  Usage: python collatz_excess_operator_resonance_scan_20261015.py out.json
"""
import sys, json
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else 'resonance_scan.json'


def build(Q):
    m = 3 ** Q
    i2 = pow(2, -1, m)
    U = np.zeros((m, m))
    V = np.zeros((m, m))
    K = np.zeros((m, m))
    for a in range(m):
        U[a, (a * i2) % m] = 1
        V[a, ((3 * a + 1) * i2) % m] = 1
        for c in range(2, m, 3):
            K[a, c] = 1.0 / m
    B = np.zeros((m, m - 1))
    for a in range(m - 1):
        B[a, a] = 1
        B[m - 1, a] = -1
    Bp = np.linalg.pinv(B)
    return Bp @ U @ B, Bp @ V @ B, Bp @ K @ B


def det_minus_one(mats, s):
    Ub, Vb, Kb = mats
    A = -(2.0 ** (-s)) * Ub - (1.5 ** s) / 3 * Vb + (1.5 ** s) * Kb
    sign, logdet = np.linalg.slogdet(A - np.eye(A.shape[0]))
    return sign, logdet


def winding(mats, path):
    ph = []
    for s in path:
        sg, ld = det_minus_one(mats, s)
        ph.append(np.angle(sg))
    ph = np.unwrap(np.array(ph))
    return (ph[-1] - ph[0]) / (2 * np.pi)


def rect(s0, s1, t0, t1, n):
    a = np.linspace(s0, s1, n) + 1j * t0
    b = s1 + 1j * np.linspace(t0, t1, n)
    c = np.linspace(s1, s0, n) + 1j * t1
    d = s0 + 1j * np.linspace(t1, t0, n)
    return np.concatenate([a, b, c, d])


res = {}
# exact spectrum check of THM-4658 statement 7 at sample points: eigenvalues of P~_Q(s) on mean-zero step densities
# are -2^(-s) zeta with zeta^(2 3^k) = 1, k = 0..Q-1 (with multiplicity: one full set of 2 3^(Q-j-1)-th roots per j)
spec_err = {}
for Q in range(1, 5):
    mats = build(Q)
    worst = 0.0
    for s in (0.7 + 3.1j, 0.3 - 11.0j, 1.2 + 0.4j, 0.05 + 40.0j):
        Ub, Vb, Kb = mats
        Mx = -(2.0 ** (-s)) * Ub - (1.5 ** s) / 3 * Vb + (1.5 ** s) * Kb
        ev = np.sort_complex(np.linalg.eigvals(Mx))
        pred = []
        for j in range(Q):
            L = 2 * 3 ** (Q - j - 1)
            pred += [-(2.0 ** (-s)) * np.exp(2j * np.pi * r / L) for r in range(L)]
        pred = np.sort_complex(np.array(pred))
        # match each predicted eigenvalue to the nearest computed one
        used = np.zeros(len(ev), dtype=bool)
        for z in pred:
            d = np.abs(ev - z) + used * 1e9
            i = int(np.argmin(d))
            used[i] = True
            worst = max(worst, float(d[i]))
    spec_err[Q] = worst
    print(f'Q = {Q}: max distance between computed and predicted spectrum (statement 7): {worst:.2e}', flush=True)
res['statement7_spectrum_max_error'] = spec_err
for Q in range(1, 5):
    mats = build(Q)
    rows = []
    for (s0, t0, t1) in [(0.05, 0.01, 50.0), (0.05, 50.0, 150.0), (0.05, 150.0, 300.0), (0.5, 0.01, 300.0)]:
        w1 = winding(mats, rect(s0, 1.5, t0, t1, 4000))
        w2 = winding(mats, rect(s0, 1.5, t0, t1, 12000))
        rows.append(dict(sigma_min=s0, t_range=[t0, t1], zeros=round(w2), winding_coarse=float(w1), winding_fine=float(w2)))
        print(f'Q = {Q}: zeros of det(P_Q(s) - I) in [{s0}, 1.5] x [{t0}, {t1}]: {w2:.3f} (coarse {w1:.3f})', flush=True)
    res[Q] = rows
json.dump(res, open(OUT, 'w'), indent=1)
