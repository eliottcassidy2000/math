#!/usr/bin/env python3
"""Multiplier-only block certificates: fully adversarial hidden digits, identity steps deleted.

Model.  The coupling multiplier evolves deterministically, a' = a * mbar_(a j + b) / mbar_j (mod d); the translation
part b' of the next coupling depends on hidden digits and on the constants r_i, so let an adversary choose it, subject
only to the next coupling being non-identity (identity steps do not move the debt and are deleted: a block is k
non-identity steps, a stopping time; identity runs are finite).  For a form Q and a unit vector u, the badness of a block
second moment S, f_Q(S) = lambda_max(Q^-1/2 S Q^-1/2) - tr(Q^-1 S)/2 = max_u <S, P_u>, P_u = Q^-1/2 (u u^T - I/2) Q^-1/2,
is linear in S for fixed u, so the adversary's best response decomposes over branches:
    g_s(a, b; u) = <C(a, b), P_u> + (1/d) sum_j max_{b' : (a'_j, b') != (1, 0)} g_(s-1)(a'_j, b'; u),   g_0 = 0.
If max_u g_k(a, b; u) < 0 for every non-identity coupling (a, b), then every possible k-block (whatever the constants r_i
and the hidden digits) is balanced by Q, and THM-4611's argument (with the time change) gives transience.
The maximum over the unit sphere is bounded rigorously by a Fibonacci grid plus a Lipschitz term:
|g(u) - g(u')| <= 2 L |u - u'| with L = k max_pi lambda_max(Q^-1/2 C_pi Q^-1/2).
Usage: python3 adversarial_cert.py"""
import math, sys
import numpy as np
from block_balance import int_coords, ratio_group, coupling_D

def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); th = math.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)

def covering_radius(pts, probe=200000, seed=0):
    """Estimate (upper bound by probing) of the covering radius of a point set on S^2 (chordal)."""
    rng = np.random.default_rng(seed); X = rng.normal(size=(probe, 3)); X /= np.linalg.norm(X, axis=1)[:, None]
    best = np.zeros(probe)
    for s in range(0, len(pts), 2000):
        best = np.maximum(best, X @ pts[s:s + 2000].T.max(axis=1) if False else np.maximum(best, (X @ pts[s:s + 2000].T).max(axis=1)))
    return float(np.sqrt(np.maximum(0, 2 - 2 * best)).max())

def adversarial_g(d, m, Q, k, U):
    rho, v = int_coords(m)
    G = ratio_group(m, d); mb = [x % d for x in m]
    w, V = np.linalg.eigh(np.array(Q, dtype=float)); Qmh = V @ np.diag(w ** -0.5) @ V.T
    states = [(a, b) for a in G for b in range(d)]
    idx = {s: n for n, s in enumerate(states)}
    # per-state linear weight <C, P_u> for all u: u^T Ct u - tr(Ct)/2 with Ct = Qmh C Qmh, C = D/d
    W = np.zeros((len(states), len(U))); lmax = 0.0
    for (a, b), n in idx.items():
        C = np.array(coupling_D(v, a, b, d, rho), dtype=float) / d
        Ct = Qmh @ C @ Qmh
        W[n] = np.einsum('ni,ij,nj->n', U, Ct, U) - np.trace(Ct) / 2
        lmax = max(lmax, float(np.linalg.eigvalsh(Ct)[-1]))
    nxt = {}
    for (a, b) in states:
        nxt[(a, b)] = [(a * mb[(a * j + b) % d] * pow(mb[j], -1, d)) % d for j in range(d)]
    g = np.zeros((len(states), len(U)))
    for s in range(1, k + 1):
        gn = np.empty_like(g)
        best_by_a = {}
        for a in G:
            rows = [idx[(a, bb)] for bb in range(d) if (a, bb) != (1, 0)]
            best_by_a[a] = g[rows].max(axis=0)
        for (a, b), n in idx.items():
            gn[n] = W[n] + sum(best_by_a[ap] for ap in nxt[(a, b)]) / d
        g = gn
    return states, g, lmax

if __name__ == '__main__':
    U = fib_sphere(200000)
    cr = covering_radius(U)
    tests = [
        (7, [1, 1, 1, 1, 2, 3, 5], [[6, -1, -1], [-1, 5, -1], [-1, -1, 6]]),
        (5, [1, 2, 3, 7, 1], [[20, -4, -8], [-4, 15, -3], [-8, -3, 20]]),
        (5, [1, 1, 2, 3, 7], [[8, -2, -3], [-2, 6, -1], [-3, -1, 8]]),
    ]
    print(f"grid: {len(U)} Fibonacci points, probed covering radius (chordal) {cr:.5f}", flush=True)
    for d, m, Q in tests:
        for k in range(1, 9):
            states, g, lmax = adversarial_g(d, m, Q, k, U)
            worst = []
            for n, (a, b) in enumerate(states):
                if (a, b) == (1, 0): continue
                worst.append((g[n].max(), (a, b)))
            mx, st = max(worst)
            lip = 2 * k * lmax * cr
            print(f"Z_{d} m = {m} Q = {Q} k = {k}: max over non-identity couplings and grid of g_k = {mx:+.5f} at {st}; "
                  f"Lipschitz slack {lip:.5f} -> {'CERTIFIED' if mx + lip < 0 else 'not certified'}", flush=True)
