#!/usr/bin/env python3
"""
u21_core_embed_lm_monad.py  (monad-explorer, 2026-06-07)

Robust Levenberg-Marquardt unit-distance embedder for the 5 stored
u(21)=57 graph6 cores.  Replaces the diverging gradient descent of
u21_core_embed_monad.py (which blew up to NaN -> spurious "residual 0").

Goal: decide REALIZABILITY of each core as a planar unit-distance graph,
and produce coordinates with all 57 edges at distance ~1 so that a faithful
embedding (no extra accidental unit distances) certifies u(21) >= 57.

Residual per edge e=(i,j):  r_e = |p_i - p_j| - 1.
Minimize sum r_e^2 by damped Gauss-Newton (LM) with multiple random starts.
"""
import math
import random
import numpy as np

N21_GRAPH6 = [
    "Tsc@IGC@GD?R?S?Wd@A_CK@HG@VM??PRKOUZ",
    "TCKx?D?OI?OMCBSA_L?ApA_gEA\\EG?PBSCPV",
    "TCTWACAG@@CDKC?e?QgQA@OMOq]F??OUcCEj",
    "TCS`?H??XIZ?K_Co`CG@JO[?EOSCpGOTSCE\\",
    "T??_`OhSCSYA@I?c?OyWBEa@c?SIU?Aa[?el",
]


def graph6_decode(code):
    data = [ord(ch) - 63 for ch in code.strip()]
    n = data[0]
    bits = []
    for byte in data[1:]:
        for k in range(5, -1, -1):
            bits.append((byte >> k) & 1)
    adj = [[0] * n for _ in range(n)]
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if idx < len(bits) and bits[idx]:
                adj[i][j] = adj[j][i] = 1
            idx += 1
    return adj


def edges_of(adj):
    n = len(adj)
    return [(i, j) for i in range(n) for j in range(i + 1, n) if adj[i][j]]


def lm_embed(n, edges, seed, iters=200):
    """Damped Gauss-Newton on edge residuals r_e=|p_i-p_j|-1.
    Returns (P (n,2) array, max_edge_residual)."""
    rng = np.random.default_rng(seed)
    P = rng.uniform(-math.sqrt(n) / 2, math.sqrt(n) / 2, size=(n, 2))
    m = len(edges)
    lam = 1e-2
    prev_cost = None
    for it in range(iters):
        # residuals + Jacobian (m x 2n)
        J = np.zeros((m, 2 * n))
        r = np.zeros(m)
        for e, (i, j) in enumerate(edges):
            d = P[i] - P[j]
            dist = math.hypot(d[0], d[1])
            if dist < 1e-12:
                dist = 1e-12
            r[e] = dist - 1.0
            gx, gy = d[0] / dist, d[1] / dist
            J[e, 2 * i] = gx
            J[e, 2 * i + 1] = gy
            J[e, 2 * j] = -gx
            J[e, 2 * j + 1] = -gy
        cost = float(r @ r)
        # LM step: (J^T J + lam I) dx = -J^T r
        JT = J.T
        H = JT @ J
        g = JT @ r
        A = H + lam * np.eye(2 * n)
        try:
            dx = np.linalg.solve(A, -g)
        except np.linalg.LinAlgError:
            lam *= 10
            continue
        Pnew = P + dx.reshape(n, 2)
        # evaluate new cost
        rn = np.array([math.hypot(*(Pnew[i] - Pnew[j])) - 1.0 for (i, j) in edges])
        newcost = float(rn @ rn)
        if newcost < cost:
            P = Pnew
            lam = max(lam * 0.5, 1e-9)
            if prev_cost is not None and abs(prev_cost - newcost) < 1e-16:
                break
            prev_cost = newcost
            if newcost < 1e-20:
                break
        else:
            lam *= 3.0
            if lam > 1e8:
                break
    maxres = max(abs(math.hypot(*(P[i] - P[j])) - 1.0) for (i, j) in edges)
    return P, maxres


def count_unit(P, tol):
    n = len(P)
    c, pairs = 0, []
    for i in range(n):
        for j in range(i + 1, n):
            d = math.hypot(*(P[i] - P[j]))
            if abs(d - 1.0) < tol:
                c += 1
                pairs.append((i, j))
    return c, pairs


def min_pair_dist(P):
    n = len(P)
    mn = float("inf")
    for i in range(n):
        for j in range(i + 1, n):
            mn = min(mn, math.hypot(*(P[i] - P[j])))
    return mn


def main():
    print("=" * 72)
    print("u(21)=57 cores: robust LM unit-distance embedding (monad-explorer)")
    print("=" * 72)
    for idx, code in enumerate(N21_GRAPH6, 1):
        adj = graph6_decode(code)
        n = len(adj)
        E = edges_of(adj)
        print(f"\n--- Core {idx}: n={n}, |E|={len(E)} ---")
        best = None
        for s in range(40):
            P, res = lm_embed(n, E, seed=7919 * idx + 31 * s, iters=160)
            if best is None or res < best[1]:
                best = (P, res)
            if res < 1e-9:
                break
        P, res = best
        print(f"  best edge residual over starts: {res:.3e}")
        if res < 1e-6:
            for tol in (1e-6, 1e-4):
                c, pairs = count_unit(P, tol)
                print(f"  total unit distances (tol={tol:.0e}): {c}")
            extra = c - len(E)
            mind = min_pair_dist(P)
            print(f"  min pairwise distance: {mind:.5f}")
            print(f"  FAITHFUL embedding: {'YES (UD graph == core, 57 units)' if c == len(E) else f'NO -- {extra} accidental extra unit pairs'}")
            print(f"  ==> REALIZABLE: u(21) >= {c}")
            # save coordinates of the first realizable core
            np.save(f"05-knowledge/results/u21_core{idx}_coords.npy", P)
            print(f"  coords saved to 05-knowledge/results/u21_core{idx}_coords.npy")
        else:
            print("  did NOT converge to a unit-distance embedding from these starts")
            print("  (inconclusive: either non-realizable or needs more starts)")


if __name__ == "__main__":
    main()
