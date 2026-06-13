#!/usr/bin/env python3
"""
u21_embed_v4_coulomb_monad.py  (monad-explorer, 2026-06-07)

Realize the genuine AMP extremal u(21)=57 graphs via Coulomb-ANNEALED stress:
minimize  E(P) = sum_{edges}(|p_i-p_j|-1)^2  +  lam * sum_{i<j} 1/|p_i-p_j|
with lam annealed 0.5 -> 0.  The 1/d Coulomb term forbids the vertex-collapse
that defeated edges-only LM, is scale-aware (no hard DMIN that excluded the
true close-vertex realization in v3), and as lam->0 the edge term dominates so
we land ON a true unit-distance realization if the basin was found.

Self-tested on the same known-realizable graphs as v3, then applied to the 5
cores with many MDS/random starts.  A success = distinct points (min dist
>1e-3) with edge residual <1e-9 and >=57 unit distances.
"""
import math
import numpy as np
from collections import deque

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


def graph_dist(adj):
    n = len(adj)
    D = np.full((n, n), float(n))
    for s in range(n):
        dist = [-1] * n
        dist[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in range(n):
                if adj[u][v] and dist[v] < 0:
                    dist[v] = dist[u] + 1
                    q.append(v)
        for v in range(n):
            if dist[v] >= 0:
                D[s, v] = dist[v]
    return D


def mds_init(adj, rng, jitter=0.1):
    n = len(adj)
    D = graph_dist(adj)
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(B)
    order = np.argsort(w)[::-1]
    L = np.sqrt(np.maximum(w[order][:2], 0))
    P = V[:, order][:, :2] * L
    return P + rng.normal(scale=jitter, size=P.shape)


def coulomb_anneal(P, edges, n, rng):
    """Gradient descent on edge stress + annealed Coulomb, with adaptive step."""
    Ee = edges
    lam_schedule = [0.5, 0.3, 0.2, 0.12, 0.07, 0.04, 0.02, 0.01,
                    5e-3, 2e-3, 1e-3, 3e-4, 1e-4, 0.0]
    for lam in lam_schedule:
        step = 0.02
        for it in range(600):
            G = np.zeros((n, 2))
            # edge stress
            for (i, j) in Ee:
                d = P[i] - P[j]
                dist = math.hypot(d[0], d[1]) or 1e-12
                f = 2.0 * (dist - 1.0) / dist
                G[i] += f * d
                G[j] -= f * d
            # Coulomb repulsion (all pairs)
            if lam > 0:
                for i in range(n):
                    for j in range(i + 1, n):
                        d = P[i] - P[j]
                        r2 = d[0] * d[0] + d[1] * d[1]
                        if r2 < 1e-9:
                            r2 = 1e-9
                        # dE/dP_i for lam/ r  = -lam * d / r^3
                        coef = -lam / (r2 ** 1.5)
                        G[i] += coef * d
                        G[j] -= coef * d
            gn = np.linalg.norm(G)
            if gn < 1e-12:
                break
            P = P - step * G / (1.0 + gn)
    return P


def gauss_newton_polish(P, edges, n, steps=60):
    m = len(edges)
    lam = 1e-6
    for _ in range(steps):
        J = np.zeros((m, 2 * n))
        r = np.zeros(m)
        for e, (i, j) in enumerate(edges):
            d = P[i] - P[j]
            dist = math.hypot(d[0], d[1]) or 1e-15
            r[e] = dist - 1.0
            gx, gy = d[0] / dist, d[1] / dist
            J[e, 2 * i], J[e, 2 * i + 1] = gx, gy
            J[e, 2 * j], J[e, 2 * j + 1] = -gx, -gy
        H = J.T @ J + lam * np.eye(2 * n)
        g = J.T @ r
        for c in (0, 1, 3):  # gauge pin v0, v1.y
            H[c, :] = 0; H[:, c] = 0; H[c, c] = 1; g[c] = 0
        try:
            dx = np.linalg.solve(H, -g)
        except np.linalg.LinAlgError:
            break
        Pn = P + dx.reshape(n, 2)
        if max(abs(math.hypot(*(Pn[i] - Pn[j])) - 1.0) for (i, j) in edges) < \
           max(abs(math.hypot(*(P[i] - P[j])) - 1.0) for (i, j) in edges):
            P = Pn
            lam = max(lam * 0.5, 1e-12)
        else:
            lam *= 5
            if lam > 1e6:
                break
    return P


def analyze(P, adj):
    n = len(P)
    E = edges_of(adj)
    er = max(abs(math.hypot(*(P[i] - P[j])) - 1.0) for (i, j) in E)
    mind = min(math.hypot(*(P[i] - P[j])) for i in range(n) for j in range(i + 1, n))
    units = [(i, j) for i in range(n) for j in range(i + 1, n)
             if abs(math.hypot(*(P[i] - P[j])) - 1.0) < 1e-7]
    extra = [(i, j) for (i, j) in units if not adj[i][j]]
    return er, mind, len(units), extra


def attempt(adj, starts, label):
    n = len(adj)
    E = edges_of(adj)
    best = None
    for s in range(starts):
        rng = np.random.default_rng(101 * s + 7)
        P0 = mds_init(adj, rng, jitter=0.05 + 0.1 * (s % 4))
        P = coulomb_anneal(P0.copy(), E, n, rng)
        P = gauss_newton_polish(P, E, n)
        er, mind, nu, extra = analyze(P, adj)
        ok = er < 1e-9 and mind > 1e-3
        key = (0 if ok else 1, er + (1.0 if mind <= 1e-3 else 0.0))
        if best is None or key < best[0]:
            best = (key, P, er, mind, nu, extra)
        if ok:
            break
    _, P, er, mind, nu, extra = best
    realizable = er < 1e-9 and mind > 1e-3
    faithful = realizable and len(extra) == 0 and nu == len(E)
    print(f"[{label}] edge_res={er:.2e} min_dist={mind:.5f} units={nu} "
          f"extra={len(extra)} REALIZABLE={realizable} FAITHFUL={faithful}")
    return P, realizable, faithful


def hexagon_center():
    adj = [[0] * 7 for _ in range(7)]
    def add(a, b): adj[a][b] = adj[b][a] = 1
    for k in range(1, 7):
        add(0, k)
    for k in range(1, 7):
        add(k, k % 6 + 1)
    return adj


def main():
    print("SELF-TEST (Coulomb-anneal):")
    attempt(hexagon_center(), 6, "hexagon+center u7=12")
    print("\nCORES:")
    for idx, code in enumerate(N21_GRAPH6, 1):
        adj = graph6_decode(code)
        P, realizable, faithful = attempt(adj, 24, f"core {idx}")
        if realizable:
            np.save(f"05-knowledge/results/u21_core{idx}_v4coords.npy", np.array(P))
            with open(f"05-knowledge/results/u21_core{idx}_coords.txt", "w") as f:
                f.write(f"# AMP extremal u(21)=57 core {idx}; graph6={code}\n")
                for i, (x, y) in enumerate(P):
                    f.write(f"{i}\t{x:.18f}\t{y:.18f}\n")
            print(f"   -> coords saved for core {idx}")
            break


if __name__ == "__main__":
    main()
