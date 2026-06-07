#!/usr/bin/env python3
"""
u21_embed_v3_monad.py  (monad-explorer, 2026-06-07)

A SELF-VALIDATED unit-distance embedder, then applied to the 5 stored
u(21)=57 "cores" to decide whether they are geometrically realizable.

Method:
  init   = classical MDS from BFS graph-distance (scaled so edges ~1)
  refine = damped Gauss-Newton (LM) on:
             - edge residuals      r = |p_i-p_j| - 1                (target 1)
             - non-edge repulsion  r = w*max(0, DMIN - |p_i-p_j|)   (keep apart)
           The repulsion forbids the vertex-collapse that fooled v2
           (which hit residual 1e-13 with two coincident vertices).
A core is REALIZABLE iff we find distinct points (min pair dist > eps) with all
edges at distance ~1.  Faithful iff additionally no non-edge is at distance ~1
(=> the unit-distance graph is exactly the core, 57 units => u(21) >= 57).

SELF-TEST first on graphs whose realizability is known:
  - hexagon+center  (u(7)=12, realizable, faithful)
  - triangular lattice 12-pt patch (realizable)
If the embedder solves those, failure on a core is evidence of NON-realizability.
"""
import math
import random
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
    D = [[0] * n for _ in range(n)]
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
            D[s][v] = dist[v] if dist[v] >= 0 else n  # disconnected -> n
    return np.array(D, dtype=float)


def mds_init(adj, rng):
    """Classical MDS on BFS graph distances -> 2D, jittered."""
    n = len(adj)
    D = graph_dist(adj)
    # mild scale: nearest-neighbour graph distance 1 ~ euclidean 1
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    w, V = np.linalg.eigh(B)
    order = np.argsort(w)[::-1]
    w = w[order]
    V = V[:, order]
    L = np.sqrt(np.maximum(w[:2], 0))
    P = V[:, :2] * L
    P = P + rng.normal(scale=0.05, size=P.shape)
    return P


def lm_refine(n, edges, nonedges, P, iters=400, DMIN=0.55, wrep=1.0):
    m = len(edges)
    mr = len(nonedges)
    lam = 1e-2
    for it in range(iters):
        rows = m + mr
        J = np.zeros((rows, 2 * n))
        r = np.zeros(rows)
        for e, (i, j) in enumerate(edges):
            d = P[i] - P[j]
            dist = math.hypot(d[0], d[1]) or 1e-12
            r[e] = dist - 1.0
            gx, gy = d[0] / dist, d[1] / dist
            J[e, 2 * i], J[e, 2 * i + 1] = gx, gy
            J[e, 2 * j], J[e, 2 * j + 1] = -gx, -gy
        for q, (i, j) in enumerate(nonedges):
            e = m + q
            d = P[i] - P[j]
            dist = math.hypot(d[0], d[1]) or 1e-12
            if dist < DMIN:
                r[e] = wrep * (DMIN - dist)
                gx, gy = d[0] / dist, d[1] / dist
                # d r/d P_i = -wrep * gx
                J[e, 2 * i], J[e, 2 * i + 1] = -wrep * gx, -wrep * gy
                J[e, 2 * j], J[e, 2 * j + 1] = wrep * gx, wrep * gy
            # else residual 0, jacobian 0
        cost = float(r @ r)
        H = J.T @ J
        g = J.T @ r
        try:
            dx = np.linalg.solve(H + lam * np.eye(2 * n), -g)
        except np.linalg.LinAlgError:
            lam *= 10
            continue
        Pnew = P + dx.reshape(n, 2)
        rn = np.zeros(rows)
        for e, (i, j) in enumerate(edges):
            rn[e] = math.hypot(*(Pnew[i] - Pnew[j])) - 1.0
        for q, (i, j) in enumerate(nonedges):
            dd = math.hypot(*(Pnew[i] - Pnew[j]))
            rn[m + q] = wrep * max(0.0, DMIN - dd)
        newcost = float(rn @ rn)
        if newcost < cost:
            P = Pnew
            lam = max(lam * 0.5, 1e-10)
            if cost - newcost < 1e-18:
                break
        else:
            lam *= 3
            if lam > 1e10:
                break
    edge_res = max(abs(math.hypot(*(P[i] - P[j])) - 1.0) for (i, j) in edges)
    return P, edge_res


def all_nonedges(adj):
    n = len(adj)
    return [(i, j) for i in range(n) for j in range(i + 1, n) if not adj[i][j]]


def analyze(P, edges, adj, tol=1e-5):
    n = len(P)
    edge_res = max(abs(math.hypot(*(P[i] - P[j])) - 1.0) for (i, j) in edges)
    mind = min(math.hypot(*(P[i] - P[j])) for i in range(n) for j in range(i + 1, n))
    units = []
    for i in range(n):
        for j in range(i + 1, n):
            if abs(math.hypot(*(P[i] - P[j])) - 1.0) < tol:
                units.append((i, j))
    nonedge_units = [(i, j) for (i, j) in units if not adj[i][j]]
    return edge_res, mind, len(units), nonedge_units


def try_embed(adj, starts=60, label=""):
    n = len(adj)
    E = edges_of(adj)
    NE = all_nonedges(adj)
    best = None
    for s in range(starts):
        rng = np.random.default_rng(13 * s + 1)
        if s % 2 == 0:
            P0 = mds_init(adj, rng)
        else:
            P0 = rng.uniform(-math.sqrt(n) / 2, math.sqrt(n) / 2, size=(n, 2))
        P, res = lm_refine(n, E, NE, P0, iters=300)
        er, mind, nunits, ne_units = analyze(P, E, adj)
        # quality: edges at 1 AND distinct points
        ok = er < 1e-5 and mind > 1e-3
        score = (0 if ok else 1, er + (1.0 if mind < 1e-3 else 0.0))
        if best is None or score < best[0]:
            best = (score, P, er, mind, nunits, ne_units)
        if ok and er < 1e-8:
            break
    _, P, er, mind, nunits, ne_units = best
    return P, er, mind, nunits, ne_units, E


def report(label, adj, starts):
    P, er, mind, nunits, ne_units, E = try_embed(adj, starts=starts, label=label)
    print(f"\n[{label}] n={len(adj)} |E|={len(E)}")
    print(f"   edge residual:   {er:.2e}")
    print(f"   min pair dist:   {mind:.5f}  ({'DISTINCT' if mind>1e-3 else 'DEGENERATE/collapsed'})")
    print(f"   total unit dist: {nunits}   (extra non-edge unit pairs: {len(ne_units)})")
    realizable = er < 1e-5 and mind > 1e-3
    faithful = realizable and len(ne_units) == 0 and nunits == len(E)
    print(f"   REALIZABLE: {realizable}   FAITHFUL: {faithful}")
    return P, realizable, faithful, nunits


# ---------- known-realizable test graphs ----------

def hexagon_center():
    # 7 vertices: center 0, hexagon 1..6
    adj = [[0] * 7 for _ in range(7)]

    def add(a, b):
        adj[a][b] = adj[b][a] = 1
    for k in range(1, 7):
        add(0, k)
    for k in range(1, 7):
        add(k, k % 6 + 1)
    return adj  # 12 edges, u(7)=12


def tri_patch(rings=1):
    # triangular lattice disk, eisenstein coords, unit^2 = 1 (nearest nbr)
    pts = []
    for x in range(-3, 4):
        for y in range(-3, 4):
            if x * x + x * y + y * y <= rings * rings + 1:
                pts.append((x, y))
    n = len(pts)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            dx = pts[i][0] - pts[j][0]
            dy = pts[i][1] - pts[j][1]
            if dx * dx + dx * dy + dy * dy == 1:
                adj[i][j] = adj[j][i] = 1
    return adj


def main():
    print("=" * 72)
    print("SELF-TEST: embedder on KNOWN-realizable unit-distance graphs")
    print("=" * 72)
    report("hexagon+center (u7=12)", hexagon_center(), starts=20)
    report("triangular patch", tri_patch(2), starts=20)

    print("\n" + "=" * 72)
    print("THE 5 STORED u(21)=57 CORES")
    print("=" * 72)
    for idx, code in enumerate(N21_GRAPH6, 1):
        adj = graph6_decode(code)
        P, realizable, faithful, nunits = report(f"core {idx}", adj, starts=80)
        if realizable:
            np.save(f"05-knowledge/results/u21_core{idx}_v3_coords.npy", P)


if __name__ == "__main__":
    main()
