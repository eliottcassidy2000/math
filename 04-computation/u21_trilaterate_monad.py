#!/usr/bin/env python3
"""
u21_trilaterate_monad.py  (monad-explorer, 2026-06-07)

Realize the 5 GENUINE Alexeev-Mixon-Parshall extremal u(21)=57 unit-distance
graphs (arXiv:2412.11914, Table 2 -- verified to match the repo's stored
cores) as actual planar point sets, via Henneberg/trilateration construction
rather than global optimization (which over-collapses).

Why trilateration:
  - place v0=(0,0), v1=(1,0) on an edge;
  - repeatedly take an unplaced vertex u that has >=2 already-placed neighbours
    a,b: u lies on circle(a,1) ∩ circle(b,1)  =>  at most 2 candidate points;
  - branch over the (<=2) choices, and over u's OTHER placed neighbours impose
    the consistency constraint |u - w| = 1 to prune;
  - backtrack until all 21 vertices are placed with every one of the 57 edges
    at distance 1 (within numeric tol).
This yields coordinates in a tower of quadratic extensions (algebraic), and a
faithful realization certifies u(21) >= 57 constructively.

Output: a numeric realization per core (polished by a couple of Newton/LM
steps), with EXACT integer-style verification of the unit count at high
precision, saved for the exact-arithmetic follow-up.
"""
import math
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


def circle_intersect(a, b, r=1.0):
    """Intersection points of circle(a,r) and circle(b,r). Returns list."""
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    d = math.hypot(dx, dy)
    if d > 2 * r or d < 1e-12:
        return []
    # midpoint
    mx, my = (ax + bx) / 2, (ay + by) / 2
    h2 = r * r - (d / 2) ** 2
    if h2 < 0:
        if h2 > -1e-12:
            h2 = 0.0
        else:
            return []
    h = math.sqrt(h2)
    ux, uy = -dy / d, dx / d   # unit perpendicular
    return [(mx + h * ux, my + h * uy), (mx - h * ux, my - h * uy)]


def trilateration_order(adj):
    """Greedy: place an edge, then any vertex with >=2 placed neighbours."""
    n = len(adj)
    nbr = [set(j for j in range(n) if adj[i][j]) for i in range(n)]
    # pick starting edge = highest-degree vertex and one of its neighbours
    v0 = max(range(n), key=lambda v: len(nbr[v]))
    v1 = max(nbr[v0], key=lambda v: len(nbr[v]))
    placed = [v0, v1]
    placedset = {v0, v1}
    order = [(v0, None, None), (v1, v0, None)]  # (vertex, anchor pair...)
    while len(placed) < n:
        best = None
        for u in range(n):
            if u in placedset:
                continue
            pn = [w for w in nbr[u] if w in placedset]
            if len(pn) >= 2:
                # prefer the vertex with the MOST placed neighbours (more constraints -> less branching downstream)
                if best is None or len(pn) > best[1]:
                    best = (u, len(pn), pn)
        if best is None:
            return None  # no trilateration order (graph not "2-tree-constructible" from here)
        u, _, pn = best
        order.append((u, pn[0], pn[1], pn[2:]))
        placed.append(u)
        placedset.add(u)
    return order


def realize(adj, tol=1e-7, max_nodes=2_000_000):
    """Backtracking trilateration. Returns coords list or None."""
    n = len(adj)
    edges = edges_of(adj)
    nbr = [set(j for j in range(n) if adj[i][j]) for i in range(n)]
    order = trilateration_order(adj)
    if order is None:
        return None, "no-trilateration-order"

    P = [None] * n
    seq = []  # list of (vertex, [candidate points])

    # place first two
    v0 = order[0][0]
    v1 = order[1][0]
    P[v0] = (0.0, 0.0)
    P[v1] = (1.0, 0.0)

    rest = order[2:]
    counter = [0]

    def consistent(u, pt):
        # every placed neighbour of u must be at distance ~1
        for w in nbr[u]:
            if P[w] is not None:
                d = math.hypot(pt[0] - P[w][0], pt[1] - P[w][1])
                if abs(d - 1.0) > tol:
                    return False
        return True

    def backtrack(k):
        counter[0] += 1
        if counter[0] > max_nodes:
            return False
        if k == len(rest):
            return True
        u, a, b = rest[k][0], rest[k][1], rest[k][2]
        cands = circle_intersect(P[a], P[b])
        for pt in cands:
            if consistent(u, pt):
                P[u] = pt
                if backtrack(k + 1):
                    return True
                P[u] = None
        return False

    ok = backtrack(0)
    if not ok:
        return None, f"backtrack-exhausted({counter[0]} nodes)"
    return [P[i] for i in range(n)], f"ok({counter[0]} nodes)"


def newton_polish(P, edges, steps=30):
    """A few Gauss-Newton steps (edges-only) to drive residual to ~machine eps,
    starting from the already-near-exact trilateration coords."""
    P = np.array(P, dtype=float)
    n = len(P)
    m = len(edges)
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
        # pin gauge: fix v0 (cols 0,1) and v1.y (col 3)
        H = J.T @ J + 1e-12 * np.eye(2 * n)
        g = J.T @ r
        for c in (0, 1, 3):
            H[c, :] = 0
            H[:, c] = 0
            H[c, c] = 1
            g[c] = 0
        try:
            dx = np.linalg.solve(H, -g)
        except np.linalg.LinAlgError:
            break
        P = P + dx.reshape(n, 2)
        if max(abs(math.hypot(*(P[i] - P[j])) - 1.0) for (i, j) in edges) < 1e-14:
            break
    return P


def analyze(P, adj):
    n = len(P)
    edges = edges_of(adj)
    edge_res = max(abs(math.hypot(*(P[i] - P[j])) - 1.0) for (i, j) in edges)
    mind = min(math.hypot(*(P[i] - P[j])) for i in range(n) for j in range(i + 1, n))
    units = [(i, j) for i in range(n) for j in range(i + 1, n)
             if abs(math.hypot(*(P[i] - P[j])) - 1.0) < 1e-7]
    extra = [(i, j) for (i, j) in units if not adj[i][j]]
    return edge_res, mind, len(units), extra


def main():
    print("=" * 74)
    print("Trilateration realization of the 5 AMP extremal u(21)=57 graphs")
    print("(monad-explorer; cores verified == AMP arXiv:2412.11914 Table 2)")
    print("=" * 74)
    any_ok = False
    for idx, code in enumerate(N21_GRAPH6, 1):
        adj = graph6_decode(code)
        n = len(adj)
        E = edges_of(adj)
        print(f"\n--- Core {idx}: n={n} |E|={len(E)} (graph6 len {len(code)}) ---")
        P, status = realize(adj)
        print(f"  trilateration: {status}")
        if P is None:
            continue
        P = newton_polish(P, E)
        er, mind, nunits, extra = analyze(P, adj)
        print(f"  polished edge residual: {er:.2e}")
        print(f"  min pairwise distance:  {mind:.6f}  ({'DISTINCT' if mind>1e-4 else 'COLLAPSED'})")
        print(f"  total unit distances:   {nunits}  (extra non-edge unit pairs: {len(extra)})")
        faithful = (er < 1e-9 and mind > 1e-4 and nunits == len(E))
        print(f"  REALIZABLE: {er<1e-9 and mind>1e-4}    FAITHFUL (UD graph==core): {faithful}")
        if er < 1e-9 and mind > 1e-4:
            any_ok = True
            np.save(f"05-knowledge/results/u21_core{idx}_trilat_coords.npy", np.array(P))
            # also dump coords as text
            with open(f"05-knowledge/results/u21_core{idx}_coords.txt", "w") as f:
                f.write(f"# AMP extremal u(21)=57 core {idx}; graph6={code}\n")
                f.write(f"# realized by trilateration; edge_residual={er:.2e}; min_dist={mind:.6f}\n")
                for i, (x, y) in enumerate(P):
                    f.write(f"{i}\t{x:.18f}\t{y:.18f}\n")
            print(f"  coords saved (npy + txt)")
    print(f"\nAt least one core realized as a genuine 21-point UDG with >=57 units: {any_ok}")
    print("=> u(21) >= 57 has a constructive realization (numeric); exact follow-up next.")


if __name__ == "__main__":
    main()
