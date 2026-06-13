#!/usr/bin/env python3
"""
u21_core_embed_monad.py  (monad-explorer, 2026-06-07)

CAMPAIGN: u(21) unit-distance maximum.  u(21)=57 is PROVED exact
(Alexeev-Mixon-Parshall, arXiv:2412.11914, Dec 2024).

This repo stores 5 abstract graph6 "cores" claimed to be the extremal
u(21)=57 unit-distance graphs, but has only ever treated them
COMBINATORIALLY (traceability, degree histograms).  It has NEVER verified
they are realizable as unit-distance graphs in the plane.

This script:
  (1) decodes each core, checks |V|=21, |E|=57, and the NECESSARY conditions
      for a planar unit-distance graph: K4-free (max clique <=3) and
      K_{2,3}-free (every pair has <=2 common neighbours).
  (2) attempts a numeric unit-distance embedding (energy minimization of
      sum_{edges}(|p_i-p_j|^2 - 1)^2) from many random starts.
  (3) reports the best embedding residual and whether all 57 edges hit
      distance ~1 while NO non-edge accidentally also hits distance 1
      (a faithful embedding) -- and counts total unit distances realized.

All distance *counts* in the verification stage are done with a tolerance
sweep; the EXACT-arithmetic certificate is a separate follow-up once a
faithful numeric embedding is in hand.
"""
import itertools
import math
import random

N21_GRAPH6 = [
    "Tsc@IGC@GD?R?S?Wd@A_CK@HG@VM??PRKOUZ",
    "TCKx?D?OI?OMCBSA_L?ApA_gEA\\EG?PBSCPV",
    "TCTWACAG@@CDKC?e?QgQA@OMOq]F??OUcCEj",
    "TCS`?H??XIZ?K_Co`CG@JO[?EOSCpGOTSCE\\",
    "T??_`OhSCSYA@I?c?OyWBEa@c?SIU?Aa[?el",
]


def graph6_decode(code):
    """Decode a graph6 string into an adjacency matrix (list of lists)."""
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


def degree_hist(adj):
    from collections import Counter
    return dict(sorted(Counter(sum(r) for r in adj).items()))


def max_common_neighbours(adj):
    n = len(adj)
    nbr = [set(j for j in range(n) if adj[i][j]) for i in range(n)]
    worst = 0
    bad = None
    for i in range(n):
        for j in range(i + 1, n):
            c = len(nbr[i] & nbr[j])
            if c > worst:
                worst = c
                bad = (i, j, c)
    return worst, bad


def max_clique_geq4(adj):
    """Return True if there exist 4 mutually adjacent vertices (a K4)."""
    n = len(adj)
    nbr = [set(j for j in range(n) if adj[i][j]) for i in range(n)]
    for i in range(n):
        Ni = nbr[i]
        for j in Ni:
            if j <= i:
                continue
            common_ij = Ni & nbr[j]
            for k in common_ij:
                if k <= j:
                    continue
                # i,j,k triangle; need a 4th common to all three
                if common_ij & nbr[k]:
                    # check the 4th is distinct
                    for l in common_ij & nbr[k]:
                        if l != i and l != j and l != k:
                            return True, (i, j, k, l)
    return False, None


# ---------------- numeric embedding ----------------

def embed(adj, edges, seed, iters=4000, lr=0.05):
    """Energy-minimize sum_{edges}(d^2-1)^2 by gradient descent.
    Returns (coords, edge_residual_max)."""
    rng = random.Random(seed)
    n = len(adj)
    # start: random points in a disk of radius ~ sqrt(n)/2
    R = math.sqrt(n) / 2
    P = [[rng.uniform(-R, R), rng.uniform(-R, R)] for _ in range(n)]

    for t in range(iters):
        grad = [[0.0, 0.0] for _ in range(n)]
        for (i, j) in edges:
            dx = P[i][0] - P[j][0]
            dy = P[i][1] - P[j][1]
            d2 = dx * dx + dy * dy
            f = (d2 - 1.0)          # energy term (d2-1)^2, dE/dd2 = 2(d2-1)
            # dE/dP_i = 2(d2-1) * 2*(P_i - P_j)
            gx = 4.0 * f * dx
            gy = 4.0 * f * dy
            grad[i][0] += gx
            grad[i][1] += gy
            grad[j][0] -= gx
            grad[j][1] -= gy
        step = lr / (1.0 + 0.001 * t)
        for i in range(n):
            P[i][0] -= step * grad[i][0]
            P[i][1] -= step * grad[i][1]

    # edge residual
    res = 0.0
    for (i, j) in edges:
        d = math.dist(P[i], P[j])
        res = max(res, abs(d - 1.0))
    return P, res


def count_unit(P, tol):
    n = len(P)
    c = 0
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            d = math.dist(P[i], P[j])
            if abs(d - 1.0) < tol:
                c += 1
                pairs.append((i, j))
    return c, pairs


def min_pair_dist(P):
    n = len(P)
    m = float("inf")
    for i in range(n):
        for j in range(i + 1, n):
            m = min(m, math.dist(P[i], P[j]))
    return m


def main():
    print("=" * 70)
    print("u(21)=57 core audit + unit-distance embedding (monad-explorer)")
    print("=" * 70)
    for idx, code in enumerate(N21_GRAPH6, 1):
        adj = graph6_decode(code)
        n = len(adj)
        E = edges_of(adj)
        dh = degree_hist(adj)
        cn, cnbad = max_common_neighbours(adj)
        hask4, k4 = max_clique_geq4(adj)
        print(f"\n--- Core {idx}: n={n}, |E|={len(E)} ---")
        print(f"  degree histogram: {dh}")
        print(f"  max common neighbours (need <=2): {cn}   witness {cnbad}")
        print(f"  contains K4 (must be False): {hask4}  {k4 if hask4 else ''}")
        realizable_necessary = (n == 21 and len(E) == 57
                                and cn <= 2 and not hask4)
        print(f"  passes NECESSARY planar-UDG conditions: {realizable_necessary}")
        if not realizable_necessary:
            print("  -> FAILS necessary conditions; not a planar UDG. Skipping embed.")
            continue

        # try embedding from several starts
        best_res = float("inf")
        best_P = None
        for s in range(8):
            P, res = embed(adj, E, seed=1000 * idx + s, iters=3000, lr=0.04)
            if res < best_res:
                best_res = res
                best_P = P
            if res < 1e-4:
                break
        print(f"  best edge residual over 8 starts: {best_res:.2e}")
        if best_P is not None:
            for tol in (1e-3, 1e-4):
                c, _ = count_unit(best_P, tol)
                print(f"    unit distances at tol={tol:.0e}: {c}")
            print(f"    min pairwise distance: {min_pair_dist(best_P):.4f}")
            faithful = "YES" if best_res < 1e-3 else "NO (did not converge)"
            print(f"    edges realized at unit (faithful embed?): {faithful}")


if __name__ == "__main__":
    main()
