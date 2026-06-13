#!/usr/bin/env python3
"""
even_graph_burnside_s260.py — Burnside formula for even graph iso classes
opus-2026-03-23-S260

An even graph on n vertices = a subgraph of K_n where every vertex has even degree.
These form the cycle space of K_n over GF(2), which has dimension C(n-1,2).

By Burnside: V_even(n) = (1/n!) * sum_{g in S_n} Fix_even(g)

Fix_even(g) = number of even graphs fixed by g.

For a permutation g with cycle type (c_1, c_2, ..., c_k):
The edge orbits of g partition the C(n,2) edges. An even graph fixed by g
must include either all or none of each edge orbit. Additionally, for the
graph to be even, every vertex must have even degree.

For each edge orbit, let its orbit structure determine degree contributions.
The constraint is: for each vertex v, the sum of contributions from included
edge orbits to deg(v) must be even.

This is a system of linear equations over GF(2):
    A * x = 0  (mod 2)
where x ∈ {0,1}^{#edge_orbits} and A[v, orbit] = |{edges in orbit incident to v}| mod 2.

Fix_even(g) = 2^{dim(null_space(A))} = 2^{#edge_orbits - rank(A)}.
"""

from math import comb, factorial, gcd
from collections import Counter
import numpy as np

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield []; return
    for f in range(min(n, mx), 0, -1):
        for r in partitions(n - f, f): yield [f] + r

def conjugacy_class_size(n, ct):
    c = Counter(ct); r = factorial(n)
    for l, k in c.items(): r //= (l ** k) * factorial(k)
    return r

def edge_orbits_and_constraints(n, ct):
    """For a permutation with cycle type ct, compute edge orbits
    and the GF(2) constraint matrix for even-degree subgraphs."""
    # Build a representative permutation
    perm = [0] * n
    pos = 0
    for c in ct:
        for i in range(c):
            perm[pos + i] = pos + (i + 1) % c
        pos += c

    # Find edge orbits
    P = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(P)
    edge_to_idx = {e: i for i, e in enumerate(P)}

    visited = [False] * m
    orbits = []
    for idx in range(m):
        if visited[idx]:
            continue
        orbit = []
        i, j = P[idx]
        while True:
            e = (min(i,j), max(i,j))
            eidx = edge_to_idx[e]
            if visited[eidx]:
                break
            visited[eidx] = True
            orbit.append(eidx)
            i, j = perm[i], perm[j]
        orbits.append(orbit)

    num_orbits = len(orbits)

    # Build constraint matrix A (n x num_orbits) over GF(2)
    # A[v, orbit_idx] = |{edges in orbit incident to v}| mod 2
    A = np.zeros((n, num_orbits), dtype=np.int8)
    for oidx, orbit in enumerate(orbits):
        for eidx in orbit:
            i, j = P[eidx]
            A[i, oidx] = (A[i, oidx] + 1) % 2
            A[j, oidx] = (A[j, oidx] + 1) % 2

    return num_orbits, A

def gf2_rank(M):
    """Compute rank of matrix M over GF(2)."""
    M = M.copy()
    rows, cols = M.shape
    rank = 0
    for col in range(cols):
        # Find pivot
        found = -1
        for r in range(rank, rows):
            if M[r, col]:
                found = r
                break
        if found == -1:
            continue
        # Swap
        M[[rank, found]] = M[[found, rank]]
        # Eliminate
        for r in range(rows):
            if r != rank and M[r, col]:
                M[r] = (M[r] + M[rank]) % 2
        rank += 1
    return rank

def burnside_even_graphs(n):
    """Compute V_even(n) = number of even graph iso classes on n vertices."""
    nf = factorial(n)
    total = 0

    for ct in partitions(n):
        ct_list = list(ct)
        ccs = conjugacy_class_size(n, ct_list)
        num_orbits, A = edge_orbits_and_constraints(n, ct_list)
        r = gf2_rank(A)
        fix = 2 ** (num_orbits - r)
        total += ccs * fix

    return total // nf

def burnside_even_graphs_detail(n):
    """V_even(n) with detailed breakdown by cycle type."""
    nf = factorial(n)
    total = 0
    details = []

    for ct in partitions(n):
        ct_list = list(ct)
        ccs = conjugacy_class_size(n, ct_list)
        num_orbits, A = edge_orbits_and_constraints(n, ct_list)
        r = gf2_rank(A)
        fix = 2 ** (num_orbits - r)
        total += ccs * fix
        details.append({
            'ct': ct_list, 'ccs': ccs,
            'edge_orbits': num_orbits, 'rank': r,
            'nullity': num_orbits - r, 'fix': fix,
            'contribution': ccs * fix
        })

    return total // nf, details

if __name__ == '__main__':
    print("=" * 70)
    print("  EVEN GRAPH ISO CLASS COUNTS VIA BURNSIDE")
    print("  opus-2026-03-23-S260")
    print("=" * 70)

    # Also compute tournament iso class counts for comparison
    from tournament_metagraph import burnside_quantities

    print(f"\n{'n':>3} {'V_even':>10} {'V_tourn':>10} {'ratio':>8} {'V_graphs':>10}")
    for n in range(3, 16):
        V_even = burnside_even_graphs(n)
        bq = burnside_quantities(n)
        V_tourn = bq['V']

        # Also compute total graph iso classes (standard Burnside)
        nf = factorial(n)
        V_graphs = 0
        for ct in partitions(n):
            ct_list = list(ct)
            ccs = conjugacy_class_size(n, ct_list)
            # Count edge orbits under this permutation
            perm = [0] * n; pos = 0
            for c in ct_list:
                for i in range(c):
                    perm[pos + i] = pos + (i + 1) % c
                pos += c
            P = [(i,j) for i in range(n) for j in range(i+1,n)]
            m = len(P)
            eidx_map = {e: i for i, e in enumerate(P)}
            visited = [False] * m
            num_edge_orbits = 0
            for idx in range(m):
                if visited[idx]: continue
                i, j = P[idx]
                while True:
                    e = (min(i,j), max(i,j))
                    ei = eidx_map[e]
                    if visited[ei]: break
                    visited[ei] = True
                    i, j = perm[i], perm[j]
                num_edge_orbits += 1
            V_graphs += ccs * 2 ** num_edge_orbits
        V_graphs //= nf

        ratio_te = V_tourn / V_even if V_even > 0 else float('inf')
        print(f"{n:>3} {V_even:>10} {V_tourn:>10} {ratio_te:>8.2f} {V_graphs:>10}")

    # Detailed analysis at small n
    print(f"\n{'='*70}")
    print("  DETAILED ANALYSIS")
    print(f"{'='*70}")

    for n in range(3, 9):
        V_even, details = burnside_even_graphs_detail(n)
        print(f"\n  n = {n}: V_even = {V_even}")
        print(f"  {'cycle type':>20} {'|class|':>8} {'edge_orb':>9} {'rank':>5} "
              f"{'null':>5} {'fix':>10} {'contrib':>12}")
        for d in details:
            print(f"  {str(d['ct']):>20} {d['ccs']:>8} {d['edge_orbits']:>9} "
                  f"{d['rank']:>5} {d['nullity']:>5} {d['fix']:>10} {d['contribution']:>12}")

    # Check: V_even = V_tourn (Royle) for all n?
    print(f"\n{'='*70}")
    print("  ROYLE EQUINUMEROSITY CHECK: V_even == V_tourn?")
    print(f"{'='*70}")
    for n in range(3, 16):
        V_even = burnside_even_graphs(n)
        bq = burnside_quantities(n)
        match = "✓" if V_even == bq['V'] else "✗"
        print(f"  n={n:>2}: V_even={V_even:>10}, V_tourn={bq['V']:>10}, {match}")

    # V_even sequence
    print(f"\n  V_even sequence: ", end="")
    seq = []
    for n in range(3, 16):
        seq.append(burnside_even_graphs(n))
    print(", ".join(str(s) for s in seq))

    print("\nDONE.")
