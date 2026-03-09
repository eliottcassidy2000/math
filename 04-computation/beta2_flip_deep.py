#!/usr/bin/env python3
"""
beta2_flip_deep.py — Flip Obstruction Analysis

Investigates why 3-cycles among "bad" vertices (those with beta_1(T\v)=1)
force beta_1(T)>=1, while generic 3-cycles don't.

Analysis steps:
1. For each tournament, classify 3-cycles by how many bad vertices they contain
2. Test: does beta_1=0 imply no all-bad 3-cycle?
3. Cocycle extension analysis for bad vertex triples
4. Flip analysis: transitive bad triples -> what if we flip to make a 3-cycle?
"""

import sys
import os
import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
import random

# Import from path_homology_v2 — but it runs validation on import.
# We'll suppress that by redirecting stdout temporarily.
import io
old_stdout = sys.stdout
sys.stdout = io.StringIO()
sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')
import path_homology_v2 as ph
sys.stdout = old_stdout

# Use ph functions
path_betti_numbers = ph.path_betti_numbers
enumerate_allowed_paths = ph.enumerate_allowed_paths
build_full_boundary_matrix = ph.build_full_boundary_matrix
compute_omega_basis = ph.compute_omega_basis

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def find_3cycles(A, n):
    """Return list of 3-cycles as frozensets of vertices."""
    cycles = []
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    cycles.append(frozenset([i,j,k]))
    return cycles

def subtournament(A, n, exclude_v):
    """Return adjacency matrix of T \ {v}."""
    verts = [u for u in range(n) if u != exclude_v]
    m = len(verts)
    B = [[0]*m for _ in range(m)]
    for i2, u in enumerate(verts):
        for j2, w in enumerate(verts):
            B[i2][j2] = A[u][w]
    return B, m, verts

def find_bad_vertices(A, n):
    """Find vertices v where beta_1(T \ {v}) = 1."""
    bad = []
    for v in range(n):
        B, m, _ = subtournament(A, n, v)
        betti = path_betti_numbers(B, m, max_dim=1)
        if len(betti) > 1 and betti[1] >= 1:
            bad.append(v)
    return bad

def flip_arc(A, n, u, v):
    """Return a copy of A with arc u->v flipped to v->u (or vice versa)."""
    B = [row[:] for row in A]
    B[u][v] = 1 - B[u][v]
    B[v][u] = 1 - B[v][u]
    return B

def is_3cycle_among(A, verts):
    """Check if the 3 vertices form a directed 3-cycle."""
    u, v, w = verts
    return ((A[u][v] and A[v][w] and A[w][u]) or
            (A[v][u] and A[w][v] and A[u][w]))

def edge_pattern(A, verts):
    """Describe edge pattern among 3 vertices."""
    u, v, w = sorted(verts)
    edges = []
    for a, b in [(u,v),(u,w),(v,w)]:
        if A[a][b]:
            edges.append(f"{a}->{b}")
        else:
            edges.append(f"{b}->{a}")
    return ", ".join(edges)

def compute_cocycle_info(A, n):
    """Compute Z^1 and B^1 dimensions (cohomology via transpose)."""
    # H^1 = ker(delta_1) / im(delta_0) where delta is coboundary
    # dim H^1 = beta_1 (by universal coefficients over Q)
    # Actually for path homology, H^p and H_p are the same over Q.
    # Let's just compute beta_1 directly and also get the cycle space info.

    allowed_0 = enumerate_allowed_paths(A, n, 0)
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_2 = enumerate_allowed_paths(A, n, 2)

    omega_1 = compute_omega_basis(A, n, 1, allowed_1, allowed_0)

    # boundary d1: Omega_1 -> A_0
    bd_1 = build_full_boundary_matrix(allowed_1, allowed_0)
    bd_1_omega = bd_1 @ omega_1

    # kernel of d1 restricted to Omega_1 = Z_1
    if bd_1_omega.shape[0] > 0 and bd_1_omega.shape[1] > 0:
        U, S, Vt = np.linalg.svd(bd_1_omega, full_matrices=True)
        S = np.atleast_1d(S)
        rank_d1 = int(np.sum(S > 1e-8))
        ker_d1 = Vt[rank_d1:].T  # in Omega_1 coords
    else:
        rank_d1 = 0
        ker_d1 = omega_1

    # boundary d2: Omega_2 -> A_1
    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    bd_2 = build_full_boundary_matrix(allowed_2, allowed_1)
    bd_2_omega = bd_2 @ omega_2

    if bd_2_omega.shape[0] > 0 and bd_2_omega.shape[1] > 0:
        S2 = np.linalg.svd(bd_2_omega, compute_uv=False)
        S2 = np.atleast_1d(S2)
        rank_d2 = int(np.sum(S2 > 1e-8))
    else:
        rank_d2 = 0

    dim_omega_1 = omega_1.shape[1] if omega_1.ndim == 2 else 0
    dim_Z1 = dim_omega_1 - rank_d1
    dim_B1 = rank_d2
    beta_1 = dim_Z1 - dim_B1

    # Get a representative cycle if beta_1 > 0
    rep_cycle = None
    if beta_1 > 0 and ker_d1 is not None and ker_d1.shape[1] > 0:
        # Z_1 basis in A_1 coords
        z1_basis = omega_1 @ ker_d1
        # B_1 basis in A_1 coords
        if bd_2_omega.shape[1] > 0:
            b1_basis = bd_2_omega
        else:
            b1_basis = np.zeros((len(allowed_1), 0))

        # Find element in Z_1 not in B_1
        # Project z1 onto complement of B1
        if b1_basis.shape[1] > 0:
            # Orthogonal complement
            Q, R = np.linalg.qr(b1_basis)
            for i in range(z1_basis.shape[1]):
                z = z1_basis[:, i]
                proj = Q @ (Q.T @ z)
                residual = z - proj
                if np.linalg.norm(residual) > 1e-8:
                    rep_cycle = residual / np.linalg.norm(residual)
                    break
        else:
            rep_cycle = z1_basis[:, 0]
            rep_cycle = rep_cycle / np.linalg.norm(rep_cycle)

    return {
        'beta_1': beta_1,
        'dim_omega_1': dim_omega_1,
        'dim_Z1': dim_Z1,
        'dim_B1': dim_B1,
        'allowed_1': allowed_1,
        'rep_cycle': rep_cycle
    }

def sample_tournaments(n, count):
    """Sample random tournaments at given n."""
    tournaments = []
    for _ in range(count):
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        tournaments.append(A)
    return tournaments


# ================================================================
# MAIN ANALYSIS
# ================================================================

print("=" * 70)
print("FLIP OBSTRUCTION ANALYSIS")
print("3-cycles among bad vertices and beta_1 forcing")
print("=" * 70)

for n in [5, 6]:
    print(f"\n{'='*70}")
    print(f"n = {n}  (exhaustive, {1 << (n*(n-1)//2)} tournaments)")
    print(f"{'='*70}")

    # Counters
    total = 0
    b1_0_count = 0
    b1_1_count = 0
    b1_0_no_bad = 0  # beta_1=0, no bad vertices at all

    # 3-cycle classification: for beta_1=0 and beta_1>=1
    # bad_in_cycle[beta_1_val][num_bad_in_cycle] = count
    bad_in_cycle = {0: Counter(), 1: Counter()}

    # Track specific patterns
    b1_0_all_bad_cycle = []  # beta_1=0 tournaments with all-bad 3-cycle
    b1_1_all_bad_cycle = []  # beta_1=1 tournaments with all-bad 3-cycle
    b1_0_has_3bad_vertices = []  # beta_1=0, has >=3 bad vertices
    b1_0_3bad_transitive = []  # beta_1=0, 3 bad vertices form transitive triple

    # Distribution of number of bad vertices
    bad_count_dist = {0: Counter(), 1: Counter()}

    for A in all_tournaments(n):
        total += 1
        t3 = count_3cycles(A, n)
        betti = path_betti_numbers(A, n, max_dim=1)
        b1 = betti[1] if len(betti) > 1 else 0

        if b1 == 0:
            b1_0_count += 1
        else:
            b1_1_count += 1

        b1_key = min(b1, 1)  # group beta_1>=1 together

        # Skip t3=0 (transitive tournaments) — no 3-cycles to analyze
        if t3 == 0:
            bad_count_dist[b1_key][0] += 1
            continue

        # Find bad vertices
        bad = find_bad_vertices(A, n)
        bad_set = set(bad)
        bad_count_dist[b1_key][len(bad)] += 1

        if len(bad) == 0:
            if b1 == 0:
                b1_0_no_bad += 1
            continue

        # Find 3-cycles and classify by bad vertex count
        cycles = find_3cycles(A, n)
        for cyc in cycles:
            num_bad = len(cyc & bad_set)
            bad_in_cycle[b1_key][num_bad] += 1

        # Check for all-bad 3-cycles
        has_all_bad_cycle = False
        for cyc in cycles:
            if len(cyc & bad_set) == 3:
                has_all_bad_cycle = True
                break

        if has_all_bad_cycle:
            if b1 == 0:
                b1_0_all_bad_cycle.append(A)
            else:
                b1_1_all_bad_cycle.append(A)

        # Check for 3 bad vertices in transitive triple
        if b1 == 0 and len(bad) >= 3:
            b1_0_has_3bad_vertices.append((A, bad))
            for triple in combinations(bad, 3):
                if not is_3cycle_among(A, triple):
                    b1_0_3bad_transitive.append((A, triple, bad))

    # --- REPORT ---
    print(f"\nTotal tournaments: {total}")
    print(f"  beta_1 = 0: {b1_0_count}")
    print(f"  beta_1 >= 1: {b1_1_count}")

    print(f"\n--- STEP 1 & 2: 3-cycle classification by bad vertex count ---")
    for b1_val, label in [(0, "beta_1=0"), (1, "beta_1>=1")]:
        print(f"\n  {label}:")
        print(f"    Distribution of #bad vertices in tournament:")
        for k in sorted(bad_count_dist[b1_val]):
            print(f"      {k} bad vertices: {bad_count_dist[b1_val][k]} tournaments")
        print(f"    3-cycles with k bad vertices:")
        for k in sorted(bad_in_cycle[b1_val]):
            print(f"      k={k}: {bad_in_cycle[b1_val][k]} 3-cycles")

    print(f"\n--- STEP 3: THE KEY TEST ---")
    print(f"  beta_1=0 with all-bad 3-cycle: {len(b1_0_all_bad_cycle)}")
    print(f"  beta_1>=1 with all-bad 3-cycle: {len(b1_1_all_bad_cycle)}")

    if len(b1_0_all_bad_cycle) > 0:
        print(f"  *** CONJECTURE FAILS: found beta_1=0 with all-bad 3-cycle ***")
        for A in b1_0_all_bad_cycle[:3]:
            bad = find_bad_vertices(A, n)
            print(f"    bad vertices: {bad}")
            for cyc in find_3cycles(A, n):
                if len(cyc & set(bad)) == 3:
                    print(f"    all-bad cycle: {sorted(cyc)}")
    else:
        print(f"  CONFIRMED: beta_1=0 => no 3-cycle has all 3 vertices bad")

    # --- STEP 4: Cocycle analysis for bad triples with beta_1>=1 ---
    print(f"\n--- STEP 4: Cocycle analysis (beta_1>=1, all-bad 3-cycles) ---")
    analyzed = 0
    for A in b1_1_all_bad_cycle[:5]:  # analyze up to 5
        bad = find_bad_vertices(A, n)
        bad_set = set(bad)
        cycles = find_3cycles(A, n)
        all_bad_cycles = [cyc for cyc in cycles if len(cyc & bad_set) == 3]

        if not all_bad_cycles:
            continue
        analyzed += 1
        cyc = all_bad_cycles[0]
        verts = sorted(cyc)

        print(f"\n  Tournament (first all-bad cycle: {verts}):")
        print(f"    Edges among cycle: {edge_pattern(A, verts)}")
        print(f"    Bad vertices: {bad}")

        # For each bad vertex, compute cocycle info on T\v
        for v in verts:
            B, m, vert_map = subtournament(A, n, v)
            info = compute_cocycle_info(B, m)
            print(f"    T\\{v}: beta_1={info['beta_1']}, dim(Omega_1)={info['dim_omega_1']}, "
                  f"dim(Z_1)={info['dim_Z1']}, dim(B_1)={info['dim_B1']}")
            if info['rep_cycle'] is not None:
                # Show which edges the representative cycle uses
                cycle_vec = info['rep_cycle']
                nonzero = [(info['allowed_1'][i], cycle_vec[i])
                           for i in range(len(cycle_vec)) if abs(cycle_vec[i]) > 1e-8]
                # Map back to original vertex labels
                mapped = [((vert_map[e[0]], vert_map[e[1]]), coeff) for e, coeff in nonzero]
                print(f"      Rep cycle (in T): {[(e, round(c,3)) for e,c in mapped[:8]]}")

        # Compute cocycle info for full T
        info_T = compute_cocycle_info(A, n)
        print(f"    Full T: beta_1={info_T['beta_1']}, dim(Z_1)={info_T['dim_Z1']}, dim(B_1)={info_T['dim_B1']}")
        if info_T['rep_cycle'] is not None:
            cycle_vec = info_T['rep_cycle']
            nonzero = [(info_T['allowed_1'][i], cycle_vec[i])
                       for i in range(len(cycle_vec)) if abs(cycle_vec[i]) > 1e-8]
            print(f"      Rep cycle: {[(e, round(c,3)) for e,c in nonzero[:8]]}")

    if analyzed == 0:
        print("  (No all-bad 3-cycles found for beta_1>=1)")

    # --- STEP 5: Transitive bad triples and flip analysis ---
    print(f"\n--- STEP 5: Transitive bad triples (beta_1=0) & flip analysis ---")
    print(f"  Tournaments with beta_1=0 and >=3 bad vertices: {len(b1_0_has_3bad_vertices)}")
    print(f"  Of those, triples forming transitive triple: {len(b1_0_3bad_transitive)}")

    flip_forces_b1 = 0
    flip_keeps_b1_0 = 0
    flip_details = []

    for A, triple, all_bad in b1_0_3bad_transitive[:20]:  # analyze up to 20
        u, v, w = triple
        # Find the edge pattern
        pattern = edge_pattern(A, triple)

        # Try all possible single arc flips among the triple to see if any creates a 3-cycle
        flips_creating_cycle = []
        for a, b in [(u,v),(u,w),(v,w),(v,u),(w,u),(w,v)]:
            if A[a][b] == 1:
                B = flip_arc(A, n, a, b)
                if is_3cycle_among(B, triple):
                    betti_B = path_betti_numbers(B, n, max_dim=1)
                    b1_B = betti_B[1] if len(betti_B) > 1 else 0
                    flips_creating_cycle.append((a, b, b1_B))
                    if b1_B >= 1:
                        flip_forces_b1 += 1
                    else:
                        flip_keeps_b1_0 += 1

        if flips_creating_cycle:
            flip_details.append((triple, pattern, all_bad, flips_creating_cycle))

    print(f"\n  Flip results (flipping arc to create 3-cycle among bad vertices):")
    print(f"    Flips that force beta_1>=1: {flip_forces_b1}")
    print(f"    Flips that keep beta_1=0:   {flip_keeps_b1_0}")

    if flip_details:
        print(f"\n  Detailed flip analysis (first {min(10, len(flip_details))}):")
        for triple, pattern, all_bad, flips in flip_details[:10]:
            print(f"    Triple {triple}, edges: {pattern}")
            print(f"      All bad vertices in T: {all_bad}")
            for a, b, b1_B in flips:
                print(f"      Flip {a}->{b}: beta_1 becomes {b1_B}")

# ================================================================
# n=7 SAMPLING
# ================================================================
print(f"\n{'='*70}")
print(f"n = 7  (sampling 500 random tournaments)")
print(f"{'='*70}")

random.seed(42)
n = 7
tournaments_7 = sample_tournaments(n, 500)

b1_0_count = 0
b1_1_count = 0
bad_in_cycle_7 = {0: Counter(), 1: Counter()}
bad_count_dist_7 = {0: Counter(), 1: Counter()}
b1_0_all_bad_cycle_7 = 0
b1_1_all_bad_cycle_7 = 0
b1_0_3bad_transitive_7 = []
flip_forces_7 = 0
flip_keeps_7 = 0

for idx, A in enumerate(tournaments_7):
    t3 = count_3cycles(A, n)
    betti = path_betti_numbers(A, n, max_dim=1)
    b1 = betti[1] if len(betti) > 1 else 0
    b1_key = min(b1, 1)

    if b1 == 0:
        b1_0_count += 1
    else:
        b1_1_count += 1

    if t3 == 0:
        bad_count_dist_7[b1_key][0] += 1
        continue

    bad = find_bad_vertices(A, n)
    bad_set = set(bad)
    bad_count_dist_7[b1_key][len(bad)] += 1

    if len(bad) == 0:
        continue

    cycles = find_3cycles(A, n)
    for cyc in cycles:
        num_bad = len(cyc & bad_set)
        bad_in_cycle_7[b1_key][num_bad] += 1

    has_all_bad_cycle = any(len(cyc & bad_set) == 3 for cyc in cycles)
    if has_all_bad_cycle:
        if b1 == 0:
            b1_0_all_bad_cycle_7 += 1
        else:
            b1_1_all_bad_cycle_7 += 1

    # Flip analysis for transitive bad triples with beta_1=0
    if b1 == 0 and len(bad) >= 3:
        for triple in combinations(bad, 3):
            if not is_3cycle_among(A, triple):
                u, v, w = triple
                for a, b_v in [(u,v),(u,w),(v,w),(v,u),(w,u),(w,v)]:
                    if A[a][b_v] == 1:
                        B = flip_arc(A, n, a, b_v)
                        if is_3cycle_among(B, triple):
                            betti_B = path_betti_numbers(B, n, max_dim=1)
                            b1_B = betti_B[1] if len(betti_B) > 1 else 0
                            if b1_B >= 1:
                                flip_forces_7 += 1
                            else:
                                flip_keeps_7 += 1

    if (idx + 1) % 100 == 0:
        print(f"  ...processed {idx+1}/500")

print(f"\nResults (n=7 sample):")
print(f"  beta_1=0: {b1_0_count}, beta_1>=1: {b1_1_count}")
print(f"\n  Bad vertex distribution:")
for b1_val, label in [(0, "beta_1=0"), (1, "beta_1>=1")]:
    print(f"    {label}:")
    for k in sorted(bad_count_dist_7[b1_val]):
        print(f"      {k} bad: {bad_count_dist_7[b1_val][k]}")

print(f"\n  3-cycles by bad vertex count:")
for b1_val, label in [(0, "beta_1=0"), (1, "beta_1>=1")]:
    print(f"    {label}:")
    for k in sorted(bad_in_cycle_7[b1_val]):
        print(f"      k={k}: {bad_in_cycle_7[b1_val][k]}")

print(f"\n  KEY TEST:")
print(f"    beta_1=0 with all-bad 3-cycle: {b1_0_all_bad_cycle_7}")
print(f"    beta_1>=1 with all-bad 3-cycle: {b1_1_all_bad_cycle_7}")

print(f"\n  FLIP ANALYSIS (transitive bad triples, beta_1=0):")
print(f"    Flips forcing beta_1>=1: {flip_forces_7}")
print(f"    Flips keeping beta_1=0:  {flip_keeps_7}")

# ================================================================
# SUMMARY
# ================================================================
print(f"\n{'='*70}")
print("SUMMARY OF FINDINGS")
print(f"{'='*70}")
print("""
Key questions answered:
1. Does beta_1=0 imply no 3-cycle has all 3 vertices bad?
2. Does having an all-bad 3-cycle force beta_1>=1?
3. If 3 bad vertices form a transitive triple, does flipping
   one arc (to create a 3-cycle) force beta_1 >= 1?

See detailed results above for each n.
""")
print("Done.")
