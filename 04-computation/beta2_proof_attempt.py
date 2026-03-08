#!/usr/bin/env python3
"""
beta2_proof_attempt.py — Attempt to prove β₂=0 for tournaments algebraically

APPROACH: The order complex of the tournament's domination poset.

For a tournament T on [n], define the "domination order" (partial order):
  For vertex sets S, T: S ≤ T iff every vertex in S beats every vertex in T.

This is NOT a well-defined poset on vertices since T is not transitive.
Instead, consider the ACYCLIC condensation: contract each strongly connected
component to a point. The SCC quotient is a DAG.

For tournaments, every SCC is a "strongly connected tournament" and the
condensation is a TOTAL ORDER on SCCs (Rédei-type). The transitive triples
are exactly the triples from this total order.

Wait — that's not quite right. A transitive triple (a,b,c) means a→b→c
and a→c. This doesn't require a,b,c to be in different SCCs.

KEY QUESTION: What does a 2-cycle in ker(∂₂) look like concretely?

Let's study small cases exhaustively to build intuition.

Author: kind-pasteur-2026-03-08-S41
"""
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import enumerate_allowed_paths

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def transitive_triples(A, n):
    """Find all transitive triples (a,b,c) with a→b→c and a→c."""
    tt = []
    for a in range(n):
        for b in range(n):
            if b == a: continue
            if not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    tt.append((a,b,c))
    return tt

def dt_4paths(A, n):
    """Find all DT 4-paths (a,b,c,d) with a→b→c→d, a→c, b→d."""
    dt = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if not A[b][c] or not A[a][c]: continue
                for d in range(n):
                    if d == a or d == b or d == c: continue
                    if A[c][d] and A[b][d]:
                        dt.append((a,b,c,d))
    return dt

def compute_ker_boundary2(tt, a1_list, n):
    """Compute ker(∂₂: span(tt) → span(a1))."""
    a1_idx = {tuple(p): i for i, p in enumerate(a1_list)}
    bd2 = np.zeros((len(a1_list), len(tt)))
    for j, (a, b, c) in enumerate(tt):
        bd2[a1_idx[(b,c)], j] += 1
        bd2[a1_idx[(a,c)], j] -= 1
        bd2[a1_idx[(a,b)], j] += 1
    rank = np.linalg.matrix_rank(bd2, tol=1e-8)
    ker_dim = len(tt) - rank

    if ker_dim > 0:
        # Actually compute the kernel
        U, S, Vt = np.linalg.svd(bd2, full_matrices=True)
        null_mask = np.abs(S) < 1e-8 if len(S) > 0 else np.array([])
        # Null space is last columns of V where s_i ≈ 0
        null_start = len(S) - sum(np.abs(S) < 1e-8) if len(S) > 0 else 0
        # Actually use the proper null space computation
        _, s_vals, Vt = np.linalg.svd(bd2)
        null_space = Vt[len(s_vals):].T if len(s_vals) < len(tt) else np.zeros((len(tt), 0))
        # Also include singular values near zero
        for i, sv in enumerate(s_vals):
            if sv < 1e-8:
                null_space = np.column_stack([null_space, Vt[i]])
        return ker_dim, null_space
    return 0, None

# ===== Study 2-cycles at n=4 and n=5 =====
print("=" * 70)
print("STRUCTURE OF ker(∂₂) — 2-CYCLES IN TOURNAMENTS")
print("=" * 70)

for n in [4, 5]:
    print(f"\n{'='*40}")
    print(f"n = {n}")
    print(f"{'='*40}")

    max_ker_dim = 0
    cycle_counts = Counter()

    for A in all_tournaments_gen(n):
        # Get allowed 1-paths
        a1 = [(i,j) for i in range(n) for j in range(n) if i != j and A[i][j]]

        tt = transitive_triples(A, n)
        if len(tt) == 0:
            continue

        ker_dim, null_space = compute_ker_boundary2(tt, a1, n)
        cycle_counts[ker_dim] += 1

        if ker_dim > max_ker_dim:
            max_ker_dim = ker_dim

        if ker_dim > 0 and n == 4:
            # Show what the 2-cycle looks like
            # Get DT paths
            dt = dt_4paths(A, n)

            # Check: does im(∂₃|DT) fill ker(∂₂)?
            if len(dt) > 0 and null_space is not None:
                tt_idx = {t: i for i, t in enumerate(tt)}
                bd3 = np.zeros((len(tt), len(dt)))
                for j, (v0,v1,v2,v3) in enumerate(dt):
                    faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
                    signs = [1, -1, 1, -1]
                    for face, sign in zip(faces, signs):
                        if face in tt_idx:
                            bd3[tt_idx[face], j] += sign

                im_rank = np.linalg.matrix_rank(bd3, tol=1e-8)
                print(f"\n  Example with ker_dim={ker_dim}:")
                print(f"    TT: {len(tt)}, DT: {len(dt)}")
                print(f"    im(∂₃|DT) rank = {im_rank}")
                print(f"    β₂ = ker_dim - im_rank = {ker_dim - im_rank}")

                # Show the 2-cycle
                if null_space is not None and null_space.shape[1] > 0:
                    cycle = null_space[:, 0]
                    nonzero = [(tt[i], cycle[i]) for i in range(len(tt)) if abs(cycle[i]) > 1e-8]
                    print(f"    2-cycle: {nonzero[:10]}")

    print(f"\n  ker(∂₂) dimension distribution:")
    for dim, cnt in sorted(cycle_counts.items()):
        print(f"    dim {dim}: {cnt} tournaments")

# ===== What does a typical 2-cycle look like at n=5? =====
print("\n" + "=" * 70)
print("DETAILED 2-CYCLE ANALYSIS AT n=5")
print("=" * 70)

n = 5
example_count = 0
for A in all_tournaments_gen(n):
    a1 = [(i,j) for i in range(n) for j in range(n) if i != j and A[i][j]]
    tt = transitive_triples(A, n)
    if len(tt) == 0:
        continue

    ker_dim, null_space = compute_ker_boundary2(tt, a1, n)

    if ker_dim >= 3 and example_count < 3:
        example_count += 1
        print(f"\n  Example {example_count}: ker_dim={ker_dim}, #TT={len(tt)}")

        scores = [sum(A[i]) for i in range(n)]
        print(f"    Scores: {scores}")

        # Count 3-cycles
        c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]))
        print(f"    c3 = {c3}")

        # Show one 2-cycle
        if null_space is not None and null_space.shape[1] > 0:
            cycle = null_space[:, 0]
            # Normalize to integers
            nonzero_vals = [abs(cycle[i]) for i in range(len(cycle)) if abs(cycle[i]) > 1e-8]
            if nonzero_vals:
                scale = min(nonzero_vals)
                cycle_int = [round(cycle[i]/scale) for i in range(len(cycle))]
                triples_in_cycle = [(tt[i], cycle_int[i]) for i in range(len(tt)) if cycle_int[i] != 0]
                print(f"    2-cycle ({len(triples_in_cycle)} triples):")
                for trip, coeff in triples_in_cycle[:20]:
                    print(f"      {coeff:+d} * {trip}")

        # DT check
        dt = dt_4paths(A, n)
        tt_idx = {t: i for i, t in enumerate(tt)}
        if len(dt) > 0:
            bd3 = np.zeros((len(tt), len(dt)))
            for j, (v0,v1,v2,v3) in enumerate(dt):
                faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
                signs = [1, -1, 1, -1]
                for face, sign in zip(faces, signs):
                    if face in tt_idx:
                        bd3[tt_idx[face], j] += sign
            im_rank = np.linalg.matrix_rank(bd3, tol=1e-8)
            print(f"    DT: {len(dt)}, im(∂₃|DT)={im_rank}, β₂={ker_dim-im_rank}")

# ===== Vertex participation analysis =====
print("\n" + "=" * 70)
print("VERTEX PARTICIPATION IN 2-CYCLES")
print("=" * 70)

n = 5
vertex_participation = Counter()
for A in all_tournaments_gen(n):
    a1 = [(i,j) for i in range(n) for j in range(n) if i != j and A[i][j]]
    tt = transitive_triples(A, n)
    if not tt:
        continue

    ker_dim, null_space = compute_ker_boundary2(tt, a1, n)
    if ker_dim == 0:
        continue

    # How many distinct vertices appear in the 2-cycle?
    if null_space is not None and null_space.shape[1] > 0:
        cycle = null_space[:, 0]
        verts_in_cycle = set()
        for i in range(len(tt)):
            if abs(cycle[i]) > 1e-8:
                verts_in_cycle.update(tt[i])
        vertex_participation[len(verts_in_cycle)] += 1

print("  Number of vertices in 2-cycles:")
for v_count, cnt in sorted(vertex_participation.items()):
    print(f"    {v_count} vertices: {cnt} tournaments")
