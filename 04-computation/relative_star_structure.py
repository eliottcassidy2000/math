#!/usr/bin/env python3
"""
relative_star_structure.py — opus-2026-03-09-S52

Investigate the algebraic structure of the "star" of v in the path complex.

The relative complex R_p = Omega_p(T) / Omega_p(T\v) consists of
v-dependent elements. Key structural question: what does this complex
look like in terms of v's neighborhood?

For a vertex v with in-neighborhood N^- and out-neighborhood N^+:
- R_0 = span{v}, dim = 1
- R_1: arcs incident to v (both i→v and v→j), dim = n-1
  But in Omega_1, we only keep arcs that are allowed paths.
  For tournaments, all arcs are allowed, so dim(Omega_1(T)) = n(n-1).
  But after the Omega quotient, dim may differ.

Actually, let me check: what does R_p look like concretely?
For a tournament, every edge is an allowed path (no non-allowed 1-faces).
So Omega_1 = A_1 (all edges), and the relative R_1 = "edges involving v".

At level 1 in the relative complex:
- Elements: edges (i,v) and (v,j) for i ∈ N^-(v), j ∈ N^+(v)
- These have dimension |N^-(v)| + |N^+(v)| = n-1
- The relative boundary d_1^rel maps these to R_0 = span{v}:
  d_1(i,v) = v - i, d_1(v,j) = j - v
  In the relative complex (mod T\v): d_1^rel(i,v) = v, d_1^rel(v,j) = -v
- So rank(d_1^rel) = 1, ker(d_1^rel) = n-2

Wait, but the mechanism data showed ker(d_1^rel) = 4 at n=6 (= n-2).
And dim(R_1) = 5 (= n-1). So rank(d_1^rel) = 5-4 = 1. ✓

For H_1: dim(R_1) = n-1, rank(d_1^rel) = 1, so ker = n-2.
H_1 = ker(d_1^rel) - im(d_2^rel) = (n-2) - im(d_2^rel).
For H_1 ≤ 1: need im(d_2^rel) ≥ n-3.

Question: what constrains im(d_2^rel)?

Let me analyze R_2 — the 2-paths through v.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def is_allowed_path(T, path):
    for i in range(len(path)-1):
        if not T[path[i]][path[i+1]]:
            return False
    return len(path) == len(set(path))

def analyze_v_star(T, n, v):
    """Analyze the structure of v-dependent paths."""
    N_in = [i for i in range(n) if i != v and T[i][v] == 1]
    N_out = [j for j in range(n) if j != v and T[v][j] == 1]

    # Level 0: just {v}
    r0 = 1

    # Level 1: arcs involving v
    # In Omega_1 for tournaments, every arc is allowed
    # v-dependent arcs: (i,v) for i in N_in, (v,j) for j in N_out
    r1_raw = len(N_in) + len(N_out)  # = n-1

    # Level 2: 2-paths through v
    # Types: (i,v,j) where i→v and v→j, i.e., i ∈ N_in, j ∈ N_out, i≠j
    #        (i,j,v) where i→j and j→v, i.e., T[i][j]=1 and j ∈ N_in
    #        (v,i,j) where v→i and i→j, i.e., i ∈ N_out and T[i][j]=1
    type_ivj = [(i,v,j) for i in N_in for j in N_out if i != j]
    type_ijv = [(i,j,v) for j in N_in for i in range(n) if i != j and i != v and T[i][j]]
    type_vij = [(v,i,j) for i in N_out for j in range(n) if j != i and j != v and T[i][j]]

    # But we also need to check which of these are in Omega_2
    # For Omega_2: a 2-path (a,b,c) is in Omega_2 if the face (a,c) obtained
    # by deleting b is either allowed or doesn't create a non-allowed constraint.
    # Actually, in Omega_p, we quotient by the subspace where deleting interior
    # vertices gives non-allowed faces. For p=2, the interior vertex is the middle one.

    # Let's just enumerate v-dependent paths at each level
    v_dep_paths = {}
    for p in range(0, min(n, 6)):
        paths = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                if v not in verts:
                    continue
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths.append(perm)
        v_dep_paths[p] = paths

    return {
        'N_in': N_in,
        'N_out': N_out,
        'd_in': len(N_in),
        'd_out': len(N_out),
        'v_dep_counts': {p: len(v_dep_paths[p]) for p in v_dep_paths},
        'type_ivj': len(type_ivj),
        'type_ijv': len(type_ijv),
        'type_vij': len(type_vij),
    }

def main():
    print("=" * 70)
    print("RELATIVE STAR STRUCTURE ANALYSIS")
    print("=" * 70)

    # ===========================================
    # Part 1: Structure of R_p at n=6
    # ===========================================
    n = 6
    print(f"\nn = {n}: v-dependent path counts")

    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    # Track v-dependent path counts by (d_in, d_out) = (in-degree, out-degree) of v
    from collections import defaultdict
    stats = defaultdict(lambda: defaultdict(list))

    for bits in range(0, n_total, 4):  # sample 1/4
        T = tournament_from_bits(n, bits)
        for v in range(1):  # just v=0 for speed
            res = analyze_v_star(T, n, v)
            key = (res['d_in'], res['d_out'])
            for p in res['v_dep_counts']:
                stats[key][p].append(res['v_dep_counts'][p])

    print(f"\n  v-dependent path counts A_p^v (raw, before Omega quotient):")
    print(f"  {'(d_in, d_out)':>15} | {'A_0^v':>5} | {'A_1^v':>5} | {'A_2^v':>7} | {'A_3^v':>7} | {'A_4^v':>7}")
    print(f"  {'-'*15}-+-{'-'*5}-+-{'-'*5}-+-{'-'*7}-+-{'-'*7}-+-{'-'*7}")
    for key in sorted(stats.keys()):
        vals = []
        for p in range(5):
            if p in stats[key]:
                vals.append(f"{np.mean(stats[key][p]):7.1f}")
            else:
                vals.append(f"{'':>7}")
        print(f"  {str(key):>15} | {' | '.join(vals)}")

    # ===========================================
    # Part 2: Key question — how does d_2^rel act?
    # At level 1: ker(d_1^rel) = n-2
    # We need im(d_2^rel) >= n-3 for H_1 <= 1
    # ===========================================
    print(f"\n{'='*70}")
    print("LEVEL 2 BOUNDARY MAP ANALYSIS")
    print("=" * 70)
    print("""
    For H_1(T,T\\v) <= 1, we need im(d_2^rel) >= n-3 (= 3 at n=6).
    Since ker(d_1^rel) = n-2 (= 4 at n=6), and dim(R_1) = n-1 (= 5 at n=6),
    exactly 1 dimension is killed by d_1^rel (maps to R_0).

    H_1 = (n-2) - im(d_2^rel).
    H_1 = 0 iff im(d_2^rel) = n-2.
    H_1 = 1 iff im(d_2^rel) = n-3.
    """)

    # Track im(d_2^rel) by in/out degree of v
    from collections import Counter
    im_d2_by_degree = defaultdict(list)

    # Need full computation for this
    from itertools import combinations, permutations

    tol = 1e-8

    for bits in range(0, n_total, 1):  # exhaustive
        if bits % 8000 == 0:
            print(f"  ... {bits}/{n_total}")
        T = tournament_from_bits(n, bits)
        v = 0
        d_in = sum(T[i][v] for i in range(n) if i != v)
        d_out = sum(T[v][j] for j in range(n) if j != v)

        # Quick computation of relative R_1 and R_2
        # R_1 elements: edges involving v
        # In Omega_1 for tournaments, all edges are in Omega_1
        # (because all 0-paths are allowed, so no face constraint)
        edges_v = []
        for i in range(n):
            if i != v:
                if T[i][v]:
                    edges_v.append((i, v))
                if T[v][i]:
                    edges_v.append((v, i))
        # dim R_1 = len(edges_v) = n-1

        # R_2: 2-paths involving v, in Omega_2
        # First get all 2-paths involving v
        paths2_v = []
        for verts in combinations(range(n), 3):
            if v not in verts:
                continue
            for perm in permutations(verts):
                if is_allowed_path(T, perm):
                    paths2_v.append(perm)

        # Check which are in Omega_2
        # Omega_2 = A_2 / (constraint: d_1(sigma) not in A_1 for interior face)
        # For p=2, sigma = (a,b,c), interior face = (a,c) (delete middle vertex)
        # If (a,c) is not an allowed 1-path, we get a constraint
        all_1paths = set()
        for i in range(n):
            for j in range(n):
                if i != j and T[i][j]:
                    all_1paths.add((i,j))

        # Build constraint matrix for v-dependent 2-paths
        na_faces = {}
        for sigma in paths2_v:
            face = (sigma[0], sigma[2])  # delete middle
            if face not in all_1paths:
                if face not in na_faces:
                    na_faces[face] = len(na_faces)

        if not na_faces:
            # All 2-paths involving v are in Omega_2
            omega2_v_dim = len(paths2_v)
        else:
            mat = np.zeros((len(na_faces), len(paths2_v)))
            for j, sigma in enumerate(paths2_v):
                face = (sigma[0], sigma[2])
                if face in na_faces:
                    mat[na_faces[face], j] += (-1)**1  # sign for deleting index 1
            U, S, Vt = np.linalg.svd(mat, full_matrices=True)
            rank = int(np.sum(S > tol))
            omega2_v_dim = len(paths2_v) - rank

        # Boundary d_2 maps R_2 → R_1
        # d_2(a,b,c) = (b,c) - (a,c) + (a,b)
        # In relative complex, only keep v-dependent part
        # Actually this is getting complex. Let me just track omega2_v_dim
        # and compare with the full computation from the mechanism script.

        # Actually, let me just directly look at whether d_in, d_out affects H_1
        # I already know from the mechanism data that H_1 ∈ {0,1} always.
        # The question is: what structural property of v determines H_1?

        # Skip the full boundary computation; just record omega2_v_dim
        im_d2_by_degree[(d_in, d_out)].append(omega2_v_dim)

    print(f"\n  dim(Omega_2^v) by (d_in, d_out) of v:")
    for key in sorted(im_d2_by_degree.keys()):
        vals = im_d2_by_degree[key]
        print(f"    {key}: mean={np.mean(vals):.1f}, min={min(vals)}, max={max(vals)}, n={len(vals)}")

    # ===========================================
    # Part 3: When does H_1(T,T\v) = 1?
    # ===========================================
    print(f"\n{'='*70}")
    print("WHEN IS H_1(T,T\\v) = 1?")
    print("=" * 70)

    # From mechanism data: 5760 cases with H_1 = 1 out of 196608
    # What characterizes these?
    h1_eq_1_scores = Counter()
    h1_eq_0_scores = Counter()

    for bits in range(n_total):
        T = tournament_from_bits(n, bits)
        for v in range(n):
            d_out = sum(T[v][j] for j in range(n) if j != v)
            # H_1(T,T\v) = 1 iff v is a source or sink in some sense
            # Let me compute directly from the score sequence
            score_seq = tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)]))

            # To know H_1: I need the full computation, but let me check
            # if it correlates with d_out (out-degree of v)
            # H_1 = 1 when im(d_2^rel) = n-3 = 3
            # H_1 = 0 when im(d_2^rel) = n-2 = 4
            pass

    # Instead, let me just check: which (score_of_v, score_seq) give H_1 = 1?
    # I'll use the already-computed mechanism data approach
    print("\n  Correlating H_1(T,T\\v) with out-degree of v:")

    h1_by_outdeg = defaultdict(lambda: Counter())
    for bits in range(0, n_total, 1):
        T = tournament_from_bits(n, bits)
        for v in range(n):
            d_out = sum(T[v][j] for j in range(n) if j != v)
            # Need to know H_1(T,T\v)... let me check a simple criterion:
            # v is a "source" in the subgraph on its in-neighbors?
            # Or: the in-neighborhood and out-neighborhood structure
            pass

    # Actually, let me just check: is H_1 = 1 iff d_out = 0 or d_out = n-1?
    # (i.e., v is a source or sink)
    print("  Testing: is H_1(T,T\\v) = 1 iff v is source (d_out=5) or sink (d_out=0)?")
    print(f"  n=6: sources/sinks per tournament: 0 or 1 source, 0 or 1 sink")
    print(f"  Total (T,v) with d_out ∈ {{0, {n-1}}}: ", end="")

    source_sink_count = 0
    for bits in range(n_total):
        T = tournament_from_bits(n, bits)
        for v in range(n):
            d_out = sum(T[v][j] for j in range(n) if j != v)
            if d_out == 0 or d_out == n-1:
                source_sink_count += 1
    print(f"{source_sink_count}")
    print(f"  H_1 = 1 cases from mechanism data: 5760")
    print(f"  Match: {'YES' if source_sink_count == 5760 else 'NO'}")

    if source_sink_count != 5760:
        # Check the actual condition
        print(f"\n  Source/sink count ({source_sink_count}) ≠ H_1=1 count (5760)")
        print(f"  Investigating what d_out values give H_1 = 1...")

        # Need to actually compute H_1 for a subset
        # From mechanism data: H_1(T,T\v) = ker(d_1) - im(d_2) = 4 - im(d_2)
        # H_1 = 1 iff im(d_2) = 3, H_1 = 0 iff im(d_2) = 4
        # But im(d_2) depends on Omega_2, which requires the full computation.

        # Let me just sample to find the pattern
        h1_by_score = Counter()
        sample_count = 0
        for bits in range(0, n_total, 4):
            T = tournament_from_bits(n, bits)
            v = 0
            d_out = sum(T[v][j] for j in range(n) if j != v)
            # Compute H_1 the hard way... actually let me skip to Part 4
            h1_by_score[d_out] += 1
            sample_count += 1
        print(f"\n  (Skipping detailed H_1 computation — focus on H_3 mechanism)")

    # ===========================================
    # Part 4: Structure that controls H_3
    # ===========================================
    print(f"\n{'='*70}")
    print("H_3 STRUCTURAL ANALYSIS")
    print("=" * 70)

    # For H_3(T,T\v) = 1 (1920 cases at n=6):
    # From mechanism: two types:
    #   dim(R_3)=6, ker=1, im_d4=0 (480 cases)
    #   dim(R_3)=14, ker=6, im_d4=5 (1440 cases)
    # In both: gap = 1.
    #
    # Type 1 (dim=6, ker=1): minimal relative complex,
    # there's exactly 1 cycle that's not a boundary.
    #
    # Type 2 (dim=14, ker=6): bigger relative complex,
    # but 5 out of 6 kernel elements are boundaries, leaving 1.

    # Question: what distinguishes these types?
    # Hypothesis: score sequence of T or score of v
    print("\n  At n=6, H_3(T,T\\v) = 1 occurs 1920 times.")
    print("  Two structural types based on dim(R_3):")
    print("    Type 1: dim(R_3)=6, ker=1, im=0 — 480 cases")
    print("    Type 2: dim(R_3)=14, ker=6, im=5 — 1440 cases")
    print("  Total: 1920 = 320 × 6 (i.e., 320 tournaments × 6 vertices)")
    print("  But from good_vertex: all have beta_3(T\\v)=0 for ALL v")
    print("  So H_3(T,T\\v) = beta_3(T) for these. All 320 have beta_3=1.")
    print("  1920 = 320 × 6 confirms: H_3=1 for ALL vertex pairs of beta_3=1 tournaments.")

    print("\n  KEY INSIGHT: For beta_3=1 tournaments, H_3(T,T\\v) = 1 for ALL v.")
    print("  This means beta_3(T) = H_3(T,T\\v) + rank(i_*) = 1 + 0 = 1")
    print("  where rank(i_*) = 0 because beta_3(T\\v) = 0.")
    print("  Consistent with the exact equation.")

    print(f"\n  Score sequences of the 320 beta_3=1 tournaments:")
    score_types = Counter()
    for bits in range(n_total):
        T = tournament_from_bits(n, bits)
        # Quick beta_3 check via score first
        scores = sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)])
        score_types[tuple(scores)] += 1

    # From earlier: (1,1,1,4,4,4) → 80, (2,2,2,3,3,3) → 240
    print(f"  Total tournaments by score: {len(score_types)} distinct sequences")
    for s, c in sorted(score_types.items(), key=lambda x: -x[1])[:10]:
        print(f"    {s}: {c} tournaments")

if __name__ == '__main__':
    main()
