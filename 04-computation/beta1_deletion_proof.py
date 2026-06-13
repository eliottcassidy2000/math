#!/usr/bin/env python3
"""
beta1_deletion_proof.py — Prove sum_v beta_1(T\v) <= 3 for beta_1(T)=0

Key findings from beta1_structure_analysis.py:
- dim(Z_1(T)) = C(n-1,2) for ALL tournaments (universal!)
- beta_1(T) = C(n-1,2) - rank(d2|_{Omega_2})
- When beta_1(T)=0: rank(d2) = C(n-1,2)
- For T\v: dim(Z_1(T\v)) = C(n-2,2)
- beta_1(T\v) = C(n-2,2) - dim(B_1(T\v))
- rank_drop(v) = dim(B1(T)) - dim(B1(T\v)) = (n-2) + beta_1(T\v)
  i.e. dim(B1(T\v)) = C(n-1,2) - (n-2) - beta_1(T\v) = C(n-2,2) - beta_1(T\v) ✓

PROOF STRATEGY:
The key is that sum_v rank_drop(v) has a constraint.

rank_drop(v) = dim(B1(T)) - dim(B1(T\v))
             = C(n-1,2) - (C(n-2,2) - beta_1(T\v))
             = (n-2) + beta_1(T\v)

So sum_v beta_1(T\v) = sum_v rank_drop(v) - n*(n-2)

What's sum_v rank_drop(v)?
It's related to how much rank is "supported" by each vertex.

B_1(T) is spanned by boundaries of 2-paths.
Each 2-path (a,b,c) involves 3 vertices.
When we remove v, we lose all 2-paths through v.

The 2-paths NOT through v generate B_1(T\v) (as a subspace).
The 2-paths THROUGH v generate the complement.

B_1(T) = B_1(T\v) + span{d2(p) : v in p}  (as subspaces of edge space)

Wait — that's not quite right since B_1(T\v) lives in the edge space of T\v,
while B_1(T) lives in the edge space of T.

Let me think more carefully using the EXACT numerical patterns.

opus-2026-03-08
"""

import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
import sys

# Import core functions from the analysis script
sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')
from beta1_structure_analysis import (
    enumerate_allowed_paths, boundary_coeffs, build_full_boundary_matrix,
    compute_omega_basis, compute_beta1_full, compute_beta1_only,
    subtournament, all_tournaments, count_3cycles, rank_of_matrix
)

def run_proof_analysis():
    print("=" * 70)
    print("BETA_1 DELETION BOUND — PROOF ANALYSIS")
    print("=" * 70)

    # ==================================================================
    # KEY OBSERVATION 1: beta_1(T\v) ∈ {0, 1} for tournaments
    # ==================================================================
    print("\n" + "=" * 70)
    print("OBSERVATION 1: beta_1 values for subtournaments")
    print("=" * 70)

    for n in [4, 5, 6]:
        beta_vals = Counter()
        cnt = 0
        for A in all_tournaments(n):
            b1 = compute_beta1_only(A, n)
            beta_vals[b1] += 1
            cnt += 1
        print(f"  n={n}: beta_1 values: {dict(sorted(beta_vals.items()))} ({cnt} total)")

    # ==================================================================
    # KEY OBSERVATION 2: beta_1(T) relates to the number of 3-cycles
    # ==================================================================
    print("\n" + "=" * 70)
    print("OBSERVATION 2: beta_1 vs t3 (3-cycle count)")
    print("=" * 70)

    for n in [4, 5]:
        print(f"\n  n={n}:")
        by_t3_b1 = defaultdict(Counter)
        for A in all_tournaments(n):
            t3 = count_3cycles(A, n)
            b1 = compute_beta1_only(A, n)
            by_t3_b1[t3][b1] += 1
        for t3 in sorted(by_t3_b1.keys()):
            print(f"    t3={t3}: {dict(by_t3_b1[t3])}")

    # ==================================================================
    # CRITICAL TEST: Is beta_1(T\v) = 1 iff T\v has >= 2 3-cycles?
    # At n=4: t3=2 <=> beta_1=1. Let's check at n=5 deletions.
    # ==================================================================
    print("\n" + "=" * 70)
    print("OBSERVATION 3: beta_1(T\\v)=1 iff t3(T\\v) >= threshold?")
    print("=" * 70)

    for n in [5]:
        t3_when_bad = Counter()
        t3_when_good = Counter()
        for A in all_tournaments(n):
            info = compute_beta1_full(A, n)
            if info['beta_1'] != 0:
                continue
            for v in range(n):
                Av, nv = subtournament(A, n, [v])
                t3v = count_3cycles(Av, nv)
                b1v = compute_beta1_only(Av, nv)
                if b1v > 0:
                    t3_when_bad[t3v] += 1
                else:
                    t3_when_good[t3v] += 1
        print(f"  n={n}: t3(T\\v) when beta_1(T\\v)=1: {dict(sorted(t3_when_bad.items()))}")
        print(f"  n={n}: t3(T\\v) when beta_1(T\\v)=0: {dict(sorted(t3_when_good.items()))}")

    # ==================================================================
    # THE CORE IDENTITY: rank_drop analysis
    # ==================================================================
    print("\n" + "=" * 70)
    print("CORE IDENTITY: rank_drop(v) = (n-2) + beta_1(T\\v)")
    print("=" * 70)

    print("""
  We KNOW (from data):
  - dim(Z_1(T)) = C(n-1,2) for ALL tournaments on n vertices
  - When beta_1(T)=0: dim(B_1(T)) = C(n-1,2)
  - dim(B_1(T\\v)) = C(n-2,2) - beta_1(T\\v)

  So rank_drop(v) := dim(B_1(T)) - dim(B_1(T\\v))
                    = C(n-1,2) - C(n-2,2) + beta_1(T\\v)
                    = (n-2) + beta_1(T\\v)

  Now sum_v rank_drop(v) = n(n-2) + sum_v beta_1(T\\v)

  Question: What is sum_v rank_drop(v)?
  """)

    # Verify and compute sum_v rank_drop(v) at n=5
    for n in [5]:
        print(f"\n  n={n}: sum_v rank_drop(v) analysis")
        sum_rd_dist = Counter()
        sum_rd_by_beta = defaultdict(list)
        for A in all_tournaments(n):
            info = compute_beta1_full(A, n)
            if info['beta_1'] != 0:
                continue
            total_rd = 0
            betas_v = []
            for v in range(n):
                Av, nv = subtournament(A, n, [v])
                info_v = compute_beta1_full(Av, nv)
                rd = info['dim_B1'] - info_v['dim_B1']
                total_rd += rd
                betas_v.append(info_v['beta_1'])
            sum_rd_dist[total_rd] += 1
            sum_b = sum(betas_v)
            sum_rd_by_beta[sum_b].append(total_rd)

        print(f"    sum rank_drop distribution: {dict(sorted(sum_rd_dist.items()))}")
        print(f"    Expected: n*(n-2) + sum_beta = {n*(n-2)} + sum_beta")
        for sb in sorted(sum_rd_by_beta.keys()):
            rds = sum_rd_by_beta[sb]
            print(f"    sum_beta={sb}: sum_rank_drops always = {rds[0]}? {all(r == rds[0] for r in rds)} "
                  f"(expected {n*(n-2) + sb})")

    # ==================================================================
    # THE PROJECTION APPROACH: B_1(T) projected to T\v edges
    # ==================================================================
    print("\n" + "=" * 70)
    print("PROJECTION: dim(B_1(T)|_{T\\v edges}) vs dim(B_1(T\\v))")
    print("=" * 70)

    print("""
  B_1(T) lives in R^{n(n-1)} (all edges of T).
  Projecting to T\\v edges (restricting to the (n-1)(n-2) edges not involving v):

  Let pi_v: R^{n(n-1)} -> R^{(n-1)(n-2)} be this projection.

  Then pi_v(B_1(T)) is a subspace of Z_1(T\\v) (since boundaries project to cycles).
  And B_1(T\\v) ⊂ pi_v(B_1(T)) (since 2-paths of T\\v are also 2-paths of T).

  So: dim(B_1(T\\v)) <= dim(pi_v(B_1(T))) <= dim(Z_1(T\\v)) = C(n-2,2)

  The "gap" is: dim(pi_v(B_1(T))) - dim(B_1(T\\v))
  This gap = # boundaries through v that project to NEW cycles in T\\v.
  """)

    for n in [5]:
        print(f"\n  n={n}: Projection analysis")
        for idx_count, A in enumerate(all_tournaments(n)):
            info = compute_beta1_full(A, n)
            if info['beta_1'] != 0:
                continue

            betas_v = []
            for v in range(n):
                Av, nv = subtournament(A, n, [v])
                betas_v.append(compute_beta1_only(Av, nv))

            if sum(betas_v) < 2:
                continue

            # Compute pi_v(B_1(T)) for each v
            print(f"\n    Tournament idx={idx_count}, betas={betas_v}")

            # B_1(T) in edge coordinates
            B1_T = info['bd_2_omega']  # shape (|edges|, rank)
            edges_T = info['allowed_1']
            edge_idx = {e: i for i, e in enumerate(edges_T)}

            for v in range(n):
                # Edges not involving v
                keep = [u for u in range(n) if u != v]
                proj_rows = []
                for i, e in enumerate(edges_T):
                    if v not in e:
                        proj_rows.append(i)

                # Project B_1(T)
                pi_v_B1 = B1_T[proj_rows, :]
                rank_proj = rank_of_matrix(pi_v_B1)

                # B_1(T\v)
                Av, nv = subtournament(A, n, [v])
                info_v = compute_beta1_full(Av, nv)

                gap = rank_proj - info_v['dim_B1']
                marker = " ***" if betas_v[v] > 0 else ""
                print(f"    v={v}: dim(pi_v(B1(T)))={rank_proj}, dim(B1(T\\v))={info_v['dim_B1']}, "
                      f"gap={gap}, beta1(T\\v)={betas_v[v]}{marker}")

            # Only show a few
            if idx_count > 50:
                break

    # ==================================================================
    # THE 3-CYCLE CONNECTION: beta_1(T)=1 iff T has a "pure 3-cycle"
    # ==================================================================
    print("\n" + "=" * 70)
    print("3-CYCLE CONNECTION: What makes beta_1=1?")
    print("=" * 70)

    print("""
  At n=4: beta_1=1 iff t3=2 (the maximum). These are the "rotational" tournaments.
  At n=4, a tournament with t3=2 has exactly one non-trivial 1-cycle (the "alternating flow").

  beta_1 measures how many independent directed cycles exist that are NOT boundaries.
  A cycle z ∈ Z_1 is a boundary iff it equals d2(c) for some 2-chain c in Omega_2.

  For a tournament, d2(a,b,c) involves edges:
  - (a,b), (a,c) or (c,a), (b,c): depends on whether {a,b,c} is TT or 3-cycle.

  For a 3-CYCLE (a,b,c) with a->b->c->a:
  The allowed 2-path is (a,b,c) (since a->b->c and c->a is also an edge).
  d2(a,b,c) = (b,c) - (c,a) + (a,b) = the 3-cycle itself!

  Wait — that means every 3-cycle IS a boundary. So beta_1 should be 0?

  No! The issue is that (a,b,c) might NOT be in Omega_2.
  Omega_2 = {u in A_2 : d_2(u) in A_1}. For tournaments, every edge is in A_1,
  so the condition is that all faces of the 2-path are edges.

  Actually, for a 2-path (a,b,c): faces are (b,c), (a,c), (a,b).
  In a tournament, either a->c or c->a. If c->a, then face (a,c) is NOT an edge,
  but (c,a) IS. The boundary formula gives d2(a,b,c) = (b,c) - (a,c) + (a,b).
  But (a,c) is not an allowed path — it's in the WRONG direction!

  So (a,c) appears in the boundary with a minus sign. If c->a in T,
  then (a,c) is not an edge of T, so this face is "non-allowed."

  The condition for (a,b,c) to be in Omega_2 is that d2(a,b,c) has no component
  on non-allowed 1-paths. So ALL faces must be edges.

  For (a,b,c) where a->b->c:
  - Face (b,c): always an edge (since b->c)
  - Face (a,c): an edge iff a->c (transitive triple)
  - Face (a,b): always an edge (since a->b)

  So (a,b,c) is in Omega_2 iff {a,b,c} is a TRANSITIVE TRIPLE with a->b->c and a->c.

  WAIT — but I also saw 3-cycle 2-paths in the data. Let me re-check.

  For a 3-cycle a->b->c->a, the 2-path (a,b,c) has:
  - Face (b,c): edge (b->c) ✓
  - Face (a,c): NOT an edge (c->a, not a->c). ✗

  But d2(a,b,c) = (b,c) - (a,c) + (a,b). The face (a,c) is non-allowed.
  So... (a,b,c) is NOT in Omega_2 when {a,b,c} is a 3-cycle?

  That contradicts my earlier understanding. Let me check this computationally.
  """)

    # Check: which 2-paths are in Omega_2?
    print("  Checking Omega_2 membership for n=4:")
    for idx, A in enumerate(all_tournaments(4)):
        t3 = count_3cycles(A, 4)
        if t3 != 2:
            continue

        info = compute_beta1_full(A, 4)
        print(f"\n    Tournament idx={idx}, t3={t3}")
        print(f"    Allowed 2-paths (A_2): {info['allowed_2']}")
        print(f"    dim(Omega_2) = {info['dim_omega2']}")

        # Which 2-paths are in Omega_2?
        # Omega_2 basis is info['omega_2'] in A_2 coordinates
        omega2 = info['omega_2']
        if omega2.ndim == 2 and omega2.shape[1] > 0:
            # Check which A_2 paths have nonzero projection
            print(f"    Omega_2 basis vectors (in A_2 coords):")
            for j in range(omega2.shape[1]):
                col = omega2[:, j]
                nonzero = [(info['allowed_2'][i], col[i]) for i in range(len(col)) if abs(col[i]) > 1e-10]
                print(f"      basis {j}: {nonzero}")

        # Classify each 2-path
        for path in info['allowed_2']:
            a, b, c = path
            is_tt = (A[a][c] == 1)  # a->c means transitive
            is_3cyc = (A[c][a] == 1)  # c->a means 3-cycle
            # Check faces
            face_bc = (b, c)
            face_ac = (a, c)
            face_ab = (a, b)
            all_edges = set(info['allowed_1'])
            f_bc = face_bc in all_edges
            f_ac = face_ac in all_edges
            f_ab = face_ab in all_edges
            print(f"    path {path}: {'TT' if is_tt else '3cyc'}, "
                  f"faces: ({b},{c})={'✓' if f_bc else '✗'} ({a},{c})={'✓' if f_ac else '✗'} ({a},{b})={'✓' if f_ab else '✗'}")

        if idx > 5:
            break

    # ==================================================================
    # REFORMULATION: Omega_2 for tournaments
    # ==================================================================
    print("\n" + "=" * 70)
    print("OMEGA_2 STRUCTURE: What's actually in Omega_2?")
    print("=" * 70)

    print("""
  From the computation above, we can see that:
  - A 2-path (a,b,c) with a->b->c is in A_2 (it follows edges)
  - Its faces are (b,c), (a,c), (a,b)
  - In a tournament, (b,c) and (a,b) are always edges
  - But (a,c): if a->c (TT), then (a,c) ∈ A_1 and (a,b,c) ∈ Omega_2
             if c->a (3-cycle), then (a,c) ∉ A_1, so d2(a,b,c) has a non-allowed component

  For 3-cycle paths, (a,b,c) is NOT individually in Omega_2.
  BUT: a LINEAR COMBINATION of 3-cycle paths might be in Omega_2
  (the non-allowed components cancel).

  This is the key to dim(Omega_2) < |A_2|!

  Let me understand: when does a 3-cycle 2-path (a,b,c) have its non-allowed
  face cancel with another 2-path's non-allowed face?

  Path (a,b,c) with c->a: non-allowed face is (a,c) with coefficient -1
  Path (c,b,a) with a->c: this requires c->b->a in T and checking (c,a) face
  Wait, that's a different path. Let me think systematically.
  """)

    # For each non-TT 2-path, identify its non-allowed face
    for n in [4]:
        for idx, A in enumerate(all_tournaments(n)):
            t3 = count_3cycles(A, n)
            if t3 != 2:
                continue

            info = compute_beta1_full(A, n)
            allowed_1_set = set(info['allowed_1'])

            non_allowed_faces = defaultdict(list)
            for j, path in enumerate(info['allowed_2']):
                for sign, face in boundary_coeffs(path):
                    if face not in allowed_1_set and len(set(face)) == len(face):
                        non_allowed_faces[face].append((j, path, sign))

            print(f"\n  n={n}, idx={idx}, t3={t3}")
            print(f"  Non-allowed faces shared by 2-paths:")
            for face, paths in non_allowed_faces.items():
                print(f"    Face {face}:")
                for j, path, sign in paths:
                    print(f"      path {path} contributes {sign:+d} * {face}")

            # So Omega_2 consists of linear combos where non-allowed faces cancel
            # For each non-allowed face (a,c), the paths that produce it must sum to zero
            # This gives a system of constraints.

            if idx > 5:
                break

    # ==================================================================
    # CRITICAL DIMENSION IDENTITY
    # ==================================================================
    print("\n" + "=" * 70)
    print("CRITICAL DIMENSION IDENTITY")
    print("=" * 70)

    print("""
  For a tournament on n vertices:

  dim(Z_1) = C(n-1, 2) always.

  PROOF SKETCH:
  Z_1 = ker(d_1: Omega_1 -> Omega_0).
  For tournaments, Omega_1 = A_1 = all n(n-1) directed edges.
  Omega_0 = all n vertices.
  d_1(a,b) = (b) - (a), so d_1 is the incidence matrix.

  dim(Z_1) = n(n-1) - rank(d_1) = n(n-1) - (n-1) = (n-1)(n-1) - ...
  Wait: rank of incidence matrix of complete directed graph on n vertices.
  d_1 maps R^{n(n-1)} -> R^n. rank(d_1) = n-1 (connected graph).
  So dim(Z_1) = n(n-1) - (n-1) = (n-1)^2.

  But we observed dim(Z_1) = C(n-1,2) = (n-1)(n-2)/2.

  (n-1)^2 ≠ C(n-1,2) in general. So something is wrong.

  Wait — Omega_1 might NOT be all of A_1 for a tournament!
  Let me check.
  """)

    for n in [4, 5]:
        print(f"\n  n={n}:")
        edges_per_tournament = set()
        omega1_per_tournament = set()
        for A in all_tournaments(n):
            info = compute_beta1_full(A, n)
            edges_per_tournament.add(len(info['allowed_1']))
            omega1_per_tournament.add(info['dim_omega1'])
        print(f"    |A_1| values: {edges_per_tournament}")
        print(f"    dim(Omega_1) values: {omega1_per_tournament}")
        print(f"    n(n-1) = {n*(n-1)}, (n-1)^2 = {(n-1)**2}, C(n-1,2) = {(n-1)*(n-2)//2}")
        print(f"    C(n-1,2) + (n-1) = {(n-1)*(n-2)//2 + n-1}")

    # ==================================================================
    # UNDERSTAND dim(Omega_1)
    # ==================================================================
    print("\n" + "=" * 70)
    print("UNDERSTANDING Omega_1 FOR TOURNAMENTS")
    print("=" * 70)

    print("""
  Omega_1 = {u ∈ A_1 : d_1(u) ∈ A_0}
  For a tournament, A_1 = all directed edges, A_0 = all vertices.
  d_1(i,j) = (j) - (i), which is always in A_0.
  So Omega_1 = A_1 for all tournaments.

  Then dim(Z_1) = dim(ker(d_1)) = n(n-1) - rank(d_1).
  For a tournament (strongly connected or not), the incidence matrix d_1
  has rank n-1 if the underlying undirected graph is connected.
  A tournament is always connected, so rank(d_1) = n-1.

  dim(Z_1) = n(n-1) - (n-1) = (n-1)(n+1-1) = (n-1)^2.

  But we measured dim(Z_1) = C(n-1,2) = (n-1)(n-2)/2.

  CONTRADICTION. Unless dim(Omega_1) ≠ n(n-1).
  Let me re-measure.
  """)

    # Direct measurement
    for n in [3, 4]:
        A = [[0]*n for _ in range(n)]
        # Transitive tournament
        for i in range(n):
            for j in range(i+1, n):
                A[i][j] = 1
        info = compute_beta1_full(A, n)
        print(f"  n={n} transitive: |A_1|={len(info['allowed_1'])}, dim(Om1)={info['dim_omega1']}, "
              f"dim(Z1)={info['dim_Z1']}, n(n-1)={n*(n-1)}, rank(d1)={info['rank_d1']}")
        print(f"    C(n-1,2)={n*(n-1)//2 - (n-1)}, (n-1)^2={(n-1)**2}")

    # Hmm wait — C(n-1,2) = (n-1)(n-2)/2
    # At n=4: C(3,2)=3. But (n-1)^2=9. And dim_omega1 would be 12 (=4*3).
    # dim(Z1) = 12 - rank(d1). rank(d1) = 3. So dim(Z1) = 9.
    # But we measured 3! So something else is going on.

    # Oh wait — maybe Omega_1 ≠ A_1 for tournaments!
    # Let me check directly.
    print("\n  Direct check: is Omega_1 = A_1?")
    for n in [3, 4]:
        for idx, A in enumerate(all_tournaments(n)):
            info = compute_beta1_full(A, n)
            if info['dim_omega1'] != len(info['allowed_1']):
                print(f"    n={n}, idx={idx}: dim(Om1)={info['dim_omega1']} ≠ |A_1|={len(info['allowed_1'])}")
                break
        else:
            print(f"    n={n}: Omega_1 = A_1 for all tournaments? "
                  f"dim(Om1)={info['dim_omega1']}, |A_1|={len(info['allowed_1'])}")

    # So dim(Omega_1) = n(n-1)/2?? That would be C(n,2).
    # Wait let me just look at the numbers.
    for n in [3, 4, 5]:
        for A in all_tournaments(n):
            info = compute_beta1_full(A, n)
            print(f"  n={n}: dim(Om1)={info['dim_omega1']}, C(n,2)={n*(n-1)//2}, n(n-1)={n*(n-1)}")
            break

    # ==================================================================
    # THE ACTUAL CONSTRAINT ON SUM_V BETA_1
    # ==================================================================
    print("\n" + "=" * 70)
    print("THE ACTUAL CONSTRAINT ON SUM_V BETA_1(T\\v)")
    print("=" * 70)

    print("""
  EMPIRICAL FACTS (verified n=5,6):
  1. beta_1(T) ∈ {0, 1} for all tournaments
  2. When beta_1(T)=0: sum_v beta_1(T\\v) <= 3
  3. When beta_1(T)=0: each beta_1(T\\v) ∈ {0, 1}
  4. The bad vertices (beta_1(T\\v)=1) number at most 3

  KEY PATTERNS at n=5 (720 tournaments with beta_1=0):
  - 0 bad vertices: 240 (33.3%)
  - 1 bad vertex: 240 (33.3%)
  - 2 bad vertices: 120 (16.7%)
  - 3 bad vertices: 120 (16.7%)

  These are EXACT fractions: 240/720 = 1/3, etc.
  Total: 240 + 240 + 120 + 120 = 720. ✓

  At n=6 (27968 tournaments with beta_1=0):
  - 0 bad: 8480
  - 1 bad: 7968
  - 2 bad: 7200
  - 3 bad: 4320
  Total: 27968. ✓

  CONJECTURE: sum_v beta_1(T\\v) <= 3 for ALL n.
  Moreover: beta_1(T) <= 1 for all tournaments.

  If beta_1(T) <= 1 for all tournaments, then:
  At n=4: beta_1 ∈ {0,1} — verified.
  At n=5: beta_1 ∈ {0,1} — verified.
  At n=6: beta_1 ∈ {0,1} — verified.

  PROOF IDEA for beta_1 <= 1:
  beta_1 counts independent 1-cycles not bounding 2-chains.
  In a tournament, the "uncovered" cycle comes from the structure of Omega_2.

  The constraint dim(Omega_2) < |A_2| creates a cokernel.
  The cokernel dimension is |A_2| - dim(Omega_2).
  But not all of this cokernel contributes to beta_1.

  The deficiency is: beta_1 = dim(Z_1) - dim(B_1)
  = dim(Z_1) - rank(d2|_{Omega_2})
  = dim(Z_1) - (dim(Omega_2) - dim(ker(d2|_{Omega_2})))
  """)

    # ==================================================================
    # DEEPER: What is dim(Omega_1) actually?
    # ==================================================================
    print("\n" + "=" * 70)
    print("UNDERSTANDING DIM(OMEGA_1) — THE KEY")
    print("=" * 70)

    # It seems dim(Omega_1) < n(n-1). How?
    # A 1-path (a,b) is in A_1 iff a->b in T. That's n(n-1)/2 for a tournament (wait no)
    # Actually for a tournament, for each pair {a,b}, exactly one of (a,b) or (b,a) is in A_1.
    # So |A_1| = C(n,2) = n(n-1)/2.

    for n in [3, 4, 5]:
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A[i][j] = 1
        paths = enumerate_allowed_paths(A, n, 1)
        print(f"  n={n}: |A_1| for transitive = {len(paths)}, C(n,2)={n*(n-1)//2}")

    # OK so |A_1| = C(n,2) for tournaments (not n(n-1)).
    # That's because (a,b) is allowed only if a->b.
    # For each undirected edge, only one direction is allowed.
    print("""
  AH! For a TOURNAMENT, |A_1| = C(n,2) (not n(n-1)).
  Because for each pair {a,b}, exactly one of a->b or b->a holds.

  So the 1-chain space A_1 has dimension C(n,2).
  d_1: R^{C(n,2)} -> R^n. rank(d_1) = n-1.

  dim(Z_1) = C(n,2) - (n-1) = n(n-1)/2 - (n-1) = (n-1)(n-2)/2 = C(n-1,2).

  THIS CONFIRMS: dim(Z_1) = C(n-1,2) for all tournaments. ✓
  """)

    # ==================================================================
    # NOW THE DELETION BOUND
    # ==================================================================
    print("\n" + "=" * 70)
    print("DELETION BOUND: THE ALGEBRA")
    print("=" * 70)

    print("""
  For T on n vertices with beta_1(T) = 0:

  dim(B_1(T)) = dim(Z_1(T)) = C(n-1,2)

  For T\\v on (n-1) vertices:
  dim(Z_1(T\\v)) = C(n-2,2)
  dim(B_1(T\\v)) = C(n-2,2) - beta_1(T\\v)

  |A_1(T)| = C(n,2)
  |A_1(T\\v)| = C(n-1,2)

  Now B_1(T) ⊂ R^{C(n,2)} (subspace of dimension C(n-1,2)).
  B_1(T\\v) ⊂ R^{C(n-1,2)} (subspace of dimension C(n-2,2) - beta_1(T\\v)).

  The edges of T partition: C(n,2) = C(n-1,2) + (n-1).
  (C(n-1,2) edges of T\\v, plus (n-1) edges incident to v.)

  Projection: pi_v: R^{C(n,2)} -> R^{C(n-1,2)} drops the v-edges.

  pi_v(B_1(T)) ⊂ Z_1(T\\v) and B_1(T\\v) ⊂ pi_v(B_1(T)).

  KEY: dim(ker(pi_v|_{B_1(T)})) = dim(B_1(T)) - dim(pi_v(B_1(T)))

  = C(n-1,2) - dim(pi_v(B_1(T)))

  The kernel consists of elements of B_1(T) supported entirely on v-edges.
  But a 1-cycle supported on v-edges must have: for each vertex u ≠ v,
  the flow in equals flow out. Since v-edges are v->u or u->v (for each u),
  a cycle on v-edges would need each u to have equal in/out flow from v.
  But the v-edges form a "star" — a bipartite structure — so the only
  cycle on v-edges is the zero cycle. Therefore ker(pi_v|_{B_1(T)}) = 0!

  Wait, not quite. An element of B_1(T) supported on v-edges:
  z = sum c_e * e where e involves v. For z ∈ Z_1: d_1(z) = 0.
  d_1(z) = sum_u c_{v,u} * (u - v) + sum_u c_{u,v} * (v - u)
  = sum_u (c_{v,u} - c_{u,v}) * (u - v)

  For z ∈ Z_1: coefficient of each vertex must be zero.
  Coefficient of u (for u ≠ v): c_{v,u} - c_{u,v} (from edges involving u)
  Wait — but these are tournament edges. For each u, exactly ONE of (v,u) or (u,v)
  exists as an edge. So c_{v,u} or c_{u,v} is defined but not both.

  Actually, the edge space A_1 for a tournament has basis {(a,b) : a->b in T}.
  An element z = sum c_{a,b} (a,b) where the sum is over edges a->b.
  d_1(z) = sum c_{a,b} ((b) - (a)) = 0
  means: for each vertex w, sum_{a: a->w} c_{a,w} - sum_{b: w->b} c_{w,b} = 0.

  For z supported on v-edges only:
  For w ≠ v: the terms involving w are either c_{v,w} (if v->w) or -c_{w,v} (if w->v).
  So for each w ≠ v: either c_{v,w} = 0 (if v->w) or c_{w,v} = 0 (if w->v).
  This means ALL coefficients are zero. So z = 0.

  THEREFORE: ker(pi_v|_{B_1(T)}) = 0, so pi_v is injective on B_1(T).

  dim(pi_v(B_1(T))) = dim(B_1(T)) = C(n-1,2).

  But pi_v(B_1(T)) ⊂ R^{C(n-1,2)} and dim(Z_1(T\\v)) = C(n-2,2).

  So we need: C(n-1,2) ≤ C(n-1,2). That's trivially true but not useful.

  Wait — pi_v(B_1(T)) has dimension C(n-1,2) as a subspace of R^{C(n-1,2)}?
  That would mean pi_v(B_1(T)) = R^{C(n-1,2)}, i.e., the projection spans everything.
  But B_1(T) ⊂ Z_1(T) (where dim(Z_1(T)) = C(n-1,2)).
  And Z_1(T) is a PROPER subspace of R^{C(n,2)}.

  pi_v(Z_1(T)): the projection of Z_1(T) into R^{C(n-1,2)}.
  Is pi_v injective on Z_1(T)?
  ker(pi_v|_{Z_1(T)}) = Z_1(T) ∩ {elements supported on v-edges} = {0} (same argument).

  So dim(pi_v(Z_1(T))) = C(n-1,2) = dim(R^{C(n-1,2)}).
  This means pi_v(Z_1(T)) = R^{C(n-1,2)}!

  So pi_v(B_1(T)) = pi_v(Z_1(T)) = R^{C(n-1,2)} ⊃ Z_1(T\\v).

  That means B_1(T\\v) ⊂ pi_v(B_1(T)) ∩ Z_1(T\\v) = Z_1(T\\v)
  (which is obvious) and dim(pi_v(B_1(T)) ∩ Z_1(T\\v)) = C(n-2,2).

  Hmm, this doesn't immediately give the constraint. Let me think differently.
  """)

    # Let me verify: is pi_v injective on Z_1(T)?
    for n in [5]:
        print(f"\n  Verifying pi_v injective on Z_1(T) at n={n}:")
        for idx, A in enumerate(all_tournaments(n)):
            info = compute_beta1_full(A, n)
            if info['beta_1'] != 0:
                continue
            edges = info['allowed_1']
            for v in range(n):
                proj_rows = [i for i, e in enumerate(edges) if v not in e]
                # Z_1(T) = ker(d_1|_{Omega_1})
                # Build d1
                bd1 = build_full_boundary_matrix(edges, [(u,) for u in range(n)])
                # Kernel of d1
                U, S, Vt = np.linalg.svd(bd1, full_matrices=True)
                r = sum(s > 1e-8 for s in S)
                Z1_basis = Vt[r:].T  # columns are Z_1 basis

                # Project Z_1 to non-v edges
                Z1_proj = Z1_basis[proj_rows, :]
                rank_proj = rank_of_matrix(Z1_proj)
                if rank_proj != Z1_basis.shape[1]:
                    print(f"    NOT injective: n={n}, v={v}, rank={rank_proj}, dim(Z1)={Z1_basis.shape[1]}")

            # Just check one
            print(f"    Tournament idx={idx}: pi_v injective on Z_1 for all v ✓")
            break

    # ==================================================================
    # THE RIGHT APPROACH: B_1(T) = B_1(T\v) + (v-contributions)
    # ==================================================================
    print("\n" + "=" * 70)
    print("THE RIGHT APPROACH: DECOMPOSING B_1(T)")
    print("=" * 70)

    print("""
  Since pi_v is injective on B_1(T), we can work in the projected space.

  pi_v(B_1(T)) has dimension C(n-1,2) in R^{C(n-1,2)}.
  So pi_v(B_1(T)) = R^{C(n-1,2)}.

  B_1(T\\v) has dimension C(n-2,2) - beta_1(T\\v) inside Z_1(T\\v).
  Z_1(T\\v) has dimension C(n-2,2) inside R^{C(n-1,2)}.

  So B_1(T\\v) is a codimension beta_1(T\\v) subspace of Z_1(T\\v).

  Now sum_v beta_1(T\\v) = sum_v (dim(Z_1(T\\v)) - dim(B_1(T\\v)))
                         = n * C(n-2,2) - sum_v dim(B_1(T\\v))

  We need: sum_v dim(B_1(T\\v)) >= n*C(n-2,2) - 3.

  Equivalently: sum_v codim_{Z_1(T\\v)}(B_1(T\\v)) <= 3.

  The codimension measures how many "missing" boundary generators there are
  when we restrict to T\\v.
  """)

    # Compute sum_v dim(B1(T\v)) vs n*C(n-2,2)
    for n in [5, 6]:
        print(f"\n  n={n}: n*C(n-2,2) = {n * (n-2)*(n-3)//2}")
        if n == 5:
            for A in all_tournaments(n):
                info = compute_beta1_full(A, n)
                if info['beta_1'] != 0:
                    continue
                total_B1 = 0
                for v in range(n):
                    Av, nv = subtournament(A, n, [v])
                    info_v = compute_beta1_full(Av, nv)
                    total_B1 += info_v['dim_B1']
                deficit = n * (n-2)*(n-3)//2 - total_B1
                # Just show a few
                if deficit > 0:
                    print(f"    deficit = {deficit}")
                    break
            else:
                print(f"    All deficits verified")

    # ==================================================================
    # FINAL: The 3-cycle argument
    # ==================================================================
    print("\n" + "=" * 70)
    print("FINAL: THE 3-CYCLE ARGUMENT")
    print("=" * 70)

    print("""
  KEY INSIGHT: At n=4, beta_1(T)=1 iff t3(T)=2 (all 4 triples are 3-cycles
  is impossible; max t3=4 at n=4 but wait, C(4,3)=4 triples, and t3=2 already
  gives beta_1=1).

  Actually at n=4: t3=0 or 1 => beta_1=0, t3=2 => beta_1=1.

  At n=5: beta_1(T\\v)=1 requires t3(T\\v) >= 2 (i.e., the sub-tournament
  on 4 vertices has at least 2 three-cycles).

  Now: T on 5 vertices, beta_1(T)=0. When does T\\v have t3>=2?

  Each 3-cycle in T\\v is also a 3-cycle in T (that avoids v).
  Total 3-cycles in T avoiding v: t3 - #(3-cycles through v).

  For beta_1(T\\v)=1: we need t3(T\\v) >= 2.

  A 3-cycle through v involves v and 2 other vertices.
  There are C(n-1,2) possible pairs for these 2 vertices.

  If all 3-cycles in T go through 3 specific vertices a,b,c,
  then T\\a, T\\b, T\\c could each still have 3-cycles (that avoid a,b,c resp.),
  but there's a limit to how many 3-cycles can be distributed.

  THE CONSTRAINT: For beta_1(T)=0, we need B_1(T) = Z_1(T).
  This requires enough transitive triples to generate all of Z_1.
  If too many triples are 3-cycles, Omega_2 shrinks and B_1 loses rank.

  Specifically: each 3-cycle triple {a,b,c} removes one transitive-triple
  2-path from A_2 and replaces it with 3-cycle 2-paths that are NOT individually
  in Omega_2. The Omega_2 condition creates dependencies.
  """)

    # Connection between t3 and beta_1 at n=4
    print("  n=4: t3 threshold for beta_1=1:")
    for A in all_tournaments(4):
        t3 = count_3cycles(A, 4)
        b1 = compute_beta1_only(A, 4)
        if t3 == 2:
            print(f"    t3=2: beta_1={b1}")
            break

    # At n=5, what are the t3 values for T with beta_1=0 that have 3 bad vertices?
    print("\n  n=5: t3 distribution for beta_1=0 tournaments by # bad vertices:")
    by_bad_count = defaultdict(Counter)
    for A in all_tournaments(5):
        info = compute_beta1_full(A, 5)
        if info['beta_1'] != 0:
            continue
        n_bad = 0
        for v in range(5):
            Av, nv = subtournament(A, 5, [v])
            if compute_beta1_only(Av, nv) > 0:
                n_bad += 1
        t3 = count_3cycles(A, 5)
        by_bad_count[n_bad][t3] += 1

    for nb in sorted(by_bad_count.keys()):
        print(f"    {nb} bad vertices: t3 distribution = {dict(sorted(by_bad_count[nb].items()))}")

if __name__ == '__main__':
    run_proof_analysis()
