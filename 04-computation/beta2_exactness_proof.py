"""
beta2_exactness_proof.py — WHY is the chain complex exact at Omega_2?

We know:
1. Constraint matrix at Omega_2 has full rank (disjoint support theorem)
2. dim(Omega_2) = |A_2| - #NA_faces

We need: ker(d_2|Omega_2) = im(d_3|Omega_3)
i.e., every 2-cycle in Omega_2 is a 2-boundary from Omega_3.

APPROACH: Study the "rank defect" of d_2 and show it equals im(d_3).

Alternative approach: show that the chain complex
  Omega_3 --d_3--> Omega_2 --d_2--> Omega_1
is exact at Omega_2 by an explicit construction or dimension count.

Key formulas (from disjoint support):
  dim(Omega_2) = |A_2| - #{(a,c) : c->a, exists b with a->b->c}
  dim(Omega_1) = C(n,2)
  rank(d_1|Omega_1) = n - 1

We need to understand:
  rank(d_2|Omega_2) and rank(d_3|Omega_3) separately,
  and show their sum = dim(Omega_2).

NEW IDEA: Study the boundary map d_2 BEFORE restricting to Omega.
On all of R^{A_2}: d_2(a,b,c) = (b,c) - (a,c)[if a->c] + (a,b)
Then Omega_2 = {c in R^{A_2} : for each NA (a,c), sum_b coeff(a,b,c) = 0}

Can we decompose R^{A_2} and study the maps?

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = int(sum(s > 1e-10 for s in S))
    ns = Vt[rank:].T
    return ns if ns.shape[1] > 0 else np.zeros((dim_Ap, 0))

def build_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx:
                M[idx[face], j] += sign
    return M


def main():
    print("=" * 70)
    print("EXACTNESS AT OMEGA_2: DEEP STRUCTURE")
    print("=" * 70)

    # =====================================================================
    # APPROACH 1: Decompose A_2 paths by their "skip type"
    #
    # A 2-path (a,b,c) has skip edge (a,c).
    # Type F: a->c (forward skip). No constraint. Path is "free" in Omega_2.
    # Type B: c->a (backward skip). Path is "constrained" in Omega_2.
    #
    # Type B paths group into constraint classes C_{a,c} = {(a,b,c) : a->b->c, c->a}
    # Each class has |C_{a,c}| = #{b : a->b, b->c} paths.
    # The constraint: sum of coefficients in each class = 0.
    # So each class contributes |C_{a,c}| - 1 free dimensions.
    #
    # dim(Omega_2) = #Type_F + sum_{NA (a,c)} (|C_{a,c}| - 1)
    #              = #Type_F + #Type_B - #NA_faces
    #              = |A_2| - #NA_faces  (confirmed)
    #
    # Now: what does d_2 do to each type?
    # d_2(a,b,c) = (b,c) - (a,c) + (a,b)  in A_1 (= Omega_1)
    #
    # Type F (a->c): d_2(a,b,c) = (b,c) - (a,c) + (a,b)
    #   All three edges are present in the tournament.
    # Type B (c->a): d_2(a,b,c) = (b,c) + (a,b)  [no (a,c) term since a->c is FALSE]
    #   Wait - GLMY boundary maps to Omega_{p-1}, not A_{p-1}.
    #   For Omega_2, the boundary d_2: Omega_2 -> Omega_1 = A_1.
    #   But the path (a,c) is NOT in A_1 when c->a.
    #   The GLMY boundary of (a,b,c) is still (b,c) - (a,c) + (a,b) as elements of R^{paths}.
    #   But in R^{A_1}, (a,c) doesn't exist (c->a), so it's projected out.
    #
    #   Actually, the GLMY construction says:
    #   d: R^{A_p} -> R^{A_{p-1}} is the boundary restricted to allowed paths.
    #   The non-allowed face just doesn't appear.
    #   So d_2(a,b,c) = sum_{face in A_1} sign * face
    #   For Type F: = (b,c) - (a,c) + (a,b) [all 3 in A_1]
    #   For Type B: = (b,c) + (a,b) [only 2 in A_1, middle face drops]
    #
    # Wait, actually I need to be more careful about the GLMY definition.
    # Let me re-read the definition.
    # =====================================================================

    # The GLMY boundary is: d_p : R^{all p-paths} -> R^{all (p-1)-paths}
    # d_p(v_0,...,v_p) = sum_{i=0}^p (-1)^i (v_0,...,hat{v_i},...,v_p)
    #
    # Omega_p = {c in R^{A_p} : d(c) has no component on non-allowed (p-1)-paths}
    #
    # So d(c) is computed using ALL faces (including non-allowed ones),
    # and the constraint is that the non-allowed components vanish.
    # Then the ACTUAL boundary d_p|Omega_p maps into Omega_{p-1} automatically.
    #
    # For our case: if c is in Omega_2, then d_2(c) = sum alpha_j * d_2(a_j,b_j,c_j)
    # where d_2(a,b,c) = (b,c) - (a,c) + (a,b) in R^{all 1-paths}.
    # The Omega_2 constraint ensures: for each non-allowed (a,c), the coefficient of (a,c) is 0.
    # So d_2(c) automatically has no non-allowed components.
    # And d_2(c) in R^{A_1} = Omega_1.
    #
    # So: d_2: Omega_2 -> Omega_1 sends
    # c = sum alpha_j e_{(a_j,b_j,c_j)} to
    # sum alpha_j [(b_j,c_j) - (a_j,c_j) + (a_j,b_j)]
    # where (a_j,c_j) only contributes if a_j -> c_j.

    print("\n--- APPROACH 1: Boundary structure by skip type ---")

    for n in [4, 5, 6]:
        N = 2**(n*(n-1)//2)
        rank_d2_equals_typeF_plus = 0  # does rank(d_2) relate to type counts?
        total = 0
        data = []

        for trial in range(N):
            A = bits_to_adj(trial, n)
            allowed = {}
            for p in range(-1, min(n+1, 6)):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            # Count types
            a1_set = set(allowed[1])
            type_F = 0  # forward skip: a->c
            type_B = 0  # backward skip: c->a
            na_faces = set()
            class_sizes = []

            # Group by NA face
            face_groups = defaultdict(list)
            for j, path in enumerate(allowed[2]):
                a, b, c = path
                if A[a][c]:  # a->c: forward
                    type_F += 1
                else:  # c->a: backward
                    type_B += 1
                    face_groups[(a,c)].append(j)
                    na_faces.add((a,c))

            for face, paths in face_groups.items():
                class_sizes.append(len(paths))

            # Compute Omega_2 and d_2
            ob2 = compute_omega_basis(A, n, 2, allowed[2], allowed[1])
            dim_om2 = ob2.shape[1] if ob2.ndim == 2 else 0

            bd2 = build_boundary_matrix(allowed[2], allowed[1])
            if dim_om2 > 0:
                d2_om = bd2 @ ob2
                sv2 = np.linalg.svd(d2_om, compute_uv=False)
                rank_d2 = int(sum(s > 1e-8 for s in sv2))
            else:
                rank_d2 = 0

            ker_d2 = dim_om2 - rank_d2

            # Also: rank of d_2 on ALL of R^{A_2} (not just Omega_2)
            full_rank_d2 = np.linalg.matrix_rank(bd2, tol=1e-8)

            # And rank of d_3
            ob3 = compute_omega_basis(A, n, 3, allowed[3], allowed[2])
            dim_om3 = ob3.shape[1] if ob3.ndim == 2 else 0
            bd3 = build_boundary_matrix(allowed[3], allowed[2])
            if dim_om3 > 0:
                d3_om = bd3 @ ob3
                sv3 = np.linalg.svd(d3_om, compute_uv=False)
                rank_d3 = int(sum(s > 1e-8 for s in sv3))
            else:
                rank_d3 = 0

            total += 1
            data.append({
                'type_F': type_F, 'type_B': type_B,
                'na_faces': len(na_faces), 'class_sizes': class_sizes,
                'dim_om2': dim_om2, 'rank_d2': rank_d2, 'ker_d2': ker_d2,
                'full_rank_d2': full_rank_d2, 'rank_d3': rank_d3,
                'dim_om3': dim_om3
            })

        print(f"\n  n={n} ({total} tournaments):")

        # Check: is rank_d2 = type_F always? Or some formula involving types?
        # rank(d_2|Omega_2) = dim(Omega_2) - ker(d_2) = (|A_2| - #NA) - rank(d_3)
        # So rank(d_2) = |A_2| - #NA - rank(d_3) IF beta_2 = 0.
        # = (type_F + type_B) - #NA - rank(d_3)
        # = type_F + (type_B - #NA) - rank(d_3)
        # = type_F + sum(|C_{a,c}| - 1) - rank(d_3)

        # What is rank(d_2) vs type_F?
        rank_vs_typeF = Counter()
        for d in data:
            diff = d['rank_d2'] - d['type_F']
            rank_vs_typeF[diff] += 1
        print(f"    rank(d_2|Omega_2) - type_F: {dict(sorted(rank_vs_typeF.items()))}")

        # What is rank(d_2) on full A_2 vs restricted to Omega_2?
        rank_diff = Counter()
        for d in data:
            rank_diff[d['full_rank_d2'] - d['rank_d2']] += 1
        print(f"    rank(d_2|A_2) - rank(d_2|Omega_2): {dict(sorted(rank_diff.items()))}")

        # What is full_rank_d2?
        full_rank_dist = Counter()
        for d in data:
            full_rank_dist[d['full_rank_d2']] += 1
        print(f"    rank(d_2|A_2) distribution: {dict(sorted(full_rank_dist.items()))}")
        print(f"    C(n,2) = {comb(n,2)}, C(n,2)-1 = {comb(n,2)-1}")

        # Key check: is full_rank_d2 always = C(n,2) - 1?
        # This would mean d_2 on all of A_2 always has rank n(n-1)/2 - 1.
        # Since d_1 * d_2 = 0, im(d_2) subset ker(d_1).
        # dim(ker d_1) = C(n,2) - (n-1).
        # So rank(d_2|A_2) <= C(n,2) - n + 1.
        max_possible = comb(n, 2) - n + 1
        print(f"    Max possible rank(d_2) = C(n,2) - n + 1 = {max_possible}")

        # CRITICAL: does d_2 on A_2 always achieve max rank?
        achieves_max = sum(1 for d in data if d['full_rank_d2'] == max_possible)
        print(f"    Achieves max rank: {achieves_max}/{total}")

    # =====================================================================
    # KEY FINDING: Let's check if d_2 on A_2 always has rank C(n,2) - n + 1.
    # If so: im(d_2|A_2) = ker(d_1) (since d_1*d_2 = 0 and rank = codim).
    # This would mean: the chain complex R^{A_p} is exact at A_1!
    # i.e., H_1(A) = 0 always for tournaments.
    #
    # Actually, this is exactly beta_1 for the UNRESTRICTED complex.
    # The Omega complex restricts to paths in the tournament, but
    # if the A complex is exact at level 1, that tells us something.
    # =====================================================================

    print("\n--- KEY: Is the A_p complex exact at A_1? ---")
    for n in [4, 5, 6]:
        N = 2**(n*(n-1)//2)
        exact_at_A1 = 0
        for trial in range(N):
            A = bits_to_adj(trial, n)
            a1 = enumerate_allowed_paths(A, n, 1)
            a2 = enumerate_allowed_paths(A, n, 2)
            bd2 = build_boundary_matrix(a2, a1)
            bd1 = build_boundary_matrix(a1, [(v,) for v in range(n)])
            # Check: im(d_2) = ker(d_1)?
            # rank(d_2) should equal dim(ker d_1) = C(n,2) - (n-1)
            rank_d2_full = np.linalg.matrix_rank(bd2, tol=1e-8)
            ker_d1_dim = comb(n, 2) - (n - 1)
            if rank_d2_full == ker_d1_dim:
                exact_at_A1 += 1
        print(f"  n={n}: A-complex exact at A_1 in {exact_at_A1}/{N} cases")
        # This is beta_1(A) = 0 for the A-level complex

    # =====================================================================
    # APPROACH 2: What formulas describe rank(d_3|Omega_3)?
    #
    # If beta_2 = 0, then ker(d_2) = im(d_3), so rank(d_3) = ker(d_2).
    # And rank(d_2) = dim(Omega_2) - ker(d_2).
    # So rank(d_2) + rank(d_3) = dim(Omega_2).
    #
    # Question: is rank(d_3|Omega_3) related to a simple combinatorial quantity?
    # =====================================================================

    print("\n--- APPROACH 2: rank(d_3|Omega_3) formula search ---")
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        formula_data = []
        for trial in range(N):
            A = bits_to_adj(trial, n)
            allowed = {}
            for p in range(-1, min(n+1, 6)):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            ob3 = compute_omega_basis(A, n, 3, allowed[3], allowed[2])
            dim_om3 = ob3.shape[1] if ob3.ndim == 2 else 0
            bd3 = build_boundary_matrix(allowed[3], allowed[2])
            if dim_om3 > 0:
                d3_om = bd3 @ ob3
                sv3 = np.linalg.svd(d3_om, compute_uv=False)
                rank_d3 = int(sum(s > 1e-8 for s in sv3))
            else:
                rank_d3 = 0

            c3 = sum(1 for i,j,k in combinations(range(n), 3)
                     if (A[i][j] and A[j][k] and A[k][i]) or
                        (A[i][k] and A[k][j] and A[j][i]))

            formula_data.append({
                'c3': c3, 'rank_d3': rank_d3, 'dim_om3': dim_om3,
                '|A_2|': len(allowed[2]), '|A_3|': len(allowed[3])
            })

        # Try: is rank_d3 a function of c3?
        c3_to_rd3 = defaultdict(set)
        for d in formula_data:
            c3_to_rd3[d['c3']].add(d['rank_d3'])
        print(f"\n  n={n}: rank(d_3) as function of c3:")
        for c3, vals in sorted(c3_to_rd3.items()):
            print(f"    c3={c3}: rank_d3 in {sorted(vals)}")

    # =====================================================================
    # APPROACH 3: Study rank(d_2|Omega_2) vs rank(d_2|A_2)
    #
    # The restriction from A_2 to Omega_2 cuts out some paths.
    # Does this REDUCE the rank of d_2, or keep it the same?
    # If rank(d_2|Omega_2) = rank(d_2|A_2) - something predictable,
    # we might get a formula.
    # =====================================================================

    print("\n--- APPROACH 3: rank(d_2|Omega_2) vs rank(d_2|A_2) ---")
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        diff_dist = Counter()
        diff_to_na = defaultdict(Counter)
        for trial in range(N):
            A = bits_to_adj(trial, n)
            allowed = {}
            for p in range(-1, min(n+1, 6)):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            a1_set = set(allowed[1])
            na_faces = set()
            for path in allowed[2]:
                for sign, face in boundary_coeffs(path):
                    if len(face) == 2 and tuple(face) not in a1_set:
                        na_faces.add(tuple(face))

            ob2 = compute_omega_basis(A, n, 2, allowed[2], allowed[1])
            dim_om2 = ob2.shape[1] if ob2.ndim == 2 else 0

            bd2 = build_boundary_matrix(allowed[2], allowed[1])
            full_rank = np.linalg.matrix_rank(bd2, tol=1e-8)
            if dim_om2 > 0:
                d2_om = bd2 @ ob2
                om_rank = np.linalg.matrix_rank(d2_om, tol=1e-8)
            else:
                om_rank = 0

            diff = full_rank - om_rank
            diff_dist[diff] += 1
            diff_to_na[diff][len(na_faces)] += 1

        print(f"\n  n={n}: rank(d_2|A_2) - rank(d_2|Omega_2) distribution:")
        for d, cnt in sorted(diff_dist.items()):
            na_vals = sorted(diff_to_na[d].items())
            print(f"    diff={d}: {cnt} cases, #NA distribution={na_vals}")

    # =====================================================================
    # APPROACH 4: Is d_2 SURJECTIVE on ker(d_1)?
    #
    # Recall: beta_1 = dim(ker d_1 / im d_2|Omega_2)
    # If im(d_2|Omega_2) = ker(d_1) then beta_1 = 0.
    # If im(d_2|Omega_2) has codim 1 in ker(d_1) then beta_1 = 1.
    #
    # We know beta_1 in {0,1} for all tournaments.
    # So rank(d_2|Omega_2) = ker_d1_dim or ker_d1_dim - 1.
    #
    # From n=6 data: rank(d_2|A_2) distribution was [10,11].
    # And ker_d1_dim = C(6,2) - 5 = 10.
    # So rank(d_2|A_2) can be 10 or 11.
    # But rank must be <= ker_d1_dim = 10 (since im(d_2) subset ker(d_1)).
    # Contradiction! Unless I'm computing wrong. Let me check.
    # =====================================================================

    print("\n--- APPROACH 4: rank(d_2|A_2) vs ker(d_1) check ---")
    for n in [4, 5, 6]:
        N = 2**(n*(n-1)//2)
        ker_d1 = comb(n, 2) - (n - 1)
        print(f"\n  n={n}: ker(d_1) dim = {ker_d1}")
        surj_count = 0
        rank_dist = Counter()
        for trial in range(N):
            A = bits_to_adj(trial, n)
            a1 = enumerate_allowed_paths(A, n, 1)
            a2 = enumerate_allowed_paths(A, n, 2)
            bd2 = build_boundary_matrix(a2, a1)
            rank_d2_full = np.linalg.matrix_rank(bd2, tol=1e-8)
            rank_dist[rank_d2_full] += 1
            if rank_d2_full == ker_d1:
                surj_count += 1
        print(f"    rank(d_2|A_2) distribution: {dict(sorted(rank_dist.items()))}")
        print(f"    Surjective (rank = ker_d1 = {ker_d1}): {surj_count}/{N}")

    # =====================================================================
    # CRUCIAL TEST: Does d_2 on A_2 always map onto ker(d_1)?
    # i.e., is H_1 of the A-complex always 0?
    #
    # If YES for tournaments, this would be a strong structural property.
    # Combined with the Omega constraint analysis, it might prove beta_2 = 0.
    # =====================================================================

    print("\n--- IS H_1(A-complex) = 0 FOR ALL TOURNAMENTS? ---")
    for n in [4, 5, 6]:
        N = 2**(n*(n-1)//2)
        ker_d1_dim = comb(n, 2) - (n - 1)
        h1_zero = 0
        h1_nonzero = 0
        for trial in range(N):
            A = bits_to_adj(trial, n)
            a1 = enumerate_allowed_paths(A, n, 1)
            a2 = enumerate_allowed_paths(A, n, 2)
            a0 = [(v,) for v in range(n)]
            bd1 = build_boundary_matrix(a1, a0)
            bd2 = build_boundary_matrix(a2, a1)

            # im(d_2) subset ker(d_1)
            # Check if they're equal
            rank_d1 = np.linalg.matrix_rank(bd1, tol=1e-8)
            rank_d2 = np.linalg.matrix_rank(bd2, tol=1e-8)
            ker_d1 = len(a1) - rank_d1

            if rank_d2 == ker_d1:
                h1_zero += 1
            else:
                h1_nonzero += 1

        print(f"  n={n}: H_1(A) = 0 in {h1_zero}/{N} cases"
              f" (H_1 != 0 in {h1_nonzero})")

    # =====================================================================
    # So H_1(A) is NOT always 0. Good — this means Omega matters.
    #
    # Let me now check: what's the relationship between H_1(A) and beta_1(Omega)?
    # And what happens at levels 0,1,2 of the A-complex vs Omega-complex?
    # =====================================================================

    print("\n--- H_p(A) vs beta_p(Omega) for p=0,1,2 ---")
    for n in [5, 6]:
        N = 2**(n*(n-1)//2)
        comparison = Counter()
        for trial in range(N):
            A = bits_to_adj(trial, n)
            allowed = {}
            for p in range(-1, min(n+1, 6)):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            a0 = [(v,) for v in range(n)]
            bd1 = build_boundary_matrix(allowed[1], a0)
            bd2 = build_boundary_matrix(allowed[2], allowed[1])
            bd3 = build_boundary_matrix(allowed[3], allowed[2])

            # H_1(A): ker(d_1)/im(d_2)
            rank_bd1 = np.linalg.matrix_rank(bd1, tol=1e-8)
            rank_bd2 = np.linalg.matrix_rank(bd2, tol=1e-8)
            ker_d1 = len(allowed[1]) - rank_bd1
            h1_A = ker_d1 - rank_bd2

            # H_2(A): ker(d_2)/im(d_3)
            rank_bd3 = np.linalg.matrix_rank(bd3, tol=1e-8)
            ker_d2_A = len(allowed[2]) - rank_bd2
            h2_A = ker_d2_A - rank_bd3

            # beta_1(Omega)
            ob2 = compute_omega_basis(A, n, 2, allowed[2], allowed[1])
            dim_om2 = ob2.shape[1] if ob2.ndim == 2 else 0
            if dim_om2 > 0:
                d2_om = bd2 @ ob2
                rank_d2_om = np.linalg.matrix_rank(d2_om, tol=1e-8)
            else:
                rank_d2_om = 0
            beta_1 = ker_d1 - rank_d2_om

            comparison[(h1_A, h2_A, beta_1)] += 1

        print(f"\n  n={n}: (H_1(A), H_2(A), beta_1(Omega)) distribution:")
        for key, cnt in sorted(comparison.items(), key=lambda x: -x[1]):
            h1, h2, b1 = key
            print(f"    H_1(A)={h1}, H_2(A)={h2}, beta_1={b1}: {cnt} cases")

    print("\nDONE.")


if __name__ == '__main__':
    main()
