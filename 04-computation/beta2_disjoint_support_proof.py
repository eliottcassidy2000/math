"""
beta2_disjoint_support_proof.py — Prove beta_2 = 0 via disjoint support structure.

Key insight: At level 2 in tournaments, each 2-path (a,b,c) has AT MOST ONE
non-allowed face: (a,c) when c->a. This means the Omega_2 constraint matrix
has disjoint row supports, hence full row rank.

But full rank of constraint matrix is NOT enough for beta_2=0.
We need: ker(d_2|Omega_2) = im(d_3|Omega_3).

This script:
1. Verifies the disjoint support property (proved below)
2. Studies ker(d_2) explicitly to understand why it equals im(d_3)
3. Compares with level 4 where the property fails
4. Looks for a structural reason for exactness at level 2

Author: kind-pasteur-S45 (2026-03-09)
"""
import sys
import numpy as np
from math import comb
from itertools import combinations
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
    print("DISJOINT SUPPORT PROOF AND EXACTNESS ANALYSIS")
    print("=" * 70)

    # =====================================================================
    # THEOREM (Disjoint Support):
    # In a tournament, every 2-path (a,b,c) has at most one non-allowed
    # 1-face. This is the face (a,c) at position 1, which is non-allowed
    # iff c->a in the tournament.
    #
    # PROOF:
    # The faces of (a,b,c) are:
    #   d_0(a,b,c) = (b,c)  — in A_1 iff b->c, which is TRUE (a->b->c is allowed)
    #   d_1(a,b,c) = (a,c)  — in A_1 iff a->c
    #   d_2(a,b,c) = (a,b)  — in A_1 iff a->b, which is TRUE (a->b->c is allowed)
    # So only d_1 = (a,c) can be non-allowed, and this happens iff c->a.
    # QED
    # =====================================================================

    # Verify exhaustively
    print("\n--- THEOREM: Disjoint support at level 2 ---")
    for n in [4, 5, 6]:
        N = 2**(n*(n-1)//2)
        violations = 0
        for trial in range(N):
            A = bits_to_adj(trial, n)
            a1 = enumerate_allowed_paths(A, n, 1)
            a2 = enumerate_allowed_paths(A, n, 2)
            a1_set = set(a1)
            for path in a2:
                na_count = sum(1 for sign, face in boundary_coeffs(path)
                              if len(face) == 2 and tuple(face) not in a1_set)
                if na_count > 1:
                    violations += 1
        print(f"  n={n}: {violations} violations in {N} tournaments (should be 0)")

    # =====================================================================
    # COROLLARY: Constraint matrix at Omega_2 has full row rank.
    # Since each path appears in at most one constraint row, and each
    # constraint has at least one path in its support (by definition of
    # NA face), the rows have disjoint non-empty supports, hence are
    # linearly independent.
    # =====================================================================

    print("\n--- COROLLARY: Full rank constraints at Omega_2 ---")
    for n in [4, 5, 6]:
        N = 2**(n*(n-1)//2)
        full_rank = 0
        for trial in range(N):
            A = bits_to_adj(trial, n)
            a1 = enumerate_allowed_paths(A, n, 1)
            a2 = enumerate_allowed_paths(A, n, 2)
            a1_set = set(a1)
            na_faces = set()
            for path in a2:
                for sign, face in boundary_coeffs(path):
                    if len(face) == 2 and tuple(face) not in a1_set:
                        na_faces.add(tuple(face))
            if len(na_faces) == 0:
                full_rank += 1
                continue
            # Build constraint matrix
            na_list = list(na_faces)
            na_idx = {f: i for i, f in enumerate(na_list)}
            P = np.zeros((len(na_list), len(a2)))
            for j, path in enumerate(a2):
                for sign, face in boundary_coeffs(path):
                    ft = tuple(face) if len(face) == 2 else None
                    if ft in na_idx:
                        P[na_idx[ft], j] += sign
            sv = np.linalg.svd(P, compute_uv=False)
            rank = int(sum(s > 1e-10 for s in sv))
            if rank == len(na_list):
                full_rank += 1
        print(f"  n={n}: full rank in {full_rank}/{N} tournaments")

    # =====================================================================
    # KEY QUESTION: Why does ker(d_2) = im(d_3) at Omega level?
    #
    # Let's track the exact dimensions:
    # dim(Omega_2), ker(d_2), im(d_3), and their relationship
    # =====================================================================

    print("\n--- EXACTNESS AT LEVEL 2: Detailed dimension tracking ---")

    for n in [4, 5, 6]:
        N = 2**(n*(n-1)//2)
        print(f"\n  n={n} (exhaustive, {N} tournaments):")

        # Track (dim_om2, ker_d2, dim_om3, im_d3, rank_d2) distributions
        from collections import Counter
        dim_dist = Counter()
        for trial in range(N):
            A = bits_to_adj(trial, n)
            allowed = {}
            for p in range(-1, min(n+1, 6)):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            ob2 = compute_omega_basis(A, n, 2, allowed[2], allowed[1])
            ob3 = compute_omega_basis(A, n, 3, allowed[3], allowed[2])
            dim_om2 = ob2.shape[1] if ob2.ndim == 2 else 0
            dim_om3 = ob3.shape[1] if ob3.ndim == 2 else 0

            # d_2 on Omega_2
            bd2 = build_boundary_matrix(allowed[2], allowed[1])
            if dim_om2 > 0:
                d2_om = bd2 @ ob2
                sv2 = np.linalg.svd(d2_om, compute_uv=False)
                rank_d2 = int(sum(s > 1e-8 for s in sv2))
            else:
                rank_d2 = 0
            ker_d2 = dim_om2 - rank_d2

            # d_3 on Omega_3
            bd3 = build_boundary_matrix(allowed[3], allowed[2])
            if dim_om3 > 0:
                d3_om = bd3 @ ob3
                sv3 = np.linalg.svd(d3_om, compute_uv=False)
                im_d3 = int(sum(s > 1e-8 for s in sv3))
            else:
                im_d3 = 0

            beta2 = ker_d2 - im_d3
            dim_dist[(dim_om2, ker_d2, dim_om3, im_d3, beta2)] += 1

        for key, cnt in sorted(dim_dist.items(), key=lambda x: -x[1]):
            do2, kd2, do3, id3, b2 = key
            print(f"    dim_om2={do2}, ker_d2={kd2}, dim_om3={do3}, im_d3={id3}, beta2={b2}: {cnt} cases")

    # =====================================================================
    # HYPOTHESIS: ker(d_2) = im(d_3) because of a STRUCTURAL identity
    # related to tournaments being complete oriented graphs.
    #
    # In a tournament, every pair of vertices has an edge.
    # This means: for any 2-cycle in Omega_2, there's always a
    # "filling" from Omega_3 that bounds it.
    #
    # Let's check: when ker(d_2) > 0, what do the kernel elements look like?
    # And can we always explicitly construct a 3-chain that bounds them?
    # =====================================================================

    print("\n--- KERNEL ELEMENT ANALYSIS ---")
    n = 5
    N = 2**(n*(n-1)//2)
    ker_examples = 0
    for trial in range(N):
        A = bits_to_adj(trial, n)
        allowed = {}
        for p in range(-1, 5):
            allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

        ob2 = compute_omega_basis(A, n, 2, allowed[2], allowed[1])
        dim_om2 = ob2.shape[1] if ob2.ndim == 2 else 0
        if dim_om2 == 0: continue

        bd2 = build_boundary_matrix(allowed[2], allowed[1])
        d2_om = bd2 @ ob2
        U, S, Vt = np.linalg.svd(d2_om, full_matrices=True)
        rank_d2 = int(sum(s > 1e-8 for s in S))
        ker_d2 = dim_om2 - rank_d2
        if ker_d2 == 0: continue

        ker_examples += 1
        if ker_examples > 5: break

        # Get kernel elements in A_2 coordinates
        ker_om2 = Vt[rank_d2:].T  # Omega_2 coords
        ker_A2 = ob2 @ ker_om2      # A_2 coords

        c3 = sum(1 for i,j,k in combinations(range(n), 3)
                 if (A[i][j] and A[j][k] and A[k][i]) or
                    (A[i][k] and A[k][j] and A[j][i]))

        print(f"\n  trial {trial}: c3={c3}, dim_om2={dim_om2}, ker_d2={ker_d2}")

        for k in range(ker_d2):
            vec = ker_A2[:, k]
            support = [(allowed[2][j], round(vec[j], 4))
                       for j in range(len(vec)) if abs(vec[j]) > 1e-8]
            print(f"    ker element {k}: {len(support)} terms")
            for path, coeff in support:
                a, b, c = path
                skip_allowed = "a->c" if A[a][c] else "c->a"
                print(f"      {coeff:+.4f} * ({a},{b},{c})  [{skip_allowed}]")

    # =====================================================================
    # NEW APPROACH: Acyclicity of the filling complex
    #
    # For a tournament, the key observation might be that the 2-cycles
    # in Omega_2 can always be "filled" because the tournament graph
    # is rich enough in 3-paths.
    #
    # At level 2: ker(d_2) consists of formal sums of 2-paths whose
    # boundary in Omega_1 is zero. Can we show every such cycle is
    # automatically a boundary?
    #
    # Approach: count dimensions.
    # dim(ker d_2) = dim(Omega_2) - rank(d_2)
    # im(d_3) = rank(d_3|Omega_3)
    #
    # We need: dim(Omega_2) - rank(d_2) = rank(d_3)
    # i.e.: dim(Omega_2) = rank(d_2) + rank(d_3)
    #
    # Using beta_1 formula: rank(d_2) = C(n,2) - n + 1 - beta_1
    # So we need: dim(Omega_2) = C(n,2) - n + 1 - beta_1 + rank(d_3)
    #
    # Since dim(Omega_2) = |A_2| - #NA_faces (full rank constraints):
    # |A_2| - #NA_faces = C(n,2) - n + 1 - beta_1 + rank(d_3)
    #
    # This is the IDENTITY we need to verify/prove.
    # =====================================================================

    print("\n--- IDENTITY: |A_2| - #NA_faces = C(n,2) - n + 1 - beta_1 + rank(d_3) ---")

    for n in [4, 5, 6]:
        N = 2**(n*(n-1)//2)
        identity_holds = 0
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

            lhs = len(allowed[2]) - len(na_faces)

            # Compute beta_1 and rank(d_3)
            ob2 = compute_omega_basis(A, n, 2, allowed[2], allowed[1])
            ob3 = compute_omega_basis(A, n, 3, allowed[3], allowed[2])
            dim_om2 = ob2.shape[1] if ob2.ndim == 2 else 0
            dim_om3 = ob3.shape[1] if ob3.ndim == 2 else 0

            bd2 = build_boundary_matrix(allowed[2], allowed[1])
            if dim_om2 > 0:
                d2_om = bd2 @ ob2
                sv2 = np.linalg.svd(d2_om, compute_uv=False)
                rank_d2 = int(sum(s > 1e-8 for s in sv2))
            else:
                rank_d2 = 0
            beta_1 = comb(n, 2) - n + 1 - rank_d2

            bd3 = build_boundary_matrix(allowed[3], allowed[2])
            if dim_om3 > 0:
                d3_om = bd3 @ ob3
                sv3 = np.linalg.svd(d3_om, compute_uv=False)
                rank_d3 = int(sum(s > 1e-8 for s in sv3))
            else:
                rank_d3 = 0

            rhs = comb(n, 2) - n + 1 - beta_1 + rank_d3

            if lhs == rhs:
                identity_holds += 1
            else:
                print(f"  FAIL at n={n}, trial {trial}: LHS={lhs}, RHS={rhs}")
                print(f"    |A_2|={len(allowed[2])}, #NA={len(na_faces)}, beta_1={beta_1}, rank_d3={rank_d3}")

        print(f"  n={n}: identity holds in {identity_holds}/{N} tournaments"
              f" {'(= beta_2=0!)' if identity_holds == N else ''}")

    # =====================================================================
    # ALTERNATIVE: Direct combinatorial approach
    #
    # #NA_faces at level 2 = number of backward pairs (a,c) with c->a
    # that participate in some 3-cycle.
    #
    # For each such pair, the number of intermediate vertices b
    # (with a->b->c) is the number of common successors of a and
    # predecessors of c. Call this d_{a,c}.
    #
    # Then: #NA_faces = number of backward pairs with d_{a,c} >= 1
    # And the constraint for that pair is: sum of d_{a,c} path coefficients = 0.
    # The number of constrained paths = sum_{(a,c) NA} d_{a,c}.
    # The number of free paths = |A_2| - sum d_{a,c}.
    #
    # But dim(Omega_2) = |A_2| - #NA_faces (each constraint kills one dim).
    # So the free paths plus the constrained-but-surviving paths = dim(Omega_2).
    # Free paths: those (a,b,c) where a->c (no constraint at all).
    # Constrained groups: for each NA face (a,c), d_{a,c} paths are
    # constrained to sum to 0, giving d_{a,c}-1 free parameters.
    #
    # Total dim(Omega_2) = #{(a,b,c) : a->c} + sum_{NA (a,c)} (d_{a,c} - 1)
    #                    = #{(a,b,c) : a->c} + sum d_{a,c} - #NA
    #                    = (|A_2| - sum d_{a,c}) + sum d_{a,c} - #NA
    #                    = |A_2| - #NA
    # Which we already know. OK, this is circular.
    # =====================================================================

    # =====================================================================
    # Let's try: show rank(d_2) + rank(d_3) = dim(Omega_2) directly.
    #
    # rank(d_2|Omega_2) = dim(im d_2) in Omega_1
    # rank(d_3|Omega_3) = dim(im d_3) in Omega_2
    #
    # We need: im(d_2) + ker(d_2) spans Omega_2
    #     AND: ker(d_2) = im(d_3)
    #
    # The first is just dim count: rank + nullity = dim.
    # The second is the content of beta_2 = 0.
    #
    # NEW IDEA: Can we show that im(d_3) generates ker(d_2)?
    # For each kernel element, we need a preimage in Omega_3.
    # The existence of this preimage depends on the tournament structure.
    #
    # For 3-paths: a 3-path (a,b,c,d) is allowed iff a->b, b->c, c->d.
    # Its boundary in A_2 is: (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
    # After projection to Omega_2 (only allowed 2-faces survive):
    # The terms that survive depend on which edges exist.
    # =====================================================================

    print("\n--- BOUNDARY OF 3-PATHS: Which 2-faces survive? ---")
    n = 5
    A = bits_to_adj(42, n)  # pick a specific tournament
    allowed = {}
    for p in range(-1, 5):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

    a2_set = set(allowed[2])
    print(f"  n={n}, trial 42: |A_3|={len(allowed[3])}")
    for path in allowed[3][:10]:
        a, b, c, d = path
        faces_info = []
        for i in range(4):
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            in_a2 = tuple(face) in a2_set
            faces_info.append(f"{'+' if sign > 0 else '-'}({','.join(map(str,face))}){'*' if in_a2 else 'X'}")
        print(f"  ({a},{b},{c},{d}): {' '.join(faces_info)}")
    print("  (* = in A_2, X = not in A_2)")

    # =====================================================================
    # COUNTING NON-ALLOWED FACES OF 3-PATHS IN TOURNAMENTS
    #
    # For 3-path (a,b,c,d): a->b, b->c, c->d
    # Faces:
    #   (b,c,d): allowed iff b->c (YES) and c->d (YES) and b->d?
    #            Actually (b,c,d) is allowed iff b->c AND c->d (consecutive edges).
    #            So (b,c,d) is ALWAYS allowed.
    #   (a,c,d): allowed iff a->c AND c->d. c->d is YES. a->c?
    #            If a->c: allowed. If c->a: NOT allowed.
    #   (a,b,d): allowed iff a->b (YES) AND b->d.
    #            If b->d: allowed. If d->b: NOT allowed.
    #   (a,b,c): allowed iff a->b (YES) AND b->c (YES). ALWAYS allowed.
    #
    # So a 3-path (a,b,c,d) has:
    #   - (a,c,d) non-allowed iff c->a  (1 non-allowed face)
    #   - (a,b,d) non-allowed iff d->b  (1 non-allowed face)
    # Both can be simultaneously non-allowed!
    # So a 3-path can have 0, 1, or 2 non-allowed faces.
    #
    # At level 2: max 1 non-allowed face
    # At level 3: max 2 non-allowed faces
    # At level p: up to p-1 non-allowed faces
    # =====================================================================

    print("\n--- NON-ALLOWED FACE COUNT PER 3-PATH ---")
    for n in [5, 6, 7]:
        if n <= 6:
            N = 2**(n*(n-1)//2)
            rng = None
        else:
            N = 200
            rng = np.random.RandomState(42)

        from collections import Counter
        dist = Counter()
        for trial in range(N):
            if n <= 6:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)
            a2 = enumerate_allowed_paths(A, n, 2)
            a3 = enumerate_allowed_paths(A, n, 3)
            a2_set = set(a2)
            for path in a3:
                na = sum(1 for sign, face in boundary_coeffs(path)
                         if len(set(face)) == len(face) and tuple(face) not in a2_set)
                dist[na] += 1

        total = sum(dist.values())
        print(f"  n={n}: NA faces per 3-path distribution:")
        for k in sorted(dist.keys()):
            print(f"    {k}: {dist[k]} ({100*dist[k]/total:.1f}%)")

    # =====================================================================
    # The critical difference between levels 2 and 3+:
    # Level 2: max 1 NA face per path => disjoint row supports => full rank
    # Level 3: up to 2 NA faces per path => overlapping rows => rank deficit possible
    #
    # But beta_3 CAN be nonzero (at n=6+), while beta_2 CANNOT.
    # So the disjoint support property at level 2 is NECESSARY
    # for beta_2=0, but the overlapping at level 3 is what ALLOWS beta_3>0.
    #
    # However, disjoint support alone doesn't prove beta_2=0.
    # We need the chain complex to be exact.
    #
    # Let me look at this from the chain complex perspective more carefully.
    # =====================================================================

    # =====================================================================
    # APPROACH: Euler characteristic argument
    #
    # chi = sum_{p>=0} (-1)^p dim(Omega_p) = sum_{p>=0} (-1)^p beta_p
    #
    # For a tournament on n vertices:
    # dim(Omega_0) = n
    # dim(Omega_1) = C(n,2)  (all edges, since all edges are allowed)
    #
    # If beta_2 = 0:
    # chi = 1 - beta_1 + 0 - beta_3 + beta_4 - ...
    # chi = 1 - beta_1 - beta_3 + beta_4 - ...
    #
    # But chi = sum (-1)^p dim(Omega_p)
    # = n - C(n,2) + dim(Omega_2) - dim(Omega_3) + ...
    #
    # Can we compute chi directly?
    # =====================================================================

    print("\n--- EULER CHARACTERISTIC ---")
    for n in [4, 5, 6, 7]:
        if n <= 6:
            N = 2**(n*(n-1)//2)
            rng = None
        else:
            N = 200
            rng = np.random.RandomState(42)

        from collections import Counter
        chi_dist = Counter()
        for trial in range(N):
            if n <= 6:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            allowed = {}
            for p in range(-1, n+1):
                allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

            # Compute all Omega dimensions
            od = {}
            for p in range(n):
                ap = allowed[p]
                apm1 = allowed[p-1]
                dim_Ap = len(ap)
                if dim_Ap == 0:
                    od[p] = 0
                    continue
                if p == 0:
                    od[p] = dim_Ap
                    continue
                apm1_set = set(apm1)
                non_allowed = {}
                na_count = 0
                for j, path in enumerate(ap):
                    for sign, face in boundary_coeffs(path):
                        if len(set(face)) == len(face) and face not in apm1_set:
                            if face not in non_allowed:
                                non_allowed[face] = na_count
                                na_count += 1
                if na_count == 0:
                    od[p] = dim_Ap
                    continue
                P_mat = np.zeros((na_count, dim_Ap))
                for j, path in enumerate(ap):
                    for sign, face in boundary_coeffs(path):
                        if face in non_allowed:
                            P_mat[non_allowed[face], j] += sign
                sv = np.linalg.svd(P_mat, compute_uv=False)
                rank = int(sum(s > 1e-10 for s in sv))
                od[p] = dim_Ap - rank

            chi = sum((-1)**p * od.get(p, 0) for p in range(n))
            chi_dist[chi] += 1

            if trial < 3 and n <= 5:
                od_list = [od.get(p, 0) for p in range(n)]
                print(f"  n={n}, trial {trial}: Omega={od_list}, chi={chi}")

        print(f"  n={n}: chi distribution: {dict(sorted(chi_dist.items()))}")

    # =====================================================================
    # Is chi always 1 for tournaments? Let me check.
    # If chi = 1 for all tournaments, then:
    # 1 = 1 - beta_1 + beta_2 - beta_3 + ...
    # => beta_1 - beta_2 + beta_3 - beta_4 + ... = 0
    #
    # But at n=8, beta = [1,0,0,0,1,0,0,0] gives chi = 2.
    # So chi is NOT always 1. It depends on the tournament.
    #
    # What IS chi_omega = sum (-1)^p dim(Omega_p) at the A_p level?
    # chi_A = sum (-1)^p |A_p|
    # This might be tournament-independent!
    # =====================================================================

    print("\n--- ALTERNATING SUM OF |A_p| ---")
    for n in [4, 5, 6, 7]:
        if n <= 6:
            N = 2**(n*(n-1)//2)
            rng = None
        else:
            N = 200
            rng = np.random.RandomState(42)

        chi_A_dist = Counter()
        for trial in range(N):
            if n <= 6:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)
            chi_A = 0
            for p in range(n):
                ap = enumerate_allowed_paths(A, n, p)
                chi_A += (-1)**p * len(ap)
            chi_A_dist[chi_A] += 1

        print(f"  n={n}: chi_A distribution: {dict(sorted(chi_A_dist.items()))}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
