"""
Analysis of beta_1+beta_5 coexistence and beta_4>0 at n=8.

Questions:
1. What does the beta_1=beta_5=1 tournament look like?
2. Is there a pattern in beta_4>0 tournaments (scores, c3, structure)?
3. Does beta_2=0 have an algebraic explanation?
4. Onset dimensions: beta_0 always, beta_1 at n=3, beta_3 at n=6,
   beta_4 at n=8, beta_5 at n=8. What's the pattern?
"""
import numpy as np
from math import comb
from itertools import combinations

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
    rank = sum(s > 1e-10 for s in S)
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

def betti_single(A, n, target_p):
    allowed = {}
    for p in [target_p - 1, target_p, target_p + 1]:
        if p < 0: allowed[p] = []
        else: allowed[p] = enumerate_allowed_paths(A, n, p)
    omega_p = compute_omega_basis(A, n, target_p, allowed[target_p], allowed[target_p-1])
    dim_om = omega_p.shape[1] if omega_p.ndim == 2 else 0
    if dim_om == 0: return 0
    bd_p = build_boundary_matrix(allowed[target_p], allowed[target_p-1])
    bd_p_om = bd_p @ omega_p
    if bd_p_om.size > 0:
        sv = np.linalg.svd(bd_p_om, compute_uv=False)
        rk = sum(s > 1e-8 for s in sv)
    else: rk = 0
    ker = dim_om - rk
    omega_p1 = compute_omega_basis(A, n, target_p+1, allowed[target_p+1], allowed[target_p])
    dim_om1 = omega_p1.shape[1] if omega_p1.ndim == 2 else 0
    if dim_om1 > 0:
        bd1 = build_boundary_matrix(allowed[target_p+1], allowed[target_p])
        bd1_om = bd1 @ omega_p1
        sv1 = np.linalg.svd(bd1_om, compute_uv=False)
        im = sum(s > 1e-8 for s in sv1)
    else: im = 0
    return ker - im

def count_directed_3cycles(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: c3 += 1
        if A[i][k] and A[k][j] and A[j][i]: c3 += 1
    return c3

def is_strongly_connected(A, n):
    # BFS from 0 using forward edges
    visited_fwd = set()
    queue = [0]
    visited_fwd.add(0)
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[v][u] and u not in visited_fwd:
                visited_fwd.add(u)
                queue.append(u)
    if len(visited_fwd) != n: return False
    # BFS from 0 using backward edges
    visited_bwd = set()
    queue = [0]
    visited_bwd.add(0)
    while queue:
        v = queue.pop(0)
        for u in range(n):
            if A[u][v] and u not in visited_bwd:
                visited_bwd.add(u)
                queue.append(u)
    return len(visited_bwd) == n

def main():
    print("=" * 70)
    print("BETA COEXISTENCE AND ONSET ANALYSIS")
    print("=" * 70)

    n = 8

    # Part 1: Find and analyze all exotic Betti profiles
    print("\n--- Part 1: Finding exotic profiles at n=8 (2000 samples) ---", flush=True)

    rng = np.random.RandomState(11111)
    exotic = []

    for trial in range(2000):
        A = random_tournament(n, rng)

        # Quick screen: compute beta_1, beta_3 first
        bettis = [1]  # beta_0 = 1 always
        for p in range(1, n):
            b = betti_single(A, n, p)
            bettis.append(b)

        nonzero = [p for p in range(1, n) if bettis[p] > 0]
        if len(nonzero) > 1 or any(bettis[p] > 0 for p in [2, 4, 5, 6]):
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            c3 = count_directed_3cycles(A, n)
            sc = is_strongly_connected(A, n)
            chi = sum((-1)**p * bettis[p] for p in range(n))
            exotic.append({
                'trial': trial, 'bettis': bettis, 'scores': scores,
                'c3': c3, 'sc': sc, 'chi': chi
            })
            print(f"  Trial {trial}: betti={bettis}, scores={scores}, c3={c3}, SC={sc}, chi={chi}", flush=True)

        if (trial+1) % 500 == 0:
            print(f"  {trial+1}/2000 done, {len(exotic)} exotic found", flush=True)

    print(f"\n  Total exotic: {len(exotic)}/2000", flush=True)

    # Part 2: Analysis of exotic cases
    print("\n--- Part 2: Exotic case analysis ---", flush=True)

    # Group by profile type
    from collections import Counter
    type_counts = Counter()
    for e in exotic:
        nonzero = tuple(p for p in range(1, n) if e['bettis'][p] > 0)
        type_counts[nonzero] += 1

    print("  Profile types:")
    for nonzero, cnt in type_counts.most_common():
        print(f"    Nonzero at p={nonzero}: {cnt} cases")

    # Score analysis for beta_4>0
    b4_cases = [e for e in exotic if e['bettis'][4] > 0]
    if b4_cases:
        print(f"\n  beta_4>0 cases ({len(b4_cases)}):")
        for e in b4_cases:
            print(f"    betti={e['bettis']}, scores={e['scores']}, c3={e['c3']}, SC={e['sc']}")

    # Score analysis for beta_1+beta_5
    b15_cases = [e for e in exotic if e['bettis'][1] > 0 and e['bettis'][5] > 0]
    if b15_cases:
        print(f"\n  beta_1+beta_5 cases ({len(b15_cases)}):")
        for e in b15_cases:
            print(f"    betti={e['bettis']}, scores={e['scores']}, c3={e['c3']}, SC={e['sc']}")

    # Part 3: Onset dimensions table
    print("\n--- Part 3: Onset dimensions ---", flush=True)
    print("""
    beta_0: n=1 (trivial, always 1)
    beta_1: n=3 (first 3-cycle)
    beta_2: NEVER (conjecture: beta_2=0 for ALL tournaments)
    beta_3: n=6 (first appearance, ~1% at n=6)
    beta_4: n=8 (first appearance, ~0.5% at n=8)
    beta_5: n=8 (first appearance, ~0.2% at n=8, WITH beta_1)

    Pattern: onset of beta_p at n = ?
    p:     0  1  2  3  4  5
    onset: 1  3  -  6  8  8

    The gap at beta_2 is striking. beta_2=0 is the DEEPEST invariant.
    The seesaw mechanism (THM-095) depends on it.
    """)

    # Part 4: Is beta_2=0 related to tournament completeness?
    # In a tournament, ALL edges exist (complete digraph).
    # Omega_1 = all directed edges = n(n-1)/2 (one per pair).
    # Actually Omega_1 = n(n-1) (both orderings of each pair contribute IF edge exists).
    # No wait - Omega_1 = A_1 = directed edges = one per pair = C(n,2) in a tournament.
    # Wait - for a tournament, (i,j) is allowed iff i->j. So |A_1| = C(n,2) = n(n-1)/2.
    # But every pair contributes exactly one direction. So |A_1| = C(n,2).
    #
    # The constraint matrix for Omega_1 involves faces of 1-paths:
    # Boundary of (i,j) = (j) - (i). These are vertices. All vertices are in A_0.
    # So there are NO non-allowed faces. Hence Omega_1 = A_1.
    # dim(Omega_1) = C(n,2).
    #
    # For Omega_2: a 2-path (a,b,c) has faces (b,c), (a,c), (a,b).
    # (b,c) and (a,b) are always in A_1 (they're edges of the tournament).
    # (a,c): in A_1 iff a->c. If c->a, then (a,c) NOT in A_1.
    # So the Omega_2 constraint is: for each (a,b,c) with c->a,
    # the coefficient of (a,b,c) in any Omega_2 element must satisfy:
    # the "non-allowed face constraint" at (a,c).
    #
    # The constraint says: for each non-allowed face f,
    # sum of (-1)^i * coefficient where deleting position i gives f = 0.
    #
    # For face (a,c) from path (a,b,c): this comes from deleting position 1 (middle).
    # Sign = (-1)^1 = -1. So the constraint is: -coefficient(a,b,c) = 0
    # for each (a,b,c) where c->a.
    #
    # BUT WAIT: multiple paths can have the same non-allowed face (a,c).
    # If c->a, then for ANY b with a->b and b->c, the path (a,b,c)
    # contributes to the constraint at face (a,c).
    # The constraint: sum_{b: a->b, b->c} (-1)^1 * coeff(a,b,c) = 0
    # i.e., sum_{b: a->b, b->c} coeff(a,b,c) = 0
    #
    # This is the constraint: the sum of coefficients of all 2-paths from a to c
    # (going through various b) must be zero, whenever c->a.

    print("\n--- Part 4: Constraint structure for Omega_2 ---", flush=True)

    for n_test in [5, 6, 7]:
        rng2 = np.random.RandomState(42)
        A = random_tournament(n_test, rng2)

        a1 = enumerate_allowed_paths(A, n_test, 1)
        a2 = enumerate_allowed_paths(A, n_test, 2)
        a1_set = set(a1)

        # Count non-allowed 1-faces of 2-paths
        na_faces = set()
        for path in a2:
            a, b, c = path
            face = (a, c)
            if face not in a1_set:
                na_faces.add(face)

        # Each non-allowed face (a,c) means c->a in the tournament
        # Number of constraints = number of (a,c) pairs with c->a
        # that appear as a skip in some 2-path
        print(f"  n={n_test}: |A_2|={len(a2)}, non-allowed faces={len(na_faces)}")
        print(f"    These are (a,c) pairs with c->a AND some b with a->b->c")

        # How many 2-paths per non-allowed face?
        paths_per_face = {}
        for path in a2:
            a, b, c = path
            face = (a, c)
            if face not in a1_set:
                paths_per_face[face] = paths_per_face.get(face, 0) + 1

        if paths_per_face:
            vals = list(paths_per_face.values())
            print(f"    Paths per NA face: min={min(vals)}, max={max(vals)}, avg={np.mean(vals):.1f}")

    # Part 5: The n=3 case for beta_2 — why impossible?
    print("\n--- Part 5: Why beta_2=0 at n=3 (exhaustive) ---", flush=True)
    for bits in range(8):
        A = np.zeros((3, 3), dtype=int)
        idx = 0
        for i in range(3):
            for j in range(i+1, 3):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        a2 = enumerate_allowed_paths(A, 3, 2)
        a1 = enumerate_allowed_paths(A, 3, 1)
        dim_om2 = len(a2)  # Quick: how big is Omega_2?

        ob2 = compute_omega_basis(A, 3, 2, a2, a1)
        dim_om2_real = ob2.shape[1] if ob2.ndim == 2 else 0

        # beta_2 needs dim(Omega_3), but n=3 so A_3 = Hamiltonian paths
        a3 = enumerate_allowed_paths(A, 3, 3)  # would need 4 vertices, impossible at n=3

        b2 = betti_single(A, 3, 2)
        scores = tuple(sorted([sum(A[i]) for i in range(3)]))
        print(f"    bits={bits}: scores={scores}, |A_2|={len(a2)}, dim(Omega_2)={dim_om2_real}, |A_3|={len(a3)}, beta_2={b2}")

if __name__ == '__main__':
    main()
