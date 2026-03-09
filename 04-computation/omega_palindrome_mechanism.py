"""
Investigate WHY Paley T_7 has palindromic Omega dims.

Key: Paley T_7 has Omega = [7,21,42,63,63,42,21].
This is EXACTLY C(7,1), C(7,2), C(7,3), C(7,4)(?), C(7,4)(?), C(7,3)(?), C(7,2)(?).

Wait: C(7,1)=7, C(7,2)=21, C(7,3)=35 ≠ 42. So it's NOT just binomial coefficients.

Actually, for a tournament on n=7 vertices:
- |A_0| = 7 = C(7,1)
- |A_1| = 21 = 3*7 (each vertex has 3 out-neighbors for regular)
  For regular tournament, each vertex has out-degree 3, so |A_1| = 7*3 = 21 = C(7,2)
- |A_2| = Hamiltonian paths on 3-vertex subtournaments
  For regular: each 3-subset has exactly 1 transitive triple or 2 cycles
  |A_2| = sum over 3-subsets of #(directed 2-paths)
  = sum over 3-subsets of (if transitive: 2, if cyclic: 3)
  Wait, for Paley: all C(7,3)=35 triples... each has 1 transitive triple giving 2 paths
  or 1 cycle giving 3 paths. With c3=14 cycles: |A_2| = 2*(35-14) + 3*14 = 42 + 42 = 84?
  But we got |A_2|=63. Hmm.

Actually wait — |A_2| is the number of directed paths of length 2 (3 vertices), not the number
of Hamiltonian paths of 3-subsets. A directed 2-path is (a,b,c) with a→b→c, all distinct.
For regular tournament with out-degree 3: each vertex has 3 out-neighbors, each of those has
3 out-neighbors, but some go back. |A_2| = sum_v sum_{u: v→u} (out-degree of u not back to v).

Let me just compute this directly and see why dim(Omega_p) = what it is.

The big question: is dim(Omega_p) = C(7,p+1) for Paley T_7?
C(7,1)=7, C(7,2)=21, C(7,3)=35, C(7,4)=35, C(7,5)=21, C(7,6)=7, C(7,7)=1
But Omega = [7,21,42,63,63,42,21] ≠ [7,21,35,35,21,7,1]

So it's 42 = 2*21, 63 = 3*21. There's a multiplier pattern:
7*1, 7*3, 7*6, 7*9, 7*9, 7*6, 7*3
= 7 * (1, 3, 6, 9, 9, 6, 3)

The sequence 1,3,6,9,9,6,3 sums to 37. Hmm.

Actually: dim(Omega_p) for a self-complementary tournament should satisfy:
the duality σ: v → n-1-v maps T to T^op, and then the path complex of T^op is
"reversed". Does this create a duality on Omega?

For the GLMY path complex: the converse operation T → T^op reverses all edges.
Under path reversal, a p-path (v_0,...,v_p) in T maps to (v_p,...,v_0) in T^op.
This gives a chain isomorphism Omega_p(T) ≅ Omega_p(T^op) but with
alternating signs on the boundary map.

If T ≅ T^op (self-converse), this gives an AUTO-isomorphism of the chain complex.
But this doesn't directly give Poincare duality (which maps Omega_p to Omega_{n-1-p}).

For Poincare duality we'd need a map Omega_p → Omega_{n-1-p}^*, i.e.,
sending a p-path to a functional on (n-1-p)-paths. The natural candidate is:
(v_0,...,v_p) → the "complementary" path on V\{v_0,...,v_p}.

For a tournament where every subset has a unique Hamiltonian path (which is NOT true
in general), this would work. But even for Paley T_7, subsets can have multiple
Hamiltonian paths.
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict

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

def main():
    n = 7
    # Paley T_7
    A = np.zeros((n, n), dtype=int)
    qr = {1, 2, 4}
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % n) in qr:
                A[i][j] = 1

    print("=" * 70)
    print("Paley T_7: detailed |A_p| and dim(Omega_p) analysis")
    print("=" * 70)

    for p in range(n):
        paths = enumerate_allowed_paths(A, n, p)
        omega = compute_omega_basis(A, n, p, paths,
                                     enumerate_allowed_paths(A, n, p-1) if p > 0 else [])
        dim_om = omega.shape[1] if omega.ndim == 2 else 0

        # Count paths by vertex set
        by_vset = defaultdict(int)
        for path in paths:
            vset = frozenset(path)
            by_vset[vset] += 1

        # How many distinct vertex sets?
        n_vsets = len(by_vset)
        # How many paths per vertex set?
        ppvs = sorted(by_vset.values())
        ppvs_dist = defaultdict(int)
        for v in ppvs:
            ppvs_dist[v] += 1

        print(f"\np={p}: |A_p|={len(paths)}, dim(Omega_p)={dim_om}")
        print(f"  Vertex sets: {n_vsets} (= C({n},{p+1})={len(list(combinations(range(n), p+1)))})")
        print(f"  Paths per vertex set: {dict(sorted(ppvs_dist.items()))}")
        print(f"  Omega constraint rank: {len(paths) - dim_om}")

    # Part 2: Self-complementary check
    print("\n" + "=" * 70)
    print("Complement map analysis")
    print("=" * 70)
    # The self-converse map for Paley is sigma: i -> 2i mod 7 (or some automorphism)
    # Actually, path reversal maps T to T^op. For self-converse T, this is an automorphism.

    # Check: does vertex complement give a duality?
    # Complement of a p-path (v_0,...,v_p) is V\{v_0,...,v_p}, a set of n-1-p vertices.
    # Does the complement always have a UNIQUE Hamiltonian path in the subtournament?
    print("\nDoes each (n-1-p)-element complement have a unique Ham path in Paley?")
    for p in range(n):
        q = n - 1 - p
        if q < 0: continue
        paths_p = enumerate_allowed_paths(A, n, p)
        # Group by vertex set
        by_vset = defaultdict(list)
        for path in paths_p:
            vset = frozenset(path)
            by_vset[vset] += [path]

        # For each vertex set, look at complement
        comp_hampath_counts = []
        for vset, these_paths in by_vset.items():
            comp = sorted(set(range(n)) - vset)
            if len(comp) == 0:
                comp_hampath_counts.append(1)  # empty set, trivial
                continue
            # Count Ham paths of A[comp]
            sub_n = len(comp)
            if sub_n == 1:
                comp_hampath_counts.append(1)
                continue
            # Build sub-adjacency
            sub_A = np.zeros((sub_n, sub_n), dtype=int)
            for i, ci in enumerate(comp):
                for j, cj in enumerate(comp):
                    sub_A[i][j] = A[ci][cj]
            hp = len(enumerate_allowed_paths(sub_A, sub_n, sub_n - 1))
            comp_hampath_counts.append(hp)

        hp_dist = defaultdict(int)
        for c in comp_hampath_counts:
            hp_dist[c] += 1
        print(f"  p={p} (comp size {q+1}): comp Ham path counts: {dict(sorted(hp_dist.items()))}")

    # Part 3: Check dim(Omega_p) = dim(Omega_{n-2-p}) (shift by 1 from palindrome)
    print("\n" + "=" * 70)
    print("Palindrome pattern: dim(Omega_p) = dim(Omega_{n-2-p})?")
    print("=" * 70)
    print(f"Paley T_7 omega = [7, 21, 42, 63, 63, 42, 21]")
    print(f"Index: p=0,1,2,3,4,5,6")
    print(f"n-2-p: 5,4,3,2,1,0,-1")
    print(f"Omega_{n-2-p}: 42,63,63,42,21,7,_")
    print(f"NOT palindromic with shift. Let's check p vs n-1-p:")
    print(f"Omega_0 vs Omega_6: 7 vs 21 - NO")
    print(f"Omega_1 vs Omega_5: 21 vs 42 - NO")
    print(f"Omega_2 vs Omega_4: 42 vs 63 - NO")
    print(f"So the palindrome is in the LATTER half only?")
    print(f"Actually [7,21,42,63|63,42,21] - symmetry at the midpoint p=3!")
    print(f"dim(Omega_p) = dim(Omega_{6-p}) for p >= 3: 63=63, 42=42, 21=21")
    print(f"But p < 3 does NOT match: 7≠21, 21≠42, 42≠63")
    print(f"\nWait — it IS palindromic IF we list from p=0 to p=6:")
    print(f"[7, 21, 42, 63, 63, 42, 21]")
    print(f"Reverse: [21, 42, 63, 63, 42, 21, 7]")
    print(f"These are NOT the same! So it's NOT palindromic!")
    print(f"\nBut it looks like Omega_p = dim for p=3,...,6 mirrors the reversed sequence.")
    print(f"Actually the sequence [21, 42, 63, 63, 42, 21] from p=1 to p=6 IS palindromic.")
    print(f"And then Omega_0 = 7 is extra.")

    # Part 4: Compare number of allowed paths |A_p| across the 3 regular classes
    print("\n" + "=" * 70)
    print("Part 4: |A_p| for all 3 regular classes")
    print("=" * 70)

    rng = np.random.RandomState(42)
    classes = {}
    for trial in range(100000):
        B = random_tournament(n, rng)
        scores = tuple(sorted(sum(B[i][j] for j in range(n)) for i in range(n)))
        if scores != (3,3,3,3,3,3,3):
            continue
        # Compute H to identify class
        dp = {}
        for v in range(n):
            dp[(1 << v, v)] = 1
        for mask in range(1, 1 << n):
            bits = bin(mask).count('1')
            if bits == 1: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                prev_mask = mask ^ (1 << v)
                for u in range(n):
                    if (prev_mask & (1 << u)) and B[u][v]:
                        if (mask, v) not in dp:
                            dp[(mask, v)] = 0
                        dp[(mask, v)] += dp.get((prev_mask, u), 0)
        H = sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

        if H not in classes:
            Ap_sizes = []
            Om_dims = []
            for p in range(n):
                paths = enumerate_allowed_paths(B, n, p)
                omega = compute_omega_basis(B, n, p, paths,
                                             enumerate_allowed_paths(B, n, p-1) if p > 0 else [])
                Ap_sizes.append(len(paths))
                Om_dims.append(omega.shape[1] if omega.ndim == 2 else 0)
            classes[H] = (Ap_sizes, Om_dims)
            print(f"  H={H}: |A_p| = {Ap_sizes}")
            print(f"         Omega  = {Om_dims}")
            print(f"         constraints = {[a-o for a,o in zip(Ap_sizes, Om_dims)]}")
        if len(classes) == 3:
            break

if __name__ == '__main__':
    main()
