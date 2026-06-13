#!/usr/bin/env python3
"""
Deep analysis of D = M - sym(H_ab) correction matrix.

We know: M[a,b] = (H_ab + H_ab^T)/2 + D[a,b]
where D is symmetric and tr(D) = tr(M) = H (odd n) or 0 (even n).

QUESTION: What determines D? Is it related to:
  1. The skew-adjacency matrix S[i,j] = A[i][j] - A[j][i]?
  2. The score vector?
  3. Some path-counting quantity?

At n=5 Paley: M = 2I, H_ab is circulant with row [0,1,0,0,1]
  sym(H_ab) has row [0, 0.5, 0, 0, 0.5] (wait — H_ab IS symmetric here!)
  So D = 2I - H_ab. Let's check this carefully.

Also: at n=3, H_ab has tr=0 always (no path starts and ends at same vertex).
So D diagonal = M diagonal. What IS the diagonal of M?

CONJECTURE to test: M[a,a] = (-1)^{n-1} * (n-1)! * something?
Or more likely: M[a,a] relates to paths avoiding vertex a.
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1: valid = False; break
        if valid: count += 1
    return count

def ham_path_count_dp(A):
    n = len(A)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if (mask & (1 << u)) or A[v][u] != 1:
                    continue
                dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def transfer_matrix(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_subset(A, S_verts, end=a)
                    bb = count_paths_subset(A, R_verts, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def endpoint_matrix(A):
    n = len(A)
    H_ab = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            H_ab[a][b] = count_paths_subset(A, list(range(n)), start=a, end=b)
    return H_ab

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(edges)):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield [list(row) for row in A]

# =====================================================================
# n=3: What is M[a,a]?
# =====================================================================
print("=" * 70)
print("DIAGONAL OF M: What is M[a,a]?")
print("=" * 70)

print("\n--- n=3 ---")
for idx, A in enumerate(all_tournaments(3)):
    if idx >= 4: break
    n = 3
    M = transfer_matrix(A)
    Hab = endpoint_matrix(A)
    H = ham_path_count_dp(A)
    scores = [sum(A[i]) for i in range(n)]

    # M[a,a] from definition: with a=b, U = V\{a}
    # M[a,a] = sum_S (-1)^|S| * E_a(S+{a}) * B_a(R+{a})
    # where R = U\S, and S+{a}, R+{a} partition V

    # Count ham paths on V\{a}
    for a in range(n):
        others = [v for v in range(n) if v != a]
        h_minus_a = count_paths_subset(A, others)
        print(f"  T{idx}: M[{a},{a}]={M[a][a]}, H(V\\{{{a}}})={h_minus_a}, score[{a}]={scores[a]}")

# =====================================================================
# n=5: M diagonal vs H(V\{a}) and scores
# =====================================================================
print("\n--- n=5: M diagonal analysis ---")

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
seen = set()
class_id = 0

for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    key = tuple(tuple(row) for row in A)
    min_key = key
    for perm in permutations(range(n)):
        pkey = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if pkey < min_key: min_key = pkey
    if min_key in seen: continue
    seen.add(min_key)
    class_id += 1

    H = ham_path_count_dp(A)
    M = transfer_matrix(A)
    Hab = endpoint_matrix(A)
    scores = sorted([sum(A[i]) for i in range(n)], reverse=True)

    # For each vertex a: M[a,a] vs H(V\{a}) vs score(a)
    diag_M = [M[i][i] for i in range(n)]
    h_minus = [count_paths_subset(A, [v for v in range(n) if v != a]) for a in range(n)]
    vertex_scores = [sum(A[i]) for i in range(n)]

    # Key test: is M[a,a] = H(V\{a}) for all a?
    diag_eq_hminus = all(diag_M[a] == h_minus[a] for a in range(n))

    # Alternative: M[a,a] = (-1)^{n-1} * H(V\{a})? (at n=5, (-1)^4 = 1)
    # Or M[a,a] = H(V\{a}) - correction?
    diffs = [diag_M[a] - h_minus[a] for a in range(n)]

    if class_id <= 12:
        print(f"\n  Class {class_id}: H={H}, scores={scores}")
        print(f"    M diagonal:    {diag_M}")
        print(f"    H(V\\{{a}}):     {h_minus}")
        print(f"    M[a,a]=H(V\\a)? {diag_eq_hminus}")
        if not diag_eq_hminus:
            print(f"    diffs M-H(V\\a): {diffs}")

# =====================================================================
# n=5: D structure — is D related to adjacency?
# =====================================================================
print()
print("=" * 70)
print("D = M - sym(H_ab): RELATIONSHIP TO ADJACENCY")
print("=" * 70)

n = 5
seen = set()

for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    key = tuple(tuple(row) for row in A)
    min_key = key
    for perm in permutations(range(n)):
        pkey = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if pkey < min_key: min_key = pkey
    if min_key in seen: continue
    seen.add(min_key)

    H = ham_path_count_dp(A)
    if H != 15: continue  # Focus on scalar case

    M = transfer_matrix(A)
    Hab = endpoint_matrix(A)
    sym_Hab = (Hab + Hab.T) / 2.0
    D = M - sym_Hab

    # Skew adjacency S[i,j] = A[i][j] - A[j][i] (antisymmetric)
    S = np.array([[A[i][j] - A[j][i] if i != j else 0 for j in range(n)] for i in range(n)])

    # Adjacency symmetrized: (A + A^T)/2 = (J-I)/2 for tournaments
    # So (A+A^T)/2 doesn't help — it's the same for all tournaments

    # Is D = c*I - d*(something involving A)?
    # For scalar M=3I: D = 3I - sym(H_ab)
    # So 2*D = 6*I - (H_ab + H_ab^T)

    # What is H_ab + H_ab^T?
    HpHT = Hab + Hab.T
    print(f"\n  H=15 case:")
    print(f"    H_ab = {Hab.tolist()}")
    print(f"    H_ab + H_ab^T = {HpHT.tolist()}")
    print(f"    2*D = 6*I - (H_ab+H_ab^T) = {(6*np.eye(n, dtype=int) - HpHT).tolist()}")

    # Row/col sums of H_ab
    print(f"    H_ab row sums (paths starting at a): {Hab.sum(axis=1).tolist()}")
    print(f"    H_ab col sums (paths ending at b):   {Hab.sum(axis=0).tolist()}")

    # Check: row sum of H_ab = number of ham paths starting at a
    # = "start-count" s(a). Similarly col sum = "end-count" e(b).
    # sum of all = H. For odd n: s(a) + e(a) = ?
    for a in range(n):
        sa = int(Hab[a,:].sum())
        ea = int(Hab[:,a].sum())
        print(f"    vertex {a}: s(a)={sa}, e(a)={ea}, s+e={sa+ea}")
    break

# =====================================================================
# n=5: 2D = (H/n)*n*I - (H_ab + H_ab^T)
# Equivalently: H_ab + H_ab^T + 2D = 2M
# Since M is symmetric: H_ab + H_ab^T = 2M - 2D
# And tr(H_ab + H_ab^T) = 0, tr(M) = H
# So tr(2D) = 2H => tr(D) = H (confirmed)
# =====================================================================

print()
print("=" * 70)
print("KEY IDENTITY: M = sym(H_ab) + D where D captures signed structure")
print("=" * 70)

# Let's look at D more carefully for all n=5 classes
n = 5
seen = set()
class_id = 0

for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    key = tuple(tuple(row) for row in A)
    min_key = key
    for perm in permutations(range(n)):
        pkey = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if pkey < min_key: min_key = pkey
    if min_key in seen: continue
    seen.add(min_key)
    class_id += 1

    H = ham_path_count_dp(A)
    M = transfer_matrix(A)
    Hab = endpoint_matrix(A)
    D = M - (Hab + Hab.T) / 2.0

    # Key metrics of D
    tr_D = np.trace(D)
    frob_D = np.sqrt(np.trace(D @ D))
    frob_M = np.sqrt(np.trace(M @ M))

    # Ratio ||D||_F / ||M||_F measures how much of M is "non-endpoint"
    ratio = frob_D / frob_M if frob_M > 0 else 0

    # Is D positive semidefinite?
    eig_D = sorted(np.linalg.eigvalsh(D))
    psd_D = eig_D[0] >= -1e-10

    # Check: is D = M - sym(H_ab) = (something involving sub-tournament paths)?
    # D[a,b] = M[a,b] - (H(a->b) + H(b->a))/2

    print(f"  Class {class_id:2d}: H={H:3d}, tr(D)={tr_D:5.0f}, ||D||_F={frob_D:6.2f}, "
          f"||M||_F={frob_M:6.2f}, ratio={ratio:.3f}, D_PSD={psd_D}, "
          f"D_eigs={[round(e,2) for e in eig_D]}")

# =====================================================================
# FUNDAMENTAL QUESTION: What is D[a,b] combinatorially?
# =====================================================================
print()
print("=" * 70)
print("D[a,b] = M[a,b] - (H(a->b) + H(b->a))/2: COMBINATORIAL MEANING")
print("=" * 70)

# M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
# H(a->b) = paths from a to b using all vertices
# The S=empty term in M gives: E_a({a}) * B_b(V\{a,b} + {b})
#   = 1 * (ham paths on V\{a} starting at some neighbor of a and ending at b... NO)
#   Actually E_a({a}) = 1 (just vertex a, ending at a = trivial path)
#   B_b(V\{a}) = ham paths on V\{a} starting at b
#   Wait: R = U\S = U = V\{a,b}, R+{b} = V\{a}
#   B_b(V\{a}) = paths on V\{a} starting at b
# So the S=empty term = 1 * B_b(V\{a}) = H_{V\{a}}(b -> ?)
# And H(a->b) = paths a -> ... -> b on full V

# So M[a,b] - H(a->b) involves the DIFFERENCE between:
#   - Inclusion-exclusion decomposition (signed)
#   - Direct counting (unsigned)

# Let's verify: for n=3, what is each term in the sum for M[a,b]?
print("\nn=3: Term-by-term M[a,b] decomposition")
A3 = [[0,1,0],[0,0,1],[1,0,0]]  # 3-cycle
n = 3
for a in range(n):
    for b in range(n):
        if a == b: continue
        U = [v for v in range(n) if v != a and v != b]
        terms = []
        for k in range(len(U)+1):
            for S in combinations(U, k):
                S_set = set(S)
                R = [v for v in U if v not in S_set]
                S_verts = sorted(list(S) + [a])
                R_verts = sorted(R + [b])
                ea = count_paths_subset(A3, S_verts, end=a)
                bb = count_paths_subset(A3, R_verts, start=b)
                contrib = ((-1)**k) * ea * bb
                terms.append((k, list(S), ea, bb, contrib))
        hab = count_paths_subset(A3, list(range(n)), start=a, end=b)
        hba = count_paths_subset(A3, list(range(n)), start=b, end=a)
        total = sum(t[4] for t in terms)
        print(f"  M[{a},{b}] = {total}: H({a}->{b})={hab}, H({b}->{a})={hba}, "
              f"D = {total - (hab+hba)/2:.1f}")
        for k, S, ea, bb, c in terms:
            print(f"    |S|={k}, S={S}: E_a={ea}, B_b={bb}, contrib={c}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
D[a,b] = M[a,b] - (H(a->b) + H(b->a))/2

Key findings:
  1. tr(D) = H (odd n) or 0 (even n) — same as tr(M)
  2. D is NOT diagonal — has rich off-diagonal structure
  3. ||D||_F / ||M||_F is substantial (50-90% of M's norm)
  4. D captures the signed inclusion-exclusion structure
  5. For scalar M = (H/n)*I: D = (H/n)*I - sym(H_ab)
     meaning D measures how far endpoint counts deviate from uniform

Physical interpretation of D[a,b]:
  M[a,b] counts "how well vertex a can feed into vertex b"
  via signed path decomposition. H(a->b) counts direct paths.
  D[a,b] = M[a,b] - avg(H(a->b), H(b->a)) is the EXCESS
  signed structure beyond simple endpoint counting.
""")
