"""
c5_exact_lambda_formula.py -- kind-pasteur-2026-03-13-S61

Key identity: tr(A^5) = 5*c5_dir for any tournament (no non-Hamiltonian
closed 5-walks exist in tournaments due to arc uniqueness).

So c5_dir = (1/5) * tr(A^5).

Now, tr(A^5) = sum_i (A^5)_{ii} = sum of products of 5 arc indicators
along directed 5-cycles.

We want to express this in terms of lambda(u,v).

APPROACH: At n=5, enumerate all tournaments and build a lookup table
from lambda graph -> c5_dir. Then look for polynomial patterns.

At n=5: lambda(u,v) ranges from 0 to 3.
The sum S1 = sum lambda(u,v) = 3 * c3_vertex_set_count.

KNOWN IDENTITY for tournaments:
c3 = C(n,3) - sum_v C(s_v, 2) where s_v = out-degree

At n=5: c3_vertex_set = 10 - sum_v C(s_v,2) = 10 - sum_v s_v(s_v-1)/2

Also: sum_v s_v = C(n,2) = 10, sum_v s_v^2 = sum_v s_v * (s_v)
And sum_v C(s_v,2) = (sum s_v^2 - sum s_v) / 2 = (sum s_v^2 - 10) / 2

For REGULAR tournaments (s_v=2 for all): sum s_v^2 = 5*4 = 20, c3 = 10 - 5 = 5.

So S1 = 3*c3 = 3*(10 - (sum s_v^2 - 10)/2) = 30 - 3*(sum s_v^2 - 10)/2
     = 30 - 3*sum(s_v^2)/2 + 15 = 45 - 3*sum(s_v^2)/2

Now, at n=5 the number of directed 5-cycles c5_dir = tr(A^5)/5.

For the score sequence (1,2,2,2,3): c5 values are {1, 2, 3}.
Different values for the SAME score! So score does NOT determine c5.
But lambda DOES. So lambda carries more info than score (at n=5).

Let me find the formula by regression with more features.
"""

import numpy as np
from itertools import combinations, permutations
from numpy.linalg import lstsq

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

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def count_directed_hc(A, vset):
    combo = tuple(sorted(vset))
    k = len(combo)
    count = 0
    for perm in permutations(combo[1:]):
        path = (combo[0],) + perm
        valid = True
        for i in range(k):
            if A[path[i]][path[(i+1) % k]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

# Build n=5 data
n = 5
total_bits = n * (n-1) // 2

print("=" * 60)
print(f"EXACT c5_dir FORMULA IN TERMS OF LAMBDA (n={n})")
print("=" * 60)

all_data = []
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    c5 = count_directed_hc(A, range(n))

    # Lambda features
    lam_pairs = []
    for i in range(n):
        for j in range(i+1, n):
            lam_pairs.append(L[i][j])

    # Scores
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))

    all_data.append({
        'A': A.copy(),
        'L': L.copy(),
        'c5': c5,
        'lam': lam_pairs,
        'scores': scores,
    })

# Symmetric functions of lambda values
# Since we want a PERMUTATION-INVARIANT formula (same c5 regardless of vertex labeling),
# the formula should be in terms of SYMMETRIC functions of the lambda matrix.

# Symmetric functions of the lambda upper triangle:
# e1 = sum lambda[i][j] (= 3*c3)
# e2 = sum lambda[i][j]^2
# e3 = sum lambda[i][j]^3
# But also VERTEX-centered functions:
# d_v = sum_u lambda[v][u] (degree in lambda graph, = 2*delta_v = 2*c3_through_v)
# q_v = sum_u lambda[v][u]^2

# And TRIANGLE functions in the lambda graph:
# T = sum_{u<v<w} lambda[u][v]*lambda[v][w]*lambda[w][u] (kind of)

# Try: build a feature matrix with many symmetric lambda features
features_list = []
target_list = []

for d in all_data:
    L = d['L']
    lam = d['lam']
    c5 = d['c5']

    # Basic symmetric functions
    S1 = sum(lam)  # = 3*c3_vertex_sets
    S2 = sum(l**2 for l in lam)
    S3 = sum(l**3 for l in lam)

    # Vertex degrees in lambda graph
    d_v = [sum(L[v][u] for u in range(n) if u != v) for v in range(n)]
    D1 = sum(d_v)  # = 2*S1
    D2 = sum(dv**2 for dv in d_v)
    D3 = sum(dv**3 for dv in d_v)

    # Lambda-squared matrix trace: sum_{u,v} lambda[u][v]^2
    L2_trace = sum(L[u][v]**2 for u in range(n) for v in range(n) if u != v)

    # Lambda matrix product trace: tr(L^2) = sum_k L[i][k]*L[k][j]
    L2 = L @ L
    trL2 = int(np.trace(L2))

    # tr(L^3)
    L3 = L2 @ L
    trL3 = int(np.trace(L3))

    features = [1, S1, S2, S3, D2, trL2, trL3]
    features_list.append(features)
    target_list.append(c5)

X = np.array(features_list, dtype=float)
y = np.array(target_list, dtype=float)

# Try fitting
coeffs, _, _, _ = lstsq(X, y, rcond=None)
err = np.max(np.abs(X @ coeffs - y))
print(f"\nFeatures: [1, S1, S2, S3, D2, trL2, trL3]")
print(f"Coefficients: {[f'{c:.6f}' for c in coeffs]}")
print(f"Max abs error: {err:.6f}")
print(f"Perfect fit? {err < 0.001}")

# Try with fewer features
for feature_names, feature_indices in [
    ("[1, S1, S2]", [0, 1, 2]),
    ("[1, S1, S2, trL2]", [0, 1, 2, 5]),
    ("[1, S1, S2, S3, trL2]", [0, 1, 2, 3, 5]),
    ("[1, S1, trL2]", [0, 1, 5]),
    ("[1, S1, D2]", [0, 1, 4]),
    ("[1, S1, trL3]", [0, 1, 6]),
    ("[1, trL2]", [0, 5]),
    ("[1, S1, S2, D2]", [0, 1, 2, 4]),
]:
    X_sub = X[:, feature_indices]
    c_sub, _, _, _ = lstsq(X_sub, y, rcond=None)
    err_sub = np.max(np.abs(X_sub @ c_sub - y))
    if err_sub < 0.001:
        # Try to find rational coefficients
        for denom in range(1, 121):
            int_coeffs = np.round(c_sub * denom)
            err_int = np.max(np.abs(X_sub @ (int_coeffs / denom) - y))
            if err_int < 0.001:
                print(f"\n  EXACT FIT: {feature_names}")
                terms = []
                for i, idx in enumerate(feature_indices):
                    name = ["1", "S1", "S2", "S3", "D2", "trL2", "trL3"][idx]
                    c = int(int_coeffs[i])
                    if c != 0:
                        terms.append(f"({c})*{name}")
                print(f"  c5 = ({' + '.join(terms)}) / {denom}")
                print(f"  Error: {err_int:.8f}")
                break
        else:
            print(f"\n  Exact fit with {feature_names}, coeffs={[f'{c:.6f}' for c in c_sub]}")

# Now verify at n=6: does the formula hold per-5-vertex-set?
print(f"\n{'='*60}")
print(f"VERIFICATION AT n=6")
print(f"{'='*60}")

n6 = 6
total_bits_6 = n6 * (n6-1) // 2

# The formula for c5_dir at n=5 should be:
# c5 = f(lambda values on 5 vertices)
# At n=6, c5_dir(T) = sum over C(6,5) = 6 vertex sets of hc(T[V])
# Each term depends on the lambda values WITHIN V (restricted lambda).

# But our formula uses lambda(u,v) which at n=5 counts witnesses among 3 vertices.
# At n=6, the RESTRICTED lambda within V has 4 witness vertices.
# These are DIFFERENT from the full lambda!

# So we need to check: is the formula in terms of RESTRICTED lambda?
# Or in terms of the FULL lambda?

# The exhaustive check at n=6 showed c5_dir IS determined by the FULL lambda.
# But the FULL lambda at n=6 has C(6,2)=15 values.

# Let me check if the per-vertex-set formula (using restricted lambda)
# gives the right answer.
np.random.seed(42)
restricted_formula_works = True

for trial in range(500):
    bits = np.random.randint(0, 1 << total_bits_6)
    A = bits_to_adj(bits, n6)
    L = lambda_graph(A, n6)

    c5_total = 0
    c5_formula_total = 0

    for combo in combinations(range(n6), 5):
        # Direct count
        hc = count_directed_hc(A, combo)
        c5_total += hc

        # Restricted lambda
        L_restricted = np.zeros((5, 5), dtype=int)
        combo_list = list(combo)
        for i in range(5):
            for j in range(5):
                if i != j:
                    L_restricted[i][j] = lambda_graph(A[np.ix_(combo_list, combo_list)], 5)[i][j]

        # Use formula from n=5
        lam_r = []
        for i in range(5):
            for j in range(i+1, 5):
                lam_r.append(L_restricted[i][j])

        S1_r = sum(lam_r)
        S2_r = sum(l**2 for l in lam_r)
        S3_r = sum(l**3 for l in lam_r)
        D2_r = sum(sum(L_restricted[v][u] for u in range(5) if u != v)**2 for v in range(5))
        trL2_r = int(np.trace(L_restricted @ L_restricted))
        trL3_r = int(np.trace(L_restricted @ L_restricted @ L_restricted))

        feat = np.array([1, S1_r, S2_r, S3_r, D2_r, trL2_r, trL3_r], dtype=float)
        c5_predicted = feat @ coeffs

        c5_formula_total += c5_predicted

    if abs(c5_total - c5_formula_total) > 0.001:
        restricted_formula_works = False
        if trial < 3:
            print(f"  Trial {trial}: c5_direct={c5_total}, c5_formula={c5_formula_total:.2f}")

print(f"  Restricted lambda formula works at n=6? {restricted_formula_works}")

# If the restricted formula doesn't work, it's because the formula
# depends on the FULL lambda, not just the restricted lambda.
# In that case, we need a formula using the full lambda(u,v) at n=6.

print("\nDone.")
