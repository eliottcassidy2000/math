"""
c5_formula_exact.py -- kind-pasteur-2026-03-13-S61
Find the exact polynomial c5 = f(S1, S2, S3) at n=5.

Known from v2 analysis:
- c5 IS determined by the lambda MULTISET (only 9 groups)
- (S1, S2, S3) suffice to determine c5 (but (S1, S2) don't)
- The ambiguous case: S1=12, S2=18 splits into c5=2 (S3=30) vs c5=3 (S3=36)

Since there are only 9 groups, we can fit c5 = a + b*S1 + c*S2 + d*S3
(or with cross terms if needed).
"""

import numpy as np
from itertools import combinations
from numpy.linalg import lstsq
from collections import defaultdict

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

# Build the 9-group lookup table
n = 5
total_bits = n * (n-1) // 2

groups = defaultdict(set)
group_data = defaultdict(list)
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A5 = np.linalg.matrix_power(A, 5)
    c5 = int(np.trace(A5)) // 5

    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    S1 = sum(lam)
    S2 = sum(l**2 for l in lam)
    S3 = sum(l**3 for l in lam)

    key = (S1, S2, S3)
    groups[key].add(c5)
    if len(group_data[key]) == 0:
        group_data[key].append({'c5': c5, 'lam_multiset': tuple(sorted(lam)),
                                'scores': tuple(sorted(int(sum(A[i])) for i in range(n)))})

print("=" * 60)
print(f"THE 9 GROUPS: (S1, S2, S3) -> c5")
print("=" * 60)

table = []
for key in sorted(groups.keys()):
    s1, s2, s3 = key
    c5_vals = groups[key]
    assert len(c5_vals) == 1, f"Ambiguous! {key} -> {c5_vals}"
    c5 = list(c5_vals)[0]
    table.append((s1, s2, s3, c5))
    print(f"  S1={s1:2d}, S2={s2:3d}, S3={s3:4d} -> c5={c5}")

# Now fit c5 = polynomial in S1, S2, S3
# Try linear first
X_linear = np.array([[1, s1, s2, s3] for s1, s2, s3, _ in table], dtype=float)
y_table = np.array([c5 for _, _, _, c5 in table], dtype=float)

coeffs_lin, _, _, _ = lstsq(X_linear, y_table, rcond=None)
err_lin = np.max(np.abs(X_linear @ coeffs_lin - y_table))
print(f"\nLinear fit: c5 = {coeffs_lin[0]:.6f} + {coeffs_lin[1]:.6f}*S1 + {coeffs_lin[2]:.6f}*S2 + {coeffs_lin[3]:.6f}*S3")
print(f"Max error: {err_lin:.6f}")

if err_lin < 0.001:
    # Find rational form
    for denom in range(1, 1201):
        int_c = np.round(coeffs_lin * denom)
        err_r = np.max(np.abs(X_linear @ (int_c / denom) - y_table))
        if err_r < 0.001:
            a, b, c, d = [int(x) for x in int_c]
            print(f"\nEXACT: c5 = ({a} + {b}*S1 + {c}*S2 + {d}*S3) / {denom}")
            # Verify on all 1024 tournaments
            verified = True
            for bits in range(1 << total_bits):
                A = bits_to_adj(bits, n)
                L = lambda_graph(A, n)
                A5 = np.linalg.matrix_power(A, 5)
                c5_true = int(np.trace(A5)) // 5
                lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
                s1 = sum(lam)
                s2 = sum(l**2 for l in lam)
                s3 = sum(l**3 for l in lam)
                c5_pred = (a + b*s1 + c*s2 + d*s3) / denom
                if abs(c5_pred - c5_true) > 0.001:
                    verified = False
                    print(f"  FAIL: bits={bits}, c5_true={c5_true}, c5_pred={c5_pred:.4f}")
                    break
            print(f"Verified on all 1024 tournaments? {verified}")
            break
else:
    # Try with cross terms
    print("\nLinear doesn't fit. Trying quadratic...")
    # Features: 1, S1, S2, S3, S1^2, S1*S2, S1*S3, S2^2, S2*S3, S3^2
    X_quad = np.array([[1, s1, s2, s3, s1**2, s1*s2, s1*s3, s2**2, s2*s3, s3**2]
                       for s1, s2, s3, _ in table], dtype=float)
    coeffs_q, _, _, _ = lstsq(X_quad, y_table, rcond=None)
    err_q = np.max(np.abs(X_quad @ coeffs_q - y_table))
    print(f"Quadratic fit max error: {err_q:.6f}")

    if err_q < 0.001:
        # Feature selection
        feat_names = ['1', 'S1', 'S2', 'S3', 'S1^2', 'S1*S2', 'S1*S3', 'S2^2', 'S2*S3', 'S3^2']
        for denom in range(1, 1201):
            int_c = np.round(coeffs_q * denom)
            err_r = np.max(np.abs(X_quad @ (int_c / denom) - y_table))
            if err_r < 0.001:
                print(f"\nEXACT QUADRATIC (denom={denom}):")
                terms = []
                for i, name in enumerate(feat_names):
                    c = int(int_c[i])
                    if c != 0:
                        terms.append(f"{c}*{name}")
                print(f"  c5 = ({' + '.join(terms)}) / {denom}")
                break

# Also check at n=5: what is the relationship between S1, S2, S3 and scores?
print(f"\n{'='*60}")
print("RELATIONSHIP TO SCORE SEQUENCE")
print(f"{'='*60}")

score_groups = defaultdict(set)
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    L = lambda_graph(A, n)
    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    S1 = sum(lam)
    S2 = sum(l**2 for l in lam)
    S3 = sum(l**3 for l in lam)
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
    score_groups[scores].add((S1, S2, S3, c5))

for score, vals in sorted(score_groups.items()):
    if len(vals) > 1:
        print(f"  Score {score}: {len(vals)} distinct (S1,S2,S3,c5) tuples")
        for v in sorted(vals):
            print(f"    S1={v[0]}, S2={v[1]}, S3={v[2]}, c5={v[3]}")
    else:
        v = list(vals)[0]
        print(f"  Score {score}: S1={v[0]}, S2={v[1]}, S3={v[2]}, c5={v[3]}")

# Now verify the formula at n=6 using RESTRICTED lambda
# For each 5-vertex subset V of [6], compute restricted lambda and c5
# Then sum over all subsets
print(f"\n{'='*60}")
print("VERIFICATION AT n=6")
print(f"{'='*60}")

from itertools import permutations

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

n6 = 6
total_bits_6 = n6 * (n6-1) // 2

# We found that c5 = f(S1, S2, S3) with the exact formula above.
# At n=6, the total c5_dir = sum over C(6,5)=6 subsets V of hc(T[V]).
# For each V, the RESTRICTED lambda uses 3 witness vertices.
# This is the same as lambda at n=5.
# So c5_per_V = f(S1_V, S2_V, S3_V) where S_k uses restricted lambda.

# Let's check if this works.
if err_lin < 0.001:
    print("Testing linear formula c5 = f(S1,S2,S3) on restricted lambda at n=6...")
    np.random.seed(42)
    success = 0
    total_tests = 500

    for trial in range(total_tests):
        bits = np.random.randint(0, 1 << total_bits_6)
        A = bits_to_adj(bits, n6)

        c5_total = 0
        c5_formula_total = 0.0

        for combo in combinations(range(n6), 5):
            # Direct count
            hc = count_directed_hc(A, combo)
            c5_total += hc

            # Restricted lambda and formula
            combo_list = list(combo)
            A_sub = A[np.ix_(combo_list, combo_list)]
            L_sub = lambda_graph(A_sub, 5)
            lam_r = [L_sub[i][j] for i in range(5) for j in range(i+1, 5)]
            s1_r = sum(lam_r)
            s2_r = sum(l**2 for l in lam_r)
            s3_r = sum(l**3 for l in lam_r)

            c5_pred = (a + b*s1_r + c*s2_r + d*s3_r) / denom
            c5_formula_total += c5_pred

        if abs(c5_total - c5_formula_total) < 0.001:
            success += 1
        elif trial < 5:
            print(f"  Trial {trial}: c5_direct={c5_total}, c5_formula={c5_formula_total:.4f}")

    print(f"  Success: {success}/{total_tests}")
    print(f"  Restricted lambda formula works at n=6? {success == total_tests}")

# Now the BIG question: does c5_dir at n=6 equal a function of the FULL lambda?
# We know from THM-172 that it does. But can we express it?
# c5_total(T) = sum_{V in C(n,5)} f(S1_V, S2_V, S3_V)
# where S_k(V) = sum_{u<v in V} lambda_V(u,v)^k
# and lambda_V(u,v) = #{w in V\{u,v} : (u,v,w) or (v,u,w) is directed 3-cycle}

# This is NOT the same as the full lambda!
# lambda_V(u,v) <= lambda(u,v) since fewer witnesses.
# lambda_V(u,v) = #{w in V\{u,v} : w contributes to lambda(u,v)}

# At n=5: lambda_V = full lambda (V is everything).
# At n=6: for V = [6]\{k}, lambda_V(u,v) = lambda(u,v) - [k is witness for (u,v)]

# So lambda_V(u,v) = lambda(u,v) - indicator(k witnesses (u,v))
# where "k witnesses (u,v)" means (u,v,k) or (v,u,k) is a directed 3-cycle.

# This connects the restricted lambda to the FULL lambda in a beautiful way!

print(f"\n{'='*60}")
print("RESTRICTED vs FULL LAMBDA RELATIONSHIP")
print(f"{'='*60}")

np.random.seed(42)
for trial in range(3):
    bits = np.random.randint(0, 1 << total_bits_6)
    A = bits_to_adj(bits, n6)
    L_full = lambda_graph(A, n6)

    print(f"\nTrial {trial}:")
    for k in range(n6):
        combo = [v for v in range(n6) if v != k]
        A_sub = A[np.ix_(combo, combo)]
        L_sub = lambda_graph(A_sub, 5)

        # Check: lambda_V(u,v) = lambda_full(u,v) - [k witnesses (u,v)]
        for i in range(5):
            for j in range(i+1, 5):
                u, v = combo[i], combo[j]
                lam_full = L_full[u][v]
                lam_restricted = L_sub[i][j]
                # Does k witness (u,v)?
                k_witnesses = int(
                    (A[u][v] and A[v][k] and A[k][u]) or
                    (A[v][u] and A[u][k] and A[k][v])
                )
                expected = lam_full - k_witnesses
                if lam_restricted != expected:
                    print(f"  FAIL: k={k}, u={u}, v={v}: restricted={lam_restricted}, full-witness={expected}")

    print(f"  lambda_V(u,v) = lambda(u,v) - witness(k,u,v): VERIFIED for all pairs")

print("\nDone.")
