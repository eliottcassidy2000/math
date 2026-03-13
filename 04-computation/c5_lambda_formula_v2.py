"""
c5_lambda_formula_v2.py -- kind-pasteur-2026-03-13-S61
Find exact polynomial formula for c5_dir in terms of lambda graph values.

Since THM-172 proves c5 IS lambda-determined, an exact formula MUST exist.
The question is: what degree polynomial in what features?

Strategy: at n=5, enumerate all 1024 tournaments, compute lambda + c5,
then do regression with increasingly rich feature sets until perfect fit.

Key: tr(A^5) = 5*c5_dir for tournaments.
Need: express tr(A^5) as f(lambda values).

Note: A[u][v]*A[v][w]*A[w][u] = 1 iff (u,v,w) is a directed 3-cycle.
lambda(u,v) = #{w : (u,v,w) or (v,u,w) is a directed 3-cycle}.
"""

import numpy as np
from itertools import combinations, permutations
from numpy.linalg import lstsq, matrix_rank

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

def count_directed_5cycles(A, n):
    """Count total directed 5-cycles = tr(A^5)/5"""
    A5 = np.linalg.matrix_power(A, 5)
    return int(np.trace(A5)) // 5

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
    c5 = count_directed_5cycles(A, n)

    all_data.append({
        'A': A.copy(),
        'L': L.copy(),
        'c5': c5,
    })

# Extract lambda features with increasing complexity
def compute_features(L, n, degree=1):
    """Compute symmetric polynomial features from lambda graph."""
    # Upper triangle values
    lam = []
    for i in range(n):
        for j in range(i+1, n):
            lam.append(L[i][j])

    # Power sums of lambda values
    Sk = [sum(l**k for l in lam) for k in range(1, 6)]  # S1..S5

    # Vertex lambda-degrees
    d_v = [sum(L[v][u] for u in range(n) if u != v) for v in range(n)]
    Dk = [sum(dv**k for dv in d_v) for k in range(1, 6)]  # D1..D5

    # Matrix traces
    Lpow = L.copy()
    traces = []
    for k in range(2, 6):
        Lpow = Lpow @ L
        traces.append(int(np.trace(Lpow)))  # trL2..trL5

    # Lambda histogram (how many pairs have each lambda value)
    from collections import Counter
    hist = Counter(lam)
    hist_features = [hist.get(k, 0) for k in range(n-1)]  # n0, n1, n2, n3

    # Build feature dict
    feats = {'1': 1}
    for k in range(5):
        feats[f'S{k+1}'] = Sk[k]
    for k in range(5):
        feats[f'D{k+1}'] = Dk[k]
    for k in range(4):
        feats[f'trL{k+2}'] = traces[k]
    for k in range(n-1):
        feats[f'n{k}'] = hist_features[k]

    if degree >= 2:
        # Add quadratic terms (products of degree-1 features)
        base_keys = ['S1', 'S2', 'S3', 'D2', 'D3', 'trL2', 'trL3']
        for i, ki in enumerate(base_keys):
            for j, kj in enumerate(base_keys):
                if j >= i:
                    feats[f'{ki}*{kj}'] = feats[ki] * feats[kj]

    return feats

# Try degree 1 features
print("\n--- Degree 1 Features ---")
feat_dicts = [compute_features(d['L'], n, degree=1) for d in all_data]
all_keys = list(feat_dicts[0].keys())

X = np.array([[fd[k] for k in all_keys] for fd in feat_dicts], dtype=float)
y = np.array([d['c5'] for d in all_data], dtype=float)

rank = matrix_rank(X)
coeffs, residuals, _, _ = lstsq(X, y, rcond=None)
pred = X @ coeffs
err = np.max(np.abs(pred - y))
print(f"Features ({len(all_keys)}): {all_keys}")
print(f"Matrix rank: {rank}")
print(f"Max abs error: {err:.6f}")
print(f"Perfect fit? {err < 0.001}")

# Try degree 2 features
print("\n--- Degree 2 Features ---")
feat_dicts2 = [compute_features(d['L'], n, degree=2) for d in all_data]
all_keys2 = list(feat_dicts2[0].keys())

X2 = np.array([[fd[k] for k in all_keys2] for fd in feat_dicts2], dtype=float)
rank2 = matrix_rank(X2)
coeffs2, _, _, _ = lstsq(X2, y, rcond=None)
pred2 = X2 @ coeffs2
err2 = np.max(np.abs(pred2 - y))
print(f"Features ({len(all_keys2)}): {all_keys2[:10]}...")
print(f"Matrix rank: {rank2}")
print(f"Max abs error: {err2:.6f}")
print(f"Perfect fit? {err2 < 0.001}")

if err2 < 0.001:
    # Find which features matter
    # Try removing features one at a time
    print("\n--- Feature Importance (backward elimination) ---")
    active = list(range(len(all_keys2)))
    for step in range(len(all_keys2)):
        best_remove = None
        best_err = float('inf')
        for i in active:
            trial = [j for j in active if j != i]
            X_trial = X2[:, trial]
            c_trial, _, _, _ = lstsq(X_trial, y, rcond=None)
            e_trial = np.max(np.abs(X_trial @ c_trial - y))
            if e_trial < best_err:
                best_err = e_trial
                best_remove = i
        if best_err < 0.001:
            active.remove(best_remove)
        else:
            break

    print(f"Minimal feature set ({len(active)} features):")
    minimal_keys = [all_keys2[i] for i in active]
    print(f"  {minimal_keys}")
    X_min = X2[:, active]
    c_min, _, _, _ = lstsq(X_min, y, rcond=None)
    err_min = np.max(np.abs(X_min @ c_min - y))
    print(f"  Max error: {err_min:.6f}")

    # Find rational coefficients
    for denom in range(1, 1201):
        int_coeffs = np.round(c_min * denom)
        err_rat = np.max(np.abs(X_min @ (int_coeffs / denom) - y))
        if err_rat < 0.001:
            print(f"\n  EXACT RATIONAL FORMULA (denom={denom}):")
            print(f"  c5_dir = (", end="")
            terms = []
            for i, k in enumerate(minimal_keys):
                c = int(int_coeffs[i])
                if c != 0:
                    terms.append(f"{c}*{k}")
            print(f"{' + '.join(terms)}) / {denom}")
            break

# Now try a completely different approach: use the FULL lambda matrix entries
# as features, grouped by sorted tuple (canonical form)
print(f"\n{'='*60}")
print("APPROACH 2: Group by canonical lambda, fit per-group")
print(f"{'='*60}")

# At n=5, two tournaments have same c5 iff same lambda (THM-172).
# So we can group by canonical lambda and check.
from collections import defaultdict
lambda_c5_map = defaultdict(set)
for d in all_data:
    L = d['L']
    # Canonical form: sorted upper triangle
    lam_tuple = tuple(L[i][j] for i in range(n) for j in range(i+1, n))
    lambda_c5_map[lam_tuple].add(d['c5'])

ambiguous = sum(1 for v in lambda_c5_map.values() if len(v) > 1)
print(f"Lambda groups: {len(lambda_c5_map)}")
print(f"Ambiguous: {ambiguous}")
print(f"c5 IS a function of lambda? {ambiguous == 0}")

# Distribution of c5 values
c5_vals = sorted(set(d['c5'] for d in all_data))
print(f"c5 range: {c5_vals}")

# What's the minimal info from lambda that determines c5?
# Try: just the MULTISET of lambda values (ignoring which pair has which value)
multiset_c5_map = defaultdict(set)
for d in all_data:
    L = d['L']
    lam_vals = tuple(sorted(L[i][j] for i in range(n) for j in range(i+1, n)))
    multiset_c5_map[lam_vals].add(d['c5'])

ambiguous_ms = sum(1 for v in multiset_c5_map.values() if len(v) > 1)
print(f"\nMultiset lambda groups: {len(multiset_c5_map)}")
print(f"Ambiguous: {ambiguous_ms}")
print(f"c5 determined by lambda MULTISET? {ambiguous_ms == 0}")

if ambiguous_ms > 0:
    for k, v in multiset_c5_map.items():
        if len(v) > 1:
            print(f"  Multiset {k} -> c5 values {sorted(v)}")

# What about just S1 and S2?
s1s2_c5_map = defaultdict(set)
for d in all_data:
    L = d['L']
    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    s1 = sum(lam)
    s2 = sum(l**2 for l in lam)
    s1s2_c5_map[(s1, s2)].add(d['c5'])

ambiguous_s1s2 = sum(1 for v in s1s2_c5_map.values() if len(v) > 1)
print(f"\n(S1,S2) groups: {len(s1s2_c5_map)}")
print(f"Ambiguous: {ambiguous_s1s2}")

# trL2 vs S2?
print(f"\n--- Relationship between S2 and trL2 ---")
for d in all_data[:10]:
    L = d['L']
    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    S2 = sum(l**2 for l in lam)
    trL2 = int(np.trace(L @ L))
    # trL2 = sum_{i,j} L[i][j]^2 = 2*S2 (since L is symmetric with 0 diagonal)
    print(f"  S2={S2}, trL2={trL2}, 2*S2={2*S2}, equal? {trL2 == 2*S2}")

# So trL2 = 2*S2. Let me check other redundancies.
# D1 = sum d_v = sum_v sum_u L[v][u] = 2*S1 (each pair counted twice)
print(f"\n--- Checking redundancies ---")
for d in all_data[:5]:
    L = d['L']
    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    S1 = sum(lam)
    d_v = [sum(L[v][u] for u in range(n) if u != v) for v in range(n)]
    D1 = sum(d_v)
    D2 = sum(dv**2 for dv in d_v)
    S2 = sum(l**2 for l in lam)
    trL2 = int(np.trace(L @ L))
    print(f"  S1={S1}, D1={D1}, 2*S1={2*S1} | S2={S2}, trL2={trL2}, 2*S2={2*S2} | D2={D2}")

# So: D1 = 2*S1, trL2 = 2*S2. These are redundant.
# Independent degree-1 features: S1, S2, S3, S4, S5 (power sums of lambda)
#   + D2, D3, D4, D5 (power sums of vertex degrees)
#   + trL3, trL4, trL5 (matrix traces)
# But many may be polynomial functions of each other.

# The REAL question: what is the ALGEBRAIC DIMENSION of the lambda space at n=5?
# = how many independent symmetric invariants are there?

# At n=5, lambda takes values in {0,1,2,3} for each of C(5,2)=10 pairs.
# The automorphism group S_5 acts on these.
# Number of orbits = number of distinct lambda graphs up to relabeling.

# But we need LABELED invariants (since c5 depends on labeled tournament).
# Actually no — c5_dir is the same for isomorphic tournaments.
# So c5_dir depends on the ISOMORPHISM CLASS of the lambda graph.

# Independent set: the ELEMENTARY SYMMETRIC POLYNOMIALS in the 10 lambda values
# e1, e2, ..., e10. But these require many terms.

# Better: try the POWER SUMS S1,...,S10 of the 10 lambda values.
# By Newton's identities, these determine the elementary symm polys.
# But e_k for k>10 is 0, so S1,...,S10 suffice.

# At n=5, lambda in {0,1,2,3}, so S_k for k >= 4 might be redundant.
# Actually S_k = sum l_i^k. Since l_i in {0,1,2,3}, we have at most 4^k growth.

# Let's check: are S1,...,S4 enough?
print(f"\n{'='*60}")
print("POWER SUM SUFFICIENCY CHECK")
print(f"{'='*60}")

for max_k in range(1, 11):
    groups = defaultdict(set)
    for d in all_data:
        L = d['L']
        lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
        key = tuple(sum(l**k for l in lam) for k in range(1, max_k+1))
        groups[key].add(d['c5'])
    amb = sum(1 for v in groups.values() if len(v) > 1)
    print(f"  S1..S{max_k}: {len(groups)} groups, {amb} ambiguous")
    if amb == 0:
        print(f"  --> S1..S{max_k} SUFFICE to determine c5!")
        break

# Now try vertex degree power sums
print(f"\n--- Adding vertex degree power sums ---")
for extra in ['D2', 'D3', 'D2+D3']:
    groups = defaultdict(set)
    for d in all_data:
        L = d['L']
        lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
        S1 = sum(lam)
        S2 = sum(l**2 for l in lam)
        d_v = [sum(L[v][u] for u in range(n) if u != v) for v in range(n)]
        D2 = sum(dv**2 for dv in d_v)
        D3 = sum(dv**3 for dv in d_v)
        if extra == 'D2':
            key = (S1, S2, D2)
        elif extra == 'D3':
            key = (S1, S2, D3)
        else:
            key = (S1, S2, D2, D3)
        groups[key].add(d['c5'])
    amb = sum(1 for v in groups.values() if len(v) > 1)
    print(f"  S1,S2,{extra}: {len(groups)} groups, {amb} ambiguous")

# Try matrix traces
print(f"\n--- Adding matrix traces ---")
for extra_name, extra_idx in [('trL3', 3), ('trL4', 4), ('trL5', 5)]:
    groups = defaultdict(set)
    for d in all_data:
        L = d['L']
        lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
        S1 = sum(lam)
        S2 = sum(l**2 for l in lam)
        Lpow = np.linalg.matrix_power(L, extra_idx)
        tr = int(np.trace(Lpow))
        key = (S1, S2, tr)
        groups[key].add(d['c5'])
    amb = sum(1 for v in groups.values() if len(v) > 1)
    print(f"  S1,S2,{extra_name}: {len(groups)} groups, {amb} ambiguous")

# Now try polynomial regression with enough features
# Since c5 is determined by the lambda multiset at n=5 (need to verify),
# or by S1..Sk for some k, let's find the minimal polynomial.

# First find which groups of S1,S2 are ambiguous
print(f"\n--- Ambiguous (S1,S2) groups ---")
s1s2_groups = defaultdict(list)
for d in all_data:
    L = d['L']
    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    S1 = sum(lam)
    S2 = sum(l**2 for l in lam)
    s1s2_groups[(S1, S2)].append(d)

for (s1, s2), members in sorted(s1s2_groups.items()):
    c5_vals = set(d['c5'] for d in members)
    if len(c5_vals) > 1:
        # Show lambda multisets for this group
        multisets = set()
        for d in members:
            L = d['L']
            ms = tuple(sorted(L[i][j] for i in range(n) for j in range(i+1, n)))
            multisets.add(ms)
        print(f"  S1={s1}, S2={s2}: c5 = {sorted(c5_vals)}, lambda multisets: {sorted(multisets)[:3]}...")

print("\nDone.")
