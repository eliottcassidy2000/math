"""
c5_minimal_formula.py -- kind-pasteur-2026-03-13-S61

Extract the EXACT minimal formula for c5_dir at n=5 in terms of lambda.
Then test whether the per-vertex-set formula extends to n=6, 7, 8.

From c5_algebraic_v2: minimal feature set is 8 features from degree <= 3
symmetric lambda invariants. Let me find the simplest exact formula
and verify it generalizes.
"""

import numpy as np
from itertools import combinations, permutations
from numpy.linalg import lstsq, matrix_rank
from collections import defaultdict, Counter
from fractions import Fraction

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

n = 5
total_bits = n * (n-1) // 2
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]

print("=" * 60)
print(f"EXACT c5 FORMULA AT n={n}")
print("=" * 60)

# The 9 groups: (n0, n1, n2, n3) -> c5
groups = [
    (10, 0, 0, 0, 0),
    (7, 3, 0, 0, 0),
    (5, 4, 1, 0, 0),
    (3, 5, 2, 0, 1),
    (3, 6, 0, 1, 1),
    (1, 6, 3, 0, 2),
    (0, 9, 0, 1, 3),
    (2, 4, 4, 0, 1),
    (0, 5, 5, 0, 2),
]

# Since there are only 9 points and 3 free variables (n1,n2,n3),
# the unique interpolating polynomial has degree at most 8.
# But with 3 variables, we need to be smart about the basis.

# Key observation: c5 = f(n0,n1,n2,n3) with n0+n1+n2+n3=10.
# Using the constraint, work in (n1,n2,n3) space.

# From the table, notice:
# c5 = 0 when c3 = S1/3 <= 2 (S1 <= 6)
# c5 >= 1 when c3 >= 3

# At c3 = S1/3: S1 = n1+2*n2+3*n3, so c3 = (n1+2*n2+3*n3)/3.

# Actually, let me just solve the 9x9 system.
# Use monomial basis: 1, n1, n2, n3, n1*n2, n1*n3, n2*n3, n2^2, n3^2
# (9 basis functions for 9 equations)

monomial_names = ['1', 'n1', 'n2', 'n3', 'n1*n2', 'n1*n3', 'n2*n3', 'n2^2', 'n3^2']

def monomial_eval(n1, n2, n3):
    return [1, n1, n2, n3, n1*n2, n1*n3, n2*n3, n2**2, n3**2]

X = np.array([monomial_eval(n1, n2, n3) for _, n1, n2, n3, _ in groups], dtype=float)
y = np.array([c5 for _, _, _, _, c5 in groups], dtype=float)

rank = matrix_rank(X)
print(f"Monomial basis rank: {rank}")

if rank < 9:
    print(f"Rank deficient! Need different basis.")
    # Try with more monomials
    monomial_names2 = ['1', 'n1', 'n2', 'n3', 'n1*n2', 'n1*n3', 'n2*n3',
                       'n2^2', 'n3^2', 'n1^2', 'n1*n2*n3']
    def monomial_eval2(n1, n2, n3):
        return [1, n1, n2, n3, n1*n2, n1*n3, n2*n3, n2**2, n3**2, n1**2, n1*n2*n3]

    X2 = np.array([monomial_eval2(n1, n2, n3) for _, n1, n2, n3, _ in groups], dtype=float)
    rank2 = matrix_rank(X2)
    print(f"Extended basis rank: {rank2}")

    coeffs2, _, _, _ = lstsq(X2, y, rcond=None)
    err2 = np.max(np.abs(X2 @ coeffs2 - y))
    print(f"Max error: {err2:.8f}")

    if err2 < 0.001:
        for denom in range(1, 5001):
            int_c = np.round(coeffs2 * denom)
            err_r = np.max(np.abs(X2 @ (int_c / denom) - y))
            if err_r < 0.001:
                print(f"\nEXACT FORMULA (denom={denom}):")
                terms = []
                for i, name in enumerate(monomial_names2):
                    c = int(int_c[i])
                    if c != 0:
                        terms.append(f"{c}*{name}")
                print(f"  c5 = ({' + '.join(terms)}) / {denom}")

                # Verify on all groups
                for n0, n1, n2, n3, c5_true in groups:
                    feat = monomial_eval2(n1, n2, n3)
                    c5_pred = sum(int_c[i]*feat[i] for i in range(len(feat))) / denom
                    status = "OK" if abs(c5_pred - c5_true) < 0.001 else "FAIL"
                    print(f"  ({n0},{n1},{n2},{n3}): c5={c5_true}, pred={c5_pred:.4f} {status}")
                break
else:
    coeffs, _, _, _ = lstsq(X, y, rcond=None)
    err = np.max(np.abs(X @ coeffs - y))
    print(f"Max error: {err:.8f}")

    if err < 0.001:
        for denom in range(1, 5001):
            int_c = np.round(coeffs * denom)
            err_r = np.max(np.abs(X @ (int_c / denom) - y))
            if err_r < 0.001:
                print(f"\nEXACT FORMULA (denom={denom}):")
                terms = []
                for i, name in enumerate(monomial_names):
                    c = int(int_c[i])
                    if c != 0:
                        terms.append(f"{c}*{name}")
                print(f"  c5 = ({' + '.join(terms)}) / {denom}")
                break

# Alternative: express as function of (c3, S2, S3) using exact rational arithmetic
print(f"\n{'='*60}")
print("EXACT FORMULA VIA (c3, S2) + S3 CORRECTION")
print(f"{'='*60}")

# c3 = S1/3 determines c5 except at c3=4, S2=18 where S3 discriminates.
# c3=0: c5=0, c3=1: c5=0, c3=2: c5=0, c3=3: c5=1
# c3=4: c5 = {1 if S2=20, 2 if (S2=18,S3=30), 3 if (S2=18,S3=36)}
# c3=5: c5=2

# This suggests a Lagrange-type formula on the discrete variable c3 plus correction.
# Let's compute the exact Lagrange polynomial in c3, plus correction at c3=4.

# At c3=4: the three sub-cases have:
#   (S2, S3) = (20, 36): c5=1
#   (S2, S3) = (18, 30): c5=2
#   (S2, S3) = (18, 36): c5=3
# The S3 values: 30 (c5=2) vs 36 (c5=3). Difference in c5 = 1, difference in S3 = 6.
# At S2=18: c5 = 2 + (S3-30)/6 = (S3-18)/6. Check: S3=30 -> c5=2, S3=36 -> c5=3. YES!
# At S2=20: c5=1. Note S3=36 here too. So S3 alone doesn't determine c5 within c3=4.

# Using (S1, S2, S3) = ((n1+2n2+3n3), (n1+4n2+9n3), (n1+8n2+27n3)):
# The Vandermonde inverse gives:
# n0 = (10*6 - 6*S1 + 3*S2 - S3) / 6
# n1 = (-6*0 + 11*S1 - 7*S2 + 2*S3) / 6  ... actually let me compute properly

# Vandermonde: [1 1 1 1] [n0]   [10]
#              [0 1 2 3] [n1] = [S1]
#              [0 1 4 9] [n2]   [S2]
#              [0 1 8 27][n3]   [S3]

V = np.array([[1,1,1,1],[0,1,2,3],[0,1,4,9],[0,1,8,27]], dtype=float)
V_inv = np.linalg.inv(V)
print("Vandermonde inverse:")
print(V_inv)

# Check:
for _, n1, n2, n3, c5 in groups:
    n0 = 10 - n1 - n2 - n3
    rhs = [10, n1+2*n2+3*n3, n1+4*n2+9*n3, n1+8*n2+27*n3]
    nk = V_inv @ rhs
    assert abs(nk[0]-n0) < 0.01 and abs(nk[1]-n1) < 0.01 and abs(nk[2]-n2) < 0.01 and abs(nk[3]-n3) < 0.01

print("Vandermonde inverse verified.")

# Express n0,...,n3 in terms of S1, S2, S3 (with 10 = sum):
# n0 = (60 - 11*S1 + 3*S2 - (1/3)*S3) / 2  ... let me use the inverse
for row in range(4):
    coeffs_row = V_inv[row]
    print(f"  n{row} = {coeffs_row[0]:.4f}*10 + {coeffs_row[1]:.4f}*S1 + {coeffs_row[2]:.4f}*S2 + {coeffs_row[3]:.4f}*S3")

# Now find exact rational Vandermonde inverse
from fractions import Fraction
V_frac = [[Fraction(int(V[i][j])) for j in range(4)] for i in range(4)]

# Invert by Gauss elimination
aug = [row + [Fraction(1) if i == j else Fraction(0) for j in range(4)] for i, row in enumerate(V_frac)]
for col in range(4):
    # Find pivot
    for row in range(col, 4):
        if aug[row][col] != 0:
            aug[col], aug[row] = aug[row], aug[col]
            break
    pivot = aug[col][col]
    for j in range(8):
        aug[col][j] /= pivot
    for row in range(4):
        if row != col:
            factor = aug[row][col]
            for j in range(8):
                aug[row][j] -= factor * aug[col][j]

V_inv_frac = [[aug[i][j+4] for j in range(4)] for i in range(4)]
print("\nExact Vandermonde inverse (rational):")
for i in range(4):
    print(f"  n{i} = " + " + ".join(f"({V_inv_frac[i][j]})*{'[10,S1,S2,S3]'[j*3+1:j*3+3] if j > 0 else '10'}" for j in range(4)))

# Now: c5 = f(n0, n1, n2, n3). We know f from the 9 groups.
# The function f can be interpolated as a polynomial in (n1, n2, n3).

# But first, let me try the SIMPLEST approach:
# express c5 directly as a rational function of S1, S2, S3.

# From the 9 groups:
# (S1, S2, S3, c5) =
# (0,  0,  0,  0)
# (3,  3,  3,  0)
# (6,  8,  12, 0)
# (9,  13, 21, 1)
# (9,  15, 33, 1)
# (12, 18, 30, 2)
# (12, 18, 36, 3)
# (12, 20, 36, 1)
# (15, 25, 45, 2)

# Fit c5 = p(S1, S2, S3) / q where p is low degree polynomial.
# We need at least cubic (quadratic didn't work).

# Try: c5 = (a*S1^3 + b*S1^2*S2 + c*S1*S2^2 + d*S2*S3 + e*S1*S3 + f*S3^2/S1 + ...) / denom
# Hmm, rational functions are hard to fit.

# Let me just use the (n1, n2, n3) approach with Lagrange interpolation.
# c5 is determined by (n1, n2, n3). There are 9 points.
# Each nk is an integer >= 0 with n0+n1+n2+n3 = 10.

# Actually the simplest way: since n3 in {0, 1} in all cases,
# split into n3=0 and n3=1 sub-formulas.

print(f"\n{'='*60}")
print("SPLIT BY n3 VALUE")
print(f"{'='*60}")

print("n3=0 cases:")
for n0, n1, n2, n3, c5 in groups:
    if n3 == 0:
        print(f"  n1={n1}, n2={n2}: c5={c5}")

print("\nn3=1 cases:")
for n0, n1, n2, n3, c5 in groups:
    if n3 == 1:
        print(f"  n1={n1}, n2={n2}: c5={c5}")

# n3=0: (n1, n2, c5) = (0,0,0), (3,0,0), (4,1,0), (5,2,1), (6,3,2), (4,4,1), (5,5,2)
# Pattern for n3=0:
# c5 = 0 when n2 <= 1
# c5 = n2-1 when n2 >= 2 and ... no, (4,4,1): n2=4 but c5=1, not 3.
# Hmm.

# Let me compute c3 = S1/3 = (n1+2*n2)/3 for n3=0:
# (0,0): c3=0, c5=0
# (3,0): c3=1, c5=0
# (4,1): c3=2, c5=0
# (5,2): c3=3, c5=1
# (6,3): c3=4, c5=2
# (4,4): c3=4, c5=1
# (5,5): c3=5, c5=2

# For n3=0, c3 < 3: c5=0
# For n3=0, c3=3: c5=1
# For n3=0, c3=4: c5 = 2 if n2=3 (star), 1 if n2=4 (distributed)
# For n3=0, c3=5: c5=2

# For c3=4, n3=0: c5 = 6-n2? Check: n2=3: c5=2=6-4? No, 6-3=3 not 2.
# c5 = n2-1? n2=3: 2 yes, n2=4: 3 no.
# c5 = (n2 choose 2) / something?

# n3=1: (n1, n2, c5) = (6,0,1), (9,0,3)
# With n3=1: n0 = 10-n1-n2-1 = 9-n1-n2
# (6,0): n0=3, c3=3, c5=1
# (9,0): n0=0, c3=4, c5=3

# c5 = n3*(c3-2)? c3=3: 1*1=1 yes. c3=4: 1*2=2 no (c5=3).
# c5 = 3*n3 when c3=4? (9,0,1): c5=3=3*1 yes. But (6,0,1): c5=1=1*1 yes.
# So c5 = n3*c3 - 2*n3? c3=3: 3-2=1 yes. c3=4: 4-2=2 no.

# Let me just try a different polynomial basis
print(f"\n{'='*60}")
print("POLYNOMIAL IN (c3, n2, n3)")
print(f"{'='*60}")

# c3 = (n1 + 2*n2 + 3*n3) / 3. Since S1 = 3*c3, c3 is always integer.
basis_names = ['1', 'c3', 'n2', 'n3', 'c3^2', 'c3*n2', 'c3*n3', 'n2^2', 'n2*n3', 'n3^2']
def basis_eval(c3, n2, n3):
    return [1, c3, n2, n3, c3**2, c3*n2, c3*n3, n2**2, n2*n3, n3**2]

X_b = np.array([basis_eval((n1+2*n2+3*n3)//3, n2, n3) for _, n1, n2, n3, _ in groups], dtype=float)
y_b = np.array([c5 for _, _, _, _, c5 in groups], dtype=float)

# Try subsets
for size in range(3, 10):
    found = False
    for combo in combinations(range(len(basis_names)), size):
        X_sub = X_b[:, list(combo)]
        if matrix_rank(X_sub) < min(len(combo), 9):
            continue
        c_sub, _, _, _ = lstsq(X_sub, y_b, rcond=None)
        e_sub = np.max(np.abs(X_sub @ c_sub - y_b))
        if e_sub < 0.001:
            # Find rational
            for denom in range(1, 601):
                int_c = np.round(c_sub * denom)
                err_r = np.max(np.abs(X_sub @ (int_c / denom) - y_b))
                if err_r < 0.001:
                    selected = [basis_names[i] for i in combo]
                    terms = [f"{int(int_c[j])}*{selected[j]}" for j in range(len(selected)) if int(int_c[j]) != 0]
                    print(f"  {size} features: c5 = ({' + '.join(terms)}) / {denom}")
                    found = True
                    break
            if found:
                break
    if found:
        break

# Verify the formula on all 1024 tournaments
print(f"\n{'='*60}")
print("VERIFICATION ON ALL 1024 TOURNAMENTS")
print(f"{'='*60}")

# Let me extract the best formula found above and verify
# Use the 9-group lookup as ground truth
group_lookup = {}
for n0, n1, n2, n3, c5 in groups:
    group_lookup[(n0, n1, n2, n3)] = c5

verified = 0
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A5 = np.linalg.matrix_power(A, 5)
    c5_true = int(np.trace(A5)) // 5

    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    hist = Counter(lam)
    nk = tuple(hist.get(k, 0) for k in range(4))

    c5_lookup = group_lookup.get(nk, -1)
    if c5_lookup == c5_true:
        verified += 1
    else:
        print(f"  FAIL: bits={bits}, hist={nk}, c5_true={c5_true}, lookup={c5_lookup}")

print(f"Lookup table verified: {verified}/1024")

# Now test: does the PER-VERTEX-SET formula work at n=6?
# c5_total(T, n=6) = sum_{V in C(6,5)} c5(T[V])
# where c5(T[V]) uses the RESTRICTED lambda histogram on V.
print(f"\n{'='*60}")
print("PER-VERTEX-SET FORMULA AT n=6")
print(f"{'='*60}")

n6 = 6
total_bits_6 = n6 * (n6-1) // 2

def count_hc_on_subset(A, vset, k):
    """Count directed k-cycles on vertex set vset."""
    combo = tuple(sorted(vset))
    assert len(combo) == k
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

np.random.seed(42)
success = 0
total_tests = 500
for trial in range(total_tests):
    bits = np.random.randint(0, 1 << total_bits_6)
    A = bits_to_adj(bits, n6)

    # Direct count
    c5_total_direct = 0
    c5_total_formula = 0
    for combo in combinations(range(n6), 5):
        # Direct
        hc = count_hc_on_subset(A, combo, 5)
        c5_total_direct += hc

        # Formula via restricted lambda histogram
        combo_list = list(combo)
        A_sub = A[np.ix_(combo_list, combo_list)]
        L_sub = lambda_graph(A_sub, 5)
        lam_sub = [L_sub[i][j] for i in range(5) for j in range(i+1, 5)]
        hist = Counter(lam_sub)
        nk = tuple(hist.get(k, 0) for k in range(4))
        c5_formula = group_lookup.get(nk, -1)
        if c5_formula == -1:
            print(f"  Unknown histogram: {nk}")
        c5_total_formula += c5_formula

    if c5_total_direct == c5_total_formula:
        success += 1
    elif trial < 3:
        print(f"  Trial {trial}: direct={c5_total_direct}, formula={c5_total_formula}")

print(f"Per-vertex-set formula at n=6: {success}/{total_tests}")

# And at n=7
print(f"\n--- n=7 ---")
n7 = 7
total_bits_7 = n7 * (n7-1) // 2

np.random.seed(42)
success7 = 0
total_tests7 = 200
for trial in range(total_tests7):
    bits = np.random.randint(0, 1 << total_bits_7)
    A = bits_to_adj(bits, n7)

    c5_total_direct = 0
    c5_total_formula = 0
    for combo in combinations(range(n7), 5):
        hc = count_hc_on_subset(A, combo, 5)
        c5_total_direct += hc

        combo_list = list(combo)
        A_sub = A[np.ix_(combo_list, combo_list)]
        L_sub = lambda_graph(A_sub, 5)
        lam_sub = [L_sub[i][j] for i in range(5) for j in range(i+1, 5)]
        hist = Counter(lam_sub)
        nk = tuple(hist.get(k, 0) for k in range(4))
        c5_formula = group_lookup.get(nk, -1)
        if c5_formula == -1:
            print(f"  Unknown histogram at n=7: {nk}")
            c5_total_formula = -999
            break
        c5_total_formula += c5_formula

    if c5_total_direct == c5_total_formula:
        success7 += 1

print(f"Per-vertex-set formula at n=7: {success7}/{total_tests7}")

print("\nDone.")
