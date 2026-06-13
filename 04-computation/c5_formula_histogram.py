"""
c5_formula_histogram.py -- kind-pasteur-2026-03-13-S61
Express c5_dir in terms of lambda histogram (n0, n1, n2, n3) at n=5.

From c5_formula_v2.py:
- c5 is determined by the lambda MULTISET (9 groups)
- The multiset is captured by (n0, n1, n2, n3) where nk = #{pairs with lambda=k}
- Constraint: n0 + n1 + n2 + n3 = C(5,2) = 10
- S_k = sum n_j * j^k

So c5 = f(n0, n1, n2, n3) where n0+n1+n2+n3=10.

Also try: express tr(A^5) directly via lambda and A entries.

Key relationship discovered:
lambda_V(u,v) = lambda(u,v) - witness(k,u,v)
where witness(k,u,v) = 1 if (u,v,k) or (v,u,k) is a directed 3-cycle.
"""

import numpy as np
from itertools import combinations, permutations
from numpy.linalg import lstsq, matrix_rank
from collections import defaultdict, Counter

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

print("=" * 60)
print(f"c5_dir FROM LAMBDA HISTOGRAM (n={n})")
print("=" * 60)

# Build the 9-group table: histogram -> c5
hist_c5 = defaultdict(set)
hist_count = Counter()
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A5 = np.linalg.matrix_power(A, 5)
    c5 = int(np.trace(A5)) // 5

    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    hist = tuple(Counter(lam).get(k, 0) for k in range(n-1))  # (n0, n1, n2, n3)

    hist_c5[hist].add(c5)
    hist_count[hist] += 1

print(f"\n{'n0':>3} {'n1':>3} {'n2':>3} {'n3':>3} | {'c5':>3} | {'count':>5} | S1  S2  S3")
print("-" * 55)
table = []
for hist in sorted(hist_c5.keys()):
    n0, n1, n2, n3 = hist
    c5_vals = hist_c5[hist]
    assert len(c5_vals) == 1
    c5 = list(c5_vals)[0]
    S1 = 0*n0 + 1*n1 + 2*n2 + 3*n3
    S2 = 0*n0 + 1*n1 + 4*n2 + 9*n3
    S3 = 0*n0 + 1*n1 + 8*n2 + 27*n3
    table.append((n0, n1, n2, n3, c5, hist_count[hist], S1, S2, S3))
    print(f"{n0:3d} {n1:3d} {n2:3d} {n3:3d} | {c5:3d} | {hist_count[hist]:5d} | {S1:3d} {S2:3d} {S3:3d}")

# With constraint n0+n1+n2+n3=10, use n1, n2, n3 as free variables
# (n0 = 10 - n1 - n2 - n3)
print(f"\n--- Fitting c5 = polynomial in (n1, n2, n3) ---")
features_and_names = [
    ('1', lambda n1, n2, n3: 1),
    ('n1', lambda n1, n2, n3: n1),
    ('n2', lambda n1, n2, n3: n2),
    ('n3', lambda n1, n2, n3: n3),
    ('n1^2', lambda n1, n2, n3: n1**2),
    ('n1*n2', lambda n1, n2, n3: n1*n2),
    ('n1*n3', lambda n1, n2, n3: n1*n3),
    ('n2^2', lambda n1, n2, n3: n2**2),
    ('n2*n3', lambda n1, n2, n3: n2*n3),
    ('n3^2', lambda n1, n2, n3: n3**2),
]

# Full quadratic
X = np.array([[f(n1, n2, n3) for _, f in features_and_names]
              for n0, n1, n2, n3, c5, _, _, _, _ in table], dtype=float)
y = np.array([c5 for _, _, _, _, c5, _, _, _, _ in table], dtype=float)

rank = matrix_rank(X)
print(f"Design matrix: {X.shape}, rank: {rank}")

coeffs, _, _, _ = lstsq(X, y, rcond=None)
pred = X @ coeffs
err = np.max(np.abs(pred - y))
print(f"Max error: {err:.8f}")

if err < 0.001:
    feat_names = [name for name, _ in features_and_names]
    for denom in range(1, 2001):
        int_c = np.round(coeffs * denom)
        err_r = np.max(np.abs(X @ (int_c / denom) - y))
        if err_r < 0.001:
            print(f"\nEXACT FORMULA (denom={denom}):")
            terms = []
            for i, name in enumerate(feat_names):
                c = int(int_c[i])
                if c != 0:
                    terms.append(f"({c:+d})*{name}")
            print(f"  c5 = ({' '.join(terms)}) / {denom}")

            # Verify on all 1024 tournaments
            success = 0
            for bits in range(1 << total_bits):
                A = bits_to_adj(bits, n)
                L = lambda_graph(A, n)
                A5 = np.linalg.matrix_power(A, 5)
                c5_true = int(np.trace(A5)) // 5
                lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
                hist = Counter(lam)
                n1v, n2v, n3v = hist.get(1, 0), hist.get(2, 0), hist.get(3, 0)
                c5_pred = sum(int_c[i] * f(n1v, n2v, n3v) for i, (_, f) in enumerate(features_and_names)) / denom
                if abs(c5_pred - c5_true) < 0.001:
                    success += 1
            print(f"  Verified: {success}/1024")
            break
else:
    print("Quadratic in (n1,n2,n3) doesn't give exact fit.")

    # Try cubic terms
    print("\n--- Adding cubic terms ---")
    cubic_features = features_and_names + [
        ('n1^3', lambda n1, n2, n3: n1**3),
        ('n1^2*n2', lambda n1, n2, n3: n1**2 * n2),
        ('n1^2*n3', lambda n1, n2, n3: n1**2 * n3),
        ('n1*n2^2', lambda n1, n2, n3: n1 * n2**2),
        ('n1*n2*n3', lambda n1, n2, n3: n1 * n2 * n3),
        ('n1*n3^2', lambda n1, n2, n3: n1 * n3**2),
        ('n2^3', lambda n1, n2, n3: n2**3),
        ('n2^2*n3', lambda n1, n2, n3: n2**2 * n3),
        ('n2*n3^2', lambda n1, n2, n3: n2 * n3**2),
        ('n3^3', lambda n1, n2, n3: n3**3),
    ]
    X3 = np.array([[f(n1, n2, n3) for _, f in cubic_features]
                   for n0, n1, n2, n3, c5, _, _, _, _ in table], dtype=float)
    rank3 = matrix_rank(X3)
    print(f"Design matrix: {X3.shape}, rank: {rank3}")
    coeffs3, _, _, _ = lstsq(X3, y, rcond=None)
    err3 = np.max(np.abs(X3 @ coeffs3 - y))
    print(f"Max error: {err3:.8f}")

# Try the SIMPLEST approach: direct lookup
# Since there are only 9 groups, just use indicator functions
print(f"\n{'='*60}")
print("DIRECT ALGEBRAIC APPROACH: tr(A^5) via lambda")
print(f"{'='*60}")

# tr(A^5) = sum_{i,j,k,l,m} A[i,j]*A[j,k]*A[k,l]*A[l,m]*A[m,i]
# where all indices distinct (for tournament) since non-Hamiltonian closed 5-walks
# cannot exist.
#
# So tr(A^5) = 5 * c5_dir = 5 * sum_{5-vertex-sets V} hc(T[V])
# where hc = number of directed Hamiltonian cycles in T[V].
#
# At n=5, V = {0,1,2,3,4} always.
# hc(T[V]) counts the directed 5-cycles.
#
# For a tournament on 5 vertices, hc depends on the score sequence:
# Transitive (0,1,2,3,4): hc=0
# Near-transitive (0,1,2,3,4) variants: hc=0-1
# Regular (2,2,2,2,2): hc=2 (always!)
#
# Actually, from the table:
# S1=15 (c3=5, regular): c5=2. So regular T_5 has 2 directed 5-cycles.
# This makes sense: T_5 has exactly 12 directed Hamiltonian paths,
# and 2*5=10 directed Hamiltonian cycle arcs... actually let me verify.

print("Score sequence -> c5 mapping:")
score_c5 = defaultdict(set)
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    scores = tuple(sorted(int(sum(A[i])) for i in range(n)))
    A5 = np.linalg.matrix_power(A, 5)
    c5 = int(np.trace(A5)) // 5
    score_c5[scores].add(c5)

for score, vals in sorted(score_c5.items()):
    print(f"  {score}: c5 in {sorted(vals)}")

# Now try to express c5 in terms of c3 and some other lambda invariant
print(f"\n--- c5 as function of c3 and S2 ---")
c3_s2_c5 = defaultdict(set)
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    S1 = sum(lam)
    S2 = sum(l**2 for l in lam)
    c3 = S1 // 3
    A5 = np.linalg.matrix_power(A, 5)
    c5 = int(np.trace(A5)) // 5
    c3_s2_c5[(c3, S2)].add(c5)

print(f"{'c3':>3} {'S2':>4} | c5")
for key in sorted(c3_s2_c5.keys()):
    vals = c3_s2_c5[key]
    marker = " <-- AMBIGUOUS" if len(vals) > 1 else ""
    print(f"  {key[0]:3d} {key[1]:4d} | {sorted(vals)}{marker}")

# c3=4, S2=18 splits into c5=2 and c5=3. What distinguishes them?
# S3 does (30 vs 36). What does S3 mean combinatorially?
# S3 = sum lambda(u,v)^3. At lambda in {0,1,2,3}: S3 counts "weighted triple-lambda"
# The difference: (0,1,1,1,1,1,1,2,2,2) has S3=30, (1,1,1,1,1,1,1,1,1,3) has S3=36.
# First: one pair not in any c3, three pairs in 2 c3s each.
# Second: one pair in ALL 3 c3s (lambda=3), rest in exactly 1.

print(f"\n--- The critical case: c3=4, S2=18 ---")
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    S1 = sum(lam)
    S2 = sum(l**2 for l in lam)
    if S1 == 12 and S2 == 18:
        S3 = sum(l**3 for l in lam)
        A5 = np.linalg.matrix_power(A, 5)
        c5 = int(np.trace(A5)) // 5
        scores = tuple(sorted(int(sum(A[i])) for i in range(n)))
        ms = tuple(sorted(lam))
        print(f"  bits={bits}: scores={scores}, lam_multiset={ms}, S3={S3}, c5={c5}")
        # Print the 3-cycles
        c3_sets = []
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if A[a][b]*A[b][c]*A[c][a] + A[c][b]*A[b][a]*A[a][c] > 0:
                        c3_sets.append((a, b, c))
        print(f"    c3 sets: {c3_sets}")
        if len(c3_sets) == 4:
            # Check: is there a vertex in all 4?
            from functools import reduce
            all_verts = [set(s) for s in c3_sets]
            common = reduce(lambda x, y: x & y, all_verts)
            print(f"    common vertex: {common}")
        break  # Just show one

print(f"\n--- The c5=3 case (one pair in all 3-cycles) ---")
found = False
for bits in range(1 << total_bits):
    if found:
        break
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    lam = [L[i][j] for i in range(n) for j in range(i+1, n)]
    S1 = sum(lam)
    S2 = sum(l**2 for l in lam)
    S3 = sum(l**3 for l in lam)
    if S1 == 12 and S2 == 18 and S3 == 36:
        A5 = np.linalg.matrix_power(A, 5)
        c5 = int(np.trace(A5)) // 5
        scores = tuple(sorted(int(sum(A[i])) for i in range(n)))
        ms = tuple(sorted(lam))
        print(f"  bits={bits}: scores={scores}, lam_multiset={ms}, c5={c5}")
        print(f"  Adjacency:")
        for i in range(n):
            print(f"    {[int(A[i][j]) for j in range(n)]}")
        c3_sets = []
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if A[a][b]*A[b][c]*A[c][a] + A[c][b]*A[b][a]*A[a][c] > 0:
                        c3_sets.append((a, b, c))
        print(f"  c3 sets ({len(c3_sets)}): {c3_sets}")
        # Which pair has lambda=3?
        for i in range(n):
            for j in range(i+1, n):
                if L[i][j] == 3:
                    print(f"  lambda({i},{j})=3: this pair is in ALL {len(c3_sets)} three-cycles")
        # Show the 5-cycles
        for perm in permutations(range(1, n)):
            path = (0,) + perm
            valid = all(A[path[k]][path[(k+1) % n]] == 1 for k in range(n))
            if valid:
                print(f"  5-cycle: {path}")
        found = True

print("\nDone.")
