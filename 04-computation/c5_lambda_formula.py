"""
c5_lambda_formula.py -- kind-pasteur-2026-03-13-S61
Find an explicit formula for c5_dir in terms of lambda(u,v).

Approach: At n=5, enumerate all tournaments and fit c5_dir as a
polynomial in the lambda values.

For n=5: lambda(u,v) for C(5,2)=10 pairs gives a 10-dim feature vector.
c5_dir is the target.

Since c5_dir on 5 vertices counts directed Hamiltonian cycles on the
FULL tournament, there's only ONE 5-vertex set.

For general n: c5_dir = sum over all C(n,5) vertex sets of hc(T[V]).
Each term depends on the 10 arcs within V.

But lambda involves THIRD-party vertices. Lambda(u,v) = #{w: {u,v,w} cyclic}.
For a 5-vertex set V = {a,b,c,d,e}:
  lambda_V(a,b) = #{w in V: {a,b,w} cyclic} (restricted lambda)

The full lambda(a,b) includes vertices outside V too.

So c5_dir on the FULL n-vertex tournament sums over V, and each term
depends on the induced tournament T[V]. The induced tournament on V
has its OWN lambda, which is a restriction of the full lambda.

But the full lambda captures MORE info (external witnesses for 3-cycles).

Let me try: at n=5, find c5_dir as polynomial in lambda.
Then at n=6, check if the same formula works per-5-vertex-set.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

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

# n=5: find c5_dir in terms of lambda
n = 5
total_bits = n * (n-1) // 2
print("=" * 60)
print(f"FINDING c5_dir FORMULA AT n={n}")
print("=" * 60)

# Collect data
data_X = []
data_y = []

for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    # Lambda features (upper triangle)
    lam_features = []
    for i in range(n):
        for j in range(i+1, n):
            lam_features.append(L[i][j])

    c5 = count_directed_hc(A, range(n))
    data_X.append(lam_features)
    data_y.append(c5)

X = np.array(data_X)
y = np.array(data_y)

print(f"Data: {len(data_X)} tournaments, {len(data_X[0])} features")
print(f"c5 values: {sorted(set(data_y))}")
print(f"Lambda value range: {X.min()} to {X.max()}")

# At n=5, lambda(u,v) counts # of 3rd vertices w making a 3-cycle.
# There are n-2 = 3 possible w values, so lambda in {0,1,2,3}.

# Try linear regression: c5 = a0 + sum a_i * lambda_i
from numpy.linalg import lstsq

# Add constant term
X_aug = np.column_stack([np.ones(len(X)), X])
coeffs, residuals, rank, sv = lstsq(X_aug, y, rcond=None)

print(f"\nLinear fit: c5 = {coeffs[0]:.4f} + sum a_i * lambda_i")
print(f"  Coefficients: {[f'{c:.4f}' for c in coeffs]}")
print(f"  Residuals: {residuals if len(residuals) > 0 else 'perfect fit'}")
print(f"  Max abs error: {np.max(np.abs(X_aug @ coeffs - y)):.6f}")

# Are all lambda values the same for each pair? (By symmetry of the labeling)
# At n=5, all edges are equivalent by S_5 symmetry? No — the specific
# tournament breaks the symmetry.

# Try simpler aggregates
# S1 = sum of all lambda values
# S2 = sum of lambda^2
# S3 = sum of lambda^3
S1 = X.sum(axis=1)
S2 = (X**2).sum(axis=1)
S3 = (X**3).sum(axis=1)

# Products: sum lambda(u,v) * lambda(u,w) for distinct pairs
# This is (sum_u (sum_{v>u} lambda(u,v))^2 - sum_{v>u} lambda(u,v)^2) / 2
# Actually let me compute vertex-centered sums
def vertex_lambda_sums(lam_features, n):
    # Reconstruct lambda matrix
    L = np.zeros((n, n))
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            L[i][j] = lam_features[idx]
            L[j][i] = lam_features[idx]
            idx += 1
    vertex_sums = L.sum(axis=1)
    vertex_sq_sums = (L**2).sum(axis=1)
    return vertex_sums, vertex_sq_sums

# Try: c5 = a + b * S1 + c * S2
X_agg = np.column_stack([np.ones(len(X)), S1, S2])
coeffs_agg, res_agg, _, _ = lstsq(X_agg, y, rcond=None)
err_agg = np.max(np.abs(X_agg @ coeffs_agg - y))
print(f"\nAggregate fit: c5 = {coeffs_agg[0]:.4f} + {coeffs_agg[1]:.4f}*S1 + {coeffs_agg[2]:.4f}*S2")
print(f"  Max abs error: {err_agg:.6f}")

# Add S3
X_agg3 = np.column_stack([np.ones(len(X)), S1, S2, S3])
coeffs_agg3, _, _, _ = lstsq(X_agg3, y, rcond=None)
err_agg3 = np.max(np.abs(X_agg3 @ coeffs_agg3 - y))
print(f"\nWith S3: c5 = {coeffs_agg3[0]:.4f} + {coeffs_agg3[1]:.4f}*S1 + {coeffs_agg3[2]:.4f}*S2 + {coeffs_agg3[3]:.4f}*S3")
print(f"  Max abs error: {err_agg3:.6f}")

# Try score sequence
scores = np.array([np.sort([sum(bits_to_adj(bits, n)[i]) for i in range(n)]) for bits in range(1 << total_bits)])
Ss1 = scores.sum(axis=1)
Ss2 = (scores**2).sum(axis=1)
X_score = np.column_stack([np.ones(len(scores)), Ss1, Ss2])
coeffs_score, _, _, _ = lstsq(X_score, y, rcond=None)
err_score = np.max(np.abs(X_score @ coeffs_score - y))
print(f"\nScore fit: c5 = {coeffs_score[0]:.4f} + {coeffs_score[1]:.4f}*Ss1 + {coeffs_score[2]:.4f}*Ss2")
print(f"  Max abs error: {err_score:.6f}")
print(f"  Score determines c5? {err_score < 0.001}")

# At n=5, the score determines the iso class (mostly), and iso
# determines c5. So score should determine c5.
# But lambda is MORE than score. Does lambda EXACTLY determine c5?
# Check: are there tournaments with same lambda but different c5?
print(f"\n  Lambda ambiguity check: {sum(1 for v in Counter(tuple(x) for x in X).values() if v > 1)} repeated lambda vectors")

# Group by lambda and check c5
lam_c5 = {}
for i in range(len(X)):
    key = tuple(X[i])
    if key not in lam_c5:
        lam_c5[key] = set()
    lam_c5[key].add(y[i])

amb = sum(1 for v in lam_c5.values() if len(v) > 1)
print(f"  Lambda groups: {len(lam_c5)}, ambiguous: {amb}")

# Since lambda determines c5 at n=5, let me find the EXACT formula.
# The issue is that the lambda values are correlated with each other.
# At n=5, we can look at it differently:

# For a 5-vertex tournament, there are exactly 2 non-isomorphic types
# with score (2,2,2,2,2) (regular): the Paley and its complement.
# Wait, at n=5 the unique regular tournament has c5 = 2.
# The other score sequences give c5 = 0 or 1.

# Let me check c5 values by score sequence
score_c5 = {}
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    sc = tuple(sorted([sum(A[i]) for i in range(n)]))
    c5 = count_directed_hc(A, range(n))
    if sc not in score_c5:
        score_c5[sc] = set()
    score_c5[sc].add(c5)

print(f"\n  Score -> c5 values:")
for sc in sorted(score_c5.keys()):
    print(f"    {sc}: {sorted(score_c5[sc])}")

# At n=5, how many c3's?
# c3 vertex sets are 3-vertex subsets that support a 3-cycle.
# lambda(u,v) = # of 3-cycles through {u,v}. At n=5, this = 0,1,2, or 3.
# c3 = (1/3) * sum delta(v) where delta(v) = #{3-cycle sets containing v}

# For the regular tournament: every 3-vertex subset is a 3-cycle iff
# the tournament is "strongly regular". Check:
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    sc = tuple(sorted([sum(A[i]) for i in range(n)]))
    if sc == (2,2,2,2,2):
        c3_count = 0
        for combo in combinations(range(n), 3):
            a, b, c = combo
            if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                c3_count += 1
        c5_v = count_directed_hc(A, range(n))
        L = lambda_graph(A, n)
        lam_upper = [L[i][j] for i in range(n) for j in range(i+1,n)]
        print(f"  Regular: c3_sets={c3_count}, c5={c5_v}, lambda={lam_upper}")
        break

# For a GENERAL n, we need c5 as a function of lambda over the FULL graph.
# c5_dir(T) = sum over all 5-vertex V of hc(T[V])
# Each hc(T[V]) depends on the 10 arcs within V.
# But lambda captures info about 3rd parties.

# Key question: is there a known identity expressing c5 in terms of
# c3 and score sequence?

# Moon's formula: c5 = C(n,5) - C(n,3)*c3 + ... ?
# Actually, for tournaments, the number of transitive subtournaments
# of order k is related to score sequence. The number of 5-vertex sets
# that support 0, 1, 2, or 3 Hamiltonian cycles can be expressed in
# terms of the adjacency matrix.

# Let me try expressing c5_dir directly. For a tournament T with
# adjacency matrix A (A[i][j] = 1 if i->j):
# c5_dir = (1/5) * sum over all directed cycles of length 5
#         = (1/5) * sum_{perm sigma, cycle type (5)} prod A[i][sigma(i)]

# This is related to tr(A^5) minus corrections for shorter cycles.
# tr(A^5) counts closed walks of length 5, which includes:
# - 5-cycles (what we want)
# - walks that visit some vertex more than once

# For tournaments: A[i][j] + A[j][i] = 1 for i != j, A[i][i] = 0.
# This constrains A^k in specific ways.

# Let me compute tr(A^5) and c5_dir to see if there's a relationship.
print(f"\n{'='*60}")
print(f"tr(A^5) vs c5_dir (n={n})")
print(f"{'='*60}")

tr_c5_data = {}
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    trA5 = int(np.trace(np.linalg.matrix_power(A, 5)))
    c5 = count_directed_hc(A, range(n))
    if trA5 not in tr_c5_data:
        tr_c5_data[trA5] = set()
    tr_c5_data[trA5].add(c5)

print(f"  tr(A^5) values: {sorted(tr_c5_data.keys())}")
for tr5 in sorted(tr_c5_data.keys()):
    if len(tr_c5_data[tr5]) > 1:
        print(f"  tr(A^5)={tr5}: c5={sorted(tr_c5_data[tr5])} (AMBIGUOUS)")

# At n=5, tr(A^5) = 5*c5_dir + corrections from non-Hamiltonian walks
# Non-Ham 5-walks: visit some vertex twice in 5 steps.
# Types: (v,w,x,v,y,v), etc.

# Actually for n=5, ALL vertices are visited in a 5-walk iff it's Hamiltonian.
# A closed 5-walk on 5 vertices that visits all vertices is a Hamiltonian cycle.
# A closed 5-walk that visits some vertex twice must miss some vertex.
# So: tr(A^5) at n=5 = 5 * c5_dir + (walks missing at least 1 vertex)

# Walks missing vertex v: closed 5-walks on {0,...,4}\{v} = 4-vertex walks.
# On 4 vertices, a closed walk of length 5 visits at least one vertex twice.
# 4-vertex closed 5-walks are enumerable.

# tr(A^5) = 5*c5 + sum_v tr(A_v^5) where A_v = adj matrix of T without vertex v
# This is a recursive formula!

print(f"\n  Verifying: tr(A^5) = 5*c5 + sum_v tr(A_v^5)")
for bits in [0, 1, 42, 100, 500, 1000]:
    if bits >= (1 << total_bits):
        continue
    A = bits_to_adj(bits, n)
    trA5 = int(np.trace(np.linalg.matrix_power(A, 5)))
    c5 = count_directed_hc(A, range(n))

    sub_traces = 0
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = A[np.ix_(remaining, remaining)]
        sub_traces += int(np.trace(np.linalg.matrix_power(A_sub, 5)))

    print(f"  bits={bits}: tr(A^5)={trA5}, 5*c5={5*c5}, sum tr(A_v^5)={sub_traces}, diff={trA5 - 5*c5 - sub_traces}")

# Hmm, the "missing vertex" decomposition may include longer corrections.
# Let me just check: tr(A^5) for 5-vertex tournament, can a closed 5-walk
# visit fewer than 5 vertices?

# Closed walk of length 5 on 5 vertices starting at v0:
# v0 -> v1 -> v2 -> v3 -> v4 -> v0
# This visits vertices {v0, v1, v2, v3, v4} which could be < 5 distinct vertices.

# For a tournament (complete directed graph with one direction per edge),
# EVERY 5-step walk exists (from any vertex there are exactly n-1 outgoing arcs).
# But the walk v0,v1,v2,v3,v4,v0 requires A[v0][v1]*A[v1][v2]*...*A[v4][v0] = 1.

# A walk visiting exactly 4 distinct vertices: one vertex appears twice.
# Possible patterns: v,w,x,y,w,v (vertex w appears at positions 1 and 4)
# Or v,w,v,x,y,v (vertex v appears 3 times)
# Etc.

# This is getting complex. Let me instead just verify the lambda formula.

# At n=5, c5_dir is determined by the lambda graph.
# The score sequence is sum_j A[i][j] for each i.
# lambda(u,v) is related to the "triangle count" through the edge (u,v).
# For a tournament: lambda(u,v) = |{w: A[u][v]*A[v][w]*A[w][u] + A[v][u]*A[u][w]*A[w][v]}|
# This can be expressed as:
# lambda(u,v) = |{w != u,v: exactly one of (u->v->w->u, v->u->w->v) holds}|
# = |{w: (u,v,w) form a directed 3-cycle}|

# At n=5, sum_{u<v} lambda(u,v) = 3 * c3_vertex_set_count

# FORMULA ATTEMPT: c5 as degree-2 polynomial in lambda
# c5 = a0 + a1*S1 + a2*S2 + a3*sum_{uv,uw: v!=w} lambda(u,v)*lambda(u,w)

# Build the feature matrix for degree-2 monomials in lambda
print(f"\n{'='*60}")
print(f"DEGREE-2 POLYNOMIAL FIT")
print(f"{'='*60}")

features = []
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    f = [1]  # constant
    lam = []
    for i in range(n):
        for j in range(i+1, n):
            lam.append(L[i][j])
    # Degree 1
    f.append(sum(lam))  # S1
    # Degree 2
    f.append(sum(l**2 for l in lam))  # S2 = sum lambda^2
    # Vertex-centered product: for each vertex u, sum_{v,w>u, v!=w} lambda(u,v)*lambda(u,w)
    vertex_prod = 0
    for u in range(n):
        neighbors = [L[u][v] for v in range(n) if v != u]
        vertex_prod += sum(neighbors)**2 - sum(l**2 for l in neighbors)
    f.append(vertex_prod // 2)  # sum_{u} C(delta_u, 2) related

    features.append(f)
    data_y_check = count_directed_hc(A, range(n))

X_poly = np.array(features, dtype=float)
y_check = np.array(data_y)
coeffs_poly, _, _, _ = lstsq(X_poly, y_check, rcond=None)
err_poly = np.max(np.abs(X_poly @ coeffs_poly - y_check))
print(f"  Coefficients: {[f'{c:.4f}' for c in coeffs_poly]}")
print(f"  Max abs error: {err_poly:.6f}")
print(f"  Perfect fit? {err_poly < 0.001}")

# Check if we can get integer coefficients
print(f"  Trying integer coefficients:")
for scale in [1, 2, 3, 4, 5, 6, 10, 12, 20, 30, 60]:
    int_coeffs = np.round(coeffs_poly * scale)
    err_int = np.max(np.abs(X_poly @ (int_coeffs / scale) - y_check))
    if err_int < 0.001:
        print(f"    Scale {scale}: {int_coeffs.astype(int)} -> c5 = ({' + '.join(f'{int(c)}*f{i}' for i,c in enumerate(int_coeffs) if c != 0)}) / {scale}")
        print(f"    Error: {err_int:.6f}")
        break

print("\nDone.")
