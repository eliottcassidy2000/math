#!/usr/bin/env python3
"""
Derive c5 = tr(A^5)/5 = tr(A² · A² · A)/5 as a function of lambda and scores.

Key facts:
- (A²)_{ij} = #{k≠i,j: i→k→j} for i≠j, (A²)_{ii} = 0
- λ(i,j) = A_{ij}·(A²)_{ji} + A_{ji}·(A²)_{ij}
- Since A_{ij} + A_{ji} = 1, we have λ(i,j) = A_{ij}·(A²)_{ji} + (1-A_{ij})·(A²)_{ij}

We want: tr(A^5) = Σ_i (A^5)_{ii} = Σ_{i,j} (A²)_{ij} · (A³)_{ji}
       = Σ_{i,j} (A²)_{ij} · Σ_k A_{jk}(A²)_{ki}

Now define P_{ij} = (A²)_{ij} for i≠j, P_{ii} = 0.

Key identity for tournaments (to be verified):
P_{ij} + P_{ji} = (n-2) - λ(i,j)  ... is this right?

No! Let me derive it properly.
P_{ij} = #{k: i→k→j, k≠i,j} (number of directed 2-paths i→k→j through k)
P_{ji} = #{k: j→k→i, k≠i,j}

For each k ≠ i,j, exactly one of 4 cases:
(1) i→k, k→j: contributes to P_{ij}
(2) i→k, j→k: neither (parallel out-arrows from i and j to k) -- wait, j→k means A_{jk}=1
    Actually: k→j means A_{kj}=1. j→k means A_{jk}=1.
    Case (2): A_{ik}=1, A_{jk}=1 (i→k and j→k). Neither P_{ij} nor P_{ji}.
(3) k→i, k→j: A_{ki}=1, A_{kj}=1. Neither P_{ij} nor P_{ji}.
(4) k→i, j→k: A_{ki}=1, A_{jk}=1. Hmm, for P_{ji} = #{k: j→k→i}, need A_{jk}=1, A_{ki}=1. Yes!
    So case (4) contributes to P_{ji}.

Wait, let me redo:
P_{ij} counts k with A_{ik}=1 AND A_{kj}=1.
P_{ji} counts k with A_{jk}=1 AND A_{ki}=1.

Case for k≠i,j: (A_{ik}, A_{ki}) = (1,0) or (0,1); (A_{kj}, A_{jk}) = (1,0) or (0,1).

(A_{ik}, A_{kj}) = (1,1): P_{ij} += 1. Also A_{ki}=0, A_{jk}=0 → P_{ji} += 0.
(A_{ik}, A_{jk}) = (1,1): P_{ij} += 0 (A_{kj}=0). P_{ji} += 0 (A_{ki}=0).
(A_{ki}, A_{kj}) = (1,1): P_{ij} += 0. P_{ji} += 0.
(A_{ki}, A_{jk}) = (1,1): P_{ij} += 0. P_{ji} += 1 (A_{jk}=1, A_{ki}=1).

So exactly: P_{ij} + P_{ji} + (cases 2 and 3) = n-2.
Cases 2: #{k: A_{ik}=1, A_{jk}=1} = #{common OUT-neighbors of i,j} (excluding each other)
         = d_i + d_j - #{k: A_{ik}=1 or A_{jk}=1} ... complicated.

Let's just count:
#{k: A_{ik}=1, A_{jk}=1, k≠i,j} = denote as co(i,j) (common out-successors)
#{k: A_{ki}=1, A_{kj}=1, k≠i,j} = denote as ci(i,j) (common in-predecessors)

Then: P_{ij} + P_{ji} + co(i,j) + ci(i,j) = n-2.

Now: co(i,j) + ci(i,j) = ?
Out-degree of i among {1,...,n}\{i,j}: d_i - A_{ij}
Out-degree of j among {1,...,n}\{i,j}: d_j - A_{ji}

Hmm, this isn't leading to a clean formula in terms of lambda.

Let me try a different approach: just verify computationally what polynomial
in lambda values gives c5 at small n.

opus-2026-03-13-S71c
"""
import sys, time
import numpy as np
from itertools import combinations
from collections import defaultdict
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

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

# At n=5, try to express c5 in terms of:
# - λ(i,j) for all pairs
# - d_i (scores)
# Since A²_{ij} is related to λ(i,j) and d_i.

n = 5
tb = n*(n-1)//2

print("=" * 60)
print(f"EXPRESSING c5 IN TERMS OF λ AND SCORES AT n={n}")
print("=" * 60)

# Build feature matrix and target
features_data = []
targets = []

for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A5 = np.linalg.matrix_power(A, 5)
    c5 = int(np.trace(A5)) // 5

    # Score sequence (ordered by vertex)
    scores = [int(sum(A[i])) for i in range(n)]

    # Lambda values (upper triangle)
    lam_vals = [L[i][j] for i in range(n) for j in range(i+1, n)]

    # Compute various lambda-based features
    S1 = sum(lam_vals)
    S2 = sum(l*l for l in lam_vals)
    S3 = sum(l*l*l for l in lam_vals)

    # Per-vertex lambda sums: S(v) = Σ_{u≠v} λ(v,u)
    per_vertex_lam = [sum(L[v][u] for u in range(n) if u != v) for v in range(n)]
    # Note: S1 = sum(per_vertex_lam) / 2 (each pair counted twice)

    # Try: Σ_{i<j} λ(i,j) * (λ(i,j) - 1)
    lam_sq_shifted = sum(l*(l-1) for l in lam_vals)

    # Try: Σ_i d_i * pv_i where pv_i = Σ_j λ(i,j)
    cross_d_lam = sum(scores[i] * per_vertex_lam[i] for i in range(n))

    # Key attempt: tr(L²)/2 where L is the lambda matrix
    trL2 = sum(L[i][j]**2 for i in range(n) for j in range(n)) // 2  # = S2
    # Note: tr(L²) = Σ_{i,j} L_{ij}² = 2 * Σ_{i<j} L_{ij}² = 2*S2

    # tr(L³)/? -- might give S3
    L_mat = L.astype(float)
    trL3 = int(np.trace(np.linalg.matrix_power(L_mat, 3)))
    # trL3 = Σ_{i,j,k} L_{ij}L_{jk}L_{ki}

    features_data.append([1, S1, S2, S3, cross_d_lam, trL3])
    targets.append(c5)

X = np.array(features_data, dtype=float)
y = np.array(targets, dtype=float)

from numpy.linalg import lstsq

# Try fitting
coeffs, _, _, _ = lstsq(X, y, rcond=None)
pred = X @ coeffs
err = np.max(np.abs(pred - y))
print(f"\nFit with [1, S1, S2, S3, cross_d_lam, trL3]: max error = {err:.6f}")

if err < 0.001:
    # Find rational form
    for denom in range(1, 10000):
        int_c = np.round(coeffs * denom)
        err_r = np.max(np.abs(X @ (int_c / denom) - y))
        if err_r < 0.001:
            feat_names = ['1', 'S1', 'S2', 'S3', 'cross_d_lam', 'trL3']
            terms = []
            for i, name in enumerate(feat_names):
                c = int(int_c[i])
                if c != 0:
                    terms.append(f"{c}*{name}")
            print(f"\n  EXACT: c5 = ({' + '.join(terms)}) / {denom}")
            break

# Now try with JUST lambda (no scores) — can we do it?
# We know from kind-pasteur that (S1, S2, S3) suffice at n=5.
# But (S1, S2) don't — the S1=12, S2=18 case is ambiguous.
# Question: is there a SIMPLER lambda-only formula?

# Try: c5 = f(S1, trL3) — because trL3 captures 3-cycle structure of the lambda graph
X_simple = np.column_stack([np.ones(len(targets)), [f[1] for f in features_data], [f[5] for f in features_data]])
coeffs_s, _, _, _ = lstsq(X_simple, y, rcond=None)
err_s = np.max(np.abs(X_simple @ coeffs_s - y))
print(f"\nFit with [1, S1, trL3]: max error = {err_s:.6f}")

# Try: c5 = f(Σ C(λ(i,j), 2)) — counting lambda pairs
features_comb = []
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A5 = np.linalg.matrix_power(A, 5)
    c5 = int(np.trace(A5)) // 5
    lam_vals = [L[i][j] for i in range(n) for j in range(i+1, n)]
    S1 = sum(lam_vals)
    CL2 = sum(l*(l-1)//2 for l in lam_vals)  # Σ C(λ,2)
    features_comb.append([1, S1, CL2])

X_comb = np.array(features_comb, dtype=float)
coeffs_c, _, _, _ = lstsq(X_comb, y, rcond=None)
err_c = np.max(np.abs(X_comb @ coeffs_c - y))
print(f"Fit with [1, S1, Σ C(λ,2)]: max error = {err_c:.6f}")

# S1 = 3*c3. So try c5 = f(c3, CL2):
# At n=5: c3 = S1/3.
# CL2 = Σ C(λ,2) = (S2 - S1)/2 (since Σ λ² - Σ λ = Σ λ(λ-1) = 2·Σ C(λ,2))
print(f"\nRelation: CL2 = (S2-S1)/2: {all(features_comb[i][2] == (features_data[i][2] - features_data[i][1])/2 for i in range(len(targets)))}")

# The winning formula (from kind-pasteur's analysis):
# c5 is determined by (S1, S2, S3) but the formula is cubic in the histogram.
# At n=5: groups = 9, so a cubic polynomial in (n0,n1,n2,n3) works.
# Histogram: n_k = #{pairs (i,j): λ(i,j) = k}

# Build histogram features
hist_features = []
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    lam_vals = [L[i][j] for i in range(n) for j in range(i+1, n)]
    hist = [0, 0, 0, 0]  # λ ∈ {0,1,2,3} at n=5
    for l in lam_vals:
        hist[l] += 1
    hist_features.append(hist)

# Since n0+n1+n2+n3 = 10, we have 3 independent vars.
# Also S1 = n1+2*n2+3*n3, S2 = n1+4*n2+9*n3, S3 = n1+8*n2+27*n3
# These are linear combinations of (n1,n2,n3). So (S1,S2,S3) and (n1,n2,n3) span same space.

# Cubic in (n1,n2,n3) has C(5,3)=10 terms (up to degree 3 with 3 vars + constant).
# But we only have 9 groups. So full cubic CANNOT be unique.
# The key is: the extra constraint n0+n1+n2+n3=10 reduces freedom.

# Let me find the exact formula using the 9 data points.
data_points = {}
for i, hf in enumerate(hist_features):
    key = tuple(hf)
    data_points[key] = targets[i]

print(f"\n{'='*60}")
print("EXACT c5 FORMULA (9 data points)")
print(f"{'='*60}")

for key in sorted(data_points.keys()):
    print(f"  n=(n0={key[0]}, n1={key[1]}, n2={key[2]}, n3={key[3]}) -> c5={data_points[key]}")

# Since S1,S2 don't suffice but S1,S2,S3 do, and these are linear in (n1,n2,n3),
# the formula must be linear in the POWER SUMS but the issue is that S3 adds
# independent info beyond S1,S2.
# At n=5, (S1,S2,S3) determine (n1,n2,n3) uniquely (since the Vandermonde
# on {1,2,3} is invertible — this is Newton's identity for 3 values).

# So: (n1,n2,n3) determined by (S1,S2,S3), and c5 is a function of (n1,n2,n3).
# The formula c5 = f(n1,n2,n3) should be low-degree polynomial.

# With 9 data points and 3 variables, let's fit:
pts = sorted(data_points.items())
n1_vals = [p[0][1] for p in pts]
n2_vals = [p[0][2] for p in pts]
n3_vals = [p[0][3] for p in pts]
c5_vals = [p[1] for p in pts]

# Linear: c5 = a + b*n1 + c*n2 + d*n3
X_hist_lin = np.column_stack([np.ones(9), n1_vals, n2_vals, n3_vals])
coeffs_hl, _, _, _ = lstsq(X_hist_lin, c5_vals, rcond=None)
err_hl = np.max(np.abs(X_hist_lin @ coeffs_hl - c5_vals))
print(f"\nLinear in (n1,n2,n3): max error = {err_hl:.6f}")

# Quadratic
feats_q = []
for i in range(9):
    n1, n2, n3 = n1_vals[i], n2_vals[i], n3_vals[i]
    feats_q.append([1, n1, n2, n3, n1**2, n1*n2, n1*n3, n2**2, n2*n3, n3**2])
X_q = np.array(feats_q, dtype=float)
coeffs_q, _, _, _ = lstsq(X_q, c5_vals, rcond=None)
err_q = np.max(np.abs(X_q @ coeffs_q - c5_vals))
print(f"Quadratic in (n1,n2,n3): max error = {err_q:.6f}")

if err_q < 0.001:
    for denom in range(1, 10000):
        int_c = np.round(np.array(coeffs_q) * denom)
        err_r = np.max(np.abs(X_q @ (int_c / denom) - c5_vals))
        if err_r < 0.001:
            feat_names = ['1', 'n1', 'n2', 'n3', 'n1²', 'n1·n2', 'n1·n3', 'n2²', 'n2·n3', 'n3²']
            print(f"\n  EXACT (denom={denom}):")
            terms = []
            for j, name in enumerate(feat_names):
                c = int(int_c[j])
                if c != 0:
                    terms.append(f"{c}*{name}")
            print(f"  c5 = ({' + '.join(terms)}) / {denom}")

            # Verify
            for i in range(9):
                n1, n2, n3 = n1_vals[i], n2_vals[i], n3_vals[i]
                pred = sum(int_c[j] * feats_q[i][j] for j in range(10)) / denom
                print(f"    n=({10-n1_vals[i]-n2_vals[i]-n3_vals[i]},{n1},{n2},{n3}): pred={pred:.1f}, actual={c5_vals[i]}")
            break
    # Also verify on ALL 1024 tournaments
    print(f"\n  Verifying on all 1024 tournaments...")
    verified = True
    for bits_check in range(1 << tb):
        A = bits_to_adj(bits_check, n)
        L = lambda_graph(A, n)
        lam_vals = [L[i][j] for i in range(n) for j in range(i+1, n)]
        hist = [0, 0, 0, 0]
        for l in lam_vals:
            hist[l] += 1
        n1, n2, n3 = hist[1], hist[2], hist[3]
        feats = [1, n1, n2, n3, n1**2, n1*n2, n1*n3, n2**2, n2*n3, n3**2]
        pred = sum(int_c[j] * feats[j] for j in range(10)) / denom
        actual = targets[bits_check]
        if abs(pred - actual) > 0.001:
            verified = False
            print(f"    FAIL: bits={bits_check}, pred={pred:.4f}, actual={actual}")
            break
    print(f"  All 1024 verified: {verified}")

# Now the BIG test: does this formula (or ANY polynomial in histogram) work at n=6?
# At n=6, λ(i,j) can range from 0 to n-2=4, so histogram has bins n0,...,n4.
# And we need to count c5 = #{5-subsets with Hamiltonian directed cycle}.

# The formula at n=5 gives c5 for a SINGLE tournament on 5 vertices.
# At n=6, the total c5 = Σ_{V ∈ C(6,5)} c5(T[V]).
# Each T[V] has its own lambda graph (RESTRICTED to V).

# Can we express Σ_V c5(T[V]) in terms of the FULL lambda graph?
# This requires summing the formula over all 5-subsets, which involves
# restricted lambda histograms.

print(f"\n{'='*60}")
print("EXTENDING TO n=6: SUM OVER 5-SUBSETS")
print(f"{'='*60}")

# At n=6: for each 5-subset V, lambda_V(u,v) = lambda(u,v) - [excluded vertex witnesses (u,v)]
# This was verified by kind-pasteur.

# So: c5(T) = Σ_V f(histogram_V)
# where histogram_V depends on restricted lambda values lambda_V(u,v).

# This IS a polynomial in the lambda values of the full graph,
# because lambda_V(u,v) = lambda(u,v) - indicator(excluded vertex)
# and the indicator is itself determined by the adjacency (which is lambda-related).

# BUT: the indicator function depends on A_{ij}, not just lambda!
# witness(k, u, v) = 1 iff (u,v,k) forms a directed 3-cycle.
# This depends on the DIRECTION of the edges, not just whether
# the 3-cycle exists (which lambda captures).

# Wait: lambda(u,v) counts the NUMBER of directed 3-cycles through {u,v}.
# But for a specific witness k, whether k witnesses {u,v} depends on the
# edge DIRECTIONS between k, u, v — not just on lambda.

# Hmm, but at n=5 the formula only uses the histogram, which is symmetric.
# At n=6, the per-vertex witness indicator might break the symmetry.

# Let me just verify: at n=6, is c5 determined by labeled lambda?
print("Checking c5 vs labeled lambda at n=6 (exhaustive)...")

n6 = 6
tb6 = n6*(n6-1)//2
lam_groups_6 = defaultdict(set)

for bits in range(1 << tb6):
    A = bits_to_adj(bits, n6)
    L = lambda_graph(A, n6)
    labeled_key = tuple(L[i][j] for i in range(n6) for j in range(i+1, n6))

    c5_total = 0
    for combo in combinations(range(n6), 5):
        A_sub = A[np.ix_(list(combo), list(combo))]
        c5_total += int(np.trace(np.linalg.matrix_power(A_sub, 5))) // 5

    lam_groups_6[labeled_key].add(c5_total)

ambig_6 = sum(1 for v in lam_groups_6.values() if len(v) > 1)
print(f"  Lambda groups: {len(lam_groups_6)}, ambiguous: {ambig_6}")

if ambig_6 == 0:
    print("  c5 IS labeled-lambda-determined at n=6!")
else:
    print(f"  c5 NOT labeled-lambda-determined at n=6: {ambig_6} ambiguities")

print(f"\nDone.")
