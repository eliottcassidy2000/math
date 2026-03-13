"""
sigma_c7_formula.py -- kind-pasteur-2026-03-13-S61

THM-183 says c7 = f(tr(Sigma^2), tr(Sigma^3), tr(Sigma^4)).
This script attempts to find the EXACT formula.

Key observations:
- tr(Sigma) = total_sigma = sum(s_i^2) - n(n-1)/2 (determined by scores)
- tr(Sigma^2) = sum_{u,v} sigma(u,v)^2
- tr(Sigma^3) = sum_{u,v,w} sigma(u,v)*sigma(v,w)*sigma(w,u) (weighted triangles)
- tr(Sigma^4) = sum of weighted 4-paths

Strategy: Fit c7 = a*tr1 + b*tr2 + c*tr3 + d*tr4 + e
with rational coefficients.
"""

import numpy as np
from itertools import combinations
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

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print(f"SIGMA-C7 FORMULA SEARCH AT n={n}")
print("=" * 60)

np.random.seed(42)

# Collect data with traces of sigma, lambda, delta, and A^k
data_points = []

for trial in range(10000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A2 = A @ A

    c3 = count_directed_k_cycles(A, n, 3)
    c5 = count_directed_k_cycles(A, n, 5)
    c7 = count_directed_k_cycles(A, n, 7)

    # Build sigma matrix
    sig_mat = np.zeros((n, n), dtype=int)
    lam_mat = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            sig = n - 2 - int(A2[u][v]) - int(A2[v][u])
            lam = int(L[u][v])
            sig_mat[u][v] = sig_mat[v][u] = sig
            lam_mat[u][v] = lam_mat[v][u] = lam

    scores = [int(sum(A[i])) for i in range(n)]

    # Traces
    tr_s1 = int(np.trace(sig_mat))  # = 0 (diagonal is 0)
    total_sigma = int(sig_mat.sum()) // 2  # = sum of upper triangle
    tr_s2 = int(np.trace(sig_mat @ sig_mat))
    tr_s3 = int(np.trace(sig_mat @ sig_mat @ sig_mat))
    tr_s4 = int(np.trace(np.linalg.matrix_power(sig_mat, 4)))

    # Also lambda traces
    tr_l2 = int(np.trace(lam_mat @ lam_mat))
    tr_l3 = int(np.trace(lam_mat @ lam_mat @ lam_mat))

    # Also A^k traces
    tr_A2 = int(np.trace(A2))
    tr_A3 = int(np.trace(A2 @ A)) // 1  # = 3*c3

    # Score-dependent quantities
    ssq = sum(s*s for s in scores)

    data_points.append({
        'c7': c7,
        'c3': c3,
        'c5': c5,
        'total_sigma': total_sigma,
        'tr_s2': tr_s2,
        'tr_s3': tr_s3,
        'tr_s4': tr_s4,
        'tr_l2': tr_l2,
        'tr_l3': tr_l3,
        'ssq': ssq,
    })

# Feature matrix: try different combinations
features_names = ['c3', 'c5', 'tr_s2', 'tr_s3', 'tr_s4', 'tr_l2', 'tr_l3', 'total_sigma', 'ssq']
X_full = np.array([[d[f] for f in features_names] for d in data_points], dtype=float)
y = np.array([d['c7'] for d in data_points], dtype=float)

# Add constant term
X_const = np.column_stack([X_full, np.ones(len(y))])

# 1. Full linear fit
print("\n--- Full Linear Fit ---")
coeffs, residuals, rank, sv = np.linalg.lstsq(X_const, y, rcond=None)
pred = X_const @ coeffs
rmse = np.sqrt(np.mean((y - pred)**2))
r2 = 1 - np.var(y - pred) / np.var(y)
print(f"Features: {features_names} + const")
for i, name in enumerate(features_names + ['const']):
    print(f"  {name}: {coeffs[i]:.6f}")
print(f"  RMSE: {rmse:.4f}, R^2: {r2:.6f}")

# 2. Try with just sigma traces
print("\n--- Sigma Traces Only ---")
sig_features = ['total_sigma', 'tr_s2', 'tr_s3', 'tr_s4']
X_sig = np.array([[d[f] for f in sig_features] for d in data_points], dtype=float)
X_sig_c = np.column_stack([X_sig, np.ones(len(y))])
coeffs_sig, _, _, _ = np.linalg.lstsq(X_sig_c, y, rcond=None)
pred_sig = X_sig_c @ coeffs_sig
rmse_sig = np.sqrt(np.mean((y - pred_sig)**2))
r2_sig = 1 - np.var(y - pred_sig) / np.var(y)
print(f"Features: {sig_features} + const")
for i, name in enumerate(sig_features + ['const']):
    print(f"  {name}: {coeffs_sig[i]:.6f}")
print(f"  RMSE: {rmse_sig:.4f}, R^2: {r2_sig:.6f}")

# 3. Try with (c3, tr_s2, tr_s3, tr_s4)
print("\n--- (c3, tr_s2, tr_s3, tr_s4) ---")
c3_sig_features = ['c3', 'tr_s2', 'tr_s3', 'tr_s4']
X_c3sig = np.array([[d[f] for f in c3_sig_features] for d in data_points], dtype=float)
X_c3sig_c = np.column_stack([X_c3sig, np.ones(len(y))])
coeffs_c3sig, _, _, _ = np.linalg.lstsq(X_c3sig_c, y, rcond=None)
pred_c3sig = X_c3sig_c @ coeffs_c3sig
rmse_c3sig = np.sqrt(np.mean((y - pred_c3sig)**2))
r2_c3sig = 1 - np.var(y - pred_c3sig) / np.var(y)
print(f"Features: {c3_sig_features} + const")
for i, name in enumerate(c3_sig_features + ['const']):
    print(f"  {name}: {coeffs_c3sig[i]:.6f}")
print(f"  RMSE: {rmse_c3sig:.4f}, R^2: {r2_c3sig:.6f}")

# 4. Try adding c3^2, c3*c5, c5^2 etc.
print("\n--- Quadratic Features ---")
quad_data = np.array([
    [d['c3'], d['c5'], d['tr_s2'], d['tr_s3'], d['tr_s4'],
     d['c3']**2, d['c5']**2, d['c3']*d['c5'],
     d['c3']*d['tr_s3'], d['c5']*d['tr_s3'],
     d['tr_s2']*d['tr_s3'], d['tr_s3']**2]
    for d in data_points
], dtype=float)
quad_names = ['c3', 'c5', 'tr_s2', 'tr_s3', 'tr_s4',
              'c3^2', 'c5^2', 'c3*c5',
              'c3*tr_s3', 'c5*tr_s3', 'tr_s2*tr_s3', 'tr_s3^2']
X_quad_c = np.column_stack([quad_data, np.ones(len(y))])
coeffs_quad, _, _, _ = np.linalg.lstsq(X_quad_c, y, rcond=None)
pred_quad = X_quad_c @ coeffs_quad
rmse_quad = np.sqrt(np.mean((y - pred_quad)**2))
r2_quad = 1 - np.var(y - pred_quad) / np.var(y)
print(f"  RMSE: {rmse_quad:.4f}, R^2: {r2_quad:.6f}")

# 5. Try with (c3, c5, c3^2, c5^2, c3*c5)
print("\n--- Just c3, c5, and products ---")
c35_data = np.array([
    [d['c3'], d['c5'], d['c3']**2, d['c5']**2, d['c3']*d['c5']]
    for d in data_points
], dtype=float)
X_c35_c = np.column_stack([c35_data, np.ones(len(y))])
coeffs_c35, _, _, _ = np.linalg.lstsq(X_c35_c, y, rcond=None)
pred_c35 = X_c35_c @ coeffs_c35
rmse_c35 = np.sqrt(np.mean((y - pred_c35)**2))
r2_c35 = 1 - np.var(y - pred_c35) / np.var(y)
c35_names = ['c3', 'c5', 'c3^2', 'c5^2', 'c3*c5', 'const']
for i, name in enumerate(c35_names):
    print(f"  {name}: {coeffs_c35[i]:.6f}")
print(f"  RMSE: {rmse_c35:.4f}, R^2: {r2_c35:.6f}")

# 6. The key: try the EXACT formula c7 = (tr(A^7))/7
#    But tr(A^k) relates to directed walks, not sigma.
#    However, sigma is defined via A^2, so maybe there's a trace identity.
print("\n\n--- Sigma-to-A Relationship ---")
# sigma(u,v) = n-2 - A^2[u,v] - A^2[v,u]
# So Sigma = (n-2)*J - diag_correction - (A^2 + (A^2)^T)
# where J is all-ones, diag_correction handles the diagonal
# Actually Sigma(u,v) = n-2 - (A^2)[u,v] - (A^2)[v,u] for u != v
# and Sigma(u,u) = 0.
# Note: (A^2)[u,v] + (A^2)[v,u] = #{common successors} + #{common predecessors}
# = #{w: u->w, v->w} + #{w: w->u, w->v}
# This is the OVERLAP between the out-neighborhoods of u,v
# plus the overlap of in-neighborhoods.
# By inclusion-exclusion: C(u,v) = (A^2 + (A^T)^2)[u,v]? No...
# (A^2)[u,v] = sum_w A[u,w]*A[w,v] = #{paths u->w->v of length 2}
# This is the number of common out-successors of u that are in-predecessors of v.
# For a tournament: A[u,w] + A[w,u] = 1 (u!=w), so
# (A^2)[u,v] = #{w: u->w AND w->v} = #{w: w is between u and v in some sense}

# Let's verify:
# (A^2 + (A^T)^2)[u,v] = sum_w (A[u,w]*A[w,v] + A[w,u]*A[v,w])
# = #{u->w->v} + #{u<-w<-v}
# = number of common "between" vertices (u->w->v) + common (v->w->u)
# These are NOT the same as common successors/predecessors.
# Common successors = #{w: u->w AND v->w} = sum_w A[u,w]*A[v,w] = (A*A^T)[u,v]
# Common predecessors = #{w: w->u AND w->v} = (A^T*A)[u,v]
# sigma(u,v) = n-2 - (A*A^T)[u,v] - (A^T*A)[u,v]... NO!
# Wait, let me re-derive.
# sigma(u,v) = (common succ) + (common pred) = (A*A^T)[u,v] + (A^T*A)[u,v]
# But wait, for u!=v:
# Common succ of {u,v} = #{w!=u,v: u->w AND v->w} = sum_{w!=u,v} A[u,w]*A[v,w]
# = (A*A^T)[u,v] - A[u,v]*A[v,v]  [but A[v,v]=0]
# Hmm, actually (A*A^T)[u,v] = sum_w A[u,w]*A[v,w] = #{w: u->w AND v->w}
# This includes w=u and w=v? No: A[u,u]=0 and A[v,v]=0, so:
# (A*A^T)[u,v] = sum_{w!=u,v} A[u,w]*A[v,w] + A[u,u]*A[v,u] + A[u,v]*A[v,v]
# = sum_{w!=u,v} A[u,w]*A[v,w] + 0 + 0
# = #{w!=u,v: u->w AND v->w} = common successors.
# Similarly (A^T*A)[u,v] = sum_w A[w,u]*A[w,v] = common predecessors.
# So sigma(u,v) = (AAT + ATA)[u,v] for u!=v.

print("Sigma = A*A^T + A^T*A (off-diagonal)")
print("Let C = A*A^T + A^T*A. Then Sigma[u,v] = C[u,v] for u!=v, Sigma[u,u]=0.")
print()

# So C[u,v] = (AAT)[u,v] + (ATA)[u,v]
# The matrix C - diag(C) = Sigma (the sigma graph)
# C = AAT + ATA
# Note: AAT + ATA = A*A^T + A^T*A. For tournament: A + A^T = J - I.
# So AAT + ATA = A(J-I-A) + (J-I-A)A = AJ - A - A^2 + JA - A - A^2
# = AJ + JA - 2A - 2A^2
# AJ = (row sums of A) = score matrix. AJ[u,v] = score(u) for all v.
# JA[u,v] = column sum of A at v = in-degree(v) = n-1-score(v) for all u.
# So C = diag(scores)*J + J*diag(n-1-scores) - 2A - 2A^2
# Hmm, this is getting complicated. Let me just verify numerically.

test_bits = np.random.randint(0, 1 << total_bits)
A_test = bits_to_adj(test_bits, n)
L_test = lambda_graph(A_test, n)
A2_test = A_test @ A_test
AAT = A_test @ A_test.T
ATA = A_test.T @ A_test
C = AAT + ATA

# Build Sigma
sig_test = np.zeros((n, n), dtype=int)
for u in range(n):
    for v in range(u+1, n):
        sig_test[u][v] = sig_test[v][u] = n - 2 - int(A2_test[u][v]) - int(A2_test[v][u])

print("Verification: Sigma = C - diag(C) ?")
C_offdiag = C.copy()
np.fill_diagonal(C_offdiag, 0)
# sigma(u,v) = n-2 - (A^2)[u,v] - (A^2)[v,u]
# C[u,v] = (AAT)[u,v] + (ATA)[u,v]
# But (A^2)[u,v] != (AAT)[u,v] in general!
# (A^2)[u,v] = sum_w A[u,w]*A[w,v]
# (AAT)[u,v] = sum_w A[u,w]*A^T[w,v] = sum_w A[u,w]*A[v,w]
# These are DIFFERENT!

# sigma = common_succ + common_pred = AAT[u,v] + ATA[u,v]
# Let me check: is Sigma == C_offdiag?
match = np.array_equal(sig_test, C_offdiag)
print(f"  Sigma == C_offdiag: {match}")

# If not, what is the relationship?
if not match:
    print(f"  First mismatch: Sigma[0,1]={sig_test[0,1]}, C[0,1]={C_offdiag[0,1]}")
    # Check: sigma = n-2 - A2[u,v] - A2[v,u]
    # and C[u,v] = AAT[u,v] + ATA[u,v]
    # Are these the same? n-2 - A2[u,v] - A2[v,u] vs AAT[u,v] + ATA[u,v]
    # A2[u,v] + A2[v,u] = #{u->w->v} + #{v->w->u} = #{transitive triples + cyclic triples?}
    # AAT[u,v] + ATA[u,v] = common_succ + common_pred
    # For tournament pair (u,v) with n-2 other vertices:
    # Each w is classified as: succ of both, pred of both, cyclic, or transitive
    # common_succ + common_pred = sigma
    # transitive = #{u->w->v when u->v} + #{v->w->u when v->u} = delta
    # Actually let me be more careful. Not right either.

    # Let's just check:
    u, v = 0, 1
    cs = sum(1 for w in range(n) if w != u and w != v and A_test[u][w] and A_test[v][w])
    cp = sum(1 for w in range(n) if w != u and w != v and A_test[w][u] and A_test[w][v])
    a2uv = int(A2_test[u][v])
    a2vu = int(A2_test[v][u])
    print(f"  u={u},v={v}: cs={cs}, cp={cp}, A2[u,v]={a2uv}, A2[v,u]={a2vu}")
    print(f"  cs+cp={cs+cp}, n-2-A2[u,v]-A2[v,u]={n-2-a2uv-a2vu}")
    print(f"  AAT[u,v]={AAT[u][v]}, ATA[u,v]={ATA[u][v]}")

    # Check if cs = AAT[u,v]
    print(f"  cs == AAT[u,v]: {cs == AAT[u][v]}")
    print(f"  cp == ATA[u,v]: {cp == ATA[u][v]}")

# 7. The true identity for sigma
# For a tournament: every w != u,v falls into one of 4 categories:
#   S: u->w AND v->w (common successor)
#   P: w->u AND w->v (common predecessor)
#   L: u->v->w->u OR v->u->w->v (cyclic, 3-cycle witness)
#   D: u->w->v (u->v case) or v->w->u (v->u case) (transitive triple)
# sigma = S + P, lambda = L, delta = D, and S+P+L+D = n-2

# Now: A^2[u,v] = #{w: u->w AND w->v} = #{w transitively between u and v}
#              = D (when u->v) or S (when v->u)? No...
# A^2[u,v] = #{u->w->v} = #{w: A[u,w]=1 and A[w,v]=1}
# If u->v:
#   w is S-type (u->w, v->w): contributes to A^2[u,v] iff w->v, but v->w, so no
#   w is P-type (w->u, w->v): contributes iff u->w, but w->u, so no
#   w is L-type (3-cycle): u->v->w->u: contributes iff u->w and w->v, which needs w->u->v->w reversed...
#     Actually L-type for pair (u,v) with u->v means either u->v->w->u or v->u->... wait
#     For u->v: L counts w where (u,v,w) form a 3-cycle.
#     That's u->v->w->u, meaning A[u,v]=1, A[v,w]=1, A[w,u]=1.
#     So A^2[u,v] contribution from L: need u->w AND w->v. But w->u (from 3-cycle), not u->w.
#     So L does NOT contribute to A^2[u,v].
#   w is D-type (transitive): u->w->v (consistent with u->v).
#     A^2[u,v]: u->w AND w->v = exactly the D-type!
# So when u->v: A^2[u,v] = D (the transitive witness count for pair u->v)

# Similarly: A^2[v,u] = #{v->w->u}.
# If u->v: this needs v->w AND w->u.
#   S-type: u->w, v->w: v->w yes, w->u? w->u iff u doesn't beat w. But u->w. No.
#   P-type: w->u yes, w->v? w->v iff v doesn't beat w. But v->w... wait P-type means w->u AND w->v.
#     So v->w->u? w->u yes. v->w? No, w->v (P-type means w beats both).
#   L-type: 3-cycle u->v->w->u: need v->w AND w->u. v->w yes (from cycle), w->u yes. YES!
#   D-type: u->w->v: need v->w? No, w->v (from D). So no.
# So when u->v: A^2[v,u] = L (the lambda count!)

# This gives us: sigma = n-2 - A^2[u,v] - A^2[v,u] = n-2 - D - L = S+P = sigma. Correct!

print("\n\n--- Key Identity ---")
print("For pair (u,v) with u->v:")
print("  A^2[u,v] = delta(u,v) (transitive witnesses)")
print("  A^2[v,u] = lambda(u,v) (cyclic witnesses = 3-cycle count)")
print("  sigma(u,v) = n-2 - delta - lambda = n-2 - A^2[u,v] - A^2[v,u]")
print()
print("This means the sigma graph encodes COMPLEMENTARY information to A^2!")
print("Sigma = (n-2)*J_off - (A^2 + (A^2)^T)_off")
print("where J_off = all-ones with zero diagonal, and _off means zero diagonal.")

# Verify
print("\nVerification:")
J_off = np.ones((n, n), dtype=int)
np.fill_diagonal(J_off, 0)
A2_sym_off = (A2_test + A2_test.T).copy()
np.fill_diagonal(A2_sym_off, 0)
Sigma_formula = (n-2) * J_off - A2_sym_off
print(f"  Sigma == (n-2)*J_off - (A^2 + (A^2)^T)_off: {np.array_equal(sig_test, Sigma_formula)}")

# So Sigma = (n-2)*J_off - sym(A^2)_off
# tr(Sigma^k) can be expressed in terms of tr((sym(A^2))^k)
# and mixed terms with J_off.

# Since J_off has known spectrum ({n-1} with mult 1, {-1} with mult n-1 for n=7),
# the sigma spectrum is a shifted version of the sym(A^2) spectrum.

# Actually: Sigma = c*J_off - M where c=n-2, M = sym(A^2)_off
# J_off has eigenvalues: n-1 (eigvec = all-1s), -1 (n-1 times)
# So for any matrix X: Sigma*X = c*J_off*X - M*X
# The spectrum of Sigma is NOT simply related to M unless they commute.

# But we DO know: tr(Sigma^k) = tr((cJ_off - M)^k)
# which expands via binomial... this is messy in general.

# Better approach: just fit c7 as a polynomial in the data
# Since c7 is determined by (tr2, tr3, tr4) and we have exact integer data,
# let's find the EXACT rational coefficients.

print("\n\n--- Exact Coefficient Search ---")

# Collect unique (tr_s2, tr_s3, tr_s4, c7) tuples
unique_data = set()
for d in data_points:
    unique_data.add((d['tr_s2'], d['tr_s3'], d['tr_s4'], d['c7']))

print(f"Unique (tr2, tr3, tr4, c7) tuples: {len(unique_data)}")

# Check if c7 is a POLYNOMIAL in (tr2, tr3, tr4)
# Try: c7 = a*tr2 + b*tr3 + c*tr4 + d
from fractions import Fraction

# Use first 10 data points to solve, then verify on rest
unique_list = sorted(unique_data)

# Linear: c7 = a*tr2 + b*tr3 + c*tr4 + d
X_lin = np.array([[t2, t3, t4, 1] for t2, t3, t4, c7 in unique_list], dtype=float)
y_lin = np.array([c7 for t2, t3, t4, c7 in unique_list], dtype=float)
coeffs_lin, _, _, _ = np.linalg.lstsq(X_lin, y_lin, rcond=None)
pred_lin = X_lin @ coeffs_lin
rmse_lin = np.sqrt(np.mean((y_lin - pred_lin)**2))
print(f"\nLinear fit c7 = a*tr2 + b*tr3 + c*tr4 + d:")
print(f"  a={coeffs_lin[0]:.8f}, b={coeffs_lin[1]:.8f}, c={coeffs_lin[2]:.8f}, d={coeffs_lin[3]:.8f}")
print(f"  RMSE: {rmse_lin:.6f}")

# Check if coefficients are rational with small denominators
for coeff in coeffs_lin:
    frac = Fraction(coeff).limit_denominator(1000)
    print(f"  {coeff:.8f} ~= {frac} = {float(frac):.8f}")

# Also include c3 (= total_lambda / 3):
print("\n\nLinear fit c7 = a*c3 + b*tr2 + c*tr3 + d*tr4 + e:")
X_c3lin = np.array([[d['c3'], d['tr_s2'], d['tr_s3'], d['tr_s4'], 1]
                      for d in data_points], dtype=float)
y_c3 = np.array([d['c7'] for d in data_points], dtype=float)
coeffs_c3lin, _, _, _ = np.linalg.lstsq(X_c3lin, y_c3, rcond=None)
pred_c3lin = X_c3lin @ coeffs_c3lin
rmse_c3lin = np.sqrt(np.mean((y_c3 - pred_c3lin)**2))
print(f"  a={coeffs_c3lin[0]:.8f}, b={coeffs_c3lin[1]:.8f}, c={coeffs_c3lin[2]:.8f}, d={coeffs_c3lin[3]:.8f}, e={coeffs_c3lin[4]:.8f}")
print(f"  RMSE: {rmse_c3lin:.6f}")
for i, coeff in enumerate(coeffs_c3lin):
    frac = Fraction(coeff).limit_denominator(10000)
    print(f"  coeff[{i}] = {coeff:.8f} ~= {frac}")

# Max absolute error
max_err = np.max(np.abs(y_c3 - pred_c3lin))
print(f"  Max error: {max_err:.4f}")

# Check roundability: if c7 is exact integer and prediction is close
near_int = np.all(np.abs(y_c3 - np.round(pred_c3lin)) < 0.01)
print(f"  Predictions round to correct c7: {near_int}")

print("\n\nDone.")
