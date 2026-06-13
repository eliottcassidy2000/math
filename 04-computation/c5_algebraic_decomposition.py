"""
c5_algebraic_decomposition.py -- kind-pasteur-2026-03-13-S61

ALGEBRAIC PROOF that c5_dir is lambda-determined.

Key identity: tr(A^5) = 5*c5_dir.
Goal: express tr(A^5) as a polynomial in lambda(u,v) values.

APPROACH: We know A[u][v]^2 = A[u][v] (binary matrix).
And: A[u][v] + A[v][u] = 1 (tournament).
And: lambda(u,v) = sum_{w != u,v} [A[u][v]*A[v][w]*A[w][u] + A[v][u]*A[u][w]*A[w][v]]

Let's define:
  a(u,v) = A[u][v] - 1/2 in {-1/2, +1/2}  (centered)
  Then A[u][v] = 1/2 + a(u,v), a(u,v) = -a(v,u).

tr(A^5) = sum_{(i,j,k,l,m) distinct} A[i][j]*A[j][k]*A[k][l]*A[l][m]*A[m][i]
         = sum_{(i,j,k,l,m) distinct} prod of 5 arc indicators

Expand A = 1/2 + a:
A[i][j]*A[j][k]*A[k][l]*A[l][m]*A[m][i]
= (1/2 + a_{ij})(1/2 + a_{jk})(1/2 + a_{kl})(1/2 + a_{lm})(1/2 + a_{mi})

This is a degree-5 polynomial in the a values. The terms of each degree:
- Degree 0: (1/2)^5 = 1/32
- Degree 1: (1/2)^4 * (a_{ij} + a_{jk} + a_{kl} + a_{lm} + a_{mi})
  Sum over all 5-tuples: each a_{uv} appears in many 5-tuples, but
  sum of a_{uv} over all distinct (u,v) is 0 (since a is antisymmetric).
  Actually need careful counting.
- Degree 2: ...
- ...

This expansion is known in the literature. Let me instead verify computationally
whether tr(A^5) = f(lambda) via a DIRECT algebraic identity.

ALTERNATIVE: Use the known identity
  sum_{w} A[u][w]*A[w][v] = number of 2-paths u->w->v
  = number of w with u->w and w->v.

For tournaments: sum_w A[u][w]*A[w][v] = (A^2)[u][v].
And lambda(u,v) relates to the 3-cycle count:
  lambda(u,v) = #{w : (u,v,w) is a 3-cycle}
  = #{w : A[u][v]*A[v][w]*A[w][u] = 1} + #{w : A[v][u]*A[u][w]*A[w][v] = 1}
  = A[u][v] * sum_w A[v][w]*A[w][u] + A[v][u] * sum_w A[u][w]*A[w][v]
  = A[u][v] * (A^2)[v][u] + A[v][u] * (A^2)[u][v]

So: lambda(u,v) = A[u][v]*(A^2)[v][u] + (1-A[u][v])*(A^2)[u][v]
                = A[u][v]*((A^2)[v][u] - (A^2)[u][v]) + (A^2)[u][v]

Let D = A^2. Then:
  lambda(u,v) = A[u][v]*(D[v][u] - D[u][v]) + D[u][v]

Now: D[u][v] + D[v][u] = #{w : u->w->v} + #{w : v->w->u}
                        = #{w != u,v} = n-2

So D[v][u] = n-2 - D[u][v], and:
  lambda(u,v) = A[u][v]*(n-2-2*D[u][v]) + D[u][v]
  = (n-2)*A[u][v] - 2*A[u][v]*D[u][v] + D[u][v]

This gives: D[u][v] = (lambda - (n-2)*A[u][v]) / (1 - 2*A[u][v])
when A[u][v] != 1/2 (which is always true for tournaments).

If A[u][v]=1: D[u][v] = (lambda - (n-2)) / (-1) = n-2 - lambda
If A[u][v]=0: D[u][v] = lambda / 1 = lambda

So: (A^2)[u][v] = lambda(u,v)         when A[u][v]=0  (u does NOT beat v)
    (A^2)[u][v] = n-2-lambda(u,v)     when A[u][v]=1  (u beats v)

Equivalently: (A^2)[u][v] = (n-2)*A[u][v] + lambda(u,v)*(1-2*A[u][v])
Wait, let me re-derive:
  If A[u][v]=1: D[u][v] = #{w: u->w->v} = #{w: A[u][w]=1 and A[w][v]=1}
  lambda(u,v) = 1*(D[v][u]) + 0*D[u][v] = D[v][u] = n-2 - D[u][v]
  So D[u][v] = n-2 - lambda(u,v)

  If A[u][v]=0 (so A[v][u]=1): D[u][v] = #{w: u->w->v}
  lambda(u,v) = 0*D[v][u] + 1*D[u][v] = D[u][v]
  So D[u][v] = lambda(u,v)

THIS IS KEY: (A^2)[u][v] depends on both A[u][v] and lambda(u,v).
So A^2 = f(A, lambda). And A^3 = A * A^2, A^4 = A * A^3, etc.
tr(A^5) = tr(A * A^4) = sum_i (A * A^4)[i][i] = sum_{i,j} A[i][j]*(A^4)[j][i]

But A^k still involves A, not just lambda. So tr(A^5) being lambda-determined
means that the A-dependence cancels out in the trace.

Let me verify this algebraically at n=5 and also check at n=7 whether
tr(A^7) = 7*c7_dir is NOT lambda-determined (it shouldn't be, since c7 isn't).
"""

import numpy as np
from itertools import combinations, permutations
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

# Verify the key identity: (A^2)[u][v] = lambda(u,v) when u does not beat v
#                           (A^2)[u][v] = n-2-lambda(u,v) when u beats v
print("=" * 60)
print("VERIFYING: A^2 = f(A, lambda)")
print("=" * 60)

n = 5
total_bits = n * (n-1) // 2
success = 0
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A2 = A @ A
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            if A[u][v] == 1:
                expected = n - 2 - L[u][v]
            else:
                expected = L[u][v]
            assert A2[u][v] == expected, f"FAIL at bits={bits}, ({u},{v})"
    success += 1
print(f"A^2 identity verified for all {success} tournaments at n={n}")

# Now: (A^2)[u][v] = (n-2)*A[u][v] + lambda(u,v) * (1 - 2*A[u][v])
# In centered variables a = A - 1/2:
# A = 1/2 + a, 1-2*A = -2a
# (A^2)[u][v] = (n-2)*(1/2+a[u][v]) + lambda[u][v]*(-2*a[u][v])
#             = (n-2)/2 + (n-2)*a[u][v] - 2*lambda[u][v]*a[u][v]
#             = (n-2)/2 + (n-2-2*lambda[u][v]) * a[u][v]

# So A^2 = (n-2)/2 * J + diag(A^2[i][i]) + [(n-2) - 2*Lambda] * a
# where * is Hadamard product and J is all-ones (off-diag).
# Wait, this is entry-wise for u != v.

# For diagonal: (A^2)[v][v] = sum_w A[v][w]^2 = sum_w A[v][w] = s_v (out-degree)

# Now tr(A^5) = tr(A * A^4) = tr(A * (A^2)^2)
# Let D = A^2. Then A^4 = D^2 and tr(A^5) = sum_i (A*D^2)[i][i]
#                                           = sum_{i,j} A[i][j] * (D^2)[j][i]
#                                           = sum_{i,j,k} A[i][j]*D[j][k]*D[k][i]

# Each D[j][k] depends on A[j][k] and lambda[j][k].
# D[j][k] = (n-2)*A[j][k] - (2*A[j][k]-1)*lambda[j][k]
#          = when A[j][k]=1: n-2-lambda  =  n-2-L[j][k]
#          = when A[j][k]=0: lambda      = L[j][k]

# So tr(A^5) = sum_{i,j,k} A[i][j] * g(A[j][k], L[j][k]) * g(A[k][i], L[k][i])
# where g(1, l) = n-2-l, g(0, l) = l.
# Note g(a, l) = l + a*(n-2-2l) = l*(1-2a) + a*(n-2)

# This triple sum still involves A[i][j]. But tr(A^5) is lambda-determined!
# So the sum over A[i][j] must cancel...

# Actually, let's decompose differently. Define:
# B[u][v] = A[u][v] * (n-2-2*lambda(u,v))  for u != v
# Then D[u][v] = (n-2)/2 + B[u][v] - lambda(u,v) + lambda(u,v)
# Hmm, this is getting messy. Let me try a direct computational approach.

# Can we express tr(A^5) purely in lambda entries?
# tr(A^5) = 5*c5. And c5 is lambda-determined.
# So tr(A^5) IS a function of the lambda matrix L.
# But it's not a POLYNOMIAL in L entries alone (it also involves A).

# The question is whether there's a POLYNOMIAL identity expressing tr(A^5) in terms of L.
# From the histogram analysis: c5 = cubic polynomial in (n0,n1,n2,n3).
# The nk are themselves polynomial functions of the L entries
# (e.g., n0 = #{(u,v): L[u][v]=0} = sum_{u<v} delta(L[u][v], 0)).
# But delta(L[u][v], 0) is NOT a polynomial!

# So the formula c5 = f(histogram) is NOT polynomial in L entries.
# But there must be SOME function (polynomial or not) of L that gives c5.

# Let me try: express c5 as a polynomial in individual L[u][v] entries.
# At n=5, there are C(5,2)=10 entries. c5 is invariant under S_5 action.
# A general polynomial in 10 variables... we need to find the right one.

print(f"\n{'='*60}")
print("EXPRESSING c5 AS POLYNOMIAL IN L[u][v] ENTRIES")
print(f"{'='*60}")

# Build feature matrix: all monomials in L[i][j] up to some degree
# Degree 1: L[0][1], L[0][2], ..., L[3][4] — 10 features
# Degree 2: L[0][1]^2, L[0][1]*L[0][2], ... — C(10+1,2) = 55 features
# Degree 3: ... — C(10+2,3) = 220 features

# First try degree 1 (linear in individual lambda entries)
n = 5
total_bits = n * (n-1) // 2

pairs = [(i, j) for i in range(n) for j in range(i+1, n)]  # 10 pairs

all_L = []
all_c5 = []
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A5 = np.linalg.matrix_power(A, 5)
    c5 = int(np.trace(A5)) // 5
    all_L.append([L[i][j] for i, j in pairs])
    all_c5.append(c5)

X_L = np.array(all_L, dtype=float)
y = np.array(all_c5, dtype=float)

# Degree 1
from numpy.linalg import lstsq, matrix_rank

X1 = np.column_stack([np.ones(len(y)), X_L])
c1, _, _, _ = lstsq(X1, y, rcond=None)
err1 = np.max(np.abs(X1 @ c1 - y))
print(f"Degree 1 (11 features): max error = {err1:.6f}")

# Degree 2
X2_parts = [np.ones(len(y))]
X2_parts.append(X_L)
# Add all products L[i][j] * L[k][l]
for i in range(10):
    for j in range(i, 10):
        X2_parts.append((X_L[:, i] * X_L[:, j]).reshape(-1, 1))
X2 = np.column_stack(X2_parts)
c2, _, _, _ = lstsq(X2, y, rcond=None)
err2 = np.max(np.abs(X2 @ c2 - y))
print(f"Degree 2 ({X2.shape[1]} features): rank={matrix_rank(X2)}, max error = {err2:.6f}")

# Degree 3 — add triple products
X3_parts = list(X2_parts)
for i in range(10):
    for j in range(i, 10):
        for k in range(j, 10):
            X3_parts.append((X_L[:, i] * X_L[:, j] * X_L[:, k]).reshape(-1, 1))
X3 = np.column_stack(X3_parts)
c3, _, _, _ = lstsq(X3, y, rcond=None)
err3 = np.max(np.abs(X3 @ c3 - y))
rank3 = matrix_rank(X3)
print(f"Degree 3 ({X3.shape[1]} features): rank={rank3}, max error = {err3:.6f}")

# The key question: does degree 2 suffice?
if err2 < 0.001:
    print("\n  *** DEGREE 2 IS EXACT! ***")
    print("  c5_dir is a QUADRATIC polynomial in lambda entries.")
elif err3 < 0.001:
    print("\n  Degree 3 needed.")

# Now check: is the formula SYMMETRIC?
# c5 should be invariant under vertex permutations.
# But the lambda entries L[i][j] transform under S_5.
# A symmetric polynomial in L[i][j] is one that's invariant under the induced S_5 action.
# The power sums S_k = sum L[i][j]^k are symmetric.
# The vertex degree sums D_k = sum_v (sum_u L[v][u])^k are symmetric.
# Products of these are symmetric.

# But a degree-2 polynomial in 10 variables might have non-symmetric terms
# that combine to give a symmetric result.

# Let me verify: is c5 a SYMMETRIC quadratic in L[i][j]?
# A general symmetric quadratic in the 10 lambda values has the form:
# c5 = a + b * sum L[i][j] + c * sum L[i][j]^2 + d * sum_{(ij),(kl) disjoint} L[ij]*L[kl]
#      + e * sum_{ij, kl sharing 1 vertex} L[ij]*L[kl]

# Actually, symmetric quadratics in lambda entries under S_5 action:
# Orbits of pairs of pairs under S_5:
# 1. Same pair: (ij, ij) — contributes sum L[ij]^2
# 2. Disjoint pairs: (ij, kl), {i,j} cap {k,l} = {} — C(5,4)*3 = 15
# 3. Sharing 1 vertex: (ij, ik), j != k — 5*C(4,2) = 30

# So symmetric quadratic has 5 free parameters: [1, sum L, sum L^2, sum_disjoint L*L, sum_sharing L*L]

print(f"\n{'='*60}")
print("SYMMETRIC QUADRATIC FIT")
print(f"{'='*60}")

# Compute the symmetric quadratic features
sym_features = []
for idx in range(len(all_L)):
    L_vals = all_L[idx]
    # 1
    f0 = 1
    # sum L[ij]
    f1 = sum(L_vals)
    # sum L[ij]^2
    f2 = sum(l**2 for l in L_vals)
    # sum over disjoint pairs
    f3 = 0
    for a in range(10):
        for b in range(a+1, 10):
            pa = pairs[a]
            pb = pairs[b]
            if len(set(pa) & set(pb)) == 0:
                f3 += L_vals[a] * L_vals[b]
    # sum over pairs sharing 1 vertex
    f4 = 0
    for a in range(10):
        for b in range(a+1, 10):
            pa = pairs[a]
            pb = pairs[b]
            if len(set(pa) & set(pb)) == 1:
                f4 += L_vals[a] * L_vals[b]
    sym_features.append([f0, f1, f2, f3, f4])

X_sym = np.array(sym_features, dtype=float)
c_sym, _, _, _ = lstsq(X_sym, y, rcond=None)
err_sym = np.max(np.abs(X_sym @ c_sym - y))
print(f"Features: [1, sum_L, sum_L2, sum_disjoint_LL, sum_sharing_LL]")
print(f"Rank: {matrix_rank(X_sym)}")
print(f"Max error: {err_sym:.6f}")
print(f"Coefficients: {[f'{c:.6f}' for c in c_sym]}")

if err_sym < 0.001:
    print("  *** SYMMETRIC QUADRATIC IS EXACT! ***")
    for denom in range(1, 1201):
        int_c = np.round(c_sym * denom)
        err_r = np.max(np.abs(X_sym @ (int_c / denom) - y))
        if err_r < 0.001:
            a0, a1, a2, a3, a4 = [int(x) for x in int_c]
            print(f"\n  EXACT FORMULA:")
            print(f"  c5 = ({a0} + {a1}*S1 + {a2}*sum(lam^2) + {a3}*sum_disj(l*l) + {a4}*sum_share(l*l)) / {denom}")
            break
else:
    # Try adding degree-3 symmetric terms
    print(f"\nNot exact. Adding degree-3 symmetric terms...")

    # Symmetric degree-3 monomials in lambda entries under S_5:
    # Triple orbits:
    # 1. Same pair cubed: sum L[ij]^3
    # 2. Product of pair with disjoint pair: sum L[ij]^2 * L[kl] where {ij} disjoint {kl}
    # 3. Product of pair with sharing pair: sum L[ij]^2 * L[ik]
    # 4. Three mutually disjoint pairs: IMPOSSIBLE (need 6 vertices, only have 5)
    # 5. Triangle: sum L[ij]*L[jk]*L[ik] (three pairs forming a triangle)
    # 6. Path: sum L[ij]*L[ik]*L[jl] (three pairs sharing vertex arrangement)
    # ... this gets complicated. Let me just use the degree-3 power sums + D sums.

    # Actually, let me try: [1, S1, S2, S3, D2, trL3, sum_disj, sum_share]
    sym3_features = []
    for idx in range(len(all_L)):
        L_mat = np.zeros((n, n), dtype=int)
        for pi, (i, j) in enumerate(pairs):
            L_mat[i][j] = all_L[idx][pi]
            L_mat[j][i] = all_L[idx][pi]

        L_vals = all_L[idx]
        f0 = 1
        f1 = sum(L_vals)
        f2 = sum(l**2 for l in L_vals)
        f3 = sum(l**3 for l in L_vals)
        d_v = [sum(L_mat[v][u] for u in range(n) if u != v) for v in range(n)]
        f4 = sum(dv**2 for dv in d_v)
        f5 = int(np.trace(L_mat @ L_mat @ L_mat))

        # Disjoint pair product
        f6 = 0
        for a in range(10):
            for b in range(a+1, 10):
                if len(set(pairs[a]) & set(pairs[b])) == 0:
                    f6 += L_vals[a] * L_vals[b]

        # Sharing pair product
        f7 = 0
        for a in range(10):
            for b in range(a+1, 10):
                if len(set(pairs[a]) & set(pairs[b])) == 1:
                    f7 += L_vals[a] * L_vals[b]

        # Triangle product
        f8 = 0
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    f8 += L_mat[a][b] * L_mat[b][c] * L_mat[a][c]

        sym3_features.append([f0, f1, f2, f3, f4, f5, f6, f7, f8])

    X_sym3 = np.array(sym3_features, dtype=float)
    c_sym3, _, _, _ = lstsq(X_sym3, y, rcond=None)
    err_sym3 = np.max(np.abs(X_sym3 @ c_sym3 - y))
    names3 = ['1', 'S1', 'S2', 'S3', 'D2', 'trL3', 'sum_disj', 'sum_share', 'triangle']
    print(f"Features: {names3}")
    print(f"Rank: {matrix_rank(X_sym3)}")
    print(f"Max error: {err_sym3:.6f}")

    if err_sym3 < 0.001:
        print("  *** SYMMETRIC CUBIC IS EXACT! ***")
        for denom in range(1, 2001):
            int_c = np.round(c_sym3 * denom)
            err_r = np.max(np.abs(X_sym3 @ (int_c / denom) - y))
            if err_r < 0.001:
                print(f"\n  EXACT FORMULA (denom={denom}):")
                terms = []
                for i, name in enumerate(names3):
                    c = int(int_c[i])
                    if c != 0:
                        terms.append(f"({c:+d})*{name}")
                print(f"  c5 = ({' '.join(terms)}) / {denom}")

                # Verify
                verified = True
                for bits in range(1 << total_bits):
                    A = bits_to_adj(bits, n)
                    L = lambda_graph(A, n)
                    A5 = np.linalg.matrix_power(A, 5)
                    c5_true = int(np.trace(A5)) // 5

                    L_vals = [L[i][j] for i, j in pairs]
                    feat = [1]
                    feat.append(sum(L_vals))
                    feat.append(sum(l**2 for l in L_vals))
                    feat.append(sum(l**3 for l in L_vals))
                    d_v = [sum(L[v][u] for u in range(n) if u != v) for v in range(n)]
                    feat.append(sum(dv**2 for dv in d_v))
                    feat.append(int(np.trace(L @ L @ L)))
                    disj = sum(L_vals[a]*L_vals[b] for a in range(10)
                               for b in range(a+1, 10)
                               if len(set(pairs[a]) & set(pairs[b])) == 0)
                    feat.append(disj)
                    share = sum(L_vals[a]*L_vals[b] for a in range(10)
                                for b in range(a+1, 10)
                                if len(set(pairs[a]) & set(pairs[b])) == 1)
                    feat.append(share)
                    tri = 0
                    for a in range(n):
                        for b in range(a+1, n):
                            for c_idx in range(b+1, n):
                                tri += L[a][b] * L[b][c_idx] * L[a][c_idx]
                    feat.append(tri)

                    c5_pred = sum(int_c[i] * feat[i] for i in range(len(feat))) / denom
                    if abs(c5_pred - c5_true) > 0.001:
                        verified = False
                        print(f"  FAIL at bits={bits}")
                        break
                print(f"  Verified on all 1024: {verified}")
                break

    # Also try: replace trL3 with sum_triangle (they may be related)
    # trL3 = sum_{i,j,k} L[i][j]*L[j][k]*L[k][i]
    # For symmetric L with 0 diagonal:
    # trL3 = sum_{i!=j!=k!=i} L[i][j]*L[j][k]*L[k][i]
    # = 6 * sum_{i<j<k} L[i][j]*L[j][k]*L[i][k]  (each triangle counted 6 times)
    # So trL3 = 6 * triangle!

    print(f"\n--- Check: trL3 = 6 * triangle? ---")
    for bits in range(5):
        A = bits_to_adj(bits, n)
        L = lambda_graph(A, n)
        trL3 = int(np.trace(L @ L @ L))
        tri = 0
        for a in range(n):
            for b in range(a+1, n):
                for c_idx in range(b+1, n):
                    tri += L[a][b] * L[b][c_idx] * L[a][c_idx]
        print(f"  bits={bits}: trL3={trL3}, 6*triangle={6*tri}, equal={trL3 == 6*tri}")

print("\nDone.")
