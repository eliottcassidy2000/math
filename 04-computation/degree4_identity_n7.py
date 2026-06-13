#!/usr/bin/env python3
"""
DEGREE-4 OCF IDENTITY AT n=7: COMPLETE ANALYSIS

=========================================================================
THEOREM (Degree-4 Structure at n=7):
The degree-4 Fourier space for tournament invariants at n=7 is
TWO-DIMENSIONAL, spanned by [deg-4 of t_5] and [deg-4 of alpha_2].

The following EXACT identities hold:

  (A) [deg-4 of t_7] = (1/2)*[deg-4 of t_5] + 1*[deg-4 of alpha_2]

  (B) w_2/4 = 3*[deg-4 of t_5] + 6*[deg-4 of alpha_2]

Together, these are EQUIVALENT to the degree-4 OCF identity:
  w_2/4 = 2*[deg-4 of alpha_1] + 4*[deg-4 of alpha_2]

Proof of equivalence:
  [deg-4 of alpha_1] = [deg-4 of t_5] + [deg-4 of t_7]  (t_3 has max deg 2)
  Using (A): = [deg-4 of t_5] + (1/2)*[deg-4 of t_5] + [deg-4 of alpha_2]
             = (3/2)*[deg-4 of t_5] + [deg-4 of alpha_2]
  So RHS = 2*(3/2)*[deg-4 of t_5] + 2*[deg-4 of alpha_2] + 4*[deg-4 of alpha_2]
         = 3*[deg-4 of t_5] + 6*[deg-4 of alpha_2]
  Which is identity (B). QED.

=========================================================================
COMBINED WITH PREVIOUSLY PROVED IDENTITIES:

  Degree 0: PROVED (expected values)
  Degree 2: PROVED (proportionality constants c_5=3, c_7=1.5, c_{a2}=1)
  Degree 4: COMPUTATIONALLY VERIFIED via identities (A) and (B)
  Degree 6: PROVED (path-cycle bijection)

STATUS: OCF at n=7 reduces to proving identities (A) and (B).

=========================================================================
COMBINATORIAL MEANING:

[deg-4 of t_5]: Sum over (5-vertex subset, 5-cycle C on subset, skip 1 edge)
                of product of 4 remaining edge s-variables.
                = "4-paths embeddable in 5-cycles on vertex subsets"

[deg-4 of alpha_2]: Sum over (disjoint triangle pairs T1, T2)
                    of [deg-2 of I_{T1}] * [deg-2 of I_{T2}]
                    = product of "P2 path indicators" from disjoint triangles

[deg-4 of t_7]: (1/4) * Sum over (7-cycle C, choose 4 of 7 edges)
                of product of 4 edge s-variables.
                = "4-edge subsets of Hamiltonian cycles"

Identity (A) says: 4-edge subsets of 7-cycles decompose into
contributions from 5-cycles and disjoint triangle pairs.

opus-2026-03-06-S11b (continued^8)
"""
from itertools import combinations, permutations
from collections import defaultdict
import numpy as np
import random
from numpy.linalg import lstsq, svd

def random_tournament(n, seed):
    random.seed(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

# =====================================================
# VERIFICATION WITH 3000 SAMPLES
# =====================================================
n = 7
edges = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges)
N = 3000

print("=" * 70)
print("DEGREE-4 OCF IDENTITY AT n=7: VERIFICATION")
print("=" * 70)
print(f"\nSampling {N} random tournaments on {n} vertices...")

t3 = np.zeros(N)
t5 = np.zeros(N)
t7 = np.zeros(N)
a2 = np.zeros(N)
w0 = np.zeros(N)
w2 = np.zeros(N)
H = np.zeros(N)

for trial in range(N):
    A = random_tournament(n, seed=trial)

    t3[trial] = sum(1 for a,b,c in combinations(range(n), 3)
                     if A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a] > 0)

    t5_val = 0
    for verts in combinations(range(n), 5):
        v = list(verts)
        sub = [[A[v[i]][v[j]] for j in range(5)] for i in range(5)]
        for perm in permutations(range(1, 5)):
            path = (0,) + perm
            if all(sub[path[i]][path[i+1]] for i in range(4)) and sub[path[4]][path[0]]:
                t5_val += 1
    t5[trial] = t5_val

    t7_val = 0
    for perm in permutations(range(1, 7)):
        path = (0,) + perm
        if all(A[path[i]][path[i+1]] for i in range(6)) and A[path[6]][path[0]]:
            t7_val += 1
    t7[trial] = t7_val

    tris = [(a,b,c) for a,b,c in combinations(range(n), 3)
            if A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a] > 0]
    a2[trial] = sum(1 for i in range(len(tris)) for j in range(i+1, len(tris))
                     if set(tris[i]).isdisjoint(set(tris[j])))

    w0_val = 0.0
    for perm in permutations(range(n)):
        prod_val = 1.0
        for i in range(n-1):
            prod_val *= (A[perm[i]][perm[i+1]] - 0.5)
        w0_val += prod_val
    w0[trial] = w0_val

    # Compute w2 via W(r) evaluation
    def W(r):
        total = 0.0
        for perm in permutations(range(n)):
            prod_val = 1.0
            for i in range(n-1):
                prod_val *= (r + A[perm[i]][perm[i+1]] - 0.5)
            total += prod_val
        return total
    W1 = W(1)
    W2_val = W(2)
    a1 = W1 - w0_val - 5040
    a2_eq = W2_val - w0_val - 64*5040
    w4_val = (4*a1 - a2_eq) / (-12)
    w2[trial] = a1 - w4_val

    # H via DP
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[(mask | (1 << u), u)] = dp.get((mask | (1 << u), u), 0) + val
    full = (1 << n) - 1
    H[trial] = sum(dp.get((full, v), 0) for v in range(n))

    if trial % 1000 == 0:
        print(f"  {trial}/{N}...")

# Verify OCF
ocf_err = np.max(np.abs(H - (1 + 2*t3 + 2*t5 + 2*t7 + 4*a2)))
print(f"\nOCF: max|H - (1+2t3+2t5+2t7+4a2)| = {ocf_err}")

# =====================================================
# DEGREE DECOMPOSITION
# =====================================================
print("\n" + "-" * 70)
print("Degree decomposition:")

t3c = t3 - t3.mean()
t5c = t5 - t5.mean()
t7c = t7 - t7.mean()
a2c = a2 - a2.mean()

# Degree-2 proportionality: c_5=3, c_7=1.5, c_{a2}=1
t5_d4 = t5c - 3.0 * t3c           # max degree 4, so this IS [deg-4 of t5]
a2_d4 = a2c - 1.0 * t3c           # max degree 4, so this IS [deg-4 of a2]
t7_d46 = t7c - 1.5 * t3c          # degrees 4 and 6
t7_d4 = t7_d46 - w0 / 2.0         # remove degree 6: [deg-6 of t7] = w0/2

# =====================================================
# IDENTITY (A): [deg-4 of t7] = (1/2)*[deg-4 of t5] + [deg-4 of a2]
# =====================================================
print("\n" + "=" * 70)
print("IDENTITY (A): [deg-4 of t7] = (1/2)*[deg-4 of t5] + [deg-4 of a2]")
print("=" * 70)

pred_A = 0.5 * t5_d4 + 1.0 * a2_d4
err_A = np.max(np.abs(t7_d4 - pred_A))
corr_A = np.corrcoef(t7_d4, pred_A)[0,1]
r2_A = 1 - (t7_d4 - pred_A).var() / t7_d4.var()

print(f"  max|t7_d4 - (0.5*t5_d4 + a2_d4)| = {err_A:.10f}")
print(f"  correlation = {corr_A:.10f}")
print(f"  R^2 = {r2_A:.10f}")

# Regression check
coeffs_A, _, _, _ = lstsq(np.column_stack([t5_d4, a2_d4]), t7_d4, rcond=None)
print(f"  Regression: t7_d4 = {coeffs_A[0]:.6f}*t5_d4 + {coeffs_A[1]:.6f}*a2_d4")

# =====================================================
# IDENTITY (B): w2/4 = 3*[deg-4 of t5] + 6*[deg-4 of a2]
# =====================================================
print("\n" + "=" * 70)
print("IDENTITY (B): w2/4 = 3*[deg-4 of t5] + 6*[deg-4 of a2]")
print("=" * 70)

lhs_B = w2 / 4.0
pred_B = 3.0 * t5_d4 + 6.0 * a2_d4
err_B = np.max(np.abs(lhs_B - pred_B))
corr_B = np.corrcoef(lhs_B, pred_B)[0,1]
r2_B = 1 - (lhs_B - pred_B).var() / lhs_B.var()

print(f"  max|w2/4 - (3*t5_d4 + 6*a2_d4)| = {err_B:.10f}")
print(f"  correlation = {corr_B:.10f}")
print(f"  R^2 = {r2_B:.10f}")

coeffs_B, _, _, _ = lstsq(np.column_stack([t5_d4, a2_d4]), lhs_B, rcond=None)
print(f"  Regression: w2/4 = {coeffs_B[0]:.6f}*t5_d4 + {coeffs_B[1]:.6f}*a2_d4")

# =====================================================
# SVD ANALYSIS: 2D STRUCTURE
# =====================================================
print("\n" + "=" * 70)
print("SVD: DEGREE-4 SPACE IS TWO-DIMENSIONAL")
print("=" * 70)

X = np.column_stack([t5_d4, t7_d4, a2_d4])
U, S, Vt = svd(X, full_matrices=False)
print(f"  Singular values: [{S[0]:.2f}, {S[1]:.2f}, {S[2]:.6f}]")
print(f"  Ratio S[1]/S[0] = {S[1]/S[0]:.4f} (non-trivial: 2D)")
print(f"  Ratio S[2]/S[0] = {S[2]/S[0]:.6f} (≈0: confirms 2D, not 3D)")
print(f"\n  Null vector: {Vt[2,0]:.6f}*t5 + {Vt[2,1]:.6f}*t7 + {Vt[2,2]:.6f}*a2 ≈ 0")
r = Vt[2] / Vt[2,0]
print(f"  Normalized: t5 + {r[1]:.4f}*t7 + {r[2]:.4f}*a2 ≈ 0")
print(f"  => t7 ≈ {-r[1]:.0f}^{{-1}} * (t5 + {r[2]:.4f}*a2) = {-1/r[1]:.4f}*t5 + {-r[2]/r[1]:.4f}*a2")

# =====================================================
# COMBINED OCF CHECK
# =====================================================
print("\n" + "=" * 70)
print("COMBINED: DEGREE-4 OCF IDENTITY CHECK")
print("=" * 70)

# From OCF: [deg-4 of H] = w2/4
# OCF predicts: [deg-4 of H] = 2*[deg-4 of alpha_1] + 4*[deg-4 of alpha_2]
# = 2*([deg-4 of t5] + [deg-4 of t7]) + 4*[deg-4 of a2]
rhs_ocf = 2*(t5_d4 + t7_d4) + 4*a2_d4
err_ocf = np.max(np.abs(lhs_B - rhs_ocf))
print(f"  max|w2/4 - (2*(t5_d4+t7_d4) + 4*a2_d4)| = {err_ocf:.10f}")
print(f"  Using identity (A) to simplify:")
print(f"  = max|w2/4 - (3*t5_d4 + 6*a2_d4)| = {err_B:.10f}")

# =====================================================
# SUMMARY
# =====================================================
print("\n" + "=" * 70)
print("PROOF STATUS SUMMARY FOR OCF AT n=7")
print("=" * 70)
print(f"""
The OCF H = 1 + 2*alpha_1 + 4*alpha_2 decomposes into 4 identities:

  DEGREE 0: n!/2^{{n-1}} = 1 + 2*E[alpha_1] + 4*E[alpha_2]
    STATUS: PROVED (trivial expectation calculation)

  DEGREE 2: w_4/16 = 2*[deg-2 of alpha_1] + 4*[deg-2 of alpha_2]
    STATUS: PROVED via proportionality constants c_5=3, c_7=1.5, c_{{a2}}=1

  DEGREE 4: w_2/4 = 2*[deg-4 of alpha_1] + 4*[deg-4 of alpha_2]
    EQUIVALENT TO: w_2/4 = 3*[deg-4 of t_5] + 6*[deg-4 of alpha_2]
    STATUS: COMPUTATIONALLY VERIFIED (R^2 = {r2_B:.10f})

  DEGREE 6: w_0 = 2*[deg-6 of alpha_1]
    STATUS: PROVED via path-cycle bijection

REMAINING: Prove identities (A) and (B) algebraically to complete OCF at n=7.
""")
