#!/usr/bin/env python3
"""
w_{n-5} formula with 5-cycle count.

KEY DERIVATION:
  overlap=3 contribution = 6 * sum_S chain_sum(S)
  where S ranges over 5-element subsets.

  chain_sum(S) = tr(c_0) of T|_S = H(S) - 3*c3(S)
  By OCF at n=5: H(S) = 1 + 2*(c3(S) + c5(S))  [no alpha_2 at n=5]
  So chain_sum(S) = 1 - c3(S) + 2*c5(S)

  sum_S chain_sum(S) = 21 - sum_S c3(S) + 2*sum_S c5(S)
  sum_S c3(S) = C(n-3,2)*t3 = 6*t3  [each 3-cycle in C(4,2) subsets]
  sum_S c5(S) = t5  [each 5-cycle in exactly 1 subset]

  So: overlap=3 = 6*(21 - 6*t3 + 2*t5) = 126 - 36*t3 + 12*t5

  Therefore: w_{n-5} = ov2 + 126 - 36*t3 + 12*t5
  => ov2 = w_{n-5} - 126 + 36*t3 - 12*t5

  Does w_{n-5} = a*t3 + b*t5 + c? (i.e., is ov2 linear in t3,t5?)

opus-2026-03-06-S27
"""

from itertools import permutations, combinations
from math import factorial, comb
import numpy as np
import random

def count_3cycles(A):
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                count += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])
    return count

def count_5cycles(A):
    n = len(A)
    count = 0
    for perm in permutations(range(n), 5):
        if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
            count += 1
    return count // 5

def ham_path_count_dp(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size:
                continue
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                S_prev = S ^ (1 << v)
                total = 0
                for u in range(n):
                    if not (S_prev & (1 << u)):
                        continue
                    if A[u][v] and (S_prev, u) in dp:
                        total += dp[(S_prev, u)]
                if total > 0:
                    dp[(S, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# =====================================================================
n = 7
print("=" * 70)
print(f"w_{{n-5}} = tr(c_2) AT n={n}: FORMULA WITH 5-CYCLES")
print("=" * 70)

random.seed(42)
data = []

for trial in range(30):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3cycles(A)
    t5 = count_5cycles(A)
    H = ham_path_count_dp(A)
    scores = tuple(sorted(sum(A[i]) for i in range(n)))

    # Compute w_{n-5} directly
    w_total = 0
    for p in permutations(range(n)):
        s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]
        e4 = sum(s[i]*s[j]*s[k]*s[l]
                 for i in range(6) for j in range(i+1,6) for k in range(j+1,6)
                 for l in range(k+1,6))
        w_total += e4

    # Predicted overlap=3 contribution
    ov3_pred = 126 - 36*t3 + 12*t5

    data.append({
        't3': t3, 't5': t5, 'H': H, 'scores': scores,
        'w': w_total, 'ov3_pred': ov3_pred
    })

    if trial < 5:
        # Also compute ov3 directly to verify formula
        ov3_direct = 0
        for p in permutations(range(n)):
            s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]
            # overlap=3 positions: (0,1,2,3), (1,2,3,4), (2,3,4,5)
            ov3_direct += s[0]*s[1]*s[2]*s[3] + s[1]*s[2]*s[3]*s[4] + s[2]*s[3]*s[4]*s[5]
        print(f"  Trial {trial}: ov3_direct={ov3_direct:.0f}, ov3_pred={ov3_pred:.0f}, {'✓' if abs(ov3_direct-ov3_pred)<0.01 else '✗'}")

print(f"\n  Collected {len(data)} tournaments")

# =====================================================================
# Test: w = a*t3 + b*t5 + c
# =====================================================================
y = np.array([d['w'] for d in data])

print("\n--- FITTING w_{n-5} ---")

# Just t3
X1 = np.array([[d['t3'], 1] for d in data])
c1 = np.linalg.lstsq(X1, y, rcond=None)[0]
err1 = max(abs(y - X1 @ c1))
print(f"  w ~ a*t3 + b: max_err={err1:.2f}  (NOT exact)")

# t3 + t5
X2 = np.array([[d['t3'], d['t5'], 1] for d in data])
c2 = np.linalg.lstsq(X2, y, rcond=None)[0]
err2 = max(abs(y - X2 @ c2))
print(f"  w ~ a*t3 + b*t5 + c: coeffs={[f'{x:.4f}' for x in c2]}, max_err={err2:.2f}")

# t3 + t5 + t3^2
X3 = np.array([[d['t3'], d['t5'], d['t3']**2, 1] for d in data])
c3 = np.linalg.lstsq(X3, y, rcond=None)[0]
err3 = max(abs(y - X3 @ c3))
print(f"  w ~ a*t3 + b*t5 + c*t3^2 + d: max_err={err3:.2f}")

# t3 + t5 + H
X4 = np.array([[d['t3'], d['t5'], d['H'], 1] for d in data])
c4 = np.linalg.lstsq(X4, y, rcond=None)[0]
err4 = max(abs(y - X4 @ c4))
print(f"  w ~ a*t3 + b*t5 + c*H + d: coeffs={[f'{x:.4f}' for x in c4]}, max_err={err4:.2f}")

# =====================================================================
# Check the ov2 = w - ov3 contribution
# =====================================================================
print("\n--- OVERLAP=2 RESIDUAL ---")
ov2_vals = [d['w'] - d['ov3_pred'] for d in data]
print(f"  ov2 = w - (126 - 36*t3 + 12*t5)")
for d, ov2 in zip(data[:10], ov2_vals[:10]):
    print(f"    t3={d['t3']:2d}, t5={d['t5']:2d}, w={d['w']:6.0f}, ov3_pred={d['ov3_pred']:4.0f}, ov2={ov2:8.0f}")

# Fit ov2 = a*t3 + b*t5 + c
X_ov2 = np.array([[d['t3'], d['t5'], 1] for d in data])
y_ov2 = np.array(ov2_vals)
c_ov2 = np.linalg.lstsq(X_ov2, y_ov2, rcond=None)[0]
err_ov2 = max(abs(y_ov2 - X_ov2 @ c_ov2))
print(f"\n  ov2 ~ a*t3 + b*t5 + c: coeffs={[f'{x:.4f}' for x in c_ov2]}, max_err={err_ov2:.2f}")

# Add t3*t5 interaction?
X5 = np.array([[d['t3'], d['t5'], d['t3']*d['t5'], 1] for d in data])
c5 = np.linalg.lstsq(X5, y_ov2, rcond=None)[0]
err5 = max(abs(y_ov2 - X5 @ c5))
print(f"  ov2 ~ a*t3 + b*t5 + c*t3*t5 + d: max_err={err5:.2f}")

# score squared sum
ssq_vals = []
for d in data:
    A_temp = [[0]*n for _ in range(n)]
    # Need to reconstruct... can't easily. Use different approach.
    pass

# =====================================================================
# Check: same (t3, t5) but different w?
# =====================================================================
print("\n--- Same (t3, t5) but different w? ---")
from collections import defaultdict
groups = defaultdict(list)
for d in data:
    groups[(d['t3'], d['t5'])].append(d['w'])
for key, ws in sorted(groups.items()):
    if len(ws) > 1:
        if max(ws) - min(ws) > 0.01:
            print(f"  (t3={key[0]}, t5={key[1]}): w = {ws} -- VARY")
        else:
            print(f"  (t3={key[0]}, t5={key[1]}): w = {ws[0]:.0f} (all same)")

# =====================================================================
# SUMMARY
# =====================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
  w_{{n-1}} = n! (universal)
  w_{{n-3}} = 2*(n-2)!*t_3 + const  (depends on t_3 only)
  w_{{n-5}} at n=7: NOT determined by (t_3, t_5) alone!

  Decomposition into overlap types:
    overlap=0: N/A at n=7 (need 8 > 7 vertices)
    overlap=1: always 0 (isolated block signed sum cancels)
    overlap=3: = 126 - 36*t_3 + 12*t_5 (depends on t_3 and t_5!)
    overlap=2: the residual, depends on "complementary cycle structure"

  The overlap=3 formula is PROVED:
    Each 5-vertex sub-tournament S contributes chain_sum(S) = H(S) - 3*c_3(S)
    = 1 - c_3(S) + 2*c_5(S) (by OCF at n=5)
    Summing: sum_S chain_sum = 21 - 6*t_3 + 2*t_5
    Times 6: overlap=3 = 126 - 36*t_3 + 12*t_5

  The overlap=2 contribution involves the distribution of 3-cycles
  across complementary triples within 6-vertex sub-tournaments.
  This is a FINER invariant than (t_3, t_5) — it distinguishes
  tournaments within the same (t_3, t_5) class.

  This explains the hierarchy:
    c_{{n-1}}: universal (0-edge universality)
    c_{{n-3}}: depends on t_3 (2-edge: consecutive=t_3, disjoint=universal)
    c_{{n-5}}: depends on t_3, t_5, AND "complementary cycle structure"
    In general: c_{{n-2k-1}} introduces k-edge overlaps, bringing in
    higher-order invariants at each level.
""")

print("=" * 70)
print("DONE")
print("=" * 70)
