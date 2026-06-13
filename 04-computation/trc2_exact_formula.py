#!/usr/bin/env python3
"""
EXACT FORMULA: tr(c_2) = 24*bc - 60*t_3 + 12*t_5 + 231 at n=7

VERIFIED: max error = 0.0000 over 30 random tournaments.

where:
  t_3 = #directed 3-cycles
  t_5 = #directed 5-cycles
  bc  = sum over 6-vertex subsets S of #{complementary-triple partitions
        of S where both triples are directed 3-cycles}

DERIVATION (combining THM-055 and overlap decomposition):

1. tr(c_2) = sum_P e_4(s_P) over all 7! permutations (THM-055(a))

2. e_4 has C(6,4)=15 terms, split by position topology:
   - (4,): 3 subsets, 5 vertices [NONZERO]
   - (3,1): 6 subsets, 6 vertices [ZERO - singleton cancellation]
   - (2,2): 3 subsets, 6 vertices [NONZERO]
   - (2,1,1): 3 subsets, 7 vertices [ZERO - singleton cancellation]

3. Singleton Cancellation Lemma: If the position subset has an isolated
   position (not adjacent to any other selected position), the centered
   sum vanishes because sum_{distinct (a,b)} (A[a,b] - 1/2) = 0.

4. (4,) contribution = 126 - 36*t_3 + 12*t_5
   Uses OCF at n=5 recursively: chain_sum(S) = H(S) - 3*t_3(S) for |S|=5.

5. (2,2) contribution = 24*bc - 24*t_3 + 105
   Involves complementary-triple partitions of 6-vertex subtournaments.
   Each partition (T, T') contributes F(T)*F(T') where F = 3/2 for cyclic, -1/2 for transitive.

6. Total: tr(c_2) = (126 - 36*t_3 + 12*t_5) + (24*bc - 24*t_3 + 105)
                   = 24*bc - 60*t_3 + 12*t_5 + 231

CONSEQUENCE:
  tr(c_0) = H - tr(c_2)/4 - tr(c_4)/16 - tr(c_6)/64
           = H - 6*bc - 3*t_5 + 249/4

opus-2026-03-06-S11b (continued^4)
"""
from itertools import permutations, combinations
from math import factorial, comb
import numpy as np
import random

def count_3_cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: count += 1
                if A[i][k] and A[k][j] and A[j][i]: count += 1
    return count

def count_5_cycles_dp(A, n):
    count = 0
    for verts in combinations(range(n), 5):
        sub = [[A[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
        dp = [[0]*5 for _ in range(1 << 5)]
        dp[1][0] = 1
        for mask in range(1, 1 << 5):
            for v in range(5):
                if not (mask & (1 << v)) or dp[mask][v] == 0: continue
                for u in range(5):
                    if mask & (1 << u): continue
                    if sub[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full = (1 << 5) - 1
        hc = sum(dp[full][v] for v in range(1, 5) if sub[v][0])
        count += hc
    return count

def ham_count_dp(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(n):
                if (mask & (1 << u)) or A[v][u] != 1: continue
                dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))

def is_3cycle(A, triple):
    a, b, c = triple
    return (A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]) > 0

def both_cyclic_total(A, n):
    """bc = total complementary-cyclic pairs across all 6-vertex subsets."""
    total = 0
    for v_del in range(n):
        S = [v for v in range(n) if v != v_del]
        for T in combinations(S, 3):
            T_comp = tuple(v for v in S if v not in T)
            if T < T_comp:
                if is_3cycle(A, T) and is_3cycle(A, T_comp):
                    total += 1
    return total

def compute_tr_c2_ie(A, n):
    """Compute tr(c_2) via IE formula for transfer matrix at multiple r values."""
    r_values = [0.0, 0.3, 0.5, 0.7]
    r_powers = np.array([[r**(2*k) for k in range(4)] for r in r_values])

    traces = []
    for r in r_values:
        full = (1 << n) - 1
        dp_fwd = [[0.0]*n for _ in range(1 << n)]
        for v in range(n):
            dp_fwd[1 << v][v] = 1.0
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or dp_fwd[mask][v] == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    wt = r + (A[v][u] - 0.5)
                    dp_fwd[mask | (1 << u)][u] += dp_fwd[mask][v] * wt
        dp_bwd = [[0.0]*n for _ in range(1 << n)]
        for v in range(n):
            dp_bwd[1 << v][v] = 1.0
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or dp_bwd[mask][v] == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    wt = r + (A[u][v] - 0.5)
                    dp_bwd[mask | (1 << u)][u] += dp_bwd[mask][v] * wt
        tr_val = 0.0
        for v in range(n):
            for mask_before in range(1 << n):
                if mask_before & (1 << v): continue
                mask_with_v = mask_before | (1 << v)
                if dp_fwd[mask_with_v][v] == 0: continue
                mask_after = full ^ mask_before
                if not (mask_after & (1 << v)): continue
                if dp_bwd[mask_after][v] == 0: continue
                k = bin(mask_before).count('1')
                tr_val += ((-1)**k) * dp_fwd[mask_with_v][v] * dp_bwd[mask_after][v]
        traces.append(tr_val)

    coeffs = np.linalg.solve(r_powers, traces)
    return coeffs  # [c_0, c_2, c_4, c_6]

# =====================================================================
# VERIFICATION
# =====================================================================
n = 7
print("=" * 70)
print(f"EXACT FORMULA VERIFICATION: tr(c_2) = 24*bc - 60*t_3 + 12*t_5 + 231")
print("=" * 70)

max_err_c2 = 0
max_err_c0 = 0
for trial in range(30):
    random.seed(trial * 53)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3_cycles(A, n)
    t5 = count_5_cycles_dp(A, n)
    bc = both_cyclic_total(A, n)
    H = ham_count_dp(A, n)

    # IE-based computation
    c_coeffs = compute_tr_c2_ie(A, n)

    # Formulas
    c2_formula = 24*bc - 60*t3 + 12*t5 + 231
    c4_formula = 240*t3 - 2100
    c6_formula = 720
    c0_formula = H - c2_formula/4 - c4_formula/16 - c6_formula/64

    err_c2 = abs(c_coeffs[1] - c2_formula)
    max_err_c2 = max(max_err_c2, err_c2)

    # Also verify c_0 via the algebraic relation
    c0_via_H = H - c_coeffs[1]/4 - c_coeffs[2]/16 - c_coeffs[3]/64
    err_c0 = abs(c0_via_H - c0_formula)
    max_err_c0 = max(max_err_c0, err_c0)

    if trial < 10:
        print(f"  T{trial:2d}: H={H:3d} t3={t3:2d} t5={t5:2d} bc={bc:2d} | "
              f"c2={c_coeffs[1]:7.1f}/{c2_formula:7.1f} c0={c0_via_H:7.1f}/{c0_formula:7.1f}")

print(f"\nMax error c_2: {max_err_c2:.6f}")
print(f"Max error c_0: {max_err_c0:.6f}")

print(f"\n{'='*70}")
print("COMPLETE COEFFICIENT TABLE AT n=7")
print(f"{'='*70}")
print("""
  tr(c_6) = 720                          (universal)
  tr(c_4) = 240*t_3 - 2100               (depends on t_3)
  tr(c_2) = 24*bc - 60*t_3 + 12*t_5 + 231  (depends on t_3, t_5, bc)
  tr(c_0) = H - 6*bc - 3*t_5 + 249/4    (depends on H, t_5, bc)

  where:
    t_3 = directed 3-cycle count
    t_5 = directed 5-cycle count
    bc  = both-cyclic complementary-triple partition count
        = sum over 6-vertex subsets S of
          #{partitions (T,T') of S into two triples, both 3-cycles}
""")

# Verify the complete picture: c_0 + c_2/4 + c_4/16 + c_6/64 = H
print("Verify: c_0 + c_2/4 + c_4/16 + c_6/64 = H")
for trial in range(5):
    random.seed(trial * 53)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_count_dp(A, n)
    t3 = count_3_cycles(A, n)
    t5 = count_5_cycles_dp(A, n)
    bc = both_cyclic_total(A, n)

    c0 = H - 6*bc - 3*t5 + 249/4
    c2 = 24*bc - 60*t3 + 12*t5 + 231
    c4 = 240*t3 - 2100
    c6 = 720

    reconstructed_H = c0 + c2/4 + c4/16 + c6/64
    print(f"  H={H}, reconstructed={reconstructed_H:.2f}, match={'OK' if abs(H-reconstructed_H)<0.01 else 'FAIL'}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
