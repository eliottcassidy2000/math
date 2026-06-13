#!/usr/bin/env python3
"""
BRIDGE between even-r polynomial and OCF.

Even-r: H = sum_k tr(c_{2k}) / 4^k
OCF:    H = sum_k alpha_k * 2^k  where alpha_k = #{indep sets of size k in Omega(T)}

Both are polynomial evaluations giving H. What's the connection?

At n=5:
  Even-r: H = tr(c_0) + tr(c_2)/4 + tr(c_4)/16
         = tr(c_0) + (12*t3-30)/4 + 120/16
         = tr(c_0) + 3*t3 - 7.5 + 7.5
         = tr(c_0) + 3*t3
  OCF: H = 1 + 2*(t3+t5)  [alpha_2=0 at n=5]
  => tr(c_0) = 1 - t3 + 2*t5

At n=7:
  Even-r: H = tr(c_0) + tr(c_2)/4 + (240*t3-2100)/16 + 5040/64
         = tr(c_0) + tr(c_2)/4 + 15*t3 - 131.25 + 78.75
         = tr(c_0) + tr(c_2)/4 + 15*t3 - 52.5
  OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3

  So: tr(c_0) + tr(c_2)/4 = H - 15*t3 + 52.5
                            = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 - 15*t3 + 52.5
                            = 53.5 + 2*alpha_1 - 15*t3 + 4*alpha_2 + 8*alpha_3

  Since alpha_1 = t3 + t5 + t7 (number of odd cycles):
  = 53.5 + 2*(t3+t5+t7) - 15*t3 + 4*alpha_2 + 8*alpha_3
  = 53.5 - 13*t3 + 2*t5 + 2*t7 + 4*alpha_2 + 8*alpha_3

  This is messy. But at n=5 it was clean: tr(c_0) = 1 - t3 + 2*t5.

  Let me just compute these quantities and look for patterns.

opus-2026-03-06-S27
"""

from itertools import permutations, combinations
from math import factorial, comb
import numpy as np
import random

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

def count_odd_cycles(A, max_len=None):
    """Count directed odd cycles of each length."""
    n = len(A)
    if max_len is None:
        max_len = n
    counts = {}
    for L in range(3, max_len+1, 2):
        count = 0
        for perm in permutations(range(n), L):
            if all(A[perm[i]][perm[(i+1)%L]] for i in range(L)):
                count += 1
        counts[L] = count // L
    return counts

def count_indep_sets(A):
    """Count independent sets in Omega(T) at each size.
    Omega vertices = directed odd cycles, edge = share vertex."""
    n = len(A)
    # Enumerate all directed odd cycles
    cycles = []
    for L in range(3, n+1, 2):
        for perm in permutations(range(n), L):
            if perm[0] == min(perm) and all(A[perm[i]][perm[(i+1)%L]] for i in range(L)):
                cycles.append(frozenset(perm))

    # Build conflict graph (adjacency in Omega)
    m = len(cycles)
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycles[i] & cycles[j]:
                adj[i][j] = adj[j][i] = True

    # Count independent sets by size (brute force for small m)
    alpha = [0] * (m+1)
    for mask in range(1 << m):
        bits = [i for i in range(m) if mask & (1 << i)]
        independent = True
        for a in range(len(bits)):
            for b in range(a+1, len(bits)):
                if adj[bits[a]][bits[b]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            alpha[len(bits)] += 1

    return alpha

# =====================================================================
print("=" * 70)
print("EVEN-R POLYNOMIAL vs OCF")
print("=" * 70)

# n=5 detailed analysis
n = 5
print(f"\n  n={n}:")
random.seed(55)
for trial in range(8):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_path_count_dp(A)
    cyc = count_odd_cycles(A)
    t3 = cyc.get(3, 0)
    t5 = cyc.get(5, 0)

    # Even-r: tr(c_0) = H - 3*t3
    tr_c0 = H - 3*t3

    # OCF prediction: tr(c_0) should = 1 - t3 + 2*t5
    tr_c0_pred = 1 - t3 + 2*t5

    # Independent sets
    if trial < 4:
        alpha = count_indep_sets(A)
        H_ocf = sum(alpha[k] * 2**k for k in range(len(alpha)))
        print(f"    Trial {trial}: H={H}, t3={t3}, t5={t5}, "
              f"tr(c_0)={tr_c0}, pred={tr_c0_pred}, "
              f"alpha={alpha[:4]}, H_ocf={H_ocf}")
    else:
        print(f"    Trial {trial}: H={H}, t3={t3}, t5={t5}, "
              f"tr(c_0)={tr_c0}, pred={tr_c0_pred}, "
              f"match={'✓' if tr_c0 == tr_c0_pred else '✗'}")

# =====================================================================
# n=7: what does tr(c_0) equal in terms of cycle counts?
# =====================================================================
n = 7
print(f"\n  n={n}: tr(c_0) decomposition")
print(f"  tr(c_0) = H - tr(c_2)/4 - (240*t3-2100)/16 - 5040/64")

# For this we need tr(c_2), which varies.
# For regular tournaments (t3=14, tr(c_2)=63):
print(f"\n  REGULAR tournaments at n=7:")
print(f"    tr(c_0) = H - 63/4 - 1260/16 - 5040/64")
print(f"           = H - 15.75 - 78.75 - 78.75")
print(f"           = H - 173.25")
print(f"    For H=189 (Paley): tr(c_0) = 15.75")
print(f"    For H=175: tr(c_0) = 1.75")
print(f"    For H=171: tr(c_0) = -2.25")

# What does OCF give? H = 1 + 2*alpha_1 + 4*alpha_2 + ...
# tr(c_0) = H - 173.25
# = 1 + 2*alpha_1 + 4*alpha_2 + ... - 173.25
# = -172.25 + 2*alpha_1 + 4*alpha_2 + ...

# For Paley T_7: let me count cycles
print(f"\n  Counting cycles for Paley T_7...")
# Paley T_7: gen = {1,2,4} (QRs mod 7)
A_paley = [[0]*7 for _ in range(7)]
QR = {1, 2, 4}  # quadratic residues mod 7
for i in range(7):
    for j in range(7):
        if i != j and (j-i) % 7 in QR:
            A_paley[i][j] = 1

H_paley = ham_path_count_dp(A_paley)
cyc_paley = count_odd_cycles(A_paley)
print(f"    H = {H_paley}")
print(f"    Cycles: {cyc_paley}")

# OCF
alpha_paley = count_indep_sets(A_paley)
H_ocf = sum(alpha_paley[k] * 2**k for k in range(len(alpha_paley)))
print(f"    alpha = {alpha_paley[:6]}")
print(f"    H via OCF = {H_ocf}")

# Bridge formula
tr_c0_paley = H_paley - 173.25
print(f"    tr(c_0) = {tr_c0_paley}")
print(f"    OCF: 1 + 2*{cyc_paley.get(3,0)+cyc_paley.get(5,0)+cyc_paley.get(7,0)} + ...")

# =====================================================================
# KEY INSIGHT: tr(c_0) = w_0 = sum_P prod(s_e)
# = sum_P prod(T(e)-1/2) = (1/2^{n-1}) * sum_P (-1)^{b_P}
# = (1/2^{n-1}) * [sum_P with b_P even - sum_P with b_P odd]
# =====================================================================
print("\n" + "=" * 70)
print("tr(c_0) AS SIGNED PATH COUNT")
print("=" * 70)

print(f"""
  tr(c_0) = sum_P prod_e (T(e) - 1/2)
           = (1/2^(n-1)) * sum_P (-1)^(backward edges)
           = (1/2^(n-1)) * [E(even backward) - E(odd backward)]

  At odd n (n-1 even): (-1)^b = (-1)^(n-1-f) = (-1)^f
  (since n-1 is even, so (-1)^(n-1) = 1)

  So tr(c_0) = (1/2^(n-1)) * sum_P (-1)^(f_P)
             = (1/2^(n-1)) * [paths with even forward - paths with odd forward]
""")

for n in [3, 5, 7]:
    random.seed(n * 33)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    # Direct tr(c_0)
    w0 = 0
    even_count = 0
    odd_count = 0
    for p in permutations(range(n)):
        f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        prod_s = 1.0
        for i in range(n-1):
            prod_s *= (A[p[i]][p[i+1]] - 0.5)
        w0 += prod_s
        if f % 2 == 0:
            even_count += 1
        else:
            odd_count += 1

    signed = (even_count - odd_count) / 2**(n-1)
    print(f"  n={n}: w0={w0:.4f}, (even-odd)/2^(n-1)={signed:.4f}, "
          f"even={even_count}, odd={odd_count}")

# =====================================================================
# AMAZING CONNECTION: at r=0, W(0) = tr(c_0) = (1/2^{n-1}) * sum (-1)^f
# And sum (-1)^f = det(J - 2A) where J is all-ones and A is adjacency?
# =====================================================================
print("\n" + "=" * 70)
print("DETERMINANT CONNECTION?")
print("=" * 70)

print(f"  tr(c_0) * 2^(n-1) = sum_P (-1)^f_P")
print(f"  This looks like a permanent with signs = determinant!")
print(f"  sum_P (-1)^f_P = sum_perm prod_{i} T(p_i, p_{i+1})^f * (-1)^f")
print(f"  Hmm, not exactly a determinant...")
print(f"")
print(f"  But: sum_P prod (2*T(v_i,v_{i+1}) - 1) = sum_P (-1)^(b_P) * 1^(f_P) = sum_P (-1)^b_P")
print(f"  Let B[i,j] = 2*T(i,j) - 1 in {{-1, +1}} for i≠j, B[i,i]=0.")
print(f"  Then sum_P prod B[v_i, v_{{i+1}}] might relate to a matrix permanent/determinant.")

# Actually, W(0) = sum_P prod(T-1/2) = sum_P prod(s_e) where s_e in {-1/2, +1/2}
# = (1/2^{n-1}) * sum_P prod(2*T-1) where 2*T-1 in {-1,+1}
# = (1/2^{n-1}) * sum_P prod B_e
# This is the weighted permanent of the skew matrix B with Hamiltonian path structure.

# For a CIRCULANT tournament: let me check if tr(c_0) has a nice form.
print(f"\n  CIRCULANT at n=7 (Paley):")
B_paley = [[2*A_paley[i][j]-1 if i!=j else 0 for j in range(7)] for i in range(7)]
w0_paley = 0
for p in permutations(range(7)):
    prod_b = 1
    for i in range(6):
        prod_b *= B_paley[p[i]][p[i+1]]
    w0_paley += prod_b
tr_c0_from_B = w0_paley / 2**6
print(f"    sum_P prod B = {w0_paley}")
print(f"    tr(c_0) = {w0_paley}/64 = {tr_c0_from_B:.4f}")
print(f"    Expected: 15.75")
print(f"    Match: {'✓' if abs(tr_c0_from_B - 15.75) < 0.01 else '✗'}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
