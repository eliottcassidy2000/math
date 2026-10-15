#!/usr/bin/env python3
"""Exact check of every step of the THM-4655 proof for strong transitive-plus-one tournaments (w_1 = 1, w_m = 0), m <= MMAX:
  K = {0, m} u {k : w_k = 0, w_(k+1) = 1};  C_1 = sum_(alpha: w_alpha = 0) 2^(alpha-2);  C_(>=3) = sum_(beta>=3, w_beta=1)
  sum_(alpha>beta, w_alpha=0) 2^(alpha-beta-1);  N0 = #{alpha >= 4 : w_alpha = 0}.
  (i) H = |K| + sum over crossed pairs; brute-force H and in2 = #{valid A containing u_2} agree with the case formulas
      w_2 = 1: in2 = |K| - 1 + C_(>=3) + C_1/2,  H = |K| + 3 C_1/2 + C_(>=3)
      w_2 = 0: in2 = |K| + C_(>=3) + (C_1 - 1)/2,  H = |K| + C_1 + C_(>=3);
  (ii) |K| <= 2 + N0 (injection of 01-boundaries into later zeros); (iii) 2 C_(>=3) <= C_1 - 1 - 2 N0 (w_2 = 0) and
      2 C_(>=3) <= C_1 - 2 N0 - [w_3 = 0]*... checked as 2 C_(>=3) <= sum_(alpha>=4, w_alpha=0) (2^(alpha-2) - 2);
  (iv) the final inequalities 1/3 < in2/H <= 2/3, with equality 2/3 exactly at w = 1 0 1^(m-3) 0."""
import sys
from itertools import product

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
viol = 0; nwords = 0
for m in range(3, MMAX + 1):
    for mid in product((0, 1), repeat=m - 2):
        w = (1,) + mid + (0,)                     # w[0] = w_1
        W = lambda i: w[i - 1]                    # 1-based
        nwords += 1
        # brute force
        H = in2 = 0
        for A in range(1 << m):
            B = ((1 << m) - 1) ^ A
            if A and W(A.bit_length()) != 0: continue
            if B and W((B & -B).bit_length()) != 1: continue
            H += 1
            if A >> 1 & 1: in2 += 1
        K = 2 + sum(1 for k in range(1, m) if W(k) == 0 and W(k + 1) == 1)
        C1 = sum(2 ** (a - 2) for a in range(2, m + 1) if W(a) == 0)
        C3 = sum(2 ** (a - b - 1) for b in range(3, m + 1) if W(b) == 1 for a in range(b + 1, m + 1) if W(a) == 0)
        N0 = sum(1 for a in range(4, m + 1) if W(a) == 0)
        if W(2) == 1:
            ok_i = (2 * in2 == 2 * (K - 1) + 2 * C3 + C1) and (2 * H == 2 * K + 3 * C1 + 2 * C3)
        else:
            ok_i = (2 * in2 == 2 * K + 2 * C3 + C1 - 1) and (H == K + C1 + C3)
        ok_ii = K <= 2 + N0
        ok_iii = 2 * C3 <= sum(2 ** (a - 2) - 2 for a in range(4, m + 1) if W(a) == 0)
        ok_iv = 3 * in2 > H and 3 * in2 <= 2 * H
        eq = (3 * in2 == 2 * H)
        eq_pred = (m >= 4 and w == (1, 0) + (1,) * (m - 3) + (0,))
        if not (ok_i and ok_ii and ok_iii and ok_iv and eq == eq_pred):
            viol += 1
            if viol < 5: print("FAIL", m, w, ok_i, ok_ii, ok_iii, ok_iv, eq, eq_pred, H, in2, K, C1, C3, N0)
print(f"checked {nwords} strong words, m = 3..{MMAX}: failures {viol}")
