#!/usr/bin/env python3
"""
GENERAL n MOMENT FORMULAS via the General Sigma Formula.

The General Sigma Formula gives:
  sigma((c1,...,cm)) = (n - sum_ci - m)! * sum prod H(T[Gj])
over ordered partitions into groups of sizes (c1+1,...,cm+1).

For the FIRST non-trivial moment (j=2):
  m2 = S(2,1)*1!*SIGMA_1 + S(2,2)*2!*SIGMA_2
     = SIGMA_1 + 2*SIGMA_2

  SIGMA_1 = (n-1) * sigma((1,)) = (n-1) * n!/2
  SIGMA_2 = C(n-1,2)*sigma((1,1)) + (n-2)*sigma((2,))
            (non-consecutive pairs + consecutive pairs)

  sigma((1,1)) = (n-4)! * sum_{2-groups} H(G1)*H(G2) = (n-4)! * n(n-1)(n-2)(n-3)/4
  sigma((2,)) = (n-3)! * sum_{3-subsets} H(T[S]) = (n-3)! * [C(n,3) + 2*t3]

  This gives m2 as a function of n and t3 ONLY. Let me derive for general n.

opus-2026-03-06-S28
"""
from itertools import permutations, combinations
from math import factorial, comb
import random

def count_3_cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]: count += 1
                if A[i][k]*A[k][j]*A[j][i]: count += 1
    return count

# =====================================================================
print("=" * 70)
print("GENERAL n: m2 FORMULA")
print("=" * 70)

# sigma((1,)) = n!/2 (each ordered pair appears (n-2)! times, half are edges)
# SIGMA_1 = (n-1) * n!/2

# sigma((1,1)) = (n-4)! * [sum over ordered 4-tuples of A[a,b]*A[c,d]]
# = (n-4)! * n(n-1)(n-2)(n-3)/4  (pair-partition universality)
# Number of non-consecutive pairs of positions: C(n-1,2) - (n-2) = C(n-2,2)

# sigma((2,)) = (n-3)! * [sum over ordered triples of A[a,b]*A[b,c]]
# = (n-3)! * [C(n,3) + 2*t3]
# Number of consecutive pairs: n-2

# SIGMA_2 = C(n-2,2) * (n-4)! * n(n-1)(n-2)(n-3)/4  +  (n-2) * (n-3)! * (C(n,3)+2*t3)

# Let's simplify:
# Term 1: C(n-2,2) * (n-4)! * n(n-1)(n-2)(n-3)/4
# = [(n-2)(n-3)/2] * (n-4)! * [n(n-1)(n-2)(n-3)/4]
# = (n-2)(n-3) * (n-4)! * n(n-1)(n-2)(n-3) / 8
# = n(n-1)(n-2)^2(n-3)^2 * (n-4)! / 8

# Term 2: (n-2) * (n-3)! * (C(n,3) + 2*t3)
# = (n-2) * (n-3)! * [n(n-1)(n-2)/6 + 2*t3]
# = (n-2)*(n-3)! * n(n-1)(n-2)/6 + 2*(n-2)*(n-3)!*t3

# m2 = SIGMA_1 + 2*SIGMA_2
# = (n-1)*n!/2 + 2*SIGMA_2

# This is getting messy. Let me just compute numerically for various n and verify.

for n in [3, 5, 7, 9]:
    print(f"\n  n={n}:")

    # SIGMA_1
    SIGMA_1 = (n-1) * factorial(n) // 2

    # Number of non-consecutive position pairs
    num_nonconsec = comb(n-1, 2) - (n-2)
    # Number of consecutive position pairs
    num_consec = n - 2

    # sigma formulas
    sig_11 = factorial(n-4) * n*(n-1)*(n-2)*(n-3) // 4  if n >= 4 else 0
    # sigma((2,)) = (n-3)! * (C(n,3) + 2*t3)

    # SIGMA_2(t3) = num_nonconsec * sig_11 + num_consec * (n-3)! * (C(n,3) + 2*t3)
    # = num_nonconsec * sig_11 + num_consec * (n-3)! * C(n,3) + 2*num_consec*(n-3)!*t3

    const_SIGMA_2 = num_nonconsec * sig_11 + num_consec * factorial(max(n-3,0)) * comb(n,3)
    t3_coeff_SIGMA_2 = 2 * num_consec * factorial(max(n-3,0))

    # m2 = SIGMA_1 + 2*SIGMA_2
    m2_const = SIGMA_1 + 2 * const_SIGMA_2
    m2_t3 = 2 * t3_coeff_SIGMA_2

    print(f"    SIGMA_1 = {SIGMA_1}")
    print(f"    sig((1,1)) = {sig_11}, count = {num_nonconsec}")
    print(f"    sig((2,)) coeff: const={(n-2)*factorial(max(n-3,0))*comb(n,3)}, "
          f"t3_coeff={2*(n-2)*factorial(max(n-3,0))}")
    print(f"    SIGMA_2 = {const_SIGMA_2} + {t3_coeff_SIGMA_2}*t3")
    print(f"    m2 = {m2_const} + {m2_t3}*t3")

    # Verify with random tournaments
    if n <= 9:
        for trial in range(3):
            random.seed(n*1000 + trial)
            A = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
            t3 = count_3_cycles(A, n)
            m2_actual = sum(
                sum(1 for i in range(n-1) if A[p[i]][p[i+1]])**2
                for p in permutations(range(n))
            )
            m2_pred = m2_const + m2_t3 * t3
            print(f"    Trial {trial}: t3={t3}, m2={m2_actual}, pred={m2_pred}, "
                  f"{'OK' if m2_actual == m2_pred else 'FAIL'}")

# =====================================================================
print(f"\n{'='*70}")
print("GENERAL n: m3 FORMULA")
print(f"{'='*70}")

# m3 = S(3,1)*SIGMA_1 + S(3,2)*2*SIGMA_2 + S(3,3)*6*SIGMA_3
# S(3,1)=1, S(3,2)=3, S(3,3)=1
# = SIGMA_1 + 6*SIGMA_2 + 6*SIGMA_3

# SIGMA_3 involves patterns:
# (1,1,1): C(n-2,2) - ... need to count exactly
# (2,1): mix of consecutive pair + isolated
# (3,): three consecutive

# Number of each pattern in size-3 position subsets of [0,n-2]:
# (1,1,1): no two adjacent. This is C(n-3, 3) at n-1 positions... tricky.
# Let me just compute for specific n.

def count_position_patterns(n, k):
    from collections import Counter
    patterns = Counter()
    for S in combinations(range(n-1), k):
        pos = sorted(S)
        comps = []
        comp = [pos[0]]
        for i in range(1, len(pos)):
            if pos[i] == comp[-1] + 1:
                comp.append(pos[i])
            else:
                comps.append(len(comp))
                comp = [pos[i]]
        comps.append(len(comp))
        pat = tuple(sorted(comps, reverse=True))
        patterns[pat] += 1
    return dict(patterns)

for n in [5, 7, 9, 11]:
    print(f"\n  n={n}: position patterns at k=3:")
    pats = count_position_patterns(n, 3)
    for pat, cnt in sorted(pats.items()):
        print(f"    {pat}: {cnt}")

    # SIGMA_3 = count(1,1,1)*sigma((1,1,1)) + count(2,1)*sigma((2,1)) + count(3,)*sigma((3,))

    # sigma((1,1,1)): 3 isolated positions using 6 vertices.
    # = (n-6)! * sum over ordered 6-tuples of A[a,b]*A[c,d]*A[e,f] (3 pair products)
    # = (n-6)! * 6!/8 * C(n,6) by pair-partition universality (k=3 case)
    # 6!/8 = 720/8 = 90. C(n,6) = ...
    sig_111 = factorial(max(n-6,0)) * 90 * comb(n, 6)

    # sigma((2,1)): consecutive pair + isolated, using 5 vertices.
    # = (n-5)! * sum_{5-subsets} sum_{(3,2) partitions} dp_2(3-set) * 1
    # = (n-5)! * C(n-3, 2) * [C(n,3) + 2*t3]  ... hmm, derived earlier
    # Actually: = (n-5)! * sum_{3-subsets T} dp_2(T) * C(n-3, 5-3)
    # = (n-5)! * C(n-3, 2) * [C(n,3) + 2*t3] / C(n,3)... no that's wrong.
    # From our n=7 derivation: sigma((2,1)) = (n-5)! * C(n-3, 2) * (C(n,3) + 2*t3)
    # At n=7: (2)! * C(4,2) * (35+2*t3) = 2*6*(35+2*t3) = 12*(35+2*t3) = 420+24*t3 ✓

    sig_21_const = factorial(max(n-5,0)) * comb(n-3, 2) * comb(n, 3)
    sig_21_t3 = factorial(max(n-5,0)) * comb(n-3, 2) * 2

    # sigma((3,)): 3 consecutive = 4 vertices, directed 3-path = H(4-set)
    # = (n-4)! * sum_{4-subsets} H(T[S]) = (n-4)! * [C(n,4) + 2*C(n-3,1)*t3]
    # = (n-4)! * [C(n,4) + 2*(n-3)*t3]
    sig_3_const = factorial(max(n-4,0)) * comb(n, 4)
    sig_3_t3 = factorial(max(n-4,0)) * 2 * (n-3)

    count_111 = pats.get((1,1,1), 0)
    count_21 = pats.get((2,1), 0)
    count_3 = pats.get((3,), 0)

    SIGMA_3_const = count_111*sig_111 + count_21*sig_21_const + count_3*sig_3_const
    SIGMA_3_t3 = count_21*sig_21_t3 + count_3*sig_3_t3

    # m3 = SIGMA_1 + 6*SIGMA_2 + 6*SIGMA_3
    m3_const = SIGMA_1 + 6*const_SIGMA_2 + 6*SIGMA_3_const  # Hmm, using vars from earlier...

    # Let me just compute directly for each n
    m3_const_n = SIGMA_1 + 6*(num_nonconsec * sig_11 + num_consec * factorial(max(n-3,0)) * comb(n,3)) + 6*SIGMA_3_const
    m3_t3_n = 6*t3_coeff_SIGMA_2 + 6*SIGMA_3_t3

    # Wait, I need to recompute SIGMA_1, SIGMA_2 for this n
    SIGMA_1_n = (n-1) * factorial(n) // 2
    num_nonconsec_n = comb(n-1, 2) - (n-2)
    num_consec_n = n - 2
    sig_11_n = factorial(max(n-4,0)) * n*(n-1)*(n-2)*(n-3) // 4 if n >= 4 else 0
    SIGMA_2_const_n = num_nonconsec_n * sig_11_n + num_consec_n * factorial(max(n-3,0)) * comb(n,3)
    SIGMA_2_t3_n = 2 * num_consec_n * factorial(max(n-3,0))

    m3_const_n = SIGMA_1_n + 6*SIGMA_2_const_n + 6*SIGMA_3_const
    m3_t3_n = 6*SIGMA_2_t3_n + 6*SIGMA_3_t3

    print(f"    m3 = {m3_const_n} + {m3_t3_n}*t3")

    # Verify
    if n <= 7:
        for trial in range(2):
            random.seed(n*1000 + trial)
            A = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
            t3 = count_3_cycles(A, n)
            m3_actual = sum(
                sum(1 for i in range(n-1) if A[p[i]][p[i+1]])**3
                for p in permutations(range(n))
            )
            m3_pred = m3_const_n + m3_t3_n * t3
            print(f"    Trial {trial}: t3={t3}, m3={m3_actual}, pred={m3_pred}, "
                  f"{'OK' if m3_actual == m3_pred else 'FAIL'}")

# =====================================================================
# Express coefficients in closed form
# =====================================================================
print(f"\n{'='*70}")
print("CLOSED FORM: m2 for general odd n")
print(f"{'='*70}")

for n in [3, 5, 7, 9, 11, 13]:
    SIGMA_1 = (n-1) * factorial(n) // 2
    nonconsec = comb(n-1, 2) - (n-2)
    consec = n - 2
    sig_11 = factorial(n-4) * n*(n-1)*(n-2)*(n-3) // 4 if n >= 4 else 0
    SIGMA_2_c = nonconsec * sig_11 + consec * factorial(n-3) * comb(n,3) if n >= 4 else 0
    SIGMA_2_t3 = 2 * consec * factorial(n-3) if n >= 3 else 0

    m2_c = SIGMA_1 + 2*SIGMA_2_c
    m2_t3 = 2*SIGMA_2_t3

    # Expected from pattern: m2_t3 = 4*(n-2)*(n-3)!
    expected_t3 = 4 * (n-2) * factorial(n-3)
    # m2_const: let's see if it has a nice form
    # m2_const / n! = ?
    ratio = m2_c / factorial(n)

    print(f"  n={n:2d}: m2 = {m2_c} + {m2_t3}*t3")
    print(f"         t3_coeff = {m2_t3} = 4*(n-2)*(n-3)! = {expected_t3} {'OK' if m2_t3 == expected_t3 else 'FAIL'}")
    print(f"         const/n! = {ratio:.6f}")

# The t3 coefficient in m2 is 4*(n-2)*(n-3)! = 4*(n-2)!
print(f"\n  RESULT: m2(t3) = 4*(n-2)! * t3 + universal_const(n)")
print(f"  => tr(c_{{n-3}}) via Newton: involves m2 only at this level")
print(f"  => tr(c_{{n-3}}) = 2*(n-2)!*t3 - (n-2)!*C(n,3)/2  [THM-054]")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
