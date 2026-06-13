#!/usr/bin/env python3
"""
ALGEBRAIC PROOF of tr(c_{n-3}) = 2*(n-2)!*t_3 - (n-2)!*C(n,3)/2

Key identity: W(r) = tr(M(r)) at odd n, where W(r) = sum_P prod_e (r + s_e).

Step 1: w_{n-3} = sum_P e_2(s_P) where e_2 = sum_{i<j} s_i*s_j
Step 2: e_2(s_P) = -(1/2)*f_P*b_P + C(n-1,2)/4
Step 3: sum_P f_P*b_P involves edge-pair counting in permutations
Step 4: Connect to 3-cycles via sum s_b^2

This script verifies each step computationally at n=5 and n=7.

opus-2026-03-06-S27
"""

from itertools import permutations, combinations
import numpy as np
from math import factorial, comb

# =====================================================================
# STEP 1: Verify e_2(s_P) = -(1/2)*f_P*b_P + C(n-1,2)/4
# =====================================================================
def verify_e2_formula(n):
    """For random tournament, check that e_2(s_P) = -(1/2)*f*b + C(n-1,2)/4."""
    import random
    random.seed(n * 137)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    max_err = 0
    for p in permutations(range(n)):
        s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]
        e2 = sum(s[i]*s[j] for i in range(n-1) for j in range(i+1, n-1))
        f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        b = n - 1 - f
        predicted = -0.5 * f * b + comb(n-1, 2) / 4
        max_err = max(max_err, abs(e2 - predicted))

    return max_err

# =====================================================================
# STEP 2: Compute sum_P f_P*b_P by edge-pair counting
# =====================================================================
def compute_sum_fb(A):
    """Direct computation of sum_P f_P * b_P."""
    n = len(A)
    total = 0
    for p in permutations(range(n)):
        f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        b = n - 1 - f
        total += f * b
    return total

def compute_sum_f2(A):
    """Direct computation of sum_P f_P^2."""
    n = len(A)
    total = 0
    for p in permutations(range(n)):
        f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        total += f * f
    return total

def score_sum_sq(A):
    n = len(A)
    scores = [sum(A[i]) for i in range(n)]
    return sum(s*s for s in scores)

def count_3cycles(A):
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                count += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])
    return count

# =====================================================================
# STEP 3: Derive the formula for sum_P f_P^2
# =====================================================================
# f_P = sum_{i=0}^{n-2} T(v_i, v_{i+1})
# f_P^2 = sum_i T_i + 2*sum_{i<j} T_i*T_j   (since T_i^2 = T_i)
#        = f_P + 2*[Sigma_consec + Sigma_non_consec]
#
# Sigma_consec = sum_P sum_i T(v_i,v_{i+1})*T(v_{i+1},v_{i+2})
#              = (n-2)! * sum_{b} s_b * indeg(b)
#              = (n-2)! * [(n-1)*C(n,2) - sum s_b^2]
#
# Sigma_non_consec: For positions (i,j) with j > i+1:
#   All four vertices are distinct (permutation).
#   For each of the C(n-1,2)-(n-2) position pairs:
#     sum_P T(v_i,v_{i+1})*T(v_j,v_{j+1}) = (n-4)! * sum_{disjoint (a,b),(c,d)} T(a,b)*T(c,d)
#
# sum_{disjoint} T(a,b)*T(c,d) = [sum T(a,b)]^2 - sum_{overlapping} T(a,b)*T(c,d)
#                               = C(n,2)^2 - Sigma_overlap

def compute_sigma_overlap(A):
    """Compute sum over overlapping directed edge pairs T(a,b)*T(c,d) where {a,b}∩{c,d}≠∅."""
    n = len(A)
    scores = [sum(A[i]) for i in range(n)]
    total = 0

    # Case 1: a=c (share first vertex)
    # sum_a sum_{b≠a} sum_{d≠a, d≠b} T(a,b)*T(a,d)
    for a in range(n):
        for b in range(n):
            if b == a: continue
            for d in range(n):
                if d == a or d == b: continue
                total += A[a][b] * A[a][d]

    # Case 2: a=d (first of one = second of other)
    for a in range(n):
        for b in range(n):
            if b == a: continue
            for c in range(n):
                if c == a or c == b: continue
                total += A[a][b] * A[c][a]

    # Case 3: b=c (second of one = first of other)
    for b in range(n):
        for a in range(n):
            if a == b: continue
            for d in range(n):
                if d == b or d == a: continue
                total += A[a][b] * A[b][d]

    # Case 4: b=d (share second vertex)
    for b in range(n):
        for a in range(n):
            if a == b: continue
            for c in range(n):
                if c == b or c == a: continue
                total += A[a][b] * A[c][b]

    return total

def compute_sigma_overlap_formula(A):
    """Same thing using closed forms."""
    n = len(A)
    scores = [sum(A[i]) for i in range(n)]
    indegs = [n - 1 - s for s in scores]

    # Case 1: a=c => sum_a s_a*(s_a-1) = sum s_a^2 - sum s_a = sum s_a^2 - C(n,2)
    case1 = sum(s*(s-1) for s in scores)

    # Case 2: a=d => sum_a s_a * (n-1-s_a - ...) — more complex
    # T(a,b)*T(c,a): b out-neighbor of a, c in-neighbor of a, c≠b
    # = sum_a sum_{b: a→b} sum_{c: c→a, c≠b} 1
    # = sum_a s_a * (indeg(a) - T(b...)) hmm
    # Actually: for each a, sum over b≠a of T(a,b) * sum over c≠a,c≠b of T(c,a)
    # = sum_a sum_{b≠a} T(a,b) * [indeg(a) - (1-T(a,b))]
    # Wait: T(c,a) where c≠b. If b is an in-neighbor of a (T(b,a)=1), then
    # removing c=b removes 1 from indeg. If b is out-neighbor (T(a,b)=1, T(b,a)=0),
    # removing c=b removes 0.
    # sum_{c≠a,c≠b} T(c,a) = indeg(a) - T(b,a) = (n-1-s_a) - (1 - T(a,b))
    # So: sum_a sum_{b≠a} T(a,b) * [(n-1-s_a) - 1 + T(a,b)]
    #   = sum_a sum_{b≠a} T(a,b) * [n-2-s_a+T(a,b)]
    #   = sum_a [s_a*(n-2-s_a) + sum_{b≠a} T(a,b)^2]
    #   = sum_a [s_a*(n-2-s_a) + s_a]
    #   = sum_a [s_a*(n-1-s_a)]
    case2 = sum(s*(n-1-s) for s in scores)

    # Case 3: b=c => sum_b indeg(b) * (s_b - 0) where we exclude a from d's range
    # Actually computed above: = sum_b (n-1-s_b)*s_b
    # Wait need to be more careful. T(a,b)*T(b,d), a≠b, d≠b, d≠a.
    # = sum_b sum_{a≠b} T(a,b) * [s_b - T(b,a)]
    # = sum_b sum_{a≠b} T(a,b) * [s_b - (1-T(a,b))]
    # = sum_b sum_{a≠b} [T(a,b)*s_b - T(a,b) + T(a,b)^2]
    # = sum_b [indeg(b)*s_b - indeg(b) + indeg(b)]
    # = sum_b indeg(b)*s_b = sum_b (n-1-s_b)*s_b
    case3 = sum((n-1-s)*s for s in scores)

    # Case 4: b=d => sum_b indeg(b)*(indeg(b)-1) = sum (n-1-s)^2 - (n-1-s)
    #              = sum (n-1-s)*((n-1-s)-1) = sum indeg*(indeg-1)
    case4 = sum((n-1-s)*((n-1-s)-1) for s in scores)

    return case1, case2, case3, case4, case1+case2+case3+case4

# =====================================================================
# VERIFICATION
# =====================================================================
import random

print("=" * 70)
print("ALGEBRAIC PROOF OF tr(c_{n-3}) FORMULA")
print("=" * 70)

# Step 1: e_2 = -(1/2)*f*b + C(n-1,2)/4
for n in [4, 5, 6, 7]:
    err = verify_e2_formula(n)
    print(f"  Step 1 at n={n}: e_2 = -(1/2)*f*b + C(n-1,2)/4, max_err={err:.2e} {'✓' if err < 1e-10 else '✗'}")

# Step 2: Verify sum_P f^2 formula
print(f"\n  Step 2: Verify sum_P f_P^2 formula")
for n in [4, 5]:
    random.seed(n * 99)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    sum_f = 0
    sum_f2 = 0
    sum_fb = 0
    for p in permutations(range(n)):
        f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        sum_f += f
        sum_f2 += f*f
        sum_fb += f*(n-1-f)

    scores = [sum(A[i]) for i in range(n)]
    ssq = sum(s*s for s in scores)

    # Theoretical sum_f = C(n,2) * (n-1)!
    theo_sum_f = comb(n,2) * factorial(n-1)

    # Sigma_consec = (n-2)! * [(n-1)*C(n,2) - ssq]
    sig_consec = factorial(n-2) * ((n-1)*comb(n,2) - ssq)

    # Sigma_overlap (direct)
    sig_overlap_direct = compute_sigma_overlap(A)
    c1, c2, c3, c4, sig_overlap_formula = compute_sigma_overlap_formula(A)

    # sum_{disjoint} T(a,b)*T(c,d)
    sig_disjoint = comb(n,2)**2 - sig_overlap_direct

    # Number of non-consecutive position pairs
    num_non_consec = comb(n-1, 2) - (n-2)

    # Sigma_non_consec = num_non_consec * (n-4)! * sig_disjoint
    sig_non_consec = num_non_consec * factorial(n-4) * sig_disjoint if n >= 4 else 0

    # Theoretical sum_f2 = sum_f + 2*(sig_consec + sig_non_consec)
    theo_sum_f2 = theo_sum_f + 2*(sig_consec + sig_non_consec)

    print(f"\n  n={n}:")
    print(f"    scores = {sorted(scores)}")
    print(f"    sum_f:  actual={sum_f}, theo={theo_sum_f}, match={sum_f==theo_sum_f}")
    print(f"    sum_f2: actual={sum_f2}, theo={theo_sum_f2}, match={sum_f2==theo_sum_f2}")
    print(f"    sig_consec = {sig_consec}")
    print(f"    sig_overlap: direct={sig_overlap_direct}, formula={sig_overlap_formula}, match={sig_overlap_direct==sig_overlap_formula}")
    print(f"    sig_disjoint = {sig_disjoint}")
    print(f"    sig_non_consec = {sig_non_consec}")
    print(f"    sum_fb: actual={sum_fb}, theo={(n-1)*theo_sum_f - theo_sum_f2}")

    # Now compute w_{n-3} = -(1/2)*sum_fb + C(n-1,2)/4*n!
    w_nm3 = -0.5 * sum_fb + comb(n-1,2)/4 * factorial(n)
    t3 = count_3cycles(A)
    expected = 2*factorial(n-2)*t3 - factorial(n-2)*comb(n,3)/2

    print(f"    w_{{n-3}}: computed={w_nm3:.1f}, expected=2*{factorial(n-2)}*{t3}-{factorial(n-2)*comb(n,3)//2}={expected:.1f}")
    print(f"    Match: {abs(w_nm3 - expected) < 0.01}")

# =====================================================================
# STEP 4: Derive closed form
# =====================================================================
print("\n" + "=" * 70)
print("DERIVATION OF CLOSED FORM")
print("=" * 70)

print("""
  Given:
    w_{n-3} = -(1/2)*sum_P f_P*b_P + C(n-1,2)/4 * n!

    sum_P f_P*b_P = (n-1)*sum_P f_P - sum_P f_P^2

    sum_P f_P = C(n,2) * (n-1)!

    sum_P f_P^2 = sum_P f_P + 2*(Sigma_consec + Sigma_non_consec)

  Now:
    Sigma_consec = (n-2)! * [(n-1)*C(n,2) - sum s_b^2]

    And sum s_b^2 = C(n,2) + 2*C(n,3) - 2*t_3
                  (from t_3 = C(n,3) - sum C(s_b,2) = C(n,3) - (ssq - C(n,2))/2)

  So Sigma_consec = (n-2)! * [(n-1)*C(n,2) - C(n,2) - 2*C(n,3) + 2*t_3]
                  = (n-2)! * [(n-2)*C(n,2) - 2*C(n,3) + 2*t_3]

  Note: (n-2)*C(n,2) = (n-2)*n*(n-1)/2
        2*C(n,3) = 2*n*(n-1)*(n-2)/6 = n*(n-1)*(n-2)/3
        (n-2)*C(n,2) - 2*C(n,3) = n*(n-1)*(n-2)/2 - n*(n-1)*(n-2)/3
                                 = n*(n-1)*(n-2)*[1/2 - 1/3]
                                 = n*(n-1)*(n-2)/6
                                 = C(n,3)

  So: Sigma_consec = (n-2)! * [C(n,3) + 2*t_3]

  The key t_3-dependent part is 2*(n-2)!*t_3.

  Now I need Sigma_non_consec. Let me verify numerically...
""")

# Verify Sigma_consec formula
for n in [4, 5, 6]:
    random.seed(n * 88)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    # Direct Sigma_consec
    sig_c_direct = 0
    for p in permutations(range(n)):
        for i in range(n-2):
            sig_c_direct += A[p[i]][p[i+1]] * A[p[i+1]][p[i+2]]

    t3 = count_3cycles(A)
    scores = [sum(A[i]) for i in range(n)]
    ssq = sum(s*s for s in scores)

    # Formula 1: (n-2)! * [(n-1)*C(n,2) - ssq]
    formula1 = factorial(n-2) * ((n-1)*comb(n,2) - ssq)

    # Formula 2: (n-2)! * [C(n,3) + 2*t3]
    formula2 = factorial(n-2) * (comb(n,3) + 2*t3)

    print(f"  n={n}: Sigma_consec: direct={sig_c_direct}, formula1={formula1}, formula2={formula2}")
    print(f"         ssq={ssq}, C(n,2)+2*C(n,3)-2*t3={comb(n,2)+2*comb(n,3)-2*t3}")
    print(f"         Match: formula1={'✓' if sig_c_direct==formula1 else '✗'}, formula2={'✓' if sig_c_direct==formula2 else '✗'}")

# =====================================================================
# Now compute Sigma_non_consec
# =====================================================================
print("\n--- Sigma_non_consec ---")
for n in [4, 5]:
    random.seed(n * 77)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    # Direct computation
    sig_nc = 0
    for p in permutations(range(n)):
        for i in range(n-1):
            for j in range(i+2, n-1):
                sig_nc += A[p[i]][p[i+1]] * A[p[j]][p[j+1]]

    t3 = count_3cycles(A)
    scores = [sum(A[i]) for i in range(n)]
    ssq = sum(s*s for s in scores)

    # Number of non-consecutive position pairs
    num_nc = comb(n-1,2) - (n-2)

    # All 4-tuples (a,b,c,d) of distinct vertices with T(a,b)*T(c,d)
    sig_disjoint_edges = 0
    for a in range(n):
        for b in range(n):
            if b == a: continue
            for c in range(n):
                if c == a or c == b: continue
                for d in range(n):
                    if d == a or d == b or d == c: continue
                    sig_disjoint_edges += A[a][b] * A[c][d]

    # Expected: num_nc * (n-4)! * sig_disjoint_edges
    expected_nc = num_nc * factorial(max(n-4,0)) * sig_disjoint_edges if n >= 4 else 0

    print(f"\n  n={n}:")
    print(f"    Sigma_non_consec: direct={sig_nc}")
    print(f"    num_nc={num_nc}, (n-4)!={factorial(max(n-4,0))}, sig_disjoint={sig_disjoint_edges}")
    print(f"    Expected: {expected_nc}")
    print(f"    Match: {'✓' if sig_nc == expected_nc else '✗'}")

    if sig_nc != expected_nc and n >= 4:
        # Check: maybe the count per position pair varies?
        # For each position pair (i,j) with j > i+1, compute the sum
        for i in range(n-1):
            for j in range(i+2, n-1):
                s = 0
                for p in permutations(range(n)):
                    s += A[p[i]][p[i+1]] * A[p[j]][p[j+1]]
                print(f"      (i,j)=({i},{j}): sum={s}, (n-4)!*disjoint={factorial(max(n-4,0))*sig_disjoint_edges}")

# =====================================================================
# FULL FORMULA: putting it all together
# =====================================================================
print("\n" + "=" * 70)
print("FULL FORMULA DERIVATION")
print("=" * 70)

for n in [5, 7]:
    print(f"\n  n={n}:")
    random.seed(n * 55)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3cycles(A)
    scores = [sum(A[i]) for i in range(n)]
    ssq = sum(s*s for s in scores)

    # Direct w_{n-3}
    w_direct = 0
    for p in permutations(range(n)):
        s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]
        e2 = sum(s[i]*s[j] for i in range(n-1) for j in range(i+1, n-1))
        w_direct += e2

    # Formula: 2*(n-2)!*t_3 - (n-2)!*C(n,3)/2
    w_formula = 2*factorial(n-2)*t3 - factorial(n-2)*comb(n,3)//2

    print(f"    t3={t3}, ssq={ssq}")
    print(f"    w_{{n-3}} direct = {w_direct:.1f}")
    print(f"    w_{{n-3}} formula = {w_formula}")
    print(f"    Match: {'✓' if abs(w_direct - w_formula) < 0.01 else '✗'}")

    # Now derive: w_{n-3} = -(1/2)*sum_fb + C(n-1,2)/4*n!
    sum_fb_direct = 0
    for p in permutations(range(n)):
        f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        sum_fb_direct += f*(n-1-f)

    w_from_fb = -0.5*sum_fb_direct + comb(n-1,2)/4*factorial(n)
    print(f"    sum_fb = {sum_fb_direct}")
    print(f"    w from fb = {w_from_fb:.1f}")

    # Theoretical: sum_fb = (n-1)*C(n,2)*(n-1)! - sum_f2
    # sum_f2 = C(n,2)*(n-1)! + 2*(Sigma_c + Sigma_nc)
    # Sigma_c = (n-2)! * [C(n,3) + 2*t3]
    sig_c = factorial(n-2) * (comb(n,3) + 2*t3)
    print(f"    Sigma_consec = {sig_c}")

    # For the full formula to work, we need:
    # w_{n-3} = -(1/2)*[(n-1)*sum_f - sum_f2] + C(n-1,2)/4*n!
    #         = -(1/2)*(n-1)*sum_f + (1/2)*sum_f2 + C(n-1,2)/4*n!
    #         = -(1/2)*(n-1)*C(n,2)*(n-1)! + (1/2)*[sum_f + 2*sig_c + 2*sig_nc] + C(n-1,2)/4*n!
    #
    # The t3-dependent part is: (1/2)*2*sig_c = sig_c = (n-2)!*[C(n,3)+2*t3]
    # The t3-coefficient is: (n-2)!*2 = 2*(n-2)!
    # The constant parts... let me just verify the final answer.

    print(f"    t3 coeff: {2*factorial(n-2)} = 2*{factorial(n-2)}")
    print(f"    constant: -{factorial(n-2)*comb(n,3)//2}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
