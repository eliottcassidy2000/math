#!/usr/bin/env python3
"""
Omega involution constraint on U_T coefficients.

For tournaments: omega(U_T) = U_{T^op}.
omega acts on power sums: omega(p_lambda) = epsilon_lambda * p_lambda
where epsilon_lambda = (-1)^{|lambda| - l(lambda)} = (-1)^{sum(lambda_i - 1)}.

Since ham(T) = ham(T^op), we get ps_1(U_T)(1) = ps_1(U_{T^op})(1).
But U_{T^op} = omega(U_T), so ps_1(omega(U_T))(1) = ps_1(U_T)(1).

This gives: sum_lambda c_lambda * epsilon_lambda = sum_lambda c_lambda

i.e., sum_lambda c_lambda * (epsilon_lambda - 1) = 0

Equivalently: the sum of c_lambda over lambda with epsilon_lambda = -1 is zero.

epsilon_lambda = -1 iff |lambda| - l(lambda) is odd
iff sum(lambda_i - 1) is odd
iff (number of even parts) + (number of odd parts >= 3) is odd
Actually: sum(lambda_i - 1) = |lambda| - l(lambda) = n - l(lambda).

For n=7: epsilon_lambda = (-1)^{7 - l(lambda)}.
epsilon = -1 when 7 - l is odd, i.e., l is even.

So the constraint is: sum of c_lambda over even-length partitions = 0.

Let me verify this and explore what other constraints omega gives.

opus-2026-03-07-S37
"""
from itertools import permutations
from collections import defaultdict

def compute_UT_coeffs(n, edges):
    """Compute c_lambda for U_T = sum c_lambda p_lambda."""
    edge_set = set(edges)
    opp_set = set((j, i) for (i, j) in edges)

    coeffs = defaultdict(int)
    for sigma in permutations(range(n)):
        visited = [False] * n
        cycles = []
        for start in range(n):
            if visited[start]:
                continue
            cycle = []
            curr = start
            while not visited[curr]:
                visited[curr] = True
                cycle.append(curr)
                curr = sigma[curr]
            cycles.append(tuple(cycle))

        ok = True
        phi = 0
        for cyc in cycles:
            if len(cyc) == 1:
                continue
            is_T = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set for i in range(len(cyc)))
            is_Top = all((cyc[i], cyc[(i+1) % len(cyc)]) in opp_set for i in range(len(cyc)))
            if not is_T and not is_Top:
                ok = False
                break
            if is_T:
                phi += len(cyc) - 1

        if ok:
            cycle_type = tuple(sorted([len(c) for c in cycles], reverse=True))
            coeffs[cycle_type] += (-1)**phi

    return dict(coeffs)

# Paley T_7
n = 7
QR = {1, 2, 4}
edges7 = [(i, j) for i in range(n) for j in range(n)
          if i != j and (j - i) % n in QR]

coeffs = compute_UT_coeffs(n, edges7)

print("=== Paley T_7: U_T p-coefficients ===")
print(f"{'lambda':<25} {'c_lambda':>10} {'l(lambda)':>10} {'epsilon':>10} {'c*epsilon':>10}")
print("-" * 70)

sum_all = 0
sum_even_l = 0
sum_odd_l = 0

for lam in sorted(coeffs.keys()):
    c = coeffs[lam]
    if c == 0:
        continue
    l = len(lam)
    epsilon = (-1)**(n - l)
    print(f"p_{str(lam):<22} {c:>10} {l:>10} {epsilon:>10} {c*epsilon:>10}")
    sum_all += c
    if l % 2 == 0:
        sum_even_l += c
    else:
        sum_odd_l += c

print(f"\nsum of all c_lambda = {sum_all}")
print(f"sum of c_lambda with even l = {sum_even_l}")
print(f"sum of c_lambda with odd l = {sum_odd_l}")
print(f"sum * epsilon = {sum_odd_l - sum_even_l} (for n=7 odd: epsilon_even_l = -1)")
print(f"Constraint: sum_even_l should be 0 (so that sum = sum*epsilon)")
print(f"Verified: {sum_even_l == 0}")

print(f"\n=== Now check OCF specialization vs ps_1 ===")
# ps_1(m=1): sum of c_lambda (all terms contribute 1)
ps1 = sum_all
# OCF: p_1->1, p_odd>=3->2, p_even->0
ocf = 0
for lam, c in coeffs.items():
    if c == 0:
        continue
    val = c
    for part in lam:
        if part == 1:
            val *= 1
        elif part % 2 == 1:
            val *= 2
        else:
            val *= 0
    ocf += val

print(f"ps_1(U_T)(1) = {ps1}")
print(f"OCF spec = {ocf}")
print(f"Both should be H(T) = 189: ps1={ps1==189}, ocf={ocf==189}")

# Now check: which partitions have ALL odd parts?
all_odd_sum = 0
has_even_sum = 0
for lam, c in coeffs.items():
    if c == 0:
        continue
    if all(p % 2 == 1 for p in lam):
        all_odd_sum += c
    else:
        has_even_sum += c

print(f"\nsum over all-odd-part lambdas = {all_odd_sum}")
print(f"sum over has-even-part lambdas = {has_even_sum}")
print(f"Check: {all_odd_sum} + {has_even_sum} = {all_odd_sum + has_even_sum} = {sum_all}")

# OCF kills even-part lambdas. For all-odd lambdas, each part >= 3 gets factor 2.
# So OCF = sum over all-odd lam of c_lam * 2^{#{parts >= 3}}
ocf_decomp = 0
for lam, c in coeffs.items():
    if c == 0:
        continue
    if not all(p % 2 == 1 for p in lam):
        continue
    nontrivial = sum(1 for p in lam if p >= 3)
    ocf_decomp += c * 2**nontrivial

print(f"\nOCF = sum over all-odd lam of c * 2^nontrivial = {ocf_decomp}")
print(f"This should be {ps1} if OCF = ps1")
print(f"Match: {ocf_decomp == ps1}")

# INSIGHT: has_even_sum = 0 means no contribution from even-part lambdas!
# That's WHY ps_1(1) = OCF for the all-odd part only... wait, OCF also
# multiplies by 2^nontrivial, not just 1. So they CAN'T be equal unless
# the weights happen to match.

# Let me check: at n=7, is ocf = ps1 = 189?
print(f"\nocf={ocf}, ps1={ps1}. Equal: {ocf == ps1}")

# They ARE both 189! But OCF applies different weights (2^k) while ps1
# applies weight 1 to each term. The only way they can be equal is if
# the "extra" weight from 2^k on all-odd terms exactly compensates for
# the "lost" contribution from even-part terms.

# Check: all-odd sum = {all_odd_sum}, has-even sum = {has_even_sum}
# OCF = all-odd weighted = {ocf_decomp}
# ps1 = all-odd + has-even = {sum_all}
# So: all-odd weighted - all-odd unweighted = has-even unweighted
# i.e., sum_odd c*(2^k - 1) = -sum_even c

print(f"\n=== Balance equation ===")
extra_from_weighting = ocf_decomp - all_odd_sum
print(f"Extra from 2^k weighting of all-odd terms: {extra_from_weighting}")
print(f"Contribution from even-part terms: {has_even_sum}")
print(f"Balance: {extra_from_weighting} + {has_even_sum} = {extra_from_weighting + has_even_sum}")
print(f"(Should be 0 if OCF = ps1)")
