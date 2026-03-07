#!/usr/bin/env python3
"""
Check the SIGNED mixed-direction sum at n=9.

The Rédei-Berge sum: sum over sigma with cycles that are directed T or T^op
cycles, weighted by (-1)^phi * p_{type(sigma)}.

phi = sum over T-cycles of (len - 1).

Under OCF spec: p_1->1, p_{odd>=3}->2, p_even->0.

For all-odd cycles: p_type -> 2^{nontrivial} (if all odd) or 0 (if any even).

The SIGNED weight is (-1)^phi * 2^{nontrivial}.

For pure T-cycle perms: phi = sum (L_i - 1) = even (each L_i odd => L_i-1 even).
So sign = +1. Contribution = +2^k.

For pure T^op-cycle perms: phi = 0 (no T-cycles). Sign = +1. Contribution = +2^k.

For mixed (some T, some T^op): phi = sum over T-cycle components of (L_i-1).
Each L_i is odd, so L_i-1 is even. phi is sum of even numbers = even.
So sign = +1 for mixed too!

This means EVERY valid permutation contributes POSITIVELY.

But we showed T-only = H(T) and total > H(T).
If total with correct signs still > H(T), then the GS formula must
be wrong about counting both directions.

Actually wait — I need to re-examine what the GS formula says.
GS says: ps_1(U_T)(1) = H(T). And U_T involves cycles of T and T_bar = T^op.
But ps_1 at m=1 gives sum of coefficients (since ps_1(p_lambda)(1) = 1).
And we showed sum of coefficients = 3 at C_5 (from gs_specialization_check.py).
H(C_5) = 15.

So ps_1(U_T)(1) ≠ H(T) in general!
Instead, the OCF comes from the SPECIFIC specialization p_1->1, p_{odd}->2, p_even->0.

And this specialization applied to U_T = sum (-1)^phi * p_type gives...
what exactly?

Let me compute it directly.

opus-2026-03-07-S37
"""
from itertools import permutations
from collections import defaultdict
import random

n = 9

def compute_all_sums(n, A):
    edge_set = set()
    opp_set = set()
    for i in range(n):
        for j in range(n):
            if i != j:
                if A[i][j]:
                    edge_set.add((i, j))
                else:
                    opp_set.add((i, j))

    # Held-Karp for H(T)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    H = sum(dp[full][v] for v in range(n))

    # Compute the SIGNED Rédei-Berge sum under OCF specialization
    # U_T = sum over sigma in S(T, T^op) of (-1)^phi * p_type
    # OCF spec: p_1->1, p_{2k+1}->2, p_{2k}->0

    signed_total = 0
    signed_T_only = 0
    signed_Top_only = 0
    signed_mixed = 0

    for sigma_tuple in permutations(range(n)):
        sigma = list(sigma_tuple)
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
        nontrivial = 0
        all_odd = True
        has_T = False
        has_Top = False

        for cyc in cycles:
            if len(cyc) == 1:
                continue
            if len(cyc) % 2 == 0:
                all_odd = False
                ok = False
                break
            nontrivial += 1

            is_T = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set
                       for i in range(len(cyc)))
            is_Top = all((cyc[i], cyc[(i+1) % len(cyc)]) in opp_set
                        for i in range(len(cyc)))

            if not is_T and not is_Top:
                ok = False
                break

            if is_T:
                phi += len(cyc) - 1
                has_T = True
            if is_Top:
                has_Top = True

        if not ok or not all_odd:
            continue

        sign = (-1)**phi
        weight = sign * (2**nontrivial)

        signed_total += weight

        if has_T and not has_Top:
            signed_T_only += weight
        elif has_Top and not has_T:
            signed_Top_only += weight
        elif has_T and has_Top:
            signed_mixed += weight
        else:
            # Identity
            signed_T_only += weight
            signed_Top_only += weight

    return H, signed_total, signed_T_only, signed_Top_only, signed_mixed

# Test
rng = random.Random(42)
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if rng.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

print("Random n=9 tournament (seed=42):")
H, total, T_only, Top_only, mixed = compute_all_sums(n, A)
print(f"  H(T) = {H}")
print(f"  Signed total:    {total}")
print(f"  Signed T-only:   {T_only}")
print(f"  Signed T^op-only: {Top_only}")
print(f"  Signed mixed:    {mixed}")
print(f"  total = H? {total == H}")
print(f"  T-only = H? {T_only == H}")
print(f"  T_only + mixed = {T_only + mixed}")

# Wait: phi = sum over T-cycles of (L-1). For odd L, L-1 is even.
# So (-1)^phi = (-1)^(sum of even numbers) = +1 always.
# So the sign is ALWAYS +1 for all-odd cycles. The signed sum = unsigned sum.
print(f"\nNote: for all-odd cycles, phi is always even, so sign is always +1.")
print(f"This means signed sum = unsigned sum.")
print(f"Therefore: the full Rédei-Berge sum with OCF spec does NOT equal H(T)")
print(f"when mixed-direction cycles are included.")
print(f"")
print(f"The OCF formula H(T) = I(Omega(T), 2) counts ONLY T-direction cycles.")
print(f"The factor of 2 per cycle in the specialization is NOT from T/T^op duality.")
