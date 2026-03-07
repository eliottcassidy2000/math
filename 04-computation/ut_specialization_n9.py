#!/usr/bin/env python3
"""
What does the OCF specialization of U_T give at n=9?

U_T = sum over sigma in S(T, T^op) of (-1)^phi * p_type(sigma)

OCF spec: p_1->1, p_{2k+1}->2, p_even->0

This should give ps_spec(U_T). Does it equal H(T)?

From gs_specialization_check.py, we CONFIRMED this at n=5 and n=7.
But those had NO mixed-direction permutations.

At n=9 with mixed perms, what happens?

Key point: sigma in S(T, T^op) means each nontrivial cycle of sigma
is a directed cycle of T OR T^op. The (-1)^phi factor uses phi =
sum over T-CYCLE components of (len-1).

For all-odd cycles: phi = sum (even) = even => (-1)^phi = +1.
For all-even-containing: p_type has an even-indexed part => specializes to 0.

So under OCF spec: ps_spec(U_T) = sum over sigma (all nontrivial odd,
each T or T^op) of 2^{nontrivial}.

From our computation: this gives 6425, not 2101 = H(T).

But at n=5,7 it DID give H(T)! Because there were no mixed perms.

So either:
(a) The GS Corollary 20 says something slightly different, or
(b) My computation has a bug.

Let me re-examine. Corollary 20 of Irving-Omar says:
  ham(D_bar) = sum_{sigma in S(D), all cycles of sigma have odd length} 2^{psi(sigma)}

where S(D) = {permutations whose nontrivial cycles are directed cycles of D}.
NOT "of D or D_bar"!

So S(D) uses ONLY D-cycles, not D_bar-cycles.

That's the key. Let me verify: the set S(D) in the theorem is
permutations whose nontrivial cycles are directed cycles of D only.
The D_bar cycles come from the omega involution on U_D, not from S(D).

So the correct formula is:
  ham(T) = ham(T^op) = sum over sigma in S_T(all nontrivial cycles are odd T-cycles) 2^{psi}
         = I(Omega(T), 2)

And this is EXACTLY the T-only sum, which equals H(T) = 2101. ✓

The confusion was that U_T involves S(T, T^op) = permutations with
cycles in T OR T^op. But the corollary specifically restricts to S(T).

Wait, but U_T as defined by Grinberg-Stanley uses BOTH T and T^op:
U_T = sum over sigma in S_V(T, T_bar) of (-1)^phi p_type

And S_V(T, T_bar) = permutations whose nontrivial cycles are
directed cycles of T OR T_bar (= T^op for tournaments).

The corollary then restricts to S(D) which is sigma with all nontrivial
cycles being D-cycles only.

So the derivation goes:
1. U_T = sum over S_V(T, T^op) of (-1)^phi p_type
2. omega(U_T) = U_{T^op} (involution)
3. ps_1(omega(f)) maps f to a specific value
4. Corollary 20 uses ham(D_bar) = ps_1(U_D)(1)... but with a TWIST.

Actually, I think the issue is that ps_1(U_T)(1) does NOT directly
use the OCF specialization. The OCF specialization comes from
a DIFFERENT evaluation.

Let me just verify: what is ps_1(U_T)(1) and ps_spec(U_T)?

opus-2026-03-07-S37
"""
from itertools import permutations
from collections import defaultdict
import random

n = 7
QR = {1, 2, 4}
paley_edges = [(i, j) for i in range(n) for j in range(n)
               if i != j and (j - i) % n in QR]
edge_set = set(paley_edges)
opp_set = set((j, i) for (i, j) in paley_edges)

# Compute U_T's p_lambda coefficients
coeffs = defaultdict(int)  # lambda -> (-1)^phi coefficient

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

print(f"Paley T_7: U_T p-coefficients:")
for lam in sorted(coeffs.keys()):
    if coeffs[lam] != 0:
        print(f"  p_{lam}: {coeffs[lam]}")

# Specialization 1: p_k -> 1 for all k (this is ps_1 at m=1)
ps1 = sum(coeffs[lam] for lam in coeffs)
print(f"\nps_1(U_T)(1) = {ps1}")

# Specialization 2: OCF spec: p_1->1, p_{odd>=3}->2, p_even->0
ocf = 0
for lam, c in coeffs.items():
    val = c
    for part in lam:
        if part == 1:
            val *= 1
        elif part % 2 == 1:
            val *= 2
        else:
            val *= 0
    ocf += val
print(f"OCF spec: {ocf}")
print(f"H(T_7) = 189")
print(f"OCF spec = H? {ocf == 189}")

# Now verify Corollary 20 directly:
# ham(D_bar) where D = T, D_bar = T^op.
# = sum over sigma in S(T) with all nontrivial cycles odd of 2^{psi(sigma)}
# where S(T) = permutations whose nontrivial cycles are directed T-cycles ONLY.

S_T_sum = 0
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
    nontrivial = 0
    for cyc in cycles:
        if len(cyc) == 1:
            continue
        if len(cyc) % 2 == 0:
            ok = False
            break
        is_T = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set for i in range(len(cyc)))
        if not is_T:
            ok = False
            break
        nontrivial += 1

    if ok:
        S_T_sum += 2**nontrivial

print(f"\nCorollary 20 S(T) sum: {S_T_sum}")
print(f"ham(T^op) = ham(T) = 189")
print(f"Match: {S_T_sum == 189}")

# The KEY question: does ps_spec(U_T) = S_T_sum?
# U_T includes both T and T^op cycles, but under OCF spec, does it reduce to S_T_sum?

# Let me decompose: which terms of U_T survive the OCF spec?
print(f"\n=== Decomposing U_T under OCF spec ===")
# Terms with all-odd type: the p_type specializes to 2^k.
# Terms with any even part: killed by p_even -> 0.

# Among all-odd-type terms: these are perms with all odd cycles.
# Some have only T-cycles, some only T^op, some mixed.
# At n=7: no mixed, so OCF spec of U_T = (T-only sum) + (T^op-only sum)
# = 2 * H(T)? No, that doesn't match.

# Let me compute more carefully.
# For all-odd sigma in S(T, T^op):
# (-1)^phi * prod p_{L_i} -> (-1)^phi * prod 2^{[L_i>=3]}
# = (+1) * 2^{nontrivial} (since phi always even for all-odd)

# So OCF spec of U_T = sum over ALL valid all-odd sigma of 2^{nontrivial}
# Including T-only, T^op-only, and mixed.

# At n=7: no mixed, and we showed T-only = T^op-only = 189.
# So OCF spec of U_T = 189 + 189 - 1 = 377? (minus 1 for identity counted in both)
# But we computed OCF spec = 189. Contradiction!

# Wait, identity is counted in T-only and Top-only. Let me be careful.
# Identity permutation: all fixed points. No nontrivial cycles.
# It IS in S(T, T^op) (trivially). Weight = 2^0 = 1.
# T-only count includes identity: yes (vacuously, all 0 nontrivial cycles are T-cycles)
# T^op-only count includes identity: yes.
# So T-only + T^op-only double-counts the identity.

# Actually, for T-only: sum over sigma with all nontrivial cycles being T-cycles.
# This includes sigma with ZERO nontrivial cycles (identity).
# Similarly for T^op-only.

# So: (all valid) = (T-only) + (Top-only) + (mixed) - (counted in both T-only and Top-only)
# The identity is the only sigma in both T-only AND Top-only (vacuously).
# So: all valid = T-only + Top-only + mixed - 1

# At n=7: mixed = 0, so all valid = 189 + 189 - 1 = 377.
# But OCF spec = 189. HOW?

# The answer must be that my understanding of S(T, T_bar) is wrong.
# Maybe S(T, T_bar) does NOT include permutations where ALL cycles go one way.

# Let me re-read: S_V(T, T_bar) from Grinberg-Stanley:
# sigma such that each nontrivial cycle is a directed cycle of T or of T_bar.
# This IS the set I've been using.

# But maybe the SIGN cancels things?
# For a T-only perm: phi = sum (L-1) for T-cycles. All even. Sign = +1.
# For a T^op-only perm: phi = 0 (no T-cycles). Sign = +1.
# For mixed: phi = sum (L-1) for T-cycle parts. All even. Sign = +1.

# So the OCF spec of U_T = sum of all valid * 2^k = T-only + Top-only + mixed - 1
# = 189 + 189 - 1 + 0 = 377 at n=7.
# But the COMPUTED value was 189!

# There's a fundamental error. Let me recompute the OCF spec at n=7.
# From the coefficients dict:
# The OCF spec maps p_lambda -> product of p_{L_i} values.
# For U_T, each p_lambda coefficient already includes the (-1)^phi sign.

print(f"\n=== Re-examining at n=7 ===")
# Let me separate the p_lambda coefficients by direction type
T_only_coeffs = defaultdict(int)
Top_only_coeffs = defaultdict(int)
mixed_coeffs = defaultdict(int)

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
    has_T = False
    has_Top = False

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
            has_T = True
        if is_Top:
            has_Top = True

    if not ok:
        continue

    cycle_type = tuple(sorted([len(c) for c in cycles], reverse=True))
    sign = (-1)**phi

    if has_T and not has_Top:
        T_only_coeffs[cycle_type] += sign
    elif has_Top and not has_T:
        Top_only_coeffs[cycle_type] += sign
    elif has_T and has_Top:
        mixed_coeffs[cycle_type] += sign
    else:
        # Identity: put in T_only
        T_only_coeffs[cycle_type] += sign

print(f"T-only p-coefficients:")
for lam in sorted(T_only_coeffs.keys()):
    if T_only_coeffs[lam] != 0:
        print(f"  p_{lam}: {T_only_coeffs[lam]}")

ocf_T = 0
for lam, c in T_only_coeffs.items():
    val = c
    for part in lam:
        if part == 1:
            val *= 1
        elif part % 2 == 1:
            val *= 2
        else:
            val *= 0
    ocf_T += val
print(f"OCF spec of T-only part: {ocf_T}")

print(f"\nTop-only p-coefficients:")
for lam in sorted(Top_only_coeffs.keys()):
    if Top_only_coeffs[lam] != 0:
        print(f"  p_{lam}: {Top_only_coeffs[lam]}")

ocf_Top = 0
for lam, c in Top_only_coeffs.items():
    val = c
    for part in lam:
        if part == 1:
            val *= 1
        elif part % 2 == 1:
            val *= 2
        else:
            val *= 0
    ocf_Top += val
print(f"OCF spec of Top-only part: {ocf_Top}")
print(f"\nOCF spec total: {ocf_T + ocf_Top}")
print(f"H(T) = 189")
