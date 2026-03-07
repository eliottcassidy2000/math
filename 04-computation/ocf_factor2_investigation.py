#!/usr/bin/env python3
"""
Investigate WHY the OCF specialization gives factor 2 per nontrivial odd cycle.

In the OCF: H(T) = sum over cycle-disjoint collections of odd T-cycles,
weighted by 2^{|collection|}. The factor 2^k comes from:

Permutation perspective:
- A permutation sigma with k nontrivial cycles (all odd, all T-cycles)
  contributes (-1)^{phi(sigma)} to U_T's p_{type(sigma)} coefficient.
- Under the specialization p_1->1, p_{2k+1}->2, this becomes
  (-1)^{phi(sigma)} * 2^k.
- phi(sigma) = sum over T-cycles of (len - 1).
- For an odd cycle of length L, len-1 is even, so (-1)^{L-1} = +1.
- Therefore every contributing permutation has (-1)^phi = +1.
- The factor 2 per cycle comes from p_{2k+1} -> 2.

But why 2 specifically? Is there a combinatorial explanation?

For a directed cycle C of odd length L in T, the "2" might count:
(a) Two orientations of the cycle? But C has a fixed direction.
(b) Two colorings? In the hard-core model, fugacity lambda = 2.
(c) Something else?

Let me investigate: in the Rédei-Berge framework, p_k specializes to
"number of colorings with k consecutive same-colors" or something.

Actually, the principal specialization ps_1(p_k)(m) = m (since p_k = x1^k + ... + xm^k,
and setting x1=...=xm=1 gives m).

So ps_1(U_T)(m) = sum of all coefficients * m^{l(lambda)}.

At m=1: every term gives 1, so sum of coefficients = ps_1(U_T)(1).
This should be the number of acyclic orderings of T, not H(T).

The OCF specialization is DIFFERENT from ps_1. It's:
p_1 -> 1, p_{2k+1} -> 2, p_{2k} -> 0.

This is NOT a principal specialization. It's a "partial" specialization
that zeroes out even power sums.

Can we understand this as evaluating at specific x-values?

p_k = x1^k + x2^k + ...
Setting x1 = 1, x2 = 1 (two variables, both = 1):
p_k = 1 + 1 = 2 for all k.
But we want p_1 = 1, not 2.

Setting x1 = 1, x2 = omega (primitive root of unity?):
p_k = 1 + omega^k. For omega = -1: p_k = 1 + (-1)^k.
So p_odd = 0, p_even = 2. That's the OPPOSITE of what we want.

Setting x1 = 1, x2 = i (fourth root of unity):
p_k = 1 + i^k. p_1 = 1+i, p_2 = 1-1=0, p_3 = 1-i, p_4 = 1+1=2.
Not quite right either.

Setting x1 = 1, x2 = zeta (primitive 2m-th root for any m):
p_k = 1 + zeta^k.

For the OCF spec: p_1 = 1, p_3 = 2, p_5 = 2, p_7 = 2, p_2 = 0, p_4 = 0.

From 1 + zeta^k = ? We need:
k=1: 1 + zeta = 1 => zeta = 0. But then all p_k = 1. No good.

Two variables x1, x2 with x1 + x2 = p_1 = 1, x1^3 + x2^3 = 2:
(x1+x2)^3 = x1^3 + 3x1^2 x2 + 3x1 x2^2 + x2^3
1 = 2 + 3 x1 x2 (x1 + x2) = 2 + 3 x1 x2
So x1 x2 = -1/3.
And x1 + x2 = 1 => x1, x2 = (1 ± sqrt(1+4/3))/2 = (1 ± sqrt(7/3))/2.
Then p_5 = x1^5 + x2^5 = (x1+x2)(x1^4 - x1^3 x2 + x1^2 x2^2 - x1 x2^3 + x2^4)
Hmm, this gets complicated. Let me just check if any finite number of
variables gives the OCF spec.

Actually, the proper framework:
For ANY symmetric function f, we can define a "ring homomorphism"
Lambda -> R by p_k -> a_k for any choice of values a_k.
The OCF spec is just: p_1 -> 1, p_{2k+1} -> 2, p_{2k} -> 0.

This doesn't need to come from substituting specific x-values.
It's a general algebra homomorphism on the ring of symmetric functions.

But the question remains: what makes this particular specialization
"natural" for tournaments?

Hypothesis: The factor 2 counts the two endpoints of each cycle.
In a directed cycle of length L, there are L vertices, each of which
could serve as the "starting point" for some construction. But that
gives L, not 2.

Alternative: The 2 comes from the BINARY nature of tournaments.
Each edge has two orientations. In a cycle of length L, reversing
ALL edges gives another directed cycle (in T^op). So cycles come
in pairs: C and C^op. But C^op is a cycle of T^op, not T.

Wait — for a tournament T, the opposite T^op has the SAME number
of Hamiltonian paths: ham(T) = ham(T^op). This is the path-reversal
bijection.

And the OCF uses cycles of T only (not T^op). The factor 2 in the
specialization applies to cycles of length >= 3.

Let me think about it from the Rédei-Berge function directly.
U_T involves permutations whose cycles are directed cycles of T OR T^op.
The sign is (-1)^{sum of (len-1) for T-cycles}.

For permutations with all odd cycles that are T-cycles:
sign = (-1)^{sum (even)} = +1.
Specialization gives 2^k (for k nontrivial cycles).

For permutations with all odd cycles that are T^op-cycles:
sign = (-1)^0 = +1.
Specialization also gives 2^k.

For MIXED permutations (some T-cycles, some T^op-cycles, all odd):
sign = (-1)^{sum (L_i - 1) for T-cycle components} = +1 (all even).
Specialization gives 2^k.

So EVERY permutation with all odd nontrivial cycles contributes +2^k,
regardless of whether cycles are T-cycles or T^op-cycles.

But in our OCF, we only count T-cycles. So what happens?

The RHS sum splits:
sum_sigma 2^k = sum over collections of disjoint odd cycles
                (each either a T-cycle or T^op-cycle) of 2^k

This counts MORE than just T-cycle collections.

Let's verify: are T-cycles and T^op-cycles in disjoint orbits?
A directed cycle C in T: v1->v2->...->vL->v1.
The reverse C^rev: v1->vL->...->v2->v1. This is a cycle of T^op.
As a permutation, C and C^rev are DIFFERENT.

So the full sum includes:
- Collections of T-cycles only: these are independent sets in Omega(T)
- Collections of T^op-cycles only: these are independent sets in Omega(T^op)
- Mixed collections

Hmm, but ham(T^op) = ham(T), and Omega(T^op) is isomorphic to Omega(T)
via the reversal map. So I(Omega(T^op), x) = I(Omega(T), x).

Actually wait. Let me re-examine. From the GS specialization check script,
we confirmed that the sum over permutations with all odd T-cycles (only T-direction),
weighted by 2^{nontrivial}, gives H(T).

Let me re-check: are T^op-cycles also counted?

opus-2026-03-07-S37
"""
from itertools import permutations
from collections import defaultdict

def verify_cycle_direction(n, edges):
    """Check: in the Rédei-Berge sum, do T^op cycles contribute separately?"""
    edge_set = set(edges)
    opp_edges = set((j, i) for i in range(n) for j in range(n)
                    if i != j and (i, j) not in edge_set)

    # Method 1: Sum over all valid permutations (T and T^op cycles mixed)
    total_mixed = 0
    # Method 2: Sum over T-only cycle permutations
    total_T_only = 0
    # Method 3: Sum over T^op-only
    total_Top_only = 0

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

        # Check each nontrivial cycle
        ok_mixed = True
        ok_T = True
        ok_Top = True
        nontrivial = 0
        all_odd = True

        for cyc in cycles:
            if len(cyc) == 1:
                continue
            if len(cyc) % 2 == 0:
                all_odd = False
                ok_mixed = ok_T = ok_Top = False
                break
            nontrivial += 1

            is_T = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set
                       for i in range(len(cyc)))
            is_Top = all((cyc[i], cyc[(i+1) % len(cyc)]) in opp_edges
                        for i in range(len(cyc)))

            if not is_T and not is_Top:
                ok_mixed = False
                ok_T = False
                ok_Top = False
                break

            if not is_T:
                ok_T = False
            if not is_Top:
                ok_Top = False

        if ok_mixed and all_odd:
            total_mixed += 2**nontrivial
        if ok_T and all_odd:
            total_T_only += 2**nontrivial
        if ok_Top and all_odd:
            total_Top_only += 2**nontrivial

    return total_mixed, total_T_only, total_Top_only

# Test on C_5
n = 5
cyc_edges = [(i, (i+1) % n) for i in range(n)] + [(i, (i+2) % n) for i in range(n)]

print(f"C_5 (H = 15):")
mixed, T_only, Top_only = verify_cycle_direction(n, cyc_edges)
print(f"  Mixed (T + T^op cycles): {mixed}")
print(f"  T-only cycles:           {T_only}")
print(f"  T^op-only cycles:        {Top_only}")
print(f"  Expected H:              15")

# Test on transitive T_5
trans_edges = [(i, j) for i in range(n) for j in range(i+1, n)]
print(f"\nTransitive T_5 (H = 1):")
mixed, T_only, Top_only = verify_cycle_direction(n, trans_edges)
print(f"  Mixed: {mixed}, T-only: {T_only}, T^op-only: {Top_only}")

# Test on Paley T_7
n = 7
QR = {1, 2, 4}
paley_edges = [(i, j) for i in range(n) for j in range(n)
               if i != j and (j - i) % n in QR]

print(f"\nPaley T_7 (H = 189):")
mixed, T_only, Top_only = verify_cycle_direction(n, paley_edges)
print(f"  Mixed: {mixed}, T-only: {T_only}, T^op-only: {Top_only}")

# Key question: does mixed = T_only always? Or does mixed > T_only?
# If mixed > T_only, then T^op cycles are separately counted and
# the factor of 2 partially comes from the T/T^op duality.

print(f"\n=== Factor 2 Analysis ===")
print(f"If mixed > T_only, the factor of 2 partially comes from T^op direction.")
print(f"If mixed = T_only, then the sum automatically excludes T^op.")

# Also check: what happens with JUST the identity permutation?
print(f"\n=== Identity contribution ===")
print(f"Identity permutation: all fixed points, 0 nontrivial cycles")
print(f"Contribution: 2^0 = 1")
print(f"This is the alpha_0 = 1 term in I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...")
