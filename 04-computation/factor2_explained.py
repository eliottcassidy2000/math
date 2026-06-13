#!/usr/bin/env python3
"""
DISCOVERY: The factor of 2 in the OCF is explained by T/T^op duality in U_T.

THEOREM: ps_1(U_T)(1) = H(T) for any tournament T.

PROOF SKETCH:
U_T = sum_{sigma in S(T, T^op)} (-1)^phi(sigma) * p_{type(sigma)}

where S(T, T^op) = permutations whose nontrivial cycles are directed cycles
of T or T^op, and phi(sigma) = sum_{gamma a T-cycle} (len(gamma) - 1).

Under ps_1(p_lambda)(1) = 1 for all lambda, we get:
  ps_1(U_T)(1) = sum_{sigma in S(T, T^op)} (-1)^phi(sigma)

This sum splits into contributions by cycle type:

1. EVEN-LENGTH CYCLES CANCEL:
   A k-cycle (k even) in T-direction: phi += k-1 (odd), sign = -1.
   Same cycle in T^op-direction: phi += 0, sign = +1.
   These two cancel: -1 + 1 = 0.
   (Verified: at n=5, (4,1)-type has +5 from T^op and -5 from T.)

2. ODD-LENGTH CYCLES BOTH CONTRIBUTE +1:
   A k-cycle (k odd) in T-direction: phi += k-1 (even), sign = +1.
   Same cycle in T^op-direction: phi += 0, sign = +1.
   These two ADD: +1 + 1 = 2.
   (Verified: at n=5, (3,1,1)-type has +5 from T and +5 from T^op = 10 total.)

3. MIXED-DIRECTION permutations (some T-cycles, some T^op-cycles):
   Each odd T-cycle contributes +1 to phi (net sign +1).
   Each odd T^op-cycle contributes 0 to phi (net sign +1).
   Each even T-cycle contributes -1 (odd phi increment).
   Each even T^op-cycle contributes +1.
   Even cycles in a mixed perm still cancel pairwise.

CONCLUSION:
  ps_1(U_T)(1) = sum over collections of disjoint ODD cycles
                  (each either T or T^op direction) of 1^|collection|
                = sum over T-direction collections of 2^|collection|
                  (each cycle counted twice: once as T, once as T^op)
                = I(Omega(T), 2) = H(T)

The factor of 2 in the OCF (2^psi) comes from EXACTLY the T/T^op duality:
each odd directed cycle of T has a twin in T^op, and both contribute
with the same sign (+1) to ps_1(U_T)(1).

VERIFIED: n=3 (cyclic, transitive), n=5 (C_5, transitive).

kind-pasteur-2026-03-07-S30
"""

# Verification (same as inline computation above, cleaned up)
from itertools import permutations
from collections import defaultdict


def verify_factor2(n, A, name, expected_H):
    edge_set = set()
    opp_set = set()
    for i in range(n):
        for j in range(n):
            if i != j:
                if A[i][j]:
                    edge_set.add((i, j))
                else:
                    opp_set.add((i, j))

    # Categorize by (cycle type, #T-cycles, #Top-cycles)
    by_direction = defaultdict(int)
    total_signed = 0

    for sigma in permutations(range(n)):
        visited = [False] * n
        cycles = []
        for s in range(n):
            if visited[s]:
                continue
            cyc = []
            c = s
            while not visited[c]:
                visited[c] = True
                cyc.append(c)
                c = sigma[c]
            cycles.append(tuple(cyc))

        valid = True
        n_T, n_Top = 0, 0
        phi = 0
        has_even = False
        for cyc in cycles:
            if len(cyc) == 1:
                continue
            is_T = all((cyc[i], cyc[(i + 1) % len(cyc)]) in edge_set for i in range(len(cyc)))
            is_Top = all((cyc[i], cyc[(i + 1) % len(cyc)]) in opp_set for i in range(len(cyc)))
            if not is_T and not is_Top:
                valid = False
                break
            if is_T:
                n_T += 1
                phi += len(cyc) - 1
            else:
                n_Top += 1
            if len(cyc) % 2 == 0:
                has_even = True

        if valid:
            sign = (-1) ** phi
            ctype = tuple(sorted([len(c) for c in cycles], reverse=True))
            by_direction[(ctype, n_T, n_Top)] += sign
            total_signed += sign

    print(f"\n=== {name} (n={n}, expected H={expected_H}) ===")
    print(f"ps_1(U_T)(1) = {total_signed}")
    print(f"Match: {total_signed == expected_H}")

    # Show even-cycle cancellation
    even_types = {}
    for (ctype, nT, nTop), cnt in by_direction.items():
        if any(p % 2 == 0 for p in ctype):
            if ctype not in even_types:
                even_types[ctype] = 0
            even_types[ctype] += cnt
    if even_types:
        print(f"Even-cycle type cancellations:")
        for ct, total in sorted(even_types.items()):
            print(f"  {ct}: net = {total} (should be 0)")

    # Show odd-cycle doubling
    odd_type_sums = defaultdict(lambda: [0, 0])
    for (ctype, nT, nTop), cnt in by_direction.items():
        if all(p % 2 == 1 for p in ctype) and cnt != 0:
            if nT > 0 and nTop == 0:
                odd_type_sums[ctype][0] += cnt
            elif nTop > 0 and nT == 0:
                odd_type_sums[ctype][1] += cnt
    for ct in sorted(odd_type_sums.keys()):
        T_cnt, Top_cnt = odd_type_sums[ct]
        if T_cnt != 0 or Top_cnt != 0:
            print(f"  {ct}: T-only={T_cnt}, T^op-only={Top_cnt}, ratio={Top_cnt}/{T_cnt if T_cnt else '?'}")


# n=3 cyclic
A3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
verify_factor2(3, A3, "Cyclic C_3", 3)

# n=3 transitive
A3t = [[0, 1, 1], [0, 0, 1], [0, 0, 0]]
verify_factor2(3, A3t, "Transitive T_3", 1)

# n=5 cyclic
n = 5
A5 = [[0] * n for _ in range(n)]
for i in range(n):
    A5[i][(i + 1) % n] = 1
    A5[i][(i + 2) % n] = 1
verify_factor2(5, A5, "Cyclic C_5", 15)

# n=5 transitive
A5t = [[0] * n for _ in range(n)]
for i in range(n):
    for j in range(i + 1, n):
        A5t[i][j] = 1
verify_factor2(5, A5t, "Transitive T_5", 1)
