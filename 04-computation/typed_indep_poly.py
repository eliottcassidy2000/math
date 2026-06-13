#!/usr/bin/env python3
"""
The TYPED independence polynomial I_typed(Omega; y_3, y_5, y_7, ...).

This is the multivariate refinement of I(Omega, x) that tracks cycle LENGTHS.
GS proved: U_T = p_1^n * I_typed(Omega; 2p_3/p_1^3, 2p_5/p_1^5, ...)

Key questions:
1. I(Omega, x) = I_typed(Omega; x, x, x, ...) -- does this hold?
2. Does I_typed factor or have special structure?
3. Does the deletion T -> T-v respect the typed structure?
4. Does Claim A have a typed version?

NEW INSIGHT: alpha_0 = 1 is the empty independent set = p_1^n term in U_T.
This is the "transitive baseline" contribution to H(T).
Every tournament "starts with" H=1 (transitive case) and adds cycles.

The OCF says: H(T) = 1 + 2*alpha_1 + 4*alpha_2 + ...
In typed form: H(T) = 1 + 2*t3 + 2*t5 + 2*t7 + ... + 4*bc33 + 4*bc35 + ...
where t_k = #(k-cycles), bc_{jk} = #(disjoint j,k-cycle pairs), etc.

CLAIM A in typed form:
  H(T) - H(T-v) = 2 * sum_{C ni v} mu_typed(C)
where mu_typed(C) accounts for cycle lengths in the sub-problem.

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import factorial, comb

def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx] == 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def find_ALL_directed_odd_cycles(A, n, vertices=None):
    if vertices is None:
        vertices = list(range(n))
    all_cycles = []
    nv = len(vertices)
    for length in range(3, nv+1, 2):
        for subset in combinations(vertices, length):
            for perm in permutations(subset):
                is_cycle = True
                for i in range(length):
                    if A[perm[i]][perm[(i+1)%length]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    min_idx = list(perm).index(min(perm))
                    canon = perm[min_idx:] + perm[:min_idx]
                    all_cycles.append((canon, length))
    seen = set()
    result = []
    for c, l in all_cycles:
        if c not in seen:
            seen.add(c)
            result.append((c, l))
    return result

def conflict_graph_adj(cycles):
    nc = len(cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(cycles[i][0]) & set(cycles[j][0]):
                adj[i][j] = adj[j][i] = True
    return adj

def enumerate_independent_sets(adj, nc):
    result = []
    for mask in range(2**nc):
        nodes = [i for i in range(nc) if (mask >> i) & 1]
        ok = True
        for a in range(len(nodes)):
            for b in range(a+1, len(nodes)):
                if adj[nodes[a]][nodes[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            result.append(nodes)
    return result

def typed_indep_poly(cycles):
    nc = len(cycles)
    if nc == 0:
        return {(): 1}
    adj = conflict_graph_adj(cycles)
    indep = enumerate_independent_sets(adj, nc)
    result = defaultdict(int)
    for s in indep:
        lengths = tuple(sorted([cycles[i][1] for i in s], reverse=True))
        result[lengths] += 1
    return dict(result)

def H_from_typed(tip):
    """Compute H(T) from typed independence polynomial."""
    return sum(2**len(lengths) * count for lengths, count in tip.items())

def main():
    print("=" * 70)
    print("TYPED INDEPENDENCE POLYNOMIAL AND VERTEX DELETION")
    print("=" * 70)

    # Test at n=5: typed Claim A
    n = 5
    m = n*(n-1)//2
    print(f"\n--- n={n}: Typed vertex deletion ---")

    claim_a_typed_ok = 0
    claim_a_typed_fail = 0

    for bits_int in range(2**m):
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)

        cycles_T = find_ALL_directed_odd_cycles(A, n)
        tip_T = typed_indep_poly(cycles_T)
        H_T = H_from_typed(tip_T)

        for v in range(n):
            # T-v
            remaining = [u for u in range(n) if u != v]
            cycles_Tv = find_ALL_directed_odd_cycles(A, n, remaining)
            tip_Tv = typed_indep_poly(cycles_Tv)
            H_Tv = H_from_typed(tip_Tv)

            diff = H_T - H_Tv

            # Cycles through v
            cycles_v = [(c, l) for c, l in cycles_T if v in c]

            # For each cycle C through v, compute mu(C)
            # mu(C) = I(Omega(T-v)|_{avoid C\{v}}, 2)
            # where "avoid C\{v}" means only cycles vertex-disjoint from C\{v}
            rhs = 0
            for cv, cv_len in cycles_v:
                forbidden = set(cv) - {v}
                # Cycles in T-v that are vertex-disjoint from forbidden
                allowed_cycles = [(c, l) for c, l in cycles_Tv
                                  if not (set(c) & forbidden)]
                tip_allowed = typed_indep_poly(allowed_cycles)
                mu = H_from_typed(tip_allowed)
                rhs += mu

            rhs *= 2

            if diff == rhs:
                claim_a_typed_ok += 1
            else:
                claim_a_typed_fail += 1
                if claim_a_typed_fail <= 3:
                    print(f"  FAIL: bits={bits_int}, v={v}, diff={diff}, rhs={rhs}")

    print(f"  Claim A verified: {claim_a_typed_ok} OK, {claim_a_typed_fail} FAIL")

    # Now investigate the TYPED structure of the deletion
    print(f"\n--- Typed deletion analysis at n=5 ---")
    # Pick a specific tournament with interesting structure
    # bits=100 has 3 three-cycles and 1 five-cycle
    for bits_int in [100, 200, 920]:
        if bits_int >= 2**m:
            continue
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)

        cycles_T = find_ALL_directed_odd_cycles(A, n)
        tip_T = typed_indep_poly(cycles_T)
        H_T = H_from_typed(tip_T)

        print(f"\nbits={bits_int}, H={H_T}")
        print(f"  Typed IP: {tip_T}")

        for v in range(n):
            remaining = [u for u in range(n) if u != v]
            cycles_Tv = find_ALL_directed_odd_cycles(A, n, remaining)
            tip_Tv = typed_indep_poly(cycles_Tv)
            H_Tv = H_from_typed(tip_Tv)

            cycles_v = [(c, l) for c, l in cycles_T if v in c]
            cycle_types_v = [l for _, l in cycles_v]

            print(f"  v={v}: H(T-v)={H_Tv}, diff={H_T-H_Tv}, "
                  f"cycles_through_v={cycle_types_v}, "
                  f"tip(T-v)={tip_Tv}")

    # KEY ANALYSIS: Does the typed deletion have a pattern?
    print(f"\n\n--- alpha_0 = 1 interpretation ---")
    print("The empty independent set contributes 1 to H(T).")
    print("For transitive T: H=1, only alpha_0=1, all others=0.")
    print("Each odd cycle ADDED increases H by 2 (if no conflicts).")
    print("")
    print("In the deletion H(T) - H(T-v):")
    print("  alpha_0 contributions cancel: 1 - 1 = 0")
    print("  Only CYCLES through v contribute to the difference!")
    print("  This is why Claim A involves sum over C ni v.")
    print("")
    print("DEEPER: The empty set cancellation is EXACTLY why")
    print("Claim A doesn't have an alpha_0 term — it's the identity")
    print("  sum_{C ni v} mu(C) = (H(T) - 1)/2  only if v is 'generic'")
    print("  More precisely: H(T) - H(T-v) = 2 * sum_{C ni v} mu(C)")

if __name__ == "__main__":
    main()
