#!/usr/bin/env python3
"""
Exact correspondence: GS U_T coefficients vs OCF independence polynomial.

At n=3: U_T = p_1^3 + alpha_1 * (2p_3)   -- VERIFIED
At n=5: U_T = p_1^5 + ? * (2p_3)*p_1^2 + ? * (2p_5)

Question: are the coefficients of (2p_3)^a * (2p_5)^b * p_1^c the independence
numbers alpha_k, counting independent sets of size a+b in Omega(T)?

More precisely: in GS's {p_1, 2p_3, 2p_5, ...} polynomial, the monomial
  (2p_3)^{a_3} (2p_5)^{a_5} ... p_1^{n - 3a_3 - 5a_5 - ...}
should correspond to independent sets of odd cycles of sizes (3^{a_3}, 5^{a_5}, ...).

opus-2026-03-07-S34
"""
from itertools import permutations, combinations
from collections import Counter, defaultdict
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

def compute_UT(A, n):
    """Compute U_T expansion in power-sum basis."""
    T_bar = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                T_bar[i][j] = 1 - A[i][j]

    result = defaultdict(int)
    for perm in permutations(range(n)):
        visited = [False]*n
        cycles = []
        for i in range(n):
            if not visited[i]:
                cycle = []
                j = i
                while not visited[j]:
                    visited[j] = True
                    cycle.append(j)
                    j = perm[j]
                cycles.append(tuple(cycle))

        valid = True
        phi = 0
        for cyc in cycles:
            k = len(cyc)
            if k == 1:
                continue
            is_T = all(A[cyc[i]][cyc[(i+1)%k]] == 1 for i in range(k))
            is_Tbar = all(T_bar[cyc[i]][cyc[(i+1)%k]] == 1 for i in range(k))
            if is_T:
                phi += k - 1
            elif is_Tbar:
                pass
            else:
                valid = False
                break

        if valid:
            cycle_type = tuple(sorted([len(c) for c in cycles], reverse=True))
            result[cycle_type] += (-1)**phi

    return dict(result)

def find_directed_odd_cycles(A, n):
    """Find all directed odd cycles."""
    all_cycles = []
    for length in range(3, n+1, 2):
        for subset in combinations(range(n), length):
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
                    break
    # Deduplicate
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

def main():
    print("=" * 70)
    print("EXACT GS <-> OCF COEFFICIENT CORRESPONDENCE")
    print("=" * 70)

    # n=5: exhaustive over all 1024 tournaments
    n = 5
    m = n*(n-1)//2
    print(f"\n--- Exhaustive n={n} ({2**m} tournaments) ---")

    # Group by (U_T partition coeffs) to find patterns
    # For each tournament, record:
    # - U_T coefficients for each partition
    # - Independence polynomial of Omega(T)
    # - Cycle type decomposition of independent sets

    patterns = defaultdict(list)

    for bits_int in range(2**m):
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)
        ut = compute_UT(A, n)

        # Get odd cycles with their lengths
        cycles = find_directed_odd_cycles(A, n)
        nc = len(cycles)

        if nc == 0:
            ip = {0: 1}
            typed_ip = {(): 1}
        else:
            adj = conflict_graph_adj(cycles)
            indep_sets = enumerate_independent_sets(adj, nc)
            ip = defaultdict(int)
            typed_ip = defaultdict(int)  # partition of cycle lengths -> count
            for s in indep_sets:
                ip[len(s)] += 1
                lengths = tuple(sorted([cycles[i][1] for i in s], reverse=True))
                typed_ip[lengths] += 1
            ip = dict(ip)
            typed_ip = dict(typed_ip)

        H = sum(2**k * ip.get(k, 0) for k in ip)

        # Extract U_T coefficients in the relevant basis
        # For n=5, relevant partitions: (1,1,1,1,1), (3,1,1), (5)
        c_11111 = ut.get((1,1,1,1,1), 0)
        c_311 = ut.get((3,1,1), 0)
        c_5 = ut.get((5,), 0)

        key = (c_11111, c_311, c_5)
        patterns[key].append({
            'bits': bits_int,
            'ip': ip,
            'typed_ip': typed_ip,
            'H': H,
            'ut': ut,
        })

    print(f"\nFound {len(patterns)} distinct U_T coefficient patterns:")
    for key in sorted(patterns.keys()):
        entries = patterns[key]
        rep = entries[0]
        c_111, c_31, c_5_val = key

        # Check: is c_31 = 2 * (number of 3-cycles)?
        # And c_5_val = 2 * (number of 5-cycles)?
        t3 = rep['typed_ip'].get((3,), 0)
        t5 = rep['typed_ip'].get((5,), 0)
        # Independent sets with two 3-cycles:
        t33 = rep['typed_ip'].get((3,3), 0)

        alpha_0 = rep['ip'].get(0, 0)
        alpha_1 = rep['ip'].get(1, 0)
        alpha_2 = rep['ip'].get(2, 0)

        print(f"\n  U_T = {c_111}*p1^5 + {c_31}*p3*p1^2 + {c_5_val}*p5")
        print(f"    Count: {len(entries)} tournaments")
        print(f"    H = {rep['H']}")
        print(f"    I(Omega,x) = {rep['ip']}")
        print(f"    Typed IP: {rep['typed_ip']}")
        print(f"    In (p1, 2p3, 2p5) basis: {c_111}*p1^5 + {c_31//2}*(2p3)*p1^2 + {c_5_val//2}*(2p5)")
        print(f"    alpha_1={alpha_1} = #3cycles({t3}) + #5cycles({t5})")
        print(f"    alpha_2={alpha_2} = #(3,3)-pairs({t33})")

        # KEY TEST: does c_31/2 = t3 and c_5/2 = t5?
        if c_31 == 2*t3 and c_5_val == 2*t5:
            print(f"    *** MATCH: coeff(2p3*p1^2) = #3cycles, coeff(2p5) = #5cycles ***")
        else:
            print(f"    MISMATCH: c_31/2={c_31//2} vs t3={t3}, c_5/2={c_5_val//2} vs t5={t5}")

    # Summary statistics
    print("\n\n--- SUMMARY ---")
    all_match = True
    for key in sorted(patterns.keys()):
        entries = patterns[key]
        rep = entries[0]
        c_111, c_31, c_5_val = key
        t3 = rep['typed_ip'].get((3,), 0)
        t5 = rep['typed_ip'].get((5,), 0)
        if c_31 != 2*t3 or c_5_val != 2*t5:
            print(f"  MISMATCH at {key}")
            all_match = False

    if all_match:
        print("ALL MATCH: U_T coefficients exactly encode cycle-type counts!")
        print("  coeff of (2p_k)*p_1^{n-k} = number of directed k-cycles")
        print("  NOT the same as alpha_k (which sums over ALL cycle lengths)")
    else:
        print("Some mismatches found")

    # Deeper: what about the (2p_3)^2 term?
    # For n=5, (2p_3)^2 would need 6 vertices, so impossible. But let's check at n=7
    print("\n\n--- INTERPRETATION ---")
    print("The GS expansion U_T = sum c_lambda * p_lambda")
    print("separates cycle-type information that our OCF merges into alpha_k.")
    print("Specifically:")
    print("  - Our alpha_1 = #(3-cycles) + #(5-cycles) + #(7-cycles) + ...")
    print("  - GS gives SEPARATE coefficients for each cycle length")
    print("  - This is FINER than the independence polynomial I(Omega, x)")
    print("  - GS's expansion is the TYPED independence polynomial!")

if __name__ == "__main__":
    main()
