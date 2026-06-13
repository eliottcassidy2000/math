#!/usr/bin/env python3
"""
FIXED: Exact GS U_T <-> OCF correspondence.
Previous version had bug: only found ONE directed cycle per vertex subset.
A tournament can have MULTIPLE directed Hamiltonian cycles on the same vertices.

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

def find_ALL_directed_odd_cycles(A, n):
    """Find ALL directed odd cycles (including multiple on same vertex set)."""
    all_cycles = []
    for length in range(3, n+1, 2):
        for subset in combinations(range(n), length):
            # Find ALL directed cycles on this subset
            for perm in permutations(subset):
                is_cycle = True
                for i in range(length):
                    if A[perm[i]][perm[(i+1)%length]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    # Canonicalize: start from smallest vertex
                    min_idx = list(perm).index(min(perm))
                    canon = perm[min_idx:] + perm[:min_idx]
                    all_cycles.append((canon, length))
    # Deduplicate
    seen = set()
    result = []
    for c, l in all_cycles:
        if c not in seen:
            seen.add(c)
            result.append((c, l))
    return result

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
    print("EXACT GS <-> OCF (FIXED CYCLE COUNTING)")
    print("=" * 70)

    # n=5: exhaustive
    n = 5
    m = n*(n-1)//2
    print(f"\n--- Exhaustive n={n} ({2**m} tournaments) ---")

    patterns = defaultdict(list)
    all_match = True

    for bits_int in range(2**m):
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)
        ut = compute_UT(A, n)

        cycles = find_ALL_directed_odd_cycles(A, n)
        nc = len(cycles)

        if nc == 0:
            typed_ip = {(): 1}
        else:
            adj = conflict_graph_adj(cycles)
            indep_sets = enumerate_independent_sets(adj, nc)
            typed_ip = defaultdict(int)
            for s in indep_sets:
                lengths = tuple(sorted([cycles[i][1] for i in s], reverse=True))
                typed_ip[lengths] += 1
            typed_ip = dict(typed_ip)

        # Count cycles by length
        t3 = sum(1 for c, l in cycles if l == 3)
        t5 = sum(1 for c, l in cycles if l == 5)

        c_311 = ut.get((3,1,1), 0)
        c_5 = ut.get((5,), 0)

        key = (c_311, c_5, t3, t5)
        patterns[key].append(bits_int)

    print(f"\nFound {len(patterns)} distinct (U_T, cycle-count) patterns:")
    for key in sorted(patterns.keys()):
        c_31, c_5_val, t3, t5 = key
        count = len(patterns[key])

        # Each directed k-cycle contributes 2 to coeff of p_k
        # (forward traversal with (-1)^{k-1} = +1 for odd k,
        #  plus backward traversal in T^op with (-1)^0 = +1)
        expected_p3 = 2 * t3
        expected_p5 = 2 * t5

        match_3 = (c_31 == expected_p3)
        match_5 = (c_5_val == expected_p5)

        if not (match_3 and match_5):
            all_match = False

        print(f"  c_p3p1^2={c_31}, c_p5={c_5_val}, #3cyc={t3}, #5cyc={t5}, "
              f"count={count} "
              f"{'OK' if match_3 and match_5 else 'MISMATCH'}"
              f"{'' if match_3 else f' [p3: expected {expected_p3}]'}"
              f"{'' if match_5 else f' [p5: expected {expected_p5}]'}")

    # For the mismatching cases, investigate: what are the 5-cycles?
    print("\n\n--- INVESTIGATING MISMATCHES ---")
    for key in sorted(patterns.keys()):
        c_31, c_5_val, t3, t5 = key
        if c_5_val != 2 * t5:
            bits_int = patterns[key][0]
            b = [(bits_int >> k) & 1 for k in range(m)]
            A = tournament_from_bits(n, b)
            cycles = find_ALL_directed_odd_cycles(A, n)

            print(f"\nbits={bits_int}, c_p5={c_5_val}, #5-cycles={t5}")
            print(f"  All cycles: {[(c,l) for c,l in cycles]}")

            # What permutations contribute to p_5?
            T_bar = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    if i != j:
                        T_bar[i][j] = 1 - A[i][j]

            p5_perms = []
            for perm in permutations(range(n)):
                # Check if single 5-cycle
                visited = [False]*n
                cycle = []
                j = 0
                while not visited[j]:
                    visited[j] = True
                    cycle.append(j)
                    j = perm[j]
                if len(cycle) != n:
                    continue
                # It's a single n-cycle
                k = n
                is_T = all(A[cycle[i]][cycle[(i+1)%k]] == 1 for i in range(k))
                is_Tbar = all(T_bar[cycle[i]][cycle[(i+1)%k]] == 1 for i in range(k))
                if is_T or is_Tbar:
                    phi_val = (k-1) if is_T else 0
                    sign = (-1)**phi_val
                    canon_start = cycle.index(min(cycle))
                    canon = tuple(cycle[canon_start:] + cycle[:canon_start])
                    p5_perms.append((canon, 'T' if is_T else 'Tbar', sign))

            print(f"  Contributing 5-cycles to p_5:")
            seen = set()
            for canon, typ, sign in p5_perms:
                if canon not in seen:
                    seen.add(canon)
                    print(f"    {canon} ({typ}, sign={sign})")

    if all_match:
        print("\n\nALL MATCH at n=5")
    else:
        print("\n\nSome mismatches — investigating T^op cycles")

    # Also check: directed 5-cycles vs Hamiltonian cycles
    print("\n\n--- HAMILTONIAN DIRECTED CYCLES AT n=5 ---")
    max_ham = 0
    for bits_int in range(2**m):
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = tournament_from_bits(n, b)
        # Count Hamiltonian directed cycles
        count = 0
        for perm in permutations(range(n)):
            if all(A[perm[i]][perm[(i+1)%n]] == 1 for i in range(n)):
                count += 1
        # Each cycle counted n times (starting positions)
        count //= n
        if count > max_ham:
            max_ham = count
            best = bits_int

    print(f"Max Hamiltonian directed cycles: {max_ham} (bits={best})")

    # Check specific case
    b = [(best >> k) & 1 for k in range(m)]
    A = tournament_from_bits(n, b)
    ut = compute_UT(A, n)
    cycles = find_ALL_directed_odd_cycles(A, n)
    t5 = sum(1 for c, l in cycles if l == 5)
    print(f"  U_T p_5 coeff = {ut.get((5,), 0)}")
    print(f"  #5-cycles found = {t5}")
    print(f"  Expected (2 * #5-cycles) = {2*t5}")
    print(f"  Match: {ut.get((5,), 0) == 2*t5}")

if __name__ == "__main__":
    main()
