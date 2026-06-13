#!/usr/bin/env python3
"""overnight_s1_fixed.py -- Fixed computation with correct cycle counting.

Session: kind-pasteur-2026-03-20-S1

BUG FIX: The original independence_polynomial function only counted ONE directed
cycle per vertex set, but a vertex set of size k can have MULTIPLE directed
Hamiltonian cycles (e.g., at n=5 the full vertex set has 2 directed cycles:
0->1->2->3->4->0 and 0->2->4->1->3->0). Each directed cycle is a separate
vertex in Omega(T), and OCF requires ALL of them.

This script:
  Part A: Verify OCF with corrected cycle counting
  Part B: Paley vs Interval with correct I(Omega,2)
  Part C: SC Maximizer orbit analysis with correct I(Omega,2)
  Part D: Deep H/|Aut| factorization
  Part E: Forbidden H values at n=8,9
"""

import sys, os
from itertools import combinations, permutations
from math import comb, factorial, gcd
from collections import defaultdict, Counter

# ================================================================
# CORE UTILITIES
# ================================================================

def adj_from_bits(bits, n):
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj

def held_karp_H(adj, n):
    """Compute H(T) via Held-Karp DP. O(n^2 * 2^n)."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

def count_directed_odd_cycles(adj, n):
    """Count ALL directed odd cycles, properly.
    Returns list of (vertex_set_frozenset, cycle_tuple) pairs.
    Each directed Hamiltonian cycle on a vertex subset is a separate entry.
    """
    cycles = []
    for size in range(3, n+1, 2):
        for verts in combinations(range(n), size):
            sub_adj = [[adj[verts[i]][verts[j]] for j in range(size)] for i in range(size)]
            # Enumerate all directed Hamiltonian cycles starting from vertex 0
            for perm in permutations(range(1, size)):
                path = [0] + list(perm)
                is_cycle = True
                for idx in range(size):
                    if not sub_adj[path[idx]][path[(idx+1) % size]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # This is a valid directed cycle
                    actual_verts = tuple(verts[path[i]] for i in range(size))
                    cycles.append((frozenset(verts), actual_verts))
    return cycles

def independence_polynomial_correct(adj, n):
    """Compute independence polynomial of Omega(T) CORRECTLY.
    Each directed odd cycle is a vertex of Omega.
    Two cycles conflict iff they share a vertex.
    Returns [alpha_0, alpha_1, ...].
    """
    all_cycles = count_directed_odd_cycles(adj, n)
    num_cycles = len(all_cycles)

    if num_cycles == 0:
        return [1]

    # Vertex sets for conflict detection
    vsets = [c[0] for c in all_cycles]

    # Build conflict adjacency
    conflict = [[False]*num_cycles for _ in range(num_cycles)]
    for i in range(num_cycles):
        for j in range(i+1, num_cycles):
            if vsets[i] & vsets[j]:  # Shared vertex
                conflict[i][j] = True
                conflict[j][i] = True

    # Count independent sets by size
    alpha = [1]  # alpha_0 = 1
    for k in range(1, num_cycles + 1):
        count = 0
        for subset in combinations(range(num_cycles), k):
            independent = True
            for a in range(len(subset)):
                for b in range(a+1, len(subset)):
                    if conflict[subset[a]][subset[b]]:
                        independent = False
                        break
                if not independent:
                    break
            if independent:
                count += 1
        if count == 0:
            break
        alpha.append(count)
    return alpha

def eval_ip(alpha, x):
    return sum(a * x**k for k, a in enumerate(alpha))

def circulant_adj(n, conn_set):
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in conn_set:
            j = (i + s) % n
            adj[i][j] = 1
    return adj

def paley_tournament(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    return circulant_adj(p, qr)

def interval_tournament(p):
    S = set(range(1, (p+1)//2))
    return circulant_adj(p, S)

def is_self_converse(adj, n):
    """Check if tournament is self-converse."""
    opp = [[1 - adj[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if adj[perm[i]][perm[j]] != opp[i][j]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True
    return False

def is_score_self_complementary(score):
    n = len(score)
    s = sorted(score)
    return all(s[i] + s[n-1-i] == n-1 for i in range(n))

def prime_factors(n):
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


# ================================================================
# PART A: VERIFY OCF WITH CORRECT COUNTING
# ================================================================

def part_a_verify_ocf():
    print("=" * 70)
    print("PART A: VERIFY OCF WITH CORRECT CYCLE COUNTING")
    print("=" * 70)

    for n in [3, 4, 5, 6, 7]:
        m = n * (n - 1) // 2
        if n <= 6:
            sample = range(1 << m)
            label = "exhaustive"
        else:
            import random
            random.seed(42)
            sample = [random.randint(0, (1 << m) - 1) for _ in range(500)]
            label = "500 samples"

        violations = 0
        total = 0
        for bits in sample:
            adj = adj_from_bits(bits, n)
            H = held_karp_H(adj, n)
            ip = independence_polynomial_correct(adj, n)
            I2 = eval_ip(ip, 2)
            if H != I2:
                violations += 1
                if violations <= 3:
                    print(f"  VIOLATION at n={n}, bits={bits}: H={H}, I(Omega,2)={I2}, IP={ip}")
            total += 1

        print(f"  n={n} ({label}): {total} tested, {violations} violations")


# ================================================================
# PART B: PALEY vs INTERVAL
# ================================================================

def part_b_paley_vs_interval():
    print("\n" + "=" * 70)
    print("PART B: PALEY vs INTERVAL COMPARISON")
    print("=" * 70)

    for p in [3, 5, 7]:
        paley_adj = paley_tournament(p) if p % 4 == 3 else None
        interval_adj = interval_tournament(p)

        H_int = held_karp_H(interval_adj, p)
        ip_int = independence_polynomial_correct(interval_adj, p)

        print(f"\n  p={p} (p mod 4 = {p%4})")

        if paley_adj:
            H_pal = held_karp_H(paley_adj, p)
            ip_pal = independence_polynomial_correct(paley_adj, p)
            qr = sorted(set((x*x)%p for x in range(1,p)))

            print(f"    QR = {qr}")
            print(f"    Interval S = {list(range(1, (p+1)//2))}")
            print(f"    H(Paley)    = {H_pal}, I(Omega,2) = {eval_ip(ip_pal, 2)}, IP = {ip_pal}")
            print(f"    H(Interval) = {H_int}, I(Omega,2) = {eval_ip(ip_int, 2)}, IP = {ip_int}")

            # Count cycles by length
            pal_cycles = count_directed_odd_cycles(paley_adj, p)
            int_cycles = count_directed_odd_cycles(interval_adj, p)
            pal_by_len = Counter(len(c[0]) for c in pal_cycles)
            int_by_len = Counter(len(c[0]) for c in int_cycles)
            print(f"    Paley cycles by length: {dict(sorted(pal_by_len.items()))}")
            print(f"    Interval cycles by length: {dict(sorted(int_by_len.items()))}")

            if H_pal != H_int:
                winner = 'PALEY' if H_pal > H_int else 'INTERVAL'
                print(f"    Winner: {winner} (ratio: {H_pal/H_int:.6f})")
        else:
            print(f"    No Paley tournament (p = 1 mod 4)")
            print(f"    H(Interval) = {H_int}, I(Omega,2) = {eval_ip(ip_int, 2)}, IP = {ip_int}")


# ================================================================
# PART C: SC MAXIMIZER WITH CORRECT IP
# ================================================================

def part_c_sc_maximizer():
    print("\n" + "=" * 70)
    print("PART C: SC MAXIMIZER WITH CORRECT INDEPENDENCE POLYNOMIAL")
    print("=" * 70)

    n = 6
    m = n * (n - 1) // 2

    # Collect all tournaments with SC score sequences
    by_score = defaultdict(list)
    for bits in range(1 << m):
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        if not is_score_self_complementary(sc):
            continue
        H = held_karp_H(adj, n)
        is_sc = is_self_converse(adj, n)
        ip = independence_polynomial_correct(adj, n)
        by_score[sc].append({'bits': bits, 'H': H, 'is_sc': is_sc, 'IP': tuple(ip)})

    print(f"\n  n={n}, SC score classes: {len(by_score)}")

    for sc in sorted(by_score.keys()):
        entries = by_score[sc]
        sc_entries = [e for e in entries if e['is_sc']]
        nsc_entries = [e for e in entries if not e['is_sc']]

        sc_max_H = max((e['H'] for e in sc_entries), default=0)
        nsc_max_H = max((e['H'] for e in nsc_entries), default=0)

        print(f"\n  Score {sc}:")
        print(f"    SC: {len(sc_entries)} tours, max H = {sc_max_H}")
        print(f"    NSC: {len(nsc_entries)} tours, max H = {nsc_max_H}")

        if sc_entries and nsc_entries:
            gap = sc_max_H - nsc_max_H
            print(f"    Gap = {gap} ({'SC wins' if gap > 0 else 'NSC wins' if gap < 0 else 'TIE'})")

            # Analyze WHY: compare IP distributions
            sc_ips = Counter(e['IP'] for e in sc_entries)
            nsc_ips = Counter(e['IP'] for e in nsc_entries)

            # What IP achieves max H?
            sc_max_ip = [e['IP'] for e in sc_entries if e['H'] == sc_max_H]
            nsc_max_ip = [e['IP'] for e in nsc_entries if e['H'] == nsc_max_H] if nsc_entries else []

            print(f"    SC max IP: {sc_max_ip[0] if sc_max_ip else 'N/A'}")
            print(f"    NSC max IP: {nsc_max_ip[0] if nsc_max_ip else 'N/A'}")

            # Key insight: where does SC win? Higher alpha_1 or alpha_2?
            if sc_max_ip and nsc_max_ip:
                sc_ip = sc_max_ip[0]
                nsc_ip = nsc_max_ip[0]
                for k in range(max(len(sc_ip), len(nsc_ip))):
                    a_sc = sc_ip[k] if k < len(sc_ip) else 0
                    a_nsc = nsc_ip[k] if k < len(nsc_ip) else 0
                    if a_sc != a_nsc:
                        print(f"    Difference at alpha_{k}: SC={a_sc} vs NSC={a_nsc} (delta={a_sc-a_nsc})")


# ================================================================
# PART D: H/|Aut| DEEP FACTORIZATION
# ================================================================

def part_d_factorization():
    print("\n" + "=" * 70)
    print("PART D: H(T_p) / |Aut(T_p)| DEEP ANALYSIS")
    print("=" * 70)

    # Compute |Aut| for p=3,7 directly
    for p in [3, 7]:
        adj = paley_tournament(p)
        H = held_karp_H(adj, p)
        # Count automorphisms
        aut_count = 0
        for perm in permutations(range(p)):
            is_aut = True
            for i in range(p):
                for j in range(i+1, p):
                    if adj[perm[i]][perm[j]] != adj[i][j]:
                        is_aut = False
                        break
                if not is_aut:
                    break
            if is_aut:
                aut_count += 1
        ratio = H // aut_count
        print(f"\n  p={p}: H={H}, |Aut|={aut_count}, H/|Aut|={ratio}")
        print(f"    H factored: {' * '.join(str(f) for f in prime_factors(H))}")
        print(f"    |Aut| factored: {' * '.join(str(f) for f in prime_factors(aut_count))}")
        print(f"    H/|Aut| factored: {' * '.join(str(f) for f in prime_factors(ratio)) if ratio > 1 else '1'}")

    # Known values for p=11,19
    known = [
        (11, 95095, 55),
        (19, 1172695746915, 171),
    ]
    for p, H, aut_size in known:
        ratio = H // aut_size
        print(f"\n  p={p}: H={H}, |Aut|={aut_size}, H/|Aut|={ratio}")
        print(f"    H factored: {' * '.join(str(f) for f in prime_factors(H))}")
        print(f"    |Aut| factored: {' * '.join(str(f) for f in prime_factors(aut_size))}")
        print(f"    H/|Aut| factored: {' * '.join(str(f) for f in prime_factors(ratio))}")

    # Look for arithmetic patterns
    print("\n  --- Arithmetic patterns ---")
    ratios = [1, 9, 1729, 6857869865]
    primes_list = [3, 7, 11, 19]
    H_vals = [3, 189, 95095, 1172695746915]
    aut_vals = [3, 21, 55, 171]

    for i, (p, H, aut, r) in enumerate(zip(primes_list, H_vals, aut_vals, ratios)):
        # |Aut(T_p)| = p * (p-1)/2 for Paley primes
        expected_aut = p * (p - 1) // 2
        print(f"\n  p={p}:")
        print(f"    |Aut| = {aut}, p*(p-1)/2 = {expected_aut}, match: {aut == expected_aut}")
        print(f"    H/(p!) = {H / factorial(p):.6e}")
        print(f"    H / (p!/2^(p-1)) = {H * 2**(p-1) / factorial(p):.6f}")
        print(f"    H = sum 2^k * alpha_k (by OCF)")

    # Check if 1729 = sum of cubes
    print(f"\n  1729 = 12^3 + 1^3 = {12**3 + 1**3}")
    print(f"  1729 = 10^3 + 9^3 = {10**3 + 9**3}")
    print(f"  1729 = 7 * 13 * 19 = {7*13*19}")
    print(f"  These primes are: p, p+2, p+6 for p=7? No, p=7 and we need mod structure.")
    print(f"  7 = (p-4)/1, 13 = (p+2)/1, 19 = (p+8)/1 for p=11")
    print(f"  Or: 7*13*19 = C(20,3) - 11? No, C(20,3) = {comb(20,3)}")
    print(f"  Actually: p=11, primes 7,13,19 are p-4, p+2, p+8 = 2k-1 for k=4,7,10 (AP!)")
    print(f"  Arithmetic progression: 4, 7, 10 (common diff 3 = (p-2)/3)")

    # Check if 6857869865 has similar structure
    r19 = 6857869865
    f19 = prime_factors(r19)
    print(f"\n  H/|Aut| for p=19: {r19} = {' * '.join(str(f) for f in f19)}")
    print(f"  Unique primes: {sorted(set(f19))}")
    # Check if primes relate to p=19 arithmetically
    up19 = sorted(set(f19))
    print(f"  Differences from p=19: {[q - 19 for q in up19]}")
    print(f"  Differences from each other: {[up19[i+1]-up19[i] for i in range(len(up19)-1)]}")


# ================================================================
# PART E: FORBIDDEN VALUES n=8
# ================================================================

def part_e_forbidden():
    print("\n" + "=" * 70)
    print("PART E: FORBIDDEN H VALUES AT n=8 (200K samples)")
    print("=" * 70)

    import random
    random.seed(2026)

    n = 8
    m = n * (n - 1) // 2  # 28
    sample_size = 200000

    H_values = Counter()
    for trial in range(sample_size):
        bits = random.randint(0, (1 << m) - 1)
        adj = adj_from_bits(bits, n)
        H = held_karp_H(adj, n)
        H_values[H] += 1

    all_H = sorted(H_values.keys())
    max_H = max(all_H)

    print(f"  {sample_size} random n={n} tournaments")
    print(f"  H range: {min(all_H)} to {max_H}")
    print(f"  Distinct H: {len(all_H)}")

    expected = set(range(1, max_H + 1, 2))
    achieved = set(all_H)
    gaps = sorted(expected - achieved)
    print(f"  Gaps in odd [1,{max_H}]: {len(gaps)}")

    # Separate permanent gaps (small) from sampling gaps (large)
    print(f"\n  Gaps below H=100: {[g for g in gaps if g < 100]}")
    print(f"  H=7: {H_values.get(7, 0)}, H=21: {H_values.get(21, 0)}")
    print(f"  H=63: {H_values.get(63, 0)}")  # known to be achievable at n=8

    # Distribution
    print(f"\n  H value counts (top 20):")
    for h, c in H_values.most_common(20):
        print(f"    H={h}: {c}")


# ================================================================
# PART F: SC MAXIMIZER — WHY DOES IT WORK?
# ================================================================

def part_f_sc_why():
    print("\n" + "=" * 70)
    print("PART F: SC MAXIMIZER — ALGEBRAIC ANALYSIS")
    print("=" * 70)

    # For n=6, analyze the SC maximizer mechanism in detail
    n = 6
    m = n * (n - 1) // 2

    # Focus on the most interesting score class: (2,2,2,3,3,3) — the regular one
    target_score = (2, 2, 2, 3, 3, 3)

    sc_tours = []
    nsc_tours = []

    for bits in range(1 << m):
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        if sc != target_score:
            continue

        H = held_karp_H(adj, n)
        is_sc = is_self_converse(adj, n)
        all_cycles = count_directed_odd_cycles(adj, n)
        by_len = Counter(len(c[0]) for c in all_cycles)

        # Count disjoint cycle pairs
        cycle_vsets = [c[0] for c in all_cycles]
        disjoint_pairs = 0
        for i in range(len(cycle_vsets)):
            for j in range(i+1, len(cycle_vsets)):
                if not (cycle_vsets[i] & cycle_vsets[j]):
                    disjoint_pairs += 1

        entry = {
            'bits': bits, 'H': H, 'is_sc': is_sc,
            'num_cycles': len(all_cycles),
            'c3': by_len.get(3, 0), 'c5': by_len.get(5, 0),
            'disjoint_pairs': disjoint_pairs,
        }

        if is_sc:
            sc_tours.append(entry)
        else:
            nsc_tours.append(entry)

    print(f"\n  Score {target_score} at n={n}:")
    print(f"    SC: {len(sc_tours)}, NSC: {len(nsc_tours)}")

    # Analyze H distribution
    sc_H = Counter(e['H'] for e in sc_tours)
    nsc_H = Counter(e['H'] for e in nsc_tours)
    print(f"    SC H values: {dict(sorted(sc_H.items()))}")
    print(f"    NSC H values: {dict(sorted(nsc_H.items()))}")

    # Compare cycle structure at max H
    sc_max = max(sc_tours, key=lambda e: e['H'])
    nsc_max = max(nsc_tours, key=lambda e: e['H']) if nsc_tours else None

    print(f"\n    SC max: H={sc_max['H']}, #cycles={sc_max['num_cycles']}, "
          f"c3={sc_max['c3']}, c5={sc_max['c5']}, disjoint_pairs={sc_max['disjoint_pairs']}")
    if nsc_max:
        print(f"    NSC max: H={nsc_max['H']}, #cycles={nsc_max['num_cycles']}, "
              f"c3={nsc_max['c3']}, c5={nsc_max['c5']}, disjoint_pairs={nsc_max['disjoint_pairs']}")

    # The KEY question: does SC always have MORE disjoint pairs?
    sc_dp = [e['disjoint_pairs'] for e in sc_tours]
    nsc_dp = [e['disjoint_pairs'] for e in nsc_tours]
    print(f"\n    SC disjoint pairs: min={min(sc_dp)}, max={max(sc_dp)}, avg={sum(sc_dp)/len(sc_dp):.1f}")
    if nsc_dp:
        print(f"    NSC disjoint pairs: min={min(nsc_dp)}, max={max(nsc_dp)}, avg={sum(nsc_dp)/len(nsc_dp):.1f}")

    # Correlation between H and disjoint_pairs
    print(f"\n    Correlation analysis (SC):")
    for h_val in sorted(sc_H.keys()):
        entries = [e for e in sc_tours if e['H'] == h_val]
        avg_dp = sum(e['disjoint_pairs'] for e in entries) / len(entries)
        avg_cycles = sum(e['num_cycles'] for e in entries) / len(entries)
        print(f"      H={h_val}: avg_disjoint_pairs={avg_dp:.1f}, avg_total_cycles={avg_cycles:.1f}")

    print(f"\n    Correlation analysis (NSC):")
    for h_val in sorted(nsc_H.keys()):
        entries = [e for e in nsc_tours if e['H'] == h_val]
        avg_dp = sum(e['disjoint_pairs'] for e in entries) / len(entries)
        avg_cycles = sum(e['num_cycles'] for e in entries) / len(entries)
        print(f"      H={h_val}: avg_disjoint_pairs={avg_dp:.1f}, avg_total_cycles={avg_cycles:.1f}")


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":
    part_a_verify_ocf()
    part_b_paley_vs_interval()
    part_c_sc_maximizer()
    part_d_factorization()
    part_e_forbidden()
    part_f_sc_why()
