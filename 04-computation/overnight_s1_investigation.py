#!/usr/bin/env python3
"""overnight_s1_investigation.py — Deep investigation of core open questions.

Session: kind-pasteur-2026-03-20-S1

Investigations:
  Part 1: SC Maximizer — algebraic structure analysis
  Part 2: Paley vs Interval H comparison at p=7,11
  Part 3: H(T_p)/|Aut| factorization
  Part 4: Forbidden H values at n=9 (sampling)
  Part 5: SC Maximizer orbit-pairing analysis
"""

import sys, os
from itertools import combinations, permutations
from math import comb, factorial, gcd
from collections import defaultdict, Counter
from functools import lru_cache

# ================================================================
# CORE TOURNAMENT UTILITIES
# ================================================================

def adj_from_bits(bits, n):
    """Create adjacency matrix from bit encoding."""
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

def bits_from_adj(adj, n):
    """Get bit encoding from adjacency matrix."""
    bits = 0
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j]:
                bits |= (1 << k)
            k += 1
    return bits

def held_karp_H(adj, n):
    """Compute H(T) via Held-Karp DP."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def score_seq(adj, n):
    """Return sorted score sequence."""
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

def count_3cycles(adj, n):
    """Count directed 3-cycles."""
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check both orientations
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1
    return c3

def count_odd_cycles(adj, n):
    """Count all directed odd cycles by vertex set, return (c3, c5, c7, ...)."""
    counts = {}
    for size in range(3, n+1, 2):
        count = 0
        for verts in combinations(range(n), size):
            # Count directed Hamiltonian cycles on this vertex set
            sub_adj = [[adj[verts[i]][verts[j]] for j in range(size)] for i in range(size)]
            # Count directed Hamiltonian cycles using permutations
            # A directed Hamiltonian cycle visits all vertices exactly once
            for perm in permutations(range(1, size)):
                path = [0] + list(perm)
                is_cycle = True
                for idx in range(size):
                    if not sub_adj[path[idx]][path[(idx+1) % size]]:
                        is_cycle = False
                        break
                if is_cycle:
                    count += 1
            # Each cycle counted size times (once per starting vertex)
        counts[size] = count // size
    return counts

def independence_polynomial(adj, n):
    """Compute independence polynomial of conflict graph Omega(T).
    Returns list [alpha_0, alpha_1, alpha_2, ...]"""
    # First find all directed odd cycles
    odd_cycles = []
    for size in range(3, n+1, 2):
        for verts in combinations(range(n), size):
            sub_adj = [[adj[verts[i]][verts[j]] for j in range(size)] for i in range(size)]
            for perm in permutations(range(1, size)):
                path = [0] + list(perm)
                is_cycle = True
                for idx in range(size):
                    if not sub_adj[path[idx]][path[(idx+1) % size]]:
                        is_cycle = False
                        break
                if is_cycle:
                    odd_cycles.append(set(verts))
                    break  # Only need one per vertex set

    # Remove duplicate vertex sets
    unique_cycles = []
    seen = set()
    for c in odd_cycles:
        key = frozenset(c)
        if key not in seen:
            seen.add(key)
            unique_cycles.append(c)

    # Build conflict graph (two cycles conflict if they share a vertex)
    num_cycles = len(unique_cycles)
    conflict = [[False]*num_cycles for _ in range(num_cycles)]
    for i in range(num_cycles):
        for j in range(i+1, num_cycles):
            if unique_cycles[i] & unique_cycles[j]:
                conflict[i][j] = True
                conflict[j][i] = True

    # Count independent sets by size using inclusion
    alpha = [1]  # alpha_0 = 1
    for k in range(1, num_cycles + 1):
        count = 0
        for subset in combinations(range(num_cycles), k):
            # Check if all pairs are non-conflicting
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

def eval_indpoly(alpha, x):
    """Evaluate independence polynomial at x."""
    return sum(a * x**k for k, a in enumerate(alpha))

def is_self_converse(adj, n):
    """Check if tournament is self-converse (isomorphic to its opposite)."""
    # Check all permutations (only feasible for small n)
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

def find_anti_auts(adj, n):
    """Find all anti-automorphisms of tournament."""
    anti_auts = []
    for perm in permutations(range(n)):
        is_anti = True
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]] != adj[j][i]:
                    is_anti = False
                    break
            if not is_anti:
                break
        if is_anti:
            anti_auts.append(perm)
    return anti_auts

def find_automorphisms(adj, n):
    """Find all automorphisms of tournament."""
    auts = []
    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]] != adj[i][j]:
                    is_aut = False
                    break
            if not is_aut:
                break
        if is_aut:
            auts.append(perm)
    return auts

def is_score_self_complementary(score):
    """Check if a score sequence is self-complementary."""
    n = len(score)
    s = sorted(score)
    return all(s[i] + s[n-1-i] == n-1 for i in range(n))

# ================================================================
# CIRCULANT TOURNAMENTS
# ================================================================

def circulant_adj(n, conn_set):
    """Create circulant tournament on Z_n with connection set S."""
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in conn_set:
            j = (i + s) % n
            adj[i][j] = 1
    return adj

def paley_tournament(p):
    """Create Paley tournament for prime p = 3 mod 4."""
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    return circulant_adj(p, qr)

def interval_tournament(p):
    """Create interval (cyclic) tournament on Z_p with S = {1, ..., (p-1)/2}."""
    S = set(range(1, (p+1)//2))
    return circulant_adj(p, S)


# ================================================================
# PART 1: SC MAXIMIZER DEEP ANALYSIS — n=5,6,7
# ================================================================

def part1_sc_maximizer():
    print("=" * 70)
    print("PART 1: SC MAXIMIZER — DEEP STRUCTURE ANALYSIS")
    print("=" * 70)

    for n in [5, 6, 7]:
        m = n * (n - 1) // 2
        total = 1 << m

        if n == 7:
            # Too many - sample
            import random
            random.seed(42)
            sample_size = 50000
            tournaments = [random.randint(0, total - 1) for _ in range(sample_size)]
        else:
            tournaments = range(total)

        # Group by score sequence
        score_groups = defaultdict(list)
        for bits in tournaments:
            adj = adj_from_bits(bits, n)
            sc = score_seq(adj, n)
            H = held_karp_H(adj, n)
            score_groups[sc].append((bits, H))

        print(f"\n--- n = {n} ---")
        print(f"Score classes: {len(score_groups)}")

        # For each SC score class, check if SC tournament maximizes H
        sc_scores = [s for s in score_groups if is_score_self_complementary(s)]
        print(f"Self-complementary score classes: {len(sc_scores)}")

        violations = 0
        for sc in sc_scores:
            entries = score_groups[sc]
            max_H = max(h for _, h in entries)

            # Check if any SC tournament achieves max_H
            sc_achieves_max = False
            nsc_max = 0
            sc_max = 0

            for bits, h in entries:
                adj = adj_from_bits(bits, n)
                if n <= 6:  # Only check SC for small n
                    is_sc = is_self_converse(adj, n)
                    if is_sc:
                        sc_max = max(sc_max, h)
                        if h == max_H:
                            sc_achieves_max = True
                    else:
                        nsc_max = max(nsc_max, h)

            if n <= 6:
                if not sc_achieves_max and sc_max > 0:
                    print(f"  VIOLATION: score {sc}: SC max={sc_max}, NSC max={nsc_max}, overall max={max_H}")
                    violations += 1
                elif sc_max > 0:
                    if nsc_max > 0:
                        gap = sc_max - nsc_max
                        print(f"  score {sc}: SC max={sc_max}, NSC max={nsc_max}, gap={gap}")
                    else:
                        print(f"  score {sc}: ALL are SC, max H={sc_max}")

        if n <= 6:
            print(f"  Violations: {violations}")


# ================================================================
# PART 2: PALEY vs INTERVAL COMPARISON
# ================================================================

def part2_paley_vs_interval():
    print("\n" + "=" * 70)
    print("PART 2: PALEY vs INTERVAL H COMPARISON")
    print("=" * 70)

    for p in [3, 5, 7, 11]:
        paley_adj = paley_tournament(p) if p % 4 == 3 else None
        interval_adj = interval_tournament(p)

        H_interval = held_karp_H(interval_adj, p)

        if paley_adj:
            H_paley = held_karp_H(paley_adj, p)
            qr = set()
            for x in range(1, p):
                qr.add((x * x) % p)
            print(f"\n  p={p} (= {p%4} mod 4)")
            print(f"    QR = {sorted(qr)}")
            print(f"    Interval = {list(range(1, (p+1)//2))}")
            print(f"    H(Paley) = {H_paley}")
            print(f"    H(Interval) = {H_interval}")
            if H_paley != H_interval:
                print(f"    Winner: {'PALEY' if H_paley > H_interval else 'INTERVAL'}")
                print(f"    Ratio: {H_paley / H_interval:.6f}")
            else:
                print(f"    EQUAL (same tournament or isomorphic)")
        else:
            print(f"\n  p={p} (≡ {p%4} mod 4, no Paley)")
            print(f"    Interval = {list(range(1, (p+1)//2))}")
            print(f"    H(Interval) = {H_interval}")

        # Also count cycles and compute independence polynomial for small p
        if p <= 7:
            adj = interval_adj
            oc = count_odd_cycles(adj, p)
            ip = independence_polynomial(adj, p)
            print(f"    Interval odd cycles: {oc}")
            print(f"    Interval I(Ω,x) = {ip}, I(Ω,2) = {eval_indpoly(ip, 2)}")

            if paley_adj and p > 3:
                oc2 = count_odd_cycles(paley_adj, p)
                ip2 = independence_polynomial(paley_adj, p)
                print(f"    Paley odd cycles: {oc2}")
                print(f"    Paley I(Ω,x) = {ip2}, I(Ω,2) = {eval_indpoly(ip2, 2)}")


# ================================================================
# PART 3: H/|Aut| FACTORIZATION FOR PALEY
# ================================================================

def part3_factorization():
    print("\n" + "=" * 70)
    print("PART 3: H(T_p)/|Aut(T_p)| FACTORIZATION")
    print("=" * 70)

    def prime_factors(n):
        """Return list of prime factors with multiplicity."""
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

    # Known values
    data = [
        (3, 3, 3),
        (7, 189, 21),
        (11, 95095, 55),
        # p=19: H=1172695746915, |Aut|=171
    ]

    for p, H, aut_size in data:
        ratio = H // aut_size
        factors = prime_factors(ratio)
        print(f"\n  p={p}: H={H}, |Aut|={aut_size}")
        print(f"    H/|Aut| = {ratio}")
        print(f"    Factorization: {' × '.join(str(f) for f in factors)}")
        print(f"    Unique primes: {sorted(set(factors))}")

    # Try to factor H for p=19
    H_19 = 1172695746915
    aut_19 = 171
    ratio_19 = H_19 // aut_19
    factors_19 = prime_factors(ratio_19)
    print(f"\n  p=19: H={H_19}, |Aut|={aut_19}")
    print(f"    H/|Aut| = {ratio_19}")
    print(f"    Factorization: {' × '.join(str(f) for f in factors_19)}")
    print(f"    Unique primes: {sorted(set(factors_19))}")

    # Also factor H itself
    for p, H, aut_size in data:
        factors = prime_factors(H)
        print(f"\n  H(T_{p}) = {H} = {' × '.join(str(f) for f in factors)}")

    factors = prime_factors(H_19)
    print(f"\n  H(T_19) = {H_19} = {' × '.join(str(f) for f in factors)}")

    # Look for patterns in H/|Aut| sequence
    print("\n  Pattern analysis of H/|Aut| sequence:")
    ratios = [1, 9, 1729, ratio_19]
    primes_list = [3, 7, 11, 19]
    for i, (p, r) in enumerate(zip(primes_list, ratios)):
        print(f"    p={p}: H/|Aut| = {r}")
        print(f"      mod p: {r % p}")
        print(f"      mod (p-1): {r % (p-1)}")
        print(f"      mod (p²): {r % (p*p)}")
        if i > 0:
            prev_r = ratios[i-1]
            print(f"      ratio to previous: {r / prev_r:.4f}")


# ================================================================
# PART 4: FORBIDDEN H VALUES AT n=8 (sampling)
# ================================================================

def part4_forbidden_values():
    print("\n" + "=" * 70)
    print("PART 4: FORBIDDEN H VALUES — n=8 SAMPLING")
    print("=" * 70)

    import random
    random.seed(2026)

    n = 8
    m = n * (n - 1) // 2  # 28
    sample_size = 100000

    H_values = Counter()
    for trial in range(sample_size):
        bits = random.randint(0, (1 << m) - 1)
        adj = adj_from_bits(bits, n)
        H = held_karp_H(adj, n)
        H_values[H] += 1

    all_H = sorted(H_values.keys())
    max_H = max(all_H)
    min_H = min(all_H)

    print(f"  Sampled {sample_size} random tournaments on {n} vertices")
    print(f"  H range: {min_H} to {max_H}")
    print(f"  Distinct H values: {len(all_H)}")

    # Find gaps
    expected = set(range(1, max_H + 1, 2))  # All odd numbers
    achieved = set(all_H)
    gaps = sorted(expected - achieved)

    print(f"  Gaps in odd values [1, {max_H}]: {len(gaps)}")
    if len(gaps) <= 50:
        print(f"  Gap values: {gaps}")
    else:
        print(f"  First 30 gaps: {gaps[:30]}")
        print(f"  Last 10 gaps: {gaps[-10:]}")

    # Check specifically for 7 and 21
    print(f"\n  H=7 count: {H_values.get(7, 0)}")
    print(f"  H=21 count: {H_values.get(21, 0)}")

    # Distribution near the gaps
    print("\n  H values near known permanent gaps:")
    for target in [7, 21]:
        nearby = [(h, H_values[h]) for h in sorted(H_values.keys()) if abs(h - target) <= 4]
        print(f"    Near H={target}: {nearby}")


# ================================================================
# PART 5: SC MAXIMIZER — ORBIT PAIRING DEEP ANALYSIS
# ================================================================

def part5_orbit_pairing():
    print("\n" + "=" * 70)
    print("PART 5: SC MAXIMIZER — ORBIT PAIRING MECHANISM")
    print("=" * 70)

    n = 6
    m = n * (n - 1) // 2

    # For n=6, analyze all SC tournaments and their anti-automorphism structure
    print(f"\n--- Detailed orbit analysis at n={n} ---")

    sc_data = []
    for bits in range(1 << m):
        adj = adj_from_bits(bits, n)
        if not is_self_converse(adj, n):
            continue

        H = held_karp_H(adj, n)
        sc = score_seq(adj, n)
        c3 = count_3cycles(adj, n)

        # Find involutory anti-automorphisms
        anti_auts = find_anti_auts(adj, n)
        invol_anti = [s for s in anti_auts if all(s[s[i]] == i for i in range(n))]

        # For each involutory anti-aut, find the orbit structure
        for sigma in invol_anti[:1]:  # Just analyze the first one
            orbits = []
            used = set()
            for i in range(n):
                if i not in used:
                    if sigma[i] == i:
                        orbits.append((i,))  # Fixed point
                    else:
                        orbits.append((i, sigma[i]))
                    used.add(i)
                    used.add(sigma[i])

            # Count how many 3-cycles come from orbit selections
            orbit_pairs = [o for o in orbits if len(o) == 2]
            fixed_pts = [o[0] for o in orbits if len(o) == 1]

            # A 3-cycle from orbit selections picks one vertex from each of 3 orbits
            orbit_3cycles = 0
            disjoint_pairs = 0
            if len(orbit_pairs) >= 3:
                for combo in combinations(range(len(orbit_pairs)), 3):
                    ops = [orbit_pairs[c] for c in combo]
                    # 2^3 = 8 vertex selections
                    for sel in range(8):
                        verts = [ops[i][1 if (sel >> i) & 1 else 0] for i in range(3)]
                        # Check if this is a 3-cycle
                        if adj[verts[0]][verts[1]] and adj[verts[1]][verts[2]] and adj[verts[2]][verts[0]]:
                            orbit_3cycles += 1
                        elif adj[verts[0]][verts[2]] and adj[verts[2]][verts[1]] and adj[verts[1]][verts[0]]:
                            orbit_3cycles += 1

            sc_data.append({
                'bits': bits, 'H': H, 'score': sc, 'c3': c3,
                'orbits': orbits, 'orbit_pairs': len(orbit_pairs),
                'fixed_pts': len(fixed_pts),
                'orbit_3cycles': orbit_3cycles,
                'num_anti_auts': len(anti_auts),
                'num_invol': len(invol_anti),
            })

    # Report
    print(f"  Total SC tournaments at n={n}: {len(sc_data)}")

    # Group by score
    by_score = defaultdict(list)
    for d in sc_data:
        by_score[d['score']].append(d)

    for sc in sorted(by_score.keys()):
        entries = by_score[sc]
        H_vals = sorted(set(d['H'] for d in entries))
        max_H = max(H_vals)
        print(f"\n  Score {sc}:")
        print(f"    SC tournaments: {len(entries)}")
        print(f"    H values: {H_vals}")

        # For max H entry
        max_entry = [d for d in entries if d['H'] == max_H][0]
        print(f"    Max H={max_H}: c3={max_entry['c3']}, orbit_pairs={max_entry['orbit_pairs']}, "
              f"fixed_pts={max_entry['fixed_pts']}, orbit_3cycles={max_entry['orbit_3cycles']}")

        # Compare with other entries
        if len(H_vals) > 1:
            for h in H_vals:
                if h != max_H:
                    entry = [d for d in entries if d['H'] == h][0]
                    print(f"    H={h}: c3={entry['c3']}, orbit_pairs={entry['orbit_pairs']}, "
                          f"orbit_3cycles={entry['orbit_3cycles']}")

    # Now for n=6, compare SC vs NSC within each SC score class
    print(f"\n--- SC vs NSC comparison at n={n} ---")
    all_by_score = defaultdict(list)
    for bits in range(1 << m):
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        if not is_score_self_complementary(sc):
            continue
        H = held_karp_H(adj, n)
        is_sc = is_self_converse(adj, n)
        c3 = count_3cycles(adj, n)
        ip = independence_polynomial(adj, n)
        alpha2 = ip[2] if len(ip) > 2 else 0
        all_by_score[sc].append({
            'bits': bits, 'H': H, 'c3': c3, 'is_sc': is_sc,
            'alpha1': ip[1] if len(ip) > 1 else 0,
            'alpha2': alpha2,
            'IP': ip,
        })

    for sc in sorted(all_by_score.keys()):
        entries = all_by_score[sc]
        sc_entries = [e for e in entries if e['is_sc']]
        nsc_entries = [e for e in entries if not e['is_sc']]

        if sc_entries and nsc_entries:
            sc_max_H = max(e['H'] for e in sc_entries)
            nsc_max_H = max(e['H'] for e in nsc_entries)
            sc_max_a2 = max(e['alpha2'] for e in sc_entries)
            nsc_max_a2 = max(e['alpha2'] for e in nsc_entries)
            sc_max_a1 = max(e['alpha1'] for e in sc_entries)
            nsc_max_a1 = max(e['alpha1'] for e in nsc_entries)

            print(f"\n  Score {sc}:")
            print(f"    SC: {len(sc_entries)} tours, max H={sc_max_H}, max a1={sc_max_a1}, max a2={sc_max_a2}")
            print(f"    NSC: {len(nsc_entries)} tours, max H={nsc_max_H}, max a1={nsc_max_a1}, max a2={nsc_max_a2}")
            print(f"    Gap: SC - NSC = {sc_max_H - nsc_max_H}")

            # Show IP distributions for SC and NSC
            sc_ips = Counter(tuple(e['IP']) for e in sc_entries)
            nsc_ips = Counter(tuple(e['IP']) for e in nsc_entries)
            print(f"    SC IP dist: {dict(sc_ips)}")
            print(f"    NSC IP dist: {dict(nsc_ips)}")


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":
    part1_sc_maximizer()
    part2_paley_vs_interval()
    part3_factorization()
    part4_forbidden_values()
    part5_orbit_pairing()
