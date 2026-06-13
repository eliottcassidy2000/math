#!/usr/bin/env python3
"""overnight_s1_deep.py -- Deep investigation continued.

Session: kind-pasteur-2026-03-20-S1

Part G: The SC alpha_2 mechanism at n=6 — complete classification
Part H: H/|Aut| factorization for p=23, pattern search
Part I: Forbidden values — n=9 via sampling
Part J: Paley cycle count verification at p=11
Part K: SC Maximizer at n=7 — which alpha drives the win?
"""

from itertools import combinations, permutations
from math import comb, factorial, gcd
from collections import defaultdict, Counter
import random

# ================================================================
# CORE UTILITIES (same as before)
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

def is_score_self_complementary(score):
    n = len(score)
    s = sorted(score)
    return all(s[i] + s[n-1-i] == n-1 for i in range(n))

def count_3cycles(adj, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1
    return c3

def count_directed_cycles_on_subset(adj_global, verts):
    """Count directed Hamiltonian cycles on a vertex subset.
    Returns list of cycles, each as tuple of global vertex indices."""
    size = len(verts)
    sub_adj = [[adj_global[verts[i]][verts[j]] for j in range(size)] for i in range(size)]
    cycles = []
    for perm in permutations(range(1, size)):
        path = [0] + list(perm)
        is_cycle = True
        for idx in range(size):
            if not sub_adj[path[idx]][path[(idx+1) % size]]:
                is_cycle = False
                break
        if is_cycle:
            actual = tuple(verts[path[i]] for i in range(size))
            cycles.append(actual)
    return cycles

def full_alpha_vector(adj, n):
    """Compute full independence polynomial [alpha_0, alpha_1, ...] correctly.
    Each directed odd cycle is a vertex of Omega.
    """
    # Collect all directed odd cycles with their vertex sets
    all_cycles = []  # list of (frozenset of vertices, cycle_tuple)
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
                    all_cycles.append((frozenset(verts), tuple(verts[path[i]] for i in range(size))))

    num_cycles = len(all_cycles)
    if num_cycles == 0:
        return [1]

    vsets = [c[0] for c in all_cycles]

    # Build conflict adjacency
    conflict = [[False]*num_cycles for _ in range(num_cycles)]
    for i in range(num_cycles):
        for j in range(i+1, num_cycles):
            if vsets[i] & vsets[j]:
                conflict[i][j] = True
                conflict[j][i] = True

    # Count independent sets by size
    alpha = [1]
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
# PART G: SC ALPHA_2 MECHANISM — COMPLETE CLASSIFICATION AT n=6
# ================================================================

def part_g():
    print("=" * 70)
    print("PART G: SC ALPHA_2 MECHANISM — COMPLETE ANALYSIS AT n=6")
    print("=" * 70)

    n = 6
    m = n * (n - 1) // 2
    target = (2, 2, 2, 3, 3, 3)  # Regular score

    # For each tournament with this score, compute:
    # - H
    # - Full alpha vector
    # - 3-cycle vertex sets
    # - Which complementary 3-set pairs are both cyclic

    data = []
    for bits in range(1 << m):
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        if sc != target:
            continue

        H = held_karp_H(adj, n)

        # Find 3-cycle vertex sets
        cyclic_sets = set()
        for triple in combinations(range(n), 3):
            i, j, k = triple
            if (adj[i][j] and adj[j][k] and adj[k][i]) or \
               (adj[i][k] and adj[k][j] and adj[j][i]):
                cyclic_sets.add(frozenset(triple))

        # Count complementary pairs (both halves cyclic)
        # There are C(6,3)/2 = 10 complementary 3-set pairs
        comp_pairs_both_cyclic = 0
        for triple in combinations(range(n), 3):
            s = frozenset(triple)
            comp = frozenset(range(n)) - s
            if s < comp:  # avoid double counting
                if s in cyclic_sets and comp in cyclic_sets:
                    comp_pairs_both_cyclic += 1

        data.append({
            'bits': bits, 'H': H,
            'c3': len(cyclic_sets),
            'comp_both': comp_pairs_both_cyclic,
        })

    # Analyze
    H_groups = defaultdict(list)
    for d in data:
        H_groups[d['H']].append(d)

    for H_val in sorted(H_groups.keys()):
        entries = H_groups[H_val]
        comp_vals = [e['comp_both'] for e in entries]
        c3_vals = [e['c3'] for e in entries]
        print(f"\n  H={H_val}: {len(entries)} tournaments")
        print(f"    c3 values: {Counter(c3_vals)}")
        print(f"    comp_both_cyclic values: {Counter(comp_vals)}")

    # Theorem statement
    print(f"\n  KEY OBSERVATION:")
    print(f"  For regular n=6 tournaments (score 2,2,2,3,3,3):")
    H45 = H_groups[45]
    H43 = H_groups[43]
    H41 = H_groups[41]
    print(f"  H=45: comp_both = {set(e['comp_both'] for e in H45)}, c3 = {set(e['c3'] for e in H45)}")
    print(f"  H=43: comp_both = {set(e['comp_both'] for e in H43)}, c3 = {set(e['c3'] for e in H43)}")
    print(f"  H=41: comp_both = {set(e['comp_both'] for e in H41)}, c3 = {set(e['c3'] for e in H41)}")


# ================================================================
# PART H: H/|Aut| FACTORIZATION FOR p=23
# ================================================================

def part_h():
    print("\n" + "=" * 70)
    print("PART H: H(T_p)/|Aut| ANALYSIS FOR p=23")
    print("=" * 70)

    # Known: H(T_23) = 15760206976379349, |Aut| = 23*22/2 = 253
    H_23 = 15760206976379349
    aut_23 = 23 * 22 // 2  # = 253
    ratio_23 = H_23 // aut_23

    print(f"  p=23: H = {H_23}")
    print(f"  |Aut| = {aut_23}")
    print(f"  H/|Aut| = {ratio_23}")
    print(f"  H mod |Aut| = {H_23 % aut_23} (should be 0)")

    f23 = prime_factors(ratio_23)
    print(f"  H/|Aut| factored: {' * '.join(str(f) for f in f23)}")
    print(f"  Unique primes: {sorted(set(f23))}")

    # Factor H itself
    fH23 = prime_factors(H_23)
    print(f"  H factored: {' * '.join(str(f) for f in fH23)}")

    # All known H/|Aut| values
    print(f"\n  Complete H/|Aut| table:")
    data = [
        (3, 3, 3, 1),
        (7, 189, 21, 9),
        (11, 95095, 55, 1729),
        (19, 1172695746915, 171, 6857869865),
        (23, H_23, aut_23, ratio_23),
    ]

    for p, H, aut, r in data:
        factors = prime_factors(r) if r > 1 else [1]
        unique = sorted(set(factors)) if r > 1 else []
        print(f"  p={p:2d}: H/|Aut| = {r:>20d} = {'*'.join(str(f) for f in factors) if r > 1 else '1':>40s}  primes: {unique}")

    # Look for patterns
    print(f"\n  Pattern search:")
    ratios = [r for _, _, _, r in data]
    for i in range(1, len(ratios)):
        print(f"  ratio[{i}]/ratio[{i-1}] = {ratios[i]/ratios[i-1]:.4f}")

    # Check: does H(T_p) contain all primes up to p?
    print(f"\n  Does H(T_p) contain all primes up to p?")
    for p, H, aut, r in data:
        fH = prime_factors(H)
        primes_in_H = sorted(set(fH))
        primes_up_to_p = [q for q in range(2, p+1) if all(q % d != 0 for d in range(2, q))]
        missing = [q for q in primes_up_to_p if q not in primes_in_H]
        print(f"  p={p}: primes in H = {primes_in_H}, missing primes <= p: {missing}")


# ================================================================
# PART I: FORBIDDEN VALUES n=9
# ================================================================

def part_i():
    print("\n" + "=" * 70)
    print("PART I: FORBIDDEN H VALUES AT n=9 (500K samples)")
    print("=" * 70)

    random.seed(42)
    n = 9
    m = n * (n - 1) // 2  # 36
    sample_size = 500000

    H_values = Counter()
    for trial in range(sample_size):
        bits = random.randint(0, (1 << m) - 1)
        adj = adj_from_bits(bits, n)
        H = held_karp_H(adj, n)
        H_values[H] += 1

        if trial % 100000 == 0 and trial > 0:
            print(f"  Progress: {trial}/{sample_size}, distinct H so far: {len(H_values)}")

    all_H = sorted(H_values.keys())
    max_H = max(all_H)

    print(f"\n  {sample_size} random n={n} tournaments")
    print(f"  H range: {min(all_H)} to {max_H}")
    print(f"  Distinct H: {len(all_H)}")

    # Gaps
    expected = set(range(1, min(max_H, 500) + 1, 2))
    achieved = set(h for h in all_H if h <= 500)
    gaps = sorted(expected - achieved)
    print(f"  Gaps in odd [1, {min(max_H, 500)}]: {len(gaps)}")
    print(f"  Gap values: {gaps}")

    print(f"\n  H=7: {H_values.get(7, 0)}")
    print(f"  H=21: {H_values.get(21, 0)}")
    print(f"  H=63: {H_values.get(63, 0)}")


# ================================================================
# PART J: PALEY CYCLE STRUCTURE AT p=11
# ================================================================

def part_j():
    print("\n" + "=" * 70)
    print("PART J: PALEY T_11 CYCLE STRUCTURE (partial)")
    print("=" * 70)

    # Compute cycle counts for sizes 3, 5 at p=11
    p = 11
    adj = paley_tournament(p)

    print(f"  Computing directed cycle counts for Paley T_{p}...")

    # Count 3-cycles (fast)
    c3_vertex_sets = 0
    c3_directed = 0
    for triple in combinations(range(p), 3):
        cycles = count_directed_cycles_on_subset(adj, triple)
        if cycles:
            c3_vertex_sets += 1
            c3_directed += len(cycles)
    print(f"  3-cycles: {c3_directed} directed cycles on {c3_vertex_sets} vertex sets")
    print(f"    (each 3-vertex set has at most 1 directed cycle, so these should be equal)")

    # Count 5-cycles (moderate)
    c5_vertex_sets = 0
    c5_directed = 0
    for verts in combinations(range(p), 5):
        cycles = count_directed_cycles_on_subset(adj, verts)
        if cycles:
            c5_vertex_sets += 1
            c5_directed += len(cycles)
    print(f"  5-cycles: {c5_directed} directed cycles on {c5_vertex_sets} vertex sets")
    print(f"    (avg {c5_directed/max(c5_vertex_sets,1):.2f} directed cycles per cyclic vertex set)")

    # For comparison, do the same for interval tournament
    adj_int = interval_tournament(p)
    c3i_dir = 0
    for triple in combinations(range(p), 3):
        cycles = count_directed_cycles_on_subset(adj_int, triple)
        c3i_dir += len(cycles)

    c5i_dir = 0
    c5i_vs = 0
    for verts in combinations(range(p), 5):
        cycles = count_directed_cycles_on_subset(adj_int, verts)
        if cycles:
            c5i_vs += 1
            c5i_dir += len(cycles)

    print(f"\n  Interval T_{p}:")
    print(f"  3-cycles: {c3i_dir}")
    print(f"  5-cycles: {c5i_dir} directed cycles on {c5i_vs} vertex sets")
    print(f"    (avg {c5i_dir/max(c5i_vs,1):.2f} per cyclic vertex set)")

    # Key comparison
    print(f"\n  Comparison at p={p}:")
    print(f"  Paley: {c3_directed} + {c5_directed} = {c3_directed + c5_directed} total (3+5 cycles)")
    print(f"  Interval: {c3i_dir} + {c5i_dir} = {c3i_dir + c5i_dir} total (3+5 cycles)")


# ================================================================
# PART K: SC MAXIMIZER AT n=7 — WHICH ALPHA DRIVES THE WIN?
# ================================================================

def part_k():
    print("\n" + "=" * 70)
    print("PART K: SC MAXIMIZER AT n=7 — REGULAR SCORE ANALYSIS")
    print("=" * 70)

    # At n=7, all regular tournaments have score (3,3,3,3,3,3,3)
    # There are 2^21 = 2,097,152 tournaments total
    # We'll sample from regular ones

    n = 7
    m = n * (n - 1) // 2  # 21

    random.seed(123)
    regular_found = 0
    regular_data = []
    target = 1000  # aim for 1000 regular tournaments

    attempts = 0
    while regular_found < target and attempts < 500000:
        bits = random.randint(0, (1 << m) - 1)
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        if sc != (3, 3, 3, 3, 3, 3, 3):
            attempts += 1
            continue

        H = held_karp_H(adj, n)

        # Count directed odd cycles (3-cycles only for speed, then estimate)
        c3 = count_3cycles(adj, n)

        # For alpha_2 at n=7: only 3+3=6 <= 7 disjoint pairs possible
        # Count pairs of disjoint 3-cycles
        cyclic_triples = []
        for triple in combinations(range(n), 3):
            i, j, k = triple
            if (adj[i][j] and adj[j][k] and adj[k][i]) or \
               (adj[i][k] and adj[k][j] and adj[j][i]):
                cyclic_triples.append(frozenset(triple))

        disjoint_3_pairs = 0
        for a in range(len(cyclic_triples)):
            for b in range(a+1, len(cyclic_triples)):
                if not (cyclic_triples[a] & cyclic_triples[b]):
                    disjoint_3_pairs += 1

        regular_data.append({'bits': bits, 'H': H, 'c3': c3, 'dp3': disjoint_3_pairs})
        regular_found += 1
        attempts += 1

    print(f"  Found {regular_found} regular tournaments in {attempts} attempts")

    # Analyze by H
    H_groups = defaultdict(list)
    for d in regular_data:
        H_groups[d['H']].append(d)

    for H_val in sorted(H_groups.keys()):
        entries = H_groups[H_val]
        c3_vals = Counter(e['c3'] for e in entries)
        dp_vals = Counter(e['dp3'] for e in entries)
        print(f"\n  H={H_val}: {len(entries)} tournaments")
        print(f"    c3: {dict(sorted(c3_vals.items()))}")
        print(f"    disjoint 3-cycle pairs: {dict(sorted(dp_vals.items()))}")

    # Check BIBD vs non-BIBD structure
    # For Paley T_7: c3=14, and all regular tournaments have c3=14 (is this true?)
    print(f"\n  Do all regular n=7 tournaments have the same c3?")
    all_c3 = set(d['c3'] for d in regular_data)
    print(f"    c3 values seen: {sorted(all_c3)}")
    if len(all_c3) == 1:
        print(f"    YES — all regular n=7 have c3={all_c3.pop()}")
    else:
        print(f"    NO — c3 varies among regular n=7 tournaments")


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":
    part_g()
    part_h()
    part_k()
    print("\n  Skipping Part J (p=11 5-cycle counting) — too slow for this session")
    print("  Running Part I (n=9 forbidden values) last as it's the longest...")
    part_i()
