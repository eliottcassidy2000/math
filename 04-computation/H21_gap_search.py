#!/usr/bin/env python3
"""
H21_gap_search.py -- Determine whether H=21 is achievable for any tournament.

H(T) = number of directed Hamiltonian paths in tournament T.
By OCF: H(T) = I(Omega(T), 2) = 1 + 2*i_1 + 4*i_2 + 8*i_3 + ...

Strategy:
1. Exhaustive search at n <= 7 (2^21 = 2M at n=7)
2. Large random sampling at n=8, 9, 10
3. Systematic edge-flip from transitive at n=8
4. OCF decomposition analysis
"""

import random
import sys
import time
from itertools import combinations, permutations
from collections import Counter, defaultdict

# ============================================================
# Core computation: fast H via array-based DP
# ============================================================

def count_H_from_bits(n, bits):
    """Count H for tournament encoded as integer (upper triangle bits)."""
    # Build adjacency as flat array for speed
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return count_H(adj, n)


def count_H(A, n):
    """Count Hamiltonian paths via bitmask DP (array-based, fast)."""
    if n <= 1:
        return 1
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not (mask & (1 << v)) or c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += c
    return sum(dp[full])


def random_tournament(n):
    """Generate a random tournament on n vertices."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def tournament_from_bits(n, bits):
    """Construct tournament from integer encoding of upper triangle."""
    A = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                A[i][j] = 1
            else:
                A[j][i] = 1
            k += 1
    return A


def transitive_tournament(n):
    """Transitive tournament: i beats j iff i < j."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A


def flip_edge(A, n, i, j):
    """Flip the arc between i and j."""
    B = [row[:] for row in A]
    B[i][j], B[j][i] = B[j][i], B[i][j]
    return B


# ============================================================
# OCF computation
# ============================================================

def find_odd_cycle_vsets(A, n):
    """Find all directed odd cycle VERTEX SETS in tournament A.
    A vertex set supports a directed odd cycle iff the sub-tournament
    on those vertices has a Hamiltonian cycle."""
    cycle_sets = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            # Check for Hamiltonian cycle on this subset
            s = size
            sub = list(subset)
            idx = {v: i for i, v in enumerate(sub)}
            # Build sub-adjacency
            sub_A = [[0]*s for _ in range(s)]
            for a in sub:
                for b in sub:
                    if a != b and A[a][b]:
                        sub_A[idx[a]][idx[b]] = 1
            # DP: Hamiltonian paths from node 0
            dp = [[0]*s for _ in range(1 << s)]
            dp[1][0] = 1
            for mask in range(1, 1 << s):
                if not (mask & 1):
                    continue
                for v in range(s):
                    c = dp[mask][v]
                    if not (mask & (1 << v)) or c == 0:
                        continue
                    for u in range(s):
                        if mask & (1 << u):
                            continue
                        if sub_A[v][u]:
                            dp[mask | (1 << u)][u] += c
            full = (1 << s) - 1
            has_cycle = any(dp[full][v] > 0 and sub_A[v][0] for v in range(1, s))
            if has_cycle:
                cycle_sets.append(frozenset(subset))
    return cycle_sets


def compute_ocf_decomposition(A, n):
    """Compute full OCF: returns (i_0, i_1, i_2, ...) and I(Omega,2)."""
    cycle_sets = find_odd_cycle_vsets(A, n)
    m = len(cycle_sets)

    if m == 0:
        return [1], 1  # just i_0=1

    # Build adjacency bitmasks for conflict graph
    adj_bits = [0] * m
    for a in range(m):
        for b in range(a+1, m):
            if cycle_sets[a] & cycle_sets[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    # Count independent sets by size
    size_counts = Counter()
    if m <= 25:
        for mask in range(1 << m):
            ok = True
            seen = 0
            temp = mask
            while temp:
                v = (temp & -temp).bit_length() - 1
                if adj_bits[v] & seen:
                    ok = False
                    break
                seen |= 1 << v
                temp &= temp - 1
            if ok:
                size_counts[bin(mask).count('1')] += 1
    else:
        return None, None

    max_k = max(size_counts.keys()) if size_counts else 0
    i_vals = [size_counts.get(k, 0) for k in range(max_k + 1)]
    I_at_2 = sum(i_vals[k] * (2**k) for k in range(len(i_vals)))
    return i_vals, I_at_2


# ============================================================
# PHASE 1: Exhaustive search n <= 7
# ============================================================

def exhaustive_search(n):
    """Exhaustively check ALL tournaments on n vertices."""
    m = n * (n - 1) // 2
    total = 1 << m
    h_counts = Counter()
    found_21 = None

    print(f"\n--- Exhaustive n={n}: {total} tournaments ---")
    t0 = time.time()

    for bits in range(total):
        A = tournament_from_bits(n, bits)
        h = count_H(A, n)
        h_counts[h] += 1
        if h == 21 and found_21 is None:
            found_21 = bits
            print(f"  *** H=21 FOUND at bits={bits} ***")

    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.1f}s")

    all_odd = sorted([h for h in h_counts if h <= 60 and h % 2 == 1])
    print(f"  Odd H values <= 60: {all_odd}")

    missing = [h for h in range(1, 52, 2) if h not in h_counts]
    print(f"  Missing odd H <= 51: {missing}")

    # Show counts near 21
    print(f"  Counts near 21: ", end="")
    for hv in range(15, 28, 2):
        print(f"H={hv}:{h_counts.get(hv, 0)} ", end="")
    print()

    if found_21 is not None:
        print(f"  H=21 ACHIEVED! bits={found_21}")
    else:
        print(f"  H=21 NOT achieved at n={n}")

    return h_counts, found_21


# ============================================================
# PHASE 2: Random sampling n=8,9,10
# ============================================================

def random_search(n, num_samples):
    """Large random sample search."""
    h_counts = Counter()
    found_21 = None

    print(f"\n--- Random search n={n}: {num_samples} samples ---")
    t0 = time.time()

    for trial in range(num_samples):
        A = random_tournament(n)
        h = count_H(A, n)
        h_counts[h] += 1
        if h == 21 and found_21 is None:
            found_21 = [row[:] for row in A]
            print(f"  *** H=21 FOUND at trial {trial} ***")

        if (trial + 1) % 100000 == 0:
            elapsed = time.time() - t0
            rate = (trial+1) / elapsed
            print(f"  {trial+1}/{num_samples} ({elapsed:.0f}s, {rate:.0f}/s)")
            sys.stdout.flush()

    elapsed = time.time() - t0
    print(f"  Total time: {elapsed:.1f}s")

    # Show near-21 counts
    print(f"  Counts near 21: ", end="")
    for hv in range(15, 28, 2):
        print(f"H={hv}:{h_counts.get(hv, 0)} ", end="")
    print()

    missing = [h for h in range(1, 52, 2) if h not in h_counts and h <= max(h_counts)]
    print(f"  Missing odd H <= 51: {missing}")

    if found_21:
        print(f"  H=21 ACHIEVED!")
    else:
        print(f"  H=21 NOT found in {num_samples} samples")

    return h_counts, found_21


# ============================================================
# PHASE 3: Systematic edge-flip from transitive at n=8
# ============================================================

def systematic_flip_search(n, max_flips=4):
    """Systematically flip 1..max_flips edges from transitive tournament."""
    A_trans = transitive_tournament(n)
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    ne = len(edges)
    achieved = set()
    achieved.add(1)  # transitive has H=1

    print(f"\n--- Systematic flip search n={n}, up to {max_flips} flips ---")
    t0 = time.time()

    # 1-flip
    for e1 in range(ne):
        A = flip_edge(A_trans, n, edges[e1][0], edges[e1][1])
        achieved.add(count_H(A, n))
    print(f"  After 1-flip ({ne} combos): H<=30 = {sorted(h for h in achieved if h <= 30)}")

    # 2-flip
    if max_flips >= 2:
        cnt = 0
        for e1 in range(ne):
            A1 = flip_edge(A_trans, n, edges[e1][0], edges[e1][1])
            for e2 in range(e1+1, ne):
                A2 = flip_edge(A1, n, edges[e2][0], edges[e2][1])
                achieved.add(count_H(A2, n))
                cnt += 1
        print(f"  After 2-flip ({cnt} combos): H<=30 = {sorted(h for h in achieved if h <= 30)}")

    # 3-flip
    if max_flips >= 3:
        cnt = 0
        for e1 in range(ne):
            A1 = flip_edge(A_trans, n, edges[e1][0], edges[e1][1])
            for e2 in range(e1+1, ne):
                A2 = flip_edge(A1, n, edges[e2][0], edges[e2][1])
                for e3 in range(e2+1, ne):
                    A3 = flip_edge(A2, n, edges[e3][0], edges[e3][1])
                    achieved.add(count_H(A3, n))
                    cnt += 1
        print(f"  After 3-flip ({cnt} combos): H<=30 = {sorted(h for h in achieved if h <= 30)}")

    # 4-flip
    if max_flips >= 4:
        cnt = 0
        for e1 in range(ne):
            A1 = flip_edge(A_trans, n, edges[e1][0], edges[e1][1])
            for e2 in range(e1+1, ne):
                A2 = flip_edge(A1, n, edges[e2][0], edges[e2][1])
                for e3 in range(e2+1, ne):
                    A3 = flip_edge(A2, n, edges[e3][0], edges[e3][1])
                    for e4 in range(e3+1, ne):
                        A4 = flip_edge(A3, n, edges[e4][0], edges[e4][1])
                        achieved.add(count_H(A4, n))
                        cnt += 1
            if (e1+1) % 5 == 0:
                print(f"    4-flip progress: {e1+1}/{ne} outer edges...")
                sys.stdout.flush()
        print(f"  After 4-flip ({cnt} combos): H<=30 = {sorted(h for h in achieved if h <= 30)}")

    elapsed = time.time() - t0
    print(f"  Total time: {elapsed:.1f}s")

    missing = [h for h in range(1, 52, 2) if h not in achieved]
    print(f"  Missing odd H <= 51: {missing}")

    if 21 in achieved:
        print(f"  H=21 ACHIEVED via edge flips!")
    else:
        print(f"  H=21 NOT achieved")

    return achieved


# ============================================================
# PHASE 4: OCF analysis of near-21 tournaments
# ============================================================

def ocf_analysis(n, num_samples=50000):
    """Collect tournaments with H near 21 and analyze OCF structure."""
    print(f"\n--- OCF analysis n={n}: {num_samples} samples ---")

    near_h = defaultdict(list)
    t0 = time.time()

    for trial in range(num_samples):
        A = random_tournament(n)
        h = count_H(A, n)
        if 15 <= h <= 27 and len(near_h[h]) < 3:
            near_h[h].append([row[:] for row in A])

    elapsed = time.time() - t0
    print(f"  Sampling done in {elapsed:.1f}s")

    print(f"  H values with examples: {sorted(near_h.keys())}")

    for hv in sorted(near_h.keys()):
        if near_h[hv]:
            A = near_h[hv][0]
            i_vals, I_at_2 = compute_ocf_decomposition(A, n)
            if i_vals is not None:
                i_str = ', '.join(f'i_{k}={i_vals[k]}' for k in range(len(i_vals)))
                print(f"    H={hv}: {i_str}, verify I(Omega,2)={I_at_2}")
            else:
                print(f"    H={hv}: too many cycles for full OCF")

    if 21 in near_h:
        print(f"\n  *** H=21 examples found! ***")
        for A in near_h[21]:
            i_vals, I_at_2 = compute_ocf_decomposition(A, n)
            print(f"    I(Omega,2)={I_at_2}, decomposition={i_vals}")
            print(f"    Adjacency: {A}")
    else:
        print(f"\n  H=21 NOT found at n={n}")


# ============================================================
# PHASE 5: Theoretical analysis
# ============================================================

def theoretical_analysis():
    """Print all (i_1, i_2, i_3, ...) that give H=21."""
    print("\n--- Theoretical: all (i_1, i_2, ...) giving H=21 ---")
    print("H = 1 + 2*i_1 + 4*i_2 + 8*i_3 + 16*i_4")
    print("Need: 2*i_1 + 4*i_2 + 8*i_3 + 16*i_4 = 20\n")

    solutions = []
    for i1 in range(11):
        r1 = 20 - 2*i1
        if r1 < 0: break
        for i2 in range(r1//4 + 1):
            r2 = r1 - 4*i2
            if r2 < 0: break
            for i3 in range(r2//8 + 1):
                r3 = r2 - 8*i3
                if r3 < 0: break
                if r3 == 0:
                    solutions.append((i1, i2, i3))
                elif r3 % 16 == 0:
                    i4 = r3 // 16
                    solutions.append((i1, i2, i3, i4))

    for s in solutions:
        parts = ', '.join(f'i_{k+1}={s[k]}' for k in range(len(s)) if s[k] > 0)
        if not parts:
            parts = "(no cycles -- impossible since H=21>1)"
        print(f"  {parts}")
    print(f"\n  Total: {len(solutions)} decompositions")
    print(f"\n  Note: i_1 = alpha_1 = total directed odd cycles (vertex sets)")
    print(f"  i_2 = number of vertex-disjoint PAIRS of odd cycles")
    print(f"  Key parity: i_1 must be EVEN (since H is odd and H=1+2*i_1+...)")


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("H=21 GAP SEARCH -- Is H(T)=21 achievable for any tournament T?")
    print("=" * 70)

    theoretical_analysis()

    # Phase 1: Exhaustive for n=3..7
    print("\n" + "=" * 70)
    print("PHASE 1: Exhaustive search (n=3 through 7)")
    print("=" * 70)
    all_achieved = set()
    for n in range(3, 8):
        hc, found = exhaustive_search(n)
        all_achieved.update(hc.keys())

    missing_global = [h for h in range(1, 52, 2) if h not in all_achieved]
    print(f"\nGlobal missing odd H <= 51 after exhaustive n<=7: {missing_global}")

    # Phase 2: Random sampling at n=8,9,10
    print("\n" + "=" * 70)
    print("PHASE 2: Random sampling (n=8, 9, 10)")
    print("=" * 70)

    for n, samples in [(8, 50000), (9, 20000), (10, 5000), (11, 2000), (12, 500)]:
        hc, found = random_search(n, samples)
        all_achieved.update(hc.keys())

    missing_global = [h for h in range(1, 52, 2) if h not in all_achieved]
    print(f"\nGlobal missing odd H <= 51 after random n<=12: {missing_global}")

    # Phase 3: Systematic edge flips at n=8
    print("\n" + "=" * 70)
    print("PHASE 3: Systematic edge-flip search (n=8)")
    print("=" * 70)
    achieved_flip = systematic_flip_search(8, max_flips=3)
    all_achieved.update(achieved_flip)

    missing_global = [h for h in range(1, 52, 2) if h not in all_achieved]
    print(f"\nGlobal missing odd H <= 51 after flip search: {missing_global}")

    # Phase 4: OCF analysis at n=7, 8
    print("\n" + "=" * 70)
    print("PHASE 4: OCF decomposition analysis")
    print("=" * 70)
    ocf_analysis(7, 30000)

    # Final summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"Missing odd H values <= 51: {[h for h in range(1,52,2) if h not in all_achieved]}")
    if 21 in all_achieved:
        print("RESULT: H=21 IS achievable.")
    else:
        print("RESULT: H=21 was NOT found in any search.")
        print("Combined evidence: exhaustive n<=7, 50K samples n=8,")
        print("20K n=9, 5K n=10, 2K n=11, 500 n=12, systematic 3-flip n=8.")
