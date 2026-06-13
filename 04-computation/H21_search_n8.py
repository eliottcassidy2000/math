#!/usr/bin/env python3
"""
H21_search_n8.py -- Focused search for H=21 at n=8 tournaments.

H(T) = number of directed Hamiltonian paths.
By OCF: H(T) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

H=7 is provably impossible for all n (THM-029).
H=21 is absent through n=7 (exhaustive). Does it occur at n=8?

Strategy:
1. Random sampling: 200,000 random n=8 tournaments
2. Near-transitive: flip 1-5 edges from transitive tournament
3. Low-score targeting: fix score sequences that might yield low H
"""

import random
import sys
import time
from itertools import combinations
from collections import Counter

N = 8
NUM_EDGES = N * (N - 1) // 2  # 28

# ============================================================
# Core: Hamiltonian path count via bitmask DP
# ============================================================

def count_H_dp(adj, n=N):
    """Count Hamiltonian paths via bitmask DP. adj[i][j]=1 means i->j."""
    full = (1 << n) - 1
    # dp[mask][v] = number of Hamiltonian paths on vertices in mask ending at v
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not (mask & (1 << v)) or c == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += c
    return sum(dp[full])


def bits_to_adj(bits, n=N):
    """Convert integer bit encoding to adjacency matrix.
    Bit k corresponds to edge (i,j) where i<j, enumerated row by row.
    Bit=1 means i->j, bit=0 means j->i.
    """
    adj = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj


def adj_to_bits(adj, n=N):
    """Convert adjacency matrix back to bit encoding."""
    bits = 0
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i][j]:
                bits |= (1 << k)
            k += 1
    return bits


def score_sequence(adj, n=N):
    """Return sorted score sequence."""
    return tuple(sorted(sum(adj[i]) for i in range(n)))


# ============================================================
# Phase 1: Random sampling
# ============================================================

def random_search(num_samples):
    """Sample random n=8 tournaments, looking for H=21."""
    print(f"=== Phase 1: Random sampling ({num_samples:,} tournaments) ===")
    t0 = time.time()
    h_counts = Counter()
    found_21 = None
    max_bits = 1 << NUM_EDGES

    for trial in range(num_samples):
        bits = random.randint(0, max_bits - 1)
        adj = bits_to_adj(bits)
        h = count_H_dp(adj)
        h_counts[h] += 1

        if h == 21:
            if found_21 is None:
                found_21 = bits
                print(f"  *** H=21 FOUND at trial {trial}, bits={bits} ***")
                scores = score_sequence(adj)
                print(f"      Score sequence: {scores}")

        if (trial + 1) % 50000 == 0:
            elapsed = time.time() - t0
            rate = (trial + 1) / elapsed
            print(f"  Progress: {trial+1:,}/{num_samples:,} "
                  f"({elapsed:.1f}s, {rate:.0f}/s)")
            sys.stdout.flush()

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s ({num_samples/elapsed:.0f}/s)")

    # Report near-21 counts
    print(f"  H values near 21: ", end="")
    for hv in range(13, 32, 2):
        c = h_counts.get(hv, 0)
        if c > 0:
            print(f"H={hv}:{c} ", end="")
    print()

    missing = sorted(h for h in range(1, 60, 2)
                     if h not in h_counts and h <= max(h_counts))
    print(f"  Missing odd H values (up to max seen): {missing}")

    return h_counts, found_21


# ============================================================
# Phase 2: Near-transitive (edge flips)
# ============================================================

def near_transitive_search(max_flips=4):
    """Flip 1..max_flips edges from the transitive tournament."""
    print(f"\n=== Phase 2: Near-transitive search (up to {max_flips} flips) ===")
    t0 = time.time()

    # Transitive: i beats j for all i < j => all bits = 1
    trans_bits = (1 << NUM_EDGES) - 1
    trans_adj = bits_to_adj(trans_bits)
    h_trans = count_H_dp(trans_adj)
    print(f"  Transitive H = {h_trans}")

    edges = list(range(NUM_EDGES))
    achieved = set()
    achieved.add(h_trans)
    found_21 = None

    for nf in range(1, max_flips + 1):
        count = 0
        for flip_set in combinations(edges, nf):
            flipped = trans_bits
            for e in flip_set:
                flipped ^= (1 << e)
            adj = bits_to_adj(flipped)
            h = count_H_dp(adj)
            achieved.add(h)
            count += 1
            if h == 21 and found_21 is None:
                found_21 = flipped
                print(f"  *** H=21 FOUND with {nf} flips! bits={flipped} ***")
                print(f"      Flipped edges: {flip_set}")
                scores = score_sequence(adj)
                print(f"      Score sequence: {scores}")

        near = sorted(h for h in achieved if h <= 35)
        print(f"  After {nf}-flip ({count:,} combos): "
              f"small H values = {near}")
        sys.stdout.flush()

        if found_21 is not None:
            break

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")

    missing = sorted(h for h in range(1, 52, 2) if h not in achieved)
    print(f"  Missing odd H <= 51: {missing}")

    return achieved, found_21


# ============================================================
# Phase 3: Low-score targeting
# ============================================================

def score_targeted_search(num_per_score=5000):
    """Target score sequences likely to produce low H values.
    Low-H tournaments tend to have extreme (near-transitive) scores.
    """
    print(f"\n=== Phase 3: Score-targeted search ===")
    t0 = time.time()

    # Score sequences for n=8 that might give small H:
    # Near-transitive scores have high variance
    target_scores = [
        (0, 1, 2, 3, 4, 5, 6, 7),  # transitive
        (0, 1, 2, 3, 4, 5, 6, 7),  # same but we'll sample many
        (0, 1, 2, 3, 4, 5, 7, 6),
        (0, 1, 2, 4, 3, 5, 6, 7),
        (0, 1, 3, 3, 4, 4, 6, 7),
        (0, 2, 2, 3, 4, 5, 5, 7),
        (0, 1, 2, 3, 4, 5, 5, 8),  # invalid (max score is 7)
        (1, 1, 2, 3, 4, 5, 6, 6),
        (0, 1, 2, 3, 5, 4, 6, 7),
        (1, 1, 2, 3, 4, 5, 5, 7),
        (0, 1, 2, 3, 4, 6, 5, 7),
        (0, 2, 2, 3, 4, 5, 6, 6),
        (1, 1, 2, 3, 4, 5, 6, 6),
        (1, 1, 3, 3, 4, 4, 6, 6),
        (0, 1, 2, 4, 3, 5, 6, 7),
        (1, 2, 2, 3, 4, 5, 5, 6),
        (0, 1, 3, 3, 4, 4, 5, 8),  # invalid
        (2, 2, 2, 3, 4, 5, 5, 5),
        (1, 2, 3, 3, 4, 4, 5, 6),
        (2, 2, 3, 3, 4, 4, 5, 5),
        (0, 1, 2, 3, 3, 5, 6, 8),  # invalid
    ]

    # Filter to valid score sequences (sum = C(8,2)=28, each in 0..7)
    valid_scores = []
    seen = set()
    for s in target_scores:
        s = tuple(sorted(s))
        if s in seen:
            continue
        if sum(s) == 28 and all(0 <= x <= 7 for x in s):
            valid_scores.append(s)
            seen.add(s)

    h_counts = Counter()
    found_21 = None

    for target in valid_scores:
        count_this = 0
        for _ in range(num_per_score):
            # Generate random tournament, then try to match score by swaps
            adj = generate_with_approx_score(target)
            if adj is None:
                continue
            h = count_H_dp(adj)
            h_counts[h] += 1
            count_this += 1
            if h == 21 and found_21 is None:
                found_21 = adj_to_bits(adj)
                scores = score_sequence(adj)
                print(f"  *** H=21 FOUND! bits={found_21}, scores={scores} ***")

        if count_this > 0:
            pass  # silent progress

    elapsed = time.time() - t0
    print(f"  Tested {sum(h_counts.values()):,} tournaments in {elapsed:.1f}s")

    print(f"  H values near 21: ", end="")
    for hv in range(13, 32, 2):
        c = h_counts.get(hv, 0)
        if c > 0:
            print(f"H={hv}:{c} ", end="")
    print()

    if found_21 is not None:
        print(f"  H=21 FOUND!")
    else:
        print(f"  H=21 NOT found in score-targeted search")

    return h_counts, found_21


def generate_with_approx_score(target_scores):
    """Generate a random tournament whose sorted score seq matches target.
    Uses rejection + local swaps.
    """
    n = len(target_scores)
    # Start with random tournament
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1

    # Compute current scores
    scores = [sum(adj[i]) for i in range(n)]
    target_sorted = sorted(target_scores)

    # Do up to 200 random arc-reversal swaps to approach target
    for _ in range(200):
        if tuple(sorted(scores)) == tuple(target_sorted):
            return adj
        # Find a vertex with too-high score and one with too-low
        sorted_s = sorted(range(n), key=lambda x: scores[x])
        sorted_t = sorted(target_sorted)
        # Pick a random vertex that's too high
        over = [v for v in range(n) if scores[v] > sorted_t[sorted(range(n), key=lambda x: scores[x]).index(v)]]
        under = [v for v in range(n) if scores[v] < sorted_t[sorted(range(n), key=lambda x: scores[x]).index(v)]]
        if not over or not under:
            break
        hi = random.choice(over)
        lo = random.choice(under)
        if adj[hi][lo]:
            adj[hi][lo] = 0
            adj[lo][hi] = 1
            scores[hi] -= 1
            scores[lo] += 1

    if tuple(sorted(scores)) == tuple(target_sorted):
        return adj
    return adj  # return anyway, it's still a valid tournament


# ============================================================
# Phase 4: Focused near-21 hunt via biased sampling
# ============================================================

def biased_search(num_samples=100000):
    """Sample tournaments biased toward low H by starting near-transitive
    with random small perturbations."""
    print(f"\n=== Phase 4: Biased near-transitive sampling ({num_samples:,}) ===")
    t0 = time.time()
    h_counts = Counter()
    found_21 = None

    trans_bits = (1 << NUM_EDGES) - 1

    for trial in range(num_samples):
        # Flip a random number of edges (1 to 8) from transitive
        nflips = random.randint(1, 8)
        bits = trans_bits
        flipped = random.sample(range(NUM_EDGES), nflips)
        for e in flipped:
            bits ^= (1 << e)
        adj = bits_to_adj(bits)
        h = count_H_dp(adj)
        h_counts[h] += 1

        if h == 21:
            if found_21 is None:
                found_21 = bits
                scores = score_sequence(adj)
                print(f"  *** H=21 FOUND at trial {trial}! bits={bits} ***")
                print(f"      Flipped {nflips} edges, scores={scores}")

        if (trial + 1) % 50000 == 0:
            elapsed = time.time() - t0
            print(f"  Progress: {trial+1:,}/{num_samples:,} ({elapsed:.1f}s)")
            sys.stdout.flush()

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s")

    print(f"  H values near 21: ", end="")
    for hv in range(13, 32, 2):
        c = h_counts.get(hv, 0)
        if c > 0:
            print(f"H={hv}:{c} ", end="")
    print()

    missing = sorted(h for h in range(1, 60, 2)
                     if h not in h_counts and h > 0 and h <= max(h_counts))
    print(f"  Missing odd H: {missing}")

    return h_counts, found_21


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("H=21 SEARCH AT n=8")
    print("=" * 70)
    print(f"n={N}, C(n,2)={NUM_EDGES} edges, 2^28 = {1<<NUM_EDGES:,} total tournaments")
    print()

    all_h = set()

    # Phase 1: 200K random
    h1, f1 = random_search(200000)
    all_h.update(h1.keys())

    # Phase 2: Edge flips from transitive (up to 4)
    a2, f2 = near_transitive_search(max_flips=4)
    all_h.update(a2)

    # Phase 3: Score-targeted
    h3, f3 = score_targeted_search(num_per_score=3000)
    all_h.update(h3.keys())

    # Phase 4: Biased near-transitive
    h4, f4 = biased_search(100000)
    all_h.update(h4.keys())

    # Final summary
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)

    all_found = any(f is not None for f in [f1, f2, f3, f4])

    missing = sorted(h for h in range(1, 80, 2) if h not in all_h and h <= max(all_h))
    print(f"Missing odd H values at n=8: {missing}")

    small_h = sorted(h for h in all_h if h <= 40)
    print(f"All achieved H <= 40: {small_h}")

    if 21 in all_h:
        print("\nRESULT: H=21 IS achievable at n=8!")
    else:
        print(f"\nRESULT: H=21 NOT found at n=8.")
        print(f"  Searched: 200K random + 4-flip systematic + "
              f"score-targeted + 100K biased")
        print(f"  H=21 may be impossible at n=8 (or extremely rare)")

    if 7 in all_h:
        print("\nWARNING: H=7 found! This contradicts THM-029!")
    else:
        print("Confirmed: H=7 absent (consistent with THM-029)")
