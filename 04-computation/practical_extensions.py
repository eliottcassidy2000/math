"""
practical_extensions.py — opus-2026-04-04-S8

Practical computational extensions using the multilinear theory.

TARGET 1: Extend Σ H sequence to n=12 via succession-free permutation formula
TARGET 2: Compute P(n=7) = nonzero multilinear coefficients
TARGET 3: Compute the vertex-addition transition for n=7 metagraph
"""
import numpy as np
from math import comb, factorial
from collections import defaultdict, Counter
from itertools import permutations
import time

# ========================================================================
# TARGET 1: Extend Σ H using the succession-free formula
# ========================================================================

def count_succession_free_by_bp(n):
    """Count succession-free permutations of {0,...,n-1} by base-path arc count.
    A succession is i immediately before i+1.
    Base-path arc = descent by exactly 1 (i immediately before i-1).
    """
    bp_count = defaultdict(int)
    total = 0
    for perm in permutations(range(n)):
        bp = 0
        valid = True
        for i in range(n-1):
            if perm[i+1] == perm[i] + 1:  # succession
                valid = False
                break
            if perm[i] == perm[i+1] + 1:  # descent by 1 = base path arc
                bp += 1
        if valid:
            bp_count[bp] += 1
            total += 1
    return bp_count, total

def compute_sum_H(n, bp_dist):
    """Compute Σ H from the bp distribution."""
    m = comb(n-1, 2)
    total = 0
    for bp, count in bp_dist.items():
        total += count * (2 ** (m - n + 1 + bp))
    return total

def extend_sum_H():
    print("TARGET 1: Extend Σ H sequence")
    print("=" * 60)

    results = []
    for n in range(3, 13):
        t0 = time.time()
        bp_dist, total_sf = count_succession_free_by_bp(n)
        sum_h = compute_sum_H(n, bp_dist)
        m = comb(n-1, 2)
        mean_h = sum_h / (2 ** m)
        elapsed = time.time() - t0

        results.append((n, total_sf, sum_h, mean_h))
        print(f"  n={n}: SF perms={total_sf}, Σ H={sum_h}, "
              f"E[H]={mean_h:.6f}, time={elapsed:.1f}s")

        if elapsed > 300:  # stop if taking too long
            print(f"  Stopping (time limit)")
            break

    print(f"\nΣ H sequence:")
    for n, sf, sh, mh in results:
        print(f"  n={n}: {sh}")

    print(f"\nE[H] sequence:")
    for n, sf, sh, mh in results:
        print(f"  n={n}: {mh}")

    return results

# ========================================================================
# TARGET 2: Compute P(n=7) nonzero multilinear coefficients
# ========================================================================

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_tournament(n, tiles, bits):
    T = np.zeros((n, n), dtype=np.int8)
    for k in range(n-1):
        T[k+1][k] = 1
    for i, (x, y) in enumerate(tiles):
        a, b = x-1, y-1
        if bits & (1 << i):
            T[b][a] = 1
        else:
            T[a][b] = 1
    return T

def count_hp(T, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    total = 0
    for mask in range(1, 1 << n):
        for v in range(n):
            cnt = dp.get((mask, v), 0)
            if cnt == 0: continue
            if mask == full: total += cnt; continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + cnt
    return total

def compute_P_7():
    """Compute P(7) = number of nonzero multilinear coefficients at n=7."""
    n = 7
    tiles = make_tiles(n)
    m = len(tiles)
    print(f"\nTARGET 2: Compute P({n})")
    print("=" * 60)
    print(f"  m={m} tiles, 2^m={1<<m} tilings")

    # Compute all H values
    t0 = time.time()
    H = np.zeros(1 << m, dtype=np.int64)
    for bits in range(1 << m):
        T = tiling_to_tournament(n, tiles, bits)
        H[bits] = count_hp(T, n)
        if bits % 10000 == 0 and bits > 0:
            elapsed = time.time() - t0
            rate = bits / elapsed
            eta = ((1 << m) - bits) / rate
            print(f"    {bits}/{1<<m} ({100*bits/(1<<m):.0f}%, "
                  f"rate={rate:.0f}/s, ETA={eta:.0f}s)", flush=True)

    t1 = time.time()
    print(f"  H computed in {t1-t0:.1f}s")
    print(f"  H range: [{H.min()}, {H.max()}], mean={H.mean():.2f}")
    print(f"  Σ H = {H.sum()}")

    # Now count nonzero Möbius coefficients
    # For P(7), we need to check all 2^15 = 32768 subsets
    # Möbius inversion: c_S = Σ_{T⊆S} (-1)^{|S\T|} H(T)
    t2 = time.time()
    P = 0
    by_order = defaultdict(int)
    for S in range(1 << m):
        val = 0
        T = S
        while True:
            sign = (-1) ** bin(S ^ T).count('1')
            val += sign * int(H[T])
            if T == 0: break
            T = (T - 1) & S
        if val != 0:
            P += 1
            order = bin(S).count('1')
            by_order[order] += 1

        if S % 10000 == 0 and S > 0:
            print(f"    Möbius: {S}/{1<<m} ({100*S/(1<<m):.0f}%)", flush=True)

    t3 = time.time()
    print(f"\n  P({n}) = {P} nonzero / {1<<m} total ({100*P/(1<<m):.1f}%)")
    print(f"  By order: {dict(sorted(by_order.items()))}")
    print(f"  Möbius computed in {t3-t2:.1f}s")

    d_max = max(by_order.keys()) if by_order else 0
    print(f"  Max degree: {d_max}")

    return P, by_order

# ========================================================================
# TARGET 3: Fast n=7 metagraph class count verification
# ========================================================================

def fast_metagraph_n7():
    """Build a fast version of the n=7 metagraph using score+c3+H as invariant."""
    n = 7
    tiles = make_tiles(n)
    m = len(tiles)
    print(f"\nTARGET 3: Fast n=7 metagraph")
    print("=" * 60)

    t0 = time.time()
    # Use (scores, c3, H) as fast invariant — should distinguish MOST classes
    inv_groups = defaultdict(list)
    for bits in range(1 << m):
        T = tiling_to_tournament(n, tiles, bits)
        h = count_hp(T, n)
        scores = tuple(sorted(int(T[i].sum()) for i in range(n)))
        c3 = 0
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if (T[a][b] and T[b][c] and T[c][a]) or \
                       (T[a][c] and T[c][b] and T[b][a]):
                        c3 += 1
        inv = (scores, c3, h)
        inv_groups[inv].append(bits)

        if bits % 5000 == 0 and bits > 0:
            elapsed = time.time() - t0
            print(f"    {bits}/{1<<m}, {len(inv_groups)} groups, "
                  f"{elapsed:.0f}s", flush=True)

    t1 = time.time()
    print(f"\n  Invariant groups: {len(inv_groups)} (expected ~456 classes)")
    print(f"  Time: {t1-t0:.0f}s")

    # Check: how many groups have size > 1?
    multi = sum(1 for g in inv_groups.values() if len(g) > 1)
    max_size = max(len(g) for g in inv_groups.values())
    print(f"  Groups with multiple tilings: {multi}")
    print(f"  Max group size: {max_size}")

    # The number of invariant groups is an UPPER BOUND on the class count.
    # If it's 456, our invariant is perfect (distinguishes all classes).
    expected = 456  # A000568(7)
    if len(inv_groups) == expected:
        print(f"  *** PERFECT: invariant (scores, c3, H) distinguishes all {expected} classes! ***")
    elif len(inv_groups) > expected:
        print(f"  The invariant OVER-distinguishes: {len(inv_groups)} > {expected}")
        print(f"  (Some groups need to be merged — invariant is too fine)")
    else:
        print(f"  The invariant UNDER-distinguishes: {len(inv_groups)} < {expected}")
        print(f"  (Need finer invariant to separate all classes)")

    # H-value statistics from the groups
    h_classes = defaultdict(int)
    for inv, members in inv_groups.items():
        h_classes[inv[2]] += 1

    print(f"\n  Classes per H value:")
    for h in sorted(h_classes.keys()):
        print(f"    H={h}: {h_classes[h]} classes")

    # Tiling-fibration check: size × |Aut| = H
    print(f"\n  Tiling-fibration identity size × |Aut| = H:")
    violations = 0
    for inv, members in list(inv_groups.items())[:20]:
        h = inv[2]
        size = len(members)
        aut = h / size if size > 0 else 0
        is_int = abs(aut - round(aut)) < 0.001
        if not is_int:
            violations += 1
            if violations <= 5:
                print(f"    scores={inv[0]}, c3={inv[1]}, H={h}, size={size}, "
                      f"|Aut|={aut:.3f} {'✓' if is_int else '✗'}")

    if violations == 0:
        print(f"    All {min(20, len(inv_groups))} checked: |Aut| is integer ✓")
    else:
        print(f"    {violations} violations (invariant may over-split)")

    return inv_groups

if __name__ == "__main__":
    # Target 1: Extend Σ H
    results = extend_sum_H()

    # Target 2: P(7)
    P7, by_order = compute_P_7()

    # Target 3: Fast n=7 metagraph
    groups = fast_metagraph_n7()
