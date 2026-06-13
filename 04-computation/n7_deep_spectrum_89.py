#!/usr/bin/env python3
"""
n7_deep_spectrum_89.py — opus-2026-03-14-S89
Deep analysis of n=7 tournament H-spectrum.

Key questions:
1. What is the EXACT H-spectrum at n=7? (all 2^21 tournaments)
2. Is H=63 permanently forbidden?
3. What is 131 (numerator of Var/Mean² = 131/504)?
4. The deviation sequence 0, 0, 1/60, 2/45, 37/504
5. Regular tournament H-values at n=7

Strategy: n=7 has 2^21 = 2097152 tournaments.
For each, compute H = number of Hamiltonian paths.
This requires checking 7! = 5040 permutations per tournament.
Total ops: ~10^10. Too slow naively.

OPTIMIZATION: Use dynamic programming for Hamiltonian path counting.
DP on subset mask + endpoint: dp[mask][v] = number of Hamiltonian paths
ending at v using vertices in mask.
Complexity: O(2^n * n * n) per tournament = O(2^7 * 7 * 7) = 6272 per tournament.
Total: 6272 * 2^21 ≈ 1.3 * 10^10. Still slow but feasible with bitwise ops.

Better: for n=7, 2^7 = 128 masks, so dp has 128*7 = 896 entries per tournament.
Each transition: check all pairs. Total per tournament: ~128*7*7 ≈ 6272 steps.
Total: 6272 * 2097152 ≈ 1.3 * 10^10 operations. ~30 min in Python.

Let's use C-like optimization or numpy. Actually, let me try multiprocessing.
"""

import sys
import time
from collections import Counter
from fractions import Fraction
import multiprocessing as mp

def compute_H_dp(adj_bits, n):
    """Compute H using Hamiltonian path DP.
    adj_bits[i] = bitmask of vertices that i beats.
    dp[mask][v] = number of Ham paths using vertices in mask, ending at v.
    """
    dp = [[0]*n for _ in range(1 << n)]

    # Base: single vertex
    for v in range(n):
        dp[1 << v][v] = 1

    # Fill
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            count = dp[mask][v]
            # Extend to vertex u not in mask, where v->u (v beats u)
            remaining = ((1 << n) - 1) & ~mask
            targets = adj_bits[v] & remaining
            u = targets
            while u:
                bit = u & (-u)  # lowest set bit
                idx = bit.bit_length() - 1
                dp[mask | bit][idx] += count
                u &= u - 1

    # Sum over all endpoints in full mask
    full = (1 << n) - 1
    return sum(dp[full])


def process_chunk(args):
    """Process a range of tournament indices."""
    start, end, n = args
    m = n * (n - 1) // 2
    h_counts = Counter()

    # Precompute edge list: edge i connects (u, v) where u < v
    edges = []
    for u in range(n):
        for v in range(u+1, n):
            edges.append((u, v))

    for idx in range(start, end):
        # Build adjacency from tournament index
        adj_bits = [0] * n
        bits = idx
        for e, (u, v) in enumerate(edges):
            if bits & (1 << e):
                adj_bits[u] |= (1 << v)  # u beats v
            else:
                adj_bits[v] |= (1 << u)  # v beats u

        h = compute_H_dp(adj_bits, n)
        h_counts[h] += 1

    return h_counts


def main():
    n = 7
    m = n * (n - 1) // 2  # 21
    total = 1 << m  # 2^21 = 2097152

    print("="*70)
    print(f"EXACT H-SPECTRUM AT n={n}")
    print(f"opus-2026-03-14-S89")
    print("="*70)
    print(f"\n  n={n}, m={m}, total tournaments = {total}")
    print(f"  Using Hamiltonian path DP (O(2^n * n^2) per tournament)")

    # Use multiprocessing
    num_workers = mp.cpu_count()
    print(f"  Using {num_workers} workers")

    chunk_size = total // num_workers
    chunks = []
    for i in range(num_workers):
        start = i * chunk_size
        end = (i + 1) * chunk_size if i < num_workers - 1 else total
        chunks.append((start, end, n))

    t0 = time.time()

    with mp.Pool(num_workers) as pool:
        results = pool.map(process_chunk, chunks)

    # Merge results
    h_counts = Counter()
    for r in results:
        h_counts.update(r)

    elapsed = time.time() - t0
    print(f"\n  Done in {elapsed:.1f}s")

    # Sort by H value
    h_values = sorted(h_counts.keys())

    print(f"\n  Total distinct H values: {len(h_values)}")
    print(f"  H values: {h_values}")

    # Check forbidden values
    print(f"\n  FORBIDDEN VALUE CHECK:")
    for v in [7, 21, 35, 39, 63, 73]:
        count = h_counts.get(v, 0)
        status = "ACHIEVABLE" if count > 0 else "FORBIDDEN"
        print(f"    H={v}: count={count} — {status}")

    # Full spectrum
    print(f"\n  FULL H-SPECTRUM (H: count):")
    for h in h_values:
        print(f"    H={h:5d}: {h_counts[h]:8d}")

    # Statistics
    print(f"\n" + "="*70)
    print(f"PART 2: STATISTICS AND THE NUMBER 131")
    print("="*70)

    sum_h = sum(h * c for h, c in h_counts.items())
    sum_h2 = sum(h**2 * c for h, c in h_counts.items())

    mean = Fraction(sum_h, total)
    mean_h2 = Fraction(sum_h2, total)
    var = mean_h2 - mean**2
    var_over_mean2 = var / mean**2

    print(f"\n  sum(H) = {sum_h}")
    print(f"  sum(H²) = {sum_h2}")
    print(f"  Mean = {mean} = {float(mean):.6f}")
    print(f"  E[H²] = {mean_h2} = {float(mean_h2):.6f}")
    print(f"  Var = {var} = {float(var):.6f}")
    print(f"  Var/Mean² = {var_over_mean2} = {float(var_over_mean2):.10f}")

    # Check if 131/504
    print(f"\n  Var/Mean² = {var_over_mean2.numerator}/{var_over_mean2.denominator}")
    print(f"  Expected: 131/504")
    print(f"  Match: {var_over_mean2 == Fraction(131, 504)}")

    # Analyze 131
    print(f"\n  WHAT IS 131?")
    print(f"    131 is prime: {all(131 % i != 0 for i in range(2, 12))}")
    print(f"    131 = 128 + 3 = 2^7 + 3")
    print(f"    131 mod 3 = {131 % 3}")
    print(f"    131 mod 7 = {131 % 7}")
    print(f"    131 mod 6 = {131 % 6}")

    # Analyze 504
    print(f"\n  WHAT IS 504?")
    print(f"    504 = 7 × 72 = 7 × 8 × 9")
    print(f"    504 = 7! / 10 = 5040 / 10")
    print(f"    504 = C(9,2) × 8 = 36 × 14")

    # Factor 504
    n504 = 504
    factors = []
    for p in [2, 3, 5, 7, 11, 13]:
        while n504 % p == 0:
            factors.append(p)
            n504 //= p
    print(f"    504 = {' × '.join(map(str, factors))}")
    print(f"    504 IS tribonacci number T(14)? Let's check")

    trib = [0, 0, 1]
    while trib[-1] < 600:
        trib.append(trib[-1] + trib[-2] + trib[-3])
    print(f"    Tribonacci: {[t for t in trib if t > 0 and t < 600]}")
    print(f"    504 in tribonacci: {504 in trib}")

    # Deviation sequence
    print(f"\n" + "="*70)
    print(f"PART 3: DEVIATION SEQUENCE")
    print("="*70)

    var_ratios = {
        3: Fraction(1, 3),
        4: Fraction(1, 3),
        5: Fraction(19, 60),
        6: Fraction(13, 45),
        7: Fraction(131, 504),
    }

    deviations = {}
    for nn, vr in var_ratios.items():
        dev = Fraction(1, 3) - vr
        deviations[nn] = dev
        print(f"\n  n={nn}: Var/Mean² = {vr} = {float(vr):.10f}")
        print(f"         Deviation from 1/3 = {dev} = {float(dev):.10f}")

    # Look for pattern in deviations
    print(f"\n  Deviation numerators/denominators:")
    for nn, dev in deviations.items():
        print(f"    n={nn}: {dev.numerator}/{dev.denominator}")

    # n=3: 0/1
    # n=4: 0/1
    # n=5: 1/60
    # n=6: 2/45 = 2/45
    # n=7: 37/504

    # Check denominators
    print(f"\n  Denominators: {[deviations[nn].denominator for nn in [3,4,5,6,7]]}")
    print(f"  60 = 3·4·5 = 5!/2")
    print(f"  45 = 9·5 = 5·9")
    print(f"  504 = 7·8·9")

    # Could denominator be related to n(n-1)(n-2)/something?
    for nn in [5, 6, 7]:
        prod = nn * (nn-1) * (nn-2)
        print(f"  n={nn}: n(n-1)(n-2) = {prod}, denom = {deviations[nn].denominator}, ratio = {Fraction(prod, deviations[nn].denominator)}")

    # n=5: 60/60 = 1
    # n=6: 120/45 = 8/3
    # n=7: 210/504 = 5/12
    # Hmm, not clean.

    # Try P(n,2) = n(n-1)
    for nn in [5, 6, 7]:
        prod = nn * (nn-1)
        print(f"  n={nn}: n(n-1) = {prod}, denom = {deviations[nn].denominator}, ratio = {Fraction(prod, deviations[nn].denominator)}")

    # Try m = n(n-1)/2
    for nn in [5, 6, 7]:
        mm = nn * (nn-1) // 2
        print(f"  n={nn}: m = {mm}, denom = {deviations[nn].denominator}, ratio = {Fraction(deviations[nn].denominator, mm)}")

    # Regularities in numerator sequence: 0, 0, 1, 2, 37
    print(f"\n  Numerator sequence: 0, 0, 1, 2, 37")
    print(f"    37 is prime")
    print(f"    37 = 36 + 1 = 6² + 1")
    print(f"    37 mod 6 = {37 % 6}")
    print(f"    37 = first irregular prime > 5")

    # H-spectrum structure
    print(f"\n" + "="*70)
    print(f"PART 4: H-SPECTRUM STRUCTURE AT n=7")
    print("="*70)

    # Which values are odd? (all should be)
    odd_h = [h for h in h_values if h % 2 == 1]
    even_h = [h for h in h_values if h % 2 == 0]
    print(f"\n  Odd H values: {len(odd_h)}")
    print(f"  Even H values: {len(even_h)} (should be 0)")

    # Gap analysis
    all_odd = list(range(1, max(h_values)+1, 2))
    missing = [h for h in all_odd if h not in h_counts]
    print(f"\n  Missing odd values in [1, {max(h_values)}]: {missing}")
    print(f"  Number of missing values: {len(missing)}")

    # Forbidden values are those missing
    print(f"\n  PERMANENTLY FORBIDDEN at n=7 (values that never appear):")
    for h in missing:
        print(f"    H={h}", end="")
        # Check cyclotomic/repunit status
        if h == 7:
            print(f" = 111₂ = Φ₃(2) = Fano", end="")
        elif h == 21:
            print(f" = 111₄ = Φ₃(4) = Baer", end="")
        elif h == 63:
            print(f" = 111111₂ = 2⁶-1 = Mersenne M₆", end="")
        print()

    # Mod 3 distribution
    print(f"\n  H mod 3 distribution:")
    mod3 = Counter()
    for h, c in h_counts.items():
        mod3[h % 3] += c
    for r in [0, 1, 2]:
        print(f"    H ≡ {r} mod 3: {mod3.get(r, 0)} tournaments ({mod3.get(r,0)/total*100:.2f}%)")

    # Mod 7 distribution
    print(f"\n  H mod 7 distribution:")
    mod7 = Counter()
    for h, c in h_counts.items():
        mod7[h % 7] += c
    for r in range(7):
        print(f"    H ≡ {r} mod 7: {mod7.get(r, 0)} tournaments ({mod7.get(r,0)/total*100:.2f}%)")

    # GCD of all counts
    from math import gcd
    from functools import reduce
    all_counts = [h_counts[h] for h in h_values]
    g = reduce(gcd, all_counts)
    print(f"\n  GCD of all H-spectrum counts: {g}")

    # Factor g
    gg = g
    gfactors = []
    for p in range(2, g+1):
        while gg % p == 0:
            gfactors.append(p)
            gg //= p
    if gfactors:
        print(f"  {g} = {' × '.join(map(str, gfactors))}")

    # Regular tournaments at n=7
    print(f"\n" + "="*70)
    print(f"PART 5: REGULAR TOURNAMENTS AT n=7")
    print("="*70)

    edges = []
    for u in range(n):
        for v in range(u+1, n):
            edges.append((u, v))

    regular_h = Counter()
    regular_count = 0
    for idx in range(total):
        # Check if regular: each vertex has out-degree 3
        out_deg = [0] * n
        bits = idx
        for e_idx, (u, v) in enumerate(edges):
            if bits & (1 << e_idx):
                out_deg[u] += 1
            else:
                out_deg[v] += 1

        if all(d == 3 for d in out_deg):
            regular_count += 1
            # Compute H
            adj_bits = [0] * n
            bits = idx
            for e_idx, (u, v) in enumerate(edges):
                if bits & (1 << e_idx):
                    adj_bits[u] |= (1 << v)
                else:
                    adj_bits[v] |= (1 << u)
            h = compute_H_dp(adj_bits, n)
            regular_h[h] += 1

    print(f"\n  Number of regular tournaments at n=7: {regular_count}")
    print(f"  Regular H values: {sorted(regular_h.keys())}")
    for h in sorted(regular_h.keys()):
        print(f"    H={h}: {regular_h[h]} regular tournaments")

    # Score sequences
    print(f"\n" + "="*70)
    print(f"PART 6: H=63 OBSTRUCTION ANALYSIS")
    print("="*70)

    # 63 = 2^6 - 1 = 111111₂
    # Is it achievable? If count is 0, it's forbidden
    c63 = h_counts.get(63, 0)
    print(f"\n  H=63 count: {c63}")
    if c63 == 0:
        print(f"  H=63 is PERMANENTLY FORBIDDEN at n=7!")
        print(f"  63 = 2⁶ - 1 = Mersenne number M₆")
        print(f"  63 = 7 × 9 = Φ₃(2) × 3²")
        print(f"  63 = 111111₂ (six ones in binary)")
        print(f"  63 in base 4: {63} = 333₄ (repdigit!)")
        # 63 = 3*16 + 3*4 + 3 = 48+12+3 = 63 ✓
        print(f"  63 in base 6: {63 // 6}{63 % 6} = 10{63 % 6}₆")
        # Actually: 63 = 10*6 + 3 = 103₆
        print(f"  Actually 63 = 10×6 + 3 = 103₆")

        # Check neighbors
        for delta in [-4, -2, 0, 2, 4]:
            h_check = 63 + delta
            if h_check > 0 and h_check % 2 == 1:
                c = h_counts.get(h_check, 0)
                print(f"    H={h_check}: count={c}")
    else:
        print(f"  H=63 IS achievable at n=7")

    # Check if 63 appears at n=8 (can't enumerate, but note the pattern)
    print(f"\n  REPUNIT/REPDIGIT FORBIDDEN VALUES:")
    print(f"    7 = 111₂ = R₃(2) — FORBIDDEN")
    print(f"    21 = 111₄ = R₃(4) — FORBIDDEN")
    print(f"    63 = 333₄ = 3·R₃(4) — {'FORBIDDEN' if c63==0 else 'ACHIEVABLE'}")
    print(f"    63 = 111111₂ = R₆(2) — {'FORBIDDEN' if c63==0 else 'ACHIEVABLE'}")
    print(f"    63 = 2⁶-1 = Mersenne — {'FORBIDDEN' if c63==0 else 'ACHIEVABLE'}")

    # Divisibility by 7
    print(f"\n  H values divisible by 7:")
    div7 = [(h, h_counts[h]) for h in h_values if h % 7 == 0]
    for h, c in div7:
        print(f"    H={h}: count={c}")

    # Divisibility by 21
    print(f"\n  H values divisible by 21:")
    div21 = [(h, h_counts[h]) for h in h_values if h % 21 == 0]
    for h, c in div21:
        print(f"    H={h}: count={c}")

    # Divisibility by 63
    print(f"\n  H values divisible by 63:")
    div63 = [(h, h_counts[h]) for h in h_values if h % 63 == 0]
    for h, c in div63:
        print(f"    H={h}: count={c}")

    # Analyze the forbidden set
    print(f"\n" + "="*70)
    print(f"PART 7: COMPLETE FORBIDDEN VALUE ANALYSIS")
    print("="*70)

    print(f"\n  All missing odd values up to max(H)={max(h_values)}:")
    for h in missing:
        # Binary representation
        binary = bin(h)[2:]
        # Base 4
        b4 = ""
        hh = h
        while hh > 0:
            b4 = str(hh % 4) + b4
            hh //= 4
        # Base 6
        b6 = ""
        hh = h
        while hh > 0:
            b6 = str(hh % 6) + b6
            hh //= 6

        # Factorization
        facts = []
        hh = h
        for p in range(2, h+1):
            while hh % p == 0:
                facts.append(p)
                hh //= p
            if hh == 1:
                break
        fact_str = " × ".join(map(str, facts)) if facts else "1"

        print(f"    H={h:5d} = {binary}₂ = {b4}₄ = {b6}₆ = {fact_str}")

    print(f"\n" + "="*70)
    print(f"DONE — n=7 EXACT H-SPECTRUM")
    print("="*70)


if __name__ == "__main__":
    main()
