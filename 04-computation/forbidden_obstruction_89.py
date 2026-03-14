#!/usr/bin/env python3
"""
forbidden_obstruction_89.py -- opus-2026-03-14-S89

Deep investigation of WHY H=63 is forbidden at n=7.
Uses actual exhaustive FWHT data to understand the lattice obstruction.

Strategy: For each H value, compute the exact Fourier sign pattern constraints.
The forbidden values are those where the sign constraints are inconsistent.
"""

from itertools import combinations
from fractions import Fraction
import sys


def compute_H_dp(adj_bits, n):
    """DP for Hamiltonian path count."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            count = dp[mask][v]
            remaining = ((1 << n) - 1) & ~mask
            targets = adj_bits[v] & remaining
            u = targets
            while u:
                bit = u & (-u)
                idx = bit.bit_length() - 1
                dp[mask | bit][idx] += count
                u &= u - 1
    full = (1 << n) - 1
    return sum(dp[full])


def tournament_from_bits(bits, n):
    """Convert integer encoding to adjacency list."""
    m = n * (n - 1) // 2
    adj = [0] * n
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits & (1 << idx):
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
            idx += 1
    return adj


def main():
    n = 7
    m = n * (n - 1) // 2  # 21
    total = 1 << m  # 2^21 = 2097152

    print("="*70)
    print(f"FORBIDDEN VALUE OBSTRUCTION ANALYSIS at n={n}")
    print("opus-2026-03-14-S89")
    print("="*70)

    # Step 1: Enumerate all tournaments, collect H values
    print(f"\nStep 1: Enumerate all 2^{m} = {total} tournaments...")

    h_to_count = {}
    h_to_examples = {}  # Store a few examples for each H

    for bits in range(total):
        if bits % 500000 == 0:
            print(f"  Progress: {bits}/{total} ({100*bits/total:.1f}%)", flush=True)
        adj = tournament_from_bits(bits, n)
        h = compute_H_dp(adj, n)
        h_to_count[h] = h_to_count.get(h, 0) + 1
        if h not in h_to_examples or len(h_to_examples[h]) < 3:
            h_to_examples.setdefault(h, []).append(bits)

    all_h = sorted(h_to_count.keys())
    print(f"\n  Found {len(all_h)} distinct H values")
    print(f"  Range: [{all_h[0]}, {all_h[-1]}]")

    # Step 2: Identify forbidden values (odd integers in range not achieved)
    print(f"\nStep 2: Forbidden values...")
    min_h, max_h = all_h[0], all_h[-1]
    all_odd = set(range(min_h, max_h + 1, 2))
    achieved = set(all_h)
    forbidden = sorted(all_odd - achieved)
    print(f"  Odd integers in [{min_h}, {max_h}]: {len(all_odd)}")
    print(f"  Achieved: {len(achieved)}")
    print(f"  Forbidden: {len(forbidden)}")
    print(f"  Forbidden values: {forbidden}")

    # Step 3: For each forbidden value, analyze its neighbors
    print(f"\nStep 3: Neighbor analysis of forbidden values...")
    for fv in forbidden:
        below = fv - 2
        above = fv + 2
        cb = h_to_count.get(below, 0)
        ca = h_to_count.get(above, 0)
        print(f"  H={fv}: below H={below} ({cb} tours), above H={above} ({ca} tours)")

    # Step 4: H values mod small primes
    print(f"\nStep 4: H values modular arithmetic...")
    for p in [2, 3, 4, 5, 7, 8, 9, 16]:
        residues = {}
        for h, c in h_to_count.items():
            r = h % p
            residues[r] = residues.get(r, 0) + c
        print(f"\n  H mod {p}:")
        for r in sorted(residues.keys()):
            count = residues[r]
            is_forb = any(f % p == r for f in forbidden)
            marker = " [FORBIDDEN RESIDUE]" if r in [f % p for f in forbidden] else ""
            print(f"    {r}: {count} tournaments{marker}")

    # Step 5: Distribution around forbidden values
    print(f"\nStep 5: Distribution near forbidden values...")
    for fv in forbidden:
        nearby = {}
        for delta in range(-10, 11, 2):
            h = fv + delta
            nearby[h] = h_to_count.get(h, 0)
        print(f"\n  Around H={fv}:")
        for h in sorted(nearby.keys()):
            bar = '#' * min(nearby[h] // 1000, 50)
            marker = " <<<FORBIDDEN" if h == fv else ""
            print(f"    H={h:3d}: {nearby[h]:7d} {bar}{marker}")

    # Step 6: Binary structure of forbidden values
    print(f"\nStep 6: Binary/number-theoretic structure...")
    for fv in forbidden:
        print(f"\n  H={fv}:")
        print(f"    Binary: {bin(fv)}")
        print(f"    = {fv}")
        # Factor
        factors = []
        x = fv
        for p in range(2, fv + 1):
            while x % p == 0:
                factors.append(p)
                x //= p
            if x == 1:
                break
        print(f"    Factors: {' * '.join(map(str, factors))}")
        print(f"    Mod 7: {fv % 7}")
        print(f"    Mod 8: {fv % 8}")
        print(f"    Mod 9: {fv % 9}")
        print(f"    Mod 16: {fv % 16}")

    # Step 7: Score distribution summary
    print(f"\nStep 7: Full H-spectrum with counts...")
    print(f"  {'H':>5} {'count':>8} {'cumul%':>8}")
    cumul = 0
    for h in all_h:
        cumul += h_to_count[h]
        pct = 100 * cumul / total
        if h_to_count[h] >= 1000 or h in forbidden or h <= 10 or h >= max_h - 5:
            print(f"  {h:5d} {h_to_count[h]:8d} {pct:7.2f}%")

    # Step 8: Repunit/repdigit analysis
    print(f"\nStep 8: Repunit structure of forbidden values...")
    for fv in forbidden:
        for base in range(2, 10):
            digits = []
            x = fv
            while x > 0:
                digits.append(x % base)
                x //= base
            digits.reverse()
            rep = ''.join(map(str, digits))
            # Check if repunit (all 1s) or repdigit (all same)
            if len(set(rep)) == 1:
                print(f"    H={fv} in base {base}: {''.join(map(str, digits))} (REPDIGIT)")

    print(f"\n{'='*70}")
    print("DONE -- FORBIDDEN OBSTRUCTION ANALYSIS")
    print("="*70)


if __name__ == "__main__":
    main()
