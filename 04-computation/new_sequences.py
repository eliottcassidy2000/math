"""
new_sequences.py — opus-2026-04-04-S5

Generate new sequences from the multilinear polynomial theory.

SEQUENCES TO COMPUTE:
1. P(n) = number of nonzero multilinear coefficients of H(t) at n
2. E[H(n)] = mean H over all 2^m tilings (= mean HP count with fixed base path)
3. V(n) = number of valid directed odd cycles (contributing to OCF)
4. Z(n) = number of zero multilinear coefficients
5. D(n) = maximum degree of multilinear polynomial (= 2*floor((n-1)/2))
6. B(n) = number of boundary-touching nonzero coefficients in n→n+2 decomposition
"""
import numpy as np
from math import comb
from collections import defaultdict
import time

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_tournament(n, tiles, bits):
    T = np.zeros((n, n), dtype=int)
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
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def compute_all_H(n):
    tiles = make_tiles(n)
    m = len(tiles)
    H = np.zeros(1 << m, dtype=np.int64)
    for bits in range(1 << m):
        T = tiling_to_tournament(n, tiles, bits)
        H[bits] = count_hp(T, n)
    return H, tiles, m

def mobius_count(H, m):
    """Count nonzero Mobius coefficients without storing all of them."""
    nonzero = 0
    for S in range(1 << m):
        val = 0
        T = S
        while True:
            sign = (-1) ** bin(S ^ T).count('1')
            val += sign * int(H[T])
            if T == 0: break
            T = (T - 1) & S
        if val != 0:
            nonzero += 1
    return nonzero

def mean_H_formula(n):
    """Compute E[H] from the multilinear linear terms.

    E[H] = Σ_S c_S (1/2)^|S| = c_0 + Σ_k c_k/2 + higher order terms
    = 1 + Σ_k 2^(s_k-1)/2 + ... = 1 + Σ_k 2^(s_k-2) + ...

    Linear contribution: 1 + Σ_{tiles k} 2^(skip_k - 2)
    = 1 + Σ_{s=2}^{n-1} (n-1-s) · 2^(s-2)

    Wait: tiles with skip s at n: for each s, tiles (x,y) with x-y=s.
    These satisfy y ≥ 1, x ≤ n, x = y+s. So y ranges from 1 to n-s.
    But also x = y+s ≤ n → y ≤ n-s. And y ≥ 1.
    Count = n - s (for s=2,...,n-1).

    Linear contribution to E[H] = Σ_{s=2}^{n-1} (n-s) · 2^(s-2)
    """
    linear = sum((n - s) * (2 ** (s - 2)) for s in range(2, n))
    return 1 + linear

def count_valid_cycles(n, max_len=None):
    """Count valid directed odd cycles with non-None indicators."""
    from itertools import combinations, permutations
    tiles = make_tiles(n)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    if max_len is None:
        max_len = n if n % 2 == 1 else n - 1

    total = 0
    by_length = {}
    for L in range(3, max_len + 1, 2):
        count = 0
        for verts in combinations(range(n), L):
            for perm in permutations(verts):
                # Canonical rotation
                rotations = [perm[i:] + perm[:i] for i in range(L)]
                if perm != min(rotations):
                    continue
                # Check validity
                valid = True
                for i in range(L):
                    u, v = perm[i], perm[(i+1) % L]
                    if u == v + 1:
                        continue  # base path ok
                    if v == u + 1:
                        valid = False; break  # reverse base path
                    high, low = max(u,v), min(u,v)
                    if (high+1, low+1) not in tile_to_idx:
                        valid = False; break
                if valid:
                    count += 1
        by_length[L] = count
        total += count

    return total, by_length

def compute_sequences():
    print("COMPUTING NEW SEQUENCES FROM MULTILINEAR THEORY")
    print("=" * 60)

    results = {
        'n': [],
        'm': [],
        'P': [],      # nonzero multilinear coefficients
        'Z': [],      # zero coefficients
        'mean_H': [], # actual mean H
        'mean_H_linear': [],  # linear approximation
        'D': [],      # max degree
        'V': [],      # valid cycles
        'H_sum': [],  # sum of all H values
    }

    for n in range(3, 9):
        tiles = make_tiles(n)
        m = len(tiles)
        d_max = 2 * ((n-1) // 2)
        mean_linear = mean_H_formula(n)

        print(f"\nn={n}, m={m} tiles")

        t0 = time.time()

        if m <= 21:  # feasible to enumerate all tilings
            H, _, _ = compute_all_H(n)
            actual_mean = float(np.mean(H))
            h_sum = int(np.sum(H))
            P = mobius_count(H, m)
            Z = (1 << m) - P
        else:
            actual_mean = None
            h_sum = None
            P = None
            Z = None

        t1 = time.time()

        # Count valid cycles (can be slow for large n)
        if n <= 7:
            V, by_len = count_valid_cycles(n)
            print(f"  Valid cycles: {V} ({by_len})")
        elif n == 8:
            # Only count 3-cycles (fast)
            V3 = 0
            for a in range(n):
                for b in range(a+1,n):
                    for c in range(b+1,n):
                        # Check both orientations for valid indicator
                        for perm in [(a,b,c), (a,c,b)]:
                            valid = True
                            for i in range(3):
                                u,v = perm[i], perm[(i+1)%3]
                                if v == u+1: valid = False; break
                            if valid:
                                V3 += 1
            V = V3  # approximate (3-cycles only)
            by_len = {3: V3}
            print(f"  Valid 3-cycles: {V3} (5+ cycles not enumerated)")
        else:
            V = None
            by_len = {}

        results['n'].append(n)
        results['m'].append(m)
        results['P'].append(P)
        results['Z'].append(Z)
        results['mean_H'].append(actual_mean)
        results['mean_H_linear'].append(mean_linear)
        results['D'].append(d_max)
        results['V'].append(V)
        results['H_sum'].append(h_sum)

        if actual_mean is not None:
            print(f"  H values: sum={h_sum}, mean={actual_mean:.2f}")
        print(f"  Linear mean approx: {mean_linear}")
        print(f"  Max degree: {d_max}")
        if P is not None:
            print(f"  Nonzero coeffs P(n): {P} / {1<<m}")
            print(f"  Zero coeffs Z(n): {Z}")
        print(f"  Time: {t1-t0:.1f}s")

    # Print sequences
    print(f"\n{'='*60}")
    print(f"SEQUENCES FOR OEIS")
    print(f"{'='*60}")

    print(f"\nn:           {results['n']}")
    print(f"m=C(n-1,2):  {results['m']}")
    print(f"P(n) nonzero: {results['P']}")
    print(f"Z(n) zero:   {results['Z']}")
    print(f"D(n) degree: {results['D']}")
    print(f"V(n) cycles: {results['V']}")

    # Sum H = Σ H(t) over all 2^m tilings
    print(f"Σ H:         {results['H_sum']}")

    # Mean H = Σ H / 2^m
    mean_exact = [h//(1<<m) if h is not None else None
                  for h, m in zip(results['H_sum'], results['m'])]
    # Actually mean is not integer. Let me compute Σ H
    print(f"Mean H:      {[f'{x:.2f}' if x else '?' for x in results['mean_H']]}")
    print(f"Linear approx: {results['mean_H_linear']}")

    # Check if Σ H is a known sequence
    print(f"\nΣ_{'{t}'} H(t) = total HP count over all tilings:")
    for i, n in enumerate(results['n']):
        if results['H_sum'][i] is not None:
            print(f"  n={n}: Σ H = {results['H_sum'][i]}")

    # The linear approximation to mean H
    print(f"\nLinear mean E_lin[H] = 1 + Σ (n-s)·2^(s-2) for s=2..n-1:")
    for n in range(3, 15):
        ml = mean_H_formula(n)
        print(f"  n={n}: E_lin = {ml}")

    # Check: the exact mean H involves higher order corrections
    # Mean H = E_lin + quadratic correction + ...
    if results['mean_H'][2] is not None:  # n=5
        for i, n in enumerate(results['n']):
            if results['mean_H'][i] is not None:
                correction = results['mean_H'][i] - results['mean_H_linear'][i]
                print(f"  n={n}: exact={results['mean_H'][i]:.4f}, "
                      f"linear={results['mean_H_linear'][i]}, "
                      f"correction={correction:.4f}")

    # P(n) sequence: 2, 6, 35, 200, 1782, ...
    print(f"\nP(n) = nonzero multilinear coefficients:")
    print(f"  {[p for p in results['P'] if p is not None]}")
    print(f"  (Check OEIS for this sequence)")

    # V(n) = valid directed odd cycles with base-path-compatible indicators
    print(f"\nV(n) = valid OCF cycles:")
    print(f"  {[v for v in results['V'] if v is not None]}")

if __name__ == "__main__":
    compute_sequences()
