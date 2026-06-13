#!/usr/bin/env python3
"""
fourier_n7_levels_89.py — opus-2026-03-14-S89
Fourier level decomposition at n=7.

At n=5 (m=10): only levels 2 and 4 are nonzero.
  E_2/E_0 = 3/10, E_4/E_0 = 1/60

At n=7 (m=21): which levels contribute?
  We know: N_2(7)=0, N_4(7)=315, N_6(7)=23604
  E_4/E_0 = 1/560, E_6/E_0 = 281/75600
  Total from E_4+E_6 = 1/560 + 281/75600 = 135/75600 + 281/75600 = 416/75600
  = 26/4725 ≈ 0.00550

  But Var/Mean² = 131/504 ≈ 0.25992
  So the N-formula is NOT giving the right result for "level-2" because
  N_2(7)=0 but the actual Fourier level-2 energy is large.

  The confusion: the "N_{2k}" master formula from level4_counting.out
  uses a DIFFERENT definition than vertex-covering edge subsets.
  Let me compute the ACTUAL Fourier level decomposition.

Strategy: For n=7, m=21, we can't enumerate all 2^21 Walsh coefficients
(there are 2^21 of them). But we CAN compute the energy by level.

E_k = Σ_{|S|=k} ĉ_S²

Using Parseval: Σ_S ĉ_S² = (1/2^m) Σ_T H(T)²
So Σ_{k≥0} E_k = (1/2^m) Σ H²

And E_0 = ĉ_∅² = Mean²

We need: for each level k, sum of ĉ_S² over all S with |S|=k.

Trick: E_k = (1/2^{2m}) Σ_{|S|=k} (Σ_T (-1)^{S·T} H(T))²
     = (1/2^{2m}) Σ_{|S|=k} Σ_{T1,T2} (-1)^{S·(T1⊕T2)} H(T1)H(T2)
     = (1/2^{2m}) Σ_{T1,T2} H(T1)H(T2) Σ_{|S|=k} (-1)^{S·(T1⊕T2)}

The inner sum: Σ_{|S|=k} (-1)^{S·D} where D = T1⊕T2
Let D have Hamming weight w. Then:
Σ_{|S|=k} (-1)^{|S∩D|} = Σ_{j=0}^{min(k,w)} C(w,j)(-1)^j C(m-w, k-j)

This is a Krawtchouk polynomial!
K_k(w; m) = Σ_j (-1)^j C(w,j) C(m-w, k-j)

So E_k = (1/2^{2m}) Σ_{T1,T2} H(T1)H(T2) K_k(d(T1,T2); m)
where d(T1,T2) = Hamming distance = popcount(T1⊕T2).

This means: E_k = (1/2^{2m}) Σ_{w=0}^{m} K_k(w; m) * C(w)
where C(w) = Σ_{d(T1,T2)=w} H(T1)H(T2) = correlation at distance w.

We can compute C(w) by:
C(w) = Σ_T1 H(T1) * (Σ_{T2: d(T1,T2)=w} H(T2))

Or more efficiently:
C(w) = Σ_T1 H(T1) * S_w(T1)
where S_w(T1) = Σ_{T2: d(T1,T2)=w} H(T2)

Computing S_w for all T1 is expensive (O(2^m * C(m,w)) per w).
Total: O(2^m * 2^m) which is 2^42 — too slow.

BETTER APPROACH: Compute the distance profile directly.
For each tournament T, we know H(T). We need:
Corr(w) = (1/2^m) Σ_T H(T) * average_{T': d(T,T')=w} H(T')
         = (1/2^{2m}) Σ_{d(T,T')=w} H(T)H(T')

Alternative: just compute the Walsh transform level by level using
the relationship between distance spectrum and Walsh levels.

Actually, the SIMPLEST approach for n=7:
1. Compute H for all 2^21 tournaments
2. For each pair of distances w, compute Corr(w)
3. Convert to level energies via Krawtchouk

Step 2 is the bottleneck. But we can use a histogram approach:
- Group tournaments by H value
- For each pair of H-classes, compute the distance distribution

This is still O(2^{42}) in the worst case. Let me try a sampling approach
or use the Fast Walsh-Hadamard Transform (FWHT).

FWHT approach:
1. Build H[T] for all T (array of length 2^21)
2. Apply FWHT: Ĥ[S] = Σ_T (-1)^{S·T} H[T]
3. E_k = (1/2^{2m}) Σ_{|S|=k} Ĥ[S]²

Step 2: FWHT is O(m * 2^m) = 21 * 2^21 ≈ 44M ops. FAST!
Step 3: Need to sum Ĥ[S]² grouped by |S|. Just popcount.

This is totally feasible! Let's do it.
"""

import numpy as np
import time
import multiprocessing as mp
from fractions import Fraction
from math import comb
from collections import Counter


def compute_H_dp(adj_bits, n):
    """Compute H using Hamiltonian path DP."""
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


def compute_chunk(args):
    """Process a range of tournament indices."""
    start, end, n = args
    edges = []
    for u in range(n):
        for v in range(u+1, n):
            edges.append((u, v))

    results = []
    for idx in range(start, end):
        adj_bits = [0] * n
        for e, (u, v) in enumerate(edges):
            if idx & (1 << e):
                adj_bits[u] |= (1 << v)
            else:
                adj_bits[v] |= (1 << u)
        h = compute_H_dp(adj_bits, n)
        results.append(h)
    return results


def fwht_inplace(a):
    """Fast Walsh-Hadamard Transform (in-place, no normalization)."""
    n = len(a)
    h = 1
    while h < n:
        for i in range(0, n, h * 2):
            for j in range(i, i + h):
                x = a[j]
                y = a[j + h]
                a[j] = x + y
                a[j + h] = x - y
        h *= 2


def main():
    n = 7
    m = n * (n - 1) // 2  # 21
    total = 1 << m  # 2097152

    print("="*70)
    print(f"FOURIER LEVEL DECOMPOSITION AT n={n}")
    print(f"opus-2026-03-14-S89")
    print("="*70)
    print(f"\n  n={n}, m={m}, total={total}")

    # Step 1: Compute H for all tournaments
    print(f"\n  Step 1: Computing H for all {total} tournaments...")
    t0 = time.time()

    num_workers = mp.cpu_count()
    chunk_size = total // num_workers
    chunks = []
    for i in range(num_workers):
        start = i * chunk_size
        end = (i + 1) * chunk_size if i < num_workers - 1 else total
        chunks.append((start, end, n))

    with mp.Pool(num_workers) as pool:
        chunk_results = pool.map(compute_chunk, chunks)

    # Flatten
    H_array = []
    for cr in chunk_results:
        H_array.extend(cr)

    elapsed = time.time() - t0
    print(f"    Done in {elapsed:.1f}s")

    # Convert to numpy for FWHT
    H = np.array(H_array, dtype=np.float64)
    mean_H = np.mean(H)
    print(f"    Mean(H) = {mean_H}")
    print(f"    sum(H) = {np.sum(H)}")

    # Step 2: Apply FWHT
    print(f"\n  Step 2: Fast Walsh-Hadamard Transform...")
    t0 = time.time()
    H_hat = H.copy()
    fwht_inplace(H_hat)
    elapsed = time.time() - t0
    print(f"    Done in {elapsed:.1f}s")

    # H_hat[S] = Σ_T (-1)^{S·T} H(T) (unnormalized)
    # Normalized: ĉ_S = H_hat[S] / 2^m
    # Energy: E_k = (1/2^{2m}) Σ_{|S|=k} H_hat[S]²

    # Step 3: Group by level and compute energy
    print(f"\n  Step 3: Computing level energies...")
    t0 = time.time()

    level_energy = np.zeros(m + 1)
    for S in range(total):
        k = bin(S).count('1')
        level_energy[k] += H_hat[S] ** 2

    # Normalize: E_k = level_energy[k] / 2^{2m}
    # E_0 = H_hat[0]² / 2^{2m} = (Σ H(T))² / 2^{2m} = Mean² * 2^{2m} / 2^{2m} = Mean²
    # Wait: H_hat[0] = Σ H(T), so H_hat[0]²/2^{2m} = (Σ H)² / (2^m)² = Mean² ✓

    E = level_energy / (total ** 2)
    elapsed = time.time() - t0
    print(f"    Done in {elapsed:.1f}s")

    E_0 = E[0]
    print(f"\n    E_0 = Mean² = {E_0:.6f} (expected {mean_H**2:.6f})")

    print(f"\n  FOURIER LEVEL ENERGIES:")
    var_sum = 0
    for k in range(m + 1):
        if E[k] > 1e-10:
            ratio = E[k] / E_0
            num_subsets = comb(m, k)
            avg_coeff_sq = E[k] / num_subsets if num_subsets > 0 else 0
            print(f"    Level {k:2d}: E_{k} = {E[k]:15.6f}, "
                  f"E_{k}/E_0 = {ratio:14.10f}, "
                  f"C({m},{k}) = {num_subsets:8d}, "
                  f"avg ĉ² = {avg_coeff_sq:14.10f}")
            if k > 0:
                var_sum += E[k]

    print(f"\n    Total Var = {var_sum:.6f}")
    print(f"    Var/Mean² = {var_sum/E_0:.10f}")
    print(f"    Expected: {float(Fraction(131, 504)):.10f}")

    # Exact fractions for key ratios
    print(f"\n  EXACT LEVEL RATIOS (as fractions):")
    # E_k/E_0 = level_energy[k] / level_energy[0]
    # Both are sums of integers squared, so rational
    for k in range(m + 1):
        if level_energy[k] > 0.5:
            # level_energy[k] is a sum of H_hat[S]² which are integers (H_hat has integer entries before normalization)
            # So level_energy[k] is an integer
            le_int = int(round(level_energy[k]))
            le0_int = int(round(level_energy[0]))
            ratio = Fraction(le_int, le0_int)
            if k > 0 and ratio > 0:
                print(f"    E_{k}/E_0 = {ratio} = {float(ratio):.10f}")

    # Sum of all E_k/E_0 for k>=1
    total_ratio = Fraction(0)
    for k in range(1, m+1):
        le_int = int(round(level_energy[k]))
        le0_int = int(round(level_energy[0]))
        total_ratio += Fraction(le_int, le0_int)
    print(f"\n    Total Var/Mean² = {total_ratio} = {float(total_ratio):.10f}")
    print(f"    Expected: 131/504 = {float(Fraction(131, 504)):.10f}")
    print(f"    Match: {total_ratio == Fraction(131, 504)}")

    # Analyze level structure
    print(f"\n" + "="*70)
    print(f"PART 2: LEVEL STRUCTURE ANALYSIS")
    print("="*70)

    # At n=5: only levels 2 and 4
    # At n=7: which levels?
    active_levels = [k for k in range(1, m+1) if level_energy[k] > 0.5]
    print(f"\n  Active Fourier levels (nonzero energy): {active_levels}")
    print(f"  Number of active levels: {len(active_levels)}")

    # Only even levels?
    print(f"  Only even levels? {all(k % 2 == 0 for k in active_levels)}")

    # Level energies as fractions of total variance
    print(f"\n  Level contributions to Var/Mean²:")
    for k in active_levels:
        le_int = int(round(level_energy[k]))
        le0_int = int(round(level_energy[0]))
        ratio = Fraction(le_int, le0_int)
        frac_of_var = ratio / total_ratio if total_ratio > 0 else 0
        print(f"    Level {k:2d}: {ratio} = {float(ratio):.10f} ({float(frac_of_var)*100:.2f}% of Var)")

    # Check: do all nonzero ĉ_S have the same magnitude at each level?
    print(f"\n" + "="*70)
    print(f"PART 3: COEFFICIENT MAGNITUDE ANALYSIS")
    print("="*70)

    for k in active_levels[:6]:  # first few levels
        # Find all S with |S|=k, compute |ĉ_S|
        magnitudes = Counter()
        count_nonzero = 0
        for S in range(total):
            if bin(S).count('1') == k:
                coeff = H_hat[S] / total  # normalized
                if abs(coeff) > 1e-12:
                    count_nonzero += 1
                    mag = abs(coeff)
                    # Round to avoid floating point noise
                    mag_rounded = round(mag * 1024) / 1024
                    magnitudes[mag_rounded] += 1

        num_subsets = comb(m, k)
        print(f"\n  Level {k}: {count_nonzero}/{num_subsets} nonzero coefficients")
        for mag, cnt in sorted(magnitudes.items()):
            print(f"    |ĉ| = {mag:.6f}: {cnt} coefficients")

    # Specific analysis of ĉ at level where 63 might be "blocked"
    print(f"\n" + "="*70)
    print(f"PART 4: THE H=63 FOURIER OBSTRUCTION")
    print("="*70)

    # 63 = 2^6 - 1 = 111111₂
    # In the Fourier expansion, H = Σ_S ĉ_S χ_S
    # where χ_S(T) = (-1)^{S·T}
    # H(T) = Σ_S ĉ_S (-1)^{S·T}

    # For H to equal 63, we need a specific combination of ±ĉ_S values.
    # 63 = 2^6 - 1 in binary. The OCF says H is always odd (since H ≡ 1 mod 2).
    # The mean is 315/4 = 78.75. So 63 is BELOW the mean.

    # Distribution of H values near 63
    H_counter = Counter(H_array)
    print(f"\n  H values near 63:")
    for h in range(55, 75, 2):
        count = H_counter.get(h, 0)
        bar = "#" * (count // 2000)
        print(f"    H={h:3d}: {count:8d} {bar}")

    # Is there a parity obstruction at the bit level?
    # In binary: 63 = 0b111111 (low 6 bits all 1)
    # H mod 64 distribution
    print(f"\n  H mod 64 distribution for H values divisible by... ")
    print(f"  H ≡ 63 mod 64 (the 63-family):")
    h63_mod64 = [h for h in H_counter if h % 64 == 63]
    total_63 = sum(H_counter[h] for h in h63_mod64)
    print(f"    Values: {sorted(h63_mod64)}")
    print(f"    Total tournaments: {total_63}")

    # H ≡ 63 mod 128
    h63_mod128 = [h for h in H_counter if h % 128 == 63]
    total_63_128 = sum(H_counter[h] for h in h63_mod128)
    print(f"  H ≡ 63 mod 128: {sorted(h63_mod128)}, total={total_63_128}")

    # H mod 7 analysis
    print(f"\n  H mod 7 for forbidden values:")
    print(f"    7 ≡ {7 % 7} mod 7")
    print(f"    21 ≡ {21 % 7} mod 7")
    print(f"    63 ≡ {63 % 7} mod 7")
    print(f"    All ≡ 0 mod 7")

    # How many H values are ≡ 0 mod 7?
    h_mod7_0 = [h for h in sorted(H_counter) if h % 7 == 0]
    print(f"\n  Achievable H values ≡ 0 mod 7: {h_mod7_0}")
    print(f"  Missing H values ≡ 0 mod 7 in range: ", end="")
    all_odd = set(range(1, max(H_counter)+1, 2))
    achievable = set(H_counter.keys())
    missing = sorted(all_odd - achievable)
    missing_mod7_0 = [h for h in missing if h % 7 == 0]
    print(f"{missing_mod7_0}")

    # Deeper: H mod 9 analysis
    print(f"\n  H mod 9 for forbidden values:")
    print(f"    7 ≡ {7 % 9} mod 9")
    print(f"    21 ≡ {21 % 9} mod 9")
    print(f"    63 ≡ {63 % 9} mod 9")

    # H mod 21
    print(f"\n  H mod 21 for forbidden values:")
    print(f"    7 ≡ {7 % 21} mod 21")
    print(f"    21 ≡ {21 % 21} mod 21")
    print(f"    63 ≡ {63 % 21} mod 21")

    print(f"\n" + "="*70)
    print(f"PART 5: COMPARISON WITH n=5 AND n=6")
    print("="*70)

    # n=5 Fourier levels
    print(f"\n  Level structure comparison:")
    print(f"    n=5 (m=10): levels 2, 4 only")
    print(f"    n=7 (m=21): levels {active_levels}")
    print(f"    ")
    print(f"    n=5: even levels only? True")
    print(f"    n=7: even levels only? {all(k%2==0 for k in active_levels)}")

    # The spectrum size sequence: 2, 3, 7, 19, 77
    print(f"\n  H-spectrum sizes: 2, 3, 7, 19, 77")
    print(f"    Differences: 1, 4, 12, 58")
    print(f"    Ratios: 1.5, 2.33, 2.71, 4.05")
    print(f"    77/19 = {Fraction(77, 19)} ≈ 4.053")
    print(f"    19/7 = {Fraction(19, 7)} ≈ 2.714")
    print(f"    7/3 = {Fraction(7, 3)} ≈ 2.333")

    # Is 2, 3, 7, 19, 77 related to anything known?
    # 2, 3, 7, 19 are primes (except 77 = 7×11)
    # Check OEIS-like patterns
    print(f"    2×3 = 6, 3×7 = 21 (!), 7×19 = 133, 19×77 = 1463")
    print(f"    Product of consecutive: 6, 21, 133, 1463")
    print(f"    21 = forbidden! 133 = appears in H-spectrum at n=7!")

    # a(n+1) = a(n)*a(n-1)/a(n-2) + correction?
    # 3*2/? = 6, 7*3/2 = 10.5, not 19
    # Different: 2, 3, 7, 19, 77 → check: a(n) = 5*a(n-1) - 5*a(n-2) + a(n-3)?
    # 5*7 - 5*3 + 2 = 35-15+2 = 22 ≠ 19
    # a(n) = 3*a(n-1) - a(n-2)?
    # 3*7 - 3 = 18 ≠ 19. Close!
    # 3*19 - 7 = 50 ≠ 77.
    # a(n) = 4*a(n-1) + a(n-2) - 4*a(n-3)?
    # 4*7 + 3 - 4*2 = 28+3-8 = 23 ≠ 19.
    # Hmm.

    # Let me check: are these spectrum sizes p_n - n + 1 where p_n = nth prime?
    # 2=2, 3=3, 7=7, 19=19 are primes. 77 = 7*11.
    # Primes: 2, 3, 5, 7, 11, 13, 17, 19, 23, ...
    # These aren't consecutive primes.

    # n=3→2, n=4→3, n=5→7, n=6→19, n=7→77
    # 2, 3, 7, 19, 77
    # OEIS search hint: this is close to A001339 or similar
    # 2·1+1=3, 3·2+1=7, 7·3-2=19, 19·4+1=77? Check: 19*4+1=77 ✓!
    # 3·2+1=7 ✓, 7·3-2=19 ✓, 19·4+1=77 ✓
    # Pattern: a(n) = a(n-1) * (n-2) + (-1)^n ?
    # a(3)=2, a(4)=3=2*1+1, a(5)=7=3*2+1, a(6)=19=7*3-2, a(7)=77=19*4+1
    # Corrections: +1, +1, -2, +1
    # Not clean. But a(n) ≈ a(n-1)*(n-2) is close.
    # Actually: 2*2-1=3, 3*3-2=7, 7*3-2=19, 19*4+1=77... nope.

    print(f"\n  Trying to find recurrence for spectrum sizes 2, 3, 7, 19, 77:")
    sizes = [2, 3, 7, 19, 77]
    for a in range(-5, 6):
        for b in range(-5, 6):
            ok = True
            for i in range(2, len(sizes)):
                pred = a * sizes[i-1] + b * sizes[i-2]
                if pred != sizes[i]:
                    ok = False
                    break
            if ok and (a != 0 or b != 0):
                print(f"    a(n) = {a}*a(n-1) + {b}*a(n-2) works!")

    # Try a(n) = c1*a(n-1)*a(n-2) + c2*a(n-1) + c3*a(n-2) + c4
    # This is getting complicated. Let me just note the sequence.

    print(f"\n" + "="*70)
    print(f"DONE — FOURIER LEVEL DECOMPOSITION n=7")
    print("="*70)


if __name__ == "__main__":
    main()
