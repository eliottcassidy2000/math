"""
sequences_exact.py — opus-2026-04-04-S5
Compute exact sequences for n=3..7, sample for n=8.
"""
import numpy as np
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
            if mask == full:
                total += cnt; continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + cnt
    return total

def mobius_coeffs(H, m):
    coeffs = {}
    for S in range(1 << m):
        val = 0
        T = S
        while True:
            sign = (-1) ** bin(S ^ T).count('1')
            val += sign * int(H[T])
            if T == 0: break
            T = (T - 1) & S
        if val != 0:
            coeffs[S] = val
    return coeffs

def mean_H_linear(n):
    return 1 + sum((n - s) * (2 ** (s - 2)) for s in range(2, n))

def main():
    print("EXACT SEQUENCE COMPUTATION n=3..7")
    print("=" * 60)

    seq_n = []
    seq_m = []
    seq_P = []        # nonzero multilinear coeffs
    seq_sumH = []      # sum of all H
    seq_meanH = []     # mean H
    seq_maxH = []      # max H
    seq_distinct = []  # distinct H values
    seq_Elin = []      # linear mean approximation

    for n in range(3, 8):
        tiles = make_tiles(n)
        m = len(tiles)
        t0 = time.time()

        # Compute all H
        H = np.zeros(1 << m, dtype=np.int64)
        for bits in range(1 << m):
            T = tiling_to_tournament(n, tiles, bits)
            H[bits] = count_hp(T, n)

        sumH = int(np.sum(H))
        meanH = sumH / (1 << m)
        maxH = int(np.max(H))
        distinct = len(set(H))

        # Multilinear coefficients
        coeffs = mobius_coeffs(H, m)
        P = len(coeffs)

        # Max degree
        d_max = max(bin(S).count('1') for S in coeffs) if coeffs else 0

        elapsed = time.time() - t0

        seq_n.append(n)
        seq_m.append(m)
        seq_P.append(P)
        seq_sumH.append(sumH)
        seq_meanH.append(meanH)
        seq_maxH.append(maxH)
        seq_distinct.append(distinct)
        seq_Elin.append(mean_H_linear(n))

        print(f"\nn={n}: m={m}, 2^m={1<<m}, time={elapsed:.1f}s")
        print(f"  Σ H = {sumH}")
        print(f"  E[H] = {meanH:.4f}")
        print(f"  max H = {maxH}")
        print(f"  Distinct H = {distinct}")
        print(f"  P(n) = {P}")
        print(f"  d_max = {d_max}")
        print(f"  E_lin = {seq_Elin[-1]}")

    print(f"\n{'='*60}")
    print(f"SEQUENCES")
    print(f"{'='*60}")

    print(f"\nn:         {seq_n}")
    print(f"m:         {seq_m}")

    print(f"\nP(n) = nonzero multilinear coefficients:")
    print(f"  {seq_P}")

    print(f"\nΣ H = sum over all 2^m tilings:")
    print(f"  {seq_sumH}")

    print(f"\n2^m · E[H] = Σ H (integer!):")
    print(f"  {seq_sumH}")

    print(f"\nE[H] = mean Hamiltonian paths with fixed base path:")
    print(f"  {seq_meanH}")

    print(f"\nmax H (= A038375 for global max):")
    print(f"  {seq_maxH}")

    print(f"\nDistinct H values:")
    print(f"  {seq_distinct}")

    print(f"\nE_lin(n) = 1 + Σ (n-s)·2^(s-2):")
    for n in range(3, 15):
        print(f"  E_lin({n}) = {mean_H_linear(n)}")

    # The sequence Σ H(n) = total HP count over all tilings at n
    # This counts: for EACH of the 2^m tiling configurations, the number of HPs.
    # Equivalently: number of (tiling, HP) pairs.
    # A (tiling, HP) pair is: a choice of m binary tile values + a Hamiltonian path
    # consistent with those tiles (and the base path).

    print(f"\nΣ H counts (tiling, HP) pairs = |{{(t, σ) : σ is HP of T(t)}}|")
    print(f"  This equals: Σ_σ 2^(m - |tiles constrained by σ|)")
    print(f"  = Σ_σ 2^(m - (n-1) + (base path arcs used by σ))")
    print(f"  = Σ_σ 2^(C(n-1,2) - n + 1 + bp(σ))")
    print(f"  where bp(σ) = number of base-path arcs used by σ.")

    # This is a beautiful identity! Let me verify.
    for n in range(3, 8):
        tiles = make_tiles(n)
        m = len(tiles)
        # Each HP σ uses n-1 arcs. Some are base-path (k+1→k), rest are tiles.
        # The tiles used by σ are constrained (each in a specific direction).
        # The tiles NOT used by σ are free (2 choices each).
        # So the number of tilings compatible with σ is 2^(m - tile_arcs(σ))
        # = 2^(m - (n-1-bp(σ))) = 2^(m - n + 1 + bp(σ))

        # Sum over all σ: Σ_σ 2^(m-n+1+bp(σ))

        # Compute from the actual data
        H = np.zeros(1 << m, dtype=np.int64)
        for bits in range(1 << m):
            T = tiling_to_tournament(n, tiles, bits)
            H[bits] = count_hp(T, n)
        actual_sum = int(np.sum(H))

        # Compute from the formula: need to know the distribution of bp(σ) over all HPs
        # Actually, we can compute: for EACH tiling, what is bp(σ) for each HP?
        # This is complex. Let me just verify the formula interpretation.

        # Total (tiling, HP) pairs = Σ H(t) = actual_sum
        # Expected from formula: Σ_σ 2^(m-n+1+bp(σ)) where σ ranges over ALL HPs
        # of ALL tournaments (with repetition: σ appears in multiple tilings)

        # Actually σ is an HP of the SPECIFIC tournament T(t). Each (t, σ) pair is counted
        # once. The formula says the number of tilings containing σ is 2^(m - tile_arcs(σ)).

        # This gives: Σ H(t) = Σ_{σ as labeled HP} 2^(m - tile_arcs(σ))
        # where tile_arcs(σ) = n-1-bp(σ) (arcs of σ that are tiles, not base path)

        # bp(σ) ranges from 0 (no base path arcs used) to n-1 (the base path itself)
        # The base path IS an HP, with bp = n-1, tile_arcs = 0. Contributes 2^m.
        # Other HPs use fewer base-path arcs.

        print(f"  n={n}: Σ H = {actual_sum} = {actual_sum // (1<<m)} × 2^{m} + {actual_sum % (1<<m)}")

    # Let me compute the base-path-arc distribution
    print(f"\n{'='*60}")
    print(f"BASE-PATH ARC DISTRIBUTION IN HPs")
    print(f"{'='*60}")

    for n in range(3, 7):
        tiles = make_tiles(n)
        m = len(tiles)

        # For EACH tiling, find all HPs and count base-path arcs
        bp_dist = defaultdict(int)  # bp_count → total occurrences across all tilings

        for bits in range(1 << m):
            T = tiling_to_tournament(n, tiles, bits)

            # Find all HPs and their bp counts
            # DP tracking bp count
            dp = {}  # (mask, v) → list of (count, bp_arcs)
            for v in range(n):
                dp[(1 << v, v)] = [(1, 0)]  # one path, 0 bp arcs

            full = (1 << n) - 1
            for mask in range(1, 1 << n):
                for v in range(n):
                    entries = dp.get((mask, v), [])
                    if not entries: continue
                    for count, bp in entries:
                        if mask == full:
                            bp_dist[bp] += count
                            continue
                        for u in range(n):
                            if mask & (1 << u): continue
                            if T[v][u]:
                                is_bp = (v == u + 1)
                                key = (mask | (1 << u), u)
                                if key not in dp:
                                    dp[key] = []
                                dp[key].append((count, bp + (1 if is_bp else 0)))

            # Merge dp entries (this is getting complex, let me simplify)
            # Actually the above double-counts. Let me use a cleaner approach.

        # Simpler: just count total HP × 2^(free tiles) for each bp value
        # Total (tiling, HP) pairs with bp base-path arcs in the HP:
        # = (number of distinct HP paths using bp base arcs) × 2^(m - (n-1-bp))
        # = N(bp) × 2^(m - n + 1 + bp)

        # I can compute N(bp) = number of distinct labeled HPs with exactly bp base-path arcs
        # by enumerating all permutations

        from itertools import permutations
        bp_count = defaultdict(int)
        for perm in permutations(range(n)):
            bp = 0
            valid = True
            for i in range(n-1):
                u, v = perm[i], perm[i+1]
                if u == v + 1:
                    bp += 1
                # Note: we're counting labeled HPs, not checking tournament arcs
            bp_count[bp] += 1

        print(f"\n  n={n}: HP distribution by base-path arcs:")
        formula_sum = 0
        for bp in sorted(bp_count.keys()):
            N_bp = bp_count[bp]
            contrib = N_bp * (2 ** (m - n + 1 + bp))
            formula_sum += contrib
            print(f"    bp={bp}: N={N_bp} HPs, contribution = {N_bp} × 2^{m-n+1+bp} = {contrib}")

        print(f"    Formula total: {formula_sum}")

        # Wait, this doesn't account for direction! A labeled HP (σ₁,...,σ_n) uses
        # arc σ_i → σ_{i+1}. This arc might be a base-path arc (σ_i = σ_{i+1}+1)
        # or a tile arc. For a tile arc (a,b) with a>b+1, the HP requires the tile
        # to be in the correct direction. For each tile used by the HP, there's only
        # 1 valid setting (not 2). So the number of tilings compatible with this HP
        # is 2^(tiles NOT used by HP) = 2^(m - (n-1-bp)).

        # But the formula_sum above doesn't match Σ H because:
        # Not all labeled HPs are valid in all tournaments!
        # A HP uses specific arc directions. Some arc pairs (i→j and j→i) are
        # exclusive. The formula works because:
        # Σ_t Σ_{σ∈HP(T(t))} 1 = Σ_σ |{t : σ∈HP(T(t))}| = Σ_σ 2^(m-tile_arcs(σ))

        # And σ ranges over ALL labeled permutations that COULD be HPs (i.e., all n! perms,
        # since for each perm, there exists a unique tiling making it an HP).

        # Wait, that's not right either. A permutation σ = (v₁,...,v_n) is an HP iff
        # T[v_i][v_{i+1}] = 1 for all i. The tiling t determines T, and for σ to be an HP,
        # certain tiles must be in specific directions. If two arcs of σ require the SAME
        # tile in different directions, the HP is impossible in ANY tiling. But each arc
        # of σ uses a DISTINCT tile (or base-path arc), so no conflicts! Each labeled
        # permutation is an HP of exactly 2^(m - tile_arcs(σ)) tilings.

        # So: Σ H = Σ_σ 2^(m - (n-1-bp(σ))) = Σ_{bp=0}^{n-1} N(bp) · 2^(m-n+1+bp)

        # Verify:
        H_arr = np.zeros(1 << m, dtype=np.int64)
        for bits in range(1 << m):
            T2 = tiling_to_tournament(n, tiles, bits)
            H_arr[bits] = count_hp(T2, n)
        actual = int(np.sum(H_arr))
        print(f"    Actual Σ H: {actual}")
        print(f"    Match: {'✓' if formula_sum == actual else '✗'}")

if __name__ == "__main__":
    main()
