#!/usr/bin/env python3
"""
Extended computations for opus-2026-05-28-S3:
- SC upward-tile distribution for n=2..10 (fast: only 2^m tilings with m=C(n-1,2))
- Good-cuts triangle extended to n=20
- Q(d,k) table extended
- Non-SC upward distribution
- King triangle extended analysis
- Prove K(n,1) = n*2^{C(n-1,2)} analytically
"""

import itertools
import time

def binomial(n, k):
    if k < 0 or k > n: return 0
    if k == 0 or k == n: return 1
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    return r

def get_tiles(n):
    return [(x, y) for y in range(n-1) for x in range(y+2, n)]

# ============================================================
# 1. SC and non-SC distributions by #upward tiles
# ============================================================

def sc_upward_distribution(n):
    """
    sc_dist[j] = #{SC tilings with exactly j upward tiles}
    nonsc_dist[j] = #{non-SC tilings with exactly j upward tiles}
    Uses brute force over 2^m tilings.
    """
    tiles = get_tiles(n)
    m = len(tiles)
    num_cuts = n - 1
    all_cuts = (1 << num_cuts) - 1

    # Precompute cut mask for each tile
    tile_cut_mask = []
    for x, y in tiles:
        mask = 0
        for k in range(y+1, x+1):
            mask |= (1 << (k-1))
        tile_cut_mask.append(mask)

    sc_dist = [0] * (m + 1)
    nonsc_dist = [0] * (m + 1)

    for bits in range(1 << m):
        j = bin(bits).count('1')
        # Check SC: all cuts covered
        covered = 0
        for i, tcm in enumerate(tile_cut_mask):
            if bits & (1 << i):
                covered |= tcm
        if covered == all_cuts:
            sc_dist[j] += 1
        else:
            nonsc_dist[j] += 1

    return sc_dist, nonsc_dist

def sc_nonsc_ie(n):
    """Fast IE formula for SC/non-SC."""
    if n <= 1: return 1, 0
    tiles = get_tiles(n)
    m = len(tiles)
    num_cuts = n - 1
    from collections import Counter
    tile_cut_mask = []
    for x, y in tiles:
        mask = 0
        for k in range(y+1, x+1):
            mask |= (1 << (k-1))
        tile_cut_mask.append(mask)
    mc = Counter(tile_cut_mask)
    nonzero_masks = list(mc.keys())
    nonzero_counts = [mc[k] for k in nonzero_masks]
    sc = 0
    for S in range(1 << num_cuts):
        fS = sum(c for mk, c in zip(nonzero_masks, nonzero_counts) if mk & S)
        sign = (-1) ** bin(S).count('1')
        sc += sign * (1 << (m - fS))
    return sc, (1 << m) - sc

# ============================================================
# 2. Good-cuts triangle extended to n=20
# ============================================================

def good_cuts_extended(n):
    """Compute good-cuts distribution for n vertices using Möbius transform."""
    if n == 2: return [1, 1]
    tiles = get_tiles(n)
    m = len(tiles)
    num_cuts = n - 1
    num_subsets = 1 << num_cuts

    from collections import Counter
    tile_cut_mask = []
    for x, y in tiles:
        mask = 0
        for k in range(y+1, x+1):
            mask |= (1 << (k-1))
        tile_cut_mask.append(mask)

    mc = Counter(tile_cut_mask)
    nzm = list(mc.keys())
    nzc = [mc[k] for k in nzm]

    g = [0] * num_subsets
    for S in range(num_subsets):
        fS = sum(c for mk, c in zip(nzm, nzc) if mk & S)
        g[S] = 1 << (m - fS)

    h = g[:]
    for bit in range(num_cuts):
        for S in range(num_subsets):
            if not (S & (1 << bit)):
                h[S] -= h[S | (1 << bit)]

    dist = [0] * n
    for B in range(num_subsets):
        bad_count = bin(B).count('1')
        good_count = num_cuts - bad_count
        dist[good_count] += h[B]
    return dist

# ============================================================
# 3. Q(d,k) table: convolution coefficients [x^d]A(x)^k
# ============================================================

def compute_Q_table(d_max, sc_vals):
    """
    Q[d][k] = sum over ordered k-tuples (b1,...,bk) with sum=d, bi>=2
              of product SC(bi+1).
    Returns as dict.
    """
    Q = {}
    for d in range(2, d_max + 1):
        Q[d] = {}
        for k in range(1, d // 2 + 1):
            dp = {0: 1}
            for _ in range(k):
                ndp = {}
                for s, val in dp.items():
                    for b in range(2, d - s + 1):
                        if b + 1 in sc_vals and sc_vals[b+1] is not None:
                            ns = s + b
                            ndp[ns] = ndp.get(ns, 0) + val * sc_vals[b+1]
                dp = ndp
            Q[d][k] = dp.get(d, 0)
    return Q

# ============================================================
# Main
# ============================================================

def main():
    print("=" * 70)
    print("EXTENDED COMPUTATIONS — opus-2026-05-28-S3")
    print("=" * 70)

    # -------------------------------------------------------
    # Part 1: SC upward-tile distribution
    # -------------------------------------------------------
    print("\n=== PART 1: SC and non-SC tilings by #upward tiles ===")
    print("sc_dist[j] = #{SC tilings with j upward tiles}")
    print()

    sc_upward_table = {}
    nonsc_upward_table = {}

    for n in range(2, 11):
        t1 = time.time()
        m = (n-1)*(n-2)//2
        if m > 22:  # 2^22 = 4M is feasible, 2^24 is getting slow
            print(f"n={n}: skipping (m={m} too large for brute force)")
            continue
        sc_d, nonsc_d = sc_upward_distribution(n)
        elapsed = time.time() - t1
        sc_upward_table[n] = sc_d
        nonsc_upward_table[n] = nonsc_d
        print(f"n={n} (m={m}): SC  ={sc_d}")
        print(f"       nonSC={nonsc_d}  [{elapsed:.2f}s]")
        print(f"       sums: SC={sum(sc_d)}, nonSC={sum(nonsc_d)}, total={sum(sc_d)+sum(nonsc_d)}=2^{m}={1<<m}")
        print()
        if elapsed > 30:
            print(f"  (stopping at n={n} for time)")
            break

    # Print SC upward triangle
    print("SC upward-tile distribution by rows (n=2..?):")
    for n in sorted(sc_upward_table):
        row = sc_upward_table[n]
        nonzero = [x for x in row if x > 0]
        print(f"  n={n}: {row}  (nonzero: {nonzero})")

    print("\nNon-SC upward-tile distribution by rows:")
    for n in sorted(nonsc_upward_table):
        row = nonsc_upward_table[n]
        nonzero = [x for x in row if x > 0]
        print(f"  n={n}: {row}  (nonzero: {nonzero})")

    # Relation to binomial: SC_j = C(m,j) - nonSC_j
    print("\n-- Check: SC_j + nonSC_j = C(m,j) --")
    for n in sorted(sc_upward_table):
        m = (n-1)*(n-2)//2
        sc_d = sc_upward_table[n]
        nonsc_d = nonsc_upward_table[n]
        ok = all(sc_d[j] + nonsc_d[j] == binomial(m, j) for j in range(m+1))
        print(f"  n={n}: {'OK' if ok else 'FAIL'}")

    # -------------------------------------------------------
    # Part 2: Good-cuts triangle to n=22
    # -------------------------------------------------------
    print("\n\n=== PART 2: Good-cuts triangle extended ===")

    print("\nExtending to n=17..22:")
    good_cuts_ext = {}
    for n in range(17, 23):
        t1 = time.time()
        dist = good_cuts_extended(n)
        elapsed = time.time() - t1
        good_cuts_ext[n] = dist
        print(f"n={n}: [...]  sc={dist[-1]}, total=2^{(n-1)*(n-2)//2}  [{elapsed:.1f}s]")
        if elapsed > 60:
            print(f"  (stopping at n={n})")
            break

    print("\nColumn sequences from extended table:")
    all_n = list(range(3, 23))
    all_good_dists = {}
    for n in range(3, 17):
        d = good_cuts_extended(n) if n not in good_cuts_ext else good_cuts_ext[n]
        all_good_dists[n] = d

    print("Column k=5 (d=5):")
    col5 = [all_good_dists[n][5] if 5 < len(all_good_dists.get(n,[])) else None
            for n in range(6, 17)]
    print(f"  {col5}")

    print("Column k=6 (d=6):")
    col6 = [all_good_dists[n][6] if 6 < len(all_good_dists.get(n,[])) else None
            for n in range(7, 17)]
    print(f"  {col6}")

    print("Column k=7 (d=7):")
    col7 = [all_good_dists[n][7] if 7 < len(all_good_dists.get(n,[])) else None
            for n in range(8, 17)]
    print(f"  {col7}")

    # -------------------------------------------------------
    # Part 3: Q(d,k) extended table
    # -------------------------------------------------------
    print("\n\n=== PART 3: Q(d,k) extended table ===")
    # Need SC values
    print("Computing SC values via IE for n=2..18...")
    sc_vals = {}
    for n in range(2, 19):
        sc, _ = sc_nonsc_ie(n)
        sc_vals[n] = sc

    Q = compute_Q_table(18, sc_vals)

    print("\nQ(d,k) table:")
    for d in range(2, 16):
        row = [Q.get(d,{}).get(k,0) for k in range(1, d//2+1)]
        print(f"  d={d:2d}: {row}")

    # Print Q triangle as OEIS flat sequence
    print("\n--- Q(d,k) triangle, by rows d=2..12, k=1..floor(d/2) ---")
    for d in range(2, 13):
        row = [Q.get(d,{}).get(k,0) for k in range(1, d//2+1)]
        print(", ".join(str(x) for x in row))

    # -------------------------------------------------------
    # Part 4: King distribution analysis — prove K(n,1) = n*2^{C(n-1,2)}
    # -------------------------------------------------------
    print("\n\n=== PART 4: K(n,1) formula verification ===")
    print("Claim: K(n,1) = n * 2^{C(n-1,2)}")
    print("Theorem: a tournament has exactly 1 king iff it has an Emperor (score n-1 vertex).")
    print("Proof: Emperor => unique king (Emperor beats all; others can't reach Emperor in 2 steps).")
    print("Converse: unique king v must have score n-1 (otherwise some vertex w beats v,")
    print("  and w can reach all that v can reach, so w is also a king — contradiction).")
    print()
    king_n1 = [2, 6, 32, 320, 6144, 229376]  # From previous computation
    for i, n in enumerate(range(2, 8)):
        predicted = n * (1 << ((n-1)*(n-2)//2))
        actual = king_n1[i]
        print(f"  n={n}: predicted={predicted}, actual={actual}, match={predicted==actual}")

    # -------------------------------------------------------
    # Part 5: New sequences from SC upward distribution
    # -------------------------------------------------------
    print("\n\n=== PART 5: New sequences from SC upward distribution ===")

    # The max j with SC_j > 0 for each n
    print("Max j with SC tilings:")
    for n in sorted(sc_upward_table):
        m = (n-1)*(n-2)//2
        sc_d = sc_upward_table[n]
        max_j = max(j for j,v in enumerate(sc_d) if v > 0)
        min_j = min(j for j,v in enumerate(sc_d) if v > 0) if any(sc_d) else 0
        print(f"  n={n} (m={m}): min upward={min_j}, max upward={max_j}")

    # The minimum number of upward tiles in any SC tiling
    print("\nMinimum upward tiles in SC tiling:")
    print("(Should be 1 for all n>=3, since tile (n-1,0) covers all cuts)")
    for n in sorted(sc_upward_table):
        sc_d = sc_upward_table[n]
        for j,v in enumerate(sc_d):
            if v > 0:
                print(f"  n={n}: min={j}")
                break

    # SC tilings with MINIMUM upward tiles (count for each n)
    print("\n#{SC tilings with minimum (=1) upward tiles}:")
    for n in sorted(sc_upward_table):
        sc_d = sc_upward_table[n]
        print(f"  n={n}: SC_1 = {sc_d[1]}")

    # The maximum upward tiles in SC tilings
    print("\n#{SC tilings with MAXIMUM upward tiles (=m)}:")
    for n in sorted(sc_upward_table):
        m = (n-1)*(n-2)//2
        sc_d = sc_upward_table[n]
        print(f"  n={n}: SC_m = SC_{m} = {sc_d[m]}")

    # -------------------------------------------------------
    # Part 6: The "d-good" anti-diagonal sequences
    # -------------------------------------------------------
    print("\n\n=== PART 6: Anti-diagonal sequences in T(n,k) ===")
    print("T(n, n-1) = SC(n) (already known)")
    print("T(n, n-2) = ?")
    print("T(n, n-3) = ?")

    for back in [1, 2, 3, 4, 5]:
        seq = []
        for n in range(2+back, 17):
            if n in all_good_dists and (n-back) < len(all_good_dists[n]):
                seq.append(all_good_dists[n][n-back])
        print(f"T(n, n-{back}): {seq}")

    # -------------------------------------------------------
    # Part 7: Non-SC asymptotics precision
    # -------------------------------------------------------
    print("\n\n=== PART 7: Non-SC asymptotics fine structure ===")
    print("non-SC(n) = 2^{m-n+3} - R(n) where R(n) is the correction.")
    for n in range(3, 20):
        sc, nonsc = sc_nonsc_ie(n)
        m = (n-1)*(n-2)//2
        dominant = 1 << (m - n + 3) if m >= n-3 else 0
        R = dominant - nonsc  # R = how much dominant term overestimates
        print(f"  n={n}: nonSC={nonsc}, 2^{{m-n+3}}={dominant}, R={R}")

    # -------------------------------------------------------
    # Part 8: Summary
    # -------------------------------------------------------
    print("\n\n=== PART 8: OEIS Candidate Sequences Summary ===")

    print("\n=== NEW OEIS sequences (not in OEIS as of search date) ===\n")

    print("SEQ1: SC tiling sequence (path-fixed SC tournaments)")
    sc_for_oeis = []
    for n in range(3, 19):
        sc, _ = sc_nonsc_ie(n)
        sc_for_oeis.append(sc)
    print("  " + ", ".join(str(x) for x in sc_for_oeis))

    print("\nSEQ2: Non-SC tiling sequence")
    nonsc_for_oeis = []
    for n in range(3, 19):
        _, nonsc = sc_nonsc_ie(n)
        nonsc_for_oeis.append(nonsc)
    print("  " + ", ".join(str(x) for x in nonsc_for_oeis))

    print("\nSEQ3: King distribution triangle K(n,k), n>=2, k=1..n")
    print("  Known values:")
    print("  n=2: 2, 0")
    print("  n=3: 6, 0, 2")
    print("  n=4: 32, 0, 32, 0")
    print("  n=5: 320, 0, 520, 120, 64")
    print("  n=6: 6144, 0, 11600, 7920, 5424, 1680")
    print("  n=7: 229376, 0, 402640, 527520, 491568, 336000, 110048")

    print("\nSEQ4: K(n,n) = #{tournaments with all n kings}")
    print("  0, 2, 0, 64, 1680, 110048")
    print("  (n=2..7; 0 for n=2,4)")

    print("\nSEQ5: K(n,3) = #{tournaments with exactly 3 kings}")
    print("  2, 32, 520, 11600, 402640")
    print("  (n=3..7)")

    print("\nSEQ6: Good-cuts triangle T(n,k)")
    print("  Rows n=2..6:")
    print("  1, 1")
    print("  1, 0, 1")
    print("  1, 0, 2, 5")
    print("  1, 0, 3, 10, 50")
    print("  1, 0, 4, 15, 101, 903")

    print("\nSEQ7: Q(d,k) triangle (SC convolution coefficients)")
    print("  Rows d=2..8:")
    for d in range(2, 9):
        row = [Q.get(d,{}).get(k,0) for k in range(1, d//2+1)]
        print(f"  d={d}: " + ", ".join(str(x) for x in row))

    print("\nSEQ8: SC tilings by #upward tiles (triangle)")
    print("  Row n=5 (j=0..6): 0, 1, 9, 18, 15, 6, 1")
    print("  Row n=6 (j=0..10): 0, 1, 17, 81, 180, 240, 208, 120, 45, 10, 1")
    for n in sorted(sc_upward_table):
        if n >= 7:
            print(f"  Row n={n} (j=0..{(n-1)*(n-2)//2}): {sc_upward_table[n]}")

    print("\nSEQ9: #{distinct achievable H values at n}")
    print("  1, 2, 3, 7, 19  (n=2..6; next term unknown)")

    print("\n=== KNOWN sequences confirmed ===")
    print("K(n,1) = A123903(n) = n * 2^{C(n-1,2)}  [Maurer 1980 King Chicken Theorem]")
    print("K(n,2) = 0 for all n  [No tournament has exactly 2 kings]")

    print("\nDone.")

if __name__ == "__main__":
    main()
