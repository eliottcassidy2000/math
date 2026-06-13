#!/usr/bin/env python3
"""
OEIS sequence exploration — session opus-2026-05-28-S3.

Covers:
  1. SC/non-SC tiling sequences to n=22+
  2. Good-cuts distribution triangle to n=16
  3. King-count distribution triangle for n<=7
  4. 3-cycle count distribution for n<=7 (plus formula for larger n)
  5. H-value distribution for n<=6
  6. Sequences from OCF/independence polynomial
  7. General d-good formula verification (THM-339)

Uses bitmask IE for SC (fast), brute force for small n, formulas for larger n.
"""

import itertools
import time

# ============================================================
# Utility
# ============================================================

def get_tiles(n):
    """Tiles (x,y) in n-vertex tiling model: x >= y+2, x < n."""
    return [(x, y) for y in range(n-1) for x in range(y+2, n)]

def binomial(n, k):
    if k < 0 or k > n: return 0
    if k == 0 or k == n: return 1
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    return r

# ============================================================
# 1. SC/Non-SC tiling sequences via bitmask IE  (fast, n up to ~22)
# ============================================================

def sc_and_nonsc(n):
    """Compute (sc, nonsc) tiling counts using inclusion-exclusion over cut subsets."""
    if n <= 1: return (1, 0)
    tiles = get_tiles(n)
    m = len(tiles)
    num_cuts = n - 1

    # Precompute: for each tile, which cuts (as bitmask) does it cross?
    tile_cut_mask = []
    for x, y in tiles:
        mask = 0
        for k in range(y+1, x+1):
            mask |= (1 << (k-1))
        tile_cut_mask.append(mask)

    # Group tiles by their cut-mask for speed
    from collections import Counter
    mask_count = Counter(tile_cut_mask)

    # For subset S (bitmask over cuts): f(S) = #{tiles crossing at least one cut in S}
    # = #{tiles whose cut_mask & S != 0} = m - #{tiles whose cut_mask & S == 0}
    nonzero_mask = list(mask_count.keys())
    counts_by_mask = [mask_count[k] for k in nonzero_mask]

    sc = 0
    for S in range(1 << num_cuts):
        # f(S) = # tiles with cut_mask & S != 0
        fS = sum(c for mk, c in zip(nonzero_mask, counts_by_mask) if mk & S)
        sign = (-1) ** bin(S).count('1')
        sc += sign * (1 << (m - fS))

    return sc, (1 << m) - sc

# ============================================================
# 2. Good-cuts distribution triangle via Möbius transform
# ============================================================

def good_cuts_distribution(n):
    """
    Return dist[j] = #{tilings with exactly j good cuts} for j=0..n-1.
    Uses Möbius inversion on the powerset of cuts.
    """
    if n == 2:
        return [1, 1]
    tiles = get_tiles(n)
    m = len(tiles)
    num_cuts = n - 1
    num_subsets = 1 << num_cuts

    # Precompute cut mask for each tile
    tile_cut_mask = []
    for x, y in tiles:
        mask = 0
        for k in range(y+1, x+1):
            mask |= (1 << (k-1))
        tile_cut_mask.append(mask)

    # g[S] = #{tilings with all cuts in S bad} = 2^{m - f(S)}
    from collections import Counter
    mc = Counter(tile_cut_mask)
    nonzero_masks = list(mc.keys())
    nonzero_counts = [mc[k] for k in nonzero_masks]

    g = [0] * num_subsets
    for S in range(num_subsets):
        fS = sum(c for mk, c in zip(nonzero_masks, nonzero_counts) if mk & S)
        g[S] = 1 << (m - fS)

    # Möbius transform: h[B] = #{tilings with exact bad-set B}
    # h[B] = sum_{S superset of B} (-1)^{|S\B|} g[S]
    # Computed by: for each bit, h[S] -= h[S | (1<<bit)] for S without bit.
    h = g[:]
    for bit in range(num_cuts):
        for S in range(num_subsets):
            if not (S & (1 << bit)):
                h[S] -= h[S | (1 << bit)]

    # Accumulate by good-cut count
    dist = [0] * n   # dist[j] for j=0..n-1 (n-1 is max good cuts)
    for B in range(num_subsets):
        bad_count = bin(B).count('1')
        good_count = num_cuts - bad_count
        dist[good_count] += h[B]

    return dist

# ============================================================
# 3. King-count distribution (brute force, n <= 7)
# ============================================================

def king_distribution(n):
    """
    K[k] = #{labeled n-vertex tournaments with exactly k kings}.
    Brute force over all 2^{C(n,2)} labeled tournaments.
    """
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    all_mask = (1 << n) - 1
    king_dist = [0] * (n+1)

    for T in range(1 << m):
        out = [0] * n
        for k, (i,j) in enumerate(pairs):
            if T & (1 << k):
                out[i] |= (1 << j)
            else:
                out[j] |= (1 << i)
        # Count kings: v is king iff reach2[v] = all vertices except v
        kings = 0
        for v in range(n):
            reach2 = out[v]
            tmp = out[v]
            while tmp:
                u = (tmp & -tmp).bit_length() - 1
                reach2 |= out[u]
                tmp &= tmp - 1
            if (reach2 | (1<<v)) == all_mask:
                kings += 1
        king_dist[kings] += 1

    return king_dist

# ============================================================
# 4. 3-cycle distribution (brute force n<=7, formula for larger)
# ============================================================

def three_cycle_distribution(n):
    """
    C3[k] = #{labeled n-vertex tournaments with exactly k directed 3-cycles}.
    Uses: #3-cycles = C(n,3) - sum_v C(out_degree(v), 2).
    """
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    cn3 = n*(n-1)*(n-2)//6

    # Max possible 3-cycles
    c3_dist = {}

    for T in range(1 << m):
        out_deg = [0] * n
        for k, (i,j) in enumerate(pairs):
            if T & (1 << k):
                out_deg[i] += 1
            else:
                out_deg[j] += 1
        c3 = cn3 - sum(d*(d-1)//2 for d in out_deg)
        c3_dist[c3] = c3_dist.get(c3, 0) + 1

    max_k = max(c3_dist)
    return [c3_dist.get(k, 0) for k in range(max_k+1)]

# ============================================================
# 5. H-value distribution (brute force DP, n<=7 — use n<=6 for speed)
# ============================================================

def hamiltonian_path_count(out, n):
    """Count Hamiltonian paths using DP with bitmask states."""
    dp = [0] * (n * (1 << n))  # dp[v * 2^n + mask] = paths ending at v using vertices in mask
    for v in range(n):
        dp[v * (1<<n) + (1<<v)] = 1
    for mask in range(1, 1 << n):
        pc = bin(mask).count('1')
        if pc < 2: continue
        for v in range(n):
            if not (mask & (1<<v)): continue
            prev = mask ^ (1<<v)
            val = 0
            # Sum over predecessors u with u->v
            tmp_mask = prev
            while tmp_mask:
                u = (tmp_mask & -tmp_mask).bit_length() - 1
                if out[u] & (1<<v):
                    val += dp[u * (1<<n) + prev]
                tmp_mask &= tmp_mask - 1
            dp[v * (1<<n) + mask] = val
    full = (1 << n) - 1
    return sum(dp[v * (1<<n) + full] for v in range(n))

def h_distribution(n):
    """H_dist[h] = #{labeled n-vertex tournaments with H=h}."""
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    h_dist = {}

    for T in range(1 << m):
        out = [0] * n
        for k, (i,j) in enumerate(pairs):
            if T & (1 << k):
                out[i] |= (1 << j)
            else:
                out[j] |= (1 << i)
        h = hamiltonian_path_count(out, n)
        h_dist[h] = h_dist.get(h, 0) + 1

    return h_dist

# ============================================================
# 6. Sequences from OCF/independence polynomial
# ============================================================

def odd_cycle_sequences(n_max):
    """
    Total directed odd cycles of each length across all labeled n-vertex tournaments.
    Formula: total directed k-cycles = 2^{C(n,2)} * C(n,k) * (k-1)! / 2^k
    For odd k only.
    """
    print("\n=== Total directed k-cycles across ALL labeled n-vertex tournaments ===")
    print("Formula: 2^{C(n,2)} * C(n,k) * (k-1)! / 2^k")
    print()

    # For each n, compute sum over odd k of total directed k-cycles
    # This equals 2^{C(n,2)} * E[alpha_1(T)]
    for n in range(3, n_max+1):
        cn2 = n*(n-1)//2
        total_two_power = 1 << cn2
        print(f"n={n}: 2^{cn2} = {total_two_power}")

        total_alpha1 = 0
        for k in range(3, n+1, 2):  # odd cycles of length 3,5,7,...
            cnk = 1
            for i in range(k): cnk = cnk * (n-i) // (i+1)
            fact_k1 = 1
            for i in range(1, k): fact_k1 *= i  # (k-1)!
            total_k_cycles = total_two_power * cnk * fact_k1 // (1 << k)
            print(f"  k={k}: total directed {k}-cycles = {total_k_cycles}")
            total_alpha1 += cnk * fact_k1

        # E[H] = n! / 2^{n-1}
        nfact = 1
        for i in range(1, n+1): nfact *= i
        e_h_num = nfact
        e_h_den = 1 << (n-1)
        print(f"  E[H] = {e_h_num}/{e_h_den} = {e_h_num/e_h_den:.6f}")
        print(f"  Σ_T H(T) = {total_two_power * e_h_num // e_h_den}")

# ============================================================
# 7. General d-good formula verification (THM-339)
# ============================================================

def q_values(d_max, sc_seq):
    """
    Compute Q(d,k) = sum over ordered k-tuples (b1,...,bk) with sum=d, bi>=2
                     of product SC(bi+1).
    sc_seq: sc_seq[n] = SC(n) for n=2,3,4,...
    Returns Q[d][k] for d=2..d_max, k=1..floor(d/2).
    """
    Q = {}
    for d in range(2, d_max+1):
        Q[d] = {}
        for k in range(1, d//2 + 1):
            # Sum over ordered k-tuples (b1,...,bk) with sum=d, each bi>=2
            total = 0
            # Use dynamic programming
            # dp[j][s] = sum of products for j-tuples summing to s
            dp_prev = {0: 1}  # 0-tuple summing to 0
            for _ in range(k):
                dp_new = {}
                for s_prev, val in dp_prev.items():
                    for b in range(2, d+1):
                        s_new = s_prev + b
                        if s_new <= d and b+1 < len(sc_seq) and sc_seq[b+1] is not None:
                            sc_val = sc_seq[b+1]
                            dp_new[s_new] = dp_new.get(s_new, 0) + val * sc_val
                dp_prev = dp_new
            Q[d][k] = dp_prev.get(d, 0)
    return Q

def verify_d_good_formula(n_max, sc_seq, good_cuts_data):
    """
    Verify THM-339: exactly-d-good(n) = sum_{k=1}^{floor(d/2)} C(n-d,k) * Q(d,k).
    """
    print("\n=== THM-339 Verification: General d-good formula ===")
    d_max = min(n_max - 1, 10)
    Q = q_values(d_max, sc_seq)

    print("Q(d,k) values (d=2..8, k=1..floor(d/2)):")
    for d in range(2, 9):
        row = [Q.get(d, {}).get(k, 0) for k in range(1, d//2+1)]
        print(f"  d={d}: {row}")

    print()
    all_ok = True
    for n in range(3, n_max+1):
        if n not in good_cuts_data: continue
        dist = good_cuts_data[n]
        for d in range(2, n):
            if d >= len(dist): continue
            actual = dist[d]
            predicted = sum(binomial(n-d, k) * Q.get(d, {}).get(k, 0)
                           for k in range(1, d//2+1))
            ok = actual == predicted
            if not ok:
                print(f"  MISMATCH n={n}, d={d}: actual={actual}, predicted={predicted}")
                all_ok = False
    if all_ok:
        print(f"  ALL MATCH for n=3..{min(n_max,16)}, d=2..n-1 ✓")

# ============================================================
# Main computation
# ============================================================

def main():
    print("=" * 70)
    print("OEIS SEQUENCE EXPLORATION — opus-2026-05-28-S3")
    print("=" * 70)

    # -------------------------------------------------------
    # Part 1: SC and non-SC tiling sequences
    # -------------------------------------------------------
    print("\n=== PART 1: SC and non-SC tiling sequences ===")
    print("(Exact counts via bitmask inclusion-exclusion)")
    print()

    sc_seq = {1: 1, 2: 1}
    nonsc_seq = {}
    t0 = time.time()

    for n in range(2, 24):
        t1 = time.time()
        sc, nonsc = sc_and_nonsc(n)
        elapsed = time.time() - t1
        sc_seq[n] = sc
        nonsc_seq[n] = nonsc
        m = (n-1)*(n-2)//2
        total = 1 << m
        print(f"n={n:2d}: SC={sc:>25}, nonSC={nonsc:>25}, "
              f"total=2^{m:3d}={total:>25}  [{elapsed:.2f}s]")
        if elapsed > 30:
            print(f"  (stopping at n={n} due to time)")
            break

    print(f"\nTotal time: {time.time()-t0:.1f}s")

    print("\n--- SC tiling sequence (n=2..N) ---")
    sc_vals = [sc_seq[n] for n in sorted(sc_seq)]
    print(", ".join(str(x) for x in sc_vals))

    print("\n--- Non-SC tiling sequence (n=2..N) ---")
    nonsc_vals = [nonsc_seq[n] for n in sorted(nonsc_seq)]
    print(", ".join(str(x) for x in nonsc_vals))

    print("\n--- Ratio non-SC / 2^{m-n+3} (should -> 1) ---")
    for n in sorted(nonsc_seq):
        m = (n-1)*(n-2)//2
        if m >= n-3:
            denom = 1 << (m - n + 3)
            ratio = nonsc_seq[n] / denom if denom > 0 else float('inf')
            print(f"  n={n}: {nonsc_seq[n]} / {denom} = {ratio:.8f}")

    print("\n--- log2(SC(n)) and m = C(n-1,2) ---")
    import math
    for n in sorted(sc_seq):
        if sc_seq[n] > 0:
            m = (n-1)*(n-2)//2
            log2sc = math.log2(sc_seq[n]) if sc_seq[n] > 0 else None
            print(f"  n={n}: m={m}, log2(SC)={log2sc:.4f}, m-log2(SC)={m-log2sc:.4f}")

    # -------------------------------------------------------
    # Part 2: Good-cuts distribution triangle
    # -------------------------------------------------------
    print("\n\n=== PART 2: Good-cuts distribution triangle T(n,k) ===")
    print("T(n,k) = #{tilings with exactly k good cuts}, k=0..n-1")
    print()

    good_cuts_data = {}
    for n in range(2, 17):
        t1 = time.time()
        dist = good_cuts_distribution(n)
        elapsed = time.time() - t1
        good_cuts_data[n] = dist
        print(f"n={n:2d}: {dist}  sum={sum(dist)}  [2^{(n-1)*(n-2)//2}={(1<<((n-1)*(n-2)//2))}]  [{elapsed:.2f}s]")
        if elapsed > 20:
            print(f"  Stopping at n={n} for time.")
            break

    # Print triangle for OEIS
    print("\n--- Triangle T(n,k) by rows ---")
    for n in sorted(good_cuts_data):
        dist = good_cuts_data[n]
        print(", ".join(str(x) for x in dist))

    # Diagonal sequences
    print("\n--- Column sequences k=0,1,2,3 ---")
    for k in range(5):
        col = [good_cuts_data[n][k] for n in sorted(good_cuts_data) if k < len(good_cuts_data[n])]
        print(f"k={k}: {col}")

    # -------------------------------------------------------
    # Part 3: SC sequence analysis — general d formula
    # -------------------------------------------------------
    sc_list = {n: sc_seq.get(n, None) for n in range(2, 20)}

    print("\n--- Q(d,k) and THM-339 formula check ---")
    verify_d_good_formula(16, sc_list, good_cuts_data)

    # Print exact d-good formulas
    print("\n--- Exact d-good formulas via THM-339 ---")
    Q = q_values(10, sc_list)
    for d in range(2, 11):
        terms = []
        for k in range(1, d//2+1):
            qdk = Q.get(d, {}).get(k, 0)
            if qdk:
                terms.append(f"C(n-{d},{k})*{qdk}")
        formula = " + ".join(terms) if terms else "0"
        print(f"  exactly-{d}-good(n) = {formula}")

    # -------------------------------------------------------
    # Part 4: King-count distribution
    # -------------------------------------------------------
    print("\n\n=== PART 4: King-count distribution K(n,k) ===")
    print("K(n,k) = #{labeled n-vertex tournaments with exactly k kings}")
    print()

    king_data = {}
    for n in range(2, 8):
        t1 = time.time()
        kdist = king_distribution(n)
        elapsed = time.time() - t1
        king_data[n] = kdist
        nonzero = [(k, v) for k,v in enumerate(kdist) if v]
        print(f"n={n}: {nonzero}  total={sum(kdist)}  [{elapsed:.1f}s]")
        if elapsed > 60:
            print(f"  (stopping king dist at n={n})")
            break

    print("\n--- King-count triangle by rows ---")
    for n in sorted(king_data):
        kdist = king_data[n]
        print(f"n={n}: {[kdist[k] for k in range(1,n+1)]}")

    # Derived sequences
    print("\n--- #{tournaments with exactly 1 king} (transitive-dominant) ---")
    one_king = [king_data[n][1] for n in sorted(king_data) if 1 < len(king_data[n])]
    print(one_king)

    print("--- #{tournaments with all n kings} (universal kings) ---")
    all_kings = [king_data[n][n] for n in sorted(king_data) if n < len(king_data[n])]
    print(all_kings)

    # -------------------------------------------------------
    # Part 5: 3-cycle count distribution
    # -------------------------------------------------------
    print("\n\n=== PART 5: 3-cycle count distribution ===")
    print("C3(n,k) = #{labeled n-vertex tournaments with exactly k directed 3-cycles}")
    print()

    c3_data = {}
    for n in range(3, 8):
        t1 = time.time()
        c3dist = three_cycle_distribution(n)
        elapsed = time.time() - t1
        c3_data[n] = c3dist
        print(f"n={n}: {c3dist}  total={sum(c3dist)}  [{elapsed:.1f}s]")

    print("\n--- 3-cycle distribution by rows ---")
    for n in sorted(c3_data):
        print(f"n={n}: {c3_data[n]}")

    # Formula verification
    print("\n--- Formula: total 3-cycles = 2^{C(n,2)-2} * C(n,3) ---")
    for n in range(3, 8):
        cn2 = n*(n-1)//2
        cn3 = n*(n-1)*(n-2)//6
        formula = (1 << (cn2-2)) * cn3
        if n in c3_data:
            actual = sum(k * c3_data[n][k] for k in range(len(c3_data[n])))
            print(f"  n={n}: formula={formula}, actual={actual}, match={formula==actual}")

    # -------------------------------------------------------
    # Part 6: H-value distribution for n<=6
    # -------------------------------------------------------
    print("\n\n=== PART 6: H-value distribution ===")
    print("H_dist(n) = {H: count} for labeled n-vertex tournaments")
    print()

    h_data = {}
    for n in range(2, 7):
        t1 = time.time()
        hd = h_distribution(n)
        elapsed = time.time() - t1
        h_data[n] = hd
        print(f"n={n}: {sorted(hd.items())}  total={sum(hd.values())}  [{elapsed:.1f}s]")

    # Summary sequences
    print("\n--- Max H by n (= A038375) ---")
    max_h = [max(h_data[n].keys()) for n in sorted(h_data)]
    print(max_h)

    print("\n--- Min H by n ---")
    min_h = [min(h_data[n].keys()) for n in sorted(h_data)]
    print(min_h)

    print("\n--- #{H-values achieved} by n ---")
    n_h = [len(h_data[n]) for n in sorted(h_data)]
    print(n_h)

    print("\n--- Sum of all H values (= 2^{C(n,2)} * n!/2^{n-1}) ---")
    for n in sorted(h_data):
        total = sum(h * cnt for h, cnt in h_data[n].items())
        cn2 = n*(n-1)//2
        nfact = 1
        for i in range(1, n+1): nfact *= i
        formula = (1 << cn2) * nfact // (1 << (n-1))
        print(f"  n={n}: sum={total}, formula={formula}, match={total==formula}")

    # -------------------------------------------------------
    # Part 7: OCF / odd cycle sequences
    # -------------------------------------------------------
    print("\n\n=== PART 7: Sequences from OCF formula ===")

    # Total directed k-cycles across all labeled n-vertex tournaments
    print("Total directed k-cycles across ALL labeled n-vertex tournaments:")
    print("Formula: 2^{C(n,2)} * C(n,k) * (k-1)! / 2^k")
    for k in [3, 5, 7]:
        vals = []
        for n in range(k, 10):
            cn2 = n*(n-1)//2
            cnk = binomial(n, k)
            fact_k1 = 1
            for i in range(1, k): fact_k1 *= i
            val = (1 << cn2) * cnk * fact_k1 // (1 << k)
            vals.append(val)
        print(f"  k={k}: {vals}")

    # Sequences related to independence polynomial
    print("\n--- #{SC tilings with exactly j upward tiles} at n=5 ---")
    # SC tilings on n=5 have m=6 tiles. Distribution by #upward tiles.
    n = 5
    tiles_5 = get_tiles(n)
    m5 = len(tiles_5)
    # For each SC tiling, count upward tiles
    tile_cut_masks = []
    for x, y in tiles_5:
        mask = 0
        for k in range(y+1, x+1):
            mask |= (1 << (k-1))
        tile_cut_masks.append(mask)

    sc_by_upward = [0] * (m5+1)
    num_cuts_5 = n - 1  # 4 cuts
    all_cuts = (1 << num_cuts_5) - 1
    for bits in range(1 << m5):
        # Check if all cuts are good (SC condition)
        covered = 0
        for i, tcm in enumerate(tile_cut_masks):
            if bits & (1 << i):
                covered |= tcm
        if covered == all_cuts:
            sc_by_upward[bin(bits).count('1')] += 1
    print(f"  n=5: {sc_by_upward}  total={sum(sc_by_upward)}  (expected SC(5)=50)")

    print("\n--- #{SC tilings with exactly j upward tiles} at n=6 ---")
    n = 6
    tiles_6 = get_tiles(n)
    m6 = len(tiles_6)
    tile_cut_masks_6 = []
    for x, y in tiles_6:
        mask = 0
        for k in range(y+1, x+1):
            mask |= (1 << (k-1))
        tile_cut_masks_6.append(mask)
    all_cuts_6 = (1 << (n-1)) - 1
    sc_by_upward_6 = [0] * (m6+1)
    for bits in range(1 << m6):
        covered = 0
        for i, tcm in enumerate(tile_cut_masks_6):
            if bits & (1 << i):
                covered |= tcm
        if covered == all_cuts_6:
            sc_by_upward_6[bin(bits).count('1')] += 1
    print(f"  n=6: {sc_by_upward_6}  total={sum(sc_by_upward_6)}  (expected SC(6)=903)")

    # -------------------------------------------------------
    # Part 8: Summary of new sequences for OEIS
    # -------------------------------------------------------
    print("\n\n=== PART 8: Summary of new sequences for OEIS ===")

    print("\n1. SC tiling sequence (path-fixed SC tournaments):")
    sc_for_oeis = [sc_seq[n] for n in range(3, max(sc_seq.keys())+1) if n in sc_seq]
    print("   " + ", ".join(str(x) for x in sc_for_oeis[:15]))

    print("\n2. Non-SC tiling sequence:")
    nonsc_for_oeis = [nonsc_seq[n] for n in range(3, max(nonsc_seq.keys())+1) if n in nonsc_seq]
    print("   " + ", ".join(str(x) for x in nonsc_for_oeis[:15]))

    print("\n3. Good-cuts triangle T(n,k), reading by rows:")
    triangle_vals = []
    for n in range(2, max(good_cuts_data.keys())+1):
        if n in good_cuts_data:
            triangle_vals.extend(good_cuts_data[n])
    print("   " + ", ".join(str(x) for x in triangle_vals[:50]))

    print("\n4. T(n,n-1) = SC tiling counts (matches seq 1):")
    sc_diag = [good_cuts_data[n][-1] for n in sorted(good_cuts_data) if n in good_cuts_data]
    print("   " + ", ".join(str(x) for x in sc_diag[:12]))

    print("\n5. T(n,2) = n-2 sequence: (confirmed formula n-2)")
    t_n_2 = [good_cuts_data[n][2] for n in sorted(good_cuts_data) if n in good_cuts_data and 2 < len(good_cuts_data[n])]
    print("   " + ", ".join(str(x) for x in t_n_2[:12]))

    print("\n6. King-count triangle K(n,k) for k=1..n, by rows:")
    for n in sorted(king_data):
        row = [king_data[n][k] for k in range(1, n+1)]
        print(f"   n={n}: {row}")

    print("\n7. #{tournaments with exactly 1 king} (= #{transitive-dominant tournaments}):")
    ok1 = [king_data[n][1] for n in sorted(king_data) if 1 < len(king_data[n])]
    print("   " + ", ".join(str(x) for x in ok1))

    print("\n8. #{tournaments with all n kings}:")
    okn = [king_data[n][n] for n in sorted(king_data) if n < len(king_data[n])]
    print("   " + ", ".join(str(x) for x in okn))

    print("\n9. H-value distribution triangle:")
    for n in sorted(h_data):
        h_vals = sorted(h_data[n].items())
        print(f"   n={n}: H-values={[h for h,c in h_vals]}")
        print(f"         counts  ={[c for h,c in h_vals]}")

    print("\n10. #{distinct H-values at n} sequence:")
    nh = [len(h_data[n]) for n in sorted(h_data)]
    print("    " + ", ".join(str(x) for x in nh))

    print("\n11. SC tiling by #upward tiles:")
    print(f"   n=5: {sc_by_upward}")
    print(f"   n=6: {sc_by_upward_6}")

    # -------------------------------------------------------
    # Part 9: Extensions using pure formulas (large n)
    # -------------------------------------------------------
    print("\n\n=== PART 9: Large-n formula extensions ===")

    Q = q_values(12, sc_list)
    print("Exact d-good counts for n=3..20 via THM-339:")
    print(f"{'n':>4}", end="")
    for d in range(2, 9):
        print(f"  d={d}:{'':>8}", end="")
    print()

    for n in range(3, 21):
        print(f"{n:4d}", end="")
        for d in range(2, 9):
            val = sum(binomial(n-d, k) * Q.get(d, {}).get(k, 0)
                     for k in range(1, d//2+1)) if n > d else 0
            print(f"  {val:>12}", end="")
        print()

    print("\nDone.")

if __name__ == "__main__":
    main()
