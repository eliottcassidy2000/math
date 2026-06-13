#!/usr/bin/env python3
"""
Focused analysis of W(i/2) vs H(T) at n=7.

KEY FORMULAS (derived from coefficient analysis):
  For n=7: W(i/2) = (2*N_5 - N_3) / 8
  where N_f = #{permutations with exactly f forward edges in T}
  and N_f = N_{6-f} (path reversal symmetry).

  H(T) = N_6, and the full generating function is sum_f N_f * x^f.

This script computes (H, W(i/2)) for all 2^21 = 2,097,152 tournaments at n=7.
Uses bitmask DP to compute N_3 and N_5 without full permutation enumeration.

The DP computes dp[mask][v][f_mod] where f_mod tracks forward-edge count
modulo what's needed (we only need N_1, N_3, N_5 for the formula).
"""

from collections import defaultdict
import sys
import time

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                T[i][j] = 1
            else:
                T[j][i] = 1
            k += 1
    return T


def compute_H_and_Nf_dp(T):
    """
    Bitmask DP to compute N_f = #{permutations with f forward edges}.
    Returns (H, [N_0, ..., N_6]).
    O(2^n * n^2 * n) = O(2^7 * 49 * 7) = ~44K operations. Very fast.
    """
    n = len(T)
    deg = n - 1  # 6
    full = (1 << n) - 1

    # dp[mask][v] = list of counts by f (forward-edge count)
    # f ranges 0..6. We store as a flat array.
    # Total entries: 128 * 7 = 896, each with 7 ints.

    # For memory efficiency, use a 3D array
    dp = [[[0] * n for _ in range(n)] for _ in range(1 << n)]

    for v in range(n):
        dp[1 << v][v][0] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            fv = dp[mask][v]
            for u in range(n):
                if mask & (1 << u):
                    continue
                new_mask = mask | (1 << u)
                if T[v][u]:  # forward edge
                    for f in range(deg):
                        if fv[f]:
                            dp[new_mask][u][f + 1] += fv[f]
                else:  # backward
                    for f in range(n):
                        if fv[f]:
                            dp[new_mask][u][f] += fv[f]

    Nf = [0] * n
    for v in range(n):
        for f in range(n):
            Nf[f] += dp[full][v][f]

    H = Nf[deg]
    return H, Nf


def count_3cycles(T):
    """Count directed 3-cycles."""
    n = len(T)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    count += 1
                elif T[i][k] and T[k][j] and T[j][i]:
                    count += 1
    return count


def main():
    n = 7
    m = n * (n - 1) // 2  # 21
    total = 1 << m  # 2,097,152

    print(f"W(i/2) vs H(T) analysis for all {total} tournaments on n={n}")
    print(f"Formula: W(i/2) = (2*N_5 - N_3) / 8")
    print("=" * 70)

    H_all = []
    W_all = []  # W(i/2) * 8 (integer)
    N3_all = []
    N5_all = []
    c3_all = []

    t0 = time.time()
    step = total // 20

    for bits in range(total):
        if bits % step == 0:
            elapsed = time.time() - t0
            print(f"  Progress: {bits}/{total} ({100*bits/total:.0f}%) [{elapsed:.1f}s]", file=sys.stderr)

        T = tournament_from_bits(n, bits)
        H, Nf = compute_H_and_Nf_dp(T)
        W8 = 2 * Nf[5] - Nf[3]  # W(i/2) * 8
        c3 = count_3cycles(T)

        H_all.append(H)
        W_all.append(W8)
        N3_all.append(Nf[3])
        N5_all.append(Nf[5])
        c3_all.append(c3)

    elapsed = time.time() - t0
    print(f"\nCompleted in {elapsed:.1f}s", file=sys.stderr)

    N = total

    # =====================================================================
    # BASIC STATISTICS
    # =====================================================================
    print(f"\n{'='*70}")
    print("BASIC STATISTICS")
    print(f"{'='*70}")

    # W(i/2) = W8 / 8. Check divisibility
    div8 = sum(1 for w in W_all if w % 8 == 0)
    print(f"  W*8 divisible by 8: {div8}/{N} ({100*div8/N:.1f}%)")
    print(f"  W*8 always divisible by 2: {all(w % 2 == 0 for w in W_all)}")
    print(f"  W*8 always divisible by 4: {all(w % 4 == 0 for w in W_all)}")

    # W(i/2) values
    W_vals = sorted(set(W_all))
    print(f"\n  Distinct W*8 values: {len(W_vals)}")
    print(f"  W*8 range: [{min(W_all)}, {max(W_all)}]")

    # Is W*8 always even?
    parity = set(w % 2 for w in W_all)
    print(f"  W*8 parity: {parity}")

    # W(i/2) sign distribution
    pos = sum(1 for w in W_all if w > 0)
    neg = sum(1 for w in W_all if w < 0)
    zero = sum(1 for w in W_all if w == 0)
    print(f"\n  W(i/2) > 0: {pos} ({100*pos/N:.1f}%)")
    print(f"  W(i/2) = 0: {zero} ({100*zero/N:.1f}%)")
    print(f"  W(i/2) < 0: {neg} ({100*neg/N:.1f}%)")

    # =====================================================================
    # CORRELATION ANALYSIS
    # =====================================================================
    print(f"\n{'='*70}")
    print("CORRELATIONS")
    print(f"{'='*70}")

    mean_H = sum(H_all) / N
    mean_W = sum(W_all) / N
    mean_absW = sum(abs(w) for w in W_all) / N
    mean_c3 = sum(c3_all) / N

    def corr(xs, ys):
        mx = sum(xs) / len(xs)
        my = sum(ys) / len(ys)
        cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / len(xs)
        vx = sum((x - mx)**2 for x in xs) / len(xs)
        vy = sum((y - my)**2 for y in ys) / len(ys)
        return cov / (vx**0.5 * vy**0.5) if vx > 0 and vy > 0 else 0

    absW_all = [abs(w) for w in W_all]

    print(f"  Corr(H, W(i/2))      = {corr(H_all, W_all):.6f}")
    print(f"  Corr(H, |W(i/2)|)    = {corr(H_all, absW_all):.6f}")
    print(f"  Corr(H, c3)          = {corr(H_all, c3_all):.6f}")
    print(f"  Corr(c3, W(i/2))     = {corr(c3_all, W_all):.6f}")
    print(f"  Corr(c3, |W(i/2)|)   = {corr(c3_all, absW_all):.6f}")
    print(f"  Corr(H, N3)          = {corr(H_all, N3_all):.6f}")
    print(f"  Corr(H, N5)          = {corr(H_all, N5_all):.6f}")
    print(f"  Corr(N3, N5)         = {corr(N3_all, N5_all):.6f}")

    # =====================================================================
    # JOINT DISTRIBUTION TABLE
    # =====================================================================
    print(f"\n{'='*70}")
    print("JOINT DISTRIBUTION: H vs W(i/2)*8")
    print(f"{'='*70}")

    H_to_data = defaultdict(list)
    for i in range(N):
        H_to_data[H_all[i]].append((W_all[i], c3_all[i], N3_all[i], N5_all[i]))

    H_set = sorted(H_to_data.keys())

    print(f"\n  {'H':>5}  {'count':>7}  {'W*8 min':>8}  {'W*8 max':>8}  {'|W|*8 mn':>8}  {'|W|*8 mx':>8}  {'W/8 min':>8}  {'W/8 max':>8}  {'c3 rng':>10}")
    for H in H_set:
        entries = H_to_data[H]
        Ws = [e[0] for e in entries]
        c3s = [e[1] for e in entries]
        mn_W, mx_W = min(Ws), max(Ws)
        mn_absW, mx_absW = min(abs(w) for w in Ws), max(abs(w) for w in Ws)

        print(f"  {H:>5}  {len(entries):>7}  {mn_W:>8}  {mx_W:>8}  {mn_absW:>8}  {mx_absW:>8}  {mn_W/8:>8.1f}  {mx_W/8:>8.1f}  [{min(c3s):>2},{max(c3s):>2}]")

    # =====================================================================
    # MONOTONICITY CHECK
    # =====================================================================
    print(f"\n{'='*70}")
    print("MONOTONICITY CHECK")
    print(f"{'='*70}")

    # For |W(i/2)| vs H: is it anti-monotone?
    absW_ranges = {}
    W_ranges = {}
    for H in H_set:
        entries = H_to_data[H]
        Ws = [e[0] for e in entries]
        absW_ranges[H] = (min(abs(w) for w in Ws), max(abs(w) for w in Ws))
        W_ranges[H] = (min(Ws), max(Ws))

    # Strict: H1 < H2 => min|W|(H1) >= max|W|(H2)
    violations = 0
    for H1 in H_set:
        for H2 in H_set:
            if H1 < H2 and absW_ranges[H2][1] > absW_ranges[H1][0]:
                violations += 1
    print(f"  |W| anti-monotonicity violations (H1<H2, max|W(H2)| > min|W(H1)|): {violations}")

    # Mean |W| by H (is the MEAN monotonically decreasing?)
    print(f"\n  H vs mean |W|*8:")
    prev_mean = float('inf')
    mean_violations = 0
    for H in H_set:
        entries = H_to_data[H]
        Ws = [abs(e[0]) for e in entries]
        mean_w = sum(Ws) / len(Ws)
        marker = " <-- VIOLATION" if mean_w > prev_mean else ""
        if mean_w > prev_mean:
            mean_violations += 1
        print(f"    H={H:>3}: mean|W|*8 = {mean_w:>10.1f}  (n={len(entries):>6}){marker}")
        prev_mean = mean_w
    print(f"\n  Mean |W| monotonicity violations: {mean_violations}")

    # =====================================================================
    # W(i/2) = 0 CHARACTERIZATION
    # =====================================================================
    print(f"\n{'='*70}")
    print("W(i/2) = 0 CHARACTERIZATION")
    print(f"{'='*70}")

    # W*8 = 2*N5 - N3 = 0 means N3 = 2*N5
    zero_entries = [(H_all[i], c3_all[i], N3_all[i], N5_all[i]) for i in range(N) if W_all[i] == 0]
    print(f"  Tournaments with W(i/2)=0: {len(zero_entries)}")

    if zero_entries:
        zero_H = sorted(set(e[0] for e in zero_entries))
        print(f"  H values when W=0: {zero_H}")
        for H in zero_H:
            cnt = sum(1 for e in zero_entries if e[0] == H)
            print(f"    H={H}: {cnt} tournaments")

    # =====================================================================
    # EXACT FORMULA FOR W(i/2)^2
    # =====================================================================
    print(f"\n{'='*70}")
    print("W(i/2)^2 ANALYSIS")
    print(f"{'='*70}")

    # W*8 = 2*N5 - N3
    # (W*8)^2 = 4*N5^2 - 4*N3*N5 + N3^2
    # W^2 = (2*N5 - N3)^2 / 64

    W2_all = [w * w for w in W_all]
    W2_set = sorted(set(W2_all))
    print(f"  Distinct (W*8)^2 values: {len(W2_set)}")

    # Check divisibility of W^2 by 64
    div64 = sum(1 for w2 in W2_all if w2 % 64 == 0)
    print(f"  (W*8)^2 div by 64: {div64}/{N}")

    # W^2 as integer
    W2_int = sorted(set(w2 // 64 for w2 in W2_all if w2 % 64 == 0))
    print(f"  W(i/2)^2 integer values (first 20): {W2_int[:20]}")

    # Check perfect square
    import math
    ps_count = 0
    non_ps = []
    for w2 in W2_int:
        sq = math.isqrt(w2)
        if sq * sq == w2:
            ps_count += 1
        else:
            non_ps.append(w2)
    print(f"  Perfect squares: {ps_count}/{len(W2_int)}")
    if non_ps:
        print(f"  Non-perfect-square W^2 values (first 10): {non_ps[:10]}")

    # =====================================================================
    # RELATIONSHIP TO c3 (3-cycle count)
    # =====================================================================
    print(f"\n{'='*70}")
    print("W(i/2) vs c3 (directed 3-cycle count)")
    print(f"{'='*70}")

    # At n=5, c3 uniquely determines W(i/2). Is this true at n=7?
    c3_to_W = defaultdict(set)
    for i in range(N):
        c3_to_W[c3_all[i]].add(W_all[i])

    c3_set = sorted(c3_to_W.keys())
    unique_determination = all(len(v) == 1 for v in c3_to_W.values())
    print(f"  c3 uniquely determines W: {unique_determination}")

    print(f"\n  {'c3':>5}  {'distinct W*8 values':>40}")
    for c3 in c3_set:
        ws = sorted(c3_to_W[c3])
        if len(ws) <= 8:
            print(f"  {c3:>5}  {ws}")
        else:
            print(f"  {c3:>5}  [{ws[0]}..{ws[-1]}] ({len(ws)} values)")

    # =====================================================================
    # LINEAR FORMULA: W*8 = a + b*c3 + c*H + d*??
    # =====================================================================
    print(f"\n{'='*70}")
    print("SEARCHING FOR LINEAR FORMULA")
    print(f"{'='*70}")

    # At n=5: c3 determines both H and W uniquely. Each c3 value gives one (H, W) pair.
    # At n=7: c3 does NOT uniquely determine H or W.

    # Try: is W*8 a linear function of c3 and H?
    # W*8 = a*H + b*c3 + c
    # Use least squares on first 10000 points
    sub = min(N, 10000)
    from_idx = list(range(0, N, N // sub))[:sub]

    # Simple linear regression: W8 = a*H + b*c3 + c
    sum_H = sum(H_all[i] for i in from_idx)
    sum_c3 = sum(c3_all[i] for i in from_idx)
    sum_W = sum(W_all[i] for i in from_idx)
    sum_HH = sum(H_all[i]**2 for i in from_idx)
    sum_cc = sum(c3_all[i]**2 for i in from_idx)
    sum_Hc = sum(H_all[i]*c3_all[i] for i in from_idx)
    sum_HW = sum(H_all[i]*W_all[i] for i in from_idx)
    sum_cW = sum(c3_all[i]*W_all[i] for i in from_idx)
    nn = len(from_idx)

    # Solve normal equations: [[HH, Hc, H], [Hc, cc, c3], [H, c3, n]] * [a,b,c] = [HW, cW, W]
    # Use numpy-free approach... just report correlations
    print(f"  (Skipping full regression, reporting key relationships)")

    # Check if W*8 = 2*N5 - N3 has a nice relationship with c3
    # Recall: c3 = #{directed 3-cycles}
    # N5 = #{permutations with exactly 5 forward edges}
    # N3 = #{permutations with exactly 3 forward edges}
    # These are related to the score sequence and higher-order invariants.

    # For n=5: W(i/2)*16 = -8*H + 4*N_2, and N_2 is determined by c3.
    # The reason: at n=5, score sequence determines everything, and c3 determines the score.

    # At n=7, the relationship is more complex.

    # Let me check: does (c3, score_sequence) determine W?
    # Computing score sequences would help but is expensive for all 2M tournaments.

    # Instead, let me check: among tournaments with same (H, c3), is W always the same?
    Hc3_to_W = defaultdict(set)
    for i in range(N):
        Hc3_to_W[(H_all[i], c3_all[i])].add(W_all[i])

    multi_W = sum(1 for ws in Hc3_to_W.values() if len(ws) > 1)
    print(f"\n  (H, c3) pairs with multiple W values: {multi_W}/{len(Hc3_to_W)}")

    # Among tournaments with same H, how much does W vary?
    for H in [189, 175, 171, 159, 1, 3]:
        if H in H_to_data:
            entries = H_to_data[H]
            Ws = sorted(set(e[0] for e in entries))
            c3s = sorted(set(e[1] for e in entries))
            print(f"  H={H}: W*8 values = {Ws}, c3 values = {c3s}")


    # =====================================================================
    # FORMULA VERIFICATION: W(i/2) in terms of the polynomial P(u)
    # =====================================================================
    print(f"\n{'='*70}")
    print("POLYNOMIAL STRUCTURE")
    print(f"{'='*70}")

    # Let's verify on a few examples
    from itertools import permutations as perms
    for bits in [0, 1, 12345]:
        T = tournament_from_bits(7, bits)
        H, Nf = compute_H_and_Nf_dp(T)
        W8 = 2*Nf[5] - Nf[3]
        c3 = count_3cycles(T)

        # Verify by explicit permutation enumeration (slow but correct)
        Nf_check = [0] * 7
        for P in perms(range(7)):
            f = sum(1 for k in range(6) if T[P[k]][P[k+1]])
            Nf_check[f] += 1

        match = Nf == Nf_check
        print(f"  bits={bits}: H={H}, Nf={Nf}, check={match}")
        print(f"    W*8 = 2*{Nf[5]} - {Nf[3]} = {W8}")
        print(f"    W(i/2) = {W8/8}")
        print(f"    c3 = {c3}")

        # Also verify the full W(i/2) from complex arithmetic
        W_complex = 0j
        for P in perms(range(7)):
            prod = 1+0j
            for k in range(6):
                s = T[P[k]][P[k+1]] - 0.5
                prod *= (0.5j + s)
            W_complex += prod
        print(f"    W(i/2) direct = {W_complex:.6f}")
        print(f"    Match: {abs(W_complex.real - W8/8) < 1e-6 and abs(W_complex.imag) < 1e-6}")


if __name__ == "__main__":
    main()
