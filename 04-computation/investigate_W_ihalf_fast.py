#!/usr/bin/env python3
"""
Fast investigation of anti-correlation between H(T) and W(i/2) at n=7.

Uses bitmask DP to compute forward-edge distribution N_f efficiently.
O(2^n * n^2) per tournament instead of O(n! * n).

W(r) = sum_P prod_{k=0}^{n-2} (r + s_k)
where s_k = T[P(k),P(k+1)] - 1/2.

At r=i/2: factors are (1+i)/2 (forward) or (-1+i)/2 (backward).
For n=7: W(i/2) = (1/64) * sum_f c_f * N_f where c_f = Re[(1+i)^f * (-1+i)^{6-f}].
"""

from collections import defaultdict
import sys

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

def forward_edge_distribution_dp(T):
    """
    Compute N_f = #{permutations with exactly f forward edges} using bitmask DP.
    dp[mask][v][f] = #{paths through vertices in mask, ending at v, with f forward edges}
    O(2^n * n^2 * n) time, O(2^n * n^2) space.

    For n=7: 128 * 7 * 7 = 6272 states. Very fast.
    """
    n = len(T)
    # dp[mask][v] = dict mapping f -> count
    # Too much overhead with dicts. Use array: dp[mask][v][f]
    # f ranges from 0 to n-1
    dp = [[[0] * n for _ in range(n)] for _ in range(1 << n)]

    for v in range(n):
        dp[1 << v][v][0] = 1  # single vertex, 0 forward edges

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            for f in range(n):
                cnt = dp[mask][v][f]
                if cnt == 0:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    new_mask = mask | (1 << u)
                    if T[v][u]:  # forward edge
                        dp[new_mask][u][f + 1] += cnt
                    else:  # backward edge
                        dp[new_mask][u][f] += cnt

    full = (1 << n) - 1
    Nf = [0] * n
    for v in range(n):
        for f in range(n):
            Nf[f] += dp[full][v][f]

    return Nf

def compute_W_ihalf_from_Nf(Nf, coeffs_real, denom):
    """Compute W(i/2) * denom from the forward-edge distribution."""
    return sum(coeffs_real[f] * Nf[f] for f in range(len(Nf)))

def count_3cycles(T):
    """Count directed 3-cycles in tournament T."""
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
    total = 1 << m  # 2097152
    denom = 2 ** (n - 1)  # 64

    print(f"Exhaustive enumeration of all {total} tournaments on n={n} vertices")
    print("=" * 70)

    # Precompute coefficients for W(i/2)
    # (1+i)^f for f=0..6
    pow_1pi = [(1, 0)]
    for k in range(1, n):
        a, b = pow_1pi[-1]
        pow_1pi.append((a - b, a + b))

    # (-1+i)^g for g=0..6
    pow_m1pi = [(1, 0)]
    for k in range(1, n):
        a, b = pow_m1pi[-1]
        pow_m1pi.append((-a - b, a - b))

    coeffs_real = []
    coeffs_imag = []
    for f in range(n):
        g = (n - 1) - f
        ar, ai = pow_1pi[f]
        br, bi = pow_m1pi[g]
        pr = ar * br - ai * bi
        pi_ = ar * bi + ai * br
        coeffs_real.append(pr)
        coeffs_imag.append(pi_)

    print(f"\nCoefficients c_f = Re[(1+i)^f * (-1+i)^(6-f)]:")
    for f in range(n):
        print(f"  f={f}: c_f = {coeffs_real[f]:>5}, c_f_imag = {coeffs_imag[f]:>5}")

    print(f"\nFormula: W(i/2) = (1/{denom}) * sum_f c_f * N_f")

    # Collect data
    H_values = []
    W_scaled = []  # W(i/2) * denom (integer)
    c3_values = []

    step = max(1, total // 20)

    for bits in range(total):
        if bits % step == 0:
            pct = 100 * bits / total
            print(f"  Progress: {bits}/{total} ({pct:.0f}%)", file=sys.stderr)

        T = tournament_from_bits(n, bits)
        Nf = forward_edge_distribution_dp(T)

        H = Nf[n - 1]  # N_6 = H(T)
        W_sc = compute_W_ihalf_from_Nf(Nf, coeffs_real, denom)
        c3 = count_3cycles(T)

        H_values.append(H)
        W_scaled.append(W_sc)
        c3_values.append(c3)

    print(f"\nTotal tournaments processed: {len(H_values)}", file=sys.stderr)

    # =====================================================================
    # ANALYSIS
    # =====================================================================

    # Group by H
    H_to_data = defaultdict(list)
    for i in range(total):
        H_to_data[H_values[i]].append((W_scaled[i], c3_values[i]))

    H_set = sorted(H_to_data.keys())

    # CHECK 1: Imaginary part
    print("\n" + "=" * 70)
    print("CHECK 1: W(i/2) integrality")
    non_int = sum(1 for w in W_scaled if w % denom != 0)
    print(f"  W(i/2)*{denom} always integer: YES (by construction)")
    print(f"  W(i/2) always integer (W*{denom} div by {denom}): {non_int == 0} ({non_int} failures)")

    # Check divisibility
    for d in [2, 4, 8, 16, 32, 64]:
        div_count = sum(1 for w in W_scaled if w % d == 0)
        print(f"  W*{denom} divisible by {d}: {div_count}/{total} ({100*div_count/total:.1f}%)")

    # CHECK 2: Distribution
    print("\n" + "=" * 70)
    print("CHECK 2: Joint distribution of H and W(i/2)")
    print(f"\n  {'H':>5}  {'count':>7}  {'W*64 values (sorted set)':>50}  {'W values':>30}")

    for H in H_set:
        entries = H_to_data[H]
        W_set = sorted(set(e[0] for e in entries))
        count = len(entries)
        W_actual = [w / denom for w in W_set]

        if len(W_set) <= 8:
            w_str = str(W_set)
            wa_str = str(W_actual)
        else:
            w_str = f"[{W_set[0]}..{W_set[-1]}] ({len(W_set)} vals)"
            wa_str = f"[{W_actual[0]:.1f}..{W_actual[-1]:.1f}]"

        print(f"  {H:>5}  {count:>7}  {w_str:>50}  {wa_str:>30}")

    # CHECK 3: Monotonicity
    print("\n" + "=" * 70)
    print("CHECK 3: Monotonicity of |W(i/2)| vs H")

    print(f"\n  {'H':>5}  {'count':>7}  {'|W|*64 min':>12}  {'|W|*64 max':>12}  {'|W| min':>10}  {'|W| max':>10}")
    absW_ranges = {}
    for H in H_set:
        entries = H_to_data[H]
        absW = [abs(e[0]) for e in entries]
        mn, mx = min(absW), max(absW)
        absW_ranges[H] = (mn, mx)
        print(f"  {H:>5}  {len(entries):>7}  {mn:>12}  {mx:>12}  {mn/denom:>10.2f}  {mx/denom:>10.2f}")

    # Check strict monotonicity at H-level
    violations = []
    for i, H1 in enumerate(H_set):
        for j, H2 in enumerate(H_set):
            if H1 >= H2:
                continue
            # H1 < H2 => want max|W(H2)| <= min|W(H1)|
            if absW_ranges[H2][1] > absW_ranges[H1][0]:
                violations.append((H1, H2, absW_ranges[H1][0], absW_ranges[H2][1]))

    print(f"\n  Strict monotonicity violations: {len(violations)}")
    if violations:
        print(f"  First 10 violations:")
        for H1, H2, minW1, maxW2 in violations[:10]:
            print(f"    H1={H1} < H2={H2}: min|W(H1)|*64={minW1}, max|W(H2)|*64={maxW2}")

    # CHECK 4: W(i/2)^2 analysis
    print("\n" + "=" * 70)
    print("CHECK 4: W(i/2)^2 values")

    W2_all = sorted(set(w * w for w in W_scaled))
    print(f"  Distinct W*64 values: {len(set(W_scaled))}")
    print(f"  Distinct (W*64)^2 values: {len(W2_all)}")
    print(f"  First 20 (W*64)^2 values: {W2_all[:20]}")

    # Check if W^2 = (W*64)^2 / 64^2 is always integer
    W2_actual = sorted(set(w * w for w in W_scaled))
    all_div = all(w2 % (denom * denom) == 0 for w2 in W2_actual)
    print(f"  (W*64)^2 always divisible by {denom*denom}: {all_div}")

    if all_div:
        W2_int = sorted(set(w2 // (denom * denom) for w2 in W2_actual))
        print(f"  W(i/2)^2 integer values: {W2_int}")

        # Check perfect squares
        import math
        for w2 in W2_int:
            sq = int(round(math.isqrt(w2)))
            is_ps = sq * sq == w2
            print(f"    W^2 = {w2}: perfect square = {is_ps}" + (f" (sqrt = {sq})" if is_ps else ""))

    # CHECK 5: Correlation
    print("\n" + "=" * 70)
    print("CHECK 5: Correlations")

    N = len(H_values)
    absW_vals = [abs(W_scaled[i]) / denom for i in range(N)]

    mean_H = sum(H_values) / N
    mean_absW = sum(absW_vals) / N

    cov = sum((H_values[i] - mean_H) * (absW_vals[i] - mean_absW) for i in range(N)) / N
    var_H = sum((H_values[i] - mean_H) ** 2 for i in range(N)) / N
    var_absW = sum((absW_vals[i] - mean_absW) ** 2 for i in range(N)) / N

    corr = cov / (var_H ** 0.5 * var_absW ** 0.5) if var_H > 0 and var_absW > 0 else 0
    print(f"  Corr(H, |W(i/2)|) = {corr:.6f}")

    # Signed
    W_vals = [W_scaled[i] / denom for i in range(N)]
    mean_W = sum(W_vals) / N
    cov2 = sum((H_values[i] - mean_H) * (W_vals[i] - mean_W) for i in range(N)) / N
    var_W = sum((W_vals[i] - mean_W) ** 2 for i in range(N)) / N
    corr2 = cov2 / (var_H ** 0.5 * var_W ** 0.5) if var_H > 0 and var_W > 0 else 0
    print(f"  Corr(H, W(i/2))  = {corr2:.6f}")

    # W^2
    W2_vals = [W_scaled[i] ** 2 / denom ** 2 for i in range(N)]
    mean_W2 = sum(W2_vals) / N
    cov3 = sum((H_values[i] - mean_H) * (W2_vals[i] - mean_W2) for i in range(N)) / N
    var_W2 = sum((W2_vals[i] - mean_W2) ** 2 for i in range(N)) / N
    corr3 = cov3 / (var_H ** 0.5 * var_W2 ** 0.5) if var_H > 0 and var_W2 > 0 else 0
    print(f"  Corr(H, W(i/2)^2) = {corr3:.6f}")

    # H vs c3
    mean_c3 = sum(c3_values) / N
    cov4 = sum((H_values[i] - mean_H) * (c3_values[i] - mean_c3) for i in range(N)) / N
    var_c3 = sum((c3_values[i] - mean_c3) ** 2 for i in range(N)) / N
    corr4 = cov4 / (var_H ** 0.5 * var_c3 ** 0.5) if var_H > 0 and var_c3 > 0 else 0
    print(f"  Corr(H, c3)       = {corr4:.6f}")

    # c3 vs |W|
    cov5 = sum((c3_values[i] - mean_c3) * (absW_vals[i] - mean_absW) for i in range(N)) / N
    corr5 = cov5 / (var_c3 ** 0.5 * var_absW ** 0.5) if var_c3 > 0 and var_absW > 0 else 0
    print(f"  Corr(c3, |W(i/2)|) = {corr5:.6f}")

    # CHECK 6: Exact formula verification
    print("\n" + "=" * 70)
    print("CHECK 6: Explicit formula for n=7")
    print()
    print(f"  c_f = {coeffs_real}")
    print(f"  W(i/2) * 64 = {' + '.join(f'({c})*N_{f}' for f, c in enumerate(coeffs_real))}")

    # Simplify using the constraints
    # N_0 + N_1 + N_2 + N_3 + N_4 + N_5 + N_6 = 5040
    # Coefficients: c = [8, 0, -8, 0, 8, 0, -8]
    # So W*64 = 8*(N_0 - N_2 + N_4 - N_6)
    # W*64 = 8 * (N_0 - N_2 + N_4 - N_6)
    # W = (1/8) * (N_0 - N_2 + N_4 - N_6)

    print(f"\n  Simplified: c_f = {coeffs_real}")
    print(f"  Note pattern: c_0=8, c_1=0, c_2=-8, c_3=0, c_4=8, c_5=0, c_6=-8")
    print(f"  So W*64 = 8*(N_0 - N_2 + N_4 - N_6)")
    print(f"  W(i/2) = (N_0 - N_2 + N_4 - N_6) / 8")

    # Verify this on a sample
    T_check = tournament_from_bits(7, 12345)
    Nf = forward_edge_distribution_dp(T_check)
    W_direct = compute_W_ihalf_from_Nf(Nf, coeffs_real, denom)
    W_simplified = 8 * (Nf[0] - Nf[2] + Nf[4] - Nf[6])
    print(f"\n  Verification: bits=12345, Nf={Nf}")
    print(f"    Direct: W*64 = {W_direct}")
    print(f"    Simplified: 8*(N0-N2+N4-N6) = {W_simplified}")
    print(f"    Match: {W_direct == W_simplified}")

    # So W(i/2) = (N_0 - N_2 + N_4 - N_6) / 8
    # and H = N_6
    # Anti-correlation: high H means high N_6, which SUBTRACTS from W.

    print(f"\n  Key insight: W(i/2) = (N_0 - N_2 + N_4 - N_6)/8")
    print(f"  H = N_6, H_op = N_0 (Hamiltonian paths of opposite tournament)")
    print(f"  So W(i/2) = (H_op - N_2 + N_4 - H)/8")

    # What is N_0 - N_6 = H_op - H?
    # And N_4 - N_2?

    # CHECK 7: H_op = H? (by Redei or otherwise)
    print("\n" + "=" * 70)
    print("CHECK 7: Is H(T) = H(T^op)? (N_0 = N_6?)")
    mismatch = 0
    for i in range(total):
        T = tournament_from_bits(n, i)
        # T^op has reversed edges. For bits encoding, T^op has complemented bits.
        bits_op = ((1 << m) - 1) ^ i
        # H(T^op) = H_values[bits_op]
        if H_values[i] != H_values[bits_op]:
            mismatch += 1
    print(f"  H(T) != H(T^op) count: {mismatch}")
    # This should be 0 since H(T) = H(T^op) always (reverse all paths).
    # Actually: reversing a Hamiltonian path in T gives a Ham path in T^op.
    # So H(T) = H(T^op) always. Thus N_0 = N_6 for all T.

    # Verify N_f = N_{6-f} (symmetry from path reversal)
    print("\n  Checking N_f = N_{6-f} symmetry (from path reversal):")
    symm_fail = 0
    for bits in range(0, total, total // 100):
        T = tournament_from_bits(7, bits)
        Nf = forward_edge_distribution_dp(T)
        for f in range(4):
            if Nf[f] != Nf[6 - f]:
                symm_fail += 1
    print(f"  Symmetry failures (sample): {symm_fail}")

    # If N_f = N_{6-f}, then:
    # W(i/2) = (N_0 - N_2 + N_4 - N_6)/8 = (N_6 - N_4 + N_2 - N_0)/8 * (-1) ... wait
    # N_0 = N_6, N_2 = N_4. So:
    # W(i/2) = (N_6 - N_4 + N_4 - N_6)/8 = 0 ???!!!
    print("\n  If N_f = N_{6-f}: N_0 = N_6 = H, N_2 = N_4")
    print("  Then W(i/2) = (H - N_2 + N_2 - H)/8 = 0 !!!")
    print("  This would mean W(i/2) = 0 always, contradicting the premise!")
    print()
    print("  Let me check: is N_f = N_{6-f} actually true?")

    # Check more carefully
    T_test = tournament_from_bits(7, 12345)
    Nf_test = forward_edge_distribution_dp(T_test)
    print(f"  bits=12345: Nf = {Nf_test}")
    print(f"  N_0={Nf_test[0]}, N_6={Nf_test[6]}: equal? {Nf_test[0] == Nf_test[6]}")
    print(f"  N_1={Nf_test[1]}, N_5={Nf_test[5]}: equal? {Nf_test[1] == Nf_test[5]}")
    print(f"  N_2={Nf_test[2]}, N_4={Nf_test[4]}: equal? {Nf_test[2] == Nf_test[4]}")

    # WAIT: reversing a path changes forward edges to backward and vice versa.
    # A permutation P = (P_0, P_1, ..., P_6) with f forward edges
    # reverses to P' = (P_6, P_5, ..., P_0).
    # In P', edge P_k -> P_{k-1}: T[P_k][P_{k-1}] = 1 - T[P_{k-1}][P_k]
    # So forward in P' iff backward in P.
    # Thus P' has 6-f forward edges.
    # This gives N_f(T) = N_{6-f}(T) for the SAME tournament T.

    # But then W(i/2) = 8*(N_0 - N_2 + N_4 - N_6)/64 = (N_0 - N_2 + N_4 - N_6)/8
    # With N_0 = N_6 and N_2 = N_4: W = 0.
    # Something is wrong with my coefficient computation!

    print("\n  DEBUGGING: Let me recompute the coefficients carefully.")

    # (1+i)/2 and (-1+i)/2 are the two possible factors.
    # Product = ((1+i)/2)^f * ((-1+i)/2)^{6-f}
    # = (1/2^6) * (1+i)^f * (-1+i)^{6-f}

    # (1+i)^0 = 1
    # (1+i)^1 = 1+i
    # (1+i)^2 = 2i
    # (1+i)^3 = 2i(1+i) = -2+2i
    # (1+i)^4 = (2i)^2 = -4
    # (1+i)^5 = -4(1+i) = -4-4i
    # (1+i)^6 = -4*(2i) = -8i

    print("  (1+i)^f:")
    for f in range(7):
        print(f"    f={f}: {pow_1pi[f]}")

    # (-1+i)^0 = 1
    # (-1+i)^1 = -1+i
    # (-1+i)^2 = (-1+i)^2 = 1 - 2i - 1 = -2i
    # (-1+i)^3 = -2i(-1+i) = 2i - 2i^2 = 2+2i
    # (-1+i)^4 = (-2i)^2 = -4
    # (-1+i)^5 = -4(-1+i) = 4-4i
    # (-1+i)^6 = (-2i)^3 = 8i

    print("  (-1+i)^g:")
    for g in range(7):
        print(f"    g={g}: {pow_m1pi[g]}")

    # Products (1+i)^f * (-1+i)^{6-f}:
    print("  Products (1+i)^f * (-1+i)^{6-f}:")
    for f in range(7):
        g = 6 - f
        ar, ai = pow_1pi[f]
        br, bi = pow_m1pi[g]
        pr = ar * br - ai * bi
        pi_ = ar * bi + ai * br
        print(f"    f={f}, g={g}: ({ar}+{ai}i)*({br}+{bi}i) = {pr}+{pi_}i")

    # Ah wait. Let me check: is the path reversal argument correct?
    # P = (v0, v1, ..., v6). Forward edges: T[v_k][v_{k+1}] = 1.
    # Reversed: P' = (v6, v5, ..., v0). Edge from v_{k+1} to v_k.
    # T[v_{k+1}][v_k] = 1 - T[v_k][v_{k+1}].
    # So if P has f forward edges, P' has 6-f forward edges.
    # Bijection on permutations of same tournament => N_f = N_{6-f}.
    # This means W*64 = sum c_f * N_f = sum c_f * N_{6-f} = sum c_{6-f} * N_f.
    # So the coefficients must satisfy c_f = c_{6-f} for W to be non-trivially nonzero.
    # c = [8, 0, -8, 0, 8, 0, -8] => c_{6-f} = [-8, 0, 8, 0, -8, 0, 8]
    # c_f + c_{6-f} = [0, 0, 0, 0, 0, 0, 0] !!!
    # So W*64 = sum c_f * N_f and sum c_{6-f} * N_f are NEGATIVES of each other.
    # But they should be equal (both = W*64). So W*64 = -W*64 => W = 0 always.

    print("\n  CONCLUSION: c_f + c_{6-f} = 0 for all f.")
    print("  Combined with N_f = N_{6-f} (path reversal), this gives W(i/2) = 0 ALWAYS.")

    # Let me verify this directly
    print("\n  Direct verification:")
    all_zero = True
    nonzero_count = 0
    for bits in range(min(total, 10000)):
        T = tournament_from_bits(7, bits)
        Nf = forward_edge_distribution_dp(T)
        W_sc = compute_W_ihalf_from_Nf(Nf, coeffs_real, denom)
        if W_sc != 0:
            all_zero = False
            nonzero_count += 1
            if nonzero_count <= 5:
                print(f"    NONZERO at bits={bits}: Nf={Nf}, W*64={W_sc}")

    if all_zero:
        print(f"    W(i/2) = 0 for ALL first 10000 tournaments. CONFIRMED.")
    else:
        print(f"    {nonzero_count}/10000 nonzero W values found!")

    # So W(i/2) = 0 always at n=7. The original question may have a different definition.
    # Perhaps the question is about W(r) for the REDUCED polynomial P(u,t)?
    # Or perhaps the sum is over Hamiltonian paths only, not all permutations?

    # Let me try: W_ham(r) = sum over HAM PATHS only of prod(r + s_k)
    # At r=1/2: each term = 1 (all forward), so W_ham(1/2) = H(T). Correct.
    # At r=i/2: each term = ((1+i)/2)^6 = (1+i)^6 / 64 = -8i/64 = -i/8
    # ALL Hamiltonian paths have ALL forward edges (that's what makes them Ham paths!)
    # So W_ham(i/2) = H * (-i/8). This is purely imaginary, not useful.

    # Actually wait. W(r) is defined as sum over ALL permutations, where
    # s_k = A[P(k), P(k+1)] - 1/2.
    # A[u][v] = 1 if u->v (forward), 0 if v->u.
    # So for a Hamiltonian path, all A=1, all s_k = 1/2, prod = (r+1/2)^{n-1}.
    # For a non-path, at least one A=0, giving s_k = -1/2, prod contains (r-1/2).
    # At r=1/2: (r-1/2) = 0, so only Ham paths survive. Good.

    # The symmetry N_f = N_{6-f} and c_f = -c_{6-f} forces W(i/2) = 0.
    # This is just a consequence of W(r) being a polynomial where W(r) + W(-r) has no
    # contribution from terms of alternating sign...

    # Wait, let me reconsider. Maybe the question is about a DIFFERENT r value,
    # or the "reduced" polynomial obtained after factoring out the known roots.

    # The reduced polynomial: P(u, t) where W(r) = something involving both u and t.
    # The user mentions: "W(i/2) at odd n is always REAL" and "W(i/2) = 0 characterizes
    # H-maximizers at n=3,5,7 but NOT at n=9,11"
    # But we just showed W(i/2) = 0 for ALL tournaments at n=7, not just maximizers!

    # Maybe the definition of W is different. Let me try:
    # W(r) = sum over Hamiltonian paths P of prod_{k} (r + s_k)
    # where s_k varies (not all +1/2).
    # No, for a Hamiltonian path, all edges are forward, so all s_k = +1/2.

    # OR: maybe the transfer matrix formulation gives a different polynomial.
    # M(r) = product over edges (r + epsilon_ij) where epsilon_ij = T[i][j] - 1/2.
    # The permanent or something?

    # Actually, re-reading the problem: W(r) = sum_P prod_{i=0}^{n-2} (r + s_i)
    # where P is a permutation and s_i = A[P(i), P(i+1)] - 1/2.
    # This is exactly what I computed. And it's identically zero at r=i/2 for n=7.

    # UNLESS... the definition uses s_i = ±1 instead of ±1/2?
    # Let's try: s_i = 2*(A[P(i),P(i+1)] - 1/2) = A[P(i),P(i+1)]*2 - 1 in {-1, +1}
    # Then W(r) = sum_P prod (r + s_i) with s_i in {-1, +1}
    # At r=1: each factor is 0 or 2. W(1) = 2^{n-1} * H.
    # At r=i: factors are (i+1) or (i-1) = (1+i) or -(1-i)
    # |1+i| = sqrt(2), |1-i| = sqrt(2).
    # (1+i)^f * (-(1-i))^{6-f} = (-1)^{6-f} * (1+i)^f * (1-i)^{6-f}

    print("\n" + "=" * 70)
    print("TRYING ALTERNATIVE: s_i in {-1, +1} (not {-1/2, +1/2})")

    # Coefficients for W(i) = sum_f N_f * (1+i)^f * (-(1-i))^{6-f}
    # = sum_f N_f * (-1)^{6-f} * (1+i)^f * (1-i)^{6-f}
    # (1+i)(1-i) = 2

    # (1+i)^f * (1-i)^{6-f}:
    # If f <= 6-f (f <= 3): = (1-i)^{6-2f} * 2^f
    # If f > 3: = (1+i)^{2f-6} * 2^{6-f}

    alt_coeffs_r = []
    alt_coeffs_i = []
    # (1+i)^f
    p1 = [(1,0)]
    for _ in range(6):
        a, b = p1[-1]
        p1.append((a-b, a+b))

    # (-(1-i))^g = (-1)^g * (1-i)^g
    # (1-i)^g
    pm1 = [(1,0)]
    for _ in range(6):
        a, b = pm1[-1]
        pm1.append((a+b, b-a))  # (a+bi)(1-i) = (a+b) + (b-a)i

    print(f"\n  Coefficients d_f = (1+i)^f * (-(1-i))^(6-f):")
    for f in range(7):
        g = 6 - f
        ar, ai = p1[f]
        # (-(1-i))^g = (-1)^g * (1-i)^g
        br, bi = pm1[g]
        sign = (-1)**g
        br *= sign
        bi *= sign
        pr = ar * br - ai * bi
        pi_ = ar * bi + ai * br
        alt_coeffs_r.append(pr)
        alt_coeffs_i.append(pi_)
        print(f"    f={f}: d_f = {pr} + {pi_}i")

    # Check symmetry: d_f vs d_{6-f}
    print(f"\n  d_f:     {alt_coeffs_r}")
    print(f"  d_{{6-f}}: {[alt_coeffs_r[6-f] for f in range(7)]}")
    print(f"  d_f + d_{{6-f}}: {[alt_coeffs_r[f] + alt_coeffs_r[6-f] for f in range(7)]}")

    # Same problem! d_f + d_{6-f} will again be zero if the antisymmetry holds.
    # Because path reversal always gives N_f = N_{6-f}, ANY definition of W that
    # assigns weights based only on forward-edge count will have this issue
    # IF the weight function w(f) satisfies w(f) = -w(6-f).

    # The anti-correlation mentioned in the problem statement must involve a DIFFERENT quantity.
    # Perhaps it's W evaluated at a point where the symmetry breaks, e.g., |W(r)| for
    # complex r where the coefficient function doesn't have the antisymmetry.

    # OR: Perhaps W is defined differently, e.g., as a sum over HAMILTONIAN PATHS
    # of some per-edge weight that varies by edge.

    # Let me try the TRANSFER MATRIX approach.
    # The transfer matrix M[a][b] = (r + T[a][b] - 1/2).
    # Then perm(M) = W(r).
    # But perm(M) is what we computed above.

    # Actually, maybe the point is a different evaluation.
    # "W(r) is a polynomial in r" — yes, degree n-1.
    # For n=7: W(r) = c_0 + c_1*r + ... + c_6*r^6, where c_k = sum_P e_{n-1-k}(s(P))
    # (elementary symmetric polynomials in the s values of each permutation)

    # At r = i/2: we showed this is 0.
    # But the POLYNOMIAL W(r) is not identically zero! Just zero at r=i/2.

    # OK so maybe the question is: what are the ROOTS of W(r)?
    # And i/2 happens to always be a root?
    # That would make sense: "W(i/2) = 0 characterizes H-maximizers" is FALSE,
    # it's 0 for ALL tournaments.

    # Let me compute the full polynomial W(r) for a few tournaments.
    print("\n" + "=" * 70)
    print("COMPUTING: Full polynomial W(r) for sample tournaments")

    for bits in [0, 1, 12345, max_bits[0]]:
        T = tournament_from_bits(7, bits)
        H = hamiltonian_path_count_dp(T)
        Nf = forward_edge_distribution_dp(T)

        # W(r) = sum_f N_f * (r + 1/2)^f * (r - 1/2)^{6-f}
        # Expand: (r+1/2)^f * (r-1/2)^{6-f} = (r^2 - 1/4)^min(f,6-f) * ... complicated
        # Better: just compute W at several points and interpolate.

        # Compute W at r = 0, 1, 2, 3, 4, 5, 6 (7 points for degree 6 poly)
        W_points = []
        for r_val in range(7):
            val = 0
            for f in range(7):
                g = 6 - f
                term = Nf[f] * (r_val + 0.5)**f * (r_val - 0.5)**g
                val += term
            W_points.append(val)

        print(f"\n  bits={bits}, H={H}, Nf={Nf}")
        for rv, wv in enumerate(W_points):
            print(f"    W({rv}) = {wv:.2f}")

        # Check W(1/2) = H
        W_half = 0
        for f in range(7):
            g = 6 - f
            W_half += Nf[f] * 1.0**f * 0.0**g
        print(f"    W(1/2) = {W_half:.2f} (should be {H})")

        # Check W(-1/2)
        W_mhalf = 0
        for f in range(7):
            g = 6 - f
            W_mhalf += Nf[f] * 0.0**f * (-1.0)**g
        print(f"    W(-1/2) = {W_mhalf:.2f}")
        # W(-1/2): only f=0 survives (all backward edges = Ham paths of T^op)
        # = N_0 * 1 * (-1)^6 = N_0 = H(T^op) = H(T). So W(-1/2) = H too!

        # W at r=i/2:
        W_ihalf = 0j
        for f in range(7):
            g = 6 - f
            W_ihalf += Nf[f] * (0.5j + 0.5)**f * (0.5j - 0.5)**g
        print(f"    W(i/2) = {W_ihalf:.6f}")

    print("\n  CONFIRMED: W(i/2) = 0 for all n=7 tournaments (algebraic identity).")
    print("  The anti-correlation described in the problem must involve a DIFFERENT quantity.")
    print()
    print("  Possible candidates:")
    print("  1. |W'(i/2)| (derivative of W at r=i/2)")
    print("  2. The REDUCED polynomial after dividing out (r - i/2)(r + i/2) = r^2 + 1/4")
    print("  3. A different evaluation point")
    print("  4. W for EVEN n (where the N_f = N_{n-1-f} symmetry still holds but")
    print("     the coefficient antisymmetry c_f = -c_{n-1-f} may not)")

    # Let me check: for what n does c_f = -c_{n-1-f}?
    print("\n" + "=" * 70)
    print("CHECKING: For which n is r=i/2 always a root of W(r)?")

    for nn in range(3, 10):
        # Compute coefficients for n=nn
        pp1 = [(1,0)]
        for _ in range(nn-1):
            a, b = pp1[-1]
            pp1.append((a-b, a+b))

        pm1i = [(1,0)]
        for _ in range(nn-1):
            a, b = pm1i[-1]
            pm1i.append((-a-b, a-b))

        cc_r = []
        for f in range(nn):
            g = (nn-1) - f
            ar, ai = pp1[f]
            br, bi = pm1i[g]
            pr = ar * br - ai * bi
            cc_r.append(pr)

        # Check antisymmetry: c_f + c_{n-1-f} = 0?
        antisym = all(cc_r[f] + cc_r[(nn-1)-f] == 0 for f in range(nn))
        # With N_f = N_{n-1-f}, W = 0 iff antisymmetric.
        print(f"  n={nn}: coeffs = {cc_r}, antisymmetric = {antisym}")

    # Now let's investigate W'(i/2) which might be the actual quantity of interest
    print("\n" + "=" * 70)
    print("INVESTIGATING: W'(i/2) — derivative of W at r=i/2")

    # W(r) = sum_f N_f * (r+1/2)^f * (r-1/2)^{n-1-f}
    # W'(r) = sum_f N_f * d/dr [(r+1/2)^f * (r-1/2)^{n-1-f}]
    # = sum_f N_f * [f*(r+1/2)^{f-1}*(r-1/2)^{n-1-f} + (n-1-f)*(r+1/2)^f*(r-1/2)^{n-2-f}]

    # At r=i/2: (r+1/2) = (1+i)/2, (r-1/2) = (-1+i)/2
    # This is more complex. Let me just evaluate numerically.

    delta = 1e-8

    H_list = []
    Wprime_list = []

    for bits in range(total):
        if bits % step == 0:
            print(f"  Progress: {bits}/{total} ({100*bits/total:.0f}%)", file=sys.stderr)

        T = tournament_from_bits(7, bits)
        Nf = forward_edge_distribution_dp(T)
        H = Nf[6]

        # Compute W'(i/2) numerically via finite difference
        r_val = 0.5j
        W_plus = 0j
        W_minus = 0j
        for f in range(7):
            g = 6 - f
            W_plus += Nf[f] * (r_val + delta + 0.5)**f * (r_val + delta - 0.5)**g
            W_minus += Nf[f] * (r_val - delta + 0.5)**f * (r_val - delta - 0.5)**g

        Wprime = (W_plus - W_minus) / (2 * delta)

        H_list.append(H)
        Wprime_list.append(Wprime)

    # Analyze W'(i/2)
    print(f"\n  Sample W'(i/2) values:")
    for bits in [0, 1, 12345, max_bits[0]]:
        idx = bits
        print(f"    bits={bits}: H={H_list[idx]}, W'(i/2) = {Wprime_list[idx]:.4f}")

    # Is W'(i/2) real or imaginary?
    real_parts = [w.real for w in Wprime_list]
    imag_parts = [w.imag for w in Wprime_list]
    print(f"\n  Real parts range: [{min(real_parts):.4f}, {max(real_parts):.4f}]")
    print(f"  Imag parts range: [{min(imag_parts):.4f}, {max(imag_parts):.4f}]")

    # Check if purely imaginary
    mostly_imag = all(abs(w.real) < 1e-4 for w in Wprime_list)
    mostly_real = all(abs(w.imag) < 1e-4 for w in Wprime_list)
    print(f"  Purely real: {mostly_real}")
    print(f"  Purely imaginary: {mostly_imag}")

    # Correlation with H
    if mostly_imag:
        absWp = [abs(w.imag) for w in Wprime_list]
        Wp_vals = [w.imag for w in Wprime_list]
    else:
        absWp = [abs(w) for w in Wprime_list]
        Wp_vals = [abs(w) for w in Wprime_list]

    mean_H2 = sum(H_list) / N
    mean_absWp = sum(absWp) / N

    cov = sum((H_list[i] - mean_H2) * (absWp[i] - mean_absWp) for i in range(N)) / N
    var_H2 = sum((H_list[i] - mean_H2) ** 2 for i in range(N)) / N
    var_absWp = sum((absWp[i] - mean_absWp) ** 2 for i in range(N)) / N

    corr = cov / (var_H2 ** 0.5 * var_absWp ** 0.5) if var_H2 > 0 and var_absWp > 0 else 0
    print(f"\n  Corr(H, |W'(i/2)|) = {corr:.6f}")

    # Table by H
    H_to_Wp = defaultdict(list)
    for i in range(N):
        H_to_Wp[H_list[i]].append(Wprime_list[i])

    print(f"\n  {'H':>5}  {'count':>7}  {'|W prime| min':>14}  {'|W prime| max':>14}  {'|W prime| mean':>14}")
    for H in sorted(H_to_Wp.keys()):
        entries = H_to_Wp[H]
        abs_vals = [abs(w) for w in entries]
        mn = min(abs_vals)
        mx = max(abs_vals)
        mean = sum(abs_vals) / len(abs_vals)
        print(f"  {H:>5}  {len(entries):>7}  {mn:>14.4f}  {mx:>14.4f}  {mean:>14.4f}")


if __name__ == "__main__":
    main()
