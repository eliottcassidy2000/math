#!/usr/bin/env python3
"""
Analyze W(i/2) and related quantities for tournaments.

Phase 1: Small n (3,5) exhaustive — verify W(i/2)=0 algebraic identity
Phase 2: n=7 exhaustive — compute W'(i/2) and reduced polynomial
Phase 3: Exact formulas

Key insight: W(r) = sum_P prod_{k=0}^{n-2} (r + s_k), s_k = T[P(k),P(k+1)] - 1/2.
Path reversal gives N_f = N_{n-1-f}.
At r=i/2, the coefficient c_f = -c_{n-1-f} (antisymmetric), so W(i/2)=0 always.

Therefore we investigate:
(a) The reduced polynomial Q(r) = W(r) / (r^2 + 1/4) and Q(i/2)
(b) W'(i/2)
(c) The relationship |W'(i/2)| vs H(T)
"""

from itertools import permutations, combinations
from collections import defaultdict
import sys
import cmath

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

def hamiltonian_path_count_dp(T):
    n = len(T)
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            cnt = dp[mask][v]
            if not (mask & (1 << v)) or cnt == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += cnt
    return sum(dp[full][v] for v in range(n))

def forward_edge_dist_perm(T):
    """Compute N_f by full permutation enumeration. Only for small n."""
    n = len(T)
    Nf = [0] * n
    for P in permutations(range(n)):
        f = sum(1 for k in range(n-1) if T[P[k]][P[k+1]])
        Nf[f] += 1
    return Nf

def W_polynomial_coeffs(Nf):
    """
    Compute the polynomial W(r) = sum_f N_f * (r+1/2)^f * (r-1/2)^{n-1-f}
    as a list of coefficients [a_0, a_1, ..., a_{n-1}] where W(r) = sum a_k r^k.
    Uses exact integer arithmetic (multiply by 2^{n-1} to clear denominators).
    """
    n = len(Nf)
    deg = n - 1

    # We'll work with 2^deg * W(r) to stay in integers.
    # 2^deg * (r+1/2)^f * (r-1/2)^g = (2r+1)^f * (2r-1)^g
    # Let u = 2r. Then = (u+1)^f * (u-1)^g.
    # This is a polynomial in u of degree f+g = deg.
    # W_scaled(u) = sum_f N_f * (u+1)^f * (u-1)^{deg-f}

    # Expand (u+1)^f = sum_j C(f,j) u^j
    # Expand (u-1)^g = sum_k C(g,k) (-1)^{g-k} u^k
    # Product: coeff of u^m = sum_{j+k=m} C(f,j) * C(g,k) * (-1)^{g-k}

    from math import comb

    # Polynomial in u (= 2r)
    poly_u = [0] * (deg + 1)
    for f in range(n):
        g = deg - f
        if Nf[f] == 0:
            continue
        for m in range(deg + 1):
            coeff_m = 0
            for j in range(max(0, m - g), min(f, m) + 1):
                k = m - j
                coeff_m += comb(f, j) * comb(g, k) * ((-1) ** (g - k))
            poly_u[m] += Nf[f] * coeff_m

    return poly_u  # poly_u[m] = coefficient of u^m in 2^deg * W(r) where u=2r

def eval_poly(coeffs, x):
    """Evaluate polynomial with given coefficients at x."""
    result = 0
    xpow = 1
    for c in coeffs:
        result += c * xpow
        xpow *= x
    return result

def poly_derivative(coeffs):
    """Derivative of polynomial."""
    return [coeffs[k] * k for k in range(1, len(coeffs))]

def poly_divide_by_u2_plus_1(coeffs):
    """
    Divide polynomial in u by (u^2 + 1). Return quotient and remainder.
    coeffs = [a_0, a_1, ..., a_d]
    """
    d = len(coeffs) - 1
    if d < 2:
        return [0], coeffs[:]

    # Synthetic division by u^2 + 1
    q = [0] * (d - 1)
    r = coeffs[:]

    for i in range(d, 1, -1):
        q[i-2] = r[i]
        r[i] = 0
        r[i-2] -= q[i-2]  # subtract q[i-2] * 1 from r[i-2]

    remainder = r[:2]
    return q, remainder

def count_3cycles(T):
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

def analyze_n(n):
    """Full analysis for tournaments of size n."""
    m = n * (n - 1) // 2
    total = 1 << m
    deg = n - 1

    print(f"\n{'='*70}")
    print(f"ANALYSIS FOR n = {n}  ({total} tournaments, degree {deg} polynomial)")
    print(f"{'='*70}")

    # Phase 1: Verify W(i/2) = 0
    print(f"\n--- Phase 1: Verify W(i/2) = 0 ---")

    H_values = []
    W_at_ihalf = []
    Q_at_ihalf = []  # reduced polynomial Q(i/2) = W(r)/(r^2+1/4) at r=i/2
    Wprime_at_ihalf = []
    c3_values = []
    poly_coeffs_list = []

    for bits in range(total):
        T = tournament_from_bits(n, bits)
        Nf = forward_edge_dist_perm(T)
        H = Nf[deg]  # N_{n-1} = H(T)

        # Polynomial in u (=2r): 2^deg * W(r) = P(u)
        poly_u = W_polynomial_coeffs(Nf)

        # Evaluate at u = 2*(i/2) = i
        W_scaled_ihalf = eval_poly(poly_u, 1j)

        # Divide by u^2 + 1
        quotient, remainder = poly_divide_by_u2_plus_1(poly_u)

        # Q(u) = P(u) / (u^2 + 1), so Q(i) = quotient evaluated at u=i
        Q_ihalf = eval_poly(quotient, 1j)

        # W'(r) in terms of u: dW/dr = (1/2^deg) * dP/du * du/dr = (1/2^deg) * 2 * P'(u)
        # Actually: W(r) = P(2r)/2^deg. W'(r) = 2*P'(2r)/2^deg = P'(2r)/2^{deg-1}
        dpoly = poly_derivative(poly_u)
        Wp_scaled = eval_poly(dpoly, 1j)  # P'(i)
        # W'(i/2) = Wp_scaled / 2^{deg-1}

        H_values.append(H)
        W_at_ihalf.append(W_scaled_ihalf)
        Q_at_ihalf.append(Q_ihalf)
        Wprime_at_ihalf.append(Wp_scaled)
        c3_values.append(count_3cycles(T))
        poly_coeffs_list.append(poly_u)

    # Verify W(i/2) = 0
    max_W = max(abs(w) for w in W_at_ihalf)
    print(f"  max |W(i/2)| across all tournaments: {max_W:.2e}")
    if max_W < 1e-6:
        print(f"  CONFIRMED: W(i/2) = 0 for ALL n={n} tournaments")

    # Verify remainder of division by u^2+1 is always 0
    all_zero_rem = True
    for bits in range(total):
        poly_u = poly_coeffs_list[bits]
        _, rem = poly_divide_by_u2_plus_1(poly_u)
        if any(abs(r) > 1e-10 for r in rem):
            all_zero_rem = False
            print(f"  NON-ZERO remainder at bits={bits}: {rem}")
            break
    if all_zero_rem:
        print(f"  CONFIRMED: (u^2+1) divides P(u) for all n={n} tournaments")
        print(f"  i.e., (r^2 + 1/4) divides W(r) always")

    # Phase 2: Analyze Q(i/2) = quotient at u=i
    print(f"\n--- Phase 2: Reduced polynomial Q at u=i ---")

    # Q(i) should be well-defined. Let's check its properties.
    Q_real = [q.real for q in Q_at_ihalf]
    Q_imag = [q.imag for q in Q_at_ihalf]

    print(f"  Q(i) real parts range: [{min(Q_real):.4f}, {max(Q_real):.4f}]")
    print(f"  Q(i) imag parts range: [{min(Q_imag):.4f}, {max(Q_imag):.4f}]")

    # Is Q(i) purely real or purely imaginary?
    all_real = all(abs(q.imag) < 0.01 for q in Q_at_ihalf)
    all_imag = all(abs(q.real) < 0.01 for q in Q_at_ihalf)
    print(f"  Q(i) purely real: {all_real}")
    print(f"  Q(i) purely imaginary: {all_imag}")

    # Phase 3: W'(i/2) analysis
    print(f"\n--- Phase 3: W'(i/2) analysis ---")
    Wp_real = [w.real for w in Wprime_at_ihalf]
    Wp_imag = [w.imag for w in Wprime_at_ihalf]

    print(f"  P'(i) real parts range: [{min(Wp_real):.4f}, {max(Wp_real):.4f}]")
    print(f"  P'(i) imag parts range: [{min(Wp_imag):.4f}, {max(Wp_imag):.4f}]")

    all_real_wp = all(abs(w.real) < 0.01 for w in Wprime_at_ihalf)
    all_imag_wp = all(abs(w.imag) < 0.01 for w in Wprime_at_ihalf)
    print(f"  P'(i) purely real: {all_real_wp}")
    print(f"  P'(i) purely imaginary: {all_imag_wp}")

    # Phase 4: Correlations
    print(f"\n--- Phase 4: Correlations ---")

    N = total
    # Use |Q(i)| or appropriate component
    if all_real:
        Q_vals = Q_real
        Q_label = "Re Q(i)"
    elif all_imag:
        Q_vals = [abs(q) for q in Q_imag]
        Q_label = "|Im Q(i)|"
    else:
        Q_vals = [abs(q) for q in Q_at_ihalf]
        Q_label = "|Q(i)|"

    mean_H = sum(H_values) / N
    mean_Q = sum(Q_vals) / N

    cov = sum((H_values[i] - mean_H) * (Q_vals[i] - mean_Q) for i in range(N)) / N
    var_H = sum((H_values[i] - mean_H) ** 2 for i in range(N)) / N
    var_Q = sum((Q_vals[i] - mean_Q) ** 2 for i in range(N)) / N

    if var_H > 0 and var_Q > 0:
        corr = cov / (var_H ** 0.5 * var_Q ** 0.5)
    else:
        corr = 0
    print(f"  Corr(H, {Q_label}) = {corr:.6f}")

    # With absolute value
    absQ_vals = [abs(q) for q in Q_at_ihalf]
    mean_absQ = sum(absQ_vals) / N
    cov_abs = sum((H_values[i] - mean_H) * (absQ_vals[i] - mean_absQ) for i in range(N)) / N
    var_absQ = sum((absQ_vals[i] - mean_absQ) ** 2 for i in range(N)) / N
    if var_H > 0 and var_absQ > 0:
        corr_abs = cov_abs / (var_H ** 0.5 * var_absQ ** 0.5)
    else:
        corr_abs = 0
    print(f"  Corr(H, |Q(i)|) = {corr_abs:.6f}")

    # |W'|
    absWp_vals = [abs(w) for w in Wprime_at_ihalf]
    mean_absWp = sum(absWp_vals) / N
    cov_wp = sum((H_values[i] - mean_H) * (absWp_vals[i] - mean_absWp) for i in range(N)) / N
    var_absWp = sum((absWp_vals[i] - mean_absWp) ** 2 for i in range(N)) / N
    if var_H > 0 and var_absWp > 0:
        corr_wp = cov_wp / (var_H ** 0.5 * var_absWp ** 0.5)
    else:
        corr_wp = 0
    print(f"  Corr(H, |P'(i)|) = {corr_wp:.6f}")

    # c3
    mean_c3 = sum(c3_values) / N
    cov_c3 = sum((H_values[i] - mean_H) * (c3_values[i] - mean_c3) for i in range(N)) / N
    var_c3 = sum((c3_values[i] - mean_c3) ** 2 for i in range(N)) / N
    if var_H > 0 and var_c3 > 0:
        corr_c3 = cov_c3 / (var_H ** 0.5 * var_c3 ** 0.5)
    else:
        corr_c3 = 0
    print(f"  Corr(H, c3) = {corr_c3:.6f}")

    # Phase 5: Table by H
    print(f"\n--- Phase 5: Table by H ---")

    H_to_data = defaultdict(list)
    for i in range(total):
        H_to_data[H_values[i]].append({
            'Q': Q_at_ihalf[i],
            'Wp': Wprime_at_ihalf[i],
            'c3': c3_values[i],
            'poly': poly_coeffs_list[i]
        })

    H_set = sorted(H_to_data.keys())
    denom_factor = 2 ** (deg - 1) if deg > 1 else 1

    print(f"\n  {'H':>5}  {'cnt':>5}  {'|Q(i)| rng':>25}  {'|Wprime| rng':>25}  {'c3 rng':>12}")
    for H in H_set:
        entries = H_to_data[H]
        absQ = [abs(e['Q']) for e in entries]
        absWp = [abs(e['Wp']) for e in entries]
        c3s = [e['c3'] for e in entries]

        q_range = f"[{min(absQ):.2f}, {max(absQ):.2f}]"
        wp_range = f"[{min(absWp):.2f}, {max(absWp):.2f}]"
        c3_range = f"[{min(c3s)}, {max(c3s)}]"

        print(f"  {H:>5}  {len(entries):>5}  {q_range:>25}  {wp_range:>25}  {c3_range:>12}")

    # Phase 6: Check monotonicity of |Q(i)| vs H
    print(f"\n--- Phase 6: Monotonicity of |Q(i)| vs H ---")
    absQ_ranges = {}
    for H in H_set:
        entries = H_to_data[H]
        absQ = [abs(e['Q']) for e in entries]
        absQ_ranges[H] = (min(absQ), max(absQ))

    violations = 0
    for i, H1 in enumerate(H_set):
        for j, H2 in enumerate(H_set):
            if H1 < H2 and absQ_ranges[H2][1] > absQ_ranges[H1][0]:
                violations += 1
    print(f"  Anti-monotonicity violations (H1<H2 but max|Q(H2)| > min|Q(H1)|): {violations}")

    # Phase 7: Exact polynomial examples
    print(f"\n--- Phase 7: Example polynomials ---")
    for H in H_set[:3] + H_set[-3:]:
        entries = H_to_data[H]
        e = entries[0]
        print(f"\n  H={H}: P(u) coefficients = {e['poly']}")
        q, rem = poly_divide_by_u2_plus_1(e['poly'])
        print(f"    Q(u) coefficients = {q}")
        print(f"    Remainder = {rem}")
        print(f"    Q(i) = {e['Q']:.4f}")
        print(f"    |Q(i)| = {abs(e['Q']):.4f}")

    # Phase 8: Is (u^2+1)^2 also a factor?
    print(f"\n--- Phase 8: Does (u^2+1)^2 divide P(u)? ---")
    double_div = 0
    for bits in range(total):
        poly_u = poly_coeffs_list[bits]
        q, rem = poly_divide_by_u2_plus_1(poly_u)
        if all(abs(r) < 1e-10 for r in rem):
            q2, rem2 = poly_divide_by_u2_plus_1(q)
            if all(abs(r) < 1e-10 for r in rem2):
                double_div += 1
    print(f"  (u^2+1)^2 divides P(u): {double_div}/{total}")

    return H_values, Q_at_ihalf, Wprime_at_ihalf


def main():
    # Start with n=3 and n=5 (quick)
    for n in [3, 5]:
        analyze_n(n)

    # n=7 is 2^21 = 2M tournaments with 5040 permutations each — too slow for full perm enum
    # Use DP-based approach for n=7
    print(f"\n{'='*70}")
    print(f"n=7 ANALYSIS (using DP for forward-edge distribution)")
    print(f"{'='*70}")
    analyze_n7_dp()


def forward_edge_dist_dp(T):
    """DP-based forward edge distribution. O(2^n * n^2)."""
    n = len(T)
    deg = n - 1
    # dp[mask][v][f] = # paths through mask ending at v with f fwd edges
    # Flatten f into dp to save memory: dp[mask*n + v] = array of size n
    full = (1 << n) - 1

    # Use dict of (mask, v) -> list of f-counts for memory efficiency
    dp = {}
    for v in range(n):
        key = (1 << v, v)
        dp[key] = [0] * n
        dp[key][0] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp:
                continue
            fv = dp[key]
            for u in range(n):
                if mask & (1 << u):
                    continue
                new_mask = mask | (1 << u)
                new_key = (new_mask, u)
                if new_key not in dp:
                    dp[new_key] = [0] * n
                if T[v][u]:  # forward
                    for f in range(deg):
                        if fv[f]:
                            dp[new_key][f + 1] += fv[f]
                else:  # backward
                    for f in range(n):
                        if fv[f]:
                            dp[new_key][f] += fv[f]

    Nf = [0] * n
    for v in range(n):
        key = (full, v)
        if key in dp:
            for f in range(n):
                Nf[f] += dp[key][f]
    return Nf


def analyze_n7_dp():
    n = 7
    m = 21
    total = 1 << m
    deg = 6

    # Precompute coefficient pattern
    # For n=7: W(i/2)*2^6 = 8*(N_0 - N_2 + N_4 - N_6) = 0 always
    # So we need to compute the REDUCED polynomial

    # Actually, let me be smarter. Instead of computing the full Nf distribution
    # for each tournament, let me compute the polynomial coefficients directly.

    # W(r) = perm of matrix M where M[P_k][P_{k+1}] = r + T[P_k][P_{k+1}] - 1/2
    # Actually W(r) = sum_P prod_k (r + T[P_k][P_{k+1}] - 1/2)
    # = sum_P prod_k of (r+1/2) if forward, (r-1/2) if backward
    # In u=2r: prod_k of (u+1)/2 or (u-1)/2 = (1/2^6) * prod of (u+1) or (u-1)

    # 2^6 * W = sum_P prod_k of (u+1) or (u-1)
    # Each factor is (u+1) or (u-1). Product = sum of u^{6-j} * elem_symm_j terms.

    # For computational speed, compute polynomial W via DP.
    # dp[mask][v] = polynomial (as list of coefficients) for paths through mask ending at v.
    # Each step: multiply by (u+1) or (u-1).

    # Polynomial multiplication by (u+1): [a0, a1, ...] -> [a0, a0+a1, a1+a2, ...]
    # Polynomial multiplication by (u-1): [a0, a1, ...] -> [-a0, a0-a1, a1-a2, ...]

    # This is O(2^n * n * n) per tournament, with n-length polynomials.
    # 128 * 7 * 7 * 2M = too many... 128*49 = 6272 * 2M = 12.5B operations.
    # About 30 minutes in Python. Let me sample instead.

    print(f"\n  Sampling {min(total, 50000)} tournaments...")

    import random
    random.seed(42)

    sample_size = min(total, 50000)
    if sample_size < total:
        sample_bits = random.sample(range(total), sample_size)
    else:
        sample_bits = list(range(total))

    H_values = []
    Q_ihalf_values = []

    step = max(1, sample_size // 20)

    for idx, bits in enumerate(sample_bits):
        if idx % step == 0:
            print(f"    Progress: {idx}/{sample_size}", file=sys.stderr)

        T = tournament_from_bits(n, bits)

        # Compute polynomial via DP
        # dp[mask][v] = polynomial coefficients of 2^6 * W restricted to paths thru mask ending at v
        full = (1 << n) - 1

        # Use flat arrays for speed
        # poly[mask][v] is a list of 7 ints (coeffs of u^0..u^6)
        # Too much memory for dict approach. Use array.

        # dp[mask * n + v] = list of deg+1 coefficients
        dp = [None] * ((1 << n) * n)

        for v in range(n):
            dp[(1 << v) * n + v] = [0] * (deg + 1)
            dp[(1 << v) * n + v][0] = 1  # constant 1

        for mask in range(1, 1 << n):
            bc = bin(mask).count('1')
            if bc >= n:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pv = dp[mask * n + v]
                if pv is None:
                    continue

                for u in range(n):
                    if mask & (1 << u):
                        continue
                    new_mask = mask | (1 << u)
                    idx_new = new_mask * n + u
                    if dp[idx_new] is None:
                        dp[idx_new] = [0] * (deg + 1)

                    if T[v][u]:  # multiply by (u+1)
                        # new[k] += pv[k] + pv[k-1]
                        for k in range(bc, -1, -1):
                            if pv[k] == 0:
                                continue
                            dp[idx_new][k] += pv[k]      # constant * pv[k]
                            if k + 1 <= deg:
                                dp[idx_new][k + 1] += pv[k]  # u * pv[k]
                    else:  # multiply by (u-1)
                        for k in range(bc, -1, -1):
                            if pv[k] == 0:
                                continue
                            dp[idx_new][k] -= pv[k]      # -1 * pv[k]
                            if k + 1 <= deg:
                                dp[idx_new][k + 1] += pv[k]  # u * pv[k]

        # Sum over all ending vertices
        poly_u = [0] * (deg + 1)
        for v in range(n):
            pv = dp[full * n + v]
            if pv is not None:
                for k in range(deg + 1):
                    poly_u[k] += pv[k]

        # H = W(1/2) = P(1) / 2^6 where P(u) = poly evaluated at u=1
        P_at_1 = sum(poly_u)
        H = P_at_1 // (2 ** deg)  # should be exact

        # Divide P(u) by (u^2 + 1) to get Q(u)
        q, rem = poly_divide_by_u2_plus_1(poly_u)

        # Evaluate Q at u=i
        Q_i = eval_poly(q, 1j)

        # Also: W'(i/2) = P'(i) / 2^{deg-1}
        dpoly = poly_derivative(poly_u)
        Wp_i = eval_poly(dpoly, 1j)

        H_values.append(H)
        Q_ihalf_values.append(Q_i)

    N = len(H_values)
    print(f"\n  Processed {N} tournaments")

    # Analyze Q(i)
    Q_real = [q.real for q in Q_ihalf_values]
    Q_imag = [q.imag for q in Q_ihalf_values]

    print(f"\n  Q(i) real parts range: [{min(Q_real):.4f}, {max(Q_real):.4f}]")
    print(f"  Q(i) imag parts range: [{min(Q_imag):.4f}, {max(Q_imag):.4f}]")

    all_real_Q = all(abs(q) < 0.5 for q in Q_imag)
    all_imag_Q = all(abs(q) < 0.5 for q in Q_real)
    print(f"  Q(i) purely real: {all_real_Q}")
    print(f"  Q(i) purely imaginary: {all_imag_Q}")

    # Check if Q(i) is always an integer or half-integer
    print(f"\n  Q(i) integrality check:")
    Q_int_real = all(abs(q - round(q)) < 0.01 for q in Q_real)
    Q_int_imag = all(abs(q - round(q)) < 0.01 for q in Q_imag)
    print(f"    Real part always integer: {Q_int_real}")
    print(f"    Imag part always integer: {Q_int_imag}")

    # Correlations
    mean_H = sum(H_values) / N
    absQ = [abs(q) for q in Q_ihalf_values]
    mean_absQ = sum(absQ) / N

    cov = sum((H_values[i] - mean_H) * (absQ[i] - mean_absQ) for i in range(N)) / N
    var_H = sum((H_values[i] - mean_H) ** 2 for i in range(N)) / N
    var_absQ = sum((absQ[i] - mean_absQ) ** 2 for i in range(N)) / N

    corr = cov / (var_H ** 0.5 * var_absQ ** 0.5) if var_H > 0 and var_absQ > 0 else 0
    print(f"\n  Corr(H, |Q(i)|) = {corr:.6f}")

    # Signed correlation
    if all_real_Q:
        mean_Qr = sum(Q_real) / N
        cov_r = sum((H_values[i] - mean_H) * (Q_real[i] - mean_Qr) for i in range(N)) / N
        var_Qr = sum((Q_real[i] - mean_Qr) ** 2 for i in range(N)) / N
        corr_r = cov_r / (var_H ** 0.5 * var_Qr ** 0.5) if var_H > 0 and var_Qr > 0 else 0
        print(f"  Corr(H, Re Q(i)) = {corr_r:.6f}")

    # Table by H
    H_to_Q = defaultdict(list)
    for i in range(N):
        H_to_Q[H_values[i]].append(Q_ihalf_values[i])

    H_set = sorted(H_to_Q.keys())

    print(f"\n  {'H':>5}  {'cnt':>5}  {'|Q(i)| min':>12}  {'|Q(i)| max':>12}  {'|Q(i)| mean':>12}  {'Re Q rng':>25}")
    for H in H_set:
        entries = H_to_Q[H]
        absQ_h = [abs(q) for q in entries]
        reQ = [q.real for q in entries]
        imQ = [q.imag for q in entries]
        mn, mx = min(absQ_h), max(absQ_h)
        mean = sum(absQ_h) / len(absQ_h)
        re_rng = f"[{min(reQ):.1f}, {max(reQ):.1f}]"
        print(f"  {H:>5}  {len(entries):>5}  {mn:>12.2f}  {mx:>12.2f}  {mean:>12.2f}  {re_rng:>25}")

    # Monotonicity
    absQ_ranges = {}
    for H in H_set:
        entries = H_to_Q[H]
        absQ_h = [abs(q) for q in entries]
        absQ_ranges[H] = (min(absQ_h), max(absQ_h))

    violations = 0
    for i, H1 in enumerate(H_set):
        for j, H2 in enumerate(H_set):
            if H1 < H2 and absQ_ranges[H2][1] > absQ_ranges[H1][0]:
                violations += 1
    print(f"\n  Anti-monotonicity violations for |Q(i)|: {violations}")


if __name__ == "__main__":
    main()
