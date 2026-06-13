#!/usr/bin/env python3
"""
W(r) at COMPLEX values of r -- investigation of complex evaluations.

W(r) = sum_P prod_{j=0}^{n-2} (r + A[P(j),P(j+1)] - 1/2)

Key known results:
  W(1/2) = H(T) = number of Hamiltonian paths
  W(-1/2) = (-1)^{n-1} * H(T)
  W has r-parity at odd n: W(-r) = (-1)^{n-1} W(r)
  W(r) = sum_k w_k r^k with only even-k terms at odd n

Investigations:
  1. W(i/2) and |W(i/2)|^2
  2. W at scaled roots of unity r = e^{2pi*i*k/m}/2
  3. Parseval integral of |W|^2 on circle |r|=1/2
  4. Resultant of W with W' (discriminant)
  5. W(1) = (1/2)^{n-1} * sum_P 3^{fwd(P)}

kind-pasteur-2026-03-07
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict
import random

# ========== Tournament construction ==========

def tournament_from_tiling(n, tiling_bits):
    """Construct tournament adjacency matrix from tiling bits."""
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def random_tournament(n, seed=42):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def paley_tournament(p):
    """Construct Paley tournament T_p for prime p = 3 mod 4."""
    qr = set()
    for x in range(1, p):
        qr.add((x*x) % p)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

# ========== W polynomial computation ==========

def compute_W_coeffs_exact(A, n):
    """Compute W(r) coefficients using DP. Returns list [w_0, w_1, ..., w_{n-1}]."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = [1.0]

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            poly = dp.get((mask, v))
            if poly is None:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                s = A[v][u] - 0.5
                new_poly_len = len(poly) + 1
                key = (mask | (1 << u), u)
                if key not in dp:
                    dp[key] = [0.0] * new_poly_len
                existing = dp[key]
                while len(existing) < new_poly_len:
                    existing.append(0.0)
                for i in range(len(poly)):
                    existing[i] += s * poly[i]
                    existing[i+1] += poly[i]

    full = (1 << n) - 1
    result = [0.0] * n
    for v in range(n):
        poly = dp.get((full, v), [])
        for i in range(min(len(poly), n)):
            result[i] += poly[i]
    return result

def poly_eval_complex(coeffs, r):
    """Evaluate polynomial at complex r. coeffs[k] = coefficient of r^k."""
    val = complex(0)
    for k, c in enumerate(coeffs):
        val += c * (r ** k)
    return val

def compute_forward_edge_dist(A, n):
    """Compute forward edge distribution: a_k = # perms with k forward edges."""
    dist = defaultdict(int)
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        dist[fwd] += 1
    return dict(dist)

# ========== Resultant / Discriminant ==========

def poly_derivative(coeffs):
    """Derivative of polynomial."""
    return [k * coeffs[k] for k in range(1, len(coeffs))]

def resultant_numpy(p, q):
    """Compute resultant of two polynomials via Sylvester matrix."""
    m = len(p) - 1  # degree of p
    n_q = len(q) - 1  # degree of q
    if m < 0 or n_q < 0:
        return 0.0
    size = m + n_q
    if size == 0:
        return 1.0
    S = np.zeros((size, size))
    # Rows for p (n_q rows)
    for i in range(n_q):
        for j in range(len(p)):
            S[i, i + j] = p[len(p) - 1 - j]
    # Rows for q (m rows)
    for i in range(m):
        for j in range(len(q)):
            S[n_q + i, i + j] = q[len(q) - 1 - j]
    return np.linalg.det(S)

# ========== Main investigation ==========

def investigate_tournament(A, n, label=""):
    """Full complex evaluation investigation for one tournament."""
    coeffs = compute_W_coeffs_exact(A, n)
    t3 = count_3cycles(A, n)
    H = int(round(poly_eval_complex(coeffs, 0.5).real))

    print(f"\n{'='*70}")
    print(f"Tournament: {label}  n={n}  t3={t3}  H={H}")
    print(f"{'='*70}")
    print(f"W(r) coefficients: {['%.6g' % c for c in coeffs]}")

    # ---- 1. W(i/2) ----
    print(f"\n--- 1. W(i/2) where i = sqrt(-1) ---")
    i_half = 0.5j
    W_ihalf = poly_eval_complex(coeffs, i_half)
    print(f"  W(i/2) = {W_ihalf.real:.10f} + {W_ihalf.imag:.10f}i")
    print(f"  |W(i/2)|^2 = {abs(W_ihalf)**2:.10f}")

    # For odd n, W has only even-power terms: w_0, w_2, w_4, ...
    # At r=i/2: (i/2)^k = i^k / 2^k
    # Even k: i^k alternates +1, -1 (real); odd k: alternates +i, -i (imaginary)
    # Since at odd n only even-k terms survive:
    #   W(i/2) = sum_{k even} w_k * (i/2)^k = sum_{j} w_{2j} * (-1)^j / 2^{2j}
    # This is PURELY REAL at odd n!
    if n % 2 == 1:
        print(f"  [Odd n: only even-power terms => W(i/2) is REAL]")
        real_check = 0.0
        for j in range(len(coeffs) // 2 + 1):
            k = 2 * j
            if k < len(coeffs):
                real_check += coeffs[k] * ((-1)**j) / (4**j)
        print(f"  Direct formula: sum w_{{2j}} * (-1)^j / 4^j = {real_check:.10f}")
        print(f"  Match: {abs(W_ihalf.real - real_check) < 1e-8}")

    # W(-i/2)
    W_mihalf = poly_eval_complex(coeffs, -0.5j)
    print(f"  W(-i/2) = {W_mihalf.real:.10f} + {W_mihalf.imag:.10f}i")
    if n % 2 == 1:
        print(f"  W(-i/2) should equal (-1)^{{n-1}} W(i/2) = W(i/2) = {W_ihalf.real:.10f}")
        print(f"  Match: {abs(W_mihalf - W_ihalf) < 1e-8}")

    # ---- 2. W at scaled roots of unity ----
    print(f"\n--- 2. W(r) at r = e^{{2*pi*i*k/m}} / 2 ---")
    for m_val in [3, 4, 5, 6, 8]:
        print(f"  m={m_val}:")
        vals = []
        for k in range(m_val):
            r = np.exp(2j * np.pi * k / m_val) / 2
            w = poly_eval_complex(coeffs, r)
            vals.append(w)
            print(f"    k={k}: r={r:.4f}, W={w.real:.6f} + {w.imag:.6f}i, |W|={abs(w):.6f}")
        # Sum of W at all m-th roots of unity (scaled by 1/2)
        total = sum(vals)
        print(f"    Sum = {total.real:.6f} + {total.imag:.6f}i")
        prod_abs = 1.0
        for v in vals:
            prod_abs *= abs(v)
        print(f"    Product of |W| = {prod_abs:.6f}")

    # ---- 3. Parseval: integral of |W|^2 on |r|=1/2 ----
    print(f"\n--- 3. Parseval on |r|=1/2 ---")
    # |W(r)|^2 integrated over the circle |r|=1/2 with r = (1/2)e^{i*theta}:
    # (1/2pi) integral_0^{2pi} |W((1/2)e^{i*theta})|^2 d(theta)
    # = sum_k |w_k|^2 / 4^k  (by Parseval / orthogonality of e^{ik*theta})
    parseval_sum = sum(abs(coeffs[k])**2 / 4**k for k in range(len(coeffs)))
    print(f"  Parseval sum = sum_k |w_k|^2 / 4^k = {parseval_sum:.10f}")

    # Numerical check via quadrature
    N_quad = 10000
    thetas = np.linspace(0, 2*np.pi, N_quad, endpoint=False)
    quad_sum = 0.0
    for theta in thetas:
        r = 0.5 * np.exp(1j * theta)
        w = poly_eval_complex(coeffs, r)
        quad_sum += abs(w)**2
    quad_avg = quad_sum / N_quad
    print(f"  Numerical quadrature = {quad_avg:.10f}")
    print(f"  Match: {abs(parseval_sum - quad_avg) < 0.01}")

    # Compare to H^2 and other invariants
    print(f"  H^2 = {H**2}")
    print(f"  Parseval / H^2 = {parseval_sum / H**2:.10f}" if H > 0 else "  H=0")

    # Individual |w_k|^2 / 4^k terms
    print(f"  Breakdown: |w_k|^2 / 4^k:")
    for k in range(len(coeffs)):
        val = abs(coeffs[k])**2 / 4**k
        if abs(val) > 1e-12:
            print(f"    k={k}: |w_{k}|^2 = {coeffs[k]**2:.4f}, /4^{k} = {val:.6f}")

    # ---- 4. Discriminant = Res(W, W') / leading_coeff ----
    print(f"\n--- 4. Resultant of W with W' ---")
    W_prime = poly_derivative(coeffs)
    if len(W_prime) > 0 and any(abs(c) > 1e-12 for c in W_prime):
        # Use numpy for resultant
        res = resultant_numpy(coeffs, W_prime)
        print(f"  Res(W, W') = {res:.6f}")
        # Discriminant is Res(W, W') / ((-1)^{d(d-1)/2} * a_d) where d = degree, a_d = leading
        deg = len(coeffs) - 1
        while deg > 0 and abs(coeffs[deg]) < 1e-12:
            deg -= 1
        if deg > 0 and abs(coeffs[deg]) > 1e-12:
            sign = (-1) ** (deg * (deg - 1) // 2)
            disc = res / (sign * coeffs[deg])
            print(f"  Disc(W) = {disc:.6f}")
            print(f"  (degree={deg}, leading coeff={coeffs[deg]:.4f})")
    else:
        print(f"  W' is zero or trivial")

    # ---- 5. W(1) = (1/2)^{n-1} * sum_P 3^{fwd(P)} ----
    print(f"\n--- 5. W(1) weighted count ---")
    W1 = poly_eval_complex(coeffs, 1.0).real
    print(f"  W(1) = {W1:.6f}")
    print(f"  W(1) * 2^{{n-1}} = {W1 * 2**(n-1):.6f}  (= sum_P 3^{{fwd}})")

    if n <= 7:
        fwd_dist = compute_forward_edge_dist(A, n)
        weighted_sum = sum(count * 3**fwd for fwd, count in fwd_dist.items())
        print(f"  Direct sum_P 3^{{fwd}} = {weighted_sum}")
        print(f"  Forward edge distribution: {dict(sorted(fwd_dist.items()))}")
        print(f"  Check: W(1) * 2^{{n-1}} = {W1 * 2**(n-1):.2f} vs {weighted_sum}")

    # W(-1) for completeness
    Wm1 = poly_eval_complex(coeffs, -1.0).real
    print(f"  W(-1) = {Wm1:.6f}")
    if n % 2 == 1:
        print(f"  W(-1) = (-1)^{{n-1}} W(1) = {W1:.6f}  Match: {abs(Wm1 - W1) < 1e-8}")

    # ---- Bonus: W on imaginary axis ----
    print(f"\n--- Bonus: |W(it)|^2 for real t ---")
    for t in [0.25, 0.5, 1.0, 2.0]:
        w = poly_eval_complex(coeffs, 1j * t)
        print(f"  W({t}i) = {w.real:.6f} + {w.imag:.6f}i, |W|^2 = {abs(w)**2:.6f}")

    return coeffs, t3, H


def run_all_n3():
    """Exhaustive analysis at n=3."""
    n = 3
    m = num_tiling_bits(n)
    print(f"\n{'#'*70}")
    print(f"# EXHAUSTIVE ANALYSIS: n={n} ({2**m} tournaments)")
    print(f"{'#'*70}")

    all_data = []
    for bits in range(2**m):
        A = tournament_from_tiling(n, bits)
        coeffs, t3, H = investigate_tournament(A, n, f"bits={bits}")
        all_data.append((coeffs, t3, H, bits))

    print(f"\n{'='*70}")
    print(f"SUMMARY n=3: |W(i/2)|^2 values:")
    for coeffs, t3, H, bits in all_data:
        w = poly_eval_complex(coeffs, 0.5j)
        print(f"  bits={bits}: t3={t3}, H={H}, W(i/2)={w.real:.6f}, |W(i/2)|^2={abs(w)**2:.6f}")


def run_all_n5():
    """Exhaustive analysis at n=5."""
    n = 5
    m = num_tiling_bits(n)
    print(f"\n{'#'*70}")
    print(f"# EXHAUSTIVE ANALYSIS: n={n} ({2**m} tournaments)")
    print(f"{'#'*70}")

    # Collect statistics
    ihalf_sq_vals = defaultdict(list)
    parseval_vals = defaultdict(list)
    W1_vals = defaultdict(list)

    for bits in range(2**m):
        A = tournament_from_tiling(n, bits)
        coeffs = compute_W_coeffs_exact(A, n)
        t3 = count_3cycles(A, n)
        H = int(round(poly_eval_complex(coeffs, 0.5).real))

        w_ihalf = poly_eval_complex(coeffs, 0.5j)
        ihalf_sq = abs(w_ihalf)**2
        parseval = sum(abs(coeffs[k])**2 / 4**k for k in range(len(coeffs)))
        W1 = poly_eval_complex(coeffs, 1.0).real

        ihalf_sq_vals[(t3, H)].append(ihalf_sq)
        parseval_vals[(t3, H)].append(parseval)
        W1_vals[(t3, H)].append(W1)

    print(f"\n{'='*70}")
    print(f"n=5: |W(i/2)|^2 by (t3, H)")
    print(f"{'='*70}")
    for key in sorted(ihalf_sq_vals.keys()):
        vals = ihalf_sq_vals[key]
        unique_vals = sorted(set(round(v, 6) for v in vals))
        print(f"  (t3={key[0]}, H={key[1]}): |W(i/2)|^2 in {unique_vals}  (count={len(vals)})")

    print(f"\n{'='*70}")
    print(f"n=5: Parseval sum by (t3, H)")
    print(f"{'='*70}")
    for key in sorted(parseval_vals.keys()):
        vals = parseval_vals[key]
        unique_vals = sorted(set(round(v, 6) for v in vals))
        print(f"  (t3={key[0]}, H={key[1]}): Parseval in {unique_vals}  (count={len(vals)})")

    print(f"\n{'='*70}")
    print(f"n=5: W(1) by (t3, H)")
    print(f"{'='*70}")
    for key in sorted(W1_vals.keys()):
        vals = W1_vals[key]
        unique_vals = sorted(set(round(v, 6) for v in vals))
        print(f"  (t3={key[0]}, H={key[1]}): W(1) in {unique_vals}  (count={len(vals)})")

    # Check if |W(i/2)|^2 = f(H, t3) for some function f
    print(f"\n{'='*70}")
    print(f"n=5: Is |W(i/2)|^2 a function of (t3, H)?")
    print(f"{'='*70}")
    all_determined = True
    for key in sorted(ihalf_sq_vals.keys()):
        vals = ihalf_sq_vals[key]
        unique = set(round(v, 4) for v in vals)
        if len(unique) > 1:
            print(f"  NO: (t3={key[0]}, H={key[1]}) has multiple values: {sorted(unique)}")
            all_determined = False
    if all_determined:
        print(f"  YES! |W(i/2)|^2 is determined by (t3, H) at n=5")

    # Show a few detailed examples
    for bits in [0, 15, 31, 63]:
        if bits < 2**m:
            A = tournament_from_tiling(n, bits)
            investigate_tournament(A, n, f"bits={bits}")


def run_sample_n7():
    """Sampled analysis at n=7."""
    n = 7
    print(f"\n{'#'*70}")
    print(f"# SAMPLED ANALYSIS: n={n}")
    print(f"{'#'*70}")

    # Paley tournament T_7
    A_paley = paley_tournament(7)
    investigate_tournament(A_paley, n, "Paley T_7")

    # Transitive tournament
    A_trans = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A_trans[i][j] = 1
    investigate_tournament(A_trans, n, "Transitive")

    # Random samples
    ihalf_data = []
    parseval_data = []
    W1_data = []

    N_samples = 200
    for seed in range(N_samples):
        A = random_tournament(n, seed)
        coeffs = compute_W_coeffs_exact(A, n)
        t3 = count_3cycles(A, n)
        H = int(round(poly_eval_complex(coeffs, 0.5).real))

        w_ihalf = poly_eval_complex(coeffs, 0.5j)
        parseval = sum(abs(coeffs[k])**2 / 4**k for k in range(len(coeffs)))
        W1 = poly_eval_complex(coeffs, 1.0).real

        ihalf_data.append((t3, H, abs(w_ihalf)**2, w_ihalf.real))
        parseval_data.append((t3, H, parseval))
        W1_data.append((t3, H, W1))

    print(f"\n{'='*70}")
    print(f"n=7 SAMPLED: |W(i/2)|^2 statistics ({N_samples} samples)")
    print(f"{'='*70}")

    # Group by H
    by_H = defaultdict(list)
    for t3, H, wsq, wreal in ihalf_data:
        by_H[H].append((t3, wsq, wreal))
    for H_val in sorted(by_H.keys()):
        entries = by_H[H_val]
        wsq_vals = [e[1] for e in entries]
        print(f"  H={H_val}: |W(i/2)|^2 in [{min(wsq_vals):.4f}, {max(wsq_vals):.4f}] "
              f"(n={len(entries)}, mean={np.mean(wsq_vals):.4f})")

    # Check formula: is |W(i/2)|^2 related to H^2 + something?
    print(f"\n{'='*70}")
    print(f"n=7: Testing |W(i/2)|^2 = H^2 hypothesis")
    print(f"{'='*70}")
    for t3, H, wsq, wreal in ihalf_data[:20]:
        ratio = wsq / H**2 if H > 0 else float('inf')
        diff = wsq - H**2
        print(f"  t3={t3}, H={H}: |W(i/2)|^2={wsq:.4f}, H^2={H**2}, "
              f"ratio={ratio:.6f}, diff={diff:.4f}")

    # At odd n, W(i/2) is real, so |W(i/2)|^2 = W(i/2)^2
    print(f"\n{'='*70}")
    print(f"n=7: W(i/2) is real at odd n -- check sign pattern")
    print(f"{'='*70}")
    pos_count = 0
    neg_count = 0
    for t3, H, wsq, wreal in ihalf_data:
        if wreal > 0:
            pos_count += 1
        elif wreal < 0:
            neg_count += 1
    print(f"  W(i/2) > 0: {pos_count} times")
    print(f"  W(i/2) < 0: {neg_count} times")

    # W(i/2) formula: sum_{j} w_{2j} * (-1)^j / 4^j
    # = w_0 - w_2/4 + w_4/16 - w_6/64 (at n=7, max degree 6)
    print(f"\n{'='*70}")
    print(f"n=7: W(i/2) = w_0 - w_2/4 + w_4/16 - w_6/64")
    print(f"{'='*70}")
    for seed in range(5):
        A = random_tournament(n, seed)
        coeffs = compute_W_coeffs_exact(A, n)
        H = int(round(poly_eval_complex(coeffs, 0.5).real))
        t3 = count_3cycles(A, n)
        w_ihalf = coeffs[0] - coeffs[2]/4 + coeffs[4]/16 - coeffs[6]/64
        print(f"  seed={seed}: w=[{coeffs[0]:.1f}, {coeffs[2]:.1f}, {coeffs[4]:.1f}, {coeffs[6]:.1f}]"
              f"  W(i/2)={w_ihalf:.4f}  H={H}  t3={t3}")

    # Parseval statistics
    print(f"\n{'='*70}")
    print(f"n=7 SAMPLED: Parseval sum statistics")
    print(f"{'='*70}")
    by_H_parseval = defaultdict(list)
    for t3, H, pval in parseval_data:
        by_H_parseval[H].append(pval)
    for H_val in sorted(by_H_parseval.keys()):
        pvals = by_H_parseval[H_val]
        print(f"  H={H_val}: Parseval in [{min(pvals):.4f}, {max(pvals):.4f}] "
              f"(n={len(pvals)}, mean={np.mean(pvals):.4f})")

    # W(1) statistics
    print(f"\n{'='*70}")
    print(f"n=7 SAMPLED: W(1) = (1/2)^6 * sum 3^fwd statistics")
    print(f"{'='*70}")
    by_H_W1 = defaultdict(list)
    for t3, H, w1 in W1_data:
        by_H_W1[H].append(w1)
    for H_val in sorted(by_H_W1.keys()):
        w1s = by_H_W1[H_val]
        w1_int = [round(v * 64) for v in w1s]  # W(1) * 2^6
        print(f"  H={H_val}: W(1)*64 in [{min(w1_int)}, {max(w1_int)}] "
              f"(n={len(w1s)})")


def cross_n_patterns():
    """Look for patterns across n=3,5,7."""
    print(f"\n{'#'*70}")
    print(f"# CROSS-n PATTERNS")
    print(f"{'#'*70}")

    for n in [3, 5]:
        m = num_tiling_bits(n)
        print(f"\n--- n={n} exhaustive ---")
        # Find: is there a formula for W(i/2) in terms of H and t3?
        data = []
        for bits in range(2**m):
            A = tournament_from_tiling(n, bits)
            coeffs = compute_W_coeffs_exact(A, n)
            t3 = count_3cycles(A, n)
            H = int(round(poly_eval_complex(coeffs, 0.5).real))
            w_ihalf = poly_eval_complex(coeffs, 0.5j).real
            parseval = sum(abs(coeffs[k])**2 / 4**k for k in range(len(coeffs)))
            W1 = poly_eval_complex(coeffs, 1.0).real
            data.append((t3, H, w_ihalf, parseval, W1, coeffs))

        # Try: W(i/2) = a*H + b*t3 + c
        print(f"  Linear regression: W(i/2) ~ a*H + b*t3 + c")
        X = np.array([[d[1], d[0], 1.0] for d in data])
        y = np.array([d[2] for d in data])
        try:
            result = np.linalg.lstsq(X, y, rcond=None)
            abc = result[0]
            residuals = y - X @ abc
            max_res = np.max(np.abs(residuals))
            print(f"    a={abc[0]:.6f}, b={abc[1]:.6f}, c={abc[2]:.6f}")
            print(f"    Max residual: {max_res:.10f}")
            if max_res < 1e-8:
                print(f"    EXACT FIT! W(i/2) = {abc[0]:.6f}*H + {abc[1]:.6f}*t3 + {abc[2]:.6f}")
        except Exception as e:
            print(f"    Regression failed: {e}")

        # Try: Parseval ~ a*H^2 + b*H*t3 + c*t3^2 + d*H + e*t3 + f
        print(f"  Quadratic regression: Parseval ~ a*H^2 + b*H*t3 + c*t3^2 + d*H + e*t3 + f")
        X2 = np.array([[d[1]**2, d[1]*d[0], d[0]**2, d[1], d[0], 1.0] for d in data])
        y2 = np.array([d[3] for d in data])
        try:
            result2 = np.linalg.lstsq(X2, y2, rcond=None)
            abc2 = result2[0]
            residuals2 = y2 - X2 @ abc2
            max_res2 = np.max(np.abs(residuals2))
            print(f"    Coefficients: {['%.6f' % x for x in abc2]}")
            print(f"    Max residual: {max_res2:.10f}")
            if max_res2 < 1e-6:
                print(f"    NEAR-EXACT FIT!")
        except Exception as e:
            print(f"    Regression failed: {e}")

        # Try: W(1)*2^{n-1} ~ a*H + b*t3 + c
        print(f"  Linear regression: W(1)*2^{{n-1}} ~ a*H + b*t3 + c")
        y3 = np.array([d[4] * 2**(n-1) for d in data])
        try:
            result3 = np.linalg.lstsq(X, y3, rcond=None)
            abc3 = result3[0]
            residuals3 = y3 - X @ abc3
            max_res3 = np.max(np.abs(residuals3))
            print(f"    a={abc3[0]:.6f}, b={abc3[1]:.6f}, c={abc3[2]:.6f}")
            print(f"    Max residual: {max_res3:.10f}")
        except Exception as e:
            print(f"    Regression failed: {e}")


def circle_profile():
    """Plot |W(r)|^2 on the circle |r|=1/2 for a few tournaments."""
    print(f"\n{'#'*70}")
    print(f"# |W(r)|^2 PROFILE on |r|=1/2")
    print(f"{'#'*70}")

    n = 5
    m = num_tiling_bits(n)
    # Pick a few contrasting tournaments
    test_cases = [
        (0, "Transitive"),
        (2**m - 1, "All-flip"),
    ]
    # Add Paley-like (regular if exists)
    # At n=5 regular: t3=2 => C(5,3)-sum d_i*(d_i-1)/2 = 10 - 2*C(2,2) - 3*C(1,2)
    # Regular n=5 has scores (2,2,2,2,2), t3 = C(5,3) - 5*C(2,2)/1 = 10 - 5 = 5? No.
    # t3 = C(5,3) - sum C(s_i, 2) = 10 - 5*1 = 5

    for bits, label in test_cases:
        A = tournament_from_tiling(n, bits)
        coeffs = compute_W_coeffs_exact(A, n)
        H = int(round(poly_eval_complex(coeffs, 0.5).real))
        t3 = count_3cycles(A, n)

        print(f"\n  {label} (bits={bits}, H={H}, t3={t3}):")
        thetas = np.linspace(0, 2*np.pi, 36, endpoint=False)
        for theta in thetas[::4]:  # Every 40 degrees
            r = 0.5 * np.exp(1j * theta)
            w = poly_eval_complex(coeffs, r)
            print(f"    theta={theta*180/np.pi:6.1f} deg: |W|^2={abs(w)**2:.6f}  "
                  f"Re={w.real:.6f}  Im={w.imag:.6f}")

        # Min and max of |W|^2 on the circle
        fine_thetas = np.linspace(0, 2*np.pi, 3600, endpoint=False)
        wsq_vals = [abs(poly_eval_complex(coeffs, 0.5 * np.exp(1j * t)))**2 for t in fine_thetas]
        print(f"    min|W|^2 = {min(wsq_vals):.6f} at theta={fine_thetas[np.argmin(wsq_vals)]*180/np.pi:.1f} deg")
        print(f"    max|W|^2 = {max(wsq_vals):.6f} at theta={fine_thetas[np.argmax(wsq_vals)]*180/np.pi:.1f} deg")
        print(f"    avg|W|^2 = {np.mean(wsq_vals):.6f} (= Parseval sum)")


if __name__ == "__main__":
    print("W(r) AT COMPLEX VALUES -- INVESTIGATION")
    print("=" * 70)

    # 1. Exhaustive n=3
    run_all_n3()

    # 2. Exhaustive n=5
    run_all_n5()

    # 3. Sampled n=7
    run_sample_n7()

    # 4. Cross-n pattern search
    cross_n_patterns()

    # 5. Circle profile
    circle_profile()

    print(f"\n{'#'*70}")
    print(f"# DONE")
    print(f"{'#'*70}")
