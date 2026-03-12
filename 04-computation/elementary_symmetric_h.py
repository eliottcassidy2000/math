#!/usr/bin/env python3
"""
Express H in the elementary symmetric polynomial basis of y².

KEY INSIGHT: Each e_k(x_1,...,x_m) is Schur-CONCAVE for k >= 1.
If H = c₀ + c₂·e₂ + c₃·e₃ + ... with ALL c_k >= 0,
then H is Schur-CONCAVE → maximized at uniform (Paley).

Also: Πx_i = e_m(x) is Schur-concave, and Spearman(H, Πy²) = 1.0.
This suggests the e_m term dominates.

Additional approach: the MONOMIAL symmetric functions m_λ.
If H = Σ a_λ · m_λ(y²) with a_λ >= 0, same conclusion.
"""

import numpy as np
from itertools import combinations
from scipy.optimize import linprog

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p-1)//2, p) == 1

def get_circulants(p):
    n = (p-1)//2
    results = []
    for S in combinations(range(1,p), n):
        S_set = set(S)
        if all((p-j)%p not in S_set for j in S_set):
            results.append(frozenset(S_set))
    return results

def adjacency_matrix(S, p):
    A = np.zeros((p,p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i!=j and (j-i)%p in S:
                A[i][j] = 1
    return A

def eigenvalues_circulant(S, p):
    omega = np.exp(2j*np.pi/p)
    return [sum(omega**(k*s) for s in S) for k in range(p)]

def count_hp_dp(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1<<v, v)] = 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)): continue
            if (mask,v) not in dp: continue
            c = dp[(mask,v)]
            for u in range(n):
                if mask & (1<<u): continue
                if A[v][u]:
                    key = (mask|(1<<u), u)
                    dp[key] = dp.get(key,0) + c
    return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

def analyze(p):
    print(f"\n{'='*70}")
    print(f"ELEMENTARY SYMMETRIC POLYNOMIAL ANALYSIS — p = {p}")
    print(f"{'='*70}")

    circulants = get_circulants(p)
    m = (p-1)//2
    paley_S = frozenset(j for j in range(1,p) if is_qr(j,p))

    data = []
    for S in circulants:
        A = adjacency_matrix(S, p)
        eigs = eigenvalues_circulant(S, p)
        H = count_hp_dp(A)
        y_half = [eigs[k].imag for k in range(1, m+1)]
        x = [y**2 for y in y_half]
        data.append({'S': S, 'H': H, 'x': x, 'paley': S == paley_S})

    by_H = {}
    for d in data:
        by_H.setdefault(d['H'], []).append(d)
    reps = [by_H[H][0] for H in sorted(by_H, reverse=True)]

    n_reps = len(reps)
    print(f"{len(circulants)} tournaments, {n_reps} distinct H values, m={m}")
    print(f"Paley: {sorted(paley_S)}")

    # Compute ALL elementary symmetric polynomials of y²
    for d in data:
        x = d['x']
        d['e'] = {0: 1.0}
        for k in range(1, m+1):
            d['e'][k] = sum(np.prod(list(combo)) for combo in combinations(x, k))

    # Show e_k values for representatives
    print(f"\n  Elementary symmetric polynomials e_k(y²) for k=1,...,{m}:")
    header = f"{'H':>10}"
    for k in range(1, m+1):
        header += f"  {'e_'+str(k):>12}"
    print(header)

    for d in reps:
        line = f"{d['H']:>10}"
        for k in range(1, m+1):
            line += f"  {d['e'][k]:>12.4f}"
        line += "  P" if d['paley'] else ""
        print(line)

    # Check: is e_1 universal?
    e1_vals = set(round(d['e'][1], 4) for d in data)
    print(f"\n  e_1 universal? {len(e1_vals) == 1} (value: {e1_vals})")

    # Fit H = c₀ + c₂·e₂ + c₃·e₃ + ... + c_m·e_m
    # Skip e_1 (universal)
    print(f"\n--- Fitting H = c₀ + Σ_{k>=2} c_k · e_k(y²) ---")

    # n_reps equations, need at most n_reps unknowns
    # Use e_2,...,e_{min(m, n_reps-1)+1}
    max_k = min(m, n_reps - 1) + 1  # +1 to include up to this index

    for n_terms in range(2, max_k + 1):
        # Use e_2,...,e_{n_terms}
        n_vars = n_terms  # constant + (n_terms-1) elementary polys
        if n_vars > n_reps:
            break

        A_mat = np.zeros((n_reps, n_vars))
        b_vec = np.array([d['H'] for d in reps], dtype=float)

        for i, d in enumerate(reps):
            A_mat[i, 0] = 1  # constant
            for j in range(1, n_vars):
                A_mat[i, j] = d['e'][j + 1]  # e_{j+1}

        rank = np.linalg.matrix_rank(A_mat)

        if rank < n_vars:
            print(f"\n  Using e_2,...,e_{n_terms}: rank {rank} < {n_vars}, skip")
            continue

        if n_vars == n_reps:
            try:
                coeffs = np.linalg.solve(A_mat, b_vec)
                residual = np.max(np.abs(A_mat @ coeffs - b_vec))
            except np.linalg.LinAlgError:
                print(f"\n  Using e_2,...,e_{n_terms}: singular")
                continue
        else:
            coeffs = np.linalg.lstsq(A_mat, b_vec, rcond=None)[0]
            residual = np.max(np.abs(A_mat @ coeffs - b_vec))

        print(f"\n  Using e_2,...,e_{n_terms} ({n_vars} params, residual={residual:.2e}):")
        print(f"    H = {coeffs[0]:.4f}", end="")
        for j in range(1, n_vars):
            sign = '+' if coeffs[j] >= 0 else ''
            print(f" {sign}{coeffs[j]:.6f}·e_{j+1}", end="")
        print()

        if residual < 1.0:
            pos_count = sum(1 for j in range(1, n_vars) if coeffs[j] > 0)
            neg_count = sum(1 for j in range(1, n_vars) if coeffs[j] < 0)
            zero_count = sum(1 for j in range(1, n_vars) if abs(coeffs[j]) < 1e-10)
            print(f"    Signs: {pos_count} positive, {neg_count} negative, {zero_count} zero")

            if neg_count == 0:
                print(f"    ★★★ ALL e_k COEFFICIENTS NON-NEGATIVE! ★★★")
                print(f"    → H is a POSITIVE linear combination of Schur-concave functions!")
                print(f"    → H is SCHUR-CONCAVE → Paley maximizes H!")
            else:
                neg_terms = [f"e_{j+1}" for j in range(1, n_vars) if coeffs[j] < 0]
                print(f"    Negative terms: {neg_terms}")

    # === Try: does e_m (= product) alone explain H? ===
    print(f"\n--- Product alone: H vs e_{m}(y²) = Π y_k² ---")
    for d in reps:
        print(f"  H={d['H']}: e_{m} = {d['e'][m]:.6f} {'(Paley)' if d['paley'] else ''}")

    # Fit H = a + b·e_m
    if n_reps >= 2:
        em_vals = np.array([d['e'][m] for d in reps])
        A_em = np.column_stack([np.ones(n_reps), em_vals])
        coeffs_em = np.linalg.lstsq(A_em, b_vec, rcond=None)[0]
        fitted = A_em @ coeffs_em
        max_err = np.max(np.abs(fitted - b_vec))
        print(f"\n  H ≈ {coeffs_em[0]:.2f} + {coeffs_em[1]:.6f}·Πy²  (max error: {max_err:.2f})")

    # === Log-product (sum of logs) ===
    print(f"\n--- Log-product: H vs Σ log(y_k²) ---")
    for d in reps:
        logsum = sum(np.log(xi) for xi in d['x'] if xi > 1e-15)
        d['logprod'] = logsum
        print(f"  H={d['H']}: Σlog(y²) = {logsum:.6f} {'(Paley)' if d['paley'] else ''}")

    # Fit H = a + b·Σlog
    logprod_vals = np.array([d['logprod'] for d in reps])
    A_log = np.column_stack([np.ones(n_reps), logprod_vals])
    coeffs_log = np.linalg.lstsq(A_log, b_vec, rcond=None)[0]
    fitted_log = A_log @ coeffs_log
    max_err_log = np.max(np.abs(fitted_log - b_vec))
    print(f"\n  H ≈ {coeffs_log[0]:.2f} + {coeffs_log[1]:.4f}·Σlog(y²)  (max error: {max_err_log:.2f})")

    # Fit H = a + b·Σlog + c·(Σlog)²
    logprod2 = logprod_vals**2
    A_log2 = np.column_stack([np.ones(n_reps), logprod_vals, logprod2])
    coeffs_log2 = np.linalg.lstsq(A_log2, b_vec, rcond=None)[0]
    fitted_log2 = A_log2 @ coeffs_log2
    max_err_log2 = np.max(np.abs(fitted_log2 - b_vec))
    print(f"  H ≈ {coeffs_log2[0]:.2f} + {coeffs_log2[1]:.4f}·Σlog + {coeffs_log2[2]:.4f}·(Σlog)²  (err: {max_err_log2:.2f})")

    # === RATIO: H / e_m ===
    print(f"\n--- H / Π(y²) ratio ---")
    for d in reps:
        ratio = d['H'] / d['e'][m] if d['e'][m] > 1e-15 else float('inf')
        print(f"  H={d['H']}: H/Πy² = {ratio:.4f} {'(Paley)' if d['paley'] else ''}")

    # === CHECK: Is Σe_k * 2^k related to OCF? ===
    print(f"\n--- OCF-like: Sum 2^k * e_k(y^2) ---")
    for d in reps:
        ocf_like = sum(2**k * d['e'][k] for k in range(m+1))
        print(f"  H={d['H']}: Σ2^k·e_k = {ocf_like:.4f} {'(Paley)' if d['paley'] else ''}")

    # === The key: e_k of y² vs α_k of OCF ===
    print(f"\n--- Connection: e_k(y²) vs OCF α_k ---")
    print(f"  OCF: H = 1 + 2α₁ + 4α₂ + 8α₃ + ... where α_k counts k-tuples of disjoint odd cycles")
    print("  If e_k(y^2) ~ alpha_k, then H = Sum 2^k*alpha_k <-> Sum 2^k*e_k(y^2)")
    print(f"  This would be the BRIDGE between spectral and combinatorial!")

    # === FULL analysis with ALL symmetric function bases ===
    print(f"\n{'='*70}")
    print(f"COMPLETE SYMMETRIC FUNCTION ANALYSIS")
    print(f"{'='*70}")

    # The Schur functions s_λ form another basis.
    # For m variables and degree d, the Schur functions indexed by
    # partitions λ of d with at most m parts form a basis.
    # s_λ(x) is Schur-convex iff λ has a single part (= power sum e_1^d).

    # Actually: ALL Schur functions s_λ(x_1,...,x_m) with x_i >= 0 are Schur-convex!
    # No, that's wrong. Schur functions s_λ for different λ can be either Schur-convex or -concave.
    # Actually: Schur-convexity of s_λ(x) requires checking.
    # The monomial symmetric functions m_λ = Σ_{σ∈S_m} x^{σ(λ)} ARE Schur-convex for x >= 0.

    # KEY FACT: e_k(x) are Schur-CONCAVE for k >= 1.
    # Proof: e_k(x) = s_{(1^k)}(x), the Schur function for the partition (1,1,...,1).
    # And s_{(1^k)} is an ALTERNATING sum of monomials.
    # Actually, the standard result is:
    #   p_k (power sum) is Schur-convex for x >= 0
    #   e_k (elementary) is Schur-concave for x >= 0
    # These are conjugate notions under the Newton identities.

    # So if H = c_0 + Σ_{k>=2} c_k · e_k with c_k >= 0, H is Schur-concave.
    # And if H = c_0 - Σ_{k>=2} d_k · p_k with d_k >= 0, H is also Schur-concave.

    # We already know the power sum version has mixed signs at p=11.
    # But the e_k version might have all positive coefficients!

    # Already done above, but let me also try: H as POLYNOMIAL in e_2,...,e_m
    # (allowing products like e_2² or e_2·e_3)
    if n_reps >= 5 and m >= 3:
        print(f"\n  Fitting H = c₀ + c₂e₂ + c₃e₃ + c₄e₄ + c₅e₅:")
        A5 = np.zeros((n_reps, 5))
        for i, d in enumerate(reps):
            A5[i] = [1, d['e'][2], d['e'][3], d['e'].get(4,0), d['e'].get(5,0)]
        rank5 = np.linalg.matrix_rank(A5)
        print(f"    Rank: {rank5}")
        if rank5 >= 5 and n_reps >= 5:
            coeffs5 = np.linalg.lstsq(A5, b_vec, rcond=None)[0]
            residual5 = np.max(np.abs(A5 @ coeffs5 - b_vec))
            print(f"    H = {coeffs5[0]:.4f} + {coeffs5[1]:.6f}e₂ + {coeffs5[2]:.6f}e₃ + {coeffs5[3]:.6f}e₄ + {coeffs5[4]:.6f}e₅")
            print(f"    Residual: {residual5:.2e}")

    # === THE ULTIMATE TEST: Verify e_k Schur-concavity claim ===
    print(f"\n--- Verifying e_k is Schur-concave on actual data ---")

    paley_x = [d['x'] for d in reps if d['paley']][0] if any(d['paley'] for d in reps) else None

    if paley_x:
        # For each other tournament, check:
        # Paley << other (Paley majorized by other)
        # e_k(Paley) >= e_k(other) for all k (Schur-concave means smaller on more spread)
        paley_sorted = sorted(paley_x, reverse=True)
        paley_e = {k: sum(np.prod(list(c)) for c in combinations(paley_x, k)) for k in range(1, m+1)}

        for d in reps:
            if d['paley']:
                continue
            other_sorted = sorted(d['x'], reverse=True)
            partial_p = np.cumsum(paley_sorted)
            partial_o = np.cumsum(other_sorted)
            majorized = all(partial_p[k] <= partial_o[k] + 1e-10 for k in range(m))

            print(f"\n  H={d['H']}: Paley << this? {majorized}")
            if majorized:
                for k in range(2, m+1):
                    ep = paley_e[k]
                    eo = d['e'][k]
                    print(f"    e_{k}: Paley={ep:.4f}, other={eo:.4f}, diff={ep-eo:+.4f} {'✓' if ep >= eo - 1e-6 else '✗ VIOLATION!'}")

    return reps

for p in [7, 11]:
    reps = analyze(p)

print(f"\n{'='*70}")
print(f"SYNTHESIS: ELEMENTARY SYMMETRIC POLYNOMIAL ROUTE TO PROOF")
print(f"{'='*70}")
print(f"""
At p=7:
  H = 198.19 - 1.0·e₂(y²)  [e₂ = Σ_{i<j} y_i²y_j²]
  (Wait, this needs checking — e₂ vs s₂ = Σy⁴)

At p=11:
  Need to express H in e_k basis and check signs.

KEY THEORETICAL QUESTION:
  Does the OCF formula H = 1 + 2α₁ + 4α₂ + ...
  translate into H = c₀ + Σ_{k>=2} c_k · e_k(y²)?

  If so, the proof is: each α_k (disjoint k-tuples of odd cycles)
  corresponds to e_k (product of k eigenvalue magnitudes),
  and 2^k · α_k = c_k · e_k with c_k > 0.

  Since e_k is Schur-concave and c_k > 0,
  H is Schur-concave → maximized at uniform → Paley wins.
""")
