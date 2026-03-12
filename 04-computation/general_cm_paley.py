#!/usr/bin/env python3
"""
general_cm_paley.py — Closed forms for c_m(Paley_p) for ALL odd m.

Extends THM-130 (c_5 closed form) to general odd m.

Method:
  λ_k = (χ(k)·g - 1)/2  for k ≠ 0,  λ_0 = (p-1)/2

  c_m = (1/m) · Σ_{k=0}^{p-1} λ_k^m = (1/m) · [λ_0^m + (1/2^m) · Σ_{k≠0} (χ(k)g - 1)^m]

  Binomial expansion: (χg - 1)^m = Σ_{j=0}^m C(m,j) (χg)^j (-1)^{m-j}
                                  = Σ_j C(m,j) χ(k)^j g^j (-1)^{m-j}

  Key reductions:
    χ(k)^j = χ(k) for j odd, 1 for j even  (since χ(k)=±1 for k≠0)
    g^{2s} = (χ(-1)·p)^s
    g^{2s+1} = (χ(-1)·p)^s · g

  Summing over k≠0:
    Σ_{k≠0} χ(k) = 0  (equal QR and NQR)
    Σ_{k≠0} 1 = p-1

  So terms with χ(k) (odd j) vanish, terms with 1 (even j) give factor (p-1).

  Result: Σ_{k≠0} (χg-1)^m = (p-1) · Σ_{j even} C(m,j) g^j (-1)^{m-j}

  For m odd: (-1)^{m-j} = (-1)^{m} · (-1)^{-j} = -(-1)^{-j}
    j even: (-1)^{m-j} = -1 · (-1)^{-j} = -1  (since j even => (-1)^j=1)
    Wait: (-1)^{m-j} where m is odd, j is even => m-j is odd => (-1)^{m-j} = -1

  So: Σ_{k≠0} (χg-1)^m = -(p-1) · Σ_{j=0,2,4,...} C(m,j) g^j

  Now g^j for j even: g^j = g^{2s} = (χ(-1)p)^s  where s = j/2
  Let ε = χ(-1) = (-1)^{(p-1)/2}.  Then g^{2s} = (εp)^s.

  S_even = Σ_{s=0}^{(m-1)/2} C(m, 2s) (εp)^s

  This is a polynomial in εp! And:
    c_m = (1/m) · [((p-1)/2)^m - (p-1)/(2^m) · S_even]

  We compute S_even symbolically and verify numerically.
"""

import numpy as np
from math import comb, gcd
from fractions import Fraction
from itertools import combinations
import sys

def legendre(a, p):
    """Legendre symbol (a/p)."""
    a = a % p
    if a == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

def paley_adjacency(p):
    """Adjacency matrix of Paley tournament T_p."""
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and legendre(j - i, p) == 1:
                A[i][j] = 1
    return A

def eigenvalues_direct(p):
    """Compute eigenvalues of Paley tournament adjacency matrix directly."""
    A = paley_adjacency(p)
    eigs = np.linalg.eigvals(A)
    # Sort by real part descending
    eigs = sorted(eigs, key=lambda x: -x.real)
    return eigs, A

def cm_direct(A, m):
    """Compute c_m = (number of directed m-cycles) = tr(A^m)/m."""
    Am = np.linalg.matrix_power(A, m)
    return int(round(np.trace(Am).real)) // m

def gauss_sum(p):
    """Compute quadratic Gauss sum g = Σ_{k=0}^{p-1} χ(k) ζ^k numerically."""
    chi_minus1 = legendre(-1, p)
    # g² = χ(-1)·p
    # g = sqrt(χ(-1)·p) with appropriate sign
    # For numerical verification, compute directly
    zeta = np.exp(2j * np.pi / p)
    g = sum(legendre(k, p) * zeta**k for k in range(1, p))
    return g

def eigenvalues_formula(p):
    """Compute eigenvalues via formula λ_k = (χ(k)·g - 1)/2."""
    g = gauss_sum(p)
    eigs = []
    eigs.append((p-1)/2)  # λ_0
    for k in range(1, p):
        chi_k = legendre(k, p)
        eigs.append((chi_k * g - 1) / 2)
    return eigs

def S_even_exact(m, p):
    """
    Compute S_even = Σ_{s=0}^{(m-1)/2} C(m, 2s) · (ε·p)^s
    where ε = χ(-1) = (-1)^{(p-1)/2}.

    Returns exact Fraction.
    """
    eps = (-1)**((p-1)//2)  # χ(-1)
    ep = eps * p  # ε·p as integer (or -p)

    S = Fraction(0)
    for s in range((m-1)//2 + 1):
        j = 2 * s
        S += Fraction(comb(m, j)) * Fraction(ep)**s
    return S

def cm_formula_exact(m, p):
    """
    Exact c_m(Paley_p) for odd m using the closed form:
    c_m = (1/m) · [((p-1)/2)^m - (p-1)/2^m · S_even]
    """
    assert m % 2 == 1, "m must be odd"

    S = S_even_exact(m, p)

    lambda0_m = Fraction(p-1, 2)**m
    correction = Fraction(p-1) * S / Fraction(2)**m

    total = lambda0_m - correction
    cm = total / m

    return cm

def cm_formula_symbolic(m):
    """
    Return a symbolic description of c_m in terms of p and ε=χ(-1).

    c_m = (1/m) · [((p-1)/2)^m - (p-1)/2^m · Σ_{s=0}^{(m-1)/2} C(m,2s)(εp)^s]

    We expand ((p-1)/2)^m by binomial theorem too, and collect powers of p.
    """
    # The key formula terms
    terms = []
    for s in range((m-1)//2 + 1):
        j = 2*s
        c = comb(m, j)
        terms.append((s, c))

    return terms

def print_separator():
    print("=" * 80)

def verify_cm(primes, ms):
    """Verify formula c_m against direct computation."""
    print_separator()
    print("VERIFICATION: c_m(Paley_p) formula vs direct computation")
    print_separator()

    all_ok = True
    for p in primes:
        _, A = eigenvalues_direct(p)
        for m in ms:
            if m > p:
                continue
            cm_dir = cm_direct(A, m)
            cm_form = cm_formula_exact(m, p)

            match = (cm_form == Fraction(cm_dir))
            status = "OK" if match else "FAIL"
            if not match:
                all_ok = False

            print(f"  p={p:2d}, m={m:2d}: direct={cm_dir:>12d}, formula={cm_form!s:>20s} [{status}]")

    print(f"\nOverall verification: {'ALL PASSED' if all_ok else 'SOME FAILED'}")
    return all_ok

def compute_cm_table(primes, ms):
    """Compute c_m(Paley_p) for all m, p combinations."""
    print_separator()
    print("TABLE: c_m(Paley_p) exact values")
    print_separator()

    # Header
    col_label = 'm\\p'
    header = f"{col_label:>6s}"
    for p in primes:
        header += f"  {p:>14d}"
    print(header)
    print("-" * len(header))

    for m in ms:
        row = f"{m:>6d}"
        for p in primes:
            if m > p:
                row += f"  {'—':>14s}"
            else:
                cm = cm_formula_exact(m, p)
                row += f"  {int(cm):>14d}"
        print(row)

def analyze_general_formula(ms):
    """Analyze the structure of the formula for general odd m."""
    print_separator()
    print("GENERAL FORMULA ANALYSIS")
    print_separator()

    print("""
For Paley tournament T_p with p prime, p ≡ 1 or 3 (mod 4):

  c_m(T_p) = (1/m) · [((p-1)/2)^m - (p-1)/2^m · S_m(p)]

where S_m(p) = Σ_{s=0}^{(m-1)/2} C(m, 2s) · (ε·p)^s

and ε = χ(-1) = (-1)^{(p-1)/2} = { +1 if p ≡ 1 (mod 4)
                                    { -1 if p ≡ 3 (mod 4)
""")

    for m in ms:
        print(f"\n--- m = {m} ---")
        terms = cm_formula_symbolic(m)

        # S_m as polynomial in (εp)
        print(f"  S_{m}(p) = ", end="")
        parts = []
        for s, c in terms:
            if s == 0:
                parts.append(f"{c}")
            else:
                parts.append(f"{c}·(εp)^{s}")
        print(" + ".join(parts))

        # Expand for ε=+1 (p ≡ 1 mod 4) and ε=-1 (p ≡ 3 mod 4)
        for eps_val, label in [(1, "p ≡ 1 (mod 4)"), (-1, "p ≡ 3 (mod 4)")]:
            print(f"  When {label} (ε={eps_val:+d}):")
            # S as polynomial in p
            poly = {}  # degree -> coefficient
            for s, c in terms:
                coeff = c * (eps_val ** s)
                poly[s] = poly.get(s, 0) + coeff

            poly_str = []
            for deg in sorted(poly.keys(), reverse=True):
                coeff = poly[deg]
                if coeff == 0:
                    continue
                if deg == 0:
                    poly_str.append(f"{coeff}")
                elif deg == 1:
                    poly_str.append(f"{coeff}p")
                else:
                    poly_str.append(f"{coeff}p^{deg}")
            print(f"    S_{m} = {' + '.join(poly_str)}")

            # Now the full formula: c_m = [((p-1)/2)^m - (p-1)·S/2^m] / m
            # Let's compute symbolically for small m
            # Actually let's just show the polynomial form

            # Numerically verify the polynomial
            test_p = 7 if eps_val == -1 else 13
            S_num = sum(c * (eps_val * test_p)**s for s, c in terms)
            cm_num = cm_formula_exact(m, test_p)
            if m <= test_p:
                print(f"    Verify p={test_p}: S={S_num}, c_{m}={cm_num} = {int(cm_num)}")

def check_maximization(primes, ms):
    """Check whether Paley maximizes c_m among ALL tournaments (not just circulant)."""
    print_separator()
    print("MAXIMIZATION CHECK: Does Paley maximize c_m among all tournaments?")
    print_separator()

    from itertools import combinations as combs

    for p in primes:
        if p > 11:  # Only feasible for small p
            continue

        # Generate all tournaments on p vertices
        edges = [(i, j) for i in range(p) for j in range(i+1, p)]
        n_edges = len(edges)

        paley_A = paley_adjacency(p)

        max_cm = {}
        paley_cm = {}
        max_achieved_by = {}

        for m in ms:
            if m > p:
                continue
            paley_cm[m] = cm_direct(paley_A, m)
            max_cm[m] = paley_cm[m]
            max_achieved_by[m] = 1  # count how many achieve max

        # Check all tournaments
        n_tournaments = 2**n_edges

        if n_tournaments > 200000:
            print(f"\n  p={p}: {n_tournaments} tournaments — sampling 50000 random")
            import random
            sample_size = 50000
            for _ in range(sample_size):
                A = np.zeros((p, p), dtype=int)
                for i, j in edges:
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                for m in ms:
                    if m > p:
                        continue
                    cm = cm_direct(A, m)
                    if cm > max_cm[m]:
                        max_cm[m] = cm
                        max_achieved_by[m] = 1
                    elif cm == max_cm[m]:
                        max_achieved_by[m] += 1
        else:
            print(f"\n  p={p}: {n_tournaments} tournaments — exhaustive check")
            for bits in range(n_tournaments):
                A = np.zeros((p, p), dtype=int)
                for idx, (i, j) in enumerate(edges):
                    if (bits >> idx) & 1:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                for m in ms:
                    if m > p:
                        continue
                    cm = cm_direct(A, m)
                    if cm > max_cm[m]:
                        max_cm[m] = cm
                        max_achieved_by[m] = 1
                    elif cm == max_cm[m]:
                        max_achieved_by[m] += 1

        for m in ms:
            if m > p:
                continue
            is_max = "YES" if max_cm[m] == paley_cm[m] else "NO"
            print(f"    m={m}: Paley c_{m}={paley_cm[m]}, max c_{m}={max_cm[m]}, "
                  f"Paley maximizes: {is_max}")

def check_hamiltonian(primes):
    """For m = p, c_p counts Hamiltonian directed cycles. Check connection to H."""
    print_separator()
    print("HAMILTONIAN CYCLES: c_p(T_p) = number of directed p-cycles")
    print_separator()

    for p in primes:
        if p > 13:
            print(f"  p={p}: skipped (too large for direct computation)")
            continue

        _, A = eigenvalues_direct(p)

        # c_p from direct computation
        cp_dir = cm_direct(A, p)

        # c_p from formula
        cp_form = cm_formula_exact(p, p)

        # H = number of Hamiltonian paths
        # For tournament, H = Σ λ_k^{n-1} ... no, H = permanent-like
        # Actually tr(A^p) counts closed walks of length p = p * (directed p-cycles)
        # But a Hamiltonian cycle visits all p vertices, so it IS a p-cycle (the only p-cycle possible)
        # Wait: a p-cycle in a tournament on p vertices IS a Hamiltonian cycle.

        # H (Hamiltonian paths) is different from Hamiltonian cycles.
        # c_p = number of directed Hamiltonian cycles (each counted once).
        # Number of Hamiltonian cycles = tr(A^p) / p.
        # But also = (p-1)!/2 * (fraction that are cycles) -- not directly.

        # Compute H (Hamiltonian paths) for comparison
        # H = number of Hamiltonian paths = permanent of A in a sense...
        # Actually use the I(Ω,2) formula or direct count

        from itertools import permutations
        H_count = 0
        for perm in permutations(range(p)):
            is_path = all(A[perm[i]][perm[i+1]] == 1 for i in range(p-1))
            if is_path:
                H_count += 1

        print(f"  p={p}:")
        print(f"    c_p = {int(cp_form)} (Hamiltonian directed cycles)")
        print(f"    H   = {H_count} (Hamiltonian paths)")
        print(f"    Ratio H/c_p = {H_count/int(cp_form):.6f}" if int(cp_form) > 0 else "    c_p = 0")
        print(f"    Formula matches direct: {cp_form == Fraction(cp_dir)}")

def generating_function_analysis(primes):
    """
    Analyze the generating function G(x) = Σ_{m≥1} c_m · x^m.

    Since c_m = tr(A^m)/m, we have:
      Σ_{m≥1} c_m x^m/m ... no, c_m = tr(A^m)/m, so
      Σ_{m≥1} tr(A^m)/m · x^m = tr(Σ_{m≥1} A^m x^m / m) = tr(-log(I - xA))
      = -log det(I - xA) = Σ_k -log(1 - x·λ_k)

    Actually: Σ_{m≥1} c_m x^m = Σ_{m≥1} (tr(A^m)/m) x^m = -tr(log(I-xA))

    For Paley: det(I - xA) = (1 - x·λ_0) · Π_{k≠0} (1 - x·λ_k)

    Since eigenvalues come in conjugate pairs (or ±1 shifted):
    λ_k = (χ(k)g - 1)/2, and there are (p-1)/2 QR and (p-1)/2 NQR,

    For QR:  λ = (g-1)/2
    For NQR: λ = (-g-1)/2

    So Π_{k: QR} (1 - xλ_k) = (1 - x(g-1)/2)^{(p-1)/2}
       Π_{k: NQR} (1 - xλ_k) = (1 - x(-g-1)/2)^{(p-1)/2}

    det(I - xA) = (1 - x(p-1)/2) · (1 - x(g-1)/2)^{(p-1)/2} · (1 + x(g+1)/2)^{(p-1)/2}

    Let α = (g-1)/2, β = -(g+1)/2 = (-g-1)/2. Note α - β = g, α + β = -1.
    αβ = (g-1)(-g-1)/4 = -(g²-1)/4 = -(χ(-1)p - 1)/4

    So: det(I - xA) = (1 - x(p-1)/2) · [(1-xα)(1-xβ)]^{(p-1)/2}
                     = (1 - x(p-1)/2) · [1 - x(α+β) + x²αβ]^{(p-1)/2}
                     = (1 - x(p-1)/2) · [1 + x + x²(1-χ(-1)p)/4]^{(p-1)/2}
    """
    print_separator()
    print("GENERATING FUNCTION: det(I - xA) for Paley tournaments")
    print_separator()

    print("""
THEOREM: For Paley tournament T_p:

  det(I - xA) = (1 - x·(p-1)/2) · [1 + x + x²·(1 - εp)/4]^{(p-1)/2}

where ε = χ(-1) = (-1)^{(p-1)/2}.

  - If p ≡ 1 (mod 4): ε=1,  quadratic = 1 + x + x²(1-p)/4
  - If p ≡ 3 (mod 4): ε=-1, quadratic = 1 + x + x²(1+p)/4

This gives: -log det(I - xA) = Σ_{m≥1} c_m x^m

Note: (1-xα)(1-xβ) = 1 + x + x²αβ where αβ = (1 - εp)/4.
""")

    for p in primes:
        eps = (-1)**((p-1)//2)
        alpha_beta = Fraction(1 - eps*p, 4)

        print(f"  p={p} (ε={eps:+d}):")
        print(f"    αβ = (1 - εp)/4 = {alpha_beta}")
        print(f"    det(I-xA) = (1 - {Fraction(p-1,2)}x) · (1 + x + {alpha_beta}x²)^{(p-1)//2}")

        # Verify numerically
        _, A = eigenvalues_direct(p)
        eigs = np.linalg.eigvals(A)

        # det(I - xA) at x = 0.1
        x_test = 0.1
        det_direct = np.prod(1 - x_test * eigs)

        # Formula
        q = float(alpha_beta)
        det_formula = (1 - x_test*(p-1)/2) * (1 + x_test + q*x_test**2)**((p-1)//2)

        print(f"    Verify at x=0.1: direct={det_direct.real:.10f}, formula={det_formula:.10f}, "
              f"match={abs(det_direct.real - det_formula) < 1e-6}")

def formula_simplification(ms, primes):
    """
    Try to find simplified closed form for c_m in terms of p.

    c_m = (1/m) · [((p-1)/2)^m - (p-1)/2^m · S_m]

    Let's expand and collect terms as polynomials in p.
    """
    print_separator()
    print("SIMPLIFIED CLOSED FORMS for c_m(Paley_p)")
    print_separator()

    for m in ms:
        print(f"\n=== m = {m} ===")

        # For each residue class of p mod 4
        for eps_val, label in [(1, "p ≡ 1 (mod 4)"), (-1, "p ≡ 3 (mod 4)")]:
            print(f"\n  {label} (ε = {eps_val:+d}):")

            # c_m = [((p-1)/2)^m - (p-1)/2^m · S_m] / m
            # S_m = Σ_{s=0}^{(m-1)/2} C(m,2s) · (εp)^s

            # Let's compute c_m as rational function of p
            # by evaluating at enough primes and interpolating

            # Get primes in this residue class
            test_primes = [p for p in primes if (-1)**((p-1)//2) == eps_val and p >= m]

            # Also add more primes for interpolation
            extra_primes = []
            for pp in range(3, 200):
                if all(pp % d != 0 for d in range(2, int(pp**0.5)+1)):
                    if (-1)**((pp-1)//2) == eps_val and pp >= m:
                        extra_primes.append(pp)

            if len(extra_primes) < m + 2:
                print("    Not enough primes for interpolation")
                continue

            # Compute c_m for these primes
            values = [(pp, cm_formula_exact(m, pp)) for pp in extra_primes[:m+2]]

            # c_m should be a polynomial in p of degree m.
            # Let's verify by checking degree: ((p-1)/2)^m is degree m in p
            # S_m has max degree (m-1)/2, times (p-1) gives degree (m+1)/2
            # So c_m is degree m in p.

            # Factor out p(p-1)/[m · 2^m] and find remaining polynomial
            print(f"    c_{m} = p(p-1) · f(p) / ({m} · 2^{m})")

            # Compute f(p) = c_m · m · 2^m / (p(p-1))
            f_values = []
            for pp, cm_val in values:
                f = cm_val * m * 2**m / (pp * (pp - 1))
                f_values.append((pp, f))

            # Check if f(p) is a polynomial with integer coefficients
            print(f"    f(p) values: ", end="")
            for pp, f in f_values[:6]:
                print(f"f({pp})={f}, ", end="")
            print()

            # Lagrange interpolation to find f(p)
            # f should be degree m-2 polynomial
            n_pts = min(len(f_values), m+1)
            pts = f_values[:n_pts]

            # Use Fraction arithmetic for exact interpolation
            def lagrange_eval(pts, x):
                result = Fraction(0)
                for i, (xi, yi) in enumerate(pts):
                    term = Fraction(yi)
                    for j, (xj, _) in enumerate(pts):
                        if i != j:
                            term *= Fraction(x - xj, xi - xj)
                    result += term
                return result

            # Find polynomial coefficients by evaluating at 0, 1, 2, ...
            # Actually, let's just print the polynomial
            # Check that f is indeed polynomial of degree m-2
            deg = m - 2
            if n_pts >= deg + 2:
                # Verify: evaluate at an extra point
                test_val = lagrange_eval(pts[:deg+1], pts[deg+1][0])
                is_poly = (test_val == pts[deg+1][1])
                if not is_poly:
                    # Try degree m-1
                    deg = m - 1
                    if n_pts >= deg + 2:
                        test_val = lagrange_eval(pts[:deg+1], pts[deg+1][0])
                        is_poly = (test_val == pts[deg+1][1])

            # Get coefficients via Newton's forward differences or direct
            # Build the polynomial from the interpolation points
            coeffs = [Fraction(0)] * (deg + 1)

            # Use the Vandermonde approach
            from functools import reduce

            # Evaluate at p = 0, 1, ..., deg to get coefficients
            eval_pts = list(range(deg + 1))
            eval_vals = [lagrange_eval(pts[:deg+1], x) for x in eval_pts]

            # Newton forward differences
            diffs = [list(eval_vals)]
            for level in range(deg):
                new_diffs = []
                for i in range(len(diffs[-1]) - 1):
                    new_diffs.append(diffs[-1][i+1] - diffs[-1][i])
                diffs.append(new_diffs)

            # Coefficients in falling factorial basis
            newton_coeffs = [d[0] for d in diffs]

            # Convert to standard polynomial basis
            # f(p) = Σ newton_coeffs[k] · C(p, k) · k!  ... actually
            # f(p) = Σ newton_coeffs[k] · (p choose k) · k! = Σ newton_coeffs[k] · p(p-1)...(p-k+1)/1
            # Better: just display as evaluated formula

            # Direct display of polynomial
            poly_str = f"    f(p) = "
            # Get coefficients by solving Vandermonde system
            # Using exact arithmetic
            n_coeffs = deg + 1
            # V[i][j] = pts[i][0]^j
            xs = [pts[i][0] for i in range(n_coeffs)]
            ys = [pts[i][1] for i in range(n_coeffs)]

            # Solve V @ c = y using Gaussian elimination
            # V is Vandermonde: V[i][j] = x_i^j
            V = [[Fraction(x**j) for j in range(n_coeffs)] for i, x in enumerate(xs)]
            y_vec = [Fraction(y) for y in ys]

            # Augmented matrix
            aug = [row + [y_vec[i]] for i, row in enumerate(V)]

            for col in range(n_coeffs):
                # Find pivot
                pivot = None
                for row in range(col, n_coeffs):
                    if aug[row][col] != 0:
                        pivot = row
                        break
                if pivot is None:
                    break
                aug[col], aug[pivot] = aug[pivot], aug[col]

                for row in range(n_coeffs):
                    if row != col and aug[row][col] != 0:
                        factor = aug[row][col] / aug[col][col]
                        for j in range(n_coeffs + 1):
                            aug[row][j] -= factor * aug[col][j]

            poly_coeffs = [aug[i][n_coeffs] / aug[i][i] for i in range(n_coeffs)]

            # Display
            terms = []
            for j in range(n_coeffs - 1, -1, -1):
                c = poly_coeffs[j]
                if c == 0:
                    continue
                if j == 0:
                    terms.append(f"{c}")
                elif j == 1:
                    terms.append(f"{c}·p")
                else:
                    terms.append(f"{c}·p^{j}")

            print(poly_str + " + ".join(terms).replace("+ -", "- "))

def paley_maximization_all_m(primes_small, ms):
    """
    Detailed check: does Paley maximize c_m for ALL odd m, or just some?
    Compare against all tournaments for small p.
    """
    print_separator()
    print("PALEY MAXIMIZATION: c_m for ALL odd m?")
    print_separator()

    for p in [5, 7]:
        edges = [(i, j) for i in range(p) for j in range(i+1, p)]
        n_edges = len(edges)
        n_tournaments = 2**n_edges

        paley_A = paley_adjacency(p)

        print(f"\n  p={p} ({n_tournaments} tournaments, exhaustive):")

        # Compute c_m for all tournaments
        results = {}  # m -> list of c_m values
        paley_vals = {}

        for m in ms:
            if m > p:
                continue
            results[m] = []
            paley_vals[m] = cm_direct(paley_A, m)

        for bits in range(n_tournaments):
            A = np.zeros((p, p), dtype=int)
            for idx, (i, j) in enumerate(edges):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

            for m in results:
                results[m].append(cm_direct(A, m))

        for m in sorted(results.keys()):
            max_cm = max(results[m])
            is_max = paley_vals[m] == max_cm
            n_at_max = results[m].count(max_cm)
            rank = sorted(set(results[m]), reverse=True).index(paley_vals[m]) + 1
            print(f"    m={m}: Paley c_{m}={paley_vals[m]:>6d}, max={max_cm:>6d}, "
                  f"Paley rank={rank}/{len(set(results[m]))}, "
                  f"max achieved by {n_at_max} tournaments, "
                  f"{'MAXIMIZER' if is_max else 'NOT MAX'}")

def eigenvalue_magnitude_analysis(primes):
    """
    All non-trivial eigenvalues have |λ_k| = sqrt((p+1)/4) for Paley.
    Analyze what this means for the generating function.
    """
    print_separator()
    print("EIGENVALUE MAGNITUDES AND CONVERGENCE RADIUS")
    print_separator()

    for p in primes:
        R = np.sqrt((p+1)/4)  # common magnitude of non-trivial eigenvalues
        lambda0 = (p-1)/2

        # Convergence radius of GF is 1/max|λ_k| = 1/λ_0 = 2/(p-1)
        # since λ_0 > R for p ≥ 5

        print(f"  p={p}: λ_0 = {lambda0:.1f}, |λ_k| = √{(p+1)/4:.2f} = {R:.4f}")
        print(f"    λ_0 > |λ_k|: {lambda0 > R}")
        print(f"    Convergence radius of GF: 1/λ_0 = {1/lambda0:.6f}")

        # For large m, c_m ~ λ_0^m / m (dominant eigenvalue)
        # The correction from non-trivial eigenvalues is O(R^m / m)
        # Ratio: R/λ_0 = √((p+1)/4) / ((p-1)/2) = √(p+1) / (p-1)

        ratio = R / lambda0
        print(f"    Asymptotic ratio |λ_k|/λ_0 = {ratio:.6f}")
        print(f"    For large m: c_m ≈ λ_0^m/m · [1 + O({ratio:.4f}^m)]")

def main():
    primes = [7, 11, 13, 19, 23]
    ms = [3, 5, 7, 9, 11]

    print("=" * 80)
    print("GENERAL c_m(Paley_p) FOR ALL ODD m")
    print("Extension of THM-130 to arbitrary odd cycle lengths")
    print("=" * 80)

    # 1. Verify formula
    verify_cm(primes, ms)

    # 2. Full table
    compute_cm_table(primes, ms)

    # 3. General formula analysis
    analyze_general_formula(ms)

    # 4. Simplified closed forms
    formula_simplification(ms, primes)

    # 5. Generating function
    generating_function_analysis(primes)

    # 6. Eigenvalue magnitude analysis
    eigenvalue_magnitude_analysis(primes)

    # 7. Hamiltonian cycles (m = p)
    check_hamiltonian([5, 7, 11, 13])

    # 8. Maximization check
    paley_maximization_all_m([5, 7], ms)

    # 9. Summary
    print_separator()
    print("SUMMARY: MASTER FORMULA")
    print_separator()
    print("""
THEOREM (General c_m for Paley tournaments):

For p prime, m odd, m ≤ p:

  c_m(T_p) = (1/m) · [((p-1)/2)^m - (p-1)/2^m · S_m(εp)]

where:
  ε = χ(-1) = (-1)^{(p-1)/2}
  S_m(t) = Σ_{s=0}^{(m-1)/2} C(m, 2s) · t^s

Equivalently, S_m(t) = Re[(1 + √t)^m] when t > 0 (i.e., cosh-type),
             S_m(t) = Re[(1 + i√|t|)^m] when t < 0.

Since (1+z)^m = Σ C(m,j) z^j and z = √(εp):
  (1+z)^m + (1-z)^m = 2·Σ_{s≥0} C(m,2s) z^{2s} = 2·S_m(εp)

So: S_m(εp) = [(1 + √(εp))^m + (1 - √(εp))^m] / 2

For p ≡ 1 (mod 4): S_m = [(1+√p)^m + (1-√p)^m] / 2  (real)
For p ≡ 3 (mod 4): S_m = [(1+i√p)^m + (1-i√p)^m] / 2 = Re[(1+i√p)^m]
                        = (1+p)^{m/2} · cos(m·arctan(√p))

GENERATING FUNCTION:
  det(I - xA)^{-1} = 1/[(1 - (p-1)x/2) · (1 + x + (1-εp)x²/4)^{(p-1)/2}]

  -log det(I - xA) = Σ_{m≥1} c_m · x^m

The quadratic factor 1 + x + (1-εp)x²/4 has roots at x = 2/α, 2/β where
α = (g-1)/2, β = (-g-1)/2 are the two distinct non-trivial eigenvalues.

|α| = |β| = √((p+1)/4), so the non-trivial eigenvalue contribution to the GF
has uniform convergence radius 2/√(p+1).
""")

if __name__ == "__main__":
    main()
