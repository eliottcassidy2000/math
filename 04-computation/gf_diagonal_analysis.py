#!/usr/bin/env python3
"""
Diagonal and cross-section analysis of G_T(t, x).

G_T(t, x) = A_n(t) + sum_I x^{parts(I)} * I(T) * A_{f_I+1}(t) * (t-1)^{n-1-f_I}

Special evaluations:
  G_T(t, 0) = A_n(t)                        (Eulerian polynomial, T-independent)
  G_T(0, x) = I(Omega(T), x)                (independence polynomial)
  G_T(1, x) = n! for all T                  (trivial)
  G_T(t, 2) = E_T(t) = sum_k a_k(T) t^k    (tournament Eulerian polynomial)

This script investigates:
  1. G_T(t, t) — the diagonal
  2. G_T(t, 2t) — scaled diagonal
  3. G_T(t, 2-t) — anti-diagonal (passes through (0,2)=H and (2,0)=A_n(2))
  4. G_T(t, x) = n! level set beyond t=1
  5. G_T(t, -1) — independence polynomial at x=-1

For n=7, the invariants are:
  t3 (f=4, parts=1), t5 (f=2, parts=1), t7 (f=0, parts=1), bc (f=2, parts=2)

opus-2026-03-07
"""

from itertools import combinations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction
import random

# ============================================================
# Helper functions
# ============================================================

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

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cl):
    if n < cl: return 0
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1
        for m in range(1, 1 << cl):
            for v in range(cl):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(cl):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cl) - 1
        total += sum(dp[full][v] for v in range(1, cl) if sub[v][0])
    return total

def count_bc(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

def forward_edge_dist_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    new_fwd = fwd + A[v][u]
                    key = (mask | (1 << u), u, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)

def eulerian_number(n, k):
    """Eulerian number A(n,k) with 0-indexed descent count."""
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

def eulerian_poly_eval(n, t):
    """Evaluate A_n(t) = sum_{k=0}^{n-1} A(n,k) t^k."""
    if n == 0:
        return 1
    return sum(eulerian_number(n, k) * t**k for k in range(n))

# ============================================================
# G_T(t, x) via the formula
# ============================================================

def compute_invariants_n7(A, n=7):
    """Compute (t3, t5, t7, bc) for a tournament on n=7 vertices."""
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t7 = count_directed_cycles(A, n, 7)
    bc = count_bc(A, n)
    return t3, t5, t7, bc

def G_T_formula(t, x, t3, t5, t7, bc, n=7):
    """
    G_T(t, x) = A_n(t) + sum_I x^{parts(I)} * I(T) * A_{f+1}(t) * (t-1)^{d-f}

    For n=7, d = n-1 = 6:
      t3: f=4, parts=1  -> x * t3 * A_5(t) * (t-1)^2
      t5: f=2, parts=1  -> x * t5 * A_3(t) * (t-1)^4
      t7: f=0, parts=1  -> x * t7 * A_1(t) * (t-1)^6
      bc: f=2, parts=2  -> x^2 * bc * A_3(t) * (t-1)^4
    """
    d = n - 1
    result = eulerian_poly_eval(n, t)
    # t3 term: f=4, parts=1
    result += x * t3 * eulerian_poly_eval(5, t) * (t - 1)**2
    # t5 term: f=2, parts=1
    result += x * t5 * eulerian_poly_eval(3, t) * (t - 1)**4
    # t7 term: f=0, parts=1
    result += x * t7 * eulerian_poly_eval(1, t) * (t - 1)**6
    # bc term: f=2, parts=2
    result += x**2 * bc * eulerian_poly_eval(3, t) * (t - 1)**4
    return result

def G_T_direct(t, dist, n=7):
    """
    Compute E_T(t) = sum_k a_k t^k directly from forward-edge distribution.
    This equals G_T(t, 2).
    """
    return sum(t**k * dist.get(k, 0) for k in range(n))

# ============================================================
# Verification: G_T(t, 2) should equal E_T(t)
# ============================================================

def verify_formula(n=7, num_tests=5):
    """Verify G_T(t, 2) = E_T(t) for random tournaments."""
    print("VERIFICATION: G_T(t, 2) = E_T(t)")
    print("=" * 70)

    for seed in range(num_tests):
        A = random_tournament(n, seed)
        t3, t5, t7, bc = compute_invariants_n7(A, n)
        dist = forward_edge_dist_dp(A, n)

        ok = True
        for t_val in [0, 1, 2, 3, -1, Fraction(1, 2)]:
            gf_val = G_T_formula(t_val, 2, t3, t5, t7, bc, n)
            et_val = G_T_direct(t_val, dist, n)
            if gf_val != et_val:
                print(f"  MISMATCH at seed={seed}, t={t_val}: G_T(t,2)={gf_val}, E_T(t)={et_val}")
                ok = False

        # Also verify G_T(0, x) = I(Omega, x)
        alpha1 = t3 + t5 + t7
        alpha2 = bc
        for x_val in [0, 1, 2, 3, -1]:
            gf_val = G_T_formula(0, x_val, t3, t5, t7, bc, n)
            iomega = 1 + alpha1 * x_val + alpha2 * x_val**2
            if gf_val != iomega:
                print(f"  MISMATCH at seed={seed}, x={x_val}: G_T(0,x)={gf_val}, I(Omega,x)={iomega}")
                ok = False

        if ok:
            print(f"  seed={seed}: ALL CHECKS PASS (t3={t3}, t5={t5}, t7={t7}, bc={bc})")

    print()

# ============================================================
# Analysis 1: Diagonal x = t
# ============================================================

def analysis_diagonal(n=7, num_tests=15):
    print("ANALYSIS 1: DIAGONAL x = t  (G_T(t, t))")
    print("=" * 70)

    nfact = factorial(n)

    for seed in range(num_tests):
        A = random_tournament(n, seed)
        t3, t5, t7, bc = compute_invariants_n7(A, n)
        dist = forward_edge_dist_dp(A, n)
        H = 1 + (t3 + t5 + t7) * 2 + bc * 4  # I(Omega, 2) = H(T)

        vals = {}
        for t_val in [0, Fraction(1,2), 1, Fraction(3,2), 2, -1, 3]:
            g = G_T_formula(t_val, t_val, t3, t5, t7, bc, n)
            vals[t_val] = g

        if seed < 5:
            print(f"\n  seed={seed}: t3={t3}, t5={t5}, t7={t7}, bc={bc}, H={H}")
            for t_val, g in vals.items():
                label = ""
                if t_val == 0:
                    label = f"  [= A_{n}(0) = 1, I(Omega,0) = 1]"
                elif t_val == 1:
                    label = f"  [= n! = {nfact}]"
                elif t_val == 2:
                    label = f"  [= E_T(2) = G_T(2,2)]"
                print(f"    G_T({t_val}, {t_val}) = {g}{label}")

    # Check: is G_T(t,t) T-independent at any t other than t=0 or t=1?
    print("\n  Checking T-independence of G_T(t,t) at various t...")
    test_vals = [Fraction(k, 4) for k in range(-4, 13)] + [Fraction(-1,1)]

    for t_val in test_vals:
        results = set()
        for seed in range(num_tests):
            A = random_tournament(n, seed)
            t3, t5, t7, bc = compute_invariants_n7(A, n)
            g = G_T_formula(t_val, t_val, t3, t5, t7, bc, n)
            results.add(g)
        if len(results) == 1:
            print(f"    t={t_val}: T-INDEPENDENT! Value = {results.pop()}")

    print()

# ============================================================
# Analysis 2: Scaled diagonal x = 2t
# ============================================================

def analysis_scaled_diagonal(n=7, num_tests=15):
    print("ANALYSIS 2: SCALED DIAGONAL x = 2t  (G_T(t, 2t))")
    print("=" * 70)

    nfact = factorial(n)

    for seed in range(min(5, num_tests)):
        A = random_tournament(n, seed)
        t3, t5, t7, bc = compute_invariants_n7(A, n)

        print(f"\n  seed={seed}: t3={t3}, t5={t5}, t7={t7}, bc={bc}")
        for t_val in [0, Fraction(1,2), 1, Fraction(3,2), 2, -1]:
            g = G_T_formula(t_val, 2*t_val, t3, t5, t7, bc, n)
            label = ""
            if t_val == 0:
                label = "  [= A_n(0) = 1]"
            elif t_val == 1:
                label = f"  [= G_T(1,2) = E_T(1) = n! = {nfact}]"
            print(f"    G_T({t_val}, {2*t_val}) = {g}{label}")

    # Check T-independence
    print("\n  Checking T-independence of G_T(t, 2t) at various t...")
    test_vals = [Fraction(k, 4) for k in range(-4, 13)]
    for t_val in test_vals:
        results = set()
        for seed in range(num_tests):
            A = random_tournament(n, seed)
            t3, t5, t7, bc = compute_invariants_n7(A, n)
            g = G_T_formula(t_val, 2*t_val, t3, t5, t7, bc, n)
            results.add(g)
        if len(results) == 1:
            print(f"    t={t_val}: T-INDEPENDENT! Value = {results.pop()}")

    print()

# ============================================================
# Analysis 3: Anti-diagonal x + t = 2  (passes through (0,2) and (2,0))
# ============================================================

def analysis_anti_diagonal(n=7, num_tests=15):
    print("ANALYSIS 3: ANTI-DIAGONAL x + t = 2  (G_T(t, 2-t))")
    print("=" * 70)
    print(f"  Passes through (t=0, x=2) -> H(T) = I(Omega,2)")
    print(f"  Passes through (t=2, x=0) -> A_n(2)")
    print(f"  Passes through (t=1, x=1) -> G_T(1,1) = n! = {factorial(n)}")

    nfact = factorial(n)

    for seed in range(min(5, num_tests)):
        A = random_tournament(n, seed)
        t3, t5, t7, bc = compute_invariants_n7(A, n)
        H = 1 + (t3 + t5 + t7) * 2 + bc * 4

        print(f"\n  seed={seed}: t3={t3}, t5={t5}, t7={t7}, bc={bc}, H={H}")
        for t_val in [0, Fraction(1,4), Fraction(1,2), Fraction(3,4), 1,
                      Fraction(5,4), Fraction(3,2), Fraction(7,4), 2]:
            x_val = 2 - t_val
            g = G_T_formula(t_val, x_val, t3, t5, t7, bc, n)
            label = ""
            if t_val == 0:
                label = f"  [= H = {H}]"
            elif t_val == 1:
                label = f"  [= n! = {nfact}]"
            elif t_val == 2:
                An2 = eulerian_poly_eval(n, 2)
                label = f"  [= A_{n}(2) = {An2}]"
            print(f"    G_T({t_val}, {x_val}) = {g}{label}")

    # Check T-independence
    print("\n  Checking T-independence of G_T(t, 2-t) at various t...")
    test_vals = [Fraction(k, 4) for k in range(-4, 13)]
    for t_val in test_vals:
        x_val = 2 - t_val
        results = set()
        for seed in range(num_tests):
            A = random_tournament(n, seed)
            t3, t5, t7, bc = compute_invariants_n7(A, n)
            g = G_T_formula(t_val, x_val, t3, t5, t7, bc, n)
            results.add(g)
        if len(results) == 1:
            print(f"    t={t_val}: T-INDEPENDENT! Value = {results.pop()}")

    print()

# ============================================================
# Analysis 4: Level set G_T(t, x) = n! beyond t=1
# ============================================================

def analysis_level_set(n=7, num_tests=10):
    print("ANALYSIS 4: LEVEL SET G_T(t, x) = n!")
    print("=" * 70)

    nfact = factorial(n)

    # At t=1, G_T(1, x) = n! for all x. Is this the ONLY t?
    # Check: for fixed x, is G_T(t, x) = n! only at t=1?

    print(f"\n  n! = {nfact}")
    print(f"\n  For t=1: G_T(1, x) = n! for ALL x (check):")

    A = random_tournament(n, 0)
    t3, t5, t7, bc = compute_invariants_n7(A, n)
    for x_val in [0, 1, 2, 3, -1, 100]:
        g = G_T_formula(1, x_val, t3, t5, t7, bc, n)
        print(f"    G_T(1, {x_val}) = {g}")

    # For fixed x=0: G_T(t, 0) = A_n(t). Does A_n(t) = n! have t=1 as root?
    # A_n(1) = n!, yes. Any other roots?
    print(f"\n  A_n(t) = n! roots (fixed x=0, T-independent):")
    An_coeffs = [eulerian_number(n, k) for k in range(n)]
    # A_n(t) - n! = 0 at t=1 by definition
    # A_n(t) = sum A(n,k) t^k, and sum A(n,k) = n!
    # So A_n(t) - n! = sum A(n,k) (t^k - 1) = (t-1) * Q(t)
    print(f"    A_{n} coeffs: {An_coeffs}")
    print(f"    A_{n}(1) = {sum(An_coeffs)} = {nfact}")

    # For fixed x=2: G_T(t, 2) = E_T(t). Does E_T(t) = n! only at t=1?
    print(f"\n  For x=2: E_T(t) = n! solutions (T-dependent)")
    for seed in range(min(5, num_tests)):
        A = random_tournament(n, seed)
        t3, t5, t7, bc = compute_invariants_n7(A, n)
        dist = forward_edge_dist_dp(A, n)

        # E_T(t) - n! vanishes at t=1. Find other roots by dividing by (t-1).
        # E_T(t) = sum a_k t^k. Coefficients:
        coeffs = [dist.get(k, 0) for k in range(n)]
        # Subtract n! distributed: E_T(t) - n! = sum a_k t^k - n!
        # = (a_0 - n!) + a_1 t + ... + a_{n-1} t^{n-1}
        shifted = list(coeffs)
        shifted[0] -= nfact

        # Divide by (t-1) using synthetic division
        # If P(t) = c_0 + c_1 t + ... + c_{n-1} t^{n-1} and P(1)=0,
        # P(t)/(t-1) = q_{n-2} t^{n-2} + ... + q_0
        # Use Horner from top: q_{n-2} = c_{n-1}, q_{k-1} = q_k + c_k
        q = [0] * (n - 1)
        q[n-2] = shifted[n-1]
        for k in range(n-2, 0, -1):
            q[k-1] = q[k] + shifted[k]
        # Now E_T(t) - n! = (t-1) * Q(t) where Q(t) = sum q_k t^k

        # Evaluate Q at various points to check if it has real roots
        print(f"\n  seed={seed}: coeffs={coeffs}")
        print(f"    Q(t) coeffs (after dividing E_T-n! by t-1): {q}")
        for t_val in [Fraction(k, 2) for k in range(-4, 9)]:
            Q_val = sum(q[k] * t_val**k for k in range(n-1))
            if Q_val == 0:
                print(f"    *** Q({t_val}) = 0 => E_T({t_val}) = n! ***")

    # More generally: for what curves (t, x(t)) is G_T = n! for ALL T?
    # G_T(t, x) = A_n(t) + sum_I x^{parts} I(T) A_{f+1}(t) (t-1)^{d-f} = n!
    # This requires A_n(t) = n! (handled by t=1) AND all invariant terms = 0,
    # OR the invariant terms sum to n! - A_n(t) in a T-independent way (impossible
    # since they depend on T linearly).
    # So for ALL T: need each coefficient of I(T) to vanish:
    #   x * A_5(t) * (t-1)^2 = 0  (coeff of t3)
    #   x * A_3(t) * (t-1)^4 = 0  (coeff of t5)
    #   x * A_1(t) * (t-1)^6 = 0  (coeff of t7)
    #   x^2 * A_3(t) * (t-1)^4 = 0 (coeff of bc)
    # These all vanish iff x=0 or t=1.

    print(f"\n  CONCLUSION: G_T(t, x) = n! for ALL T iff (t=1) or (x=0 and A_n(t)=n!)")
    print(f"  Since A_n(t)=n! iff t=1 (Eulerian poly has all positive coeffs),")
    print(f"  the ONLY universal solution is t=1.")

    # But wait: A_n(t) = n! could have other solutions. Let's check.
    print(f"\n  Checking A_{n}(t) = n! for real t:")
    for t_val in [Fraction(k, 10) for k in range(-20, 31)]:
        if eulerian_poly_eval(n, t_val) == nfact:
            print(f"    A_{n}({t_val}) = {nfact} = n!")

    print()

# ============================================================
# Analysis 5: G_T(t, -1) — independence polynomial at x=-1
# ============================================================

def analysis_x_neg1(n=7, num_tests=15):
    print("ANALYSIS 5: G_T(t, -1)")
    print("=" * 70)

    nfact = factorial(n)

    # G_T(t, -1) = A_n(t) - t3*A_5(t)*(t-1)^2 - t5*A_3(t)*(t-1)^4
    #              - t7*A_1(t)*(t-1)^6 + bc*A_3(t)*(t-1)^4
    # Note: x=-1, so x^2 = 1, and x^1 = -1

    # At t=0: G_T(0, -1) = I(Omega, -1) = 1 - alpha_1 + alpha_2
    print("  G_T(0, -1) = I(Omega(T), -1) = 1 - (t3+t5+t7) + bc")

    for seed in range(min(8, num_tests)):
        A = random_tournament(n, seed)
        t3, t5, t7, bc = compute_invariants_n7(A, n)
        alpha1 = t3 + t5 + t7
        alpha2 = bc

        I_neg1 = 1 - alpha1 + alpha2

        print(f"\n  seed={seed}: t3={t3}, t5={t5}, t7={t7}, bc={bc}")
        print(f"    I(Omega, -1) = {I_neg1}")

        for t_val in [0, Fraction(1,2), 1, 2, -1]:
            g = G_T_formula(t_val, -1, t3, t5, t7, bc, n)
            label = ""
            if t_val == 0:
                label = f"  [= I(Omega,-1) = {I_neg1}]"
            elif t_val == 1:
                label = f"  [= n! = {nfact}]"
            print(f"    G_T({t_val}, -1) = {g}{label}")

    # Check: is I(Omega, -1) related to any known quantity?
    # I(Omega, -1) = chromatic polynomial at 1? (No, that's for chromatic poly.)
    # For interval graphs, I(G, -1) = (-1)^alpha(G) if G is well-covered? Not sure.

    # T-independence check
    print("\n  Checking T-independence of G_T(t, -1) at various t...")
    test_vals = [Fraction(k, 4) for k in range(-4, 13)]
    for t_val in test_vals:
        results = set()
        for seed in range(num_tests):
            A = random_tournament(n, seed)
            t3, t5, t7, bc = compute_invariants_n7(A, n)
            g = G_T_formula(t_val, -1, t3, t5, t7, bc, n)
            results.add(g)
        if len(results) == 1:
            print(f"    t={t_val}: T-INDEPENDENT! Value = {results.pop()}")

    print()

# ============================================================
# Analysis 6: Curve (t-1)*x = constant
# ============================================================

def analysis_hyperbola(n=7, num_tests=15):
    print("ANALYSIS 6: HYPERBOLA (t-1)*x = c")
    print("=" * 70)

    nfact = factorial(n)

    # At (t-1)*x = c, we have x = c/(t-1).
    # G_T(t, c/(t-1)) = A_n(t) + sum_I [c/(t-1)]^{parts(I)} I(T) A_{f+1}(t) (t-1)^{d-f}
    # For parts=1 terms: [c/(t-1)] * A_{f+1}(t) * (t-1)^{d-f} = c * A_{f+1}(t) * (t-1)^{d-f-1}
    # For parts=2 terms: [c/(t-1)]^2 * A_{f+1}(t) * (t-1)^{d-f} = c^2 * A_{f+1}(t) * (t-1)^{d-f-2}

    # Special case c=2: x = 2/(t-1)
    # At t=0: x = -2 -> G_T(0, -2) = I(Omega, -2) = 1 - 2*alpha1 + 4*alpha2
    # At t=-1: x = -1 -> G_T(-1, -1)
    # At t=2: x = 2 -> G_T(2, 2)
    # At t=3: x = 1 -> G_T(3, 1)

    for c_val in [Fraction(1), Fraction(2), Fraction(-1), Fraction(-2)]:
        print(f"\n  (t-1)*x = {c_val}:")

        for seed in range(min(3, num_tests)):
            A = random_tournament(n, seed)
            t3, t5, t7, bc = compute_invariants_n7(A, n)

            print(f"    seed={seed}: t3={t3}, t5={t5}, t7={t7}, bc={bc}")
            for t_val in [Fraction(-1), Fraction(0), Fraction(2), Fraction(3)]:
                if t_val == 1:
                    continue  # x = c/0, undefined
                x_val = c_val / (t_val - 1)
                g = G_T_formula(t_val, x_val, t3, t5, t7, bc, n)
                print(f"      G_T({t_val}, {x_val}) = {g}")

    # T-independence check for (t-1)*x = c at various t
    print("\n  Checking T-independence along (t-1)*x = c ...")
    for c_val in [Fraction(1), Fraction(2), Fraction(-1), Fraction(-2)]:
        for t_val in [Fraction(k, 4) for k in range(-4, 13) if k != 4]:  # skip t=1
            x_val = c_val / (t_val - 1)
            results = set()
            for seed in range(num_tests):
                A = random_tournament(n, seed)
                t3, t5, t7, bc = compute_invariants_n7(A, n)
                g = G_T_formula(t_val, x_val, t3, t5, t7, bc, n)
                results.add(g)
            if len(results) == 1:
                print(f"    c={c_val}, t={t_val}: T-INDEPENDENT! Value = {results.pop()}")

    print()

# ============================================================
# Analysis 7: Summary — linear dependence structure
# ============================================================

def analysis_summary(n=7, num_tests=15):
    print("ANALYSIS 7: STRUCTURE SUMMARY")
    print("=" * 70)

    nfact = factorial(n)

    # G_T(t,x) depends on T only through (t3, t5, t7, bc).
    # The coefficients are:
    # coeff of t3:  x * A_5(t) * (t-1)^2
    # coeff of t5:  x * A_3(t) * (t-1)^4
    # coeff of t7:  x * A_1(t) * (t-1)^6   [A_1(t) = 1]
    # coeff of bc:  x^2 * A_3(t) * (t-1)^4

    # Note t5 and bc share A_3(t)*(t-1)^4 factor!
    # So G_T(t,x) = A_n(t) + x*(t-1)^2*A_5(t)*t3 + x*(t-1)^4*A_3(t)*t5
    #             + x*(t-1)^6*t7 + x^2*(t-1)^4*A_3(t)*bc
    # = A_n(t) + x*(t-1)^2*A_5(t)*t3 + (t-1)^4*A_3(t)*(x*t5 + x^2*bc) + x*(t-1)^6*t7

    # When does the coefficient of t5 equal the coefficient of bc?
    # x * A_3(t) * (t-1)^4 = x^2 * A_3(t) * (t-1)^4
    # => x = x^2 => x(x-1) = 0 => x=0 or x=1
    # At x=1: t5 and bc have the SAME coefficient!

    print("  Coefficient structure of G_T(t, x) for n=7:")
    print(f"    coeff(t3) = x * A_5(t) * (t-1)^2")
    print(f"    coeff(t5) = x * A_3(t) * (t-1)^4")
    print(f"    coeff(t7) = x * (t-1)^6")
    print(f"    coeff(bc) = x^2 * A_3(t) * (t-1)^4")
    print()
    print(f"  Note: coeff(t5)/coeff(bc) = 1/x, so at x=1 they are equal!")
    print(f"  At x=1: G_T(t,1) depends on t5+bc (not separately).")
    print()

    # Values at special points on the diagonal
    print("  G_T(t, t) at t=0: G_T(0,0) = A_n(0) = 1 for all T")
    print(f"  G_T(t, t) at t=1: G_T(1,1) = n! = {nfact} for all T")
    print()

    # G_T(t, 1) — at x=1, the "reduced" generating function
    print("  G_T(t, 1): the 'reduced' generating function")
    for seed in range(min(5, num_tests)):
        A = random_tournament(n, seed)
        t3, t5, t7, bc = compute_invariants_n7(A, n)
        H = 1 + (t3+t5+t7)*2 + bc*4

        vals = []
        for t_val in [0, Fraction(1, 2), 1, 2]:
            g = G_T_formula(t_val, 1, t3, t5, t7, bc, n)
            vals.append((t_val, g))

        print(f"    seed={seed}: t3={t3}, t5+bc={t5+bc}, t7={t7}, H={H}")
        for t_val, g in vals:
            print(f"      G_T({t_val}, 1) = {g}")

    print()

# ============================================================
# BONUS: Check palindromy / symmetry under t <-> 1/t or similar
# ============================================================

def analysis_symmetry(n=7, num_tests=10):
    print("ANALYSIS 8: SYMMETRY CHECKS")
    print("=" * 70)

    nfact = factorial(n)

    # E_T(t) = sum a_k t^k has palindromic coefficients: a_k = a_{n-1-k}.
    # So E_T(t) = t^{n-1} E_T(1/t).
    # Does G_T(t, x) have any functional equation?
    # G_T(t, x) = A_n(t) + sum ...
    # A_n(t) = t^{n-1} A_n(1/t).
    # A_{f+1}(t) = t^f A_{f+1}(1/t).
    # (t-1)^{d-f} = t^{d-f} (1-1/t)^{d-f} = t^{d-f} (-(1/t-1))^{d-f} = (-1)^{d-f} t^{d-f} (1/t-1)^{d-f}

    # So each term:
    # x^p * I * A_{f+1}(t) * (t-1)^{d-f}
    # = x^p * I * t^f * A_{f+1}(1/t) * (-1)^{d-f} * t^{d-f} * (1/t-1)^{d-f}
    # = x^p * I * (-1)^{d-f} * t^{f+d-f} * A_{f+1}(1/t) * (1/t-1)^{d-f}
    # = x^p * I * (-1)^{d-f} * t^d * A_{f+1}(1/t) * (1/t-1)^{d-f}

    # So G_T(t, x) = t^d * G_T(1/t, (-1)^?? * x)?
    # The sign factor (-1)^{d-f} depends on f, so no uniform substitution works.
    # Unless x -> -x absorbs it... let's check.

    # For n=7, d=6:
    #   t3: d-f = 2 (even), factor = +1
    #   t5: d-f = 4 (even), factor = +1
    #   t7: d-f = 6 (even), factor = +1
    #   bc: d-f = 4 (even), factor = +1
    # All even! So for n=7:
    # G_T(t, x) = t^6 * G_T(1/t, x)

    print(f"  For n={n}, d={n-1}:")
    print(f"  d-f values: t3: {n-1-4}=2, t5: {n-1-2}=4, t7: {n-1-0}=6, bc: {n-1-2}=4")
    print(f"  ALL EVEN! Testing G_T(t, x) = t^{n-1} * G_T(1/t, x)...")

    for seed in range(min(5, num_tests)):
        A = random_tournament(n, seed)
        t3, t5, t7, bc = compute_invariants_n7(A, n)

        ok = True
        for t_val in [Fraction(1, 2), Fraction(1, 3), Fraction(2, 3), 2, 3, Fraction(-1,1)]:
            for x_val in [0, 1, 2, -1, Fraction(1, 2)]:
                g1 = G_T_formula(t_val, x_val, t3, t5, t7, bc, n)
                g2 = t_val**(n-1) * G_T_formula(Fraction(1, t_val), x_val, t3, t5, t7, bc, n)
                if g1 != g2:
                    ok = False
                    print(f"    FAIL: seed={seed}, t={t_val}, x={x_val}: {g1} != {g2}")
        if ok:
            print(f"    seed={seed}: G_T(t,x) = t^{n-1} G_T(1/t, x) CONFIRMED")

    # What about general n? d-f values:
    # For a k-cycle: f = n-1-k+1 = n-k, so d-f = (n-1)-(n-k) = k-1.
    # k is odd (3,5,7,...), so k-1 is even!
    # For bc (pair of 3-cycles): f = n-1-6+2 = n-5, d-f = 4 (even for any n).
    # More generally, parts(I)=p, with vertex set of size s, f = n-1-s+p.
    # d-f = (n-1)-(n-1-s+p) = s-p.
    # s = sum of cycle lengths (all odd), p = number of cycles.
    # s-p = sum(L_i - 1) = sum of even numbers = EVEN.
    # So (-1)^{d-f} = 1 ALWAYS!

    print(f"\n  GENERAL ARGUMENT: s-p = sum(L_i - 1) is always even since each L_i is odd.")
    print(f"  Therefore: G_T(t, x) = t^{{n-1}} * G_T(1/t, x) FOR ALL n!")
    print(f"  This is a PALINDROMIC symmetry in t (with x fixed).")
    print()

    # Implication: G_T(t, x) as a polynomial in t of degree n-1 is palindromic.
    # So its roots come in pairs r, 1/r (and t=-1 or t=1 could be self-paired).
    # At x=0: this is the Eulerian palindromy. At x=2: this is the known a_k = a_{d-k}.
    # This generalizes to ALL x!

    print(f"  IMPLICATION: For every fixed x, G_T(t, x) is a palindromic polynomial in t.")
    print(f"  The forward-edge-like distribution is symmetric for ALL x, not just x=2.")
    print()

# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    verify_formula()
    analysis_diagonal()
    analysis_scaled_diagonal()
    analysis_anti_diagonal()
    analysis_level_set()
    analysis_x_neg1()
    analysis_hyperbola()
    analysis_summary()
    analysis_symmetry()
