#!/usr/bin/env python3
"""
Roots of W(r) — the weighted path polynomial.

W(r) = sum_P prod_{i=0}^{n-2} (r + s_i) where s_i in {+1/2, -1/2}.

For odd n: W has only odd powers of r (from r-parity).
For even n: W has only even powers of r.

So we can study the "reduced" polynomial:
  Odd n: W(r) = r * V(r^2) where V is degree (n-2)/2
  Even n: W(r) = U(r^2) where U is degree (n-1)/2

Question: does V(x) or U(x) have all real negative roots?

Also: what is the relationship between W(r)'s roots and the OCF invariants?

opus-2026-03-07-S32
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
import random

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

def compute_W_coeffs(A, n):
    """Compute W(r) as polynomial in r, return coefficient list [w_0, w_1, ..., w_{n-1}]."""
    # DP approach: state = (mask, last_vertex), value = polynomial coefficients
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = [1.0]  # constant polynomial = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            poly = dp.get((mask, v))
            if poly is None: continue
            for u in range(n):
                if mask & (1 << u): continue
                s = A[v][u] - 0.5  # +0.5 or -0.5
                # Multiply poly by (r + s): new_poly = r*poly + s*poly
                new_poly_len = len(poly) + 1
                key = (mask | (1 << u), u)
                if key not in dp:
                    dp[key] = [0.0] * new_poly_len
                existing = dp[key]
                while len(existing) < new_poly_len:
                    existing.append(0.0)
                for i in range(len(poly)):
                    existing[i] += s * poly[i]      # s * r^i term
                    existing[i+1] += poly[i]         # r * r^i = r^{i+1} term

    full = (1 << n) - 1
    result = [0.0] * n
    for v in range(n):
        poly = dp.get((full, v), [])
        for i in range(min(len(poly), n)):
            result[i] += poly[i]
    return result

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

print("ROOTS OF W(r) FOR TOURNAMENTS")
print("=" * 70)

for n in [5, 7, 9]:
    print(f"\nn={n} ({'odd' if n % 2 else 'even'} — W has only {'odd' if n % 2 else 'even'} powers):")

    all_real = 0
    complex_ct = 0
    all_neg_sq = 0  # all roots of V(x) are real negative
    tested = 50 if n <= 7 else 20

    for trial in range(tested):
        A = random_tournament(n, n * 444 + trial)
        coeffs = compute_W_coeffs(A, n)
        t3 = count_t3(A, n)

        # Parity: W(-r) = (-1)^{n-1} W(r)
        # Odd n (n-1 even): W is EVEN, only even powers. W(r) = U(r^2).
        # Even n (n-1 odd): W is ODD, only odd powers. W(r) = r * V(r^2).
        if n % 2 == 1:
            # W is even: extract even coefficients
            v_coeffs = [coeffs[2*i] for i in range((n+1)//2)]
        else:
            # W is odd: extract odd coefficients
            v_coeffs = [coeffs[2*i+1] for i in range(n//2)]

        # Find roots of V(x) = v_coeffs[0] + v_coeffs[1]*x + ...
        if len(v_coeffs) > 1:
            np_coeffs = list(reversed(v_coeffs))  # highest degree first for np.roots
            roots = np.roots(np_coeffs)

            is_all_real = all(abs(r.imag) < 1e-4 for r in roots)
            is_all_neg = is_all_real and all(r.real < 1e-4 for r in roots)

            if is_all_real:
                all_real += 1
            else:
                complex_ct += 1

            if is_all_neg:
                all_neg_sq += 1

            if trial < 5:
                root_str = ", ".join(f"{r.real:.3f}" if abs(r.imag) < 1e-4
                                     else f"{r.real:.3f}+{r.imag:.3f}i" for r in sorted(roots, key=lambda r: r.real))
                print(f"  T{trial}: t3={t3:3d}, V-roots=[{root_str}], real={is_all_real}, neg={is_all_neg}")

    print(f"  Summary: {all_real}/{tested} V real-rooted, {all_neg_sq}/{tested} all-negative")

# Now check: is W(r) itself real-rooted for all r?
print(f"\n{'=' * 70}")
print("W(r) ROOTS (full polynomial)")
print("=" * 70)

for n in [5, 7]:
    print(f"\nn={n}:")
    all_pure_imag = 0
    tested = 50

    for trial in range(tested):
        A = random_tournament(n, n * 555 + trial)
        coeffs = compute_W_coeffs(A, n)
        np_coeffs = list(reversed(coeffs))
        roots = np.roots(np_coeffs)

        # For odd n with only odd powers: roots should be 0 and pairs ±ri (pure imaginary)
        # Actually: W(r) = r*V(r^2), so roots are 0 and ±sqrt(v_root) for each root of V
        # If V has all negative roots v_k < 0, then sqrt(v_k) = ±i*sqrt(|v_k|) — pure imaginary!
        # So W real-rooted iff V has all real negative roots.

        is_pure_imag = all(abs(r.real) < 1e-4 or abs(r.imag) < 1e-4 for r in roots)

        if is_pure_imag:
            all_pure_imag += 1

        if trial < 3:
            print(f"  T{trial}: roots={sorted(roots, key=lambda r: (abs(r.imag), r.real))}")

    print(f"  {all_pure_imag}/{tested} have all roots on real or imaginary axis")

# Key connection: W(r) roots are on imaginary axis iff V(x) has all real negative roots
# V(x) negative roots => W(r) = r * V(r^2) has roots at r=0 and r=±i*sqrt(|v_k|)
# These are pure imaginary! So the question is about V(r^2) being real-rooted with negative roots.

print(f"\n{'=' * 70}")
print("CONCLUSION")
print("=" * 70)
print("""
For odd n, W(r) = r * V(r^2) where V is the "reduced" polynomial.
V(x) having all real negative roots is equivalent to W(r) having all roots
on the imaginary axis (plus the root at 0).

This is NOT always true — V sometimes has complex roots.
The deformation from Eulerian numbers (via OCF invariants) can push
V's roots off the real line.
""")
