#!/usr/bin/env python3
"""
gk_deep_cascade_89c.py — Deep cascading extraction using polynomial assumptions
opus-2026-03-14-S89c

Strategy: assume g_k is a polynomial of degree d_k.
Use the assumption to predict values, then extract the next g_k.
Integer consistency at each step validates the assumed degree.

Known:
  g_1(m) = m                                    (deg 1)
  g_2(m) = m²                                   (deg 2)
  g_3(m) = m(2m²+1)/3                           (deg 3)
  g_4(m) = (10m³-33m²+50m-24)/3                 (deg 3)
  g_5(m) = m(74m³-352m²+550m-269)/3             (deg 4)
  g_6(m) = m(9110m³-48080m²+80485m-41512)/3     (deg 4, needs verification)
"""

from fractions import Fraction
from math import factorial
from functools import reduce
from sympy import symbols, interpolate, Rational, factor, expand
import time

x = symbols('x')

def falling_factorial(n, k):
    return reduce(lambda a, b: a * b, range(n, n - k, -1), 1)

def compute_W_dp(n):
    full = (1 << n) - 1
    dp = {}
    for v in range(n):
        dp[((1 << v), v)] = 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp.get((mask, v), 0)
            if cnt == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if u == v - 1:
                    continue
                weight = 2 * cnt if u == v + 1 else cnt
                new_key = (mask | (1 << u), u)
                dp[new_key] = dp.get(new_key, 0) + weight
    return sum(dp.get((full, v), 0) for v in range(n))

# Known polynomial formulas (exact)
g_polys = {}
g_polys[1] = interpolate([(Rational(0), 0), (Rational(1), 1)], x)
g_polys[2] = interpolate([(Rational(0), 0), (Rational(1), 1), (Rational(2), 4)], x)
g_polys[3] = interpolate([(Rational(0), 0), (Rational(1), 1), (Rational(2), 6),
                           (Rational(3), 19)], x)
g_polys[4] = interpolate([(Rational(0), 0), (Rational(1), 1), (Rational(2), 8),
                           (Rational(3), 33), (Rational(4), 96)], x)
g_polys[5] = interpolate([(Rational(0), 0), (Rational(1), 1), (Rational(2), 10),
                           (Rational(3), 211), (Rational(4), 1380), (Rational(5), 4885)], x)

def g_eval(k, m):
    """Evaluate g_k(m) using polynomial formula."""
    if k in g_polys:
        return Fraction(int(g_polys[k].subs(x, m)))
    return None

# Print known polynomials
print("="*70)
print("KNOWN g_k POLYNOMIALS")
print("="*70)
for k in sorted(g_polys):
    p = g_polys[k]
    print(f"g_{k}(m) = {factor(p)} = {expand(p)}  (degree {p.as_poly().degree()})")

# Compute W(n)
print("\nComputing W(n)...")
W = {}
for n in range(3, 22):
    t0 = time.time()
    W[n] = compute_W_dp(n)
    t1 = time.time()
    print(f"  W({n}) = {W[n]}  ({t1-t0:.1f}s)")

# Deep cascade extraction
print("\n" + "="*70)
print("DEEP CASCADE: extracting g_6, g_7, g_8, ...")
print("="*70)

gk_data = {k: {} for k in range(6, 12)}
# Boundary conditions
for k in range(6, 12):
    gk_data[k][0] = Fraction(0)
    gk_data[k][1] = Fraction(1)
    gk_data[k][2] = Fraction(2 * k)  # g_k(2) = 2k

max_known_k = 5  # highest k with polynomial formula

for n in sorted(W):
    nf = factorial(n)
    cv2 = Fraction(W[n], nf) - 1
    remaining = cv2

    # Subtract all known polynomial g_k (k=1..max_known_k)
    for k in range(1, max_known_k + 1):
        mk = n - 2*k
        if mk < 1:
            break
        remaining -= Fraction(2 * g_eval(k, mk), falling_factorial(n, 2*k))

    # Now remaining = Σ_{k≥max_known_k+1} 2·g_k(n-2k)/(n)_{2k}
    # For each k ≥ max_known_k+1 with n-2k ≥ 1
    terms = []
    for k in range(max_known_k + 1, (n-1)//2 + 1):
        mk = n - 2*k
        if mk < 0:
            break
        terms.append((k, mk))

    if not terms:
        continue

    # Subtract known g_k values (from boundary conditions)
    unknowns = []
    for k, mk in terms:
        if mk in gk_data.get(k, {}):
            remaining -= Fraction(2 * gk_data[k][mk], falling_factorial(n, 2*k))
        else:
            unknowns.append((k, mk))

    if not unknowns:
        if remaining != 0:
            print(f"n={n}: RESIDUAL = {remaining} = {float(remaining):.6e}")
        continue

    if len(unknowns) == 1:
        k, mk = unknowns[0]
        val = remaining * falling_factorial(n, 2*k) / 2
        is_int = val.denominator == 1
        print(f"n={n}: g_{k}({mk}) = {val}{'' if is_int else f' ({float(val):.4f}) NOT INT!'}")
        gk_data[k][mk] = val
    else:
        # Multiple unknowns — record but can't resolve
        print(f"n={n}: {len(unknowns)} unknowns: {unknowns}")

# Try to fit g_6 polynomial
print("\n" + "="*70)
print("g_6 POLYNOMIAL FIT")
print("="*70)

pts6 = [(Rational(m), Rational(gk_data[6][m])) for m in sorted(gk_data[6])
        if gk_data[6][m].denominator == 1]
print(f"Integer data points: {len(pts6)}")
for m, v in pts6:
    print(f"  g_6({m}) = {v}")

if len(pts6) >= 5:
    for deg in range(1, len(pts6)):
        poly6 = interpolate(pts6[:deg+1], x)
        all_match = True
        for m_val, g_val in pts6[deg+1:]:
            if poly6.subs(x, m_val) != g_val:
                all_match = False
                break
        if all_match and len(pts6) > deg + 1:
            print(f"\n→ g_6(m) = {factor(poly6)} = {expand(poly6)} (EXACT degree {deg})")
            g_polys[6] = poly6
            break
    else:
        poly6 = interpolate(pts6, x)
        print(f"\n→ g_6(m) = {factor(poly6)} (through {len(pts6)} points)")
        g_polys[6] = poly6

# If we got g_6, redo extraction for g_7
if 6 in g_polys:
    print("\n" + "="*70)
    print("SECOND PASS: extracting g_7 using g_6 polynomial")
    print("="*70)

    max_known_k = 6
    gk_data[7] = {0: Fraction(0), 1: Fraction(1), 2: Fraction(14)}

    for n in sorted(W):
        nf = factorial(n)
        cv2 = Fraction(W[n], nf) - 1
        remaining = cv2

        for k in range(1, max_known_k + 1):
            mk = n - 2*k
            if mk < 1:
                break
            remaining -= Fraction(2 * g_eval(k, mk), falling_factorial(n, 2*k))

        terms = []
        for k in range(max_known_k + 1, (n-1)//2 + 1):
            mk = n - 2*k
            if mk < 0:
                break
            terms.append((k, mk))

        if not terms:
            continue

        unknowns = []
        for k, mk in terms:
            if mk in gk_data.get(k, {}):
                remaining -= Fraction(2 * gk_data[k][mk], falling_factorial(n, 2*k))
            else:
                unknowns.append((k, mk))

        if not unknowns:
            if remaining != 0:
                print(f"n={n}: RESIDUAL = {remaining} = {float(remaining):.6e}")
            continue

        if len(unknowns) == 1:
            k, mk = unknowns[0]
            val = remaining * falling_factorial(n, 2*k) / 2
            is_int = val.denominator == 1
            print(f"n={n}: g_{k}({mk}) = {val}{'' if is_int else f' ({float(val):.4f}) NOT INT!'}")
            gk_data[k][mk] = val

    # Fit g_7
    pts7 = [(Rational(m), Rational(gk_data[7][m])) for m in sorted(gk_data[7])
            if gk_data[7][m].denominator == 1]
    print(f"\ng_7 integer data: {len(pts7)} points")
    for m, v in pts7:
        print(f"  g_7({m}) = {v}")

    if len(pts7) >= 5:
        for deg in range(1, len(pts7)):
            poly7 = interpolate(pts7[:deg+1], x)
            all_match = True
            for m_val, g_val in pts7[deg+1:]:
                if poly7.subs(x, m_val) != g_val:
                    all_match = False
                    break
            if all_match and len(pts7) > deg + 1:
                print(f"\n→ g_7(m) = {factor(poly7)} = {expand(poly7)} (EXACT degree {deg})")
                g_polys[7] = poly7
                break
        else:
            poly7 = interpolate(pts7, x)
            print(f"\n→ g_7(m) = {factor(poly7)} (through {len(pts7)} points)")
            g_polys[7] = poly7

# Third pass for g_8
if 7 in g_polys:
    print("\n" + "="*70)
    print("THIRD PASS: extracting g_8 using g_7 polynomial")
    print("="*70)

    max_known_k = 7
    gk_data[8] = {0: Fraction(0), 1: Fraction(1), 2: Fraction(16)}

    for n in sorted(W):
        nf = factorial(n)
        cv2 = Fraction(W[n], nf) - 1
        remaining = cv2

        for k in range(1, max_known_k + 1):
            mk = n - 2*k
            if mk < 1:
                break
            remaining -= Fraction(2 * g_eval(k, mk), falling_factorial(n, 2*k))

        terms = []
        for k in range(max_known_k + 1, (n-1)//2 + 1):
            mk = n - 2*k
            if mk < 0:
                break
            terms.append((k, mk))

        if not terms:
            continue

        unknowns = []
        for k, mk in terms:
            if mk in gk_data.get(k, {}):
                remaining -= Fraction(2 * gk_data[k][mk], falling_factorial(n, 2*k))
            else:
                unknowns.append((k, mk))

        if not unknowns:
            continue
        if len(unknowns) == 1:
            k, mk = unknowns[0]
            val = remaining * falling_factorial(n, 2*k) / 2
            is_int = val.denominator == 1
            print(f"n={n}: g_{k}({mk}) = {val}{'' if is_int else f' ({float(val):.4f}) NOT INT!'}")
            gk_data[k][mk] = val

    pts8 = [(Rational(m), Rational(gk_data[8][m])) for m in sorted(gk_data[8])
            if gk_data[8][m].denominator == 1]
    print(f"\ng_8 integer data: {len(pts8)} points")
    for m, v in pts8:
        print(f"  g_8({m}) = {v}")

    if len(pts8) >= 5:
        poly8 = interpolate(pts8, x)
        print(f"\n→ g_8(m) = {factor(poly8)} (through {len(pts8)} points)")

# Final summary
print("\n" + "="*70)
print("SUMMARY: g_k POLYNOMIALS")
print("="*70)
for k in sorted(g_polys):
    p = g_polys[k]
    deg = p.as_poly().degree()
    print(f"g_{k}(m) = {factor(p)}  (degree {deg})")

# Degree sequence
print("\nDegree sequence:")
for k in sorted(g_polys):
    deg = g_polys[k].as_poly().degree()
    print(f"  k={k}: degree {deg}")

# Leading coefficient sequence
print("\nLeading coefficients:")
for k in sorted(g_polys):
    p = g_polys[k].as_poly()
    lc = p.LC()
    print(f"  k={k}: leading coeff = {lc}")

# g_k(1) = 1 check
print("\ng_k(1) check:")
for k in sorted(g_polys):
    val = g_polys[k].subs(x, 1)
    print(f"  g_{k}(1) = {val} {'✓' if val == 1 else '✗'}")

# g_k(2) = 2k check
print("\ng_k(2) check:")
for k in sorted(g_polys):
    val = g_polys[k].subs(x, 2)
    print(f"  g_{k}(2) = {val}, expected {2*k} {'✓' if val == 2*k else '✗'}")

print("\nDone!")
