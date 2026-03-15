#!/usr/bin/env python3
"""
gk_vs_thm201_s112.py — Compare THM-201 (original + S111d correction) vs THM-216 (exact g_k)

Purpose:
  1. Compare g_k(m) vs m^k (THM-201 original) and m^k*[1-((m-1)/m)^m] (S111d correction)
  2. Use exact g_k coefficients from THM-216 to compute g_k(m) for k=1..9
  3. Compute ratio g_k(m)/m^k for various m, k
  4. Search for a BETTER correction factor that works exactly

kind-pasteur-2026-03-15-S112
"""

from fractions import Fraction
from math import factorial, comb, exp, log
from functools import reduce


def falling_factorial(n, k):
    return reduce(lambda a, b: a * b, range(n, n - k, -1), 1)


# --------------------------------------------------------------------------
# Exact g_k polynomials from THM-216
# For k=1,2: g_k(m) directly.
# For k>=3: 3*g_k(m) = a*m^3 + b*m^2 + c*m + d
# --------------------------------------------------------------------------
gk_coeffs = {
    1: (0, 0, 3, 0),       # g_1(m) = m
    2: (0, 3, 0, 0),       # g_2(m) = m^2
    3: (2, 0, 1, 0),       # 3*g_3 = 2m^3 + m
    4: (10, -33, 50, -24),  # 3*g_4 = 10m^3 - 33m^2 + 50m - 24
    5: (388, -2040, 3431, -1776),
    6: (69660, -380445, 653748, -342960),
    7: (19826270, -109486152, 189674605, -100014720),
    8: (7309726742, -40641958545, 70757788486, -37425556680),
    9: (3262687720240, -18232387983408, 31858349908595, -16888649645424),
}


def g_exact(k, m):
    """Exact g_k(m) as a Fraction."""
    a, b, c, d = gk_coeffs[k]
    return Fraction(a * m**3 + b * m**2 + c * m + d, 3)


# --------------------------------------------------------------------------
# THM-201 formulas
# --------------------------------------------------------------------------
def thm201_original(k, m):
    """THM-201 original: g_k(m) ~ m^k"""
    return Fraction(m**k)


def thm201_s111d(k, m):
    """THM-201 S111d correction: m^k * [1 - ((m-1)/m)^m] for k>=3, m>=2.
    For k<=2 or m=1: just m^k."""
    if k <= 2 or m <= 1:
        return Fraction(m**k)
    correction = 1 - Fraction(m - 1, m)**m
    return Fraction(m**k) * correction


# --------------------------------------------------------------------------
# PART 1: Direct comparison tables
# --------------------------------------------------------------------------
print("=" * 90)
print("PART 1: g_k(m) vs m^k [THM-201 original] vs m^k*[1-((m-1)/m)^m] [S111d correction]")
print("=" * 90)

for k in range(1, 10):
    print(f"\n--- k = {k} ---")
    if k <= 2:
        print(f"  g_{k}(m) = {'m' if k == 1 else 'm^2'}")
    else:
        a, b, c, d = gk_coeffs[k]
        print(f"  3*g_{k}(m) = {a}*m^3 + ({b})*m^2 + ({c})*m + ({d})")

    print(f"  {'m':>4}  {'g_k(m) exact':>20}  {'m^k':>15}  {'s111d':>20}  "
          f"{'g/m^k':>12}  {'g/s111d':>12}  {'err_orig':>12}  {'err_s111d':>12}")
    print(f"  {'-'*4}  {'-'*20}  {'-'*15}  {'-'*20}  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*12}")

    for m in range(1, 11):
        g = g_exact(k, m)
        orig = thm201_original(k, m)
        s111 = thm201_s111d(k, m)

        ratio_orig = float(g) / float(orig) if orig != 0 else float('inf')
        ratio_s111 = float(g) / float(s111) if s111 != 0 else float('inf')
        err_orig = float(g - orig) / float(orig) if orig != 0 else float('inf')
        err_s111 = float(g - s111) / float(s111) if s111 != 0 else float('inf')

        print(f"  {m:4d}  {str(g):>20s}  {str(orig):>15s}  {str(s111):>20s}  "
              f"{ratio_orig:12.8f}  {ratio_s111:12.8f}  {err_orig:12.6e}  {err_s111:12.6e}")


# --------------------------------------------------------------------------
# PART 2: The exact ratio g_k(m)/m^k — does it have a nice closed form?
# --------------------------------------------------------------------------
print("\n\n" + "=" * 90)
print("PART 2: Exact ratio g_k(m)/m^k as a fraction (looking for patterns)")
print("=" * 90)

for k in range(3, 10):
    print(f"\n--- k = {k} ---")
    for m in range(1, 11):
        g = g_exact(k, m)
        mk = Fraction(m**k)
        ratio = g / mk if mk != 0 else None
        # Also compute 1 - ((m-1)/m)^m for comparison
        if m >= 2:
            correction = 1 - Fraction(m - 1, m)**m
        else:
            correction = Fraction(1)
        print(f"  m={m:2d}: g_{k}/m^{k} = {str(ratio):>30s} = {float(ratio):12.10f}   "
              f"  s111d_factor = {str(correction):>20s} = {float(correction):12.10f}"
              f"  ratio/factor = {float(ratio)/float(correction):12.10f}" if correction != 0 else "")


# --------------------------------------------------------------------------
# PART 3: What is g_k(m)/m^k actually?
# Factor it as a polynomial in 1/m
# --------------------------------------------------------------------------
print("\n\n" + "=" * 90)
print("PART 3: g_k(m)/m^k expanded as polynomial in 1/m")
print("For k>=3: g_k(m) = (a*m^3+b*m^2+c*m+d)/3, so g_k/m^k = (a/3)/m^{k-3} + (b/3)/m^{k-2} + ...")
print("=" * 90)

for k in range(3, 10):
    a, b, c, d = gk_coeffs[k]
    # g_k(m) = (a*m^3 + b*m^2 + c*m + d) / 3
    # g_k(m)/m^k = a/(3*m^{k-3}) + b/(3*m^{k-2}) + c/(3*m^{k-1}) + d/(3*m^k)
    # For large m: ~ (a/3) * m^{3-k}
    print(f"\nk={k}: g_{k}(m)/m^{k} = ({a}/3)/m^{{{k-3}}} + ({b}/3)/m^{{{k-2}}} + ({c}/3)/m^{{{k-1}}} + ({d}/3)/m^{{{k}}}")
    print(f"  Leading term for large m: {Fraction(a, 3)} * m^{3-k}")
    print(f"  As m->inf: g_{k}/m^{k} -> {'0' if k > 3 else str(Fraction(a, 3))}")


# --------------------------------------------------------------------------
# PART 4: Compare g_k(m) / m^k with known combinatorial sequences
# --------------------------------------------------------------------------
print("\n\n" + "=" * 90)
print("PART 4: The 'true correction factor' r_k(m) = g_k(m)/m^k")
print("=" * 90)

# For k=3: g_3(m) = m(2m^2+1)/3
# g_3/m^3 = (2m^2+1)/(3m^2) = 2/3 + 1/(3m^2)
print("\nk=3: g_3(m)/m^3 = (2m^2+1)/(3m^2) = 2/3 + 1/(3m^2)")
for m in range(1, 11):
    exact = Fraction(2 * m**2 + 1, 3 * m**2)
    s111 = 1 - Fraction(m - 1, m)**m if m >= 2 else Fraction(1)
    print(f"  m={m:2d}: exact = {str(exact):>15s} = {float(exact):.10f}   "
          f"s111d = {float(s111):.10f}   diff = {float(exact - s111):.2e}")

# For k=4: g_4(m) = (10m^3-33m^2+50m-24)/3
# g_4/m^4 = (10m^3-33m^2+50m-24)/(3m^4)
print("\nk=4: g_4(m)/m^4 = (10m^3-33m^2+50m-24)/(3m^4)")
for m in range(1, 11):
    exact = Fraction(10 * m**3 - 33 * m**2 + 50 * m - 24, 3 * m**4)
    s111 = 1 - Fraction(m - 1, m)**m if m >= 2 else Fraction(1)
    print(f"  m={m:2d}: exact = {str(exact):>15s} = {float(exact):.10f}   "
          f"s111d = {float(s111):.10f}   diff = {float(exact - s111):.2e}")


# --------------------------------------------------------------------------
# PART 5: Verify THM-201 original and S111d against actual E_{2k}/E_0
# using CV^2 = sum 2*g_k(n-2k)/(n)_{2k}
# --------------------------------------------------------------------------
print("\n\n" + "=" * 90)
print("PART 5: E_{2k}/E_0 = 2*g_k(m)/(n)_{2k}  where m=n-2k")
print("Comparing: exact (THM-216) vs original (THM-201) vs S111d correction")
print("=" * 90)

for n in range(3, 20):
    max_k = (n - 1) // 2
    print(f"\nn = {n}:")
    cv2_exact = Fraction(0)
    cv2_orig = Fraction(0)
    cv2_s111 = Fraction(0)

    for k in range(1, max_k + 1):
        m = n - 2 * k
        if m < 1:
            break
        ff = falling_factorial(n, 2 * k)

        g = g_exact(k, m)
        e_exact = Fraction(2 * g, ff)

        o = thm201_original(k, m)
        e_orig = Fraction(2 * o, ff)

        s = thm201_s111d(k, m)
        e_s111 = Fraction(2 * s, ff)

        err_o = float(e_exact - e_orig)
        err_s = float(e_exact - e_s111)

        print(f"  k={k}: m={m:2d}  E_exact={float(e_exact):.10f}  "
              f"E_orig={float(e_orig):.10f} (err={err_o:+.2e})  "
              f"E_s111={float(e_s111):.10f} (err={err_s:+.2e})")

        cv2_exact += e_exact
        cv2_orig += e_orig
        cv2_s111 += e_s111

    print(f"  CV^2: exact={float(cv2_exact):.10f}  "
          f"orig={float(cv2_orig):.10f} (err={float(cv2_exact - cv2_orig):+.2e})  "
          f"s111={float(cv2_s111):.10f} (err={float(cv2_exact - cv2_s111):+.2e})")


# --------------------------------------------------------------------------
# PART 6: Search for a better correction factor
# --------------------------------------------------------------------------
print("\n\n" + "=" * 90)
print("PART 6: Searching for a BETTER correction factor")
print("=" * 90)

# The exact correction is r_k(m) = g_k(m)/m^k.
# For k=3: (2m^2+1)/(3m^2) -- this is NOT 1-((m-1)/m)^m
# Let's see if r_k(m) can be expressed as some function of m alone (independent of k)
# or if it must depend on k.

print("\nQ: Does r_k(m) depend on k for fixed m? (If not, universal correction exists)")
print(f"  {'m':>4}  {'r_3':>12}  {'r_4':>12}  {'r_5':>12}  {'r_6':>12}  {'s111d':>12}")
for m in range(1, 11):
    vals = []
    for k in range(3, 7):
        g = g_exact(k, m)
        r = float(g) / (m**k) if m > 0 else 0
        vals.append(r)
    s111 = float(1 - Fraction(m - 1, m)**m) if m >= 2 else 1.0
    print(f"  {m:4d}  {vals[0]:12.8f}  {vals[1]:12.8f}  {vals[2]:12.8f}  {vals[3]:12.8f}  {s111:12.8f}")

print("\nConclusion: r_k(m) DEPENDS on k, so no universal correction factor exists.")

# --------------------------------------------------------------------------
# PART 7: Leading-order correction for k>=3
# --------------------------------------------------------------------------
print("\n\n" + "=" * 90)
print("PART 7: Leading coefficient a_k/3 (the m->inf limit of m^{k-3}*r_k(m))")
print("=" * 90)

print(f"  k   a_k/3 (exact)       a_k/3 (float)       a_k/3 / a_{'{k-1}'}/3")
prev = None
for k in range(3, 10):
    a = Fraction(gk_coeffs[k][0], 3)
    ratio = float(a) / float(prev) if prev and prev != 0 else None
    print(f"  {k}   {str(a):>20s}   {float(a):18.6f}   "
          f"{'---' if ratio is None else f'{ratio:18.6f}'}")
    prev = a


# --------------------------------------------------------------------------
# PART 8: Check if g_k(m) = m^k * sum_{j=0}^{k-3} c_j(k) / m^j  for some c_j(k)
# i.e., g_k(m) is a polynomial of degree k in m with leading coeff 1?
# --------------------------------------------------------------------------
print("\n\n" + "=" * 90)
print("PART 8: Is g_k(m) close to any standard polynomial sequence?")
print("=" * 90)

# For k>=3: g_k(m) is degree 3. So g_k(m)/m^k ~ (a_k/3)/m^{k-3} -> 0 for k>3.
# This means THM-201 (g_k ~ m^k) GROSSLY overestimates for large k, small m.

print("\nKey insight: For k>=3, g_k is ALWAYS degree 3 regardless of k.")
print("THM-201 assumes degree k. The discrepancy grows as k increases.")
print("S111d correction [1-((m-1)/m)^m] is a function of m only, not k.")
print("")

# Check S111d accuracy: how close is m^k*[1-((m-1)/m)^m] to g_k(m)?
print("Relative error of S111d correction for each (k, m):")
print(f"  {'k\\m':>5}", end="")
for m in range(2, 11):
    print(f"  {'m='+str(m):>10}", end="")
print()

for k in range(3, 10):
    print(f"  k={k:2d}", end="")
    for m in range(2, 11):
        g = g_exact(k, m)
        s = thm201_s111d(k, m)
        if s != 0:
            rel_err = float((g - s) / s)
            print(f"  {rel_err:10.4e}", end="")
        else:
            print(f"  {'inf':>10}", end="")
    print()


# --------------------------------------------------------------------------
# PART 9: Check specific values: does g_k(m) match the ACTUAL Fourier energy
# computed from W(n)? (Cross-validate that THM-216 g_k are correct)
# --------------------------------------------------------------------------
print("\n\n" + "=" * 90)
print("PART 9: Cross-validation — W(n) vs sum of exact g_k terms")
print("=" * 90)

def compute_W_dp(n):
    """Compute W(n) = sum_{sigma in S_n} 2^{adj1(sigma)} via DP."""
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

print(f"  {'n':>3}  {'CV2 exact':>14}  {'CV2 gk-sum':>14}  {'CV2 THM201':>14}  {'CV2 S111d':>14}  {'match?':>8}")
for n in range(3, 16):
    W = compute_W_dp(n)
    nf = factorial(n)
    cv2_exact = Fraction(W, nf) - 1

    cv2_gk = Fraction(0)
    cv2_orig = Fraction(0)
    cv2_s111 = Fraction(0)
    max_k = (n - 1) // 2

    for k in range(1, max_k + 1):
        m = n - 2 * k
        if m < 1:
            break
        ff = falling_factorial(n, 2 * k)

        if k in gk_coeffs:
            cv2_gk += Fraction(2 * g_exact(k, m), ff)

        cv2_orig += Fraction(2 * thm201_original(k, m), ff)
        cv2_s111 += Fraction(2 * thm201_s111d(k, m), ff)

    match = "YES" if cv2_exact == cv2_gk else "NO"
    print(f"  {n:3d}  {float(cv2_exact):14.10f}  {float(cv2_gk):14.10f}  "
          f"{float(cv2_orig):14.10f}  {float(cv2_s111):14.10f}  {match:>8}")


# --------------------------------------------------------------------------
# PART 10: Summary — quantify how BAD the S111d correction is
# --------------------------------------------------------------------------
print("\n\n" + "=" * 90)
print("PART 10: SUMMARY — S111d correction quality")
print("=" * 90)

print("\nFor each (k,m), the exact ratio g_k(m)/m^k vs the S111d factor [1-((m-1)/m)^m]:")
print("(Values closer to 1.0 mean the two agree)")
print()
print(f"  {'k\\m':>5}", end="")
for m in range(1, 11):
    print(f"  {'m='+str(m):>10}", end="")
print()

for k in range(1, 10):
    print(f"  k={k:2d}", end="")
    for m in range(1, 11):
        g = g_exact(k, m)
        mk = Fraction(m**k)
        exact_ratio = float(g / mk)
        if k <= 2 or m <= 1:
            s111_factor = 1.0
        else:
            s111_factor = float(1 - Fraction(m - 1, m)**m)
        if s111_factor > 0:
            agreement = exact_ratio / s111_factor
        else:
            agreement = float('inf')
        print(f"  {agreement:10.6f}", end="")
    print()

print("\nIf the ratio is exactly 1.000000, S111d is exact for that (k,m).")
print("Deviations show the S111d correction is WRONG in general.")
print()
print("KEY FINDING: The correction factor depends on BOTH k and m,")
print("not just m. No single-variable correction can fix THM-201 for all k.")
