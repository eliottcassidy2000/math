#!/usr/bin/env python3
"""
Euler number ratio analysis.

Investigates patterns in E_{n-2m}/E_n ratios for odd tangent numbers,
and connections to Bernoulli numbers and the deformed Euler number E_T(-1).
"""

from fractions import Fraction
from math import comb
from collections import Counter

# ─── Compute tangent numbers via boustrophedon/zigzag ────────────────────────

def compute_euler_zigzag(N):
    """Compute unsigned zigzag (tangent) numbers for odd indices 1,3,...,2N-1.
    Returns dict {odd_n: unsigned_tangent_number}."""
    max_idx = 2 * N
    a = [Fraction(0)] * (max_idx + 2)
    a[0] = Fraction(1)

    results = {}
    for n in range(1, max_idx + 1):
        if n % 2 == 1:
            for k in range(1, n + 1):
                a[k] = a[k] + a[k - 1]
            results[n] = int(a[n])
        else:
            for k in range(n - 1, -1, -1):
                a[k] = a[k] + a[k + 1]
    return results


def signed_euler(n, unsigned):
    """E_n (signed) = (-1)^{(n-1)/2} * |E_n|."""
    k = (n - 1) // 2
    return ((-1) ** k) * unsigned


def bernoulli_numbers(N):
    """Compute B_0, B_1, ..., B_N."""
    B = [Fraction(0)] * (N + 1)
    B[0] = Fraction(1)
    for m in range(1, N + 1):
        s = Fraction(0)
        for k in range(m):
            s += Fraction(comb(m + 1, k)) * B[k]
        B[m] = -s / Fraction(m + 1)
    return B


def factorize(n):
    """Simple trial division factorization."""
    if n <= 1:
        return Counter()
    fc = Counter()
    temp = n
    p = 2
    while p * p <= temp:
        while temp % p == 0:
            fc[p] += 1
            temp //= p
        p += 1
    if temp > 1:
        fc[temp] += 1
    return fc


def format_factors(fc):
    return " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(fc.items())) if fc else "1"


def main():
    print("=" * 70)
    print("EULER / TANGENT NUMBER RATIO ANALYSIS")
    print("=" * 70)

    # 1. Compute tangent numbers
    T = compute_euler_zigzag(8)  # n=1,3,...,15
    print("\n--- 1. Euler numbers for odd n ---")
    print(f"  n   |E_n|          E_n(signed)     factorization")
    for n in sorted(T):
        s = signed_euler(n, T[n])
        fc = factorize(T[n])
        print(f"  {n:>2}  {T[n]:>12}  {s:>14}     {format_factors(fc)}")

    # Verify against known values
    known = {1: 1, 3: 2, 5: 16, 7: 272, 9: 7936}
    for n, v in known.items():
        assert T[n] == v, f"Mismatch at n={n}: got {T[n]}, expected {v}"
    print("\n  [Verified against known values up to n=9]")

    # ─── 2. Ratios E_{n-2m}/E_n ───────────────────────────────────────────
    print("\n--- 2. Unsigned ratios |E_{n-2m}|/|E_n| ---")
    print(f"  n   m  n-2m   |E_{{n-2m}}|/|E_n|               decimal")
    for n in sorted(T):
        if n < 3:
            continue
        for m in range(1, n // 2 + 1):
            n2 = n - 2 * m
            if n2 >= 1 and n2 in T:
                ratio = Fraction(T[n2], T[n])
                print(f"  {n:>2}  {m:>2}  {n2:>4}   {str(ratio):>25}  {float(ratio):>15.8f}")

    # ─── 3. Check 2^{a} * |E_{n-2m}| / |E_n| for various a ────────────────
    for a_label, a_func in [("2S", lambda S: 2*S),
                             ("2S+1", lambda S: 2*S+1),
                             ("2S+2", lambda S: 2*S+2)]:
        print(f"\n--- 3. Pattern: 2^{{{a_label}}} * |E_{{n-2S}}| / |E_n| ---")
        print(f"  n   S  n-2S  2^a   ratio                        decimal")
        for n in sorted(T):
            if n < 3:
                continue
            for S in range(1, n // 2 + 1):
                n2 = n - 2 * S
                if n2 >= 1 and n2 in T:
                    a_val = a_func(S)
                    power = 2 ** a_val
                    val = Fraction(power * T[n2], T[n])
                    print(f"  {n:>2}  {S:>2}  {n2:>4}  {power:>4}   {str(val):>25}  {float(val):>15.8f}")

    # ─── 4. Bernoulli connection ────────────────────────────────────────────
    print("\n--- 4. Bernoulli number verification ---")
    print("  |E_{2k-1}| = 2^{2k}(2^{2k}-1)|B_{2k}|/(2k)")
    B = bernoulli_numbers(18)
    print(f"  k   2k  |B_{{2k}}|                formula_val  |E_{{2k-1}}|  match?")
    for k in range(1, 9):
        b2k = B[2 * k]
        formula = Fraction(4**k) * Fraction(4**k - 1) * abs(b2k) / Fraction(2 * k)
        n_odd = 2 * k - 1
        actual = T.get(n_odd, "?")
        match = "YES" if formula == actual else "NO"
        print(f"  {k:>2}  {2*k:>2}  {str(abs(b2k)):>20}  {str(formula):>12}  {str(actual):>12}  {match}")

    # ─── 4b. Bernoulli decomposition of ratios ─────────────────────────────
    print("\n--- 4b. Ratio decomposition via Bernoulli ---")
    print("  |E_{2j-1}|/|E_{2k-1}| = 2^{-2m} * (4^j-1)/(4^k-1) * (k/j) * |B_{2j}|/|B_{2k}|")
    print(f"  k  j  m  2^{{-2m}}      (4^j-1)/(4^k-1)   k/j      |B2j/B2k|            product    actual   ok?")
    for n in sorted(T):
        if n < 5:
            continue
        k = (n + 1) // 2
        for m in range(1, k):
            j = k - m
            if j < 1:
                continue
            n2 = 2 * j - 1
            if n2 not in T:
                continue

            f1 = Fraction(1, 4**m)
            f2 = Fraction(4**j - 1, 4**k - 1)
            f3 = Fraction(k, j)
            f4 = abs(B[2*j]) / abs(B[2*k])

            product = f1 * f2 * f3 * f4
            actual = Fraction(T[n2], T[n])
            ok = "YES" if product == actual else "NO"
            print(f"  {k:>2} {j:>2} {m:>2}  {str(f1):>10}  {str(f2):>18}  {str(f3):>6}  {str(f4):>20}  {str(product):>10} {str(actual):>10} {ok}")

    # ─── 5. The n=7 case in detail ─────────────────────────────────────────
    print("\n--- 5. Detailed n=7 analysis ---")
    print(f"  |E_7| = {T[7]} = 16 * 17 = 2^4 * 17")
    print(f"  |E_5| = {T[5]} = 2^4")
    print(f"  |E_3| = {T[3]} = 2")
    print(f"  |E_1| = {T[1]} = 1")
    print()
    print("  Known n=7 ratios (coeff / E_7):")
    cases = [("t3", 128, -272), ("t5", -64, -272), ("t7", 128, -272), ("bc", -128, -272)]
    for label, num, den in cases:
        r = Fraction(num, den)
        print(f"    {label}: {num}/{den} = {r}")
    print()

    print("  Key observation: denominator of all ratios divides 17.")
    print(f"    8/17, 4/17, 8/17, 8/17")
    print()

    # What are these in terms of tangent numbers?
    print("  Expressing coefficients as 2^a * E_{n-2S}:")
    # coeff for t3: involves 3-cycles, so S related to n-2S
    # In the formula, t3 at n=7 means we're looking at 3-cycles
    # 128 = 2^7, |E_5| = 16 = 2^4, so 128 = 2^3 * 16 = 2^3 * |E_5|
    # 64 = 2^6, |E_3| = 2, so 64 = 2^5 * 2 = 2^5 * |E_3|
    # 128 = 2^7, |E_1| = 1, so 128 = 2^7 * 1 = 2^7 * |E_1|
    decomps = [("t3", 128, 5, "S=1"), ("t5", 64, 3, "S=2"), ("t7", 128, 1, "S=3")]
    for label, coeff, e_idx, S_label in decomps:
        import math
        a = int(math.log2(coeff // T[e_idx]))
        print(f"    {label} ({S_label}): |coeff|={coeff} = 2^{a} * |E_{e_idx}| = 2^{a} * {T[e_idx]}")
        print(f"       -> 2^{a} * |E_{e_idx}| / |E_7| = 2^{a} * {T[e_idx]} / {T[7]} = {Fraction(coeff, T[7])}")

    print()

    # ─── 5b. Factor structure: 2^{2k}-1 ────────────────────────────────────
    print("--- 5b. Factor structure of 2^{2k}-1 ---")
    for k in range(1, 9):
        val = 4**k - 1
        fc = factorize(val)
        print(f"  4^{k}-1 = {val:>8} = {format_factors(fc)}")

    # ─── 5c. The odd part of |E_n| ──────────────────────────────────────────
    print("\n--- 5c. Odd part of |E_n| (removing all factors of 2) ---")
    for n in sorted(T):
        val = T[n]
        v2 = 0
        temp = val
        while temp > 0 and temp % 2 == 0:
            v2 += 1
            temp //= 2
        odd_part = temp
        k = (n + 1) // 2
        # Also compute (4^k - 1) / 3
        q = (4**k - 1) // 3
        print(f"  |E_{n:>2}| = 2^{v2} * {odd_part:>10}    (4^{k}-1)/3 = {q:>10}    odd_part / ((4^k-1)/3) = {Fraction(odd_part, q) if q > 0 else '?'}")

    # ─── 5d. The ratio |E_{n-2}|/|E_n| and (4^k-1) ────────────────────────
    print("\n--- 5d. Consecutive ratios and 4^k structure ---")
    print("  For |E_{n-2}|/|E_n| where n=2k-1, n-2=2(k-1)-1:")
    print("  Ratio = (1/4) * (4^{k-1}-1)/(4^k-1) * k/(k-1) * |B_{2k-2}|/|B_{2k}|")
    print()
    print(f"  n   k   ratio                  (4^k-1)  ratio*(4^k-1)         simplified")
    for n in sorted(T):
        if n < 3:
            continue
        n2 = n - 2
        if n2 >= 1 and n2 in T:
            k = (n + 1) // 2
            ratio = Fraction(T[n2], T[n])
            q = 4**k - 1
            prod = ratio * q
            print(f"  {n:>2}  {k:>2}   {str(ratio):>20}  {q:>8}  {str(prod):>20}   {float(prod):>15.8f}")

    # ─── 6. Summary of key patterns ────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print("""
Key findings:

1. TANGENT NUMBERS: |E_1|=1, |E_3|=2, |E_5|=16, |E_7|=272, |E_9|=7936,
   |E_11|=353792, |E_13|=22368256, |E_15|=1903757312

2. BERNOULLI FORMULA VERIFIED:
   |E_{2k-1}| = 2^{2k} * (2^{2k}-1) * |B_{2k}| / (2k)

3. RATIO DECOMPOSITION:
   |E_{2j-1}|/|E_{2k-1}| = (1/4^m) * (4^j-1)/(4^k-1) * (k/j) * |B_{2j}|/|B_{2k}|
   where m = k-j.  This factors the ratio into:
   - A pure power-of-2 part: 1/4^m
   - A "geometric" part: (4^j-1)/(4^k-1)
   - An index ratio: k/j
   - A Bernoulli ratio: |B_{2j}|/|B_{2k}|

4. THE n=7 DENOMINATORS:
   All ratios coeff/E_7 have denominator 17 because:
   - |E_7| = 2^4 * 17
   - The 17 comes from (4^4-1) = 255 = 3*5*17
   - Coefficients are powers of 2 times smaller tangent numbers
   - The odd prime 17 is the "new" prime entering at k=4

5. ODD PARTS OF TANGENT NUMBERS:
   The 2-adic structure and the (4^k-1) factor control which odd primes
   appear. The ratio 2^{a} * |E_{n-2S}| / |E_n| always has its odd part
   determined by (4^k-1) / (4^j-1).

6. FOR E_T(-1)/E_n SIMPLIFICATION:
   The coefficient 2^{p+2S} * E_{n-2S} / E_n in the deformed Euler number
   has the form:
     +/- 2^{p+2S} * (1/4^S) * (4^j-1)/(4^k-1) * (k/j) * |B_{2j}|/|B_{2k}|
   = +/- 2^p * (4^j-1)/(4^k-1) * (k/j) * |B_{2j}|/|B_{2k}|

   The power-of-2 part 2^{p+2S}/4^S = 2^p cancels!
   So the correction terms simplify to:
     2^{parts(I)} * (4^j-1)/(4^k-1) * (k/j) * |B_{2j}|/|B_{2k}| * I(T)
   where j = k - S_I.
""")


if __name__ == "__main__":
    main()
