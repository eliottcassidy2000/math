"""
tournament_count_formulas.py
kind-pasteur-2026-03-07-S39

Compute and analyze T(n) = number of non-isomorphic tournaments on n vertices.
OEIS A000568: 1, 1, 1, 2, 4, 12, 56, 456, 6880, ...

========================================================================
MAIN FORMULA (Davis 1954, Burnside's lemma):

  T(n) = sum_{lambda} 2^{t(lambda)} / z(lambda)

where:
  - lambda runs over ALL partitions of n into ODD parts
  - If lambda has parts l_1^{a_1}, l_3^{a_3}, l_5^{a_5}, ... then:
    t(lambda) = (1/2) * [sum_{r,s odd} a_r * a_s * gcd(r,s) - sum_r a_r]
              = sum_r a_r*(r-1)/2 + sum_r C(a_r,2)*r + sum_{r<s} a_r*a_s*gcd(r,s)
    z(lambda) = prod_r r^{a_r} * a_r!

Key fact: a permutation sigma in S_n fixes a tournament iff ALL cycle lengths are odd.

========================================================================
ASYMPTOTIC:
  T(n) ~ 2^{C(n,2)} / n!   as n -> infinity
  (The identity permutation dominates; correction ratio -> 1 exponentially fast)

========================================================================
CONNECTIONS:
  - T(n) = #{even graphs on n vertices} (Glasby-Praeger-Royle, 2023)
  - Strongly connected: A051337, related by OGF B(x) with g.f. 2 - 1/B(x)
  - Self-converse tournaments: A002785 (SC counts: 1,1,2,2,8,12,88,...)
  - Labeled tournaments: 2^{C(n,2)} (A006125)
========================================================================
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
from math import gcd, factorial, comb
from fractions import Fraction


def partitions_into_odd_parts(n, max_part=None):
    """Generate all partitions of n into odd parts, as lists of (part, multiplicity)."""
    if max_part is None:
        max_part = n
    if max_part % 2 == 0:
        max_part -= 1
    if n == 0:
        yield []
        return
    if max_part <= 0:
        return
    for count in range(n // max_part, 0, -1):
        remainder = n - count * max_part
        for rest in partitions_into_odd_parts(remainder, max_part - 2):
            yield [(max_part, count)] + rest
    yield from partitions_into_odd_parts(n, max_part - 2)


def t_value(partition):
    """Compute t(lambda) = (1/2)[sum_{r,s} a_r a_s gcd(r,s) - sum_r a_r].
    Equivalently: sum_r a_r*(r-1)/2 + C(a_r,2)*r + sum_{r<s} a_r*a_s*gcd(r,s)."""
    d = 0
    for size, mult in partition:
        d += mult * (size - 1) // 2
        d += comb(mult, 2) * size
    for i in range(len(partition)):
        for j in range(i + 1, len(partition)):
            s1, m1 = partition[i]
            s2, m2 = partition[j]
            d += m1 * m2 * gcd(s1, s2)
    return d


def z_value(partition):
    """Compute z_lambda = prod(r^{a_r} * a_r!) for the cycle index."""
    z = 1
    for size, mult in partition:
        z *= (size ** mult) * factorial(mult)
    return z


def tournament_count(n):
    """T(n) via Davis/Burnside formula."""
    if n <= 1:
        return 1
    total = Fraction(0)
    for partition in partitions_into_odd_parts(n):
        t = t_value(partition)
        z = z_value(partition)
        total += Fraction(2**t, z)
    return int(total)


# ============================================================
# 1. Compute the sequence T(1) ... T(20)
# ============================================================
print("=" * 70)
print("T(n) = number of non-isomorphic tournaments on n vertices (A000568)")
print("=" * 70)
print()
print("FORMULA: T(n) = sum_{lambda |- n, all parts odd} 2^t(lambda) / z(lambda)")
print()

T = []
for n in range(1, 21):
    t = tournament_count(n)
    T.append(t)
    print(f"  T({n:2d}) = {t}")

# Verify against OEIS
# OEIS A000568 starts at n=0: 1, 1, 1, 2, 4, 12, 56, ...
# Our T starts at n=1, so compare from index 1 onwards
OEIS_from_1 = [1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 903753248,
               154108311168, 48542114686912, 28401423719122304,
               31021002160355166848, 63530415842308265100288,
               244912778438520759443245824, 1783398846284777975419600287232,
               24605641171260376770598003978281472,
               645022068557873570931850526424042500096]
assert T == OEIS_from_1, "Mismatch with OEIS!"
print("\n  Matches OEIS A000568: CONFIRMED")


# ============================================================
# 2. Asymptotic formula: T(n) ~ 2^{C(n,2)} / n!
# ============================================================
print("\n" + "=" * 70)
print("ASYMPTOTIC: T(n) ~ 2^C(n,2) / n!")
print("=" * 70)
print()
print("  The identity permutation (cycle type 1^n) contributes 2^C(n,2)/n!.")
print("  All other odd-cycle permutations contribute exponentially less.")
print()

for n in range(1, 16):
    asymp = Fraction(2**comb(n, 2), factorial(n))
    exact = T[n-1]
    ratio = Fraction(exact, 1) / asymp
    print(f"  n={n:2d}: T(n) = {exact:>15d},  2^C(n,2)/n! = {float(asymp):>18.2f},  "
          f"T/asymp = {float(ratio):.8f}")


# ============================================================
# 3. Second-order asymptotic
# ============================================================
print("\n" + "=" * 70)
print("REFINED ASYMPTOTIC: T(n) = (2^C(n,2)/n!) * [1 + eps(n)]")
print("=" * 70)
print()
print("  Main term: partition (1^n) gives 2^{C(n,2)} / n!")
print("  Next term: partition (3, 1^{n-3}) gives ratio n(n-1)(n-2)/(3*4^{n-2})")
print("  [Derivation: delta_t = t(3,1^{n-3}) - C(n,2) = -2(n-2), so")
print("   ratio = 2^{-2(n-2)} * n!/(3*(n-3)!) = n(n-1)(n-2)/(3*4^{n-2})]")
print()

for n in range(3, 16):
    main = Fraction(2**comb(n, 2), factorial(n))
    total_correction = Fraction(0)
    for partition in partitions_into_odd_parts(n):
        if len(partition) == 1 and partition[0] == (1, n):
            continue
        t = t_value(partition)
        z = z_value(partition)
        total_correction += Fraction(2**t, z)

    eps = total_correction / main
    eps_31 = Fraction(n * (n-1) * (n-2), 3 * 4**(n-2))

    print(f"  n={n:2d}: eps(n) = {float(eps):.8f},  "
          f"leading term = {float(eps_31):.8f}  "
          f"({float(eps_31/eps)*100:.1f}% of correction)")


# ============================================================
# 4. Burnside decomposition by cycle type
# ============================================================
print("\n" + "=" * 70)
print("BURNSIDE DECOMPOSITION BY CYCLE TYPE")
print("=" * 70)

for n in range(2, 10):
    print(f"\n  n={n}:")
    for partition in partitions_into_odd_parts(n):
        t = t_value(partition)
        z = z_value(partition)
        n_perms = factorial(n) // z
        contrib = Fraction(2**t, z)
        cycle_desc = " + ".join(f"{s}^{m}" for s, m in partition)
        print(f"    lambda=({cycle_desc}): {n_perms:>6d} perms, Fix=2^{t}, "
              f"contrib = {float(contrib):.6f}")
    print(f"    => T({n}) = {tournament_count(n)}")


# ============================================================
# 5. Number of odd-part partitions (= number of terms in Burnside sum)
# ============================================================
print("\n" + "=" * 70)
print("NUMBER OF ODD-PART PARTITIONS OF n (= # terms in Burnside sum)")
print("=" * 70)
print()
for n in range(1, 25):
    count = sum(1 for _ in partitions_into_odd_parts(n))
    print(f"  n={n:2d}: {count} partitions into odd parts")


# ============================================================
# 6. Growth analysis: ratios
# ============================================================
print("\n" + "=" * 70)
print("GROWTH: T(n+1)/T(n) ~ 2^n/(n+1)")
print("=" * 70)
print()
print("  Since T(n) ~ 2^{C(n,2)}/n!, the ratio T(n+1)/T(n) ~ 2^n/(n+1).")
print()

for i in range(len(T) - 1):
    n = i + 1
    ratio_actual = T[i+1] / T[i]
    ratio_predicted = 2**n / (n + 1)
    rel = ratio_actual / ratio_predicted if ratio_predicted > 0 else 0
    print(f"  T({n+1})/T({n}) = {ratio_actual:>14.4f},  "
          f"2^n/(n+1) = {ratio_predicted:>14.4f},  ratio = {rel:.6f}")


# ============================================================
# 7. Small cases by hand
# ============================================================
print("\n" + "=" * 70)
print("SMALL CASES WORKED OUT")
print("=" * 70)

print("""
  n=1: Only partition: (1). t=0, z=1. T(1) = 2^0/1 = 1.
       (One tournament: the single vertex.)

  n=2: Only partition: (1^2). t=0+1*1=1, z=2. T(2) = 2^1/2 = 1.
       (One tournament: one arc between two vertices.)

  n=3: Partitions: (3), (1^3).
       (3): t=(3-1)/2=1, z=3.   => 2^1/3 = 2/3.
       (1^3): t=0+3*1=3, z=6.   => 2^3/6 = 8/6 = 4/3.
       T(3) = 2/3 + 4/3 = 2.
       (The 3-cycle and the transitive tournament.)

  n=4: Partitions: (3+1), (1^4).
       (3+1): t=1+0+gcd(3,1)=1+1=2, z=3*1=3.  => 2^2/3 = 4/3.
       (1^4): t=0+C(4,2)=6, z=24.  => 2^6/24 = 64/24 = 8/3.
       T(4) = 4/3 + 8/3 = 4.

  n=5: Partitions: (5), (3+1^2), (1^5).
       (5): t=2, z=5.  => 4/5.
       (3+1^2): t=1+0+1+2*gcd(3,1)=1+1+2=4, z=3*2=6.  => 16/6 = 8/3.
       (1^5): t=0+C(5,2)=10, z=120.  => 1024/120 = 128/15.
       T(5) = 4/5 + 8/3 + 128/15 = 12/15 + 40/15 + 128/15 = 180/15 = 12.

  n=6: Partitions: (5+1), (3^2), (3+1^3), (1^6).
       (5+1): t=2+0+gcd(5,1)=3, z=5*1=5.  => 8/5.
       (3^2): t=2*1+C(2,2)*3=2+3=5, z=9*2=18.  => 32/18 = 16/9.
       (3+1^3): t=1+0+C(3,2)+3*gcd(3,1)=1+3+3=7, z=3*6=18.  => 128/18 = 64/9.
       (1^6): t=C(6,2)=15, z=720.  => 32768/720 = 4096/90.
       T(6) = 8/5 + 16/9 + 64/9 + 4096/90 = 144/90 + 160/90 + 640/90 + 4096/90 = 5040/90 = 56.
""")


# ============================================================
# 8. Connection to project: what fraction are regular?
# ============================================================
print("=" * 70)
print("CONNECTION: fraction of regular tournaments among all non-iso tournaments")
print("=" * 70)
print()
print("  At odd n, regular tournaments have all out-degrees = (n-1)/2.")
print("  Non-iso regular tournament counts (A000571):")

# A000571: Regular tournaments
# 1, 1, 1, 3, 15, ...
# We know: n=3: 1 (3-cycle), n=5: 3 (Paley + 2 others), n=7: 15
# For now just note the connection
print("    n=3: 1 regular out of T(3)=2 (50%)")
print("    n=5: 3 regular out of T(5)=12 (25%)")
print("    n=7: 15 regular out of T(7)=456 (3.3%)")
print("    n=9: 122 regular out of T(9)=191536 (0.064%)")
print()
print("  Regular tournaments become vanishingly rare as n grows.")
print("  But they contain the H-maximizers (Paley).")


# ============================================================
# 9. Strongly connected decomposition
# ============================================================
print("\n" + "=" * 70)
print("STRONGLY CONNECTED: T(n) = sum over compositions of prod SC(k_i)")
print("=" * 70)
print()
print("  Every tournament decomposes uniquely into strongly connected components")
print("  in a transitive order. So: OGF B(x) = sum T(n)x^n satisfies")
print("  B(x) = 1/(1 - C(x)) where C(x) = sum SC(n)x^n = strongly connected OGF")
print("  => C(x) = 1 - 1/B(x), or equivalently SC(n) = IET(T(n)) via Moebius.")
print()

# Compute SC(n) from B(x) = sum T(n) x^n
# C(x) = 1 - 1/B(x)
# B = 1 + x + x^2 + 2x^3 + 4x^4 + 12x^5 + ...
# 1/B = 1 - x - x^2 + 2x^3 - ... (need to compute)

# Actually: 1/B(x) = sum b_n x^n where b_0=1, sum_{k=0}^n T(k)*b(n-k) = 0 for n>=1
T_ext = [1] + T  # T_ext[0]=1, T_ext[1]=T(1)=1, etc.
b = [Fraction(0)] * (len(T_ext))
b[0] = Fraction(1)
for n in range(1, len(T_ext)):
    s = Fraction(0)
    for k in range(n + 1):
        if k < len(T_ext):
            s += T_ext[k] * b[n - k] if n - k >= 0 and n - k < len(b) else Fraction(0)
    # T_ext[0]*b[n] + sum_{k=1}^n T_ext[k]*b[n-k] = 0
    # b[n] = -sum_{k=1}^n T_ext[k]*b[n-k]
    s2 = Fraction(0)
    for k in range(1, n + 1):
        if k < len(T_ext) and n - k >= 0:
            s2 += T_ext[k] * b[n - k]
    b[n] = -s2

# SC(n) = coefficient of x^n in 1 - 1/B(x) = -b[n] for n >= 1
print("  SC(n) = strongly connected non-iso tournaments:")
for n in range(1, 13):
    sc = -int(b[n])
    t = T_ext[n]
    frac = sc / t * 100 if t > 0 else 0
    print(f"    SC({n:2d}) = {sc:>15d}  ({frac:>6.2f}% of T({n})={t})")


# ============================================================
# 10. 2-adic valuation pattern
# ============================================================
print("\n" + "=" * 70)
print("2-ADIC VALUATIONS of T(n)")
print("=" * 70)
print()
for i, t in enumerate(T[:15]):
    n = i + 1
    v = 0
    x = t
    while x > 0 and x % 2 == 0:
        v += 1; x //= 2
    print(f"  v_2(T({n:2d})) = {v:2d}   (odd part = {x})")

print()
print("  The 2-adic valuation grows roughly as n/2, reflecting the")
print("  powers of 2 dividing n! in the denominator of the Burnside sum.")


print("\n\nDONE")
