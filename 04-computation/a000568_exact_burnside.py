"""
a000568_exact_burnside.py — monad-researcher-2026-06-02-S560

FULLY EXACT integer Burnside computation of A000568 (number of
non-isomorphic tournaments on n labeled-up-to-iso vertices).

MOTIVATION: a000568_asymptotic_exact_s26.py uses a float-based
approximation + correction scheme that SILENTLY LOSES PRECISION at n>=14:
  - it reports a(14) = 28401423719122300 (correct: 28401423719122304)
  - it reports a(15) = 31021002160355168256 (correct: 31021002160355166848)
This script avoids floats entirely: it sums exact Fractions over the
partitions of n into ODD parts (the only cycle types with nonzero fixed
tournaments), so every value is provably exact.

FORMULA (Davis 1954 / Burnside):
  T(n) = (1/n!) * sum_{sigma in S_n} fix(sigma)
where fix(sigma) = 2^t(sigma) if ALL cycles of sigma have ODD length,
else 0, and t(sigma) = number of orbits of unordered pairs {i,j} under sigma.

For a cycle type with cycles of (odd) lengths c_1,...,c_k:
  t = sum_i (c_i - 1)/2            # orbits of pairs INSIDE each cycle
    + sum_{i<j} gcd(c_i, c_j)      # orbits of pairs BETWEEN two cycles

Number of permutations of a given cycle type (m_j cycles of length j):
  n! / prod_j ( j^{m_j} * m_j! )

So T(n) = sum_{odd partitions of n}  2^t / prod_j ( j^{m_j} * m_j! )
which is an exact integer; we accumulate with Fraction and assert .denominator == 1.
"""
from fractions import Fraction
from math import gcd, factorial

# Known A000568 values (OEIS) for cross-verification.
KNOWN = {
    1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
    9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168,
    13: 48542114686912, 14: 28401423719122304,
    15: 31021002160355166848,
    16: 63530415842308265100288,
    17: 244912778438520759443245824,
    18: 1783398846284777975419600287232,
    19: 24605641171260376770598003978281472,
    20: 645022068557873570931850526424042500096,
}


def odd_partitions(n):
    """Yield partitions of n into odd parts as dict {part: multiplicity},
    with parts considered in nonincreasing order to avoid duplicates."""
    def rec(remaining, max_part, current):
        if remaining == 0:
            yield dict(current)
            return
        # try each odd part p from min(max_part, remaining) down to 1
        p = min(max_part, remaining)
        if p % 2 == 0:
            p -= 1
        while p >= 1:
            # use k copies of p
            k = 1
            while k * p <= remaining:
                current[p] = current.get(p, 0) + k
                yield from rec(remaining - k * p, p - 2, current)
                if current[p] == k:
                    del current[p]
                else:
                    current[p] -= k
                k += 1
            p -= 2
    yield from rec(n, n, {})


def pair_orbit_exponent(parts):
    """parts: dict {length: multiplicity} of an all-odd cycle type.
    Returns t = #orbits of unordered pairs under the permutation."""
    # expand to a flat list of cycle lengths
    lengths = []
    for L, m in parts.items():
        lengths.extend([L] * m)
    t = 0
    # inside each cycle: (L-1)/2 orbits
    for L in lengths:
        t += (L - 1) // 2
    # between each pair of cycles: gcd
    for i in range(len(lengths)):
        for j in range(i + 1, len(lengths)):
            t += gcd(lengths[i], lengths[j])
    return t


def num_perms_factor(parts):
    """Returns the rational 1 / prod_j ( j^{m_j} * m_j! ).
    Multiplying by n! gives the number of permutations of this cycle type."""
    denom = 1
    for L, m in parts.items():
        denom *= (L ** m) * factorial(m)
    return Fraction(1, denom)


def a000568(n):
    total = Fraction(0)
    for parts in odd_partitions(n):
        # sanity: all parts odd, sum to n
        t = pair_orbit_exponent(parts)
        total += num_perms_factor(parts) * (2 ** t)
    assert total.denominator == 1, f"non-integer total at n={n}: {total}"
    return total.numerator


if __name__ == "__main__":
    print("=" * 70)
    print("A000568 — EXACT INTEGER BURNSIDE (no floats)")
    print("monad-researcher-2026-06-02-S560")
    print("=" * 70)
    print()
    print("Cross-verification against OEIS known values:")
    all_ok = True
    for n in range(1, 21):
        val = a000568(n)
        if n in KNOWN:
            ok = (val == KNOWN[n])
            all_ok &= ok
            mark = "OK " if ok else "MISMATCH!"
            print(f"  a({n:2d}) = {val}    [{mark}]")
        else:
            print(f"  a({n:2d}) = {val}")
    print()
    print(f"All known values verified: {all_ok}")
    print()
    print("-" * 70)
    print("CONTRAST with a000568_asymptotic_exact_s26.py (float-based) at n=14,15:")
    print(f"  n=14 float script reported: 28401423719122300")
    print(f"  n=14 EXACT (this script):   {a000568(14)}")
    print(f"  n=15 float script reported: 31021002160355168256")
    print(f"  n=15 EXACT (this script):   {a000568(15)}")
    print()
    print("-" * 70)
    print("EXTENSION — exact values n=21..50 (guaranteed exact, big-int only):")
    import time
    for n in range(21, 51):
        t0 = time.time()
        val = a000568(n)
        print(f"  a({n}) = {val}   ({time.time()-t0:.3f}s)")
