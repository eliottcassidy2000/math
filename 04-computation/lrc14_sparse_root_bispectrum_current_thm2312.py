#!/usr/bin/env python3
"""Exact checks for THM-2312 (sparse-root bispectrum positivity).

All algebraic Fourier checks are performed in Q[zeta_p] by reducing
coefficient vectors modulo 1+x+...+x^(p-1).  No floating-point arithmetic
or Python ``assert`` statements are used.
"""

from fractions import Fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def poly_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def poly_mul(a, b):
    p = len(a)
    out = [Fraction(0) for _ in range(p)]
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y != 0:
                out[(i + j) % p] += x * y
    return tuple(out)


def poly_conj(a):
    p = len(a)
    out = [Fraction(0) for _ in range(p)]
    for i, x in enumerate(a):
        out[(-i) % p] += x
    return tuple(out)


def cyclotomic_reduce(a):
    """Canonical degree <=p-2 representative for prime p."""
    tail = a[-1]
    return tuple(a[i] - tail for i in range(len(a) - 1))


def fourier_value(weights, k):
    p = len(weights)
    out = [Fraction(0) for _ in range(p)]
    for r, value in enumerate(weights):
        out[(-k * r) % p] += value
    return tuple(out)


def direct_bispectrum_aggregate(weights):
    p = len(weights)
    hats = [fourier_value(weights, k) for k in range(p)]
    total = tuple(Fraction(0) for _ in range(p))
    pair_count = 0
    for k in range(1, p):
        for ell in range(1, p):
            if (k + ell) % p == 0:
                continue
            term = poly_mul(
                poly_mul(hats[k], hats[ell]),
                poly_conj(hats[(k + ell) % p]),
            )
            total = poly_add(total, term)
            pair_count += 1
    return cyclotomic_reduce(total), pair_count


def moment_formula(weights):
    p = len(weights)
    s1 = sum(weights, Fraction(0))
    s2 = sum((x * x for x in weights), Fraction(0))
    s3 = sum((x * x * x for x in weights), Fraction(0))
    return p * p * s3 - 3 * p * s1 * s2 + 2 * s1 * s1 * s1


def sharp_two_sheet_floor(p, mass):
    return Fraction((p - 2) * (p - 4), 4) * mass**3


def check_direct_case(weights, label):
    reduced, pair_count = direct_bispectrum_aggregate(weights)
    expected = moment_formula(weights)
    require(pair_count == (len(weights) - 1) * (len(weights) - 2),
            f"{label}: wrong pair count")
    require(reduced[0] == expected and all(x == 0 for x in reduced[1:]),
            f"{label}: cyclotomic aggregate mismatch")
    return pair_count


def main():
    print("THM-2312 SPARSE-ROOT BISPECTRUM AUDIT")

    direct_cases = 0
    for p in (5, 7, 13):
        for gap in range(1, p):
            for a in range(1, 5):
                for b in range(1, 5):
                    weights = [Fraction(0) for _ in range(p)]
                    weights[0] = Fraction(a)
                    weights[gap] = Fraction(b)
                    check_direct_case(weights, f"p={p},gap={gap},a={a},b={b}")
                    value = moment_formula(weights)
                    floor = sharp_two_sheet_floor(p, Fraction(a + b))
                    require(value >= floor, "sharp two-sheet floor failed")
                    require((value == floor) == (a == b),
                            "sharp equality classification failed")
                    direct_cases += 1

        for position in range(p):
            weights = [Fraction(0) for _ in range(p)]
            weights[position] = Fraction(position + 1)
            check_direct_case(weights, f"p={p},one-sheet={position}")
            require(
                moment_formula(weights)
                > sharp_two_sheet_floor(p, sum(weights, Fraction(0))),
                "one-sheet case should be strict",
            )
            direct_cases += 1

    print(f"direct cyclotomic cases: {direct_cases}")
    print("two-sheet formula:")
    print("  B_p=(p-2)*((p-1)*(a^3+b^3)-3*a*b*(a+b))")
    print("  B_p/(a+b)^3=(p-2)*((p-1)-3*p*x), x=ab/(a+b)^2")
    print("  x<=1/4 gives B_p>=(p-2)(p-4)(a+b)^3/4")

    p = 13
    pair_count = (p - 1) * (p - 2)
    fibre_mass_constant = Fraction((p - 2) * (p - 4), 4) * p**3
    word_measure_gain = 49
    whole_word_constant = word_measure_gain * fibre_mass_constant
    one_pair_constant = whole_word_constant / pair_count
    require(pair_count == 132, "C13 pair count changed")
    require(fibre_mass_constant == Fraction(217503, 4),
            "C13 fibre-to-mass constant changed")
    require(whole_word_constant == Fraction(10657647, 4),
            "C13 word-measure constant changed")
    require(one_pair_constant == Fraction(322959, 16),
            "C13 pigeonhole constant changed")
    print("C13 constants:")
    print(f"  nonzero character pairs: {pair_count}")
    print(f"  fibre-to-mass coefficient floor: {fibre_mass_constant}")
    print(f"  word measure gain from mu(Q)<=1/7: {word_measure_gain}")
    print(f"  whole-word coefficient floor: {whole_word_constant}")
    print(f"  one-pair real-part floor: {one_pair_constant}")

    strict_rho = Fraction(39002430583, 160481782761300)
    repeated_rho = Fraction(13560199813, 160481782761300)
    strict_floor = one_pair_constant * strict_rho**3
    repeated_floor = one_pair_constant * repeated_rho**3
    require(strict_floor > 0 and repeated_floor > 0,
            "word-stratum floors must be positive")
    print("THM-2305 word-stratum corollary:")
    print(f"  strict pair floor = {strict_floor}")
    print(f"  repeated-first pair floor = {repeated_floor}")

    p3_equal = moment_formula(
        [Fraction(1), Fraction(1), Fraction(0)]
    )
    p4_equal_formula = (
        4 * 4 * 2 - 3 * 4 * 2 * 2 + 2 * 2**3
    )
    support7 = [Fraction(1) if i < 7 else Fraction(0) for i in range(13)]
    signed_pair = [Fraction(1), Fraction(-1)] + [Fraction(0)] * 11
    require(p3_equal == -2, "p=3 hostile control changed")
    require(p4_equal_formula == 0, "p=4 boundary control changed")
    require(moment_formula(support7) == -42,
            "seven-sheet hostile control changed")
    require(moment_formula(signed_pair) == 0,
            "signed two-sheet hostile control changed")
    check_direct_case(support7, "C13 seven-sheet hostile")
    check_direct_case(signed_pair, "C13 signed hostile")
    print("sharp hostile boundaries:")
    print("  p=3, two equal sheets: -2")
    print("  p=4, two equal sheets: 0")
    print("  p=13, seven equal sheets: -42")
    print("  p=13, signed weights (1,-1): 0")
    print("AUDIT PASS")


if __name__ == "__main__":
    main()
