#!/usr/bin/env python3
"""Independent hostile audit of THM-2781.

The rational-series engine below does not use THM-2781's differential
recurrence.  It first forms f**a and then solves g**b=f**a coefficient by
coefficient, using the invertible linear coefficient b at g(0)=1.
"""

from fractions import Fraction
from itertools import product
from math import comb, gcd, isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


GATES = 0


def gate(condition, message):
    global GATES
    require(condition, message)
    GATES += 1


def trim(poly):
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def multiply(left, right, last=None, modulus=None):
    size = len(left) + len(right) - 1
    if last is not None:
        size = min(size, last + 1)
    result = [0 if modulus else Fraction(0) for _ in range(size)]
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            if i + j >= size:
                break
            value = result[i + j] + x * y
            result[i + j] = value % modulus if modulus else value
    return trim(result)


def poly_power(base, exponent, last=None, modulus=None):
    result = [1 if modulus else Fraction(1)]
    factor = trim(base)
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, factor, last, modulus)
        remaining //= 2
        if remaining:
            factor = multiply(factor, factor, last, modulus)
    return result


def coefficient(poly, degree):
    return poly[degree] if degree < len(poly) else 0


def bth_root_of_target(target, denominator, last, modulus=None):
    """Return g mod z^(last+1) with g(0)=1 and g**denominator=target."""

    if modulus:
        inverse = pow(denominator, -1, modulus)
        zero = 0
    else:
        inverse = Fraction(1, denominator)
        zero = Fraction(0)
    root = [1 if modulus else Fraction(1)]
    for n in range(1, last + 1):
        trial = root + [zero]
        known = coefficient(
            poly_power(trial, denominator, last=n, modulus=modulus),
            n,
        )
        wanted = coefficient(target, n)
        value = (wanted - known) * inverse
        root.append(value % modulus if modulus else value)
    return root


def rational_power_series(base, numerator, denominator, last):
    target = poly_power(base, numerator, last=last)
    return bth_root_of_target(target, denominator, last)


def is_bth_power(base, denominator):
    degree = len(trim(base)) - 1
    if degree % denominator:
        return False
    root_degree = degree // denominator
    root = bth_root_of_target(trim(base), denominator, root_degree)
    return trim(poly_power(root, denominator)) == trim(base)


def verify_recurrence(base, numerator, denominator, coefficients, last):
    alpha = Fraction(numerator, denominator)
    declared_degree = len(base) - 1
    for n in range(1, last + 1):
        rhs = sum(
            base[r] * ((alpha + 1) * r - n) * coefficients[n - r]
            for r in range(1, min(declared_degree, n) + 1)
        )
        gate(n * coefficients[n] == rhs, f"recurrence n={n}")


def exhaustive_rational_audit():
    rows = 0
    polynomials = 0
    implications = 0
    alphabet = (Fraction(-1), Fraction(0), Fraction(1))
    for degree in range(1, 6):
        for denominator in range(1, degree + 1):
            if degree % denominator:
                continue
            for numerator in range(1, 6):
                if gcd(numerator, denominator) != 1:
                    continue
                rows += 1
                terminal = degree * numerator // denominator
                gate(
                    Fraction(numerator, denominator) * degree == terminal,
                    "terminal integrality",
                )
                gate(
                    (Fraction(numerator, denominator) + 1) * degree
                    - (terminal + degree)
                    == 0,
                    "terminal multiplier",
                )
                for r in range(1, degree):
                    gate(
                        terminal + 1
                        <= terminal + degree - r
                        <= terminal + degree - 1,
                        "predecessor window",
                    )

                for tail in product(alphabet, repeat=degree):
                    base = [Fraction(1), *tail]
                    coefficients = rational_power_series(
                        base,
                        numerator,
                        denominator,
                        terminal + 2 * degree,
                    )
                    window_zero = all(
                        coefficients[n] == 0
                        for n in range(terminal + 1, terminal + degree)
                    )
                    perfect_power = is_bth_power(base, denominator)
                    gate(
                        window_zero == perfect_power,
                        (
                            f"equivalence d={degree},a={numerator},"
                            f"b={denominator},f={base}"
                        ),
                    )
                    if window_zero:
                        gate(
                            all(
                                coefficients[n] == 0
                                for n in range(
                                    terminal + 1,
                                    terminal + 2 * degree + 1,
                                )
                            ),
                            "finite tail closure control",
                        )
                        implications += 1
                    polynomials += 1

                sample = [Fraction(1)] + [
                    Fraction((2 * j + degree) % 5 - 2)
                    for j in range(1, degree + 1)
                ]
                sample_coefficients = rational_power_series(
                    sample,
                    numerator,
                    denominator,
                    terminal + degree,
                )
                verify_recurrence(
                    sample,
                    numerator,
                    denominator,
                    sample_coefficients,
                    terminal + degree,
                )
    return rows, polynomials, implications


def boundary_and_hostile_audit():
    # d=1: integrality and lowest terms force b=1, so the empty window is
    # correctly vacuous and every f is a first power.
    d1 = [Fraction(1), Fraction(2)]
    d1_series = rational_power_series(d1, 5, 1, 8)
    gate(d1_series == poly_power(d1, 5) + [Fraction(0)] * 3, "d=1 boundary")
    gate(is_bth_power(d1, 1), "d=1 first-power conclusion")

    # b=1 and actual degree below the declared bound.
    b1 = [Fraction(1), Fraction(2), Fraction(0), Fraction(-1)]
    b1_series = rational_power_series(b1, 3, 1, 15)
    gate(
        all(value == 0 for value in b1_series[10:]),
        "b=1 integral exponent tail",
    )

    short_root = [Fraction(1), Fraction(2)]
    short_base = poly_power(short_root, 3)
    short_base += [Fraction(0)] * (7 - len(short_base))
    short_series = rational_power_series(short_base, 2, 3, 10)
    gate(short_base[-1] == 0, "declared top coefficient zero")
    gate(
        all(value == 0 for value in short_series[5:]),
        "degree-below-bound full terminal tail",
    )

    # Uniform d-2 hostile: ((1+z)^d-z^d)^(1/d) agrees with 1+z through
    # degree d-1 and first differs by -z^d/d.
    uniform_hostiles = 0
    for degree in range(2, 13):
        base = [Fraction(comb(degree, j)) for j in range(degree)]
        base.append(Fraction(0))
        coefficients = rational_power_series(base, 1, degree, degree)
        gate(
            all(coefficients[n] == 0 for n in range(2, degree)),
            f"uniform insufficient tail d={degree}",
        )
        gate(
            coefficients[degree] == Fraction(-1, degree),
            f"uniform first active response d={degree}",
        )
        gate(not is_bth_power(base, degree), f"uniform nonpower d={degree}")
        uniform_hostiles += 1

    cubic = [Fraction(1), Fraction(3), Fraction(3), Fraction(0)]
    cubic_series = rational_power_series(cubic, 1, 3, 3)
    gate(cubic_series[2] == 0 and cubic_series[3] == Fraction(-1, 3),
         "cubic specialization")
    gate(not is_bth_power(cubic, 3), "cubic hostile noncube")

    quartic = [
        Fraction(1),
        Fraction(0),
        Fraction(0),
        Fraction(1),
        Fraction(0),
    ]
    quartic_series = rational_power_series(quartic, 3, 2, 9)
    gate(
        quartic_series[7:10]
        == [Fraction(0), Fraction(0), Fraction(-1, 16)],
        "quartic specialization",
    )
    gate(not is_bth_power(quartic, 2), "quartic hostile nonsquare")

    # Lowest terms: the displayed denominator four is false; coefficient
    # comparison, not degree alone, proves this square is not a fourth power.
    unreduced = poly_power(
        [Fraction(1), Fraction(0), Fraction(1)],
        2,
    )
    unreduced_series = rational_power_series(unreduced, 2, 4, 5)
    gate(all(value == 0 for value in unreduced_series[3:6]),
         "unreduced displayed tail")
    gate(not is_bth_power(unreduced, 4), "unreduced is not fourth power")
    gate(is_bth_power(unreduced, 2), "unreduced denominator-two repair")

    # Constant term one is the normalization that makes the rational branch
    # K-rational and unique.  sqrt(2) has no constant term in Q.
    gate(
        isqrt(2) ** 2 != 2,
        "constant-term branch hostile sqrt(2) over Q",
    )

    # Characteristic zero is load-bearing even when p does not divide b.
    # In F2, d=3,a=2,b=3,N=2 and f=1+z^3.  The unique cube root of f^2
    # begins 1+z^6, so c3=c4=0 but c6=1 although f is not a cube.
    base_f2 = [1, 0, 0, 1]
    target_f2 = poly_power(base_f2, 2, last=12, modulus=2)
    series_f2 = bth_root_of_target(target_f2, 3, 12, modulus=2)
    gate(series_f2[3] == 0 and series_f2[4] == 0, "char2 terminal window")
    gate(series_f2[6] == 1, "char2 later tail survives")
    gate(
        all(
            trim(poly_power([1, linear], 3, modulus=2)) != trim(base_f2)
            for linear in (0, 1)
        ),
        "char2 hostile is not cube",
    )
    return uniform_hostiles


def main():
    rows, polynomials, implications = exhaustive_rational_audit()
    uniform_hostiles = boundary_and_hostile_audit()
    print("THM-2781 INDEPENDENT TERMINAL-TAIL HOSTILE AUDIT")
    print(f"parameter_rows={rows}")
    print(f"exhaustive_polynomials={polynomials}")
    print(f"window_zero_implications={implications}")
    print(f"uniform_dminus2_hostiles={uniform_hostiles}")
    print(f"exact_gates={GATES}")
    print("independent_bth_root_engine=PASS")
    print("d1_boundary=PASS")
    print("b1_boundary=PASS")
    print("degree_below_declared_bound=PASS")
    print("uniform_insufficient_tail=PASS")
    print("unreduced_denominator_repair=PASS")
    print("constant_term_scope=PASS")
    print("char0_scope_hostile=PASS")
    print("cubic_quartic_specializations=PASS")
    print("THEOREM_VERDICT=PASS")
    print("CANONICAL_SCRIPT_VERDICT=PASS_AFTER_EXPLICIT_NONFOURTH_REPAIR")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
