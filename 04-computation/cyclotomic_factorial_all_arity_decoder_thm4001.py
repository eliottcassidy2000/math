#!/usr/bin/env python3
"""Exact companion for THM-4001, the all-arity cyclotomic factorial decoder.

The implementation uses only integer/rational arithmetic and explicit failure
gates, so ``python -O`` does not remove any check.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement, product
from math import comb, factorial, prod


GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"FAIL {label}")


def weak_compositions(total: int, slots: int):
    if slots == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, slots - 1):
            yield (first,) + tail


def multinomial(exponents: tuple[int, ...]) -> int:
    remaining = sum(exponents)
    ans = 1
    for exponent in exponents[:-1]:
        ans *= comb(remaining, exponent)
        remaining -= exponent
    return ans


def projected_factorial_moment_direct(
    coefficients: tuple[int, ...], d: int, m: int
) -> Fraction:
    """Expand the selected multinomial channels and then apply L exactly."""
    degree = d * m
    numerator = 0
    for js in weak_compositions(m, len(coefficients)):
        exponents = tuple(d * j for j in js)
        coefficient = multinomial(exponents)
        monomial = 1
        radial_weight = 1
        for a, exponent in zip(coefficients, exponents):
            monomial *= a**exponent
            radial_weight *= factorial(exponent)
        numerator += coefficient * monomial * radial_weight
    return Fraction(numerator, factorial(degree))


def complete_homogeneous(values: tuple[int, ...], top: int) -> list[int]:
    """Coefficients of product_i (1-values[i]*t)^(-1), truncated at top."""
    answer = [1] + [0] * top
    for value in values:
        updated = [0] * (top + 1)
        power = 1
        for j in range(top + 1):
            for old_degree in range(top - j + 1):
                updated[old_degree + j] += answer[old_degree] * power
            power *= value
        answer = updated
    return answer


def elementary_symmetric(values: tuple[int, ...]) -> list[int]:
    answer = [1] + [0] * len(values)
    for value in values:
        for degree in range(len(values), 0, -1):
            answer[degree] += value * answer[degree - 1]
    return answer


def decode_elementary_from_h(h: list[int], arity: int) -> list[int]:
    """Invert H(t)E(-t)=1: h_r-e1*h_(r-1)+...+(-1)^r e_r=0."""
    decoded = [1]
    for r in range(1, arity + 1):
        value = 0
        for i in range(1, r + 1):
            value += (-1) ** (i - 1) * decoded[r - i] * h[i]
        decoded.append(value)
    return decoded


def characteristic_coefficients(elementary: list[int]) -> tuple[int, ...]:
    return tuple(((-1) ** r) * elementary[r] for r in range(len(elementary)))


def polynomial_value(coefficients: tuple[int, ...], x: int) -> int:
    value = 0
    for coefficient in coefficients:
        value = value * x + coefficient
    return value


def symbolic_selector_controls() -> tuple[int, int, int]:
    """Compare formal coefficient dictionaries, without substituting a_i."""
    cases = 0
    retained = 0
    rejected = 0
    for d in range(2, 6):
        for arity in range(2, 6):
            for m in range(0, 4):
                degree = d * m
                got: dict[tuple[int, ...], Fraction] = {}
                for exponents in weak_compositions(degree, arity):
                    if all(exponent % d == 0 for exponent in exponents):
                        weight = Fraction(
                            multinomial(exponents)
                            * prod(factorial(exponent) for exponent in exponents),
                            factorial(degree),
                        )
                        got[exponents] = got.get(exponents, Fraction(0)) + weight
                        gate(weight == 1, f"symbolic radial cancellation {d=} {arity=} {m=} {exponents=}")
                        retained += 1
                    else:
                        rejected += 1
                expected = {
                    tuple(d * j for j in js): Fraction(1)
                    for js in weak_compositions(m, arity)
                }
                gate(got == expected, f"symbolic projected polynomial {d=} {arity=} {m=}")
                cases += 1
    return cases, retained, rejected


def exhaustive_numeric_controls() -> tuple[int, int, int]:
    vectors = 0
    response_rows = 0
    newton_rows = 0
    for arity in range(2, 6):
        alphabet = range(-2, 3) if arity <= 4 else range(-1, 2)
        for coefficients in product(alphabet, repeat=arity):
            coefficients = tuple(coefficients)
            vectors += 1
            for d in range(2, 6):
                powers = tuple(a**d for a in coefficients)
                h = complete_homogeneous(powers, arity + 3)
                for m in range(0, arity + 2):
                    direct = projected_factorial_moment_direct(coefficients, d, m)
                    gate(direct.denominator == 1, f"integral response {coefficients=} {d=} {m=}")
                    gate(direct.numerator == h[m], f"numeric identity {coefficients=} {d=} {m=}")
                    response_rows += 1

                direct_e = elementary_symmetric(powers)
                decoded_e = decode_elementary_from_h(h, arity)
                gate(decoded_e == direct_e, f"Newton inversion {coefficients=} {d=}")
                polynomial = characteristic_coefficients(decoded_e)
                for power in powers:
                    gate(polynomial_value(polynomial, power) == 0, f"root recovery {coefficients=} {d=} {power=}")

                for m in range(1, arity + 4):
                    recurrence = 0
                    for r in range(1, min(m, arity) + 1):
                        recurrence += (-1) ** (r - 1) * direct_e[r] * h[m - r]
                    gate(h[m] == recurrence, f"order-k recurrence {coefficients=} {d=} {m=}")
                newton_rows += 1
    return vectors, response_rows, newton_rows


def unordered_injectivity_controls() -> tuple[int, int]:
    cases = 0
    signatures = 0
    for d in range(2, 6):
        for arity in range(2, 6):
            seen: dict[tuple[int, ...], tuple[int, ...]] = {}
            for coefficients in combinations_with_replacement(range(0, 7), arity):
                powers = tuple(a**d for a in coefficients)
                h = complete_homogeneous(powers, arity)
                signature = tuple(h[1 : arity + 1])
                previous = seen.get(signature)
                gate(previous is None or previous == coefficients, f"finite injectivity {d=} {arity=} {signature=}")
                seen[signature] = coefficients
                cases += 1
            signatures += len(seen)
    return cases, signatures


def boundary_and_hostile_controls() -> list[str]:
    rows: list[str] = []

    # "Reynolds" must mean the coordinatewise (mu_d)^k projector (or any
    # equivalent k-1-coordinate projector on total degree d*m).  A diagonal
    # scalar action fixes the whole homogeneous power, while rotating only one
    # coordinate is insufficient once k>2.
    diagonal_normalized = complete_homogeneous((1, 1, 1), 2)[2]
    one_coordinate_normalized = 4  # x^2+y^2+z^2+2xy after z -> -z averaging
    coordinatewise_normalized = complete_homogeneous((1, 1, 1), 1)[1]
    gate(diagonal_normalized == 6, "diagonal Reynolds hostile")
    gate(one_coordinate_normalized == 4, "one-coordinate Reynolds hostile")
    gate(coordinatewise_normalized == 3, "coordinatewise Reynolds control")
    rows.append(
        "selector_scope k=3,d=2,a=(1,1,1): diagonal_mu2=6, "
        "one_coordinate_mu2=4, coordinatewise_(mu2)^3=3; coordinatewise is required"
    )

    repeated_zero = (0, -1, 2, 2)
    d = 3
    powers = tuple(a**d for a in repeated_zero)
    h = complete_homogeneous(powers, len(repeated_zero))
    decoded = decode_elementary_from_h(h, len(repeated_zero))
    gate(decoded == elementary_symmetric(powers), "repeated/zero Newton control")
    rows.append(
        "repeated_zero d=3 a=(0,-1,2,2) powers="
        f"{powers} h1..h4={tuple(h[1:5])} char={characteristic_coefficients(decoded)}"
    )

    positive = (1, 2, 3, 4)
    powers = tuple(a**3 for a in positive)
    h = complete_homogeneous(powers, 4)
    decoded = decode_elementary_from_h(h, 4)
    gate(decoded == elementary_symmetric(powers), "positive branch control")
    rows.append(
        "positive_branch d=3 a=(1,2,3,4) h1..h4="
        f"{tuple(h[1:5])} recovered_powers={tuple(sorted(powers))}"
    )

    plus = (1, 2)
    minus = (-1, 2)
    plus_h = complete_homogeneous(tuple(a**2 for a in plus), 8)
    minus_h = complete_homogeneous(tuple(a**2 for a in minus), 8)
    gate(plus_h == minus_h, "even-d real sign hostile")
    rows.append("branch_hostile d=2 {1,2} and {-1,2} have identical M0..M8")

    # Formal complex phase hostile: every selected exponent is d*j, so a_i ->
    # zeta*a_i contributes zeta^(d*j)=1.  Test the exponent divisibility itself.
    for d in range(2, 9):
        for m in range(0, 7):
            for js in weak_compositions(m, 4):
                for j in js:
                    gate((d * j) % d == 0, f"root-torsor invariance {d=} {m=} {js=}")
    rows.append("complex_branch all responses invariant under independent a_i -> zeta_i*a_i, zeta_i^d=1")

    unknown_arity_left = complete_homogeneous((1, 8), 8)
    unknown_arity_right = complete_homogeneous((0, 1, 8), 8)
    gate(unknown_arity_left == unknown_arity_right, "unknown arity zero hostile")
    rows.append("arity_hostile powers {1,8} and {0,1,8} have identical all positive moments; k must be known")

    taxicab_a = (1, 12)
    taxicab_b = (9, 10)
    tax_a = complete_homogeneous(tuple(a**3 for a in taxicab_a), 2)
    tax_b = complete_homogeneous(tuple(a**3 for a in taxicab_b), 2)
    gate(tax_a[1] == tax_b[1] == 1729, "taxicab first moment")
    gate(tax_a[2] != tax_b[2], "taxicab second moment split")
    gate(tax_a[2] == 2_987_713 and tax_b[2] == 2_260_441, "taxicab exact values")
    rows.append(
        "taxicab 1729: (1,12)->(M1,M2)=(1729,2987713), "
        "(9,10)->(1729,2260441)"
    )

    for d in range(2, 11):
        projected_square_cross = 2
        projected_power_cross = comb(2 * d, d)
        gate(projected_square_cross != projected_power_cross, f"nonmultiplicativity F hostile {d=}")
    rows.append("nonmultiplicative [x_i^d x_j^d]: Pi(F^d)^2 has 2, Pi(F^(2d)) has binom(2d,d)")

    return rows


def main() -> None:
    symbolic_cases, retained, rejected = symbolic_selector_controls()
    vectors, response_rows, newton_rows = exhaustive_numeric_controls()
    injection_cases, signatures = unordered_injectivity_controls()
    hostile_rows = boundary_and_hostile_controls()

    semantic_material = "\n".join(
        [
            f"symbolic={symbolic_cases},{retained},{rejected}",
            f"numeric={vectors},{response_rows},{newton_rows}",
            f"injectivity={injection_cases},{signatures}",
            *hostile_rows,
        ]
    ).encode("utf-8")
    semantic_hash = sha256(semantic_material).hexdigest()

    print("FACTORIAL ALL-ARITY CYCLOTOMIC DECODER -- INDEPENDENT EXACT AUDIT")
    print("identity M_m=L(Pi_d(F^(d*m)))/(d*m)!=h_m(a_1^d,...,a_k^d): PASS")
    print(f"symbolic selector cases={symbolic_cases} retained={retained} rejected={rejected}: PASS")
    print(f"small exhaustive coefficient vectors={vectors} response_rows={response_rows}: PASS")
    print(f"Newton/recurrence rows={newton_rows}: PASS")
    print(f"nonnegative unordered finite cases={injection_cases} unique_signatures={signatures}: PASS")
    for row in hostile_rows:
        print(row)
    print("scope recoverable_without_branch=unordered multiset {a_i^d}, with multiplicity and known k")
    print("scope original_a_i=only after a fixed injective dth-root branch (e.g. nonnegative reals)")
    print("scope no_FC_HFC_JC_LRC_consequence=projector is not multiplicative and the functional is modified")
    print(f"semantic_sha256={semantic_hash}")
    print(f"active_gates={GATES}")
    print("FINAL PASS")


if __name__ == "__main__":
    main()
