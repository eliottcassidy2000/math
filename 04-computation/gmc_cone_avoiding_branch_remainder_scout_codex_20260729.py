#!/usr/bin/env python3
"""Exact branch-reduced scout for the consecutive cone-avoiding GMC chart.

This script reconstructs the cubic eliminant P(n,Y), the linear selector

    A(n,Y) X - N(n,Y) = 0,

and the endpoint-holonomy pseudo-remainder K(n,Y) from THM-2879.  It then
tests branch-aware Bernstein certificates on the two rational negative-Y
boxes suggested by exact fixed-depth roots:

    right: a=-Y in [1/3,4/11],
    left:  a=-Y in [4/11,1/2].

The result is a scoped scout, not an all-n endpoint-exit theorem.  All
truth-bearing checks use ``require``.
"""

from __future__ import annotations

import importlib.util
import os
from hashlib import sha256
from math import comb
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def canonical_digest(polynomial: sp.Poly) -> str:
    records = "\n".join(
        f"{','.join(str(exponent) for exponent in monomial)}:{coefficient}"
        for monomial, coefficient in polynomial.terms()
    )
    return sha256((records + "\n").encode()).hexdigest()


def load_thm2879():
    dependency = Path(__file__).with_name(
        "gmc_all_shift_cubic_null_endpoint_holonomy_thm2879.py"
    )
    dependency_bytes = dependency.read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(dependency_bytes).hexdigest()
        == "44012d84c88a22f246ef350f7f9a364116ac1fc839347361dee64c0a9c4a6e27",
        "THM-2879 exact dependency hash changed",
    )
    spec = importlib.util.spec_from_file_location("thm2879_exact", dependency)
    require(spec is not None and spec.loader is not None, "cannot load THM-2879")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bernstein_coefficients_on_interval(
    power_coefficients: list[sp.Rational | sp.Integer],
    lower: sp.Rational,
    upper: sp.Rational,
) -> list[sp.Expr]:
    """Return degree-d Bernstein coefficients on [lower,upper]."""

    degree = len(power_coefficients) - 1
    width = upper - lower
    # First convert p(a) to q(z)=p(lower+width*z).
    affine = [sp.Integer(0)] * (degree + 1)
    for exponent, coefficient in enumerate(power_coefficients):
        for z_exponent in range(exponent + 1):
            affine[z_exponent] += (
                coefficient
                * comb(exponent, z_exponent)
                * lower ** (exponent - z_exponent)
                * width**z_exponent
            )
    # If q(z)=sum_i c_i z^i=sum_k b_k B_(k,d)(z), then
    # b_k=sum_(i<=k)c_i*binom(k,i)/binom(d,i).
    return [
        sp.factor(
            sum(
                affine[index]
                * sp.Rational(comb(order, index), comb(degree, index))
                for index in range(order + 1)
            )
        )
        for order in range(degree + 1)
    ]


def interval_profile(
    polynomial: sp.Poly,
    lower: sp.Rational,
    upper: sp.Rational,
) -> tuple[int, int, int, int]:
    """Profile n-monomial coefficient polynomials after Y=-a.

    Returns ``(negative, positive, mixed, zero)`` according to whether all
    Bernstein coefficients on the a-interval have one strict sign.
    """

    require(len(polynomial.gens) == 2, "interval profile needs two variables")
    n_degree = polynomial.degree(0)
    y_degree = polynomial.degree(1)
    dictionary = polynomial.as_dict()
    counts = {"negative": 0, "positive": 0, "mixed": 0, "zero": 0}
    for n_exponent in range(n_degree + 1):
        # Y=-a contributes (-1)^j.
        coefficients = [
            dictionary.get((n_exponent, y_exponent), sp.Integer(0))
            * (-1) ** y_exponent
            for y_exponent in range(y_degree + 1)
        ]
        bernstein = bernstein_coefficients_on_interval(
            coefficients, lower, upper
        )
        if all(value < 0 for value in bernstein):
            counts["negative"] += 1
        elif all(value > 0 for value in bernstein):
            counts["positive"] += 1
        elif all(value == 0 for value in bernstein):
            counts["zero"] += 1
        else:
            counts["mixed"] += 1
    return (
        counts["negative"],
        counts["positive"],
        counts["mixed"],
        counts["zero"],
    )


def newton_coefficient_vectors(
    polynomial: sp.Poly,
    base: int = 1,
) -> list[list[sp.Expr]]:
    """Return the coefficient vectors of Delta_n^k P(base,Y).

    The outer index is k and the inner ascending index is the Y exponent.
    Constructing the difference tables coefficientwise avoids expanding a
    degree-900 bivariate expression at every step.
    """

    require(len(polynomial.gens) == 2, "Newton vectors need two variables")
    n, _ = polynomial.gens
    n_degree = polynomial.degree(0)
    y_degree = polynomial.degree(1)
    dictionary = polynomial.as_dict()
    vectors = [
        [sp.Integer(0)] * (y_degree + 1)
        for _ in range(n_degree + 1)
    ]
    for y_exponent in range(y_degree + 1):
        coefficient_poly = sp.Poly.from_dict(
            {
                (n_exponent,): coefficient
                for (n_exponent, y_exponent_key), coefficient
                in dictionary.items()
                if y_exponent_key == y_exponent
            },
            (n,),
            domain=sp.QQ,
        )
        differences = [
            coefficient_poly.eval(base + offset)
            for offset in range(n_degree + 1)
        ]
        for order in range(n_degree + 1):
            vectors[order][y_exponent] = differences[0]
            differences = [
                differences[index + 1] - differences[index]
                for index in range(len(differences) - 1)
            ]
    return vectors


def newton_bernstein_profile(
    vectors: list[list[sp.Expr]],
    lower: sp.Rational,
    upper: sp.Rational,
) -> tuple[int, int, int, int]:
    """Sign profile for Delta_n^k P(1,-a) on an a-interval."""

    counts = {"negative": 0, "positive": 0, "mixed": 0, "zero": 0}
    for vector in vectors:
        coefficients = [
            coefficient * (-1) ** exponent
            for exponent, coefficient in enumerate(vector)
        ]
        bernstein = bernstein_coefficients_on_interval(
            coefficients,
            lower,
            upper,
        )
        if all(value < 0 for value in bernstein):
            counts["negative"] += 1
        elif all(value > 0 for value in bernstein):
            counts["positive"] += 1
        elif all(value == 0 for value in bernstein):
            counts["zero"] += 1
        else:
            counts["mixed"] += 1
    return (
        counts["negative"],
        counts["positive"],
        counts["mixed"],
        counts["zero"],
    )


def newton_values(polynomial: sp.Expr, variable: sp.Symbol, base: int) -> list[sp.Expr]:
    """Finite differences Delta^k polynomial(base), through its degree."""

    current = sp.expand(polynomial)
    values: list[sp.Expr] = []
    for _ in range(sp.degree(current, variable) + 1):
        values.append(sp.expand(current.subs(variable, base)))
        current = sp.expand(
            current.subs(variable, variable + 1) - current
        )
    return values


def newton_values_poly(polynomial: sp.Poly, base: int) -> list[sp.Expr]:
    """Finite differences using exact evaluations rather than repeated expansion."""

    degree = polynomial.degree()
    differences = [
        polynomial.eval(base + offset)
        for offset in range(degree + 1)
    ]
    values: list[sp.Expr] = []
    for _ in range(degree + 1):
        values.append(differences[0])
        differences = [
            differences[index + 1] - differences[index]
            for index in range(len(differences) - 1)
        ]
    return values


def univariate_factor_record(
    factor,
    exponent: int,
    variable: sp.Symbol,
) -> tuple[str, str | None]:
    """Audit an exact FLINT factor as a univariate integer polynomial."""

    require(factor.degrees()[1] == 0, "resultant factor retained y")
    factor_poly = sp.Poly.from_dict(
        {
            (monomial[0],): int(coefficient)
            for monomial, coefficient in factor.to_dict().items()
        },
        (variable,),
        domain=sp.ZZ,
    )
    if factor_poly.LC() < 0:
        factor_poly = -factor_poly
    differences = newton_values_poly(factor_poly, 1)
    if all(value > 0 for value in differences):
        sign_status = "GN+"
        positive_root_count = 0
    elif all(value < 0 for value in differences):
        sign_status = "GN-"
        positive_root_count = 0
    else:
        sign_status = "mixed"
        positive_root_count = factor_poly.count_roots(1, sp.oo)
    shifted_coefficients = factor_poly.shift(1).all_coeffs()
    shifted_variations = sum(
        (left > 0) != (right > 0)
        for left, right in zip(
            [value for value in shifted_coefficients if value],
            [value for value in shifted_coefficients if value][1:],
        )
    )
    root_free_prime = next(
        (
            prime
            for prime in sp.primerange(5, 1000)
            if int(factor_poly.LC()) % prime
            and all(
                int(factor_poly.eval(residue)) % prime
                for residue in range(prime)
            )
        ),
        None,
    )
    record = (
        f"deg{factor_poly.degree()}^{exponent}:{sign_status}:"
        f"positive_real_roots={positive_root_count}:"
        f"shift_variations={shifted_variations}:"
        f"root_free_prime={root_free_prime}:"
        f"digest={canonical_digest(factor_poly)}"
    )
    linear = (
        f"({factor_poly.as_expr()})^{exponent}"
        if factor_poly.degree() == 1
        else None
    )
    return record, linear


def coefficient_mod_prime(coefficient: sp.Expr, prime: int) -> int:
    rational = sp.Rational(coefficient)
    numerator = int(rational.p) % prime
    denominator = int(rational.q) % prime
    require(denominator != 0, f"denominator vanished modulo {prime}")
    return numerator * pow(denominator, -1, prime) % prime


def specialize_mod_prime(
    dictionary: dict[tuple[int, int], sp.Expr],
    n_value: int,
    y_degree: int,
    prime: int,
    y: sp.Symbol,
) -> sp.Poly:
    coefficients = [0] * (y_degree + 1)
    powers = [pow(n_value, exponent, prime) for exponent in range(
        max(key[0] for key in dictionary) + 1
    )]
    for (n_exponent, y_exponent), coefficient in dictionary.items():
        coefficients[y_exponent] = (
            coefficients[y_exponent]
            + coefficient_mod_prime(coefficient, prime) * powers[n_exponent]
        ) % prime
    expression = sum(
        coefficient * y**exponent
        for exponent, coefficient in enumerate(coefficients)
    )
    return sp.Poly(expression, y, modulus=prime)


def modular_coprimality_census(
    first: sp.Poly,
    second: sp.Poly,
    primes: tuple[int, ...],
) -> list[tuple[int, int, int, int]]:
    """Return (prime, degree_drops, common_residues, max_gcd_degree)."""

    require(first.gens == second.gens and len(first.gens) == 2, "bad modular pair")
    _, y = first.gens
    first_dictionary = first.as_dict()
    second_dictionary = second.as_dict()
    first_y_degree = first.degree(y)
    second_y_degree = second.degree(y)
    rows: list[tuple[int, int, int, int]] = []
    for prime in primes:
        degree_drops = 0
        common_residues = 0
        max_gcd_degree = 0
        for residue in range(prime):
            first_fixed = specialize_mod_prime(
                first_dictionary, residue, first_y_degree, prime, y
            )
            second_fixed = specialize_mod_prime(
                second_dictionary, residue, second_y_degree, prime, y
            )
            if (
                first_fixed.degree() != first_y_degree
                or second_fixed.degree() != second_y_degree
            ):
                degree_drops += 1
                continue
            gcd_degree = sp.gcd(first_fixed, second_fixed).degree()
            max_gcd_degree = max(max_gcd_degree, gcd_degree)
            if gcd_degree > 0:
                common_residues += 1
        rows.append((prime, degree_drops, common_residues, max_gcd_degree))
    return rows


def determinant_mod_prime(matrix: list[list[int]], prime: int) -> int:
    """Return the determinant of a square integer matrix modulo ``prime``."""

    size = len(matrix)
    require(all(len(row) == size for row in matrix), "matrix is not square")
    work = [[entry % prime for entry in row] for row in matrix]
    determinant = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant = determinant * pivot_value % prime
        pivot_inverse = pow(pivot_value, -1, prime)
        for row in range(column + 1, size):
            if work[row][column]:
                multiplier = work[row][column] * pivot_inverse % prime
                for offset in range(column, size):
                    work[row][offset] = (
                        work[row][offset]
                        - multiplier * work[column][offset]
                    ) % prime
    return determinant % prime


def padded_resultant_mod_prime(
    first: sp.Poly,
    second: sp.Poly,
    first_degree: int,
    second_degree: int,
    prime: int,
) -> int:
    """Evaluate the fixed-degree Sylvester determinant modulo ``prime``.

    Keeping the generic degrees is essential: after specialization, a leading
    coefficient may vanish modulo p.  The padded determinant is still exactly
    the reduction of the generic resultant, whereas recomputing a resultant
    after silently lowering degrees is not.
    """

    require(first.gens == second.gens and len(first.gens) == 1, "bad univariate pair")
    variable = first.gens[0]
    first_coefficients = [
        int(first.nth(exponent)) % prime
        for exponent in range(first_degree, -1, -1)
    ]
    second_coefficients = [
        int(second.nth(exponent)) % prime
        for exponent in range(second_degree, -1, -1)
    ]
    size = first_degree + second_degree
    matrix: list[list[int]] = []
    for shift in range(second_degree):
        matrix.append(
            [0] * shift
            + first_coefficients
            + [0] * (size - shift - len(first_coefficients))
        )
    for shift in range(first_degree):
        matrix.append(
            [0] * shift
            + second_coefficients
            + [0] * (size - shift - len(second_coefficients))
        )
    return determinant_mod_prime(matrix, prime)


def modular_resultant_census(
    first: sp.Poly,
    second: sp.Poly,
    primes: tuple[int, ...],
) -> list[tuple[int, tuple[int, ...]]]:
    """Return each prime and its zero fixed-degree resultant residues."""

    require(first.gens == second.gens and len(first.gens) == 2, "bad modular pair")
    _, y = first.gens
    first_dictionary = first.as_dict()
    second_dictionary = second.as_dict()
    first_y_degree = first.degree(y)
    second_y_degree = second.degree(y)
    rows: list[tuple[int, tuple[int, ...]]] = []
    for prime in primes:
        zero_residues: list[int] = []
        for residue in range(prime):
            first_fixed = specialize_mod_prime(
                first_dictionary, residue, first_y_degree, prime, y
            )
            second_fixed = specialize_mod_prime(
                second_dictionary, residue, second_y_degree, prime, y
            )
            if padded_resultant_mod_prime(
                first_fixed,
                second_fixed,
                first_y_degree,
                second_y_degree,
                prime,
            ) == 0:
                zero_residues.append(residue)
        rows.append((prime, tuple(zero_residues)))
    return rows


def main() -> None:
    source = load_thm2879()
    n, x, y = source.n, source.x, source.y

    invariant_numerators = tuple(
        sp.together(invariant).as_numer_denom()[0]
        for invariant in (source.cubic_one, source.cubic_two)
    )
    subresultants = sp.subresultants(*invariant_numerators, x)
    require(
        [sp.degree(item, x) for item in subresultants] == [4, 3, 2, 1, 0],
        "cubic subresultant profile changed",
    )
    eliminant = next(
        factor
        for factor, _ in sp.factor_list(subresultants[-1])[1]
        if sp.degree(factor, y) == 15
    )
    if sp.Poly(eliminant, y).LC().subs(n, 1) < 0:
        eliminant = -eliminant
    eliminant_poly = sp.Poly(eliminant, n, y, domain=sp.QQ)

    linear = sp.Poly(subresultants[-2], x)
    common_content = sp.gcd(linear.nth(1), linear.nth(0))
    selector_a = sp.cancel(linear.nth(1) / common_content)
    selector_n = sp.cancel(-linear.nth(0) / common_content)
    if selector_a.subs({n: 1, y: 1}) < 0:
        selector_a = -selector_a
        selector_n = -selector_n
    selector_a_poly = sp.Poly(selector_a, n, y, domain=sp.QQ)
    selector_n_poly = sp.Poly(selector_n, n, y, domain=sp.QQ)

    endpoint_numerator, endpoint_denominator = sp.together(
        source.endpoint_holonomy
    ).as_numer_denom()
    require(
        all(coefficient > 0 for coefficient in sp.Poly(endpoint_denominator, n).coeffs()),
        "endpoint denominator lost positivity",
    )
    endpoint_polynomial = sp.Poly(endpoint_numerator, x)
    endpoint_degree = endpoint_polynomial.degree()
    require(endpoint_degree == 5, "endpoint x-degree changed")

    powers_a = [sp.Poly(1, n, y, domain=sp.QQ)]
    powers_n = [sp.Poly(1, n, y, domain=sp.QQ)]
    for _ in range(endpoint_degree):
        powers_a.append(powers_a[-1] * selector_a_poly)
        powers_n.append(powers_n[-1] * selector_n_poly)
    substituted_endpoint = sum(
        (
            sp.Poly(endpoint_polynomial.nth(index), n, y, domain=sp.QQ)
            * powers_n[index]
            * powers_a[endpoint_degree - index]
            for index in range(endpoint_degree + 1)
        ),
        sp.Poly(0, n, y, domain=sp.QQ),
    )
    require(
        (
            substituted_endpoint.degree(n),
            substituted_endpoint.degree(y),
            len(substituted_endpoint.terms()),
        )
        == (121, 55, 6820),
        "selector-cleared endpoint profile changed",
    )

    remainder = sp.Poly(
        sp.prem(substituted_endpoint.as_expr(), eliminant, y),
        n,
        y,
        domain=sp.QQ,
    )
    require(
        (
            remainder.degree(n),
            remainder.degree(y),
            len(remainder.terms()),
            canonical_digest(remainder),
        )
        == (
            982,
            14,
            14745,
            "c1ed1aa0ff682a7226f67c752aceb7bb"
            "4924708a2973126fe903c62c686d67a2",
        )
        and all(coefficient < 0 for coefficient in remainder.coeffs()),
        "branch-reduced endpoint remainder changed",
    )

    # The pseudo-remainder above is the canonical THM-2879 object, but it
    # carries large powers of LC_Y(P).  The primitive part in QQ[n][Y] is,
    # by Gauss's lemma, also the primitive numerator of the exact Euclidean
    # remainder in QQ(n)[Y].  Extract it directly from the sparse dictionary;
    # this avoids an expensive rational-function reconstruction.
    remainder_dictionary = remainder.as_dict()
    remainder_coefficients_n = [
        sp.Poly.from_dict(
            {
                (n_exponent,): coefficient
                for (n_exponent, y_exponent_key), coefficient
                in remainder_dictionary.items()
                if y_exponent_key == y_exponent
            },
            (n,),
            domain=sp.QQ,
        )
        for y_exponent in range(remainder.degree(y) + 1)
    ]
    field_remainder_content_poly = remainder_coefficients_n[0]
    for coefficient_poly in remainder_coefficients_n[1:]:
        field_remainder_content_poly = sp.gcd(
            field_remainder_content_poly,
            coefficient_poly,
        )
    field_remainder_n_content = field_remainder_content_poly.as_expr()
    primitive_dictionary: dict[tuple[int, int], sp.Expr] = {}
    for y_exponent, coefficient_poly in enumerate(remainder_coefficients_n):
        quotient = coefficient_poly.exquo(field_remainder_content_poly)
        for (n_exponent,), coefficient in quotient.as_dict().items():
            primitive_dictionary[(n_exponent, y_exponent)] = coefficient
    field_remainder_primitive = sp.Poly.from_dict(
        primitive_dictionary,
        (n, y),
        domain=sp.QQ,
    )
    field_remainder_content_factors = sp.factor_list(
        field_remainder_n_content
    )
    require(
        field_remainder_n_content != 0
        and all(
            all(value > 0 for value in newton_values(factor, n, 1))
            for factor, _ in field_remainder_content_factors[1]
        ),
        "a field-remainder content factor can vanish at positive integer depth",
    )

    # Exact branch separators.  At a=-Y, the cubic-eliminant signs are
    # negative, positive, negative at 1/3,4/11,1/2 for every integer n>=1.
    a = sp.symbols("a")
    separator_values: list[tuple[str, int, int]] = []
    for label, value, expected_sign in (
        ("one_third", sp.Rational(1, 3), -1),
        ("four_elevenths", sp.Rational(4, 11), 1),
        ("one_half", sp.Rational(1, 2), -1),
    ):
        specialized = sp.together(eliminant.subs(y, -value)).as_numer_denom()[0]
        differences = newton_values(specialized, n, 1)
        require(
            all(expected_sign * item > 0 for item in differences),
            f"separator {label} lost discrete-Newton sign",
        )
        separator_values.append((label, expected_sign, len(differences)))

    boxes = {
        "right": (sp.Rational(1, 3), sp.Rational(4, 11)),
        "left": (sp.Rational(4, 11), sp.Rational(1, 2)),
    }
    interval_profiles = {
        (name, object_name): interval_profile(polynomial, *interval)
        for name, interval in boxes.items()
        for object_name, polynomial in (
            ("P", eliminant_poly),
            ("A", selector_a_poly),
            ("N", selector_n_poly),
            ("K", field_remainder_primitive),
        )
    }
    newton_interval_profiles: dict[
        tuple[str, str],
        tuple[int, int, int, int],
    ] = {}
    if os.environ.get("GMC_SKIP_NEWTON") != "1":
        newton_vectors = {
            object_name: newton_coefficient_vectors(polynomial)
            for object_name, polynomial in (
                ("P", eliminant_poly),
                ("A", selector_a_poly),
                ("N", selector_n_poly),
                ("K", field_remainder_primitive),
            )
        }
        newton_interval_profiles = {
            (name, object_name): newton_bernstein_profile(
                newton_vectors[object_name],
                *interval,
            )
            for name, interval in boxes.items()
            for object_name in ("P", "A", "N", "K")
        }

    # Content/gcd controls: K is the actual branch-reduced consequence
    # object and shares no polynomial factor with P over Q[n,Y].
    remainder_rational_content = sp.Poly(remainder, n, y).content()
    require(
        sp.factorint(sp.Integer(remainder_rational_content)) == {2: 163, 3: 74},
        "endpoint remainder content changed",
    )
    primitive_remainder = field_remainder_primitive
    eliminant_remainder_gcd = sp.gcd(eliminant_poly, remainder)
    require(
        remainder_rational_content != 0
        and eliminant_remainder_gcd.total_degree() == 0,
        "endpoint remainder acquired a global cubic factor",
    )
    modular_primes = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
    modular_rows = modular_coprimality_census(
        eliminant_poly,
        primitive_remainder,
        modular_primes,
    )
    modular_resultant_rows = modular_resultant_census(
        eliminant_poly,
        primitive_remainder,
        modular_primes,
    )
    flint_resultant_profile = ""
    flint_prime_text = os.environ.get("GMC_FLINT_PRIME")
    if flint_prime_text:
        from flint import fmpz_mod_mpoly_ctx

        flint_prime = int(flint_prime_text)
        flint_context = fmpz_mod_mpoly_ctx.get(["n", "y"], flint_prime)
        flint_eliminant = flint_context.from_dict(
            {
                monomial: coefficient_mod_prime(coefficient, flint_prime)
                for monomial, coefficient in eliminant_poly.as_dict().items()
            }
        )
        flint_endpoint = flint_context.from_dict(
            {
                monomial: coefficient_mod_prime(coefficient, flint_prime)
                for monomial, coefficient
                in substituted_endpoint.as_dict().items()
            }
        )
        flint_resultant = flint_eliminant.resultant(flint_endpoint, 1)
        flint_unit, flint_factors = flint_resultant.factor()
        flint_resultant_profile = (
            f"p{flint_prime}:degree={flint_resultant.degrees()[0]}:"
            f"terms={len(flint_resultant)}:unit={flint_unit}:factors="
            + ",".join(
                f"deg{factor.degrees()[0]}^{exponent}"
                for factor, exponent in flint_factors
            )
        )
        print(f"flint_resultant_profile={flint_resultant_profile}")
        print(
            "flint_linear_factors="
            + ";".join(
                f"{factor}^{exponent}"
                for factor, exponent in flint_factors
                if factor.degrees()[0] == 1
            )
        )
    if os.environ.get("GMC_FLINT_EXACT") == "1":
        from flint import fmpz_mpoly_ctx

        integer_context = fmpz_mpoly_ctx.get(["n", "y"])
        _, integer_eliminant_poly = eliminant_poly.clear_denoms(convert=True)
        _, integer_endpoint_poly = substituted_endpoint.clear_denoms(convert=True)
        integer_eliminant = integer_context.from_dict(
            {
                monomial: int(coefficient)
                for monomial, coefficient in integer_eliminant_poly.as_dict().items()
            }
        )
        integer_endpoint = integer_context.from_dict(
            {
                monomial: int(coefficient)
                for monomial, coefficient in integer_endpoint_poly.as_dict().items()
            }
        )
        integer_resultant = integer_eliminant.resultant(integer_endpoint, 1)
        integer_unit, integer_factors = integer_resultant.factor()
        exact_records = [
            univariate_factor_record(factor, exponent, n)
            for factor, exponent in integer_factors
        ]
        exact_factor_profiles = [record for record, _ in exact_records]
        exact_linear_factors = [
            linear for _, linear in exact_records if linear is not None
        ]
        print(
            "flint_exact_resultant="
            f"degree={integer_resultant.degrees()[0]}:"
            f"terms={len(integer_resultant)}:unit={integer_unit}"
        )
        print("flint_exact_factor_profiles=" + ";".join(exact_factor_profiles))
        print("flint_exact_linear_factors=" + "*".join(exact_linear_factors))
        for label, auxiliary_poly in (
            ("selector_a", selector_a_poly),
            (
                "selector_content",
                sp.Poly(common_content, n, y, domain=sp.QQ),
            ),
        ):
            _, integer_auxiliary_poly = auxiliary_poly.clear_denoms(convert=True)
            integer_auxiliary = integer_context.from_dict(
                {
                    monomial: int(coefficient)
                    for monomial, coefficient
                    in integer_auxiliary_poly.as_dict().items()
                }
            )
            auxiliary_resultant = integer_eliminant.resultant(
                integer_auxiliary,
                1,
            )
            auxiliary_unit, auxiliary_factors = auxiliary_resultant.factor()
            auxiliary_records = [
                univariate_factor_record(factor, exponent, n)
                for factor, exponent in auxiliary_factors
            ]
            print(
                f"flint_exact_{label}_resultant="
                f"degree={auxiliary_resultant.degrees()[0]}:"
                f"terms={len(auxiliary_resultant)}:"
                f"unit={auxiliary_unit}"
            )
            print(
                f"flint_exact_{label}_factor_profiles="
                + ";".join(record for record, _ in auxiliary_records)
            )

    print("GMC CONE-AVOIDING BRANCH-REDUCED ENDPOINT SCOUT")
    print("status=FINITE-EXACT SCOUT / NO ALL-N EXIT")
    print(
        "profiles="
        f"P:{eliminant_poly.degree(n)},{eliminant_poly.degree(y)},"
        f"{len(eliminant_poly.terms())};"
        f"A:{selector_a_poly.degree(n)},{selector_a_poly.degree(y)},"
        f"{len(selector_a_poly.terms())};"
        f"N:{selector_n_poly.degree(n)},{selector_n_poly.degree(y)},"
        f"{len(selector_n_poly.terms())};"
        f"F:{substituted_endpoint.degree(n)},{substituted_endpoint.degree(y)},"
        f"{len(substituted_endpoint.terms())};"
        f"K:{remainder.degree(n)},{remainder.degree(y)},"
        f"{len(remainder.terms())}"
    )
    print(f"P_digest={canonical_digest(eliminant_poly)}")
    print(f"A_digest={canonical_digest(selector_a_poly)}")
    print(f"N_digest={canonical_digest(selector_n_poly)}")
    print(f"K_digest={canonical_digest(remainder)}")
    print(f"K_rational_content={remainder_rational_content}")
    print("K_rational_content_factorization=2^163*3^74")
    eliminant_lead = sp.Poly(eliminant, y).LC()
    pseudo_exponent = substituted_endpoint.degree(y) - sp.degree(eliminant, y) + 1
    print(f"pseudo_multiplier_exponent={pseudo_exponent}")
    print(
        "pseudo_multiplier_degree="
        f"{pseudo_exponent * sp.degree(eliminant_lead, n)}"
    )
    print(
        "field_remainder_n_content_degree="
        f"{sp.degree(field_remainder_n_content, n)}"
    )
    print(
        "field_remainder_n_content_factors="
        + "*".join(
            f"deg{sp.degree(factor, n)}^{exponent}"
            for factor, exponent in field_remainder_content_factors[1]
        )
    )
    print(
        "field_remainder_primitive_profile="
        f"{primitive_remainder.degree(n)},{primitive_remainder.degree(y)},"
        f"{len(primitive_remainder.terms())}"
    )
    print(
        "field_remainder_primitive_digest="
        f"{canonical_digest(primitive_remainder)}"
    )
    print(f"gcd(P,K)_total_degree={eliminant_remainder_gcd.total_degree()}")
    print(
        "modular_coprimality="
        + ";".join(
            f"p{prime}:degree_drops={degree_drops}:"
            f"common_residues={common_residues}:max_gcd_degree={max_gcd_degree}"
            for prime, degree_drops, common_residues, max_gcd_degree in modular_rows
        )
    )
    print(
        "modular_uniform_primes="
        + ",".join(
            str(prime)
            for prime, degree_drops, common_residues, _ in modular_rows
            if degree_drops == 0 and common_residues == 0
        )
    )
    print(
        "modular_fixed_degree_resultants="
        + ";".join(
            f"p{prime}:zero_residues={len(zero_residues)}"
            for prime, zero_residues in modular_resultant_rows
        )
    )
    print(
        "modular_resultant_zero_sets="
        + ";".join(
            f"p{prime}:{','.join(str(value) for value in zero_residues)}"
            for prime, zero_residues in modular_resultant_rows
        )
    )
    print(
        "modular_resultant_uniform_primes="
        + ",".join(
            str(prime)
            for prime, zero_residues in modular_resultant_rows
            if not zero_residues
        )
    )
    print(
        "separators="
        + ",".join(
            f"{label}:sign={sign}:newton_terms={terms}"
            for label, sign, terms in separator_values
        )
    )
    print(
        "bernstein_profiles="
        + ";".join(
            f"{box}:{object_name}:negative={profile[0]},positive={profile[1]},"
            f"mixed={profile[2]},zero={profile[3]}"
            for (box, object_name), profile in interval_profiles.items()
        )
    )
    if newton_interval_profiles:
        print(
            "newton_bernstein_profiles="
            + ";".join(
                f"{box}:{object_name}:negative={profile[0]},positive={profile[1]},"
                f"mixed={profile[2]},zero={profile[3]}"
                for (box, object_name), profile
                in newton_interval_profiles.items()
            )
        )
    print(
        "interpretation=the two rational boxes retain branch order; "
        "a boxwise sign proof exists only when every n-monomial coefficient "
        "has one strict Bernstein sign"
    )
    print(
        "scope=profiles and sufficient Bernstein tests only; no parametric "
        "root uniqueness and no all-n endpoint nonvanishing claimed"
    )


if __name__ == "__main__":
    main()
