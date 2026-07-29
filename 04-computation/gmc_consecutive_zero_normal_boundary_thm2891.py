#!/usr/bin/env python3
"""Exact companion for THM-2891.

For

    L(s^q)=q!,  f_j=s^j/j!,  d_j=f_(j+1)-f_j,

put e_1=d_n, e_2=d_(n+1), e_3=d_(n+2).  A plane through the
origin whose normal has a zero coordinate and which meets the positive
e-cone in a two-dimensional cone belongs to one of three projective
families:

    A_t = span(e_1, e_2+t e_3),
    B_t = span(e_1+t e_3, e_2),
    C_t = span(t e_1+e_2, e_3),                 t>0,

together with their three coordinate-face endpoints.

This script proves that the binary quadratic factorial moment never
divides the binary cubic factorial moment on any of these planes, for
every integer n>=0.  Families A and C have a coefficientwise-negative
division invariant.  The interlaced family B needs an exact resultant;
its sole sign-indefinite depth factor is negative at n=0 and becomes
coefficientwise positive after the integer shift n=m+1.

All truth-bearing checks use ``require`` so normal and optimized Python
execute the same audit.
"""

from __future__ import annotations

import importlib.util
from hashlib import sha256
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


def invariants(source, first, second) -> tuple[sp.Expr, sp.Expr]:
    g_0 = source.multilinear((first, first))
    g_1 = source.multilinear((first, second))
    g_2 = source.multilinear((second, second))
    c_0 = source.multilinear((first, first, first))
    c_1 = source.multilinear((first, first, second))
    c_2 = source.multilinear((first, second, second))
    c_3 = source.multilinear((second, second, second))
    return (
        sp.factor(3 * c_1 * g_0 * g_2 - c_3 * g_0**2 - 2 * c_0 * g_1 * g_2),
        sp.factor(3 * c_2 * g_0 * g_2 - 2 * c_3 * g_1 * g_0 - c_0 * g_2**2),
    )


def numerator_denominator(expression: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    numerator, denominator = sp.together(expression).as_numer_denom()
    return sp.factor(numerator), sp.factor(denominator)


def fixed_literal_resultant(source, depth: int) -> sp.Integer:
    """Independent literal-factorial check of the interlaced B boundary."""

    coordinate = sp.symbols("coordinate")
    first = {depth: sp.Integer(1), depth + 2: coordinate}
    second = {depth + 1: sp.Integer(1)}

    g_0 = source.direct_moment((first, first))
    g_1 = source.direct_moment((first, second))
    g_2 = source.direct_moment((second, second))
    c_0 = source.direct_moment((first, first, first))
    c_1 = source.direct_moment((first, first, second))
    c_2 = source.direct_moment((first, second, second))
    c_3 = source.direct_moment((second, second, second))
    first_invariant = sp.expand(
        3 * c_1 * g_0 * g_2 - c_3 * g_0**2 - 2 * c_0 * g_1 * g_2
    )
    second_invariant = sp.expand(
        3 * c_2 * g_0 * g_2
        - 2 * c_3 * g_1 * g_0
        - c_0 * g_2**2
    )
    require(
        (sp.degree(first_invariant, coordinate), sp.degree(second_invariant, coordinate))
        == (4, 3),
        f"literal depth {depth}: boundary invariant degree changed",
    )
    resultant = sp.Integer(
        sp.resultant(first_invariant, second_invariant, coordinate)
    )
    require(resultant != 0, f"literal depth {depth}: boundary resultant vanished")
    return resultant


def main() -> None:
    source = load_thm2879()
    n = source.n
    t, m = sp.symbols("t m")

    family_forms = {
        "A": ({0: sp.Integer(1)}, {1: sp.Integer(1), 2: t}),
        "B": ({0: sp.Integer(1), 2: t}, {1: sp.Integer(1)}),
        "C": ({0: t, 1: sp.Integer(1)}, {2: sp.Integer(1)}),
        "face_13": ({0: sp.Integer(1)}, {2: sp.Integer(1)}),
    }
    family_invariants = {
        label: invariants(source, first, second)
        for label, (first, second) in family_forms.items()
    }

    # Every clearing denominator is positive for n>=0.  It is independent
    # of t; coefficientwise positivity in n is an exact sufficient check.
    denominators: list[sp.Expr] = []
    for pair in family_invariants.values():
        for invariant in pair:
            _, denominator = numerator_denominator(invariant)
            denominators.append(denominator)
            require(
                not denominator.has(t)
                and all(
                    bool(coefficient > 0)
                    for coefficient in sp.Poly(denominator, n).coeffs()
                ),
                "a boundary clearing denominator can vanish at n>=0",
            )

    # The two transport-ordered boundary families are even simpler than
    # the general THM-2830 argument: one literal cubic remainder already
    # has a strict coefficient sign throughout the closed t>=0 quadrant.
    signed_profiles: dict[str, tuple[int, int, str]] = {}
    for label in ("A", "C", "face_13"):
        numerator, _ = numerator_denominator(family_invariants[label][0])
        polynomial = sp.Poly(numerator, n, t)
        require(
            all(bool(coefficient < 0) for coefficient in polynomial.coeffs()),
            f"family {label}: first invariant lost strict negative sign",
        )
        signed_profiles[label] = (
            sp.degree(numerator, t),
            len(polynomial.terms()),
            canonical_digest(polynomial),
        )

    # The interlaced B family is the finite high-chart boundary y=0.
    # Eliminate t from the two exact division invariants.
    b_numerators = tuple(
        numerator_denominator(invariant)[0]
        for invariant in family_invariants["B"]
    )
    require(
        (
            sp.degree(b_numerators[0], t),
            sp.degree(b_numerators[1], t),
            len(sp.Poly(b_numerators[0], n, t).terms()),
            len(sp.Poly(b_numerators[1], n, t).terms()),
        )
        == (4, 3, 45, 32),
        "interlaced boundary invariant profile changed",
    )
    resultant = sp.factor(sp.resultant(*b_numerators, t))
    scalar, factors = sp.factor_list(resultant)
    factor_profile = sorted(
        (
            sp.degree(factor, n),
            exponent,
            len(sp.Poly(factor, n).terms()),
        )
        for factor, exponent in factors
    )
    require(
        scalar == 1728
        and factor_profile
        == [
            (1, 2, 2),
            (1, 2, 2),
            (1, 3, 2),
            (1, 4, 2),
            (1, 17, 2),
            (10, 1, 11),
            (11, 1, 12),
        ],
        "interlaced boundary resultant factor profile changed",
    )
    expected_linear_factors = (
        (n + 1, 2),
        (n + 2, 17),
        (2 * n + 1, 4),
        (3 * n + 1, 3),
        (3 * n + 2, 2),
    )
    actual_linear_factors = tuple(
        (factor, exponent)
        for factor, exponent in factors
        if sp.degree(factor, n) == 1
    )
    require(
        len(actual_linear_factors) == len(expected_linear_factors)
        and all(
            any(
                sp.expand(actual - expected) == 0
                and actual_exponent == expected_exponent
                for actual, actual_exponent in actual_linear_factors
            )
            for expected, expected_exponent in expected_linear_factors
        )
        and all(
            factor.subs(n, 0) > 0
            and all(
                coefficient > 0
                for coefficient in sp.Poly(factor, n).all_coeffs()
            )
            for factor, _ in actual_linear_factors
        ),
        "exact positive linear resultant factors changed",
    )
    degree_ten = next(factor for factor, _ in factors if sp.degree(factor, n) == 10)
    degree_eleven = next(
        factor for factor, _ in factors if sp.degree(factor, n) == 11
    )
    if sp.Poly(degree_eleven, n).LC() < 0:
        degree_eleven = -degree_eleven
    require(
        all(
            bool(coefficient > 0)
            for coefficient in sp.Poly(degree_ten, n).coeffs()
        ),
        "degree-ten boundary factor lost positivity",
    )

    # The degree-eleven factor is genuinely not positive on the continuous
    # half-line: it has one root between 7/10 and 71/100.  Integer depth
    # sees only n=0 and the shifted tail n=m+1.
    p_at_zero = degree_eleven.subs(n, 0)
    shifted = sp.Poly(sp.expand(degree_eleven.subs(n, m + 1)), m)
    require(
        p_at_zero == -159876
        and len(shifted.terms()) == 12
        and all(bool(coefficient > 0) for coefficient in shifted.coeffs()),
        "integer-depth sign split for degree-eleven factor changed",
    )
    require(
        degree_eleven.subs(n, sp.Rational(7, 10)) < 0
        and degree_eleven.subs(n, sp.Rational(71, 100)) > 0
        and sp.count_roots(
            degree_eleven, sp.Rational(7, 10), sp.Rational(71, 100)
        )
        == 1,
        "continuous-depth hostile root changed",
    )

    # Reconstruct the exact consequence object and verify it is nonzero
    # at n=0 and positive on every integer n>=1.
    reconstructed = sp.Integer(scalar)
    for factor, exponent in factors:
        reconstructed *= factor**exponent
    require(sp.expand(reconstructed - resultant) == 0, "resultant reconstruction failed")
    require(
        resultant.subs(n, 0) != 0,
        "interlaced boundary resultant vanished at n=0",
    )

    literal_controls = {
        depth: fixed_literal_resultant(source, depth) for depth in (0, 1, 8)
    }

    print("THM-2891 ZERO-NORMAL CONSECUTIVE CONE-BOUNDARY AUDIT")
    print("status=PROVISIONAL EXACT CERTIFICATE / UNDER INDEPENDENT AUDIT")
    print("depth_domain=integer n>=0")
    print(
        "families=A_t=span(e1,e2+t*e3);"
        "B_t=span(e1+t*e3,e2);"
        "C_t=span(t*e1+e2,e3);t>0"
    )
    print("coordinate_faces=all three projective endpoints included")
    print(f"positive_denominator_checks={len(denominators)}")
    print(
        "signed_invariant_profiles="
        + ",".join(
            f"{label}:degree_t={profile[0]}:terms={profile[1]}"
            for label, profile in signed_profiles.items()
        )
    )
    print(
        "signed_invariant_digests="
        + ",".join(
            f"{label}:{profile[2]}" for label, profile in signed_profiles.items()
        )
    )
    print("interlaced_degrees=4,3")
    print("interlaced_terms=45,32")
    print(f"interlaced_resultant_profile={factor_profile}")
    print(
        "interlaced_linear_factors="
        "(n+1)^2,(n+2)^17,(2n+1)^4,(3n+1)^3,(3n+2)^2"
    )
    print(f"degree10_digest={canonical_digest(sp.Poly(degree_ten, n))}")
    print(f"degree11_digest={canonical_digest(sp.Poly(degree_eleven, n))}")
    print(f"degree11_at_n0={p_at_zero}")
    print(f"degree11_shift_terms={len(shifted.terms())}")
    print("continuous_hostile_interval=(7/10,71/100);root_count=1")
    print(
        "literal_controls="
        + ",".join(
            f"n{depth}:sign={sp.sign(value)}:digits={len(str(abs(value)))}"
            for depth, value in literal_controls.items()
        )
    )
    print(
        "consequence=no zero-normal-coordinate positive-cone boundary plane "
        "has Q dividing C"
    )
    print(
        "scope=no strict-wedge reproof; no cone-avoiding or general mixed "
        "four-slot claim; no new GMC2/DvdK/JC consequence"
    )


if __name__ == "__main__":
    main()
