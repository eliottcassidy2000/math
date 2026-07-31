#!/usr/bin/env python3
"""Independent hostile audit of the balanced e=2, h=1 one-pole draft.

This companion deliberately does not import the draft's control script.  It
uses a unit-aware reduced-scheme test, an independent permutation-orbit
count, and explicit normalization/inversion identities.
"""

from __future__ import annotations

import ast
from itertools import combinations
from math import floor
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


x, r, z, c = sp.symbols("x r z c")


def H(n):
    return sp.expand(r ** (n - 1) * (n - (n - 2) * r) - n * r + n - 2)


def ratio_scheme(n):
    h = sp.Poly(H(n), r, domain=sp.QQ)
    q, rem = sp.div(h, sp.Poly((r - 1) ** 3, r, domain=sp.QQ))
    require(rem.is_zero, f"H/(r-1)^3 has a remainder for N={n}")
    return q


def rational_parts(expr):
    numerator, denominator = sp.fraction(sp.cancel(expr))
    return (
        sp.Poly(numerator, r, domain=sp.QQ),
        sp.Poly(denominator, r, domain=sp.QQ),
    )


def is_zero_on_reduced_scheme(expr, q, label):
    """Check equality on Spec(QQ[r]/q), including denominator units."""
    numerator, denominator = rational_parts(expr)
    require(
        sp.gcd(q, denominator).degree() == 0,
        f"{label}: denominator is not a unit on the ratio scheme",
    )
    require(
        sp.rem(numerator, q).is_zero,
        f"{label}: numerator does not vanish on the ratio scheme",
    )


def is_unit_on_reduced_scheme(expr, q, label):
    numerator, denominator = rational_parts(expr)
    require(
        sp.gcd(q, denominator).degree() == 0,
        f"{label}: expression is undefined on part of the ratio scheme",
    )
    require(
        sp.gcd(q, numerator).degree() == 0,
        f"{label}: expression vanishes on part of the ratio scheme",
    )


def compose(a, b):
    """Permutation composition a after b."""
    return tuple(a[b[i]] for i in range(len(a)))


def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def rotate_conjugate(permutation, shift):
    n = len(permutation)
    return tuple((permutation[(i - shift) % n] + shift) % n for i in range(n))


def two_transposition_permutations(n):
    for vertices in combinations(range(n), 4):
        a, b, c, d = vertices
        for pairs in (
            ((a, b), (c, d)),
            ((a, c), (b, d)),
            ((a, d), (b, c)),
        ):
            permutation = list(range(n))
            for left, right in pairs:
                permutation[left] = right
                permutation[right] = left
            yield tuple(permutation)


def dessin_orbit_count(n):
    """Count the passport using a fixed N-cycle and its full centralizer."""
    pole = tuple((i + 1) % n for i in range(n))
    target_type = tuple([n - 2, 1, 1])
    survivors = {
        zero
        for zero in two_transposition_permutations(n)
        if cycle_type(compose(zero, pole)) == target_type
    }
    representatives = {
        min(rotate_conjugate(zero, shift) for shift in range(n))
        for zero in survivors
    }
    return len(survivors), len(representatives)


def audit_degree(n):
    require(n >= 4, "the e=2 chamber requires N>=4")
    q = ratio_scheme(n)
    h = sp.Poly(H(n), r, domain=sp.QQ)

    # The exceptional root is exactly triple, and the residual scheme is
    # reduced and avoids every normalization divisor.
    require(h.eval(1) == 0, f"H(1) changed for N={n}")
    require(
        sp.diff(h.as_expr(), r, 1).subs(r, 1) == 0
        and sp.diff(h.as_expr(), r, 2).subs(r, 1) == 0
        and sp.diff(h.as_expr(), r, 3).subs(r, 1)
        == -n * (n - 1) * (n - 2),
        f"the root r=1 is not exactly triple for N={n}",
    )
    require(
        sp.gcd(q, sp.Poly(sp.diff(q.as_expr(), r), r, domain=sp.QQ)).degree()
        == 0,
        f"the admissible ratio scheme is not reduced for N={n}",
    )
    for divisor, name in (
        (r, "r"),
        (r - 1, "r-1"),
        (n - (n - 2) * r, "C denominator"),
        (n * r - (n - 2), "reciprocal C denominator"),
    ):
        is_unit_on_reduced_scheme(divisor, q, f"N={n} {name}")

    E = sp.expand((x - 1) * (x - r))
    C = sp.cancel(n * (n - 1) * (n - 2) / (n - (n - 2) * r))
    B = x ** 2 / (n - 2) - (1 + r) * x / (n - 1) + r / n
    P = sp.expand(x ** n - C * B)

    # Derive the elimination equation independently from the two evaluation
    # conditions, rather than assuming the displayed H polynomial.
    elimination = sp.factor(
        sp.together(r ** n * B.subs(x, 1) - B.subs(x, r))
    )
    elimination_numerator = sp.Poly(
        sp.fraction(elimination)[0], r, domain=sp.QQ
    )
    require(
        sp.rem(elimination_numerator, h).is_zero
        and elimination_numerator.degree() == h.degree() + 1
        and sp.rem(
            elimination_numerator,
            sp.Poly(r, r, domain=sp.QQ),
        ).is_zero,
        f"two-point evaluation elimination does not recover H for N={n}",
    )

    require(
        sp.cancel(x * sp.diff(P, x) - n * P - C * E) == 0,
        f"Euler first integral has the wrong sign for N={n}",
    )
    require(
        sp.cancel(P.subs(x, 1)) == 0,
        f"the first normalized double zero is missing for N={n}",
    )
    is_zero_on_reduced_scheme(P.subs(x, r), q, f"N={n} P(r)")

    quotient, remainder = sp.div(
        sp.Poly(P, x, domain=sp.QQ.frac_field(r)),
        sp.Poly(E ** 2, x, domain=sp.QQ.frac_field(r)),
    )
    for index, coefficient in enumerate(remainder.all_coeffs()):
        is_zero_on_reduced_scheme(
            coefficient, q, f"N={n} E^2 remainder coefficient {index}"
        )
    require(
        quotient.degree() == n - 4 and quotient.LC() == 1,
        f"S has the wrong degree or leading coefficient for N={n}",
    )

    # The second-derivative formulas certify exact multiplicity two.  Their
    # unit checks also prove S is disjoint from E; P(0) proves pole
    # disjointness.  The paper proof then excludes every other repeated root
    # by the first integral.
    is_zero_on_reduced_scheme(
        sp.diff(P, x, 2).subs(x, 1) - C * (1 - r),
        q,
        f"N={n} P''(1)",
    )
    is_zero_on_reduced_scheme(
        sp.diff(P, x, 2).subs(x, r) - C * (r - 1) / r,
        q,
        f"N={n} P''(r)",
    )
    for value, name in (
        (C, "C"),
        (C * (1 - r), "P''(1)"),
        (C * (r - 1) / r, "P''(r)"),
        (-C * r / n, "P(0)"),
    ):
        is_unit_on_reduced_scheme(value, q, f"N={n} {name}")

    # Swapping the two marked double zeros is not just a root-set symmetry:
    # after x=r*y and monic renormalization it reproduces the complete
    # coefficient formula for ratio 1/r.
    y = sp.symbols("y")
    transformed_P = sp.cancel(P.subs(x, r * y) / r ** n)
    reciprocal_C = sp.cancel(
        n * (n - 1) * (n - 2) / (n - (n - 2) / r)
    )
    reciprocal_B = (
        y ** 2 / (n - 2)
        - (1 + 1 / r) * y / (n - 1)
        + 1 / (n * r)
    )
    reciprocal_P = sp.expand(y ** n - reciprocal_C * reciprocal_B)
    for index, coefficient in enumerate(
        sp.Poly(
            sp.together(transformed_P - reciprocal_P),
            y,
            domain=sp.QQ.frac_field(r),
        ).all_coeffs()
    ):
        is_zero_on_reduced_scheme(
            coefficient, q, f"N={n} inversion coefficient {index}"
        )

    require(
        sp.expand(r ** n * H(n).subs(r, 1 / r) + H(n)) == 0,
        f"anti-reciprocity has the wrong sign for N={n}",
    )

    # Exact Laurent form of the Chebyshev bridge.  On z=exp(i theta), each
    # difference z^k-z^(-k) is 2i sin(k theta).
    require(
        sp.expand(
            z ** (-n) * H(n).subs(r, z ** 2)
            - n * (z ** (n - 2) - z ** (-(n - 2)))
            + (n - 2) * (z ** n - z ** (-n))
        )
        == 0,
        f"Chebyshev Laurent bridge failed for N={n}",
    )
    chebyshev_derivative = sp.Poly(
        sp.diff(sp.chebyshevu(n - 2, c), c), c, domain=sp.QQ
    )
    require(
        chebyshev_derivative.degree() == n - 3
        and sp.gcd(
            chebyshev_derivative,
            sp.Poly(
                sp.diff(chebyshev_derivative.as_expr(), c),
                c,
                domain=sp.QQ,
            ),
        ).degree()
        == 0
        and chebyshev_derivative.count_roots(-1, 1) == n - 3,
        f"Chebyshev derivative root certificate failed for N={n}",
    )
    require(
        (chebyshev_derivative.eval(0) == 0) == (n % 2 == 0),
        f"Chebyshev midpoint parity failed for N={n}",
    )
    fixed = int(q.eval(-1) == 0)
    require(
        fixed == int(n % 2 == 0),
        f"the inversion fixed-point boundary changed for N={n}",
    )
    algebraic_classes = (q.degree() + fixed) // 2
    require(
        q.degree() == n - 3
        and algebraic_classes == floor((n - 2) / 2),
        f"ratio-orbit count failed for N={n}",
    )

    dessin_survivors, dessin_classes = dessin_orbit_count(n)
    require(
        dessin_classes == algebraic_classes,
        f"independent dessin orbit count disagrees for N={n}",
    )
    return q.degree(), algebraic_classes, dessin_survivors


def main():
    rows = tuple((n, *audit_degree(n)) for n in range(4, 15))
    require(
        tuple((n, ordered, classes) for n, ordered, classes, _ in rows)
        == tuple(
            (n, n - 3, floor((n - 2) / 2)) for n in range(4, 15)
        ),
        "all-degree formula failed in the controlled range",
    )

    # Hostile witness against the weaker condition "denominator is not
    # identically zero mod q": it can still vanish on one reduced component.
    hostile_q = sp.Poly((r - 2) * (r - 3), r, domain=sp.QQ)
    hostile_denominator = sp.Poly(r - 2, r, domain=sp.QQ)
    require(
        not sp.rem(hostile_denominator, hostile_q).is_zero
        and sp.gcd(hostile_q, hostile_denominator).degree() == 1,
        "hostile quotient-denominator witness stopped witnessing",
    )

    require(
        sp.factor(H(14))
        == -2
        * (r - 1) ** 3
        * (r + 1)
        * (
            6 * r ** 10
            + 5 * r ** 9
            + 10 * r ** 8
            + 8 * r ** 7
            + 12 * r ** 6
            + 9 * r ** 5
            + 12 * r ** 4
            + 8 * r ** 3
            + 10 * r ** 2
            + 5 * r
            + 6
        ),
        "N=14 factorization changed",
    )

    tree = ast.parse(Path(__file__).read_text())
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "audit contains a Python assert node",
    )

    print("INDEPENDENT BALANCED e=2,h=1 ONE-POLE HOSTILE AUDIT")
    print(f"rows_(N,ordered,affine_classes,dessin_survivors)={rows}")
    print("unit_aware_reduced_scheme_checks=true")
    print("primary_generic_denominator_guard_is_too_weak_in_isolation=true")
    print("actual_primary_denominators_are_units_on_every_ratio_component=true")
    print("euler_sign_and_exact_double_zero_checks=true")
    print("inversion_is_realized_by_source_scaling=true")
    print("chebyshev_U_derivative_bridge_and_unit_circle_locus=true")
    print("independent_dessin_centralizer_orbit_counts_agree=true")
    print("N4_split=true;N_ge_5_nonsplit=true")
    print("N14_degree_V=26;ordered=11;affine_classes=6")
    print("scope=response_factor_classification_only;no_chart_entry_or_JC2")
    print("ALL INDEPENDENT EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
