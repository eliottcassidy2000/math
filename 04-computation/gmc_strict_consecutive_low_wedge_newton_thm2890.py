#!/usr/bin/env python3
"""Exact companion for THM-2890.

The strict low consecutive wedge in THM-2879 is written in high coordinates

    X = 1/b > 0,                 Y = -a/b < 0,

with integer depth n >= 1 and a,b > 0.  This script clears the two exact
cubic-divisibility invariants, eliminates b, saturates the sole boundary
factor, and proves positivity of the remaining degree-fifteen eliminant in
the Gregory--Newton basis in the *integer* variable n.

All truth-bearing checks use ``require`` so ordinary and optimized Python run
the same proof audit.
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


def main() -> None:
    source = load_thm2879()
    n, x, y = source.n, source.x, source.y
    a, b, m = sp.symbols("a b m")

    # Replay the high-chart real-branch gate.  Its quadratic Gram factor is
    # also exactly the content removed from the linear selector, and is
    # positive on the entire real line, not only at Y>0.
    high_numerators = tuple(
        sp.together(invariant).as_numer_denom()[0]
        for invariant in (source.cubic_one, source.cubic_two)
    )
    high_subresultants = sp.subresultants(
        high_numerators[0], high_numerators[1], x
    )
    require(
        [sp.degree(item, x) for item in high_subresultants]
        == [4, 3, 2, 1, 0],
        "high selector subresultant profile changed",
    )
    high_y_factors = [
        (factor, exponent)
        for factor, exponent in sp.factor_list(high_subresultants[-1])[1]
        if sp.degree(factor, y) > 0
    ]
    gram = (
        (n + 2)
        + 2 * (2 * n + 3) * y
        + 2 * (2 * n + 3) * y**2
    )
    gram_factor = next(
        factor for factor, _ in high_y_factors if sp.degree(factor, y) == 2
    )
    gram_ratio = sp.cancel(
        sp.Poly(gram_factor, y).LC() / sp.Poly(gram, y).LC()
    )
    require(
        sp.expand(gram_factor - gram_ratio * gram) == 0
        and sp.discriminant(gram, y) == -4 * (2 * n + 3),
        "the globally positive Gram factor changed",
    )
    high_linear = sp.Poly(high_subresultants[-2], x)
    selector_content = sp.gcd(
        high_linear.nth(1), high_linear.nth(0)
    )
    selector_content_ratio = sp.cancel(selector_content / gram)
    require(
        not selector_content_ratio.has(y)
        and all(
            bool(coefficient > 0)
            for coefficient in sp.Poly(
                selector_content_ratio, n
            ).coeffs()
        ),
        "the selector content is not a positive multiple of Gram",
    )

    # The high-to-low projective rechart from THM-2879.
    low_pairs = tuple(
        sp.cancel(
            sp.together(
                invariant.subs({x: 1 / b, y: -a / b})
            )
        ).as_numer_denom()
        for invariant in (source.cubic_one, source.cubic_two)
    )
    low = tuple(pair[0] for pair in low_pairs)
    low_denominators = tuple(pair[1] for pair in low_pairs)
    require(
        [
            (
                sp.degree(polynomial, b),
                sp.degree(polynomial, a),
                sp.degree(polynomial, n),
                len(sp.Poly(polynomial, n, a, b).terms()),
            )
            for polynomial in low
        ]
        == [(5, 3, 8, 153), (5, 4, 7, 136)],
        "low-chart invariant profile changed",
    )
    require(
        all(
            all(bool(coefficient > 0) for coefficient in sp.Poly(
                denominator, n, a, b
            ).coeffs())
            for denominator in low_denominators
        ),
        "a low-chart clearing denominator can lose positivity",
    )

    # The affine b-resultant is the actual consequence object.
    resultant = sp.resultant(low[0], low[1], b)
    scalar, factors = sp.factor_list(resultant)
    factor_profile = sorted(
        (
            sp.degree(factor, a),
            sp.degree(factor, n),
            exponent,
            len(sp.Poly(factor, n, a).terms()),
        )
        for factor, exponent in factors
    )
    require(
        scalar == 432
        and factor_profile
        == [
            (0, 1, 2, 2),
            (0, 1, 3, 2),
            (0, 1, 3, 2),
            (0, 1, 18, 2),
            (0, 5, 1, 6),
            (0, 10, 1, 11),
            (1, 0, 2, 1),
            (2, 5, 1, 18),
            (4, 4, 1, 25),
            (15, 22, 1, 368),
        ],
        "low-chart resultant factor profile changed",
    )
    n_only = [
        factor for factor, _ in factors if sp.degree(factor, a) == 0
    ]
    require(
        all(
            all(bool(coefficient > 0) for coefficient in sp.Poly(
                factor, n
            ).coeffs())
            for factor in n_only
        ),
        "an n-only resultant factor can vanish at positive depth",
    )
    a_factor = next(
        factor for factor, _ in factors if sp.degree(factor, a) == 1
    )
    require(
        sp.cancel(a_factor / a).is_number,
        "the low boundary factor is not a scalar multiple of a",
    )
    q2 = next(
        factor for factor, _ in factors if sp.degree(factor, a) == 2
    )
    q4 = next(
        factor for factor, _ in factors if sp.degree(factor, a) == 4
    )
    r15 = next(
        factor for factor, _ in factors if sp.degree(factor, a) == 15
    )
    if sp.Poly(r15.subs(n, 1), a).LC() < 0:
        r15 = -r15
    require(
        all(
            bool(coefficient > 0)
            for coefficient in sp.Poly(q4, n, a).coeffs()
        ),
        "the quartic Gram projection lost coefficientwise positivity",
    )

    # Q2 is exactly the b=0 boundary factor.  For n>=2 it is already
    # coefficientwise positive after the integer shift n=m+2.
    constants = tuple(sp.Poly(polynomial, b).nth(0) for polynomial in low)
    constant_gcd = sp.factor(sp.gcd(*constants))
    require(
        sp.rem(constant_gcd, q2, a) == 0,
        "Q2 stopped dividing the common b=0 boundary",
    )
    q2_shift = sp.Poly(sp.expand(q2.subs(n, m + 2)), m, a)
    require(
        len(q2_shift.terms()) == 18
        and all(bool(coefficient > 0) for coefficient in q2_shift.coeffs()),
        "Q2 is not positive on the n>=2 integer tail",
    )

    # At n=1 Q2 has one positive root.  Work in Q[a]/(Q2), divide the
    # forced common b factor, and prove that the two residual quartics have
    # nonzero b-resultant.  The nonzero norm makes that resultant a unit in
    # the quadratic quotient, so b=0 is the only common root on Q2.
    q2_at_one = sp.Poly(q2.subs(n, 1), a, domain=sp.QQ)
    if q2_at_one.LC() < 0:
        q2_at_one = -q2_at_one
    q2_at_one = q2_at_one.monic()
    require(
        sp.expand(
            q2_at_one.as_expr()
            - (109 * a**2 - 435 * a - 595) / 109
        )
        == 0
        and q2_at_one.is_irreducible,
        "the exceptional n=1 quadratic changed",
    )
    fixed_constant_gcd = sp.factor(
        sp.gcd(*(constant.subs(n, 1) for constant in constants))
    )
    fixed_boundary_ratio = sp.cancel(
        fixed_constant_gcd / (a * q2_at_one.as_expr())
    )
    require(
        fixed_boundary_ratio != 0
        and not fixed_boundary_ratio.has(a),
        "the n=1 common constant is not exactly scalar*a*Q2",
    )
    residual_quartics: list[sp.Expr] = []
    for polynomial in low:
        fixed = sp.Poly(polynomial.subs(n, 1), b)
        reduced = sp.Integer(0)
        quotient_reduction = sp.Integer(0)
        for exponent in range(fixed.degree() + 1):
            coefficient = sp.Poly(
                fixed.nth(exponent), a, domain=sp.QQ
            ).rem(q2_at_one)
            quotient_reduction += coefficient.as_expr() * b**exponent
        for exponent in range(1, fixed.degree() + 1):
            coefficient = sp.Poly(
                fixed.nth(exponent), a, domain=sp.QQ
            ).rem(q2_at_one)
            reduced += coefficient.as_expr() * b ** (exponent - 1)
        require(
            sp.expand(quotient_reduction - b * reduced) == 0,
            "the n=1 quotient reduction is not exactly b times a quartic",
        )
        residual_quartics.append(sp.expand(reduced))
    require(
        [
            (sp.degree(item, b), sp.degree(item, a))
            for item in residual_quartics
        ]
        == [(4, 1), (4, 1)],
        "the n=1 saturated quotient profile changed",
    )
    saturated_resultant = sp.resultant(
        residual_quartics[0], residual_quartics[1], b
    )
    saturated_remainder = sp.Poly(
        saturated_resultant, a, domain=sp.QQ
    ).rem(q2_at_one)
    saturated_norm = sp.factor(
        sp.resultant(
            q2_at_one.as_expr(),
            saturated_remainder.as_expr(),
            a,
        )
    )
    require(
        saturated_remainder.degree() == 1
        and saturated_norm != 0,
        "the n=1 Q2 boundary acquired a finite common b root",
    )

    # Expand the genuine degree-fifteen factor in the Gregory--Newton basis:
    # R(n,a)=sum D_k(a) binom(n-1,k).  This is the coordinate natural for
    # integer depth and is strictly stronger here than monomial positivity.
    r15_poly = sp.Poly(r15, n, a)
    require(
        r15_poly.degree(n) == 22
        and r15_poly.degree(a) == 15
        and len(r15_poly.terms()) == 368,
        "R15 profile changed",
    )
    raw_signs = (
        sum(int(bool(coefficient > 0)) for coefficient in r15_poly.coeffs()),
        sum(int(bool(coefficient < 0)) for coefficient in r15_poly.coeffs()),
    )
    shifted_r15 = sp.Poly(sp.expand(r15.subs(n, m + 1)), m, a)
    shifted_signs = (
        sum(int(bool(coefficient > 0)) for coefficient in shifted_r15.coeffs()),
        sum(int(bool(coefficient < 0)) for coefficient in shifted_r15.coeffs()),
    )
    require(
        raw_signs == (326, 42) and shifted_signs == (364, 4),
        "the monomial-basis hostile profile changed",
    )
    rows = [
        sp.Poly(r15.subs(n, depth), a, domain=sp.QQ)
        for depth in range(1, 24)
    ]
    newton: list[sp.Poly] = []
    while rows:
        newton.append(rows[0])
        rows = [
            right - left for left, right in zip(rows, rows[1:])
        ]
    require(len(newton) == 23, "Newton depth changed")
    for order, polynomial in enumerate(newton[:22]):
        require(
            polynomial.degree() == 15
            and polynomial.eval(0) > 0
            and polynomial.LC() > 0
            and polynomial.count_roots(0, sp.oo) == 0,
            f"Newton coefficient D_{order} lost strict positivity",
        )
    terminal_shape = (
        (a + 2) ** 6
        * (a + 3) ** 3
        * (a**2 - 12 * a - 36) ** 2
        * (7 * a**2 + 24 * a + 18)
    )
    terminal_ratio = sp.cancel(newton[22].as_expr() / terminal_shape)
    require(
        terminal_ratio == 1024 * sp.factorial(22),
        "terminal Newton coefficient lost its nonnegative square shape",
    )

    # Explicitly replay the Newton identity as a polynomial identity.
    basis = sp.Integer(1)
    reconstruction = sp.Integer(0)
    for order, polynomial in enumerate(newton):
        if order > 0:
            basis = sp.cancel(basis * (n - order) / order)
        reconstruction += polynomial.as_expr() * basis
    require(
        sp.Poly(sp.expand(reconstruction - r15), n, a).is_zero,
        "Gregory--Newton reconstruction failed",
    )
    r15_digest = canonical_digest(r15_poly)
    newton_records = []
    for order, polynomial in enumerate(newton):
        for monomial, coefficient in polynomial.terms():
            newton_records.append(
                f"{order}:{monomial[0]}:{coefficient}"
            )
    newton_digest = sha256(
        ("\n".join(newton_records) + "\n").encode()
    ).hexdigest()
    saturation_digest = canonical_digest(
        sp.Poly(saturated_remainder, a, domain=sp.QQ)
    )

    # Fixed-depth hostile controls replay the load-bearing R15 factor.
    fixed_positive_root_counts = {
        depth: sp.Poly(r15.subs(n, depth), a).count_roots(0, sp.oo)
        for depth in (1, 2, 8, 64)
    }
    require(
        fixed_positive_root_counts == {1: 0, 2: 0, 8: 0, 64: 0},
        "a fixed-depth low branch reappeared",
    )

    print("THM-2890 DISCRETE-NEWTON STRICT LOW-WEDGE CLOSURE")
    print("status=PROVED+VERIFIED-EXACT")
    print("high_selector_content=positive(n)*Gram; discriminant=-4*(2n+3)")
    print("low_invariant_profiles=(b5,a3,n8,153);(b5,a4,n7,136)")
    print(f"resultant_factor_profile={factor_profile}")
    print("positive_factors=n-only,a-boundary,Q4")
    print("Q2_n_ge_2=18/18 shifted coefficients positive")
    print(
        "Q2_n1=(109*a^2-435*a-595)/109;"
        "saturated_residual_degrees=(4,4);norm_nonzero"
    )
    print(f"Q2_n1_saturation_sha256={saturation_digest}")
    print("R15=degree_n:22,degree_a:15,terms:368")
    print("R15_monomial_hostile=raw:326+/42-;shift_n=m+1:364+/4-")
    print(f"R15_sha256={r15_digest}")
    print("Newton_rows=23;D0..D21=strictly_positive_on_a>0")
    print(
        "D22=1024*22!*(a+2)^6*(a+3)^3*"
        "(a^2-12a-36)^2*(7a^2+24a+18)>=0"
    )
    print(f"Newton_sha256={newton_digest}")
    print(f"fixed_positive_root_counts={fixed_positive_root_counts}")
    print("consequence=no cubic-divisible strict low wedge for any integer n>=1")
    print(
        "composition=with THM-2879, all three strict consecutive "
        "cone-cutting wedge charts miss a shared cubic-quartic line"
    )
    print(
        "scope=strict consecutive wedges only; boundaries and the general "
        "mixed four-slot branch remain open"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
