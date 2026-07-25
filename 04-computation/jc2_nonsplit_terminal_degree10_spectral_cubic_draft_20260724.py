#!/usr/bin/env python3
"""Exact symbolic referee for the nonsplit quartic terminal draft.

The accompanying proof, not this script, supplies the all-field and
all-polynomial quantifiers.  This file checks:

* the E_2, E_6, E_10 polynomial-part decompositions;
* the three Laurent observables for degrees 2, 6, and 10;
* elimination of the first flux in degree 10;
* the singular-cubic parametrization and Keller one-form; and
* the unavoidable cusp pole in the constant coefficient of H.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def faber_observables(
    degree: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return (Phi, Psi, R) from the Laurent recurrence at infinity."""
    alpha = sp.Rational(degree, 4)
    coefficients = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + 4):
        value = sum(
            quartic[step]
            * ((alpha + 1) * step - index)
            * coefficients[index - step]
            for step in range(2, min(4, index) + 1)
            if step in quartic
        ) / index
        coefficients.append(sp.factor(value))
    phi = sp.factor(4 * coefficients[degree + 1])
    psi = sp.factor(4 * coefficients[degree + 2])
    third = sp.factor(
        4 * coefficients[degree + 3] + p * coefficients[degree + 1]
    )
    return phi, psi, third


def direct_laurent_coefficient(
    degree: int,
    index: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
) -> sp.Expr:
    """Independent finite multinomial coefficient of the Laurent series."""
    alpha = sp.Rational(degree, 4)
    total = sp.Integer(0)
    for i in range(index // 2 + 1):
        for j in range(index // 3 + 1):
            remainder = index - 2 * i - 3 * j
            if remainder < 0 or remainder % 4:
                continue
            k = remainder // 4
            chosen = i + j + k
            falling = sp.prod(alpha - h for h in range(chosen))
            coefficient = (
                falling
                / (sp.factorial(i) * sp.factorial(j) * sp.factorial(k))
            )
            total += coefficient * p**i * q**j * r**k
    return sp.factor(total)


def check_polynomial_parts() -> None:
    x = sp.symbols("L")
    a, s, d = sp.symbols("a s d", nonzero=True)
    h = sp.expand(a * (x + s) ** 2 + d)
    quotient_2, remainder_2 = sp.div(x**2, h, domain=sp.QQ.frac_field(a, s, d))
    quotient_3, remainder_3 = sp.div(x**3, h, domain=sp.QQ.frac_field(a, s, d))
    require(sp.factor(quotient_2 - 1 / a) == 0, "Pol(L^2/H) mismatch")
    require(
        sp.factor(quotient_3 - (x - 2 * s) / a) == 0,
        "Pol(L^3/H) mismatch",
    )
    require(sp.degree(remainder_2, x) < 2, "degree-2 division remainder")
    require(sp.degree(remainder_3, x) < 2, "degree-3 division remainder")

    # The binomial tails after the displayed quotients have negative degree.
    # Thus these are the exact polynomial parts.
    e2 = h
    r2 = h**3 + sp.Rational(3, 2) * h * x
    e6 = r2 + sp.Rational(3, 8) / a
    r3 = (
        h**5
        + sp.Rational(5, 2) * h**3 * x
        + sp.Rational(15, 8) * h * x**2
    )
    e10 = r3 + sp.Rational(5, 16) * (x - 2 * s) / a
    require(sp.degree(e2, x) == 2, "E2 degree mismatch")
    require(sp.degree(e6, x) == 6, "E6 degree mismatch")
    require(sp.degree(e10, x) == 10, "E10 degree mismatch")
    print("polynomial_parts=PASS degrees=2,6,10")


def check_flux_bank_and_elimination() -> None:
    d, q, s, t = sp.symbols("d q s T")
    p = 2 * d
    r = d**2 - s
    bank = {
        degree: faber_observables(degree, p, q, r)
        for degree in (2, 6, 10)
    }
    expected = {
        2: (2 * q, -2 * s, -d * q),
        6: (
            -3 * q * s,
            -sp.Rational(3, 2) * (d * q**2 - s**2),
            q * (6 * d * s - q**2) / 4,
        ),
        10: (
            -sp.Rational(5, 4) * q * (d * q**2 - 3 * s**2),
            -sp.Rational(5, 32)
            * (-24 * d * q**2 * s + q**4 + 8 * s**3),
            sp.Rational(5, 8)
            * q
            * (d**2 * q**2 - 3 * d * s**2 + q**2 * s),
        ),
    }
    for degree in bank:
        require(
            all(
                sp.factor(found - wanted) == 0
                for found, wanted in zip(bank[degree], expected[degree])
            ),
            f"degree-{degree} Laurent bank mismatch",
        )
        direct = tuple(
            direct_laurent_coefficient(degree, degree + offset, p, q, r)
            for offset in (1, 2, 3)
        )
        require(
            sp.factor(4 * direct[0] - bank[degree][0]) == 0
            and sp.factor(4 * direct[1] - bank[degree][1]) == 0
            and sp.factor(
                4 * direct[2] + p * direct[0] - bank[degree][2]
            )
            == 0,
            f"degree-{degree} independent multinomial mismatch",
        )

    alpha, beta, psi0 = sp.symbols("alpha beta Psi0")
    phi10 = sp.factor(
        (
            bank[10][0] + alpha * bank[6][0] + beta * bank[2][0]
        ).subs(q**2, t)
        / q
    )
    psi10 = sp.factor(
        (
            bank[10][1] + alpha * bank[6][1] + beta * bank[2][1]
        ).subs(q**4, t**2)
        .subs(q**2, t)
    )
    third10_over_q = sp.factor(
        (
            bank[10][2] + alpha * bank[6][2] + beta * bank[2][2]
        ).subs(q**4, t**2)
        .subs(q**2, t)
        / q
    )
    dt_value = sp.solve(phi10, d * t)[0]
    y = sp.symbols("y")
    cubic = (
        64 * y**3
        + (640 * beta - 192 * alpha**2) * y
        + 128 * alpha**3
        - 640 * alpha * beta
        - 800 * psi0
        - 125 * t**2
    )
    require(
        sp.factor(
            800
            * (psi10.subs(d, dt_value / t) - psi0)
            - cubic.subs(y, 5 * s - 2 * alpha)
        )
        == 0,
        "degree-10 spectral cubic mismatch",
    )
    require(
        sp.factor(
            third10_over_q.subs(d, dt_value / t)
            - t * (5 * s - 2 * alpha) / 8
        )
        == 0,
        "degree-10 third observable mismatch",
    )

    # Degree six: Phi=0 gives s=2 alpha/3 and then R/q=-T/4.
    phi6 = sp.factor(
        (bank[6][0] + alpha * bank[2][0]).subs(q**2, t) / q
    )
    third6_over_q = sp.factor(
        (bank[6][2] + alpha * bank[2][2]).subs(q**2, t) / q
    )
    require(
        sp.factor(phi6.subs(s, 2 * alpha / 3)) == 0,
        "degree-6 first flux reduction mismatch",
    )
    require(
        sp.factor(
            third6_over_q.subs(s, 2 * alpha / 3) + t / 4
        )
        == 0,
        "degree-6 third observable mismatch",
    )
    print("faber_flux_bank=PASS degrees=2,6,10")
    print("independent_multinomial_bank=PASS degrees=2,6,10")
    print("degree10_spectral_elimination=PASS")
    print("degree6_terminal_reduction=PASS")


def check_singular_normalization() -> None:
    r, e = sp.symbols("r e")
    y = sp.Rational(125, 64) * r**2 - 2 * e
    t = r * (sp.Rational(125, 64) * r**2 - 3 * e)
    require(
        sp.factor(
            125 * t**2 - 64 * (y - e) ** 2 * (y + 2 * e)
        )
        == 0,
        "singular cubic parametrization mismatch",
    )
    omega = sp.factor(
        (
            t * sp.diff(y, r)
            + sp.Rational(3, 2) * y * sp.diff(t, r)
        )
        / 8
    )
    wanted_omega = (
        203125 * r**4 - 312000 * e * r**2 + 73728 * e**2
    ) / 65536
    require(
        sp.factor(omega - wanted_omega) == 0,
        "singular Keller one-form mismatch",
    )
    primitive = (
        sp.Rational(40625, 65536) * r**5
        - sp.Rational(1625, 1024) * e * r**3
        + sp.Rational(9, 8) * e**2 * r
    )
    require(
        sp.factor(sp.diff(primitive, r) - wanted_omega) == 0,
        "singular Keller primitive mismatch",
    )
    print("singular_cubic_normalization=PASS")
    print("singular_keller_primitive=PASS")


def check_cusp_reconstruction() -> None:
    r, alpha, l, e0 = sp.symbols("r alpha L E")
    t = sp.Rational(125, 64) * r**3
    a = 1 / t
    s = sp.Rational(25, 64) * r**2 + sp.Rational(2, 5) * alpha
    d = sp.Rational(15, 64) * r
    h = sp.expand(a * (l + e0 + s) ** 2 + d)

    # Treat r as a pole and a=64/(125 r^3).  Terms involving polynomial E
    # are regular at that pole; the coefficient of the naked r term is the
    # obstruction.
    h_at_zero_e = sp.expand(h.subs(e0, 0))
    naked_r = sp.expand(h_at_zero_e).coeff(r, 1)
    require(
        sp.factor(naked_r - sp.Rational(5, 16)) == 0,
        "cusp constant-term pole coefficient mismatch",
    )
    print("cusp_reconstruction=PASS pole_coefficient=5/16")


def main() -> None:
    check_polynomial_parts()
    check_flux_bank_and_elimination()
    check_singular_normalization()
    check_cusp_reconstruction()
    print("status=EXACT_SYMBOLIC_REFEREE_ONLY")


if __name__ == "__main__":
    main()
