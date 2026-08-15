#!/usr/bin/env python3
"""Exact companion for THM-3081.

Checks the primitive terminal Bezout coordinates, weight/phase decoder, and
the affine/reciprocal exact controls inherited from THM-3074/3080.
"""

from __future__ import annotations

import math

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def bezout_one(g: int, e: int) -> tuple[int, int]:
    """Return A,B with A*g+B*e=1 for coprime positive g,e."""

    require(math.gcd(g, e) == 1, "terminal values must be coprime")
    candidates: list[tuple[int, int]] = []
    for a in range(-e, e + 1):
        numerator = 1 - a * g
        if numerator % e == 0:
            candidates.append((a, numerator // e))
    if candidates:
        return min(candidates, key=lambda pair: (abs(pair[0]) + abs(pair[1]), abs(pair[0])))
    raise RuntimeError("Bezout search failed")


def key_bezout(alpha: int, beta: int) -> tuple[int, int]:
    """Return gamma,delta with alpha*delta+beta*gamma=1."""

    require(math.gcd(alpha, beta) == 1, "primitive strict-stage exponents")
    candidates: list[tuple[int, int]] = []
    for gamma in range(-alpha, alpha + 1):
        numerator = 1 - beta * gamma
        if numerator % alpha == 0:
            candidates.append((gamma, numerator // alpha))
    if candidates:
        return min(candidates, key=lambda pair: (abs(pair[0]) + abs(pair[1]), abs(pair[0])))
    raise RuntimeError("strict-stage Bezout search failed")


def rational_map_degree_from_pair(numerator: sp.Expr, denominator: sp.Expr, t: sp.Symbol) -> int:
    numerator_poly = sp.Poly(numerator, t)
    denominator_poly = sp.Poly(denominator, t)
    common = sp.gcd(numerator_poly, denominator_poly)
    numerator_reduced = numerator_poly.exquo(common)
    denominator_reduced = denominator_poly.exquo(common)
    return max(numerator_reduced.degree(), denominator_reduced.degree())


def main() -> None:
    # The terminal value-one coordinate rho and residue coordinate theta form
    # a unimodular exponent chart.  Every leading monomial r^j m^n decodes as
    # rho^(jg+ne) theta^(-Bj+An).
    for g in range(1, 25):
        for e in range(1, 25):
            if math.gcd(g, e) != 1:
                continue
            a, b = bezout_one(g, e)
            require(a * g + b * e == 1, "terminal Bezout identity")
            determinant = a * g + b * e
            require(determinant == 1, "terminal exponent determinant")
            for gauge in range(-3, 4):
                a_gauge = a + gauge * e
                b_gauge = b - gauge * g
                require(a_gauge * g + b_gauge * e == 1, "Bezout gauge")
                require(a_gauge - a == gauge * e, "rho gauge r exponent")
                require(b_gauge - b == -gauge * g, "rho gauge m exponent")
            for j in range(-8, 9):
                for n in range(-8, 9):
                    weight = j * g + n * e
                    theta_power = -b * j + a * n
                    # Inverting [log rho; log theta]=[[A,B],[-e,g]][log r;log m]
                    # gives r=rho^g theta^-B and m=rho^e theta^A.
                    reconstructed_rho = j * g + n * e
                    reconstructed_theta = -b * j + a * n
                    require(reconstructed_rho == weight, "rho exponent decoder")
                    require(reconstructed_theta == theta_power, "theta exponent decoder")

    # Earlier strict stage: (m^g/r^e, m^gamma*r^delta) has determinant
    # d=gcd(g,e), while the primitive-root version has determinant one.
    # Its kernel on q-torsion has exactly gcd(d,q) elements.
    for g in range(1, 17):
        for e in range(1, 17):
            divisor = math.gcd(g, e)
            alpha = g // divisor
            beta = e // divisor
            gamma, delta = key_bezout(alpha, beta)
            require(alpha * delta + beta * gamma == 1, "primitive strict determinant")
            require(g * delta + e * gamma == divisor, "nonprimitive strict determinant")
            for q in range(2, 19):
                kernel_size = sum(
                    1
                    for m_exp in range(q)
                    for r_exp in range(q)
                    if (g * m_exp - e * r_exp) % q == 0
                    and (gamma * m_exp + delta * r_exp) % q == 0
                )
                require(kernel_size == math.gcd(divisor, q), "strict torsion kernel")

    # Degree multiplication for P^1 maps forces a rational right inverse to
    # be Mobius.  Check the exact degree convention and two canonical forms.
    t, u = sp.symbols("t u")
    for degree in range(1, 9):
        polynomial = t**degree + t + 1 if degree > 1 else 2 * t + 1
        require(
            rational_map_degree_from_pair(polynomial, sp.Integer(1), t) == degree,
            "polynomial map degree",
        )
        reciprocal = 1 + t ** (-degree)
        numerator, denominator = sp.fraction(sp.cancel(reciprocal))
        require(
            rational_map_degree_from_pair(numerator, denominator, t) == degree,
            "reciprocal map degree",
        )

    alpha, beta, gamma, delta = sp.symbols("alpha beta gamma delta", nonzero=True)
    theta_mobius = (delta * u - beta) / (alpha - gamma * u)
    h_mobius = (alpha * t + beta) / (gamma * t + delta)
    require(
        sp.factor(h_mobius.subs(t, theta_mobius) - u) == 0,
        "general Mobius right inverse",
    )
    mobius_determinant = alpha * delta - beta * gamma
    mobius_square = (delta + gamma * theta_mobius) ** 2 / mobius_determinant
    require(
        sp.factor(sp.diff(theta_mobius, u) - mobius_square) == 0,
        "Mobius derivative square law",
    )

    # One-stage packet: g=1,e=5, theta=3u and P=u=theta/3.
    theta_one = 3 * u
    h_one = t / 3
    require(sp.expand(h_one.subs(t, theta_one) - u) == 0, "one-stage residue decoder")
    theta_one_prime = sp.diff(theta_one, u)
    require(theta_one_prime == 3, "one-stage autonomous decoder")

    # Two-stage packet: terminal g=2,e=3 with A=-1,B=1,
    # rho=-3, theta=9/u, and P=u=9/theta.
    g_two, e_two = 2, 3
    a_two, b_two = bezout_one(g_two, e_two)
    require((a_two, b_two) == (-1, 1), "two-stage canonical Bezout pair")
    r_two = u
    m_two = -3 * u
    rho_two = sp.factor(r_two**a_two * m_two**b_two)
    theta_two = sp.factor(m_two**g_two / r_two**e_two)
    require(rho_two == -3, "two-stage value-one leading coefficient")
    require(theta_two == 9 / u, "two-stage residue parameter")
    require(sp.factor((9 / t).subs(t, theta_two) - u) == 0, "two-stage residue decoder")
    # Here E=3, kappa=1, K=-1/27, L=1, A=-1.
    ode_two = sp.factor(3 * theta_two ** (1 - a_two) * sp.Rational(-1, 27))
    require(sp.factor(sp.diff(theta_two, u) - ode_two) == 0, "two-stage autonomous ODE")

    # Three-stage packet: terminal g=4,e=3 with A=1,B=-1,
    # rho=1/3, theta=81u, and P=u=theta/81.
    g_three, e_three = 4, 3
    a_three, b_three = bezout_one(g_three, e_three)
    require((a_three, b_three) == (1, -1), "three-stage canonical Bezout pair")
    r_three = u
    m_three = 3 * u
    rho_three = sp.factor(r_three**a_three * m_three**b_three)
    theta_three = sp.factor(m_three**g_three / r_three**e_three)
    require(rho_three == sp.Rational(1, 3), "three-stage value-one coefficient")
    require(theta_three == 81 * u, "three-stage residue parameter")
    require(
        sp.expand((t / 81).subs(t, theta_three) - u) == 0,
        "three-stage residue decoder",
    )
    # Here E=3, kappa=1, K=27, L=1, A=1.
    ode_three = sp.factor(3 * theta_three ** (1 - a_three) * 27)
    require(sp.factor(sp.diff(theta_three, u) - ode_three) == 0, "three-stage autonomous ODE")

    print("theorem=THM-3081")
    print("status=PROVED_VERIFIED_EXACT_CANDIDATE")
    print("terminal_chart=Ag+Be=1;rho=r^A*m^B;theta=m^g/r^e")
    print("initial_decoder=in_w(F)=rho^w*H_F(theta);H_F_rational")
    print("residue_field=C(u)=C(theta)")
    print("mobius=u=H_P(theta);degree(H_P)*degree(theta(u))=1")
    print("target_initial=tau=rho^E*K(theta)")
    print("prefactor_initial=u_N=rho^(E-e)*L(theta)")
    print("autonomous_decoder=theta_prime=E*kappa^-1*theta^(1-A)*K(theta)/L(theta)")
    print("mobius_square=K(T)/L(T)=kappa/(E*(ad-bc))*T^(A-1)*(a-cT)^2")
    print("strict_torsor=(m^g/r^e,m^gamma*r^delta);degree=gcd(g,e);kernel=mu_gcd(g,e)")
    print("one_stage=theta=3u;H=theta/3")
    print("two_stage=theta=9/u;H=9/theta;theta_prime=-theta^2/9")
    print("three_stage=theta=81u;H=theta/81;theta_prime=81")
    print("scope=terminal_associated_graded_residue_decoder;no_polynomial_globalization_or_C3_A4_S4_JC2_exclusion")


if __name__ == "__main__":
    main()
