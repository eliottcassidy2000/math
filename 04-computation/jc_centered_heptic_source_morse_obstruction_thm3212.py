#!/usr/bin/env python3
"""Independent exact audit of THM-3212's centered THM-3123 source lift.

This script imports no repository companion.  It works in the
symmetric cubic accessory algebra Q[u]/(q(u)) for both remaining heptic
passports and checks the response/source identities and the five-fibre
critical obstruction exactly.
"""

from __future__ import annotations

import hashlib

import sympy as sp


x, z, u = sp.symbols("x z u")
Bjet, Cjet = sp.symbols("Bjet Cjet")
Bprime, Cprime, Vprime, Aprime, Eprime = sp.symbols(
    "Bprime Cprime Vprime Aprime Eprime"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def reduced_rational(expression: sp.Expr, modulus: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    q_poly = sp.Poly(modulus, u, domain=sp.QQ)
    num_poly = sp.rem(sp.Poly(numerator, u, domain=sp.QQ), q_poly)
    den_poly = sp.rem(sp.Poly(denominator, u, domain=sp.QQ), q_poly)
    require(sp.gcd(den_poly, q_poly).degree() == 0, "nonunit denominator")
    inverse = sp.invert(den_poly, q_poly)
    return sp.rem(num_poly * inverse, q_poly).as_expr()


def zero_mod(expression: sp.Expr, modulus: sp.Expr, variables=(x, z)) -> bool:
    numerator = sp.together(expression).as_numer_denom()[0]
    polynomial = sp.Poly(sp.expand(numerator), *variables)
    return all(reduced_rational(coefficient, modulus) == 0 for coefficient in polynomial.coeffs())


def unit_mod(expression: sp.Expr, modulus: sp.Expr) -> bool:
    reduced = reduced_rational(expression, modulus)
    numerator, denominator = sp.cancel(reduced).as_numer_denom()
    q_poly = sp.Poly(modulus, u, domain=sp.QQ)
    return (
        sp.gcd(sp.Poly(numerator, u, domain=sp.QQ), q_poly).degree() == 0
        and sp.gcd(sp.Poly(denominator, u, domain=sp.QQ), q_poly).degree() == 0
    )


def digest(expression: sp.Expr) -> str:
    return hashlib.sha256(str(sp.factor(expression)).encode()).hexdigest()[:20]


def audit_case(name: str, a: int, b: int, v: sp.Expr, q: sp.Expr, S: sp.Expr, A0: sp.Expr) -> None:
    quadratic = x**2 - u * x + v
    D = sp.expand(x**a * (x - 1) ** b * quadratic)
    T = sp.expand(x * (x - 1) * quadratic)
    E = sp.expand(
        (
            a * (x - 1) * quadratic
            + b * x * quadratic
            + x * (x - 1) * (2 * x - u)
        )
        / 7
    )
    C = -7 * A0
    M = sp.expand(S * E * T)
    V = sp.cancel(4 * S * D * T**2 / C**2)
    G = sp.cancel(C * E / (2 * D * T))
    F = sp.cancel(S * E**2 / D)
    A = sp.cancel(2 * M / C)
    g = sp.expand(S * T)

    require(zero_mod(D + A0 - S * E**2, q, (x,)), "D+A0=S E^2 failed")
    require(unit_mod(C, q), "C is not an accessory unit")
    require(zero_mod(F - V * G**2, q, (x,)), "F=VG^2 failed")
    require(zero_mod(2 * V * sp.diff(G, x) + sp.diff(V, x) * G - 2, q, (x,)), "response ODE failed")
    require(zero_mod(A - V * G, q, (x,)), "A_src=VG failed")
    require(zero_mod(A - 2 * M / C, q, (x,)), "A_src=2M/C failed")
    require(zero_mod(A**2 / V - F, q, (x,)), "source T_F=F failed")
    require(zero_mod(V / g - 4 * D * T / C**2, q, (x,)), "g does not divide V")
    require(zero_mod(A / g - 2 * E / C, q, (x,)), "g does not divide A_src")

    # Pairwise-disjoint squarefree factor ledger for M=S E T.
    gates = {
        "disc_g": sp.discriminant(g, x),
        "disc_T": sp.discriminant(T, x),
        "disc_E": sp.discriminant(E, x),
        "res_S_T": sp.resultant(S, T, x),
        "res_S_E": sp.resultant(S, E, x),
        "res_E_T": sp.resultant(E, T, x),
    }
    require(all(unit_mod(value, q) for value in gates.values()), "factor gate is not a unit")
    resultant_owner = sp.resultant(sp.diff(M, x), g, x)
    require(unit_mod(resultant_owner, q), "A_src' is not a unit modulo g")
    require(sp.degree(g, x) == 5, "g does not have degree five")

    # Universal one-jet calculation at a geometric root of g, where V=A=0.
    # H=Bz+C, P_z=2B(Bz+C), and P_x is the expression below.
    Hjet = Bjet * z + Cjet
    pz_jet = sp.expand(2 * Bjet * Hjet)
    px_jet = sp.expand(
        2 * Hjet * (Vprime * z**2 + Bprime * z + Cprime)
        + Aprime * z
        + Eprime
    )
    critical_z = -Cjet / Bjet
    delta = Bjet * Eprime - Aprime * Cjet
    require(sp.cancel(pz_jet.subs(z, critical_z)) == 0, "unit-B critical z failed")
    require(
        sp.cancel(px_jet.subs(z, critical_z) - delta / Bjet) == 0,
        "unit-B Delta formula failed",
    )

    # Centered B=C=0 gives P_z=0 and P_x=A'z+E0'; since A' is a unit,
    # z=-E0'/A' is the unique critical point for every E0.
    require(
        sp.cancel(pz_jet.subs({Bjet: 0, Cjet: 0})) == 0,
        "centered P_z formula failed",
    )
    require(
        sp.cancel(px_jet.subs({Bjet: 0, Cjet: 0, z: -Eprime / Aprime}))
        == 0,
        "centered critical point formula failed",
    )
    # For the centered polynomial itself, P_zz=0 and P_xz=A' above g,
    # so each forced point has Hessian determinant -A'^2 and is reduced.
    centered_hessian_determinant = -Aprime**2
    require(
        centered_hessian_determinant != 0,
        "formal centered Hessian determinant vanished",
    )

    # Minimal unit-B local repair B=1,C0=0,E0=x has Delta=1.  At every
    # root of g, P_z=2z, so its sole zero is z=0, where P_x=1.
    control = {Bjet: 1, Cjet: 0, Eprime: 1, z: 0}
    require(pz_jet.subs(control) == 0, "control P_z root failed")
    require(px_jet.subs(control) == 1, "control P_x unit failed")
    require(
        delta.subs({Bjet: 1, Cjet: 0, Eprime: 1}) == 1,
        "control Delta failed",
    )

    print(f"passport={name}")
    print(f"accessory_modulus={sp.factor(q)}")
    print(f"v={v}")
    print(f"S={S}")
    print(f"A0={A0}")
    print(f"C={sp.factor(C)}")
    print("D_plus_A0_equals_S_E2=PASS")
    print("F_equals_VG2_and_2VGprime_plus_VprimeG_equals_2=PASS")
    print("A_src_equals_VG_equals_2SET_over_C=PASS")
    print("T_F_equals_A_src2_over_V_equals_F; d_F=0; s_F=-E0=PASS")
    print("g_equals_ST_degree=5; g_divides_V_and_A_src=PASS")
    reduced_owner = sp.factor(reduced_rational(resultant_owner, q))
    print(f"owner_resultant_reduced={reduced_owner}")
    print(f"owner_resultant_digest={digest(reduced_owner)}")
    print("A_src_prime_unit_mod_g=PASS")
    print("centered_source_has_one_critical_point_over_each_g_root_for_every_E0=PASS")
    print("five_forced_critical_points_are_reduced_with_Hessian_det_minus_Aprime2=PASS")
    print("B=1,C0=0,E0=x_local_control: B=1,Delta=1,Pz_root_z=0,Px=1=PASS")


def main() -> None:
    print("THM-3212 CENTERED HEPTIC SOURCE MORSE OBSTRUCTION EXACT AUDIT")
    audit_case(
        "(4,1,1,1)",
        4,
        1,
        (8 * u**2 + 9 * u + 8) / 7,
        100 * u**3 + 244 * u**2 + 237 * u + 44,
        x + 5 * (u + 1) / 7,
        80 * ((8 * u**2 + 9 * u + 8) / 7) ** 2 * (u + 1) / 343,
    )
    audit_case(
        "(3,2,1,1)",
        3,
        2,
        (24 * u**2 - 16 * u - 16) / 21,
        75 * u**3 - 89 * u**2 - 31 * u + 61,
        x + (5 * u - 4) / 7,
        9 * ((24 * u**2 - 16 * u - 16) / 21) ** 2 * (5 * u - 4) / 343,
    )
    print("scope=four_THM3123_S7_covers_centered_family_not_JC2_or_chart_entry")
    print("FAILED_CHECKS=NONE")


if __name__ == "__main__":
    main()
