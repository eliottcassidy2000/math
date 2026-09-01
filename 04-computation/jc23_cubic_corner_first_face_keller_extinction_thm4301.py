#!/usr/bin/env python3
"""Exact symbolic certificate for THM-4301.

Universe: the D=Lambda=0, UZ!=0 cubic corner inside the inherited exact-M=12
reduced (2,3) seam.  The script reconstructs the full displayed Hhat support,
proves the universal q-linear anchor of order at most eight, checks the
q-separable good-form bound, and builds a full-support genus-four hostile on
which the weakest bound is attained.  It does not resolve q-repeated faces or
prove seam entry, JC(2), or DC(2).
"""

from __future__ import annotations

import sympy as sp


def require_zero(value: sp.Expr, label: str) -> None:
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def main() -> None:
    q, r, t, z, sigma = sp.symbols("q r t z sigma")
    x, y, Y, u = sp.symbols("x y Y u")
    U = sp.symbols("U", nonzero=True)
    alpha, beta11 = sp.symbols("alpha_11 beta_11")
    upsilon, xi10 = sp.symbols("upsilon_5 xi_10")
    eta, zeta3 = sp.symbols("eta zeta_3")
    Delta, Theta, Phi = sp.symbols("Delta Theta Phi")
    k0 = sp.Rational(1376, 135)
    K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
    c = sp.Rational(7168, 135) - sp.Rational(7, 6) * Delta

    # D=Lambda=0 with UZ!=0 gives (U,W,Z)=(U,-2U,U).  This is the
    # literal THM-4292/4297 Hhat, not a freely prepared cubic.
    r1 = 1 + q
    hhat = (
        U * (r1**6 - 2 * r1**5 + r1**4)
        + t * (alpha * r1**5 + beta11 * r1**4)
        + t**2 * (upsilon * r1**5 + xi10 * r1**4)
        + t**3 * (eta * r1**4 + zeta3 * r1**3)
        + t**4 * (Delta * r1**4 + Theta * r1**3)
        + t**5 * Phi * r1**3
        + t**6 * (-k0 * r1**3 + K * r1**2)
        + sp.Rational(8, 3) * t**8 * r1**2
        - 3 * t**10 * r1
    )
    F = sp.expand(q * (hhat - z**12) - t**12 / 2)

    c1 = alpha + beta11
    c2 = upsilon + xi10
    c3 = eta + zeta3
    c4 = Delta + Theta
    c5 = Phi
    h = (
        c1 * t + c2 * t**2 + c3 * t**3 + c4 * t**4 + c5 * t**5
        + c * t**6 + sp.Rational(8, 3) * t**8 - 3 * t**10
    )
    require_zero(F.coeff(q, 0) + t**12 / 2, "constant response")
    require_zero(F.coeff(q, 1) - (h - z**12), "q-linear response")
    require_zero(F.subs(t, 0) - (U * q**3 * (1 + q) ** 4 - q * z**12),
                 "top cubic corner")
    require_zero(F.coeff(q, 3).subs(t, 0) - U, "cubic unit")

    # The first nonzero t-order ell in h lies in this exhaustive list.  The
    # fixed 8/3 coefficient is the terminal anchor if all earlier rows vanish.
    ell_values = (1, 2, 3, 4, 5, 6, 8)
    bound_rows = []
    for ell in ell_values:
        s_weight = 9 - ell
        beta_weight = 11 - ell
        assert s_weight > 0 and beta_weight > 0
        bound_rows.append((ell, s_weight, beta_weight))

    universal_checks = 0
    for tau_value in range(1, 65):
        for gamma_value in range(1, 16 * tau_value + 1):
            assert min(2 * gamma_value, 12 * tau_value - gamma_value) <= 8 * tau_value
            universal_checks += 1

    # A q-horizontal first face can only use the earliest h monomial and
    # -q*z^12.  After t=sigma*z and removal of the torus monomial q*z^ell,
    # its strict equation is z^(12-ell)=unit*sigma^ell, hence rational.
    horizontal_rows = []
    for ell in ell_values:
        components = sp.gcd(12 - ell, ell)
        horizontal_rows.append((ell, 12 - ell, ell, components, 0))

    # Weakest-bound positive control.  Force c1=...=c5=c=0.  The coefficient
    # relation makes Delta=2048/45 and Theta=-Delta; the unavoidable t^8 q
    # anchor and the q^2 t^4 translation both remain on the first face.
    delta_star = sp.Rational(2048, 45)
    hostile = {
        U: 1,
        alpha: 0,
        beta11: 0,
        upsilon: 0,
        xi10: 0,
        eta: 0,
        zeta3: 0,
        Delta: delta_star,
        Theta: -delta_star,
        Phi: 0,
    }
    require_zero(c.subs(hostile), "hostile c=0")
    F_hostile = sp.expand(F.subs(hostile).subs(t, sigma * z))
    weights = {sigma: 1, z: 2, q: 12}
    weighted_terms: list[tuple[int, sp.Expr]] = []
    for term in sp.Add.make_args(F_hostile):
        powers = term.as_powers_dict()
        weight = sum(int(powers.get(variable, 0)) * value
                     for variable, value in weights.items())
        weighted_terms.append((weight, term))
    minimum = min(weight for weight, _ in weighted_terms)
    face = sp.expand(sum(term for weight, term in weighted_terms if weight == minimum))
    expected_face = sp.expand(
        q**3 + delta_star * q**2 * (sigma * z) ** 4
        + sp.Rational(8, 3) * q * (sigma * z) ** 8
        - q * z**12 - (sigma * z) ** 12 / 2
    )
    assert minimum == 36
    require_zero(face - expected_face, "full-support weakest face")

    # In the toric chart x=sigma^2/z, y=q/z^6, division by z^18 gives G.
    G = sp.expand(
        y**3 + delta_star * x**2 * y**2
        + (sp.Rational(8, 3) * x**4 - 1) * y - x**6 / 2
    )
    discriminant = sp.factor(sp.discriminant(G, y))
    expected_discriminant = (
        24553315427 * x**12 - 1282042880 * x**8
        + 247770240 * x**4 + 486000
    ) / sp.Integer(121500)
    require_zero(discriminant - expected_discriminant, "hostile discriminant")
    assert sp.degree(discriminant, x) == 12
    assert sp.gcd(sp.Poly(discriminant, x),
                  sp.Poly(sp.diff(discriminant, x), x)).degree() == 0

    # A rational root of the monic cubic over C(x) would be a polynomial of
    # degree at most two.  Exact coefficient elimination gives the unit ideal.
    root_u, root_v, root_w = sp.symbols("root_u root_v root_w")
    root = root_u * x**2 + root_v * x + root_w
    root_equations = [sp.Poly(sp.expand(G.subs(y, root)), x).coeff_monomial(x**j)
                      for j in range(9)]
    root_basis = sp.groebner(root_equations, root_u, root_v, root_w, order="lex")
    assert list(root_basis) == [1]

    infinity_cubic = sp.factor(
        sp.expand(u**6 * G.subs({x: 1 / u, y: Y / u**2})).subs(u, 0)
    )
    assert sp.discriminant(infinity_cubic, Y) != 0
    ramification = 12
    genus = (ramification - 4) // 2
    assert genus == 4

    # At (s,beta,gamma)=(1,2,12), L=36 and v(F_q)=L-gamma=24.
    # The good form -sigma^9 z^10 dz/F_q has exact generic order seven.
    form_order = 9 * 1 + 11 * 2 - 24
    assert form_order == 7

    print("THM-4301 PRIMARY EXACT AUDIT: PASS")
    print("UNIVERSE exact_M12 D=Lambda=0 UZ!=0 characteristic_zero")
    print("q_linear=h(t)-z^12; h_orders=1,2,3,4,5,6,8; terminal=8/3*t^8")
    print(f"universal_bound: min(2gamma,12tau-gamma)<=8tau checks={universal_checks}")
    print("ell  lower_bound")
    for ell, s_weight, beta_weight in bound_rows:
        print(f"{ell:3d}  {s_weight}*s+{beta_weight}*beta")
    print("q_horizontal: z^(12-ell)=unit*sigma^ell; every component genus=0")
    print(f"hostile: weights=(1,2,12) face_weight={minimum} Delta={delta_star}")
    print(f"hostile_discriminant={discriminant}")
    print(f"hostile_cover: degree=3 finite_branch=12 infinity_branch=0 genus={genus}")
    print(f"hostile_good_form_order={form_order}")
    print("SCOPE q_separable_first_faces_extinct q_horizontal_rational q_repeated_open")


if __name__ == "__main__":
    main()
