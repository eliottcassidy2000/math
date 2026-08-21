#!/usr/bin/env python3
"""Exact controls for the nodal-cylinder Keller-projection scaffold.

The source map

    (u,t) -> (u, X=t^2-1, Y=t(t^2-1))

is the normalization of a nodal-cubic cylinder.  It is differential-rank two
and identifies (u,1) with (u,-1).  This companion verifies that its descended
normal-multiplier gate and even a closed ambient two-form gate pass, records
why the displayed closed certificate is not a wedge of two polynomial
potentials, and checks the first degree-108 leading-face scaffold.

The final scaffold is not a Keller pair.  Its purpose is to isolate the first
globally degree-compatible correction cell, not to assert a counterexample.
All truth gates remain active under ``python -O``.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


def require(condition: bool, label: str) -> None:
    """Raise on a failed truth gate, including under optimized execution."""
    if not bool(condition):
        raise ArithmeticError(f"FAILED: {label}")


def top_form(poly: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    """Return the ordinary top homogeneous form in the named variables."""
    terms = sp.Poly(sp.expand(poly), *variables).terms()
    degree = max(sum(exponent) for exponent, _ in terms)
    return sp.expand(sum(
        coefficient * sp.prod(variable**power for variable, power
                              in zip(variables, exponent))
        for exponent, coefficient in terms if sum(exponent) == degree
    ))


def source_degree_support(target_cap: int, source_degree: int) -> list[int]:
    """Possible u-exponents in a fixed source-degree top row.

    Target weights are deg(u,X,Y)=(1,2,3) after pullback.
    """
    exponents: set[int] = set()
    for a, b, c in product(range(target_cap + 1), repeat=3):
        if a + b + c <= target_cap and a + 2 * b + 3 * c == source_degree:
            exponents.add(a)
    return sorted(exponents)


def main() -> None:
    u, t, lam = sp.symbols("u t lambda")
    U, X0, Y0 = sp.symbols("U X Y")

    X = t**2 - 1
    Y = t * X
    H = Y0**2 - X0**2 * (X0 + 1)
    require(sp.expand(H.subs({X0: X, Y0: Y})) == 0,
            "nodal-cylinder image equation")

    # Tangent minors in target coordinate order (U,X,Y).
    m_ux = sp.diff(X, t)
    m_uy = sp.diff(Y, t)
    m_xy = sp.Integer(0)
    require(m_ux == 2 * t and m_uy == 3 * X + 2,
            "tangent minors")
    require(sp.expand(((3 * X + 1) / 2) * m_uy
                      - sp.Rational(9, 4) * Y * m_ux) == 1,
            "descended unit tangent-minor certificate")

    # The image normal is X times the tangent normal.
    h_x = sp.diff(H, X0)
    h_y = sp.diff(H, Y0)
    require(sp.expand(h_x.subs({X0: X, Y0: Y}) + X * m_uy) == 0,
            "H_X pullback has multiplier X")
    require(sp.expand(h_y.subs({X0: X, Y0: Y}) - X * m_ux) == 0,
            "H_Y pullback has multiplier X")
    require(sp.rem(
        sp.expand(
            X0
            + ((3 * X0 + 1) / 2) * h_x
            + sp.Rational(9, 4) * Y0 * h_y
        ),
        H,
        Y0,
    ) == 0, "normal multiplier lies in the image Jacobian ideal")

    # A closed ambient two-form with the same unit contraction.
    q_ux = -sp.Rational(9, 4) * Y0
    q_uy = (3 * X0 + 1) / 2
    q_xy = sp.Rational(15, 4) * U
    closure = (sp.diff(q_ux, Y0) - sp.diff(q_uy, X0)
               + sp.diff(q_xy, U))
    require(closure == 0, "ambient coefficient two-form is closed")
    require(sp.expand(
        q_ux.subs({X0: X, Y0: Y}) * m_ux
        + q_uy.subs({X0: X, Y0: Y}) * m_uy
        + q_xy.subs(U, u) * m_xy
    ) == 1, "closed certificate still contracts to one")

    # Its dual divergence-free vector field becomes diagonal after V=3X+1:
    # 4D = 15 U d_U - 6 V d_V - 9 Y d_Y.  Every invariant monomial
    # U^a V^b Y^c satisfies 5a=2b+3c.  At U=0 this forces b=c=0, so two
    # invariant gradients are parallel there; they cannot have cross product D,
    # which is generically nonzero on U=0.  We check the weight implication in
    # a finite hostile box; the proof in the accompanying note is all-degree.
    invariant_without_u = [
        (b, c) for b in range(12) for c in range(12)
        if -6 * b - 9 * c == 0
    ]
    require(invariant_without_u == [(0, 0)],
            "diagonal invariant restriction on U=0")

    # Degree-108 top-row support at the first four target caps.
    cap_support = {
        cap: source_degree_support(cap, 108) for cap in range(35, 39)
    }
    require(cap_support == {35: [], 36: [0], 37: [0, 1],
                            38: [0, 1, 2, 3]},
            "degree-108 target-cap support")

    # A cap-38 lift of the first admissible common leading base
    # K=t^35(t+lambda*u).  The last term uses X^3=Y^2-X^2 on the nodal image;
    # omitting the second summand would create a spurious u-width defect.
    p_cap38 = (
        Y**24
        + 2 * lam * u * X * Y**23
        + lam**2 * u**2 * X**2 * Y**22
    )
    q_cap38 = (
        Y**36
        + 3 * lam * u * X * Y**35
        + 3 * lam**2 * u**2 * X**2 * Y**34
        + lam**3 * u**3 * (Y**35 - X**2 * Y**33)
    )
    common_base = t**35 * (t + lam * u)
    require(top_form(p_cap38, (u, t)) == sp.expand(common_base**2),
            "degree-72 common-base square")
    require(top_form(q_cap38, (u, t)) == sp.expand(common_base**3),
            "degree-108 common-base cube")
    require(sp.Poly(p_cap38, u, t).total_degree() == 72,
            "square scaffold source degree")
    require(sp.Poly(q_cap38, u, t).total_degree() == 108,
            "cube scaffold source degree")

    a2 = sp.Poly(p_cap38, u).coeff_monomial(u**2)
    b3 = sp.Poly(q_cap38, u).coeff_monomial(u**3)
    width_defect = sp.expand(2 * a2 * sp.diff(b3, t)
                             - 3 * sp.diff(a2, t) * b3)
    require(width_defect == 0, "cap-38 u-width resonance")

    # The cap-38 representative is exactly C^2,C^3 after pullback, even though
    # the unreduced ambient cube would have target degree 39.
    c_resonant = Y**12 + lam * u * X * Y**11
    p_skeleton = p_cap38
    q_skeleton = q_cap38
    require(sp.expand(p_skeleton - c_resonant**2) == 0,
            "cap-38 square equals literal square")
    require(sp.expand(q_skeleton - c_resonant**3) == 0,
            "cap-38 relation-reduced cube equals literal cube")

    c_target = Y0**12 + lam * U * X0 * Y0**11
    p_target = sp.expand(c_target**2)
    q_target = sp.expand(
        Y0**36
        + 3 * lam * U * X0 * Y0**35
        + 3 * lam**2 * U**2 * X0**2 * Y0**34
        + lam**3 * U**3 * (Y0**35 - X0**2 * Y0**33)
    )
    require(sp.expand(q_target - c_target**3
                      - lam**3 * U**3 * Y0**33 * H) == 0,
            "ambient cap reduction uses exactly the image relation")
    target_caps = tuple(sp.Poly(poly, U, X0, Y0).total_degree()
                        for poly in (p_target, q_target,
                                     c_target, c_target**3))
    require(target_caps == (26, 38, 13, 39),
            "target cap ledger before and after relation reduction")
    require(top_form(p_skeleton, (u, t)) == sp.expand(common_base**2),
            "cap-38 resonant square top")
    require(top_form(q_skeleton, (u, t)) == sp.expand(common_base**3),
            "cap-38 resonant cube top")
    require(sp.expand(sp.diff(p_skeleton, u) * sp.diff(q_skeleton, t)
                      - sp.diff(p_skeleton, t) * sp.diff(q_skeleton, u)) == 0,
            "literal square/cube bracket vanishes")

    # Gwozdziewicz's injective-line theorem forbids the tempting boundary
    # curve (u,0).  Use instead an immersed nodal boundary, so u=+/-1 and
    # t=+/-1 give a four-point fibre.  The Y-rows repair both transverse jets.
    a0 = u**2 - 1 + Y / 2 + p_skeleton
    b0 = u * (u**2 - 1) + sp.Rational(3, 4) * u * Y + q_skeleton
    jacobian = sp.expand(sp.diff(a0, u) * sp.diff(b0, t)
                         - sp.diff(a0, t) * sp.diff(b0, u))
    require(sp.expand(jacobian.subs(t, 1)) == 1
            and sp.expand(jacobian.subs(t, -1)) == 1,
            "unit Jacobian on both collision sheets")
    require(sp.expand(a0.subs(t, 1) - a0.subs(t, -1)) == 0
            and sp.expand(b0.subs(t, 1) - b0.subs(t, -1)) == 0,
            "collision is retained")
    four_point_values = {
        (sp.expand(a0.subs({u: sign_u, t: sign_t})),
         sp.expand(b0.subs({u: sign_u, t: sign_t})))
        for sign_u in (-1, 1) for sign_t in (-1, 1)
    }
    require(four_point_values == {(0, 0)}, "explicit four-point fibre")
    require(sp.expand(jacobian.subs({u: 0, t: 0}) - 1) != 0,
            "displayed scaffold is not Keller")

    # The displayed scaffold still fails Gwozdziewicz globally.  On t=c,
    # write C=beta(c)+alpha(c)u.  Since alpha is nonconstant, alpha(c)=i has
    # a complex solution.  At every root of alpha^2+1, beta and 2alpha beta
    # are nonzero, so A0 becomes a nonconstant affine-linear function of u
    # and the restriction of the pair to that horizontal line is injective.
    alpha_line = sp.expand(lam * X * Y**11)
    beta_line = sp.expand(Y**12)
    coefficient_field = sp.QQ.frac_field(lam)
    hostile_gcd = sp.gcd(
        sp.Poly(alpha_line**2 + 1, t, domain=coefficient_field),
        sp.Poly(alpha_line * beta_line, t, domain=coefficient_field),
    )
    require(hostile_gcd.degree() == 0,
            "injective-horizontal-line hostile has nonzero linear term")

    print("NODAL-CYLINDER KELLER-PROJECTION SCAFFOLD")
    print("source_packet=(u,X=t^2-1,Y=t(t^2-1))")
    print("image_relation=Y^2-X^2(X+1)")
    print("collision=(u,1)~(u,-1); differential_rank_two=True")
    print("normal_multiplier=X; multiplier_in_image_Jacobian_ideal=True")
    print("descended_unit_minor_certificate=True")
    print("closed_ambient_two_form_certificate=True")
    print("displayed_closed_certificate_polynomial_potential_realizable=False")
    print(f"degree108_u_support_by_target_cap={cap_support}")
    print("first_face_admissible_cap=38")
    print("first_exact_resonant_one_scale_cap=38")
    print("reduced_source_degree_pair=(72,108)")
    print("leading_forms=(t^35(t+lambda*u))^2,(t^35(t+lambda*u))^3")
    print("collision_first_jet_repaired=True")
    print("boundary_line_restriction_noninjective=True")
    print("four_point_fibre=((u,t)=(+/-1,+/-1))->(0,0)")
    print("displayed_scaffold_has_injective_horizontal_line=True")
    print("displayed_scaffold_Keller=False")
    print("VERDICT=typed open correction cell, not a counterexample")


if __name__ == "__main__":
    main()
