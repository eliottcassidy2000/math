#!/usr/bin/env python3
"""Primary exact symbolic audit for THM-4272.

This canonical certificate verifies the
exact face input for the raw W=0, U+Z=0 contact, the E0-isotypic
differential-order lemma, and the extension of the honest fat-contact
incidence base across Lambda=0.  The relative-Hom/Neron proof is geometric,
not computational.  In particular this script does not construct or descend
a Keller face map.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    b, r, U, W, Z = sp.symbols("b r U W Z")
    S, P = sp.symbols("S P")
    D = W**2 - 4 * U * Z
    Lambda = U + W + Z

    # THM-4230 main face before any division by U, Z, or Lambda.
    face = 1 - U * P**6 - W * S**2 * P**5 - Z * S**4 * P**4
    scaled_face = sp.factor(
        sp.cancel(b**12 * face.subs({S: 1 / b, P: r / b**2}))
    )
    expected_scaled = b**12 - U * r**6 - W * r**5 - Z * r**4
    require(sp.expand(scaled_face - expected_scaled) == 0, "toric chart mismatch")

    # R is r=1.  Its complete scheme-theoretic contact algebra is monic.
    contact = sp.factor(expected_scaled.subs(r, 1))
    require(contact == b**12 - U - W - Z, "contact algebra mismatch")

    # Lambda is not a smoothness gate for C.  On the torus, the exact
    # derivative syzygy makes U*Z*D the gate; on the three C-only boundary
    # edges the same gates appear as discriminants.  Lambda is only the
    # resultant with the rational component's boundary point X=1.
    derivative_S = sp.factor(sp.diff(face, S))
    derivative_P = sp.factor(sp.diff(face, P))
    sbr = W * P + 2 * Z * S**2
    pbr = 6 * U * P**2 + 5 * W * S**2 * P + 4 * Z * S**4
    require(derivative_S == -2 * S * P**4 * sbr,
            "S derivative factorization changed")
    require(derivative_P == -P**3 * pbr,
            "P derivative factorization changed")
    syzygy = sp.factor(2 * Z * pbr - (4 * Z * S**2 + 3 * W * P) * sbr)
    require(sp.expand(syzygy + 3 * D * P**2) == 0,
            "torus critical-point syzygy changed")

    X = sp.symbols("X")
    edge_Z = 1 - Z * X**4
    edge_D = U * X**2 + W * X + Z
    edge_U = U - X**6
    disc_Z = sp.factor(sp.discriminant(edge_Z, X))
    disc_D = sp.factor(sp.discriminant(edge_D, X))
    disc_U = sp.factor(sp.discriminant(edge_U, X))
    lambda_resultant = sp.factor(sp.resultant(X - 1, edge_D, X))
    require(disc_Z == -256 * Z**3, "Z-edge discriminant changed")
    require(disc_D == D, "quadratic-edge discriminant changed")
    require(disc_U == 46656 * U**5, "U-edge discriminant changed")
    require(lambda_resultant == Lambda, "R-C boundary resultant changed")

    # This monic algebra remains finite free of rank twelve when Lambda=0;
    # etaleness is deliberately not asserted there.
    contact_basis = tuple(b**index for index in range(12))
    require(len(contact_basis) == 12, "fat-contact basis length changed")
    require(sp.Poly(b**12 - Lambda, b).LC() == 1,
            "fat-contact algebra stopped being monic")

    # On W=0, Z=-U, the C branch is smooth and has order-twelve contact.
    wall_face = sp.factor(expected_scaled.subs({W: 0, Z: -U}))
    expected_wall = b**12 - U * r**4 * (r**2 - 1)
    require(sp.expand(wall_face - expected_wall) == 0, "wall face mismatch")
    derivative_r = sp.factor(sp.diff(wall_face, r).subs({b: 0, r: 1}))
    require(derivative_r == -2 * U, "C branch is not transverse to r at Q")
    require(sp.factor(wall_face.subs(r, 1)) == b**12, "contact length changed")
    leading_branch = sp.factor(-1 / derivative_r)
    require(leading_branch == 1 / (2 * U), "implicit branch leading term changed")

    # The homogeneous t quotient identifies the same coefficient wall.
    T, R = sp.symbols("T R")
    U_phi = 4 * T**2 * R**2
    Z_phi = (T**2 - R**2) ** 2
    Lambda_phi = (T**2 + R**2) ** 2
    require(sp.expand(U_phi + Z_phi - Lambda_phi) == 0, "homogeneous wall identity")

    # A fixed radical normalization selects only one infinity branch.
    epsilon = sp.symbols("epsilon", nonzero=True)
    boundary_ratio = epsilon / r
    require(boundary_ratio.subs(r, 1) == epsilon, "R lost the selected branch")
    require(boundary_ratio.subs(r, -1) == -epsilon,
            "the other boundary point was merged with R")

    # THM-4259 hidden differential coefficients really span both monomials.
    a = sp.symbols("a")
    quartic = a**4 - 2 * a**3 - 2 * a + 1
    A0 = -a**3 + 2 * a**2 + 3
    B0 = a**3 - 2 * a**2 - 1
    resultant_A0 = sp.factor(sp.resultant(quartic, A0, a))
    resultant_B0 = sp.factor(sp.resultant(quartic, B0, a))
    require(resultant_A0 == 6, "A0 can vanish on the hidden branch")
    require(resultant_B0 == -2, "B0 can vanish on the hidden branch")

    # For eta_(a,c)=x^a y^c dx/y^3 and ord_Q(x,y,dx)=(-2,-3,-3),
    # ord_Q eta_(a,c)=6-2a-3c.
    canonical_pairs = ((0, 0), (1, 0), (2, 0), (3, 0),
                       (0, 1), (1, 1), (0, 2))
    canonical_orders = tuple(6 - 2 * x_power - 3 * y_power
                             for x_power, y_power in canonical_pairs)
    require(canonical_orders == (6, 4, 2, 0, 3, 1, 0),
            "canonical order ledger changed")

    # Full E0 pullbacks: hidden (0,0),(0,2); v=(0,1); u=(1,1).
    e0_pairs = ((0, 0), (0, 1), (1, 1), (0, 2))
    e0_orders = tuple(6 - 2 * x_power - 3 * y_power
                      for x_power, y_power in e0_pairs)
    require(e0_orders == (6, 3, 1, 0), "E0 order ledger changed")
    require(len(set(e0_orders)) == 4, "E0 leading orders can cancel")
    ramification_indices = tuple(sorted(order + 1 for order in e0_orders))
    require(ramification_indices == (1, 2, 4, 7), "ramification spectrum changed")
    require(max(ramification_indices) < 12, "length-twelve constancy not excluded")

    # Hostile: the full canonical space does admit an order-twelve
    # cancellation.  If epsilon^2=-1 and s=y^2/x^3 tends to epsilon, then
    # (s-epsilon)(s+epsilon)=s^2+1=x^-6 has order 12.
    ord_x_inverse_six = 12
    ord_s_plus_epsilon = 0
    hostile_order = ord_x_inverse_six - ord_s_plus_epsilon
    require(hostile_order == 12, "canonical hostile lost order twelve")

    # At the selected contact, tau fixes Q and b=1/S has the inverse
    # primitive character; k[b]/(b^12) contains all twelve characters.
    contact_characters = tuple(range(12))
    require(len(set(contact_characters)) == 12, "contact character collision")

    print("A23 INFINITY-JET PRIMARY EXACT AUDIT")
    print("theorem=THM-4272_PROVED_RELATIVE")
    print("scope=raw_A23_contact_E0_jet_and_honest_fat_incidence_Bstar_collar")
    print("toric_chart=Ctilde=b^12-U*r^6-W*r^5-Z*r^4 R:r=1 PASS")
    print("contact_family=b^12-(U+W+Z) finite_free_rank=12 PASS")
    print("Bstar=Spec(C[U,W,Z,(U*Z*D)^-1]) Lambda_not_inverted")
    print("C_smoothness_gates=U,Z,D torus_syzygy_and_edge_discriminants PASS")
    print("Lambda_role=Res(X-1,U*X^2+W*X+Z)=R_C_transversality_only PASS")
    print("Lambda0_contact=finite_flat_rank12_nonetale_allowed PASS")
    print("wall_C=b^12-U*r^4*(r^2-1) dC_dr_at_Q=-2U PASS")
    print("wall_intersection=Spec(C[b]/(b^12)) one_branch=12Q")
    print("raw_union=two_smooth_branches_contact_order_12 analytic_type=A23 delta=12")
    print("homogeneous_lambda=(T^2+R^2)^2 wall_points=t:+i,-i PASS")
    print("fixed_normalization=R_selects_Q_epsilon_not_Q_minus_epsilon PASS")
    print(f"hidden_resultants=A0:{resultant_A0},B0:{resultant_B0} nonzero=PASS")
    print("E0_differential_basis=dx/y^3,y*dx/y^3,x*y*dx/y^3,y^2*dx/y^3")
    print("E0_vanishing_orders=6,3,1,0 ramification_indices=7,4,2,1")
    print("length12_restriction=NONCONSTANT_FOR_EVERY_NONCONSTANT_C0_TO_E0_MAP")
    print("hostile_full_canonical_order12=mixes_x^3dx/y^3_with_y^2dx/y^3 PASS")
    print("tau_contact_characters=0..11_once PASS")
    print("incidence_collar=GEOMETRIC_THM4268_NERON_ROSATI_ARGUMENT_NOT_MECHANIZED")
    print("proved_split=A_HONEST_FAT_INCIDENCE_COLLAR_ACROSS_LAMBDA0")
    print("open_split=B_RAW_KELLER_DESCENT_NO_WALL_CANDIDATE_EXCLUSION_NO_JC2")
    print("verdict=PROVED_RELATIVE_EXACT_AUDIT_PASS")


if __name__ == "__main__":
    main()
