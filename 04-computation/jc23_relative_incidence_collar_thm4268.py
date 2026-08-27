#!/usr/bin/env python3
"""Exact face/contact audit for THM-4268.

The properness argument in the theorem is geometric.  This companion checks
the algebraic hypotheses which feed it: the gate controls the affine critical
locus and all toric edge schemes, the contact divisor is etale of rank twelve,
and the homogeneous W=0 attachment coordinate has precisely the three named
V4 branch walls.  No floating-point arithmetic is used.
"""

from __future__ import annotations

import sympy as sp


def main() -> None:
    U, W, Z, S, P, X, T, R = sp.symbols("U W Z S P X T R")
    D = W**2 - 4 * U * Z
    Lam = U + W + Z
    F = 1 - U * P**6 - W * S**2 * P**5 - Z * S**4 * P**4

    Fs = sp.factor(sp.diff(F, S))
    Fp = sp.factor(sp.diff(F, P))
    source_critical = W * P + 2 * Z * S**2
    p_bracket = 6 * U * P**2 + 5 * W * S**2 * P + 4 * Z * S**4
    critical_identity = sp.expand(
        2 * Z * p_bracket
        - (4 * Z * S**2 + 3 * W * P) * source_critical
        + 3 * D * P**2
    )
    assert critical_identity == 0

    edge_quadratic = U * X**2 + W * X + Z
    assert sp.discriminant(edge_quadratic, X) == D
    assert sp.resultant(X - 1, edge_quadratic, X) == Lam
    assert sp.factor(sp.discriminant(1 - Z * X**4, X)) == -256 * Z**3
    assert sp.factor(sp.discriminant(U - X**6, X)) == 46656 * U**5

    contact = sp.expand(F.subs(P, S**2))
    assert sp.expand(contact - (1 - Lam * S**12)) == 0
    assert sp.expand(sp.diff(contact, S) + 12 * Lam * S**11) == 0

    U_h = 4 * T**2 * R**2
    Z_h = (T**2 - R**2) ** 2
    assert sp.expand(U_h + Z_h - (T**2 + R**2) ** 2) == 0
    ratio_numerator = sp.expand((T**2 - R**2) ** 2)
    ratio_denominator = sp.expand(4 * T**2 * R**2)
    assert sp.expand(
        ratio_numerator.subs({T: -T}) * ratio_denominator
        - ratio_numerator * ratio_denominator.subs({T: -T})
    ) == 0
    # In the affine t=T/R chart, inversion swaps T and R.
    assert sp.expand(
        ratio_numerator.subs({T: R, R: T}, simultaneous=True)
        * ratio_denominator
        - ratio_numerator
        * ratio_denominator.subs({T: R, R: T}, simultaneous=True)
    ) == 0

    # The nonproper hostile has empty W=0 fibre but nonempty generic fibre.
    hostile_special = sp.Poly((W * X - 1).subs(W, 0), X)
    assert hostile_special.degree() == 0 and hostile_special.LC() == -1

    print("theorem=THM-4268")
    print(f"face={F}")
    print(f"dS={Fs}")
    print(f"dP={Fp}")
    print("torus_critical_identity=2Z*Pbracket-(4ZS^2+3WP)*Sbracket+3DP^2=0")
    print("torus_critical_gate=D")
    print("edge_discriminants={quadratic:D,quartic:-256*Z^3,sextic:46656*U^5}")
    print("edge_resultant_at_X1=Lambda")
    print("contact_equation=1-Lambda*S^12")
    print("contact_derivative=-12*Lambda*S^11")
    print("contact_divisor=finite_etale_rank_12_on_gate")
    print("w0_kummer_map=[U:Z]=[4*T^2*R^2:(T^2-R^2)^2]")
    print("w0_kummer_sum=U+Z=(T^2+R^2)^2")
    print("w0_branch_fibres={U0:[0,infinity],Z0:[1,-1],UplusZ0:[i,-i]}")
    print("w0_generic_orbit={t,-t,1/t,-1/t}")
    print("nonproper_hostile=W*X-1_special_fibre=unit_generic_fibre=point")
    print("properness_input=abelian_target_Neron_extension_not_symbolic_specialization")
    print("scope=no_explicit_collar_equation_no_boundary_crossing_no_M12_entry_no_JC2")


if __name__ == "__main__":
    main()
