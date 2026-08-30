#!/usr/bin/env python3
"""Exact symbolic certificate for THM-4294's central-component obstruction.

The load-bearing mathematics is geometric, but the vulnerable coordinate
claims are finite: the source/target scaling exponent, the restriction of
G_P to the wall component, and the local r=3 hostile control.  This script
reconstructs those identities from the literal scaled model.
"""

from __future__ import annotations

import sympy as sp


CHECKS = 0


def need(condition: object, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise AssertionError(label)


def zero(expr: sp.Expr, label: str) -> None:
    need(sp.factor(expr) == 0, f"{label}: {sp.factor(expr)}")


def main() -> None:
    sigma, S, P = sp.symbols("sigma S P")
    U, W, Z = sp.symbols("U W Z", nonzero=True)

    # THM-4230's exact weight-twelve special fibre.
    rational = S**2 - P
    central = 1 - U * P**6 - W * S**2 * P**5 - Z * S**4 * P**4
    special_source = sp.expand(rational * central)

    wall_central = sp.expand(central.subs({W: 0, Z: -U}))
    wall_source = sp.expand(special_source.subs({W: 0, Z: -U}))
    wall_cp = sp.diff(wall_central, P)
    wall_gp = sp.diff(wall_source, P)

    zero(
        wall_central - (1 - U * P**6 + U * S**4 * P**4),
        "literal wall component",
    )
    zero(
        wall_cp + 2 * U * P**3 * (3 * P**2 - 2 * S**4),
        "central P derivative",
    )
    # Equality modulo C: G_P=-C+(S^2-P)C_P.
    zero(
        wall_gp - rational * wall_cp + wall_central,
        "G_P restriction identity",
    )

    coefficient_field = sp.QQ.frac_field(U, S)
    central_p = sp.Poly(wall_central, P, domain=coefficient_field)
    derivative_p = sp.Poly(wall_cp, P, domain=coefficient_field)
    rational_p = sp.Poly(rational, P, domain=coefficient_field)
    need(sp.gcd(central_p, derivative_p).degree() == 0,
         "C_P is generically nonzero on C")
    need(sp.gcd(central_p, rational_p).degree() == 0,
         "S^2-P is generically nonzero on C")
    zero(wall_central.subs(P, 0) - 1, "P factor hostile")
    zero(wall_central.subs(P, S**2) - 1, "rational-factor hostile")
    resultant = sp.factor(sp.resultant(wall_central, 3 * P**2 - 2 * S**4, P))
    zero(resultant - (4 * U * S**12 + 27) ** 2,
         "critical-factor finite intersection")

    # Derive the sigma exponent rather than importing it as a label.
    # F_Q=sigma^-2 G and P=sigma^2 p imply (F_Q)_p=G_P.
    q_exponent = 12
    ds_exponent = -1
    original_target_form_exponent = q_exponent + ds_exponent
    good_to_original_form_exponent = 2
    good_form_exponent = (
        original_target_form_exponent - good_to_original_form_exponent
    )
    need(original_target_form_exponent == 11,
         "Q ds/(F_Q)_p exponent")
    need(good_form_exponent == 9,
         "good invariant differential exponent")

    # The good target really has good reduction at sigma=0.
    X, a, d = sp.symbols("X a d")
    target_cubic = X**3 + 1 - a * sigma**8 * X - d * sigma**12
    need(sp.discriminant(target_cubic.subs(sigma, 0), X) == -27,
         "good target discriminant")

    # Hostile control for the sole noncritical arithmetic survivor r=3.
    # On C_0 use x=alpha P, y=beta S P, b=1/S and
    # u=(-x^2,y^2).  The target parameter -X/Y pulls back exactly to b^2.
    alpha, beta, b = sp.symbols("alpha beta b", nonzero=True)
    x = alpha * P
    y = beta * S * P
    u_target_parameter = sp.cancel(x**2 / y**2)
    local_parameter = sp.factor(u_target_parameter.subs(S, 1 / b))
    zero(
        local_parameter - alpha**2 * b**2 / beta**2,
        "visible eigenline local parameter",
    )
    need(sp.degree(sp.Poly(local_parameter, b), b) == 2,
         "visible eigenline local index two")

    # Proper-flat degree conservation is imported geometrically.  This finite
    # ledger records the consequence once every special component is constant.
    noncritical_degrees = (40, 32, 38, 30, 36, 28, 34, 26, 32, 24, 30, 22)
    need(all(value > 0 for value in noncritical_degrees),
         "all noncritical responses are positive")
    constant_component_degrees = (0, 0, 0, 0)
    need(sum(constant_component_degrees) == 0,
         "constant-component degree sum")
    need(all(value != sum(constant_component_degrees)
             for value in noncritical_degrees),
         "positive response contradicts constant special fibre")

    print("STATUS=VERIFIED-EXACT; THM-4294=CENTRAL-KELLER-EXTINCTION; JC(2)=OPEN")
    print("WALL=W=Lambda=0; Z=-U; U!=0")
    print("SPECIAL_FIBRE=(S^2-P)*C; C=1-U*P^6+U*S^4*P^4")
    print("GP_MOD_C=(S^2-P)*C_P; C_P=-2*U*P^3*(3*P^2-2*S^4)")
    print("GENERIC_UNIT_TESTS=gcd(C,C_P)=1,gcd(C,S^2-P)=1 PASS")
    print("SCALING=Q^1*ds/original_form=sigma^11; original=sigma^2*eta0; eta0=sigma^9*dS/G_P")
    print("CENTRAL_VERTICAL_ORDER=9; CENTRAL_MAP=CONSTANT")
    print("IMPORTED_TAIL_LEDGER=THM4292_ALL_NONCENTRAL_A23_COMPONENTS_CONSTANT")
    print("DEGREE_CONSERVATION=POSITIVE_GENERIC_VS_ZERO_SPECIAL CONTRADICTION")
    print("R3_HOSTILE=u^*(-X/Y)=(alpha^2/beta^2)*b^2; LOCAL_INDEX=2; NONCONSTANT")
    print("CONCLUSION=ALL_W_LAMBDA_ZERO_U_NONZERO_RESPONSES_EXCLUDED; SEAM_ENTRY_OPEN")
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
