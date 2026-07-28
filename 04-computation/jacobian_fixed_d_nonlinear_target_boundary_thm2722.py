#!/usr/bin/env python3
"""Exact symbolic companion for THM-2722.

The theorem factors every graph target in C[A,d] through the planar pair
(A,d).  The all-degree unit and integration arguments are mathematical; this
companion independently checks the formal chain rule, the graph-coordinate
factor, the complete fixed-d normal form and inverse, a nontrivial planar
Keller composition, and sharp non-Keller/zero-response controls.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def jacobian(first: sp.Expr, second: sp.Expr, left: sp.Symbol, right: sp.Symbol) -> sp.Expr:
    return sp.factor(
        sp.diff(first, left) * sp.diff(second, right)
        - sp.diff(first, right) * sp.diff(second, left)
    )


def main() -> None:
    # Universal two-by-two chain rule with no specialization of P or Q.
    PA, Pd, QA, Qd = sp.symbols("P_A P_d Q_A Q_d")
    Ax, At, dx, dt = sp.symbols("A_x A_t d_x d_t")
    composed = sp.det(
        sp.Matrix(
            [
                [PA * Ax + Pd * dx, PA * At + Pd * dt],
                [QA * Ax + Qd * dx, QA * At + Qd * dt],
            ]
        )
    )
    target_det = PA * Qd - Pd * QA
    base_det = Ax * dt - At * dx
    require(
        sp.expand(composed - target_det * base_det) == 0,
        "formal target/base Jacobian factorization changed",
    )
    graph_factor = sp.factor(composed.subs({Ax: 0, At: -2}))
    require(
        graph_factor == 2 * dx * target_det,
        "A=-2t graph factor changed",
    )

    x, t = sp.symbols("x t")
    A, d = sp.symbols("A d")
    a, alpha = sp.symbols("a alpha", nonzero=True)
    r0, r1, r2 = sp.symbols("r0 r1 r2")
    k0, k1, k2 = sp.symbols("k0 k1 k2")
    # Quadratic controls keep both the graph and target genuinely nonlinear
    # while making the two-sided inverse replay inexpensive.  The theorem's
    # arbitrary-degree statement is the unit/integration argument, not this
    # finite symbolic sample.
    R = r0 + r1 * t + r2 * t**2
    K = k0 + k1 * d + k2 * d**2
    F = a * x + R
    A_graph = -2 * t

    # The coordinate pair (A,d)=(-2t,ax+R(t)) and its inverse.
    A0, d0 = sp.symbols("A0 d0")
    t_from_base = -A0 / 2
    x_from_base = sp.factor((d0 - R.subs(t, t_from_base)) / a)
    require(
        sp.factor(A_graph.subs(t, t_from_base) - A0) == 0,
        "base inverse lost A",
    )
    require(
        sp.factor(F.subs({x: x_from_base, t: t_from_base}, simultaneous=True) - d0)
        == 0,
        "base inverse lost d",
    )
    require(
        jacobian(A_graph, F, x, t) == 2 * a,
        "base coordinate Jacobian changed",
    )

    # Complete nonlinear fixed-d normal form P=alpha*A+K(d).
    U = sp.expand(alpha * A_graph + K.subs(d, F))
    V = F
    require(
        jacobian(U, V, x, t) == 2 * a * alpha,
        "fixed-d nonlinear target Jacobian changed",
    )

    u, v = sp.symbols("u v")
    t_back = sp.factor((K.subs(d, v) - u) / (2 * alpha))
    x_back = sp.factor((v - R.subs(t, t_back)) / a)
    require(
        sp.factor(t_back.subs({u: U, v: V}, simultaneous=True) - t) == 0,
        "fixed-d inverse failed on t",
    )
    require(
        sp.factor(x_back.subs({u: U, v: V}, simultaneous=True) - x) == 0,
        "fixed-d inverse failed on x",
    )
    require(
        sp.factor(V.subs({x: x_back, t: t_back}, simultaneous=True) - v) == 0,
        "fixed-d forward map failed after inverse on v",
    )
    require(
        sp.factor(U.subs({x: x_back, t: t_back}, simultaneous=True) - u) == 0,
        "fixed-d forward map failed after inverse on u",
    )

    # A nontrivial planar Keller pair in (A,d), composed with the graph
    # coordinate automorphism.  This is a positive control for the exact
    # equivalence boundary, not evidence for arbitrary JC(2).
    P_henon = d
    Q_henon = -A + d**2
    require(
        jacobian(P_henon, Q_henon, A, d) == 1,
        "planar Hénon control determinant changed",
    )
    U_henon = sp.expand(P_henon.subs({A: A_graph, d: F}))
    V_henon = sp.expand(Q_henon.subs({A: A_graph, d: F}))
    require(
        jacobian(U_henon, V_henon, x, t) == 2 * a,
        "composed planar Keller determinant changed",
    )
    t_henon_back = sp.factor((v - u**2) / 2)
    x_henon_back = sp.factor((u - R.subs(t, t_henon_back)) / a)
    require(
        sp.factor(U_henon.subs({x: x_henon_back, t: t_henon_back}, simultaneous=True) - u)
        == 0,
        "Hénon composition inverse failed on first target",
    )
    require(
        sp.factor(V_henon.subs({x: x_henon_back, t: t_henon_back}, simultaneous=True) - v)
        == 0,
        "Hénon composition inverse failed on second target",
    )

    # Sharp hostiles: nonlinear A-dependence has variable response, and each
    # factor of the product can independently kill the zero-response branch.
    P_variable = A**2 + d
    Q_fixed = d
    variable_graph = jacobian(
        P_variable.subs({A: A_graph, d: F}),
        Q_fixed.subs(d, F),
        x,
        t,
    )
    require(variable_graph == -8 * a * t, "variable-Jacobian hostile changed")
    require(
        jacobian(A_graph, R, x, t) == 0,
        "F_x=0 zero-response boundary changed",
    )
    require(
        jacobian((d**2).subs(d, F), F, x, t) == 0,
        "P_A=0 zero-response boundary changed",
    )

    print("THM-2722 FIXED-d NONLINEAR TARGET CLASSIFICATION")
    print("graph_coordinates=A:-2t,d:F(x,t)")
    print("chain_rule=Jac(P(A,d),Q(A,d))=2*D_(P,Q)(-2t,F)*F_x")
    print("nonzero_constant_forces=F:a*x+R(t),D_(P,Q):delta_constant")
    print("base_map=(-2t,a*x+R(t)):polynomial_automorphism")
    print("fixed_d_iff=P:alpha*A+K(d),a*alpha_nonzero")
    print("fixed_d_inverse=t:(K(V)-U)/(2alpha),x:(V-R(t))/a")
    print("full_C[A,d]^2=ordinary_planar_Keller_pair_after_base_change")
    print("counterexample_in_full_subalgebra_iff=JC2_false")
    print("zero_response_boundaries=F_x_zero_or_target_determinant_zero")
    print("scope=FIXED_d_FACE_CLOSED_FULL_SUBALGEBRA_IS_JC2_NOT_A_JC2_PROOF")
    print("PASS")


if __name__ == "__main__":
    main()
