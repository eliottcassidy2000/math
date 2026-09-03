#!/usr/bin/env python3
"""Independent exact audit of Long's Hamiltonian primitive descent obstruction.

This certificate copies the four core polynomials directly from the TeX source
of arXiv:2608.23777.  It does not import the THM-4397 implementation or the
exceptional-quartic scout.  Besides checking the two-form sign and the exact
three-point fibre, it identifies the correction's residue as a nontrivial
idempotent direction and records the precise boundary of a core-descended
compensator elimination.
"""

from __future__ import annotations

import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def exact(value: sp.Expr) -> sp.Expr:
    """Normalize an exact rational expression."""

    return sp.cancel(sp.together(sp.expand(value)))


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def wedge(f: sp.Expr, g: sp.Expr, variables: tuple[sp.Symbol, ...]) -> tuple[sp.Expr, ...]:
    """Coefficients of df wedge dg in lexicographic basis pairs."""

    answer: list[sp.Expr] = []
    for i in range(len(variables)):
        for j in range(i + 1, len(variables)):
            vi, vj = variables[i], variables[j]
            answer.append(exact(sp.diff(f, vi) * sp.diff(g, vj) - sp.diff(f, vj) * sp.diff(g, vi)))
    return tuple(answer)


def add_forms(left: tuple[sp.Expr, ...], right: tuple[sp.Expr, ...]) -> tuple[sp.Expr, ...]:
    return tuple(exact(a + b) for a, b in zip(left, right, strict=True))


def sub_forms(left: tuple[sp.Expr, ...], right: tuple[sp.Expr, ...]) -> tuple[sp.Expr, ...]:
    return tuple(exact(a - b) for a, b in zip(left, right, strict=True))


def reduce_x(poly: sp.Expr, x: sp.Symbol) -> sp.Expr:
    """Reduce on the exact reduced fibre Q[x]/(x^3-x)."""

    return sp.rem(sp.Poly(sp.expand(poly), x, domain=sp.QQ), sp.Poly(x**3 - x, x, domain=sp.QQ)).as_expr()


def main() -> None:
    x, y, beta = sp.symbols("x y beta")
    variables = (x, y, beta)
    q0 = y + x * beta / 3
    u = x * y

    # Long's formulas, copied directly from equations (R-def)--(H-def).
    R = 2 * x - 3 * x**2 * y - x**3 * beta
    S = y + 3 * x * (1 + u) ** 2 * beta + 3 * x * y**2 * (4 + 3 * u)
    T = -sp.Rational(1, 2) * (
        (1 + u) ** 3 * beta + y**2 * (1 + u) * (4 + 3 * u)
    )
    H = (
        y**4 * (18 * u**2 + 78 * u + 125) / 20
        + 3 * beta * y**2 * (u**3 + 5 * u**2 + 10 * u - 5) / 10
        - beta**2 * (9 * u + 2) / 6
        - x**2 * beta**3 / 6
    )

    check(exact(R - x * (2 - 3 * x * q0)) == 0, "R factorization in q0")
    check(exact((q0 - x * beta / 3) - y) == 0, "polynomial adapted-coordinate inverse")
    adapted_jacobian = sp.det(
        sp.Matrix([[sp.diff(f, v) for v in variables] for f in (x, q0, beta)])
    )
    check(adapted_jacobian == 1, "adapted-coordinate Jacobian")

    core_jacobian = sp.factor(
        sp.det(sp.Matrix([[sp.diff(f, v) for v in variables] for f in (R, T, S)]))
    )
    check(core_jacobian == 1, "oriented core Keller identity")

    # Paper convention: omega=dR wedge dD0+Theta and
    # Theta=dR wedge dH+dT wedge dS.  Hence Xi=Theta-dT wedge dS.
    theta = (
        exact(sp.Rational(9, 2) * q0 * (beta + 21 * q0**2)),
        exact(
            3 * q0**2
            + beta * (1 + 3 * x * q0) / 6
            + sp.Rational(3, 2) * x * q0 * (beta + 21 * q0**2)
        ),
        exact((1 + 3 * x * q0) / 2),
    )
    dR_dH = wedge(R, H, variables)
    dT_dS = wedge(T, S, variables)
    check(add_forms(dR_dH, dT_dS) == theta, "Long two-form coefficient identity")
    Xi = sub_forms(theta, dT_dS)
    check(Xi == dR_dH, "Xi equals dR wedge dH with paper sign")
    check(any(coefficient != 0 for coefficient in Xi), "primitive residual is nonzero")

    # Exact all-degree kernel certificate.  In coordinates (x,q0,beta), the
    # beta wedge coefficient forces h_beta=0 because R_q0=-3*x^2.  The
    # remaining Jacobian equation is the derivative at fixed R.
    X, Q, B, rho = sp.symbols("X Q B rho")
    hx, hq = sp.symbols("h_x h_q")
    R_adapted = 2 * X - 3 * X**2 * Q
    R_X = sp.diff(R_adapted, X)
    R_Q = sp.diff(R_adapted, Q)
    check(sp.diff(R_adapted, B) == 0, "R independent of beta")
    check(R_Q == -3 * X**2, "beta-forcing q0 derivative")
    Q_at_fixed_R = (2 * X - rho) / (3 * X**2)
    fixed_R_Q_derivative = sp.diff(Q_at_fixed_R, X)
    check(
        exact(fixed_R_Q_derivative - 2 * (rho - X) / (3 * X**3)) == 0,
        "fixed-R q0 derivative",
    )
    jacobian_equation = R_X * hq - R_Q * hx
    fixed_R_derivative = hx + fixed_R_Q_derivative * hq
    check(
        exact(
            (X * jacobian_equation).subs(Q, Q_at_fixed_R)
            - 3 * X**3 * fixed_R_derivative
        )
        == 0,
        "Jacobian equation is fixed-R derivative",
    )
    check(exact(Q_at_fixed_R.subs(rho, R_adapted) - Q) == 0, "function-field inverse")

    # The paper's exact core fibre and the primitive's residue on it.
    fibre_groebner = sp.groebner(
        [R, S, T - sp.Rational(1, 8)], beta, y, x, order="lex", domain=sp.QQ
    )
    expected_basis = (
        beta - sp.Rational(27, 4) * x**2 + sp.Rational(1, 4),
        y + sp.Rational(3, 2) * x,
        x**3 - x,
    )
    actual_basis = tuple(sp.expand(poly.as_expr()) for poly in fibre_groebner.polys)
    check(actual_basis == expected_basis, "reduced lexicographic core fibre basis")
    check(sp.gcd(x**3 - x, sp.diff(x**3 - x, x)) == 1, "core fibre is reduced")

    H_residue = sp.expand(fibre_groebner.reduce(H)[1])
    expected_residue = -sp.Rational(1, 48) - sp.Rational(1093, 192) * x**2
    check(exact(H_residue - expected_residue) == 0, "Hamiltonian fibre residue")

    fibre_points = (
        (sp.Integer(0), sp.Integer(0), -sp.Rational(1, 4)),
        (sp.Integer(1), -sp.Rational(3, 2), sp.Rational(13, 2)),
        (-sp.Integer(1), sp.Rational(3, 2), sp.Rational(13, 2)),
    )
    expected_values = (
        -sp.Rational(1, 48),
        -sp.Rational(1097, 192),
        -sp.Rational(1097, 192),
    )
    values: list[sp.Expr] = []
    original_points: list[tuple[sp.Expr, ...]] = []
    for point, expected in zip(fibre_points, expected_values, strict=True):
        substitution = dict(zip(variables, point, strict=True))
        check(exact(R.subs(substitution)) == 0, "fibre point R")
        check(exact(T.subs(substitution) - sp.Rational(1, 8)) == 0, "fibre point T")
        check(exact(S.subs(substitution)) == 0, "fibre point S")
        value = exact(H.subs(substitution))
        check(value == expected, "fibre point H")
        values.append(value)

        # Paper inverse (Psi_0) with D=0, hence D0=-H.  This also audits the
        # sign of the correction against the advertised original variables.
        q_value = exact(q0.subs(substitution))
        M_value = exact(beta.subs(substitution) + 9 * q_value**2)
        A_value = exact(1 - 3 * x.subs(substitution) * q_value)
        C_value = exact((1 + 3 * x.subs(substitution) * q_value) / 2)
        D0_value = -value
        p_value = exact(3 * q_value**2 * M_value + 2 * A_value * D0_value)
        z_value = exact(C_value * M_value - 3 * x.subs(substitution) ** 2 * D0_value)
        check(
            exact(C_value * p_value - 3 * q_value**2 * z_value - D0_value) == 0,
            "original D0 recovery",
        )
        check(
            exact(
                3 * x.subs(substitution) ** 2 * p_value
                + 2 * A_value * z_value
                - 9 * q_value**2
                - beta.subs(substitution)
            )
            == 0,
            "original beta recovery",
        )
        original_points.append((x.subs(substitution), q_value, p_value, z_value))

    expected_original_points = (
        (sp.Integer(0), sp.Integer(0), sp.Rational(1, 24), -sp.Rational(1, 8)),
        (sp.Integer(1), sp.Rational(2, 3), sp.Rational(247, 96), -sp.Rational(89, 64)),
        (-sp.Integer(1), -sp.Rational(2, 3), sp.Rational(247, 96), -sp.Rational(89, 64)),
    )
    check(tuple(original_points) == expected_original_points, "paper original three-point fibre")

    h0, h1 = expected_values[0], expected_values[1]
    separation = exact(h1 - h0)
    check(separation == -sp.Rational(1093, 192), "nonzero branch separation")
    check(values[1] == values[2], "hostile pair not separated by H")

    # The normalized residue is the idempotent x^2 in Q[x]/(x^3-x).
    idempotent = exact((H_residue - h0) / separation)
    check(idempotent == x**2, "normalized primitive residue")
    check(reduce_x(idempotent**2 - idempotent, x) == 0, "primitive residue idempotent")
    Z = sp.symbols("Z")
    check(
        reduce_x((H_residue - h0) * (H_residue - h1), x) == 0,
        "two-valued minimal-polynomial relation",
    )
    check(reduce_x(H_residue - h0, x) != 0, "not first scalar value")
    check(reduce_x(H_residue - h1, x) != 0, "not second scalar value")

    # On a constant/descended source graph D0=c, D=c+H therefore has the
    # same nonconstant residue.  Conversely D=d forces D0=d-H, which carries
    # the identical non-descended idempotent.  These are fibre certificates
    # for the general intersection lemma A[D0+H] cap N=A.
    c, d = sp.symbols("c d")
    D_on_descended_graph = c + H_residue
    D0_on_target_graph = d - H_residue
    check(
        reduce_x(D_on_descended_graph.subs(x, 1) - D_on_descended_graph.subs(x, 0), x)
        == separation,
        "descended D0 graph splits the core fibre",
    )
    check(
        exact(D0_on_target_graph.subs(x, 1) - D0_on_target_graph.subs(x, 0))
        == -separation,
        "constant D graph needs branch-sensitive D0",
    )

    # Cusp hostile: reduced geometric fibre constancy is not a sufficient
    # descent test.  Q[t^2,t^3] subset Q[t] is pointwise bijective on the
    # cusp, but t cannot occur because 1 is absent from <2,3>.
    semigroup_degrees = {2 * a + 3 * b for a in range(8) for b in range(8)}
    check(1 not in semigroup_degrees, "cusp normalization coordinate does not descend")
    check(all(n in semigroup_degrees for n in range(2, 15)), "cusp numerical semigroup check")

    print("LONG HAMILTONIAN PRIMITIVE DESCENT -- INDEPENDENT EXACT AUDIT")
    print("source_form_convention: omega=dR^dD0+Theta; Theta=dR^dH+dT^dS")
    print("Xi_convention: Xi=Theta-dT^dS=dR^dH")
    print(f"core_jacobian_det_RTS: {core_jacobian}")
    print("kernel_certificate:")
    print("  R=2*x-3*x^2*q0, R_q0=-3*x^2")
    print("  x*Jac(R,h)=3*x^3*(partial h/partial x at fixed R)")
    print("  ker(h -> dR^dh)=Q[R] (function-field constancy + q0=0 specialization)")
    print("additive_exact_wedge_quotient: (dR^dN)/(dR^dA) ~= N/A")
    print("primitive_torsor: H+Q[R]")
    print("core_fibre_groebner_basis:")
    for polynomial in expected_basis:
        print(f"  {polynomial}")
    print(f"H_fibre_residue: {H_residue}")
    print("H_fibre_values: " + ", ".join(str(value) for value in values))
    print(f"branch_separation: {separation}")
    print("normalized_residue: x^2; idempotent modulo x^3-x")
    minimal_polynomial = sp.factor((Z - h0) * (Z - h1))
    print(f"minimal_polynomial_on_fibre: {minimal_polynomial}")
    print("paper_original_points:")
    for point in original_points:
        print("  " + ", ".join(str(coordinate) for coordinate in point))
    print("destabilization_boundary:")
    print("  A-descended D0 graphs cannot make D descend; D=constant needs D0=d-H")
    print("  exact algebra lemma: A[D0+H] intersect N=A, so D0 is recoverable iff H in A")
    print("hostile_controls:")
    print("  H does not separate x=+1 from x=-1")
    print("  reduced-fibre constancy is not sufficient: Q[t^2,t^3] subset Q[t]")
    print("scope: fixed Long gauge/core; no obstruction to arbitrary source symplectic regauging")
    print(f"checks: {CHECKS}")


if __name__ == "__main__":
    main()
