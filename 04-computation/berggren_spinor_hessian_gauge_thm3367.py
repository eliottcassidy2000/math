#!/usr/bin/env python3
"""Exact companion for THM-3367.

The universal statements in the theorem are proved algebraically there.  This
companion checks the load-bearing polynomial identities, generic rank-two
ideal, positive controls, and hostile controls using exact arithmetic only.
"""

from __future__ import annotations

import sympy as sp


def require(label: str, condition: object, checks: list[str]) -> None:
    """Record one check, with behavior unchanged by ``python -O``."""

    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"FAIL: {label}: {condition!r}")
    checks.append(label)


def matrix_is_zero(matrix: sp.MatrixBase) -> bool:
    return all(sp.simplify(entry) == 0 for entry in matrix)


def vector_substitute(
    vector: sp.MatrixBase, substitutions: dict[sp.Symbol, sp.Expr]
) -> sp.Matrix:
    return sp.Matrix(
        [sp.expand(entry.subs(substitutions, simultaneous=True)) for entry in vector]
    )


def main() -> None:
    checks: list[str] = []

    # The Berggren spinor pencil and its symmetric/Clifford target gauge.
    a, b, c = sp.symbols("a b c")
    L = sp.Matrix([[0, 1], [-1, 2]])
    M = sp.Matrix([[0, 1], [1, 2]])
    R = sp.Matrix([[1, 0], [2, 1]])
    D = sp.diag(-1, 1)
    J = sp.Matrix([[0, 1], [1, 0]])
    pencil = a * L + b * M + c * R
    symmetric_pencil = sp.Matrix([[b - a, c], [c, a + b]])

    require("L=M*D", L == M * D, checks)
    require("R=M*J", R == M * J, checks)
    require("M^{-1} pencil is symmetric", M.inv() * pencil == symmetric_pencil, checks)
    require("Clifford involutions", D**2 == sp.eye(2) and J**2 == sp.eye(2), checks)
    require("Clifford anticommutation", D * J + J * D == sp.zeros(2), checks)
    require(
        "determinant sign",
        sp.factor(pencil.det()) == a**2 - b**2 + c**2
        and sp.factor(symmetric_pencil.det()) == -(a**2 - b**2 + c**2),
        checks,
    )

    # A generic exact polynomial potential verifies both directions of the
    # displayed construction without assuming the Keller equation.
    s, t = sp.symbols("s t")
    H = (
        3 * s**4
        - 2 * s**3 * t
        + 5 * s**2 * t**2
        + 7 * s * t**3
        - 11 * t**4
        + 13 * s**2
        - 17 * s * t
        + 19 * t**2
        + 23 * s
        - 29 * t
    )
    Hss = sp.diff(H, s, 2)
    Hst = sp.diff(H, s, t)
    Htt = sp.diff(H, t, 2)
    weights = {
        a: sp.expand((Htt - Hss) / 2),
        b: sp.expand((Htt + Hss) / 2),
        c: Hst,
    }
    pulled_pencil = pencil.subs(weights)
    phi = sp.Matrix([sp.diff(H, t), sp.diff(H, s) + 2 * sp.diff(H, t)])
    require(
        "generic potential gives the pulled-back Jacobian",
        matrix_is_zero(pulled_pencil - phi.jacobian([s, t])),
        checks,
    )
    require(
        "generic determinant equals minus Hessian determinant",
        sp.simplify(pulled_pencil.det() + sp.det(sp.hessian(H, (s, t)))) == 0,
        checks,
    )

    # The affine-rank no-go has a short dimension proof.  Algebraically, the
    # three squared 2x2 minors reduce to zero modulo the generic equations
    # saying that two columns are isotropic and mutually orthogonal.
    u1, u2, u3, v1, v2, v3 = sp.symbols("u1 u2 u3 v1 v2 v3")
    q_u = u1**2 - u2**2 + u3**2
    q_v = v1**2 - v2**2 + v3**2
    polar_uv = u1 * v1 - u2 * v2 + u3 * v3
    rank_ideal = sp.groebner(
        [q_u, q_v, polar_uv],
        u1,
        u2,
        u3,
        v1,
        v2,
        v3,
        order="grevlex",
    )
    minors = [
        u1 * v2 - u2 * v1,
        u1 * v3 - u3 * v1,
        u2 * v3 - u3 * v2,
    ]
    require(
        "generic affine rank-two ideal kills every squared minor",
        all(rank_ideal.reduce(minor**2)[1] == 0 for minor in minors),
        checks,
    )

    # Rational affine-line Keller control and its inverse.
    line_H = s * t + s**6 / sp.Integer(30)
    line_phi = sp.Matrix(
        [sp.diff(line_H, t), sp.diff(line_H, s) + 2 * sp.diff(line_H, t)]
    )
    require(
        "rational affine-line control has determinant one",
        sp.factor(line_phi.jacobian([s, t]).det()) == 1,
        checks,
    )
    U, V = sp.symbols("U V")
    line_inverse = sp.Matrix([U, V - U**5 / 5 - 2 * U])
    require(
        "rational affine-line inverse: forward after inverse",
        matrix_is_zero(
            vector_substitute(line_phi, {s: line_inverse[0], t: line_inverse[1]})
            - sp.Matrix([U, V])
        ),
        checks,
    )
    require(
        "rational affine-line inverse: inverse after forward",
        matrix_is_zero(
            vector_substitute(line_inverse, {U: line_phi[0], V: line_phi[1]})
            - sp.Matrix([s, t])
        ),
        checks,
    )

    # A second exact control over Q(i) exercises a nontrivial isotropic
    # direction g=(1,i) for C=I and the inverse after the M target gauge.
    ell = s + sp.I * t
    complex_H = (s**2 + t**2) / 2 + ell**4 / 12
    complex_gradient = sp.Matrix([sp.diff(complex_H, s), sp.diff(complex_H, t)])
    complex_phi = M * complex_gradient
    require(
        "Q(i) affine-line control has determinant minus one",
        sp.simplify(complex_phi.jacobian([s, t]).det()) == -1,
        checks,
    )
    target1, target2 = sp.symbols("target1 target2")
    target = sp.Matrix([target1, target2])
    ungauged = M.inv() * target
    sigma = sp.expand(ungauged[0] + sp.I * ungauged[1])
    complex_inverse = sp.Matrix(
        [ungauged[0] - sigma**3 / 3, ungauged[1] - sp.I * sigma**3 / 3]
    )
    require(
        "Q(i) affine-line inverse: forward after inverse",
        matrix_is_zero(
            vector_substitute(
                complex_phi, {s: complex_inverse[0], t: complex_inverse[1]}
            )
            - target
        ),
        checks,
    )
    require(
        "Q(i) affine-line inverse: inverse after forward",
        matrix_is_zero(
            vector_substitute(
                complex_inverse,
                {target1: complex_phi[0], target2: complex_phi[1]},
            )
            - sp.Matrix([s, t])
        ),
        checks,
    )

    # Constant and zero-determinant boundary controls.
    constant_H = (s**2 + t**2) / 2
    constant_phi = M * sp.Matrix([sp.diff(constant_H, s), sp.diff(constant_H, t)])
    require(
        "constant coefficient boundary is affine Keller",
        constant_phi.jacobian([s, t]).det() == -1,
        checks,
    )
    degenerate_H = s**2 / 2 + s**4
    degenerate_phi = M * sp.Matrix(
        [sp.diff(degenerate_H, s), sp.diff(degenerate_H, t)]
    )
    require(
        "zero-determinant affine-line boundary is not Keller",
        degenerate_phi.jacobian([s, t]).det() == 0,
        checks,
    )

    # Hostile: a nonlinear two-parameter map into q=1 which fails both
    # Hessian curl equations.  Constant determinant alone is not integrability.
    p = s
    r = -t * (2 + s * t)
    hostile_c = 1 + s * t
    hostile_a = sp.expand((p + r) / 2)
    hostile_b = sp.expand((r - p) / 2)
    hostile_q = sp.expand(hostile_a**2 - hostile_b**2 + hostile_c**2)
    hostile_S = sp.Matrix(
        [[hostile_b - hostile_a, hostile_c], [hostile_c, hostile_a + hostile_b]]
    )
    hostile_curls = (
        sp.expand(sp.diff(hostile_S[0, 0], t) - sp.diff(hostile_S[0, 1], s)),
        sp.expand(sp.diff(hostile_S[0, 1], t) - sp.diff(hostile_S[1, 1], s)),
    )
    require("constant-q hostile lies on q=1", hostile_q == 1, checks)
    require(
        "constant-q hostile fails integrability",
        hostile_curls == (-t, s + t**2),
        checks,
    )

    # Exact dehomogenized binary-Hessian formula, checked with generic
    # coefficients through degree eight, plus pure-power and near controls.
    z = sp.symbols("z")
    for degree in range(2, 9):
        coefficients = sp.symbols(f"h{degree}_0:{degree + 1}")
        polynomial = sum(
            coefficients[index] * z**index for index in range(degree + 1)
        )
        binary_form = sp.expand(t**degree * polynomial.subs(z, s / t))
        formula = sp.expand(
            (degree - 1)
            * t ** (2 * degree - 4)
            * (
                degree * polynomial * sp.diff(polynomial, z, 2)
                - (degree - 1) * sp.diff(polynomial, z) ** 2
            ).subs(z, s / t)
        )
        require(
            f"binary Hessian formula degree {degree}",
            sp.expand(sp.det(sp.hessian(binary_form, (s, t))) - formula) == 0,
            checks,
        )
        require(
            f"pure-power/near-hostile degree {degree}",
            sp.det(sp.hessian((s - 2 * t) ** degree, (s, t))) == 0
            and sp.det(sp.hessian(s ** (degree - 1) * t, (s, t))) != 0,
            checks,
        )

    # The FC-facing cubic determinant is a different, non-functorial linear
    # mixing: Sym^2 agrees branchwise but not after adding the branches.
    T_L = sp.Matrix([[1, -2, 2], [2, -1, 2], [2, -2, 3]])
    T_M = sp.Matrix([[1, 2, 2], [2, 1, 2], [2, 2, 3]])
    T_R = sp.Matrix([[-1, 2, 2], [-2, 1, 2], [-2, 2, 3]])
    triple_pencil = a * T_L + b * T_M + c * T_R
    require(
        "triple-transfer determinant factorization",
        sp.factor(triple_pencil.det())
        == (a - b - c) * (a + b - c) * (a + b + c),
        checks,
    )
    require(
        "2x2/3x3 determinant-mixing hostile at (1,1,1)",
        (L + M + R).det() == 1 and (T_L + T_M + T_R).det() == -3,
        checks,
    )

    print("THM-3367 exact companion")
    print("universe=C[s,t]; symbolic rings over Z, Q, and Q(i)")
    print("binary_hessian_formula_degrees=2..8")
    print("affine_rank_audit=generic squared-minor ideal")
    print(f"hostile_hessian_curls={hostile_curls}")
    print("determinant_mixing_hostile=(spinor:1,triple:-3)")
    print(f"checks={len(checks)}")
    print("PASS")


if __name__ == "__main__":
    main()
