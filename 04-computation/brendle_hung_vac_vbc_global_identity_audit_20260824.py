#!/usr/bin/env python3
"""Exact global audit of the Brendle--Hung ``V_ac`` and ``V_bc`` identities.

Source under audit
------------------
Simon Brendle and Pei-Ken Hung, *A metric on S^2 x S^2 with positive
sectional curvature*, arXiv:2608.19068v1, Proposition 5.2.  The saved
companion notebook constructs the intended ``Vbc`` expression but then
overwrites it with ``FullSimplify[Vac,...]`` before asserting ``Vbc==0``.

This script reuses the repository's independent reconstruction of the paper's
moving-frame metric, connection, curvature, linearized-curvature operator,
and quadratic-curvature operator, then independently performs the global
normal-Hessian and mixed-coefficient reductions.  Its independence is from
the saved Mathematica state, not a second implementation of the reconstructed
geometry.  It makes the exact substitution ``x=tan(theta)>0``.  All
denominators are nonzero on the generic locus ``0<theta<pi/2``.  The three
summands of both ``V_ac`` and ``V_bc`` are reduced separately and are each
nonzero as formal expressions; their sums vanish identically before imposing
any ``SO(3)`` relations.

The script also reconciles two sign conventions for the quarter-turn
symmetry ``F(p1,p2)=(R p1,-R p2)``.  Pullback by ``F`` has signs ``(+,+,+)``
on ``(h_a,h_c,h_ac)->(h_b,h_c,h_bc)``, whereas pushforward by ``F`` (the
pullback by ``F^{-1}``) has signs ``(-,+,-)``.  Thus the source ledger's
pushforward formula and the pullback-plus computation are equivalent, not a
discrepancy.

Status and scope
----------------
VERIFIED-EXACT GLOBAL SYMBOLIC IDENTITY relative to the reconstructed source
formulas.  This closes the overwritten ``V_bc=0`` calculation, not the other
quadratic identities, cubic integral, global smoothness, or headline metric.

Reproduce from the repository root with

    python3 04-computation/brendle_hung_vac_vbc_global_identity_audit_20260824.py
    python3 -O 04-computation/brendle_hung_vac_vbc_global_identity_audit_20260824.py
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import sympy as sp

import brendle_hung_lemma54_independent_audit_20260824 as bh
import brendle_hung_vbc_exact_point_audit_20260824 as point_audit


X = sp.symbols("x", positive=True)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def exact_zero(expression: sp.Expr) -> bool:
    return sp.factor(sp.cancel(expression)) == 0


def rationalize_theta(expression: sp.Expr) -> sp.Expr:
    """Rewrite exactly using ``x=tan(theta)>0`` on the generic interval."""

    expression = sp.expand_trig(expression)
    expression = expression.subs(
        {
            sp.sin(bh.THETA): X / sp.sqrt(1 + X**2),
            sp.cos(bh.THETA): 1 / sp.sqrt(1 + X**2),
            sp.tan(bh.THETA): X,
        },
        simultaneous=True,
    )
    return sp.factor(sp.cancel(sp.powsimp(expression, force=True)))


def tensor_family() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    ex = sp.Matrix([[1, 0, 0]])
    ey = sp.Matrix([[0, 1, 0]])
    ez = sp.Matrix([[0, 0, 1]])
    k_ex, k_ey, k_ez = (bh.k_covector(vector) for vector in (ex, ey, ez))
    j_ey = bh.j_covector(ey)

    h_b = (
        sp.Rational(6, 5)
        * (2 + sp.cos(2 * bh.THETA))
        * (bh.K_MINUS.T * k_ey + k_ey.T * bh.K_MINUS)
        / 2
        - 2
        * (1 + sp.cos(2 * bh.THETA))
        * (2 - sp.cos(2 * bh.THETA))
        * (bh.K_PLUS.T * j_ey + j_ey.T * bh.K_PLUS)
        / 2
    ).applyfunc(sp.trigsimp)
    h_c = (k_ez.T * k_ez).applyfunc(sp.trigsimp)
    h_ac = (
        sp.Rational(1, 10)
        * sp.cos(bh.THETA)
        * bh.Q1Z
        * (-5 - 14 * sp.cos(2 * bh.THETA) + 3 * sp.cos(4 * bh.THETA))
        * (k_ex.T * k_ez + k_ez.T * k_ex)
        / 2
    ).applyfunc(sp.trigsimp)
    h_bc = (
        sp.Rational(1, 10)
        * sp.sin(bh.THETA)
        * bh.Q2Z
        * (-5 + 14 * sp.cos(2 * bh.THETA) + 3 * sp.cos(4 * bh.THETA))
        * (k_ey.T * k_ez + k_ez.T * k_ey)
        / 2
    ).applyfunc(sp.trigsimp)
    return h_b, h_c, h_ac, h_bc


def matrix_equal(first: sp.Matrix, second: sp.Matrix) -> bool:
    for entry in first - second:
        normalized = sp.trigsimp(sp.expand_trig(entry))
        if normalized != 0:
            return False
    return True


def symmetry_audit(
    h_b: sp.Matrix, h_c: sp.Matrix, h_ac: sp.Matrix, h_bc: sp.Matrix
) -> None:
    """Check pullback-plus and pushforward-minus conventions entrywise."""

    # R(x,y,z)=(-y,x,z).  At F(p),
    # q1'=-R q2, q2'=-R q1, q3'=-R q3, theta'=pi/2-theta.
    forward_substitution = {
        bh.THETA: sp.pi / 2 - bh.THETA,
        bh.Q1X: bh.Q2Y,
        bh.Q1Y: -bh.Q2X,
        bh.Q1Z: -bh.Q2Z,
        bh.Q2X: bh.Q1Y,
        bh.Q2Y: -bh.Q1X,
        bh.Q2Z: -bh.Q1Z,
        bh.Q3X: bh.Q3Y,
        bh.Q3Y: -bh.Q3X,
        bh.Q3Z: -bh.Q3Z,
    }
    # At F^{-1}(p), replace R by R^{-1}.
    inverse_substitution = {
        bh.THETA: sp.pi / 2 - bh.THETA,
        bh.Q1X: -bh.Q2Y,
        bh.Q1Y: bh.Q2X,
        bh.Q1Z: -bh.Q2Z,
        bh.Q2X: -bh.Q1Y,
        bh.Q2Y: bh.Q1X,
        bh.Q2Z: -bh.Q1Z,
        bh.Q3X: -bh.Q3Y,
        bh.Q3Y: bh.Q3X,
        bh.Q3Z: -bh.Q3Z,
    }
    # dF(m1,m2,m3,m4)=(-m2',-m1',-m3',-m4').
    frame_matrix = sp.Matrix(
        [[0, -1, 0, 0], [-1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]
    )

    def pullback(tensor: sp.Matrix, substitution: dict[sp.Symbol, sp.Expr]):
        return frame_matrix.T * tensor.xreplace(substitution) * frame_matrix

    require(
        matrix_equal(pullback(bh.METRIC, forward_substitution), bh.METRIC),
        "F failed the background-isometry check",
    )
    require(
        matrix_equal(pullback(bh.H_A, forward_substitution), h_b),
        "F^* h_a != h_b",
    )
    require(
        matrix_equal(pullback(h_c, forward_substitution), h_c),
        "F^* h_c != h_c",
    )
    require(
        matrix_equal(pullback(h_ac, forward_substitution), h_bc),
        "F^* h_ac != h_bc",
    )
    require(
        matrix_equal(pullback(bh.H_A, inverse_substitution), -h_b),
        "F_* h_a != -h_b",
    )
    require(
        matrix_equal(pullback(h_c, inverse_substitution), h_c),
        "F_* h_c != h_c",
    )
    require(
        matrix_equal(pullback(h_ac, inverse_substitution), -h_bc),
        "F_* h_ac != -h_bc",
    )


def audit() -> list[str]:
    h_b, h_c, h_ac, h_bc = tensor_family()

    t_a, l_a = point_audit.tensor_operators(bh.H_A)
    t_b, l_b = point_audit.tensor_operators(h_b)
    t_c, l_c = point_audit.tensor_operators(h_c)
    _, l_ac = point_audit.tensor_operators(h_ac)
    _, l_bc = point_audit.tensor_operators(h_bc)

    r_a = sp.Matrix(
        [rationalize_theta(2 * l_a(*index)) for index in bh.P1_INDICES]
    )
    r_b = sp.Matrix(
        [rationalize_theta(2 * l_b(*index)) for index in bh.P1_INDICES]
    )
    r_c = sp.Matrix(
        [rationalize_theta(2 * l_c(*index)) for index in bh.P1_INDICES]
    )

    hessian = point_audit.p2_background().applyfunc(rationalize_theta)
    hessian_determinant = sp.factor(sp.cancel(hessian.det()))
    expected_determinant = (
        256
        * X**2
        * (X**2 + 1) ** 6
        * (3 * X**2 + 5)
        * (5 * X**2 + 3)
        / (81 * (X**2 + 3) ** 5 * (3 * X**2 + 1) ** 5)
    )
    require(
        exact_zero(hessian_determinant - expected_determinant),
        "normal-Hessian determinant changed",
    )

    z_c = (-hessian.inv() * r_c).applyfunc(
        lambda entry: sp.factor(sp.cancel(entry))
    )
    expected_z_c = sp.Matrix(
        [
            -bh.Q1Z * bh.Q3Z * X / (3 * sp.sqrt(X**2 + 1)),
            0,
            -bh.Q2Z * bh.Q3Z / (3 * sp.sqrt(X**2 + 1)),
            0,
        ]
    )
    require(
        all(exact_zero(entry) for entry in z_c - expected_z_c),
        f"z_c mismatch: {tuple(z_c)}",
    )

    linear_ac = rationalize_theta(l_ac(2, 3, 2, 3))
    quadratic_ac = rationalize_theta(point_audit.q_p0(t_a, t_c))
    minimizer_ac = sp.factor(sp.cancel(sp.Rational(1, 2) * r_a.dot(z_c)))
    linear_bc = rationalize_theta(l_bc(2, 3, 2, 3))
    quadratic_bc = rationalize_theta(point_audit.q_p0(t_b, t_c))
    minimizer_bc = sp.factor(sp.cancel(sp.Rational(1, 2) * r_b.dot(z_c)))

    denominator = (X**2 + 1) ** sp.Rational(5, 2)
    factor_ac = bh.Q1Z * bh.Q3X * bh.Q3Z
    factor_bc = bh.Q2Z * bh.Q3Y * bh.Q3Z
    expected_terms = {
        "linear_ac": factor_ac
        * (X**2 - 1)
        * (27 * X**2 + 2)
        / (15 * denominator),
        "quadratic_ac": -factor_ac
        * (203 * X**6 - 150 * X**4 - 159 * X**2 - 6)
        / (45 * denominator * (3 * X**2 + 1)),
        "minimizer_ac": -2
        * factor_ac
        * X**2
        * (20 * X**4 + 3 * X**2 + 33)
        / (45 * denominator * (3 * X**2 + 1)),
        "linear_bc": -factor_bc
        * X
        * (X**2 - 1)
        * (2 * X**2 + 27)
        / (15 * denominator),
        "quadratic_bc": factor_bc
        * X
        * (6 * X**6 + 159 * X**4 + 150 * X**2 - 203)
        / (45 * denominator * (X**2 + 3)),
        "minimizer_bc": -2
        * factor_bc
        * X
        * (33 * X**4 + 3 * X**2 + 20)
        / (45 * denominator * (X**2 + 3)),
    }
    actual_terms = {
        "linear_ac": linear_ac,
        "quadratic_ac": quadratic_ac,
        "minimizer_ac": minimizer_ac,
        "linear_bc": linear_bc,
        "quadratic_bc": quadratic_bc,
        "minimizer_bc": minimizer_bc,
    }
    for name, expected in expected_terms.items():
        require(
            exact_zero(actual_terms[name] - expected),
            f"{name} mismatch: {actual_terms[name]} != {expected}",
        )

    v_ac = sp.factor(sp.cancel(linear_ac + quadratic_ac + minimizer_ac))
    v_bc = sp.factor(sp.cancel(linear_bc + quadratic_bc + minimizer_bc))
    require(v_ac == 0, f"V_ac did not vanish: {v_ac}")
    require(v_bc == 0, f"V_bc did not vanish: {v_bc}")

    hostile_drop_bc = sp.factor(sp.cancel(linear_bc + quadratic_bc))
    hostile_flip_bc = sp.factor(
        sp.cancel(linear_bc + quadratic_bc - minimizer_bc)
    )
    require(hostile_drop_bc != 0, "dropping the minimizer term was invisible")
    require(hostile_flip_bc != 0, "flipping the minimizer sign was invisible")

    symmetry_audit(h_b, h_c, h_ac, h_bc)

    script_path = Path(__file__)
    point_path = Path(point_audit.__file__)
    geometry_path = Path(bh.__file__)
    rows = [
        "status=VERIFIED-EXACT",
        "source=arXiv:2608.19068v1 Proposition 5.2 equations (14)-(15)",
        "universe=0<theta<pi/2; x=tan(theta)>0; q-components formal",
        "method=independent moving-frame reconstruction followed by exact rational reduction",
        f"normal_hessian_determinant={hessian_determinant}",
        "normal_hessian_conclusion=determinant positive and matrix invertible for every x>0; positive-definiteness not asserted",
        f"z_c={tuple(z_c)}",
        f"V_ac_terms={(linear_ac, quadratic_ac, minimizer_ac)}",
        f"V_ac_sum={v_ac}",
        f"V_bc_terms={(linear_bc, quadratic_bc, minimizer_bc)}",
        f"V_bc_sum={v_bc}",
        "SO3_algebraic_q_component_constraints_used_in_final_reduction=none",
        "symmetry_pullback=F^*(h_a,h_c,h_ac)=(h_b,h_c,h_bc)",
        "symmetry_pushforward=F_*(h_a,h_c,h_ac)=(-h_b,h_c,-h_bc)",
        "plane_coordinate_transport=(z1,z2,z3,z4)->(z3,z4,z1,z2)",
        f"hostile_drop_Vbc_minimizer={hostile_drop_bc}",
        f"hostile_flip_Vbc_minimizer={hostile_flip_bc}",
        "conclusion=V_ac and V_bc vanish globally on M_generic relative to the reconstructed formulas",
        "nonconsequence=the remaining quadratic identities, cubic integral, smooth extension, and headline theorem remain unaudited",
        f"script_sha256={sha256(script_path.read_bytes()).hexdigest()}",
        f"geometry_dependency_sha256={sha256(geometry_path.read_bytes()).hexdigest()}",
        f"point_operator_dependency_sha256={sha256(point_path.read_bytes()).hexdigest()}",
    ]
    return rows


def main() -> None:
    rows = audit()
    print("\n".join(rows))
    print(f"transcript_sha256={sha256(chr(10).join(rows).encode()).hexdigest()}")


if __name__ == "__main__":
    main()
