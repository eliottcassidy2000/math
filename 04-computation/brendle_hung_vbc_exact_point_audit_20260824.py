#!/usr/bin/env python3
"""Exact independent point audit of the overwritten Brendle--Hung ``V_bc``.

The companion notebook for arXiv:2608.19068v1 forms a raw expression for
``Vbc`` and then accidentally replaces it by ``FullSimplify[Vac,...]`` before
the saved ``Vbc==0`` assertion.  This script does not reuse that cached
output.  It rebuilds the paper's moving-frame operators and checks the three
terms of ``V_bc`` at a genuinely generic exact point.  All three terms are
nonzero and cancel exactly, so the check is sensitive to dropping or
duplicating a summand.

Status: FINITE-EXACT point control.  The separate global symbolic audit now
proves ``V_ac=V_bc=0`` relative to the reconstructed formulas.  It also
reconciles the two symmetry conventions: pullback by the quarter-turn plus
second-factor antipodal map sends ``(h_a,h_c,h_ac)`` to
``(h_b,h_c,h_bc)``, while pushforward sends it to
``(-h_b,h_c,-h_bc)``.  This script remains the sensitive generic-point
hostile for a dropped, duplicated, or sign-flipped summand.

Reproduce from the repository root with

    python3 04-computation/brendle_hung_vbc_exact_point_audit_20260824.py
    python3 -O 04-computation/brendle_hung_vbc_exact_point_audit_20260824.py
"""

from __future__ import annotations

from functools import lru_cache
from hashlib import sha256

import sympy as sp

import brendle_hung_lemma54_independent_audit_20260824 as bh


def tensor_operators(tensor: sp.Matrix):
    """Return exact covariant, ``T``, and linearized-curvature operators."""

    @lru_cache(maxsize=None)
    def covariant(i: int, j: int, k: int) -> sp.Expr:
        return bh.frame_derivative(i, tensor[j, k]) - sum(
            bh.CONNECTION[i][j][ell] * tensor[ell, k]
            + bh.CONNECTION[i][k][ell] * tensor[j, ell]
            for ell in range(4)
        )

    @lru_cache(maxsize=None)
    def covariant2(i: int, j: int, k: int, ell: int) -> sp.Expr:
        return bh.frame_derivative(i, covariant(j, k, ell)) - sum(
            bh.CONNECTION[i][j][m] * covariant(m, k, ell)
            + bh.CONNECTION[i][k][m] * covariant(j, m, ell)
            + bh.CONNECTION[i][ell][m] * covariant(j, k, m)
            for m in range(4)
        )

    @lru_cache(maxsize=None)
    def t_tensor(i: int, j: int, k: int) -> sp.Expr:
        return sp.Rational(1, 2) * (
            covariant(i, j, k) + covariant(j, k, i) - covariant(k, i, j)
        )

    @lru_cache(maxsize=None)
    def linearized(i: int, j: int, k: int, ell: int) -> sp.Expr:
        entry = sp.Rational(1, 2) * (
            -covariant2(i, k, j, ell)
            + covariant2(i, ell, j, k)
            + covariant2(j, k, i, ell)
            - covariant2(j, ell, i, k)
        )
        entry += sum(
            sp.Rational(1, 2)
            * bh.INVERSE_METRIC[m, m]
            * bh.riemann(i, j, k, m)
            * tensor[m, ell]
            - sp.Rational(1, 2)
            * bh.INVERSE_METRIC[m, m]
            * bh.riemann(i, j, ell, m)
            * tensor[k, m]
            for m in range(4)
        )
        return entry

    return t_tensor, linearized


def p2_background() -> sp.Matrix:
    """The Hessian ``P2(Rm)`` in the notebook's plane coordinates."""

    r = bh.riemann
    return sp.Matrix(
        [
            [
                2 * r(0, 3, 0, 3),
                2 * r(0, 1, 2, 3) + 2 * r(0, 3, 2, 1),
                2 * r(0, 3, 1, 3),
                2 * r(0, 3, 2, 0),
            ],
            [
                2 * r(0, 1, 2, 3) + 2 * r(0, 3, 2, 1),
                2 * r(2, 1, 2, 1),
                2 * r(1, 3, 2, 1),
                2 * r(2, 0, 2, 1),
            ],
            [
                2 * r(0, 3, 1, 3),
                2 * r(1, 3, 2, 1),
                2 * r(1, 3, 1, 3),
                2 * r(1, 0, 2, 3) + 2 * r(1, 3, 2, 0),
            ],
            [
                2 * r(0, 3, 2, 0),
                2 * r(2, 0, 2, 1),
                2 * r(1, 0, 2, 3) + 2 * r(1, 3, 2, 0),
                2 * r(2, 0, 2, 0),
            ],
        ]
    )


def q_p0(t_first, t_second) -> sp.Expr:
    """Return ``P0(Q(first,second))=Q_3434`` in zero-based indices."""

    i, j, k, ell = 2, 3, 2, 3
    return sum(
        bh.INVERSE_METRIC[m, m]
        * (
            t_first(j, k, m) * t_second(i, ell, m) / 2
            - t_first(i, k, m) * t_second(j, ell, m) / 2
            + t_second(j, k, m) * t_first(i, ell, m) / 2
            - t_second(i, k, m) * t_first(j, ell, m) / 2
        )
        for m in range(4)
    )


def audit() -> list[str]:
    ey = sp.Matrix([[0, 1, 0]])
    ez = sp.Matrix([[0, 0, 1]])
    k_ey = bh.k_covector(ey)
    j_ey = bh.j_covector(ey)
    k_ez = bh.k_covector(ez)

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
    h_bc = (
        sp.Rational(1, 10)
        * sp.sin(bh.THETA)
        * bh.Q2Z
        * (
            -5
            + 14 * sp.cos(2 * bh.THETA)
            + 3 * sp.cos(4 * bh.THETA)
        )
        * (k_ey.T * k_ez + k_ez.T * k_ey)
        / 2
    ).applyfunc(sp.trigsimp)

    t_b, l_b = tensor_operators(h_b)
    t_c, l_c = tensor_operators(h_c)
    _, l_bc = tensor_operators(h_bc)

    # Columns of this rational SO(3) matrix are q1,q2,q3.  It comes from the
    # quaternion (1,2,3,4), so orthogonality and orientation are exact.
    point = {
        bh.THETA: sp.pi / 6,
        bh.Q1X: -sp.Rational(2, 3),
        bh.Q1Y: sp.Rational(2, 3),
        bh.Q1Z: sp.Rational(1, 3),
        bh.Q2X: sp.Rational(2, 15),
        bh.Q2Y: -sp.Rational(1, 3),
        bh.Q2Z: sp.Rational(14, 15),
        bh.Q3X: sp.Rational(11, 15),
        bh.Q3Y: sp.Rational(2, 3),
        bh.Q3Z: sp.Rational(2, 15),
    }

    hessian = p2_background().applyfunc(
        lambda entry: sp.simplify(sp.trigsimp(entry.subs(point)))
    )
    determinant = sp.factor(hessian.det())
    r_b = sp.Matrix(
        [
            sp.simplify(sp.trigsimp(2 * l_b(*index).subs(point)))
            for index in bh.P1_INDICES
        ]
    )
    r_c = sp.Matrix(
        [
            sp.simplify(sp.trigsimp(2 * l_c(*index).subs(point)))
            for index in bh.P1_INDICES
        ]
    )
    z_c = -hessian.inv() * r_c

    linear_term = sp.simplify(sp.trigsimp(l_bc(2, 3, 2, 3).subs(point)))
    quadratic_term = sp.simplify(sp.trigsimp(q_p0(t_b, t_c).subs(point)))
    minimizer_term = sp.simplify(sp.Rational(1, 2) * r_b.dot(z_c))
    total = sp.simplify(linear_term + quadratic_term + minimizer_term)

    expected = {
        "determinant": sp.Rational(28672, 2278125),
        "r_b": sp.Matrix(
            [0, 0, 37 * sp.sqrt(3) / 150, 7 * sp.sqrt(3) / 27]
        ),
        "r_c": sp.Matrix(
            [
                sp.Rational(32, 1215),
                sp.Rational(56, 6075),
                896 * sp.sqrt(3) / 84375,
                1232 * sp.sqrt(3) / 151875,
            ]
        ),
        "linear": sp.Rational(581, 20250),
        "quadratic": -sp.Rational(1064, 50625),
        "minimizer": -sp.Rational(259, 33750),
        "total": sp.S.Zero,
    }
    actual = {
        "determinant": determinant,
        "r_b": r_b,
        "r_c": r_c,
        "linear": linear_term,
        "quadratic": quadratic_term,
        "minimizer": minimizer_term,
        "total": total,
    }
    for key, expected_value in expected.items():
        if actual[key] != expected_value:
            raise RuntimeError(f"{key} mismatch: {actual[key]} != {expected_value}")

    return [
        "status=FINITE-EXACT",
        "source=arXiv:2608.19068v1 Proposition 5.2 V_bc identity",
        "universe=theta=pi/6; q-frame=SO(3) quaternion (1,2,3,4)",
        f"hessian_determinant={determinant}",
        f"r_b={tuple(r_b)}",
        f"r_c={tuple(r_c)}",
        f"P0_L_hbc={linear_term}",
        f"P0_Q_hb_hc={quadratic_term}",
        f"half_r_b_dot_z_c={minimizer_term}",
        f"V_bc_sum={total}",
        "hostile_control=all three displayed summands are nonzero",
        "conclusion=the overwritten identity passes one independent exact generic point",
        "nonconsequence=this point control alone does not prove V_bc identically zero",
    ]


def main() -> None:
    rows = audit()
    print("\n".join(rows))
    print(f"transcript_sha256={sha256(chr(10).join(rows).encode()).hexdigest()}")


if __name__ == "__main__":
    main()
