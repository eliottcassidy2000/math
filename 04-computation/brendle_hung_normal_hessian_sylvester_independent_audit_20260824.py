#!/usr/bin/env python3
"""Independent Sylvester audit of the Brendle--Hung normal Hessian.

Source under audit
------------------
Simon Brendle and Pei-Ken Hung, *A metric on S^2 x S^2 with positive
sectional curvature*, arXiv:2608.19068v1.  In the notation of the paper,
``H=P_2 R^(0)`` is the Hessian in the ordered normal-plane coordinates
``(z^1,z^2,z^3,z^4)`` determined by

    m_3(z)=m_3+z^1 m_1+z^3 m_2,
    m_4(z)=m_4+z^4 m_1+z^2 m_2.

This script reconstructs the background metric, noncoordinate-frame brackets,
Levi--Civita connection, Riemann tensor, and ``P_2`` directly after the exact
change of variable ``x=tan(theta)>0``.  It imports no repository geometry
implementation.  It verifies the complete Hessian matrix and all four leading
principal minors.  Their factorizations prove positive definiteness by
Sylvester's criterion on the generic locus.

Scope
-----
VERIFIED-EXACT relative to this independent reconstruction of the displayed
background formulas.  The endpoint hostile records degeneration as
``x->0+`` and ``x->infinity``; no endpoint extension or uniform coercivity is
claimed.  A determinant-positive indefinite matrix is included to show why
the four-minor gate is genuinely stronger than the determinant check alone.

Reproduce from the repository root with

    python 04-computation/brendle_hung_normal_hessian_sylvester_independent_audit_20260824.py
    python -O 04-computation/brendle_hung_normal_hessian_sylvester_independent_audit_20260824.py
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import sympy as sp


X = sp.symbols("x", positive=True)
ZERO = sp.S.Zero
ONE = sp.S.One


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def normalized(expression: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(expression))


def exact_zero(expression: sp.Expr) -> bool:
    return normalized(expression) == 0


def zero_cube() -> list[list[list[sp.Expr]]]:
    return [
        [[ZERO for _ in range(4)] for _ in range(4)] for _ in range(4)
    ]


def frame_derivative(index: int, expression: sp.Expr) -> sp.Expr:
    """Differentiate a background coefficient in the ``m_i`` frame."""

    if index == 3:
        return sp.expand((1 + X**2) * sp.diff(expression, X))
    return ZERO


def reconstruct() -> tuple[sp.Matrix, tuple[sp.Expr, ...]]:
    """Rebuild ``H=P_2 R^(0)`` and its leading principal minors."""

    metric_entries = (
        (1 + X**2) / (1 + 3 * X**2),
        (1 + X**2) / (3 + X**2),
        sp.Rational(1, 3),
        ONE,
    )
    metric = sp.diag(*metric_entries)
    inverse_metric = metric.inv()

    # These are exactly the source brackets after x=tan(theta):
    # [m1,m2], [m1,m3], [m1,m4], [m2,m3], [m2,m4].
    structure = zero_cube()
    bracket_entries = {
        (0, 1, 2): (1 + X**2) / X,
        (0, 2, 1): -1 / X,
        (0, 3, 0): 1 / X,
        (1, 2, 0): X,
        (1, 3, 1): -X,
    }
    for (i, j, k), value in bracket_entries.items():
        structure[i][j][k] = value
        structure[j][i][k] = -value

    def lower_structure(i: int, j: int, k: int) -> sp.Expr:
        return structure[i][j][k] * metric[k, k]

    connection = zero_cube()
    for i in range(4):
        for j in range(4):
            for k in range(4):
                entry = sp.Rational(1, 2) * (
                    frame_derivative(i, metric[j, k])
                    + frame_derivative(j, metric[i, k])
                    - frame_derivative(k, metric[i, j])
                    + lower_structure(i, j, k)
                    - lower_structure(j, k, i)
                    + lower_structure(k, i, j)
                )
                connection[i][j][k] = normalized(
                    entry * inverse_metric[k, k]
                )

    def riemann(i: int, j: int, k: int, ell: int) -> sp.Expr:
        entry = frame_derivative(i, connection[j][ell][k])
        entry -= frame_derivative(j, connection[i][ell][k])
        entry += sum(
            connection[j][ell][r] * connection[i][r][k]
            - connection[i][ell][r] * connection[j][r][k]
            - structure[i][j][r] * connection[r][ell][k]
            for r in range(4)
        )
        return normalized(entry * metric[k, k])

    # This is the source P_2 ordering (z^1,z^2,z^3,z^4), with paper
    # indices shifted from 1..4 to 0..3.
    hessian = sp.Matrix(
        [
            [
                2 * riemann(0, 3, 0, 3),
                2 * riemann(0, 1, 2, 3) + 2 * riemann(0, 3, 2, 1),
                2 * riemann(0, 3, 1, 3),
                2 * riemann(0, 3, 2, 0),
            ],
            [
                2 * riemann(0, 1, 2, 3) + 2 * riemann(0, 3, 2, 1),
                2 * riemann(2, 1, 2, 1),
                2 * riemann(1, 3, 2, 1),
                2 * riemann(2, 0, 2, 1),
            ],
            [
                2 * riemann(0, 3, 1, 3),
                2 * riemann(1, 3, 2, 1),
                2 * riemann(1, 3, 1, 3),
                2 * riemann(1, 0, 2, 3) + 2 * riemann(1, 3, 2, 0),
            ],
            [
                2 * riemann(0, 3, 2, 0),
                2 * riemann(2, 0, 2, 1),
                2 * riemann(1, 0, 2, 3) + 2 * riemann(1, 3, 2, 0),
                2 * riemann(2, 0, 2, 0),
            ],
        ]
    ).applyfunc(normalized)

    minors = tuple(
        normalized(hessian[:size, :size].det()) for size in range(1, 5)
    )
    return hessian, minors


def audit() -> list[str]:
    hessian, minors = reconstruct()

    expected_hessian = sp.Matrix(
        [
            [
                2 * (X**2 + 1) ** 2 * (3 * X**2 + 7)
                / (3 * X**2 + 1) ** 3,
                2 * (X**2 + 1) ** 2 * (3 * X**2 + 13)
                / (3 * (X**2 + 3) * (3 * X**2 + 1) ** 2),
                0,
                0,
            ],
            [
                2 * (X**2 + 1) ** 2 * (3 * X**2 + 13)
                / (3 * (X**2 + 3) * (3 * X**2 + 1) ** 2),
                2 * (X**2 + 1) * (3 * X**4 + 34 * X**2 + 27)
                / (9 * (X**2 + 3) ** 2 * (3 * X**2 + 1)),
                0,
                0,
            ],
            [
                0,
                0,
                2 * (X**2 + 1) ** 2 * (7 * X**2 + 3)
                / (X**2 + 3) ** 3,
                2 * (X**2 + 1) ** 2 * (13 * X**2 + 3)
                / (3 * (X**2 + 3) ** 2 * (3 * X**2 + 1)),
            ],
            [
                0,
                0,
                2 * (X**2 + 1) ** 2 * (13 * X**2 + 3)
                / (3 * (X**2 + 3) ** 2 * (3 * X**2 + 1)),
                2 * (X**2 + 1) * (27 * X**4 + 34 * X**2 + 3)
                / (9 * (X**2 + 3) * (3 * X**2 + 1) ** 2),
            ],
        ]
    ).applyfunc(normalized)

    expected_minors = (
        2 * (X**2 + 1) ** 2 * (3 * X**2 + 7)
        / (3 * X**2 + 1) ** 3,
        16 * (X**2 + 1) ** 3 * (3 * X**2 + 5)
        / (9 * (X**2 + 3) ** 2 * (3 * X**2 + 1) ** 3),
        32 * (X**2 + 1) ** 5 * (3 * X**2 + 5) * (7 * X**2 + 3)
        / (9 * (X**2 + 3) ** 5 * (3 * X**2 + 1) ** 3),
        256
        * X**2
        * (X**2 + 1) ** 6
        * (3 * X**2 + 5)
        * (5 * X**2 + 3)
        / (81 * (X**2 + 3) ** 5 * (3 * X**2 + 1) ** 5),
    )
    expected_minors = tuple(normalized(value) for value in expected_minors)

    require(hessian == expected_hessian, "reconstructed Hessian changed")
    require(hessian == hessian.T, "reconstructed Hessian is not symmetric")
    require(
        all(exact_zero(actual - expected) for actual, expected in zip(minors, expected_minors)),
        "leading principal minor factorization changed",
    )
    require(
        all(value.is_positive is True for value in expected_minors),
        "SymPy did not certify every factorized minor positive for x>0",
    )

    # Exact positive control in the interior of the generic interval.
    point_minors = tuple(sp.factor(value.subs(X, 1)) for value in minors)
    require(
        all(value.is_positive is True for value in point_minors),
        "x=1 positive control failed",
    )

    # Boundary hostile: the full determinant degenerates at both deleted ends.
    determinant = minors[-1]
    endpoint_zero = (
        sp.limit(determinant, X, 0, dir="+"),
        sp.limit(determinant, X, sp.oo),
    )
    require(endpoint_zero == (0, 0), "generic-boundary degeneration changed")

    # Determinant positivity by itself does not imply positive definiteness.
    indefinite = sp.diag(1, -1, -1, 1)
    indefinite_det = indefinite.det()
    indefinite_delta2 = indefinite[:2, :2].det()
    require(
        (indefinite_det, indefinite_delta2) == (1, -1),
        "determinant-only indefinite hostile failed",
    )

    # A single source-entry sign mutation must be visible to the exact gate.
    mutated = hessian.copy()
    mutated[0, 0] = -mutated[0, 0]
    mutated_delta1 = normalized(mutated[0, 0])
    mutation_matches = exact_zero(mutated_delta1 - expected_minors[0])
    require(not mutation_matches, "diagonal sign mutation escaped the gate")
    require(
        exact_zero(mutated_delta1 + expected_minors[0]),
        "diagonal sign mutation did not negate Delta_1",
    )

    script_path = Path(__file__)
    rows = [
        "status=VERIFIED-EXACT",
        "source=arXiv:2608.19068v1 background metric and P2 normal Hessian",
        "universe=0<theta<pi/2; x=tan(theta)>0",
        "method=independent x-coordinate moving-frame reconstruction; no repository geometry imports",
        "basis=(z1,z2,z3,z4) in the source P2 order",
        f"hessian={hessian}",
        "hessian_symmetric=True",
        "hessian_block_structure=(z1,z2)|(z3,z4)",
        f"leading_principal_minors={minors}",
        "sylvester_conclusion=H is positive definite for every x>0",
        f"x_equals_1_minors={point_minors}",
        f"endpoint_det_limits_x0_xinf={endpoint_zero}",
        f"det_positive_indefinite_hostile=(det={indefinite_det},delta2={indefinite_delta2})",
        f"diagonal_sign_mutation_delta1={mutated_delta1}",
        f"diagonal_sign_mutation_matches_expected={mutation_matches}",
        "conclusion=the reconstructed generic normal Hessian is positive definite, hence invertible",
        "nonconsequence=no endpoint extension, uniform coercivity, remaining curvature identity, cubic integral, or headline theorem",
        f"script_sha256={sha256(script_path.read_bytes()).hexdigest()}",
    ]
    return rows


def main() -> None:
    rows = audit()
    print("\n".join(rows))
    print(f"transcript_sha256={sha256(chr(10).join(rows).encode()).hexdigest()}")


if __name__ == "__main__":
    main()
