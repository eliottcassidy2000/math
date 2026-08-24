#!/usr/bin/env python3
"""Independent exact audit of the missing ``h_a`` branch in BH Lemma 5.4.

Source under audit
------------------
Simon Brendle and Pei-Ken Hung, *A metric on S^2 x S^2 with positive
sectional curvature*, arXiv:2608.19068v1, especially Sections 3--5 and
Lemma 5.4.  The v1 PDF downloaded on 2026-08-24 has SHA-256

    46bb66824fa239b95dbfa6a661de675a2d15d0a1e38332433827bf9f31f4e4f7.

The attached Mathematica notebook's cell advertised as checking
``z(h_a)=z(h_b)=z(h_c)=0`` on the two components of Sigma instead assigns
``ra=P1[L[hc]]``.  Thus that cell checks the ``h_c`` branch twice.  This
script independently reconstructs the moving-frame metric, connection,
curvature, linearized-curvature operator ``L``, and ``h_a`` from the paper.
Since the normal Hessian ``H`` is invertible on the generic locus,
``z(h_a)=0`` is equivalent to ``P1(L(h_a))=0``.

Status and scope
----------------
FINITE-EXACT.  SymPy performs exact symbolic differentiation and arithmetic.
The script checks both parametrized components of the residual torus Sigma
and exact point controls.  It verifies only the omitted Lemma 5.4 branch; it
does not independently replay the whole Mathematica notebook or prove the
paper's positive-curvature theorem.

Reproduce from the repository root with

    python3 04-computation/brendle_hung_lemma54_independent_audit_20260824.py
    python3 -O 04-computation/brendle_hung_lemma54_independent_audit_20260824.py
"""

from __future__ import annotations

from functools import lru_cache
from hashlib import sha256

import sympy as sp


THETA, PHI = sp.symbols("theta phi", real=True)
(
    Q1X,
    Q1Y,
    Q1Z,
    Q2X,
    Q2Y,
    Q2Z,
    Q3X,
    Q3Y,
    Q3Z,
) = sp.symbols("q1x q1y q1z q2x q2y q2z q3x q3y q3z", real=True)

Q1 = sp.Matrix([[Q1X, Q1Y, Q1Z]])
Q2 = sp.Matrix([[Q2X, Q2Y, Q2Z]])
Q3 = sp.Matrix([[Q3X, Q3Y, Q3Z]])

METRIC = sp.diag(
    1 / (2 - sp.cos(2 * THETA)),
    1 / (2 + sp.cos(2 * THETA)),
    sp.Rational(1, 3),
    1,
)
INVERSE_METRIC = METRIC.inv()


def _zero_cube() -> list[list[list[sp.Expr]]]:
    return [[[sp.S.Zero for _ in range(4)] for _ in range(4)] for _ in range(4)]


STRUCTURE = _zero_cube()
for index, value in {
    (0, 1, 2): 1 / (sp.sin(THETA) * sp.cos(THETA)),
    (1, 0, 2): -1 / (sp.sin(THETA) * sp.cos(THETA)),
    (0, 2, 1): -sp.cos(THETA) / sp.sin(THETA),
    (2, 0, 1): sp.cos(THETA) / sp.sin(THETA),
    (0, 3, 0): sp.cos(THETA) / sp.sin(THETA),
    (3, 0, 0): -sp.cos(THETA) / sp.sin(THETA),
    (1, 2, 0): sp.sin(THETA) / sp.cos(THETA),
    (2, 1, 0): -sp.sin(THETA) / sp.cos(THETA),
    (1, 3, 1): -sp.sin(THETA) / sp.cos(THETA),
    (3, 1, 1): sp.sin(THETA) / sp.cos(THETA),
}.items():
    STRUCTURE[index[0]][index[1]][index[2]] = value

STRUCTURE_LOWER = [
    [
        [STRUCTURE[i][j][k] * METRIC[k, k] for k in range(4)]
        for j in range(4)
    ]
    for i in range(4)
]


@lru_cache(maxsize=None)
def frame_derivative(i: int, expression: sp.Expr) -> sp.Expr:
    """Differentiate in the paper's moving frame ``m_1,...,m_4``."""

    if i == 0:
        return (
            Q3X * sp.diff(expression, Q2X)
            + Q3Y * sp.diff(expression, Q2Y)
            + Q3Z * sp.diff(expression, Q2Z)
            - Q2X * sp.diff(expression, Q3X)
            - Q2Y * sp.diff(expression, Q3Y)
            - Q2Z * sp.diff(expression, Q3Z)
        ) / sp.sin(THETA)
    if i == 1:
        return (
            Q1X * sp.diff(expression, Q3X)
            + Q1Y * sp.diff(expression, Q3Y)
            + Q1Z * sp.diff(expression, Q3Z)
            - Q3X * sp.diff(expression, Q1X)
            - Q3Y * sp.diff(expression, Q1Y)
            - Q3Z * sp.diff(expression, Q1Z)
        ) / sp.cos(THETA)
    if i == 2:
        return (
            Q2X * sp.diff(expression, Q1X)
            + Q2Y * sp.diff(expression, Q1Y)
            + Q2Z * sp.diff(expression, Q1Z)
            - Q1X * sp.diff(expression, Q2X)
            - Q1Y * sp.diff(expression, Q2Y)
            - Q1Z * sp.diff(expression, Q2Z)
        )
    return sp.diff(expression, THETA)


CONNECTION = _zero_cube()
for i in range(4):
    for j in range(4):
        for k in range(4):
            entry = sp.Rational(1, 2) * (
                frame_derivative(i, METRIC[j, k])
                + frame_derivative(j, METRIC[i, k])
                - frame_derivative(k, METRIC[i, j])
                + STRUCTURE_LOWER[i][j][k]
                - STRUCTURE_LOWER[j][k][i]
                + STRUCTURE_LOWER[k][i][j]
            ) * INVERSE_METRIC[k, k]
            CONNECTION[i][j][k] = sp.trigsimp(entry)


@lru_cache(maxsize=None)
def riemann(i: int, j: int, k: int, ell: int) -> sp.Expr:
    entry = frame_derivative(i, CONNECTION[j][ell][k]) - frame_derivative(
        j, CONNECTION[i][ell][k]
    )
    entry += sum(
        CONNECTION[j][ell][r] * CONNECTION[i][r][k]
        - CONNECTION[i][ell][r] * CONNECTION[j][r][k]
        - STRUCTURE[i][j][r] * CONNECTION[r][ell][k]
        for r in range(4)
    )
    return sp.trigsimp(entry * METRIC[k, k])


def k_covector(vector: sp.Matrix) -> sp.Matrix:
    return sp.Matrix(
        [[
            sp.sin(THETA) * (Q1 * vector.T)[0],
            sp.cos(THETA) * (Q2 * vector.T)[0],
            (Q3 * vector.T)[0],
            0,
        ]]
    ) * METRIC


def j_covector(vector: sp.Matrix) -> sp.Matrix:
    return sp.Matrix(
        [[
            sp.cos(THETA) * (Q2 * vector.T)[0],
            sp.sin(THETA) * (Q1 * vector.T)[0],
            0,
            -(Q3 * vector.T)[0],
        ]]
    ) * METRIC


EX = sp.Matrix([[1, 0, 0]])
K_EX = k_covector(EX)
J_EX = j_covector(EX)
K_PLUS = sp.Matrix(
    [[2 * sp.cos(THETA) * sp.sin(THETA), 0, 0, 0]]
) * METRIC
K_MINUS = sp.Matrix(
    [[0, 2 * sp.cos(THETA) * sp.sin(THETA), 0, 0]]
) * METRIC

H_A = (
    sp.Rational(6, 5)
    * (2 - sp.cos(2 * THETA))
    * (K_PLUS.T * K_EX + K_EX.T * K_PLUS)
    / 2
    - 2
    * (1 - sp.cos(2 * THETA))
    * (2 + sp.cos(2 * THETA))
    * (K_MINUS.T * J_EX + J_EX.T * K_MINUS)
    / 2
).applyfunc(sp.trigsimp)

# A deliberately non-special symmetric tensor.  It is a hostile control for
# the moving-frame and P1 pipelines: unlike h_a, it has nonzero P1(L(h)) on
# Sigma at the exact point used below.
H_HOSTILE = sp.zeros(4)
H_HOSTILE[0, 3] = H_HOSTILE[3, 0] = 1


@lru_cache(maxsize=None)
def covariant_h(i: int, j: int, k: int) -> sp.Expr:
    return frame_derivative(i, H_A[j, k]) - sum(
        CONNECTION[i][j][ell] * H_A[ell, k]
        + CONNECTION[i][k][ell] * H_A[j, ell]
        for ell in range(4)
    )


@lru_cache(maxsize=None)
def covariant2_h(i: int, j: int, k: int, ell: int) -> sp.Expr:
    return frame_derivative(i, covariant_h(j, k, ell)) - sum(
        CONNECTION[i][j][m] * covariant_h(m, k, ell)
        + CONNECTION[i][k][m] * covariant_h(j, m, ell)
        + CONNECTION[i][ell][m] * covariant_h(j, k, m)
        for m in range(4)
    )


@lru_cache(maxsize=None)
def linearized_curvature(i: int, j: int, k: int, ell: int) -> sp.Expr:
    entry = sp.Rational(1, 2) * (
        -covariant2_h(i, k, j, ell)
        + covariant2_h(i, ell, j, k)
        + covariant2_h(j, k, i, ell)
        - covariant2_h(j, ell, i, k)
    )
    entry += sum(
        sp.Rational(1, 2)
        * INVERSE_METRIC[m, m]
        * riemann(i, j, k, m)
        * H_A[m, ell]
        - sp.Rational(1, 2)
        * INVERSE_METRIC[m, m]
        * riemann(i, j, ell, m)
        * H_A[k, m]
        for m in range(4)
    )
    return entry


@lru_cache(maxsize=None)
def hostile_covariant_h(i: int, j: int, k: int) -> sp.Expr:
    return frame_derivative(i, H_HOSTILE[j, k]) - sum(
        CONNECTION[i][j][ell] * H_HOSTILE[ell, k]
        + CONNECTION[i][k][ell] * H_HOSTILE[j, ell]
        for ell in range(4)
    )


@lru_cache(maxsize=None)
def hostile_covariant2_h(i: int, j: int, k: int, ell: int) -> sp.Expr:
    return frame_derivative(i, hostile_covariant_h(j, k, ell)) - sum(
        CONNECTION[i][j][m] * hostile_covariant_h(m, k, ell)
        + CONNECTION[i][k][m] * hostile_covariant_h(j, m, ell)
        + CONNECTION[i][ell][m] * hostile_covariant_h(j, k, m)
        for m in range(4)
    )


@lru_cache(maxsize=None)
def hostile_linearized_curvature(i: int, j: int, k: int, ell: int) -> sp.Expr:
    entry = sp.Rational(1, 2) * (
        -hostile_covariant2_h(i, k, j, ell)
        + hostile_covariant2_h(i, ell, j, k)
        + hostile_covariant2_h(j, k, i, ell)
        - hostile_covariant2_h(j, ell, i, k)
    )
    entry += sum(
        sp.Rational(1, 2)
        * INVERSE_METRIC[m, m]
        * riemann(i, j, k, m)
        * H_HOSTILE[m, ell]
        - sp.Rational(1, 2)
        * INVERSE_METRIC[m, m]
        * riemann(i, j, ell, m)
        * H_HOSTILE[k, m]
        for m in range(4)
    )
    return entry


P1_INDICES = (
    (0, 3, 2, 3),
    (2, 1, 2, 3),
    (1, 3, 2, 3),
    (2, 0, 2, 3),
)

PLUS_TORUS = {
    Q1X: sp.cos(PHI),
    Q1Y: sp.sin(PHI),
    Q1Z: 0,
    Q2X: -sp.sin(PHI),
    Q2Y: sp.cos(PHI),
    Q2Z: 0,
    Q3X: 0,
    Q3Y: 0,
    Q3Z: 1,
}
MINUS_TORUS = {
    Q1X: sp.cos(PHI),
    Q1Y: -sp.sin(PHI),
    Q1Z: 0,
    Q2X: -sp.sin(PHI),
    Q2Y: -sp.cos(PHI),
    Q2Z: 0,
    Q3X: 0,
    Q3Y: 0,
    Q3Z: -1,
}


def exact_zero(expression: sp.Expr) -> bool:
    """Normalize trig/rational structure and test exact zero."""

    return sp.trigsimp(sp.cancel(sp.expand_trig(expression))) == 0


def audit() -> list[str]:
    p1_lha = [2 * linearized_curvature(*index) for index in P1_INDICES]
    plus = [sp.trigsimp(entry.subs(PLUS_TORUS)) for entry in p1_lha]
    minus = [sp.trigsimp(entry.subs(MINUS_TORUS)) for entry in p1_lha]

    plus_zero = [exact_zero(entry) for entry in plus]
    minus_zero = [exact_zero(entry) for entry in minus]

    exact_points = (
        ("plus_theta_pi6_phi0", PLUS_TORUS | {THETA: sp.pi / 6, PHI: 0}),
        (
            "plus_theta_pi5_phi_pi4",
            PLUS_TORUS | {THETA: sp.pi / 5, PHI: sp.pi / 4},
        ),
        ("minus_theta_pi6_phi0", MINUS_TORUS | {THETA: sp.pi / 6, PHI: 0}),
        (
            "minus_theta_pi5_phi_pi4",
            MINUS_TORUS | {THETA: sp.pi / 5, PHI: sp.pi / 4},
        ),
    )
    point_rows: list[str] = []
    for label, substitution in exact_points:
        values = tuple(
            sp.simplify(sp.trigsimp(entry.subs(substitution))) for entry in p1_lha
        )
        point_rows.append(f"control_{label}={values}")
        if any(value != 0 for value in values):
            raise RuntimeError(f"exact point control failed: {label}: {values}")

    if not all(plus_zero + minus_zero):
        raise RuntimeError(
            "symbolic torus identity failed: "
            f"plus={plus_zero}, minus={minus_zero}, values={plus!r}/{minus!r}"
        )

    hostile_substitution = PLUS_TORUS | {
        THETA: sp.pi / 6,
        PHI: sp.pi / 4,
    }
    hostile_values = tuple(
        sp.simplify(
            sp.trigsimp(
                2
                * hostile_linearized_curvature(*index).subs(hostile_substitution)
            )
        )
        for index in P1_INDICES
    )
    hostile_expected = (0, 0, -sp.Rational(1, 3), sp.Rational(2, 9))
    if hostile_values != hostile_expected:
        raise RuntimeError(
            f"hostile tensor control failed: {hostile_values} != {hostile_expected}"
        )

    return [
        "status=FINITE-EXACT",
        "source=arXiv:2608.19068v1 Lemma 5.4",
        "target=P1(L(h_a)) on both components of Sigma",
        f"plus_component_zero={tuple(plus_zero)}",
        f"minus_component_zero={tuple(minus_zero)}",
        *point_rows,
        f"hostile_constant_h14={hostile_values}",
        "conclusion=the omitted h_a branch is independently verified",
        "nonconsequence=the rest of the notebook and Theorem 1.1 remain unaudited",
    ]


def main() -> None:
    rows = audit()
    print("\n".join(rows))
    print(f"transcript_sha256={sha256(chr(10).join(rows).encode()).hexdigest()}")


if __name__ == "__main__":
    main()
