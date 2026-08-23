#!/usr/bin/env python3
"""Exact leading-curl companion for THM-3709's Cohn decoration gate."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


x, y = sp.symbols("x y")
A = 1 + x * y
B = x**2
G = -y**2
D = 1 - x * y
C = sp.Matrix(((A, B), (G, D)))


def eplus(h: sp.Expr) -> sp.Matrix:
    return sp.Matrix(((1, h), (0, 1)))


def eminus(h: sp.Expr) -> sp.Matrix:
    return sp.Matrix(((1, 0), (h, 1)))


def curl(row: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(sp.diff(row[0], y) - sp.diff(row[1], x))


def homogeneous_part(expr: sp.Expr, degree: int) -> sp.Expr:
    poly = sp.Poly(sp.expand(expr), x, y)
    return sp.expand(
        sum(
            coefficient * x**powers[0] * y**powers[1]
            for powers, coefficient in poly.terms()
            if sum(powers) == degree
        )
    )


def row(matrix: sp.Matrix, index: int) -> tuple[sp.Expr, sp.Expr]:
    return sp.expand(matrix[index, 0]), sp.expand(matrix[index, 1])


def homogeneous_form(prefix: str, degree: int) -> sp.Expr:
    coefficients = sp.symbols(f"{prefix}0:{degree + 1}")
    return sum(
        coefficients[index] * x**index * y ** (degree - index)
        for index in range(degree + 1)
    )


gate(sp.expand(C.det()) == 1, "Cohn determinant")

c = sp.symbols("c")
semantic_rows: list[str] = []
degree_pairs = tuple((p, q) for p in range(1, 5) for q in range(1, 5))

for p, q in degree_pairs:
    up = homogeneous_form(f"u{p}_{q}_", p)
    wq = homogeneous_form(f"w{p}_{q}_", q)
    # Constants are hostile lower-degree controls.  The proof permits every
    # collection of lower homogeneous terms; none can reach the extracted top.
    u0, w0 = sp.symbols(f"u0_{p}_{q} w0_{p}_{q}")
    u = up + u0
    w = wq + w0

    right_a = eplus(w) * eminus(u)
    right_b = eminus(u) * eplus(w)
    N_a = sp.expand(C * right_a)
    N_b = sp.expand(C * right_b)

    s = B + A * w
    r = D + G * w
    t = A + B * u
    z = G + D * u
    expected_a = sp.Matrix(((A + u * s, s), (G + u * r, r)))
    expected_b = sp.Matrix(((t, B + w * t), (z, D + w * z)))
    gate(all(sp.expand(value) == 0 for value in N_a - expected_a), "right order A")
    gate(all(sp.expand(value) == 0 for value in N_b - expected_b), "right order B")

    rows_a = (row(N_a, 0), row(N_a, 1))
    rows_b = (row(N_b, 0), row(N_b, 1))
    top_degree = p + q + 2
    H1 = up * x * y * wq
    H2 = -up * y**2 * wq
    K1 = x**2 * up * wq
    K2 = -x * y * up * wq

    gate(sp.expand(homogeneous_part(rows_a[0][0], top_degree) - H1) == 0,
         "A row1 first lead")
    gate(homogeneous_part(rows_a[0][1], top_degree) == 0,
         "A row1 second lower")
    gate(sp.expand(homogeneous_part(rows_a[1][0], top_degree) - H2) == 0,
         "A row2 first lead")
    gate(homogeneous_part(rows_a[1][1], top_degree) == 0,
         "A row2 second lower")
    gate(homogeneous_part(rows_b[0][0], top_degree) == 0,
         "B row1 first lower")
    gate(sp.expand(homogeneous_part(rows_b[0][1], top_degree) - K1) == 0,
         "B row1 second lead")
    gate(homogeneous_part(rows_b[1][0], top_degree) == 0,
         "B row2 first lower")
    gate(sp.expand(homogeneous_part(rows_b[1][1], top_degree) - K2) == 0,
         "B row2 second lead")

    constant_expected = (
        sp.diff(up * y * wq * (-y + c * x), y),
        sp.diff(up * y * wq * (x - c * y), y),
        -sp.diff(x * up * wq * (-y + c * x), x),
        -sp.diff(x * up * wq * (x - c * y), x),
    )
    constant_actual = (
        homogeneous_part(curl((rows_a[1][0] + c * rows_a[0][0], rows_a[1][1] + c * rows_a[0][1])), p + q + 1),
        homogeneous_part(curl((rows_a[0][0] + c * rows_a[1][0], rows_a[0][1] + c * rows_a[1][1])), p + q + 1),
        homogeneous_part(curl((rows_b[1][0] + c * rows_b[0][0], rows_b[1][1] + c * rows_b[0][1])), p + q + 1),
        homogeneous_part(curl((rows_b[0][0] + c * rows_b[1][0], rows_b[0][1] + c * rows_b[1][1])), p + q + 1),
    )
    for actual, expected in zip(constant_actual, constant_expected):
        gate(sp.expand(actual - expected) == 0, "constant leading curl")
        gate(sp.expand(expected) != 0, "constant generic lead vanished")

    for m in range(1, 5):
        hm = homogeneous_form(f"h{p}_{q}_{m}_", m)
        actual = (
            homogeneous_part(curl((rows_a[1][0] + hm * rows_a[0][0], rows_a[1][1] + hm * rows_a[0][1])), m + p + q + 1),
            homogeneous_part(curl((rows_a[0][0] + hm * rows_a[1][0], rows_a[0][1] + hm * rows_a[1][1])), m + p + q + 1),
            homogeneous_part(curl((rows_b[1][0] + hm * rows_b[0][0], rows_b[1][1] + hm * rows_b[0][1])), m + p + q + 1),
            homogeneous_part(curl((rows_b[0][0] + hm * rows_b[1][0], rows_b[0][1] + hm * rows_b[1][1])), m + p + q + 1),
        )
        expected = (
            sp.diff(H1 * hm, y),
            sp.diff(H2 * hm, y),
            -sp.diff(K1 * hm, x),
            -sp.diff(K2 * hm, x),
        )
        for left, right in zip(actual, expected):
            gate(sp.expand(left - right) == 0, "positive-degree leading curl")
            gate(sp.expand(right) != 0, "positive-degree generic lead vanished")

    semantic_rows.append(
        f"{p},{q}:" + hashlib.sha256(
            "|".join(sp.srepr(value) for value in constant_actual).encode()
        ).hexdigest()
    )

# The mixed boundary has exactly one nonconstant right parameter and one
# nonzero scalar parameter.  In the two middle cases the leading row is a
# scalar multiple pair, so the curl becomes a directional derivative.
cc, dd = sp.symbols("cc dd", nonzero=True)
for degree in range(1, 5):
    lead = homogeneous_form(f"v{degree}_", degree)
    lower = sp.symbols(f"vlow_{degree}")
    value = lead + lower

    mixed_cases = (
        (
            "plus-minus-w-constant",
            sp.expand(C * eplus(cc) * eminus(value)),
            (
                (lead * x * (x + cc * y), 0),
                (-lead * y * (x + cc * y), 0),
            ),
            lambda F: sp.diff(F, y),
        ),
        (
            "plus-minus-u-constant",
            sp.expand(C * eplus(value) * eminus(cc)),
            (
                (cc * x * y * lead, x * y * lead),
                (-cc * y**2 * lead, -y**2 * lead),
            ),
            lambda F: cc * sp.diff(F, y) - sp.diff(F, x),
        ),
        (
            "minus-plus-w-constant",
            sp.expand(C * eminus(value) * eplus(cc)),
            (
                (x**2 * lead, cc * x**2 * lead),
                (-x * y * lead, -cc * x * y * lead),
            ),
            lambda F: sp.diff(F, y) - cc * sp.diff(F, x),
        ),
        (
            "minus-plus-u-constant",
            sp.expand(C * eminus(cc) * eplus(value)),
            (
                (0, x * lead * (y + cc * x)),
                (0, -y * lead * (y + cc * x)),
            ),
            lambda F: -sp.diff(F, x),
        ),
    )

    for label, matrix, leading_rows, directional in mixed_cases:
        matrix_rows = (row(matrix, 0), row(matrix, 1))
        top_degree = degree + 2
        for actual_row, expected_row in zip(matrix_rows, leading_rows):
            for actual_entry, expected_entry in zip(actual_row, expected_row):
                gate(
                    sp.expand(homogeneous_part(actual_entry, top_degree) - expected_entry) == 0,
                    f"{label} row lead",
                )

        # The first and fourth cases expose only one leading component.  The
        # middle cases expose (cc*F,F) and (F,cc*F), respectively.  In all
        # cases `scalar_part` is the F to which the displayed operator applies.
        scalar_index = 1 if label in {
            "plus-minus-u-constant",
            "minus-plus-u-constant",
        } else 0
        scalar_parts = (
            leading_rows[0][scalar_index],
            leading_rows[1][scalar_index],
        )
        constant_scalars = (
            sp.expand(scalar_parts[1] + dd * scalar_parts[0]),
            sp.expand(scalar_parts[0] + dd * scalar_parts[1]),
        )
        constant_rows = (
            (
                sp.expand(matrix_rows[1][0] + dd * matrix_rows[0][0]),
                sp.expand(matrix_rows[1][1] + dd * matrix_rows[0][1]),
            ),
            (
                sp.expand(matrix_rows[0][0] + dd * matrix_rows[1][0]),
                sp.expand(matrix_rows[0][1] + dd * matrix_rows[1][1]),
            ),
        )
        for actual_row, scalar in zip(constant_rows, constant_scalars):
            actual = homogeneous_part(curl(actual_row), degree + 1)
            expected = sp.expand(directional(scalar))
            gate(sp.expand(actual - expected) == 0, f"{label} constant-left curl")
            gate(expected != 0, f"{label} constant-left generic lead vanished")

        for left_degree in range(1, 5):
            hm = homogeneous_form(f"hmix_{label}_{degree}_{left_degree}_", left_degree)
            positive_rows = (
                (
                    sp.expand(matrix_rows[1][0] + hm * matrix_rows[0][0]),
                    sp.expand(matrix_rows[1][1] + hm * matrix_rows[0][1]),
                ),
                (
                    sp.expand(matrix_rows[0][0] + hm * matrix_rows[1][0]),
                    sp.expand(matrix_rows[0][1] + hm * matrix_rows[1][1]),
                ),
            )
            positive_scalars = (hm * scalar_parts[0], hm * scalar_parts[1])
            for actual_row, scalar in zip(positive_rows, positive_scalars):
                actual = homogeneous_part(curl(actual_row), degree + left_degree + 1)
                expected = sp.expand(directional(scalar))
                gate(sp.expand(actual - expected) == 0, f"{label} positive-left curl")
                gate(expected != 0, f"{label} positive-left generic lead vanished")

        semantic_rows.append(
            f"mixed-{degree}-{label}:" + hashlib.sha256(
                "|".join(
                    sp.srepr(directional(scalar)) for scalar in constant_scalars
                ).encode()
            ).hexdigest()
        )

# Reconstruct both alternating left words on a nonlinear hostile pair and
# verify the exposed row and preserved determinant directly.
u_probe = x**3 + x * y + 2 * y + 1
w_probe = x**2 * y - y**3 + x - 2
f_probe = x**3 - 2 * x * y + y + 1
g_probe = x**2 * y + x - y**2
for R in (eplus(w_probe) * eminus(u_probe), eminus(u_probe) * eplus(w_probe)):
    N = sp.expand(C * R)
    left_a = sp.expand(eplus(f_probe) * eminus(g_probe) * N)
    left_b = sp.expand(eminus(g_probe) * eplus(f_probe) * N)
    gate(row(left_a, 1) == (
        sp.expand(N[1, 0] + g_probe * N[0, 0]),
        sp.expand(N[1, 1] + g_probe * N[0, 1]),
    ), "left order A exposed row")
    gate(row(left_b, 0) == (
        sp.expand(N[0, 0] + f_probe * N[1, 0]),
        sp.expand(N[0, 1] + f_probe * N[1, 1]),
    ), "left order B exposed row")
    gate(sp.expand(left_a.det()) == 1, "decorated determinant A")
    gate(sp.expand(left_b.det()) == 1, "decorated determinant B")

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

print("theorem=THM-3709-Cohn-alternating-two-by-two-decoration-nonentry")
print("right_orders=E+(w)E-(u),E-(u)E+(w);left_orders=E+(f)E-(g),E-(g)E+(f)")
print("right_parameters=nonzero_not_both_constant;left_parameters=arbitrary_polynomial")
print("leading_obstructions=ordinary_or_directional_derivative_with_forced_linear_factor")
print("hostile_degree_grid=both_nonconstant:1..4x1..4;mixed_nonzero_constant:1..4;deg(left_lead):0..4;four_orders=PASS")
print("decorated_determinants=1;exposed_row_identities=PASS")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
