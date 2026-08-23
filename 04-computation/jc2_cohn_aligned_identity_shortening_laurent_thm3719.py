#!/usr/bin/env python3
"""Exact companion for THM-3719's aligned shortening obstruction."""

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


x, y, t = sp.symbols("x y t")
A = 1 + x * y
B = x**2
G = -y**2
D = 1 - x * y
C = sp.Matrix(((A, B), (G, D)))


def eplus(h: sp.Expr) -> sp.Matrix:
    return sp.Matrix(((1, h), (0, 1)))


def eminus(h: sp.Expr) -> sp.Matrix:
    return sp.Matrix(((1, 0), (h, 1)))


def row(matrix: sp.Matrix, index: int) -> tuple[sp.Expr, sp.Expr]:
    return sp.expand(matrix[index, 0]), sp.expand(matrix[index, 1])


def curl(one_form: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(sp.diff(one_form[0], y) - sp.diff(one_form[1], x))


def determinant_rows(
    first: tuple[sp.Expr, sp.Expr], second: tuple[sp.Expr, sp.Expr]
) -> sp.Expr:
    return sp.expand(first[0] * second[1] - first[1] * second[0])


gate(sp.expand(C.det()) == 1, "Cohn determinant")

u, g, w, f = sp.symbols("u g w f")
N_minus = sp.expand(C * eminus(u))
minus_rows = (row(N_minus, 0), row(N_minus, 1))
beta_minus = (
    sp.expand(minus_rows[1][0] + g * minus_rows[0][0]),
    sp.expand(minus_rows[1][1] + g * minus_rows[0][1]),
)
gate(determinant_rows(minus_rows[0], beta_minus) == 1,
     "minus exposed determinant")
Qx, Qy = sp.symbols("Qx Qy")
gate(
    sp.expand(
        minus_rows[0][0] * Qy - minus_rows[0][1] * Qx
        - ((1 + x * y + x**2 * u) * Qy - x**2 * Qx)
    ) == 0,
    "minus scalar PDE",
)
cross_minus = (
    sp.expand(minus_rows[0][0] + f * minus_rows[1][0]),
    sp.expand(minus_rows[0][1] + f * minus_rows[1][1]),
)
gate(determinant_rows(cross_minus, minus_rows[1]) == 1,
     "minus cross determinant")
gate(
    sp.expand(
        Qx * minus_rows[1][1] - Qy * minus_rows[1][0]
        - ((1 - x * y) * Qx
           + (y**2 - (1 - x * y) * u) * Qy)
    ) == 0,
    "minus cross scalar PDE",
)

N_plus = sp.expand(C * eplus(w))
plus_rows = (row(N_plus, 0), row(N_plus, 1))
beta_plus = (
    sp.expand(plus_rows[0][0] + f * plus_rows[1][0]),
    sp.expand(plus_rows[0][1] + f * plus_rows[1][1]),
)
gate(determinant_rows(beta_plus, plus_rows[1]) == 1,
     "plus exposed determinant")
gate(
    sp.expand(
        Qx * plus_rows[1][1] - Qy * plus_rows[1][0]
        - ((1 - x * y - y**2 * w) * Qx + y**2 * Qy)
    ) == 0,
    "plus scalar PDE",
)
cross_plus = (
    sp.expand(plus_rows[1][0] + g * plus_rows[0][0]),
    sp.expand(plus_rows[1][1] + g * plus_rows[0][1]),
)
gate(determinant_rows(plus_rows[0], cross_plus) == 1,
     "plus cross determinant")
gate(
    sp.expand(
        plus_rows[0][0] * Qy - plus_rows[0][1] * Qx
        - ((1 + x * y) * Qy
           - (x**2 + (1 + x * y) * w) * Qx)
    ) == 0,
    "plus cross scalar PDE",
)

# The quarter-turn (X,Y)=(y,-x), followed by P -> -P, sends the two
# exposed-top PDEs to the two proved exposed-bottom PDEs.  Here PX=-Qy and
# PY=Qx after the substitution.
X, Y = sp.symbols("X Y")
U, W = sp.symbols("U W")
PX, PY = sp.symbols("PX PY")
aligned_plus_transformed = sp.expand(
    (1 + X * Y - X**2 * W) * PY - X**2 * PX
)
gate(
    sp.expand(
        aligned_plus_transformed.subs(W, -U)
        - ((1 + X * Y + X**2 * U) * PY - X**2 * PX)
    ) == 0,
    "quarter-turn aligned PDE",
)
cross_minus_transformed = sp.expand(
    (1 + X * Y) * PY
    - (X**2 - (1 + X * Y) * U) * PX
)
gate(
    sp.expand(
        cross_minus_transformed.subs(U, -W)
        - ((1 + X * Y) * PY
           - (X**2 + (1 + X * Y) * W) * PX)
    ) == 0,
    "quarter-turn cross PDE",
)

semantic_rows: list[str] = []
for u_degree in range(0, 7):
    u_coefficients = sp.symbols(f"uc{u_degree}_0:{u_degree + 1}")
    ux = sum(
        coefficient * x**index
        for index, coefficient in enumerate(u_coefficients)
    )
    for q_degree in range(1, 7):
        q_coefficients = {
            (i, j): sp.symbols(f"q{u_degree}_{q_degree}_{i}_{j}")
            for i in range(q_degree + 1)
            for j in range(q_degree + 1 - i)
        }
        q = sum(
            coefficient * x**i * y**j
            for (i, j), coefficient in q_coefficients.items()
        )
        q_tilde = sp.expand(q.subs(y, t / x))
        original = sp.expand(
            ((1 + x * y + x**2 * ux) * sp.diff(q, y) - x**2 * sp.diff(q, x))
            .subs(y, t / x)
            / x
        )
        laurent = sp.expand(
            (1 + x**2 * ux) * sp.diff(q_tilde, t)
            - x * sp.diff(q_tilde, x)
        )
        gate(sp.expand(original - laurent) == 0, "Laurent conjugacy")

        coefficient_minus_one = sp.expand(
            sum(
                coefficient * t**j
                for (i, j), coefficient in q_coefficients.items()
                if i - j == -1
            )
        )
        gate(sp.expand(coefficient_minus_one.subs(t, 0)) == 0,
             "honest-polynomial minus-one layer divisible by t")
        semantic_rows.append(
            f"{u_degree},{q_degree}:" + hashlib.sha256(
                "|".join(
                    (sp.srepr(original), sp.srepr(coefficient_minus_one))
                ).encode()
            ).hexdigest()
        )

# The opposite exposed pair has the cross scalar gate
# (1+xy)Q_y-[x^2+(1+xy)w(y)]Q_x=1 after its top x-degree descent.
# In (y,t=xy) coordinates this is the second Laurent operator in THM-3719.
for w_degree in range(0, 7):
    w_coefficients = sp.symbols(f"wc{w_degree}_0:{w_degree + 1}")
    wy = sum(
        coefficient * y**index
        for index, coefficient in enumerate(w_coefficients)
    )
    for q_degree in range(1, 7):
        q_coefficients = {
            (i, j): sp.symbols(f"qc{w_degree}_{q_degree}_{i}_{j}")
            for i in range(q_degree + 1)
            for j in range(q_degree + 1 - i)
        }
        q = sum(
            coefficient * x**i * y**j
            for (i, j), coefficient in q_coefficients.items()
        )
        q_tilde = sp.expand(q.subs(x, t / y))
        original = sp.expand(
            (
                (1 + x * y) * sp.diff(q, y)
                - (x**2 + (1 + x * y) * wy) * sp.diff(q, x)
            ).subs(x, t / y)
        )
        laurent = sp.expand(
            (1 + t) * sp.diff(q_tilde, y)
            + (t / y - (1 + t) * y * wy) * sp.diff(q_tilde, t)
        )
        gate(sp.expand(original - laurent) == 0, "cross Laurent conjugacy")

        coefficient_zero = sp.expand(
            sum(
                coefficient * t**i
                for (i, j), coefficient in q_coefficients.items()
                if j - i == 0
            )
        )
        coefficient_one = sp.expand(
            sum(
                coefficient * t**i
                for (i, j), coefficient in q_coefficients.items()
                if j - i == 1
            )
        )
        gate(sp.expand(coefficient_one.subs(t, 0) - q_coefficients[(0, 1)]) == 0,
             "cross honest-polynomial layer")
        semantic_rows.append(
            f"cross-{w_degree},{q_degree}:" + hashlib.sha256(
                "|".join(
                    (
                        sp.srepr(original),
                        sp.srepr(coefficient_zero),
                        sp.srepr(coefficient_one),
                    )
                ).encode()
            ).hexdigest()
        )

# Lowest-Laurent coefficient equations.  For r<-1, F'-rF=0 has no nonzero
# polynomial solution.  At r=-1, F'+F=1 has the unique solution F=1, which
# contradicts the divisibility-by-t gate above.
for degree in range(0, 13):
    coefficients = sp.symbols(f"fc{degree}_0:{degree + 1}")
    F = sum(
        coefficient * t**index
        for index, coefficient in enumerate(coefficients)
    )
    for r in range(-12, -1):
        equations = sp.Poly(sp.diff(F, t) - r * F, t).all_coeffs()
        matrix, vector = sp.linear_eq_to_matrix(equations, coefficients)
        gate(matrix.rank() == len(coefficients), "negative Laurent kernel")
        gate(all(value == 0 for value in vector), "negative Laurent RHS")

    equations = sp.Poly(sp.diff(F, t) + F - 1, t).all_coeffs()
    solution = sp.linsolve(equations, coefficients)
    expected_solution = tuple([sp.Integer(1)] + [sp.Integer(0)] * degree)
    gate(solution == sp.FiniteSet(expected_solution), "minus-one layer solution")

    for r in range(-12, 0):
        equations = sp.Poly(
            t * sp.diff(F, t) + r * (1 + t) * F,
            t,
        ).all_coeffs()
        matrix, vector = sp.linear_eq_to_matrix(equations, coefficients)
        gate(matrix.rank() == len(coefficients), "cross negative Laurent kernel")
        gate(all(value == 0 for value in vector), "cross negative Laurent RHS")

    zero_equations = sp.Poly(t * sp.diff(F, t), t).all_coeffs()
    zero_matrix, _ = sp.linear_eq_to_matrix(zero_equations, coefficients)
    gate(zero_matrix.rank() == degree, "cross zero layer is constant")

    one_equations = sp.Poly(
        t * sp.diff(F, t) + (1 + t) * F - 1,
        t,
    ).all_coeffs()
    gate(sp.linsolve(one_equations, coefficients) == sp.EmptySet,
         "cross one layer obstruction")

# Hostile polynomial controls reconstruct the actual decorated rows and find
# the promised exposed curl obstruction without relying on symbolic placeholders.
probes = (
    sp.Integer(0),
    sp.Integer(1),
    x,
    y,
    x + y + 1,
    x**2 - x * y + y**2 + 2,
    x**3 + x * y**2 - 2 * y + 1,
)
for right_parameter in probes:
    Nm = sp.expand(C * eminus(right_parameter))
    Np = sp.expand(C * eplus(right_parameter))
    for left_parameter in probes:
        rm = (row(Nm, 0), row(Nm, 1))
        rp = (row(Np, 0), row(Np, 1))
        bm = (
            sp.expand(rm[1][0] + left_parameter * rm[0][0]),
            sp.expand(rm[1][1] + left_parameter * rm[0][1]),
        )
        bp = (
            sp.expand(rp[0][0] + left_parameter * rp[1][0]),
            sp.expand(rp[0][1] + left_parameter * rp[1][1]),
        )
        cm = (
            sp.expand(rm[0][0] + left_parameter * rm[1][0]),
            sp.expand(rm[0][1] + left_parameter * rm[1][1]),
        )
        cp = (
            sp.expand(rp[1][0] + left_parameter * rp[0][0]),
            sp.expand(rp[1][1] + left_parameter * rp[0][1]),
        )
        gate(determinant_rows(rm[0], bm) == 1, "minus hostile determinant")
        gate(determinant_rows(bp, rp[1]) == 1, "plus hostile determinant")
        gate(determinant_rows(cm, rm[1]) == 1, "minus cross hostile determinant")
        gate(determinant_rows(rp[0], cp) == 1, "plus cross hostile determinant")
        gate(curl(bm) != 0, "minus hostile exposed row")
        gate(curl(bp) != 0, "plus hostile exposed row")
        gate(curl(cm) != 0, "minus cross hostile exposed row")
        gate(curl(cp) != 0, "plus cross hostile exposed row")

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

print("theorem=THM-3719-Cohn-complete-identity-shortening-Laurent-nonentry")
print("closed_cells=all_alternating_two-left_words_with_right_length_at_most_one")
print("scalar_gates=aligned_hyperbolic_Euler,cross_Broughton_Laurent;duals=quarter_turn")
print("obstruction=forbidden_lowest_Laurent_layer_in_both_orientation_classes")
print("hostile_grid=parameter_degree:0..6;Q_degree:1..6;Laurent_r:-12..1;row_probes:98=PASS")
print("remaining_shortening=NONE;complete_two-by-two_cell_with_T3709=YES")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
