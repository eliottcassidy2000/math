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

# All-constant nonzero right parameters reduce to a general constant
# R=[[aa,bb],[rrc,dd0]] in SL_2 with bb,rrc nonzero.  The only possible
# constant exposed-row repair is a Broughton-type resonance.
aa, bb, rrc, dd0, gg, ff = sp.symbols(
    "aa bb rrc dd0 gg ff", nonzero=True
)
R0 = sp.Matrix(((aa, bb), (rrc, dd0)))
N0 = sp.expand(C * R0)
rows0 = (row(N0, 0), row(N0, 1))
ell1 = aa * y + rrc * x
ell2 = bb * y + dd0 * x
gate(
    all(
        sp.expand(actual - expected) == 0
        for actual, expected in zip(
            (homogeneous_part(entry, 2) for entry in rows0[0]),
            (x * ell1, x * ell2),
        )
    ),
    "constant R top row",
)
gate(
    all(
        sp.expand(actual - expected) == 0
        for actual, expected in zip(
            (homogeneous_part(entry, 2) for entry in rows0[1]),
            (-y * ell1, -y * ell2),
        )
    ),
    "constant R bottom row",
)
Y_lower = (aa - 2 * dd0) * x - bb * y
Y_upper = rrc * x + (2 * aa - dd0) * y
gate(sp.expand(curl(rows0[0]) - Y_lower) == 0, "constant R top curl")
gate(sp.expand(curl(rows0[1]) + Y_upper) == 0, "constant R bottom curl")

lower_constant_row = (
    sp.expand(rows0[1][0] + gg * rows0[0][0]),
    sp.expand(rows0[1][1] + gg * rows0[0][1]),
)
upper_constant_row = (
    sp.expand(rows0[0][0] + ff * rows0[1][0]),
    sp.expand(rows0[0][1] + ff * rows0[1][1]),
)
gate(
    sp.expand(
        homogeneous_part(curl(lower_constant_row), 1)
        - ((gg * (aa - 2 * dd0) - rrc) * x
           + (dd0 - 2 * aa - gg * bb) * y)
    ) == 0,
    "lower constant closure equations",
)
gate(
    sp.expand(
        homogeneous_part(curl(upper_constant_row), 1)
        - ((aa - 2 * dd0 - ff * rrc) * x
           + (ff * (dd0 - 2 * aa) - bb) * y)
    ) == 0,
    "upper constant closure equations",
)

# Check the two slope equations used to exclude every positive-degree
# exposed coefficient when bb,rrc are nonzero.
tt = sp.symbols("tt")
for degree in range(1, 7):
    pcoeff = sp.symbols(f"ptop_{degree}_0:{degree + 1}")
    ptt = sum(pcoeff[index] * tt**index for index in range(degree + 1))
    lower_slope = sp.expand(
        (bb * tt**2 + (aa + dd0) * tt + rrc) * sp.diff(ptt, tt)
        - ((degree + 1) * bb * tt + (degree + 2) * dd0 - aa) * ptt
    )
    upper_slope = sp.expand(
        (rrc * tt**2 + (aa + dd0) * tt + bb) * sp.diff(ptt, tt)
        - ((degree + 1) * rrc * tt + (degree + 2) * aa - dd0) * ptt
    )
    gate(
        sp.expand(sp.Poly(lower_slope, tt).coeff_monomial(tt ** (degree + 1))
                  + bb * pcoeff[degree]) == 0,
        "lower slope top coefficient",
    )
    gate(
        sp.expand(sp.Poly(upper_slope, tt).coeff_monomial(tt ** (degree + 1))
                  + rrc * pcoeff[degree]) == 0,
        "upper slope top coefficient",
    )

s0 = aa - dd0

# Lower exposed-row resonance: substitute its two closure equations.
lower_sub = {
    bb: (dd0 - 2 * aa) / gg,
    rrc: gg * (aa - 2 * dd0),
}
R_lower = sp.simplify(R0.subs(lower_sub))
N_lower = sp.simplify(C * R_lower)
lower_rows = (row(N_lower, 0), row(N_lower, 1))
beta_lower = (
    sp.expand(lower_rows[1][0] + gg * lower_rows[0][0]),
    sp.expand(lower_rows[1][1] + gg * lower_rows[0][1]),
)
X_lower = gg * x - y
Y_lower_res = (aa - 2 * dd0) * x - lower_sub[bb] * y
Q_lower = sp.expand(2 * s0 * X_lower + X_lower**2 * Y_lower_res / 3)
gate(sp.simplify(R_lower.det() - 2 * s0**2) == 0, "lower resonance determinant")
gate(sp.expand(sp.diff(X_lower, x) * sp.diff(Y_lower_res, y)
               - sp.diff(X_lower, y) * sp.diff(Y_lower_res, x) - 3 * s0) == 0,
     "lower resonance coordinates")
gate(sp.expand(beta_lower[0] - sp.diff(Q_lower, x)) == 0,
     "lower Broughton potential x")
gate(sp.expand(beta_lower[1] - sp.diff(Q_lower, y)) == 0,
     "lower Broughton potential y")
gate(sp.expand(curl(lower_rows[0]) - Y_lower_res) == 0,
     "lower complementary curl")

# Upper exposed-row resonance is the source/target-dual calculation.
upper_sub = {
    rrc: (aa - 2 * dd0) / ff,
    bb: ff * (dd0 - 2 * aa),
}
R_upper = sp.simplify(R0.subs(upper_sub))
N_upper = sp.simplify(C * R_upper)
upper_rows = (row(N_upper, 0), row(N_upper, 1))
beta_upper = (
    sp.expand(upper_rows[0][0] + ff * upper_rows[1][0]),
    sp.expand(upper_rows[0][1] + ff * upper_rows[1][1]),
)
X_upper = x - ff * y
Y_upper_res = upper_sub[rrc] * x + (2 * aa - dd0) * y
Q_upper = sp.expand(2 * s0 * X_upper + X_upper**2 * Y_upper_res / 3)
gate(sp.simplify(R_upper.det() - 2 * s0**2) == 0, "upper resonance determinant")
gate(sp.expand(sp.diff(X_upper, x) * sp.diff(Y_upper_res, y)
               - sp.diff(X_upper, y) * sp.diff(Y_upper_res, x) - 3 * s0) == 0,
     "upper resonance coordinates")
gate(sp.expand(beta_upper[0] - sp.diff(Q_upper, x)) == 0,
     "upper Broughton potential x")
gate(sp.expand(beta_upper[1] - sp.diff(Q_upper, y)) == 0,
     "upper Broughton potential y")
gate(sp.expand(curl(upper_rows[1]) + Y_upper_res) == 0,
     "upper complementary curl")

# In normalized Broughton coordinates Q/(2s)=X+X^2T.  Its Hamiltonian
# derivation kills every Q-power.  The homogeneous cubic kernel is exactly
# the span of (X^2T)^j in each degree, which powers the degree descent.
XX, TT = sp.symbols("XX TT")
q_broughton = XX + XX**2 * TT


def broughton_derivation(expr: sp.Expr) -> sp.Expr:
    return sp.expand(
        XX**2 * sp.diff(expr, XX)
        - (1 + 2 * XX * TT) * sp.diff(expr, TT)
    )


for power in range(0, 8):
    gate(broughton_derivation(q_broughton**power) == 0,
         "Broughton kernel power")
for total_degree in range(0, 16):
    for j in range(total_degree + 1):
        i = total_degree - j
        monomial = XX**i * TT**j
        cubic_bracket = sp.expand(
            sp.diff(monomial, XX) * sp.diff(XX**2 * TT, TT)
            - sp.diff(monomial, TT) * sp.diff(XX**2 * TT, XX)
        )
        gate(
            sp.expand(cubic_bracket - (i - 2 * j) * XX ** (i + 1) * TT**j) == 0,
            "homogeneous cubic kernel weight",
        )
        gate((cubic_bracket == 0) == (i == 2 * j),
             "homogeneous cubic kernel support")

semantic_rows.append(
    "constant-resonance:" + hashlib.sha256(
        "|".join(
            sp.srepr(value)
            for value in (
                Q_lower,
                Q_upper,
                Y_lower_res,
                Y_upper_res,
                broughton_derivation(6 * TT),
            )
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
print("right_parameters=both_nonzero;left_parameters=arbitrary_polynomial")
print("leading_obstructions=ordinary_or_directional_derivative_with_forced_linear_factor")
print("constant_right_resonance=broughton_X_plus_X2T;hamiltonian_cokernel=PASS")
print("hostile_degree_grid=both_nonconstant:1..4x1..4;mixed:1..4;constant_slope:1..6;hamiltonian:0..15;four_orders=PASS")
print("decorated_determinants=1;exposed_row_identities=PASS")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
