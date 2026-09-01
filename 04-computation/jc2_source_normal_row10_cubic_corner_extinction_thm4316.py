#!/usr/bin/env python3
"""Primary exact SymPy certificate for THM-4316.

This certificate starts from the audited THM-4308 row-eight machinery and
THM-4315 row-nine gate.  It expands the literal residual row G10, transports
the unique row-eight tangent through G9, derives the m=10 Student--Stein
cokernel, and reduces its obstruction on the cubic-corner k=1 row-nine locus.

The conclusion is deliberately local: the fixed normalized residual-weight
at-most-twelve cubic-corner k=1 lane has no degree-capped bracket lift through
row ten.  No row-ten depth iff is needed for that exclusion, and no claim is
made about the other M=12 walls, seam entry, JC(2), or DC(2).
"""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_student_stein_row9_thm4315 as R9  # noqa: E402


R8 = R9.R8
CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def is_zero(expression: sp.Expr) -> bool:
    return sp.cancel(sp.together(sp.expand(expression))) == 0


def require_zero(expression: sp.Expr, label: str) -> None:
    require(is_zero(expression), label)


x = R8.x
Phi = R8.Phi
X = R9.X


ROW10_COKERNEL = [
    46189,
    0,
    14586,
    0,
    15444,
    0,
    30888,
    0,
    99792,
    0,
    489888,
]


P10_EXPECTED = sp.Poly(
    68842335386673891964107421875 * X**3
    - 114199708075156991490870528000000 * X**2
    - 49516799750570385919467992383488000 * X
    - 2955996966894649715961849999793382752256,
    X,
)


def corner_substitutions() -> tuple[dict[sp.Symbol, sp.Expr], sp.Expr, sp.Expr, sp.Expr]:
    """Return the cubic-corner response after solving the row-nine gate."""

    xi_corner = (4343625 * Phi**2 + 124805668864) / 12798000
    eta_corner = (
        sp.Rational(2091705253888, 107983125)
        - sp.Rational(2839, 1185) * Phi**2
    ) / Phi
    alpha_corner = (
        6971519208442078125 * Phi**4
        - 14082869793796263936000 * Phi**2
        - 74378924775425263164981248
    ) / (396452079682031250 * Phi**3)
    substitutions = {
        R8.xi10: xi_corner,
        R8.eta: eta_corner,
        R8.alpha11: alpha_corner,
        R8.beta11: -alpha_corner,
    }
    require_zero(R9.E9_GATE.subs(substitutions), "corner response solves E9")
    return substitutions, xi_corner, eta_corner, alpha_corner


def derive_row_ten_obstruction() -> tuple[sp.Expr, sp.Poly, sp.Poly]:
    """Derive G10, the Student obstruction, Q, and the row-ten cubic."""

    arows, crows, bracket_subs = R8.build_bracket_rows()
    amatrix8, cmatrix8, theta8_symbols, terminal_subs = R8.hasse_audit(
        arows, crows, bracket_subs
    )
    require(
        (amatrix8.rows, amatrix8.cols, amatrix8.rank()) == (63, 131, 51),
        "inherited row-eight P2 projection",
    )
    require(
        (cmatrix8.rows, cmatrix8.cols, cmatrix8.rank()) == (72, 204, 63),
        "inherited row-eight P3 projection",
    )

    # Expand the literal source rows.  G9 is used to select the unique point
    # in the row-eight terminal cloud which can have a row-nine successor.
    g9 = sp.expand(R8.tcoeff(R8.H, 9))
    expected_g9 = (
        (20 * R8.U + 10 * R8.W + 4 * R8.Z) * x**6
        + (10 * R8.alpha11 + 6 * R8.beta11) * x**7
        + (5 * R8.upsilon5 + 4 * R8.xi10) * x**8
        + (R8.eta + R8.zeta3) * x**9
    )
    require_zero(g9 - expected_g9, "literal source G9")

    g10 = sp.expand(R8.tcoeff(R8.H, 10))
    expected_g10 = (
        (15 * R8.U + 10 * R8.W + 6 * R8.Z) * x**8
        + (5 * R8.alpha11 + 4 * R8.beta11) * x**9
        + (R8.upsilon5 + R8.xi10) * x**10
    )
    require_zero(g10 - expected_g10, "literal source G10")

    fixed_a, fixed_c, _, _ = R9.solve_row_eight_for_g9(
        arows,
        crows,
        bracket_subs,
        terminal_subs,
        theta8_symbols,
        g9,
    )
    corner_subs, _, _, alpha_corner = corner_substitutions()
    fixed_a = [
        sp.expand(row.subs(bracket_subs).subs(terminal_subs).subs(corner_subs))
        for row in fixed_a
    ]
    fixed_c = [
        sp.expand(row.subs(bracket_subs).subs(terminal_subs).subs(corner_subs))
        for row in fixed_c
    ]
    g9_corner = sp.expand(
        g9.subs(bracket_subs).subs(terminal_subs).subs(corner_subs)
    )
    require_zero(
        R8.predicted_G(9, fixed_a, fixed_c) - g9_corner,
        "unique row-eight transport matches G9",
    )

    # Solve the row-nine determinant equation.  Its ten-dimensional tangent
    # is deliberately retained long enough to verify that the Student row
    # kills its entire effect on the next compatibility.
    b9 = R8.B_row(9, fixed_a, fixed_c)
    a9base, c9base = R8.particular_row(9, b9)
    theta9_symbols = list(sp.symbols("theta9_0:10"))
    theta9 = sum(theta9_symbols[j] * x**j for j in range(10))
    a9 = sp.expand(a9base + theta9 * sp.diff(R8.A0, x))
    c9 = sp.expand(c9base + theta9 * sp.diff(R8.C0, x))

    base_a = fixed_a + [a9base]
    base_c = fixed_c + [c9base]
    tangent_a = fixed_a + [a9]
    tangent_c = fixed_c + [c9]
    tangent_change = sp.expand(
        R8.predicted_G(10, tangent_a, tangent_c)
        - R8.predicted_G(10, base_a, base_c)
    )
    expected_tangent_change = sp.expand(
        -x * theta9 + (x**2 + 6) * sp.diff(theta9, x) / 20
    )
    require_zero(
        tangent_change - expected_tangent_change,
        "row-nine tangent gives m=10 Student operator",
    )

    matrix10 = R9.tangent_matrix(10)
    student10 = R9.primitive_student_row(10)
    require(student10 == ROW10_COKERNEL, "primitive m=10 Student cokernel")
    require(
        sp.Matrix([student10]) * matrix10 == sp.zeros(1, 10),
        "m=10 Student row annihilates tangent image",
    )

    g10_corner = sp.expand(
        g10.subs(bracket_subs).subs(terminal_subs).subs(corner_subs)
    )
    difference10 = sp.expand(g10_corner - R8.predicted_G(10, base_a, base_c))
    obstruction10 = sp.factor(
        sum(student10[j] * R8.xcoeff(difference10, j) for j in range(11))
    )
    expected_obstruction = sp.Rational(
        143, 2518243204518514160156250
    ) * P10_EXPECTED.as_expr().subs(X, Phi**2) / Phi**4
    require_zero(obstruction10 - expected_obstruction, "exact row-ten obstruction")

    primitive_content, primitive_p10 = P10_EXPECTED.primitive()
    require(primitive_content == 1, "P10 primitive content")
    require(primitive_p10.LC() > 0 and primitive_p10.degree() == 3, "P10 normalization")

    # Reconstruct rather than merely import the row-nine k=1 equation Q.
    u_corner = -sp.Rational(13, 57591000) * (
        820125 * Phi**2 + 13056802816
    )
    c2_corner = sp.Rational(11, 474000) * (
        14625 * Phi**2 + 404652032
    )
    square_relation = sp.together(
        (alpha_corner**2 - 4 * u_corner * c2_corner) * Phi**2
    )
    raw_q = sp.Poly(square_relation.as_numer_denom()[0], Phi)
    _, primitive_q_phi = raw_q.primitive()
    quintic = R9.even_polynomial_in_x(primitive_q_phi.as_expr())
    require(quintic == R9.Q_EXPECTED, "reconstructed row-nine quintic Q")
    require(sp.gcd(quintic, P10_EXPECTED) == sp.Poly(1, X), "exact gcd(Q,P10)=1")

    # Compact independent modular witness for the exact gcd.
    qmod7 = sp.Poly(quintic.as_expr(), X, modulus=7)
    pmod7 = sp.Poly(P10_EXPECTED.as_expr(), X, modulus=7)
    expected_qmod7 = sp.Poly(
        3 * X**5 + 2 * X**4 + X**3 - 2 * X**2 - 3,
        X,
        modulus=7,
    )
    expected_pmod7 = sp.Poly(3 * X**3 - X**2 + 1, X, modulus=7)
    require(qmod7 == expected_qmod7, "Q reduction modulo 7")
    require(pmod7 == expected_pmod7, "P10 reduction modulo 7")
    bezout7 = sp.Poly(
        (X**2 + 3 * X + 3) * qmod7.as_expr()
        + (-X**4 + 3 * X**3 - 2 * X**2 + 2 * X + 3)
        * pmod7.as_expr()
        - 1,
        X,
        modulus=7,
    )
    require(bezout7.is_zero, "mod-7 Bezout certificate")

    # The row-nine mod-19 positive control is killed at row ten as expected.
    require(P10_EXPECTED.eval(7) % 19 == 15, "mod-19 row-nine control fails P10")

    return g10, obstruction10, quintic


def main() -> None:
    g10, obstruction10, quintic = derive_row_ten_obstruction()
    print("THM-4316 SOURCE-NORMAL ROW-TEN CUBIC-CORNER EXTINCTION: PASS")
    print(f"G10={sp.sstr(g10)}")
    print(f"COKERNEL_m10={ROW10_COKERNEL}")
    print("STUDENT_m10=((x^2+6)*theta'-20*x*theta)/20; rank=10 cokernel=1")
    print(f"P10={sp.sstr(P10_EXPECTED.as_expr())}")
    print(
        "OBSTRUCTION10=143*P10(Phi^2)/(2518243204518514160156250*Phi^4)"
    )
    print(f"ROW9_Q={sp.sstr(quintic.as_expr())}")
    print("MOD7 Q=3X^5+2X^4+X^3-2X^2-3; P10=3X^3-X^2+1")
    print(
        "MOD7_BEZOUT (X^2+3X+3)Q+(-X^4+3X^3-2X^2+2X+3)P10=1"
    )
    print("CONSEQUENCE gcd(Q,P10)=1; every cubic-corner k1 row-nine point dies before row ten")
    print("DEPTH row-ten depth iff not needed; bracket compatibility is prior necessary")
    print("MOD19_CONTROL X=7 survives row9 and has P10(7)=15!=0")
    print(f"CHECKS={CHECKS}")
    print(
        "SCOPE fixed normalized characteristic-zero residual_weight<=12 cubic-corner k1; "
        "no other M12 wall, seam entry, JC2, or DC2 conclusion"
    )


if __name__ == "__main__":
    main()
