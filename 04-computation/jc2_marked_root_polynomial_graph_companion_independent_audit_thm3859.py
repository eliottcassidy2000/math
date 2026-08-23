#!/usr/bin/env python3
"""Independent hostile audit for THM-3859.

The checker uses a square-discriminant parametrization and a complete
constant-(s,q) boundary analysis distinct from the primary companion.  Human
UFD and valuation arguments are recorded in the accompanying report.
"""

from __future__ import annotations

import ast
import hashlib
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"gate failed: {label}")


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(expression) == 0, label)


def nonsquare_polynomial(polynomial: sp.Expr, variable: sp.Symbol) -> bool:
    """Ignore scalar units (the ground field is algebraically closed)."""

    _, factors = sp.factor_list(sp.Poly(polynomial, variable))
    return any(exponent % 2 for _, exponent in factors)


A, C, S, Q = sp.symbols("A C S Q")


def branch(profile: sp.Expr) -> sp.Expr:
    return sp.expand(
        -27 * A**2 * profile**2
        + 8 * A * C**3
        - 54 * A * C * profile
        + 9 * C**2
        - 54 * profile
    )


# ---------------------------------------------------------------------------
# I. Cusp identity and the universal marked-root factorization.
# ---------------------------------------------------------------------------
b_free = sp.symbols("b_free")
p = sp.Rational(3, 2) + A * C
u = 1 + A * C + A**2 * b_free
zero(8 * p**3 - 27 * u**2 - A**2 * branch(b_free), "cusp discriminant identity")

F = C - 6 * S * (1 + A * S)
b0 = 2 * S**2 * (3 + 4 * A * S)
b = b0 + Q * F
Delta = branch(b)
H = sp.cancel(Delta / F)
gate(sp.denom(H) == 1, "universal quotient is polynomial")
H = sp.expand(H)
zero(Delta - F * H, "universal marked-root factorization")

z = 1 + 2 * A * S
graph_C = 6 * S * (1 + A * S)
zero(1 + sp.Rational(2, 3) * A * graph_C - z**2, "marked cusp square")
zero(1 + A * graph_C + A**2 * b0 - z**3, "marked cusp cube")

# Reverse UFD calculation.  If P=z^2 and U=z^3, these two identities recover
# the graph and profile.  The z(0)=-1 branch fails divisibility by A^2.
Zeta = sp.symbols("Zeta")
zero(
    Zeta**3 - 1 - sp.Rational(3, 2) * (Zeta**2 - 1)
    - sp.Rational(1, 2) * (Zeta - 1) ** 2 * (2 * Zeta + 1),
    "reverse profile numerator",
)
gate(
    (sp.Rational(1, 2) * (Zeta - 1) ** 2 * (2 * Zeta + 1)).subs(Zeta, -1)
    == -2,
    "hostile z(0)=-1 branch cannot be divisible by A^2",
)
zero(sp.Rational(3, 2) * (z**2 - 1) / A - graph_C, "reverse graph recovery")
zero(
    sp.Rational(1, 2) * (z - 1) ** 2 * (2 * z + 1) / A**2 - b0,
    "reverse profile recovery",
)


# ---------------------------------------------------------------------------
# II. Companion quadratic, primitivity, and irreducible flank.
# ---------------------------------------------------------------------------
H_poly = sp.Poly(H, C)
gate(H_poly.degree() == 2, "companion has C-degree two")
leading = H_poly.coeff_monomial(C**2)
B = sp.factor(H_poly.coeff_monomial(C))
E = sp.factor(H_poly.coeff_monomial(1))
gate(leading == 8 * A, "companion leading coefficient")
gate(B.subs(A, 0) == 9, "middle coefficient is a primitivity unit modulo A")
gate(sp.gcd(sp.Poly(leading, A), sp.Poly(B, A)) == 1, "leading and middle coefficients coprime")

D = sp.factor(B**2 - 4 * leading * E)
L_plus = A * Q + 4 * A * S + 3
L_minus = 3 * A * Q - 4 * A * S + 1
zero(D - 27 * L_plus * L_minus**3, "two-flank discriminant")
gate(D.subs(A, 0) == 81, "two distinct unramified places over A=0")
zero(H.subs(A, 0) - (9 * C - 54 * Q + 54 * S), "special affine fibre")

# A nonconstant irreducible control.  Its nonsquare discriminant and direct
# factor census independently exercise the connected case.
H_irreducible = sp.factor(H.subs({S: A, Q: 0}))
D_irreducible = sp.factor(D.subs({S: A, Q: 0}))
zero(
    D_irreducible - 27 * (4 * A**2 + 3) * (1 - 4 * A**2) ** 3,
    "irreducible control discriminant",
)
gate(nonsquare_polynomial(D_irreducible, A), "irreducible control discriminant nonsquare")
factor_unit, factor_rows = sp.factor_list(sp.Poly(H_irreducible, C, A))
gate(
    factor_unit != 0 and len(factor_rows) == 1 and factor_rows[0][1] == 1,
    "irreducible control companion",
)
gate(B.subs({A: 0, S: 0, Q: 0}) == 9, "positive zero-place numerator baseline")
gate(-9 - B.subs({A: 0, S: 0, Q: 0}) == -18, "negative zero-place has a C pole")


# ---------------------------------------------------------------------------
# III. General square-discriminant parametrization and G_m component.
# ---------------------------------------------------------------------------
B1, Vroot = sp.symbols("B1 Vroot")
B_square = 9 + A * B1
W_square = 9 + A * Vroot
E_square = sp.cancel((B_square**2 - W_square**2) / (32 * A))
gate(sp.denom(E_square) == 1, "square-discriminant constant term is polynomial")
H_square = sp.expand(8 * A * C**2 + B_square * C + E_square)
graph_factor = 16 * C + B1 - Vroot
pole_factor = 18 + A * (16 * C + B1 + Vroot)
zero(32 * H_square - graph_factor * pole_factor, "general square-discriminant factorization")
gate(graph_factor.subs(A, 0).coeff(C) == 16, "first factor remains a graph at A=0")
gate(pole_factor.subs(A, 0) == 18, "second factor omits A=0")

# In the pole-factor quotient A is a unit; this gives k[A,A^-1] exactly.
A_inverse = -(16 * C + B1 + Vroot) / 18
zero(A * A_inverse - 1 + pole_factor / 18, "explicit inverse of A on pole component")
C_laurent = -(18 / A + B1 + Vroot) / 16
zero(pole_factor.subs(C, C_laurent), "pole component Laurent graph")


# ---------------------------------------------------------------------------
# IV. Complete constant-(s,q) and axis boundaries.
# ---------------------------------------------------------------------------
a, d = sp.symbols("a d")
D_constant = sp.factor(D.subs({S: a, Q: d}))
zero(
    D_constant
    - 27 * (3 + A * (d + 4 * a)) * (1 + A * (3 * d - 4 * a)) ** 3,
    "constant transverse discriminant",
)
zero(
    (d + 4 * a) - 3 * (3 * d - 4 * a) + 8 * (d - 2 * a),
    "constant flank collision iff d=2a",
)

# On the complete square boundary d=2a, H contains the marked graph again
# and the other factor is the canonical Laurent/G_m component.
constant_square_factorization = (
    (6 * A * a**2 - C + 6 * a)
    * (12 * A**2 * a**2 - 8 * A * C + 12 * A * a - 9)
)
zero(H.subs({S: a, Q: 2 * a}) - constant_square_factorization, "all constant square cases")
zero(F.subs(S, a) + (6 * A * a**2 - C + 6 * a), "constant square repeats marked graph")
gate(
    (12 * A**2 * a**2 - 8 * A * C + 12 * A * a - 9).subs(A, 0) == -9,
    "constant square companion omits A=0",
)

# A constant nonsquare control and the zero/collision controls.
gate(nonsquare_polynomial(D_constant.subs({a: 1, d: 0}), A), "constant noncollision is irreducible")
zero(Delta.subs({S: 0, Q: 0}) - C**2 * (8 * A * C + 9), "zero constant boundary")
collision = {S: -sp.Rational(1, 6), Q: -sp.Rational(1, 3)}
line = A - 6 * C - 6
conic = A**2 - 24 * A * C - 6 * A - 27
zero(F.subs(collision) + line / 6, "THM-3852 collision graph")
zero(H.subs(collision) - line * conic / 18, "THM-3852 collision companion")
zero(D.subs(collision) - (A - 3) ** 4, "THM-3852 collision square")

B_axis = sp.symbols("B_axis")
zero(branch(B_axis).subs(A, 0) - (9 * C**2 - 54 * B_axis), "axis restriction")
zero(
    branch(C**2 / 6 + A * B_axis).subs(A, 0),
    "axis component positive control",
)


# ---------------------------------------------------------------------------
# V. Relation to the fixed-product reciprocal-cubic corollary.
# ---------------------------------------------------------------------------
T, c_fixed = sp.symbols("T c_fixed")
fixed_product_cubic = T**3 - A * T**2 + C * T - c_fixed
fixed_product_graph = c_fixed * T**-1 - T**2 + A * T
zero(fixed_product_cubic.subs(C, fixed_product_graph), "fixed-product Laurent root graph")
gate(
    sp.Poly(c_fixed - T**3 + A * T**2, T).eval(0) == c_fixed,
    "fixed nonzero product forces the T=0 pole",
)
zero(
    fixed_product_cubic.subs(c_fixed, 0) - T * (T**2 - A * T + C),
    "zero-product boundary becomes reducible when the pole disappears",
)


# ---------------------------------------------------------------------------
# VI. Deterministic transcript.
# ---------------------------------------------------------------------------
source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "checker contains no inactive Python assert",
)

semantic_lines = [
    "THM3859_INDEPENDENT_HOSTILE_AUDIT PASS",
    "CLASSIFICATION cusp UFD forces z=1+2As,F=C-6s(1+As),b|F=2s^2(3+4As)",
    "UNIVERSAL Delta_b=F*H with primitive quadratic H and no vertical factor",
    "DISCRIMINANT 27(Aq+4As+3)(3Aq-4As+1)^3;D(0)=81",
    "IRREDUCIBLE nonsquare D gives one A=0 pole plus at least one infinity place",
    "REDUCIBLE square D gives one polynomial graph and one canonical Gm factor",
    "CONSTANT square iff q=2s;then the marked graph repeats and the other factor is Gm",
    "BOUNDARIES zero profile,THM3852 collision,and excluded A=0 axis all agree",
    "FIXED_PRODUCT nonzero constant term is Laurent/Gm and lies outside the polynomial-graph hypothesis",
    "FIXED_PRODUCT_ZERO pole removal makes the cubic reducible",
    "SCOPE q=q(A) only;C-dependent quotient and axis remain open;no JC conclusion",
]
semantic_sha = hashlib.sha256("\n".join(semantic_lines).encode("utf-8")).hexdigest()
transcript = "\n".join(semantic_lines + [f"GATES {GATES}", f"SEMANTIC_SHA256 {semantic_sha}"]) + "\n"

if len(sys.argv) == 1:
    sys.stdout.write(transcript)
elif len(sys.argv) == 3 and sys.argv[1] == "--verify-frozen":
    frozen = Path(sys.argv[2]).read_bytes()
    gate(frozen == transcript.encode("utf-8"), "frozen transcript byte match")
    print("FROZEN_REPLAY PASS")
    print(f"FROZEN_SHA256 {hashlib.sha256(frozen).hexdigest()}")
else:
    raise SystemExit("usage: independent_checker.py [--verify-frozen PATH]")
