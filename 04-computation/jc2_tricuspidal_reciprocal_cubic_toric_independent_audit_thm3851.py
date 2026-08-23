#!/usr/bin/env python3
"""Independent hostile audit of THM-3851's reciprocal-cubic/toric bridge.

This checker deliberately does not reuse the tricuspidal lattice census.  It
audits only the post-promotion bridge from the canonical (A,B,C) quartic to
the reciprocal cubic and its ordered-root torus.  The quotient-ring and unit
conclusions in the accompanying report still use the standard invariant-ring
and Laurent-unit lemmas as human algebraic inputs.
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


def cyclotomic_reduce(expression: sp.Expr, omega: sp.Symbol) -> sp.Expr:
    """Reduce a rational expression modulo omega^2+omega+1 exactly."""

    numerator, denominator = sp.together(expression).as_numer_denom()
    modulus = sp.Poly(omega**2 + omega + 1, omega)
    numerator_reduced = sp.rem(sp.Poly(sp.expand(numerator), omega), modulus).as_expr()
    denominator_reduced = sp.rem(sp.Poly(sp.expand(denominator), omega), modulus).as_expr()
    gate(denominator_reduced != 0, "cyclotomic denominator remains nonzero")
    return sp.cancel(numerator_reduced / denominator_reduced)


# ---------------------------------------------------------------------------
# I. Fourier equivalence with the promoted canonical quartic.
# ---------------------------------------------------------------------------
A, B, C, U, V, Z, omega = sp.symbols("A B C U V Z omega")

delta_three = sp.expand((A * B + A * C + B * C) ** 2 - 4 * A * B * C * (A + B + C))
D = sp.expand(U**2 * V**2 - 4 * (U**3 + V**3) * Z + 18 * U * V * Z**2 - 27 * Z**4)

U_fourier = A + omega * B + omega**2 * C
V_fourier = A + omega**2 * B + omega * C
Z_fourier = (A + B + C) / 3

zero(
    cyclotomic_reduce(
        D.subs({U: U_fourier, V: V_fourier, Z: Z_fourier}) - 9 * delta_three,
        omega,
    ),
    "Fourier quartic identity",
)

forward_matrix = sp.Matrix(
    [
        [1, omega, omega**2],
        [1, omega**2, omega],
        [sp.Rational(1, 3), sp.Rational(1, 3), sp.Rational(1, 3)],
    ]
)
gate(cyclotomic_reduce(forward_matrix.det(), omega) != 0, "Fourier matrix invertible")

A_inverse = Z + (U + V) / 3
B_inverse = Z + (omega**2 * U + omega * V) / 3
C_inverse = Z + (omega * U + omega**2 * V) / 3
for original, recovered, label in (
    (A, A_inverse.subs({U: U_fourier, V: V_fourier, Z: Z_fourier}, simultaneous=True), "A"),
    (B, B_inverse.subs({U: U_fourier, V: V_fourier, Z: Z_fourier}, simultaneous=True), "B"),
    (C, C_inverse.subs({U: U_fourier, V: V_fourier, Z: Z_fourier}, simultaneous=True), "C"),
):
    zero(cyclotomic_reduce(recovered - original, omega), f"Fourier inverse {label}")
zero(3 * Z_fourier - (A + B + C), "canonical bitangent becomes Z=0")

# The report's reciprocal normalization is not merely another quartic model.
# An explicit Mobius change identifies it with the Fourier transform of the
# promoted normalization.  On the canonical t=1 chart put p=s/t and
# q=omega^2(p+omega)/(p+omega^2).
p, q = sp.symbols("p q")
A_p = p**2
B_p = (p - 1) ** 2
C_p = p**2 * (p - 1) ** 2
U_p = U_fourier.subs({A: A_p, B: B_p, C: C_p})
V_p = V_fourier.subs({A: A_p, B: B_p, C: C_p})
Z_p = Z_fourier.subs({A: A_p, B: B_p, C: C_p})
q_p = omega**2 * (p + omega) / (p + omega**2)
zero(
    cyclotomic_reduce(U_p / Z_p - (2 * q_p + q_p**-2), omega),
    "Fourier normalization U coordinate",
)
zero(
    cyclotomic_reduce(V_p / Z_p - (q_p**2 + 2 * q_p**-1), omega),
    "Fourier normalization V coordinate",
)
for address, expected, label in (
    (0, omega, "canonical cusp p=0"),
    (1, 1, "canonical cusp p=1"),
):
    zero(cyclotomic_reduce(q_p.subs(p, address) - expected, omega), label)
zero(cyclotomic_reduce(sp.limit(q_p, p, sp.oo) - omega**2, omega), "canonical cusp p=infinity")


# ---------------------------------------------------------------------------
# II. Reciprocal cubic, root sheet, branch normalization, and places.
# ---------------------------------------------------------------------------
T = sp.symbols("T")
f = T**3 - U * T**2 + V * T - 1
f_prime = sp.diff(f, T)

# Independent discriminant path: use the Sylvester determinant, not
# sympy.discriminant.  For a monic cubic disc(f)=-Res(f,f').
sylvester = sp.Matrix(
    [
        [1, -U, V, -1, 0],
        [0, 1, -U, V, -1],
        [3, -2 * U, V, 0, 0],
        [0, 3, -2 * U, V, 0],
        [0, 0, 3, -2 * U, V],
    ]
)
zero(-sylvester.det() - D.subs(Z, 1), "reciprocal cubic discriminant by Sylvester determinant")

# In the quotient by f, T(T^2-UT+V)=1.  Solving the relation for V gives
# the mutually inverse Laurent presentation k[U,T,T^-1].
T_inverse_in_quotient = T**2 - U * T + V
zero(T * T_inverse_in_quotient - 1 - f, "T is a unit in the cubic algebra")
V_laurent = T**-1 - T**2 + U * T
zero(f.subs(V, V_laurent), "Laurent presentation annihilates the cubic relation")
zero(T_inverse_in_quotient - T**2 + U * T - V, "Laurent and quotient presentations compose")

# The source map (U,T)->(U,V) has the same ramification divisor as the
# monogenic different, up to the already-invertible factor -T.
source_jacobian = sp.diff(V_laurent, T)
different_on_source = f_prime.subs(V, V_laurent)
zero(source_jacobian + different_on_source / T, "Jacobian equals different up to a unit")

# A repeated root q gives roots (q,q,q^-2), hence the exact normalization of
# the discriminant branch.  This also exposes the three triple-root cusps.
U_q = 2 * q + q**-2
V_q = q**2 + 2 * q**-1
zero(
    f.subs({U: U_q, V: V_q}) - (T - q) ** 2 * (T - q**-2),
    "double-root factorization",
)
zero(D.subs({U: U_q, V: V_q, Z: 1}), "branch normalization lies on D")
zero(sp.diff(U_q, q) - 2 * (1 - q**-3), "normalization U derivative")
zero(sp.diff(V_q, q) - 2 * q * (1 - q**-3), "normalization V derivative")
derivative_gcd = sp.monic(
    sp.gcd(
        sp.together(sp.diff(U_q, q)).as_numer_denom()[0],
        sp.together(sp.diff(V_q, q)).as_numer_denom()[0],
    ),
    q,
)
gate(derivative_gcd == q**3 - 1, "only the three triple roots ramify the normalization")

# Hostile injectivity test.  If nonzero q,r have the same U and V and q!=r,
# their difference equations reduce to the following saturated unit ideal.
r, inverse = sp.symbols("r inverse")
U_difference_numerator = 2 * q**2 * r**2 - q - r
V_difference_numerator = q * r * (q + r) - 2
injectivity_groebner = sp.groebner(
    [
        U_difference_numerator,
        V_difference_numerator,
        inverse * q * r * (q - r) - 1,
    ],
    inverse,
    q,
    r,
    order="grevlex",
)
gate(
    len(injectivity_groebner.polys) == 1
    and injectivity_groebner.polys[0].as_expr() == 1,
    "branch normalization is injective",
)

# Homogenization shows that q=0 and q=infinity are the only missing places,
# and they land at distinct smooth points of the projective quartic.
s, t = sp.symbols("s t")
U_st = t * (2 * s**3 + t**3)
V_st = s * (s**3 + 2 * t**3)
Z_st = s**2 * t**2
zero(D.subs({U: U_st, V: V_st, Z: Z_st}), "projective reciprocal normalization")
gate((U_st.subs({s: 0, t: 1}), V_st.subs({s: 0, t: 1}), Z_st.subs({s: 0, t: 1})) == (1, 0, 0), "q=0 place")
gate((U_st.subs({s: 1, t: 0}), V_st.subs({s: 1, t: 0}), Z_st.subs({s: 1, t: 0})) == (0, 1, 0), "q=infinity place")
gate(sp.diff(D, Z).subs({U: 1, V: 0, Z: 0}) == -4, "q=0 image smooth")
gate(sp.diff(D, Z).subs({U: 0, V: 1, Z: 0}) == -4, "q=infinity image smooth")


# ---------------------------------------------------------------------------
# III. Ordered-root torus and quotient dictionary.
# ---------------------------------------------------------------------------
x, y = sp.symbols("x y", nonzero=True)
z = 1 / (x * y)
U_xy = x + y + z
V_xy = x * y + x**-1 + y**-1
W_xy = (x - y) * (y - z) * (z - x)

zero(
    f.subs({U: U_xy, V: V_xy}) - (T - x) * (T - y) * (T - z),
    "ordered-root factorization",
)
zero(D.subs({U: U_xy, V: V_xy, Z: 1}) - W_xy**2, "Vandermonde square identity")

logarithmic_jacobian = sp.det(
    sp.Matrix(
        [
            [x * sp.diff(U_xy, x), y * sp.diff(U_xy, y)],
            [x * sp.diff(V_xy, x), y * sp.diff(V_xy, y)],
        ]
    )
)
zero(logarithmic_jacobian + W_xy, "logarithmic Jacobian is minus Vandermonde")

# The generators of S3 have the expected invariant/alternating action.
transposition = {x: y, y: x}
cycle = {x: y, y: z}
for expression, sign, label in (
    (U_xy, 1, "transposition fixes U"),
    (V_xy, 1, "transposition fixes V"),
    (W_xy, -1, "transposition negates W"),
):
    zero(expression.subs(transposition, simultaneous=True) - sign * expression, label)
for expression, label in (
    (U_xy, "three-cycle fixes U"),
    (V_xy, "three-cycle fixes V"),
    (W_xy, "three-cycle fixes W"),
):
    zero(expression.subs(cycle, simultaneous=True) - expression, label)

# A nontrivial 3-cycle fixes only x=y=z with x^3=1, a finite codimension-two
# set in the torus.  Thus the A3 cover is etale in codimension one; this is
# exactly the quasi-etale Kummer direction and no stronger completion claim.
fixed_cycle_groebner = sp.groebner([x - y, x**2 * y - 1], y, x, order="lex")
expected_fixed_cycle = sp.groebner([y - x, x**3 - 1], y, x, order="lex")
gate(fixed_cycle_groebner == expected_fixed_cycle, "A3 fixed locus is three torus points")


# ---------------------------------------------------------------------------
# IV. Deterministic transcript and frozen replay.
# ---------------------------------------------------------------------------
source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "checker has no inactive Python assert",
)

semantic_lines = [
    "THM3851_TORIC_BRIDGE_INDEPENDENT_AUDIT PASS",
    "FOURIER D(A+wB+w2C,A+w2B+wC,(A+B+C)/3)=9Delta_3;map invertible",
    "PARAMETER q=w2(p+w)/(p+w2) identifies the promoted normalization with (2q+q^-2,q^2+2q^-1)",
    "DISCRIMINANT disc(T^3-UT^2+VT-1)=D(U,V,1) by Sylvester determinant",
    "ROOT_SHEET k[U,V,T]/(f)=k[U,T,T^-1];finite free rank 3;T is a nonconstant unit",
    "DIFFERENT source Jacobian equals -f_T/T",
    "BRANCH normalization q->(2q+q^-2,q^2+2q^-1) is injective Gm;ramified only at q^3=1",
    "PLACES projective closure adds q=0 and q=infinity at two distinct smooth points",
    "ORDERED_ROOTS (x,y,(xy)^-1) in Gm^2;Vandermonde^2=D;log Jacobian=-Vandermonde",
    "QUOTIENT standard invariant-ring lemma gives Gm^2/S3=A2 and Gm^2/A3=(W^2=D)",
    "QUASI_ETALE A3 fixed locus consists of the three triple-root points",
    "OBSTRUCTION explicit root sheet and ordered-root cover are unit-rich, so neither admits a dominant A2 atlas",
    "SCOPE one explicit cyclic Kummer direction;no plane S3 completion,no Keller pair,no JC conclusion",
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
