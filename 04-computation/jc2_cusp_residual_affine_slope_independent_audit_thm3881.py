#!/usr/bin/env python3
"""Import-independent exact audit of THM-3881 Section 7 and its affine-slope extension."""

from __future__ import annotations

import hashlib
import sys
from collections import Counter

import sympy as sp


sys.stdout.reconfigure(newline="\n")
COUNTS: Counter[str] = Counter()


def check(section: str, condition: bool, label: str) -> None:
    COUNTS[section] += 1
    if not condition:
        raise RuntimeError(f"{section}:{label}")


def zero(section: str, expression: sp.Expr, label: str) -> None:
    check(section, sp.factor(expression) == 0, label)


def same(section: str, left: object, right: object, label: str) -> None:
    check(section, left == right, f"{label}: {left!r} != {right!r}")


def homogeneous_part(expression: sp.Expr, degree: int) -> sp.Expr:
    scale = sp.symbols("scale")
    scaled = sp.expand(expression.subs({x: scale * x, y: scale * y}))
    return sp.expand(scaled.coeff(scale, degree))


x, y = sp.symbols("x y")
a = x + 1
L = 9 * x + 4
K = y**2 - 15 * x**2 - 15 * x - 4
P = a * L**2
Q = K * L**2
Delta = sp.expand(a**3 * L**2 - K**2)
F = sp.symbols("F")


# ---------------------------------------------------------------------------
# 1. Reconstruct S(0,f) from the rank-two residual, then factor it twice.
# ---------------------------------------------------------------------------
r0 = K * F
A0 = a * P * F
S0 = sp.expand(
    L**4
    + 2 * (3 * A0 + 3 * P + r0**2) * L**2 * F
    + (8 * A0 + 6 * P + 3 * r0**2) * P * F**2
)
H_first = sp.expand(
    L**2 * (1 + a * F) ** 3 * (1 + 3 * a * F)
    - Delta * F**3 * (2 + 3 * a * F)
)
H_degree = sp.expand(
    L**2 * (1 + 2 * a * F) ** 3
    + K**2 * F**3 * (2 + 3 * a * F)
)
zero("factor", S0 - L**2 * H_first, "S0=L2*H first")
zero("factor", H_first - H_degree, "equivalent H forms")
zero(
    "factor",
    (1 + a * F) ** 3 * (1 + 3 * a * F)
    - (a * F) ** 3 * (2 + 3 * a * F)
    - (1 + 2 * a * F) ** 3,
    "one-variable factor identity",
)

# Direct discriminant reconstruction is independent of the displayed
# factorization.
P_shift = P + 2 * A0 + r0**2
Q_shift = Q + 3 * r0 * P + 3 * r0 * A0 + r0**3
zero("factor", P_shift**3 - Q_shift**2 - Delta * S0, "direct residual reconstruction")


# ---------------------------------------------------------------------------
# 2. The T=0 top degree, including the constant and zero edges.
# ---------------------------------------------------------------------------
K2 = y**2 - 15 * x**2
c = sp.symbols("c", nonzero=True)
H_constant = sp.expand(H_degree.subs(F, c))
H_constant_top = homogeneous_part(H_constant, 5)
expected_constant_top = sp.expand(648 * c**3 * x**5 + 3 * c**4 * x * K2**2)
zero("T0", H_constant_top - expected_constant_top, "constant f degree-five form")
same(
    "T0",
    sp.Poly(H_constant_top, x, y).coeff_monomial(x * y**4),
    3 * c**4,
    "constant f unique xy4 coefficient",
)
same("T0", sp.Poly(H_constant_top, x, y).total_degree(), 5, "constant f odd degree")

zero("T0", H_degree.subs(F, 0) - L**2, "zero f H square")
zero("T0", S0.subs(F, 0) - L**4, "zero f residual square")

# For a nonconstant f of degree n, its leading form f_n controls the unique
# top term 3*x*K2^2*f_n^4 in H.  Generic homogeneous controls check the exact
# all-coefficient formula through several degrees; the proof uses only degree
# additivity in a domain.
u, v = sp.symbols("u v")
for n in range(1, 7):
    fn = u * x**n + v * y**n
    specialized = sp.expand(H_degree.subs(F, fn))
    top = homogeneous_part(specialized, 4 * n + 5)
    expected = sp.expand(3 * x * K2**2 * fn**4)
    zero("T0", top - expected, f"nonconstant leading form n={n}")
    same("T0", sp.Poly(expected, x, y).total_degree(), 4 * n + 5, f"odd degree n={n}")


# ---------------------------------------------------------------------------
# 3. Strict extension: every affine-slope lane T=h*f is empty off f=0.
# ---------------------------------------------------------------------------
alpha, beta, gamma = sp.symbols("alpha beta gamma")
h = alpha * x + beta * y + gamma
R = sp.expand(a * h + K)
M = sp.expand(K * h + a * P)
N = sp.expand(P - h**2)

r_h = R * F
A_h = M * F
T_h = h * F
S_h_direct = sp.expand(
    L**4
    + 2 * (3 * A_h + 3 * P + r_h**2) * L**2 * F
    + (8 * A_h + 6 * P + 3 * r_h**2) * (P * F**2 - T_h**2)
)
S_h_packet = sp.expand(
    L**4
    + 6 * P * L**2 * F
    + 6 * (M * L**2 + P * N) * F**2
    + 2 * (R**2 * L**2 + 4 * M * N) * F**3
    + 3 * R**2 * N * F**4
)
zero("affine_h", S_h_direct - S_h_packet, "affine-slope coefficient packet")

h1 = alpha * x + beta * y
R2 = sp.expand(K2 + x * h1)
same("affine_h", homogeneous_part(R, 2), R2, "R degree-two form")
same("affine_h", sp.Poly(R2, x, y).coeff_monomial(y**2), 1, "R2 cannot cancel")
same("affine_h", homogeneous_part(N, 3), 81 * x**3, "N degree-three form")
same("affine_h", homogeneous_part(M, 4), 81 * x**4, "M degree-four form")

coef1 = sp.expand(6 * P * L**2)
coef2 = sp.expand(6 * (M * L**2 + P * N))
coef3 = sp.expand(2 * (R**2 * L**2 + 4 * M * N))
coef4 = sp.expand(3 * R**2 * N)
same("affine_h", sp.Poly(coef1, x, y).total_degree(), 5, "f coefficient degree")
same("affine_h", sp.Poly(coef2, x, y).total_degree(), 6, "f2 coefficient degree")
same("affine_h", sp.Poly(coef3, x, y).total_degree(), 7, "f3 coefficient degree")
same("affine_h", sp.Poly(coef4, x, y).total_degree(), 7, "f4 coefficient degree")
zero("affine_h", homogeneous_part(coef4, 7) - 243 * x**3 * R2**2, "f4 degree-seven form")
zero("affine_h", homogeneous_part(coef3, 7) - 52488 * x**7, "f3 degree-seven form")

# For deg(f)>=1, the f^4 row has degree 4n+7 while all lower rows cap at
# 3n+7, 2n+6, n+5, 4.  Generic leading-form controls verify no hidden sign or
# cancellation in the top row.
for n in range(1, 6):
    fn = u * x**n + v * y**n
    specialized = sp.expand(S_h_packet.subs(F, fn))
    top = homogeneous_part(specialized, 4 * n + 7)
    expected = sp.expand(243 * x**3 * R2**2 * fn**4)
    zero("affine_h", top - expected, f"affine h leading form n={n}")
    same("affine_h", sp.Poly(expected, x, y).total_degree(), 4 * n + 7, f"affine h odd degree n={n}")

# For a nonzero constant f=c, f^3 and f^4 tie at degree seven.  Their sum is
# still nonzero, detected by the x^3*y^4 coefficient 243*c^4.
S_h_constant_top = homogeneous_part(sp.expand(S_h_packet.subs(F, c)), 7)
expected_h_constant_top = sp.expand(
    81 * x**3 * c**3 * (648 * x**4 + 3 * c * R2**2)
)
zero("affine_h", S_h_constant_top - expected_h_constant_top, "affine h constant top form")
same(
    "affine_h",
    sp.Poly(S_h_constant_top, x, y).coeff_monomial(x**3 * y**4),
    243 * c**4,
    "affine h constant x3y4 coefficient",
)
same("affine_h", sp.Poly(S_h_constant_top, x, y).total_degree(), 7, "affine h constant odd degree")
zero("affine_h", S_h_packet.subs(F, 0) - L**4, "affine h zero survivor")


# ---------------------------------------------------------------------------
# 4. Characteristic boundary.
# ---------------------------------------------------------------------------
# The proof works over every field of characteristic !=3.  In characteristic
# three, both 3 and the degree-three leading coefficient 81 of P vanish.
# Already f=x^2 makes the T=0 degree even, so the odd-degree proof genuinely
# stops there (this is a proof-boundary witness, not a square counterexample).
H_char3 = sp.Poly(H_degree.subs(F, x**2), x, y, modulus=3)
same("scope", H_char3.total_degree(), 10, "char3 parity boundary witness")
same("scope", H_char3.total_degree() % 2, 0, "char3 degree is even")
for prime in (2, 5, 7, 11):
    check("scope", 3 % prime != 0 and 81 % prime != 0, f"top coefficients survive char {prime}")


semantic = (
    "THM-3881 Section 7 hostile audit",
    "S(0,f)=L^2 H_f with two exact H factorizations",
    "nonconstant f gives unique H degree 4n+5",
    "nonzero constant f gives unique xy4 coefficient 3c4 and H degree five",
    "address condition unnecessary for the T zero square iff f zero",
    "strict affine-slope extension T=h f for degree h at most one",
    "nonconstant f gives unique S degree 4n+7",
    "constant f gives x3y4 coefficient 243c4 and S degree seven",
    "valid over any field of characteristic not three; characteristic three remains outside this proof",
)
semantic_hash = hashlib.sha256(repr(semantic).encode()).hexdigest()

print("AUDIT=THM-3881-SECTION-7")
print("BASE=de787fd3a")
print("T_ZERO_FACTOR=PASS;S=L^2*H;two-H-forms-agree")
print("T_ZERO_STRONG=PASS;square-iff-f=0;address-condition-unnecessary")
print("AFFINE_SLOPE=PASS;T=h*f;deg(h)<=1;square-iff-f=0")
print("NONCONSTANT_TOP=H:4n+5;S_affine_h:4n+7")
print("CONSTANT_TOP=H:xy4-coeff-3c4;S_affine_h:x3y4-coeff-243c4")
print("FIELD_SCOPE=any-field-char-not-3;algebraic-closure-unused")
print("CHAR3=degree-parity-proof-fails;no-counterexample-claimed")
print("SECTION_GATES=" + ";".join(f"{name}:{COUNTS[name]}" for name in sorted(COUNTS)))
print(f"SEMANTIC_SHA256={semantic_hash}")
print(f"CHECKS={sum(COUNTS.values())}")
print("RESULT=PASS")
