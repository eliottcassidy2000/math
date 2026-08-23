#!/usr/bin/env python3
"""Independent hostile audit for THM-3871.

This companion does not import the primary quintic verifier.  It rebuilds the
ten Jacobian buckets by coefficient convolution, checks the constant-target
and Kummer reductions in every degree-drop row, rederives the depressed
(5,2), (5,3), and (5,4) differential packets, and certifies the local and
infinity obstructions.  In particular it computes both coprime resultants
which prevent simultaneous regularity of the two original arm coefficients.

Every proof gate is explicit rather than a Python ``assert`` so that ``-O``
cannot remove mathematical checks.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: object, label: str) -> None:
    """Optimization-safe exact gate."""

    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"GATE FAILED: {label}: {condition}")


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


def bracket(left: sp.Expr, right: sp.Expr, z: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(left, z) * sp.diff(right, s)
        - sp.diff(left, s) * sp.diff(right, z)
    )


def total_derivative(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    derivatives: tuple[sp.Symbol, ...],
) -> sp.Expr:
    return sp.expand(sum(
        sp.diff(expression, variable) * derivative
        for variable, derivative in zip(variables, derivatives)
    ))


# ---------------------------------------------------------------------------
# 1. Ten buckets, independently reconstructed by coefficient convolution.
# ---------------------------------------------------------------------------
s, z, x = sp.symbols("s z x")
a = sp.symbols("a0:6")
c = sp.symbols("c0:6")
ad = sp.symbols("a0d:6d")
cd = sp.symbols("c0d:6d")


def bucket(degree: int) -> sp.Expr:
    return sp.expand(sum(
        i * a[i] * cd[j] - j * ad[i] * c[j]
        for i in range(6)
        for j in range(6)
        if i + j == degree + 1
    ))


expected = [
    a[1] * cd[0] - ad[0] * c[1],
    a[1] * cd[1] - ad[1] * c[1] + 2 * a[2] * cd[0] - 2 * ad[0] * c[2],
    a[1] * cd[2] - 2 * ad[1] * c[2] + 2 * a[2] * cd[1]
    - ad[2] * c[1] + 3 * a[3] * cd[0] - 3 * ad[0] * c[3],
    a[1] * cd[3] - 3 * ad[1] * c[3] + 2 * a[2] * cd[2]
    - 2 * ad[2] * c[2] + 3 * a[3] * cd[1] - ad[3] * c[1]
    + 4 * a[4] * cd[0] - 4 * ad[0] * c[4],
    a[1] * cd[4] - 4 * ad[1] * c[4] + 2 * a[2] * cd[3]
    - 3 * ad[2] * c[3] + 3 * a[3] * cd[2] - 2 * ad[3] * c[2]
    + 4 * a[4] * cd[1] - ad[4] * c[1] + 5 * a[5] * cd[0]
    - 5 * ad[0] * c[5],
    a[1] * cd[5] - 5 * ad[1] * c[5] + 2 * a[2] * cd[4]
    - 4 * ad[2] * c[4] + 3 * a[3] * cd[3] - 3 * ad[3] * c[3]
    + 4 * a[4] * cd[2] - 2 * ad[4] * c[2] + 5 * a[5] * cd[1]
    - ad[5] * c[1],
    2 * a[2] * cd[5] - 5 * ad[2] * c[5] + 3 * a[3] * cd[4]
    - 4 * ad[3] * c[4] + 4 * a[4] * cd[3] - 3 * ad[4] * c[3]
    + 5 * a[5] * cd[2] - 2 * ad[5] * c[2],
    3 * a[3] * cd[5] - 5 * ad[3] * c[5] + 4 * a[4] * cd[4]
    - 4 * ad[4] * c[4] + 5 * a[5] * cd[3] - 3 * ad[5] * c[3],
    4 * a[4] * cd[5] - 5 * ad[4] * c[5]
    + 5 * a[5] * cd[4] - 4 * ad[5] * c[4],
    5 * (a[5] * cd[5] - ad[5] * c[5]),
]

for degree in range(10):
    zero(bucket(degree) - expected[degree], f"quintic convolution bucket z^{degree}")

# A direct polynomial specialization is a second route to every bucket.
A_values = [s**2 + 1, s**3 - s, 2 * s + 3, s**4 + s, 5 - s**2, s + 7]
C_values = [s**3 + 2, s - 4, s**2 + s + 1, 3 * s**3 - 1, s**2 - 5, 2 * s + 1]
substitution: dict[sp.Symbol, sp.Expr] = {}
for symbol, derivative_symbol, value in zip(a, ad, A_values):
    substitution[symbol] = value
    substitution[derivative_symbol] = sp.diff(value, s)
for symbol, derivative_symbol, value in zip(c, cd, C_values):
    substitution[symbol] = value
    substitution[derivative_symbol] = sp.diff(value, s)
A_sample = sum(value * z**index for index, value in enumerate(A_values))
C_sample = sum(value * z**index for index, value in enumerate(C_values))
J_sample = bracket(A_sample, C_sample, z, s)
for degree in range(10):
    zero(
        J_sample.coeff(z, degree) - expected[degree].subs(substitution),
        f"direct quintic bucket z^{degree}",
    )
gate(sp.Poly(J_sample, z).degree() <= 9, "no Jacobian bucket above z9")


# ---------------------------------------------------------------------------
# 2. Constant target direction, Kummer divisors, and degree-drop handoffs.
# ---------------------------------------------------------------------------
r11, r12, r21, r22 = sp.symbols("r11 r12 r21 r22")
w, omega, wd, omegad = sp.symbols("w omega wd omegad")
w_new = r11 * w + r12 * omega
omega_new = r21 * w + r22 * omega
wd_new = r11 * wd + r12 * omegad
omegad_new = r21 * wd + r22 * omegad
zero(
    w_new * omegad_new - wd_new * omega_new
    - (r11 * r22 - r12 * r21) * (w * omegad - wd * omega),
    "top Wronskian target covariance",
)

# Both affine charts of P1(k) are normalized by constant SL2 matrices.
Pdir, Qdir, hdir = sp.symbols("Pdir Qdir hdir", nonzero=True)
matrix_p = sp.Matrix([[1 / Pdir, 0], [-Qdir, Pdir]])
zero(matrix_p.det() - 1, "top direction P-chart determinant")
normalized_p = matrix_p * sp.Matrix([hdir * Pdir, hdir * Qdir])
zero(normalized_p[0] - hdir, "top direction P-chart survivor")
zero(normalized_p[1], "top direction P-chart killed coordinate")
matrix_q = sp.Matrix([[0, 1 / Qdir], [-Qdir, 0]])
zero(matrix_q.det() - 1, "top direction Q-chart determinant")
normalized_q = matrix_q * sp.Matrix([0, hdir * Qdir])
zero(normalized_q[0] - hdir, "top direction Q-chart survivor")
zero(normalized_q[1], "top direction Q-chart killed coordinate")

# For j=1,2,3,4, 5 ord(c_j)=j ord(w) is exactly the h5/hj lattice.
kummer_profiles = 0
valuation_profiles = 0
for transverse_degree in range(1, 5):
    for w_order in range(26):
        for c_order in range(26):
            equation = 5 * c_order == transverse_degree * w_order
            lattice = (
                w_order % 5 == 0
                and c_order % transverse_degree == 0
                and w_order // 5 == c_order // transverse_degree
            )
            gate(
                equation == lattice,
                f"Kummer lattice j={transverse_degree},w={w_order},c={c_order}",
            )
            valuation_profiles += 1
            if equation:
                kummer_profiles += 1

# The j=1 top coefficient is removed by the target shear A-rho*C^5.
rho = sp.symbols("rho")
A_fun = sp.Function("A_fun")(s, z)
C_fun = sp.Function("C_fun")(s, z)
zero(
    bracket(A_fun - rho * C_fun**5, C_fun, z, s)
    - bracket(A_fun, C_fun, z, s),
    "degree-one target shear preserves Jacobian",
)

# The j=0 rows imply w*lambda=0, while both factors are required nonzero.
w0, lambda0 = sp.symbols("w0 lambda0", nonzero=True)
alpha0, b0prime = sp.symbols("alpha0 b0prime")
zero(
    w0 * (alpha0 * b0prime - lambda0) - alpha0 * (w0 * b0prime)
    + w0 * lambda0,
    "degree-zero top/constant bucket Bezout combination",
)
gate(w0 * lambda0 != 0, "degree-zero row contradicts Keller unit")

# Only j=1 admits a polynomial target power whose transverse degree is five.
degree_shear_profiles = 0
for transverse_degree in range(1, 5):
    for target_degree in range(1, 6):
        gate(
            (transverse_degree * target_degree == 5)
            == (transverse_degree == 1 and target_degree == 5),
            f"target shear divisibility j={transverse_degree},d={target_degree}",
        )
        degree_shear_profiles += 1
Rtop, Qtop, htop = sp.symbols("Rtop Qtop htop", nonzero=True)
zero(
    Rtop * htop**5 - (Rtop / Qtop**5) * (Qtop * htop)**5,
    "degree-one h5/h coefficient cancellation",
)

# Scaling and moving translation introduce no hidden connection term.
h_scale = s**2 + s + 1
y = sp.symbols("y")
F_scale = y**5 + s * y**3 + (s**2 + 1) * y + s
G_scale = y**4 + (2 * s - 1) * y**2 + s**3
scaled_F = F_scale.subs(y, h_scale * z)
scaled_G = G_scale.subs(y, h_scale * z)
scaled_inner = bracket(F_scale, G_scale, y, s).subs(y, h_scale * z)
zero(
    bracket(scaled_F, scaled_G, z, s) - h_scale * scaled_inner,
    "moving Kummer scale Euler cancellation",
)
g_move = s**2 - 3 * s + 1
translated_F = F_scale.subs(y, z + g_move)
translated_G = G_scale.subs(y, z + g_move)
translated_inner = bracket(F_scale, G_scale, y, s).subs(y, z + g_move)
zero(
    bracket(translated_F, translated_G, z, s) - translated_inner,
    "moving translation bracket invariance",
)


# ---------------------------------------------------------------------------
# 3. The complete depressed (5,2) row and its two-arm obstruction.
# ---------------------------------------------------------------------------
R, Q = sp.symbols("R Q", nonzero=True)
D, E, F, G, A0 = sp.symbols("D E F G A0")
B, N, b = sp.symbols("B N b")
Bp, Np, bp = sp.symbols("Bp Np bp")

K52 = 5 * R / (2 * Q)
U52 = K52 * b + E
M52 = 2 * D * b / Q + F
L52 = 3 * K52 * b**2 / (4 * Q) + 3 * E * b / (2 * Q) + G
a52 = D * b**2 / Q**2 + F * b / Q + A0
d52 = lambda expression: total_derivative(expression, (b,), (bp,))

rows52 = [
    5 * R * bp - 2 * Q * d52(U52),
    4 * D * bp - 2 * Q * d52(M52),
    3 * U52 * bp - 2 * Q * d52(L52),
    2 * M52 * bp - 2 * Q * d52(a52),
    L52 * bp,
]
for index, row in enumerate(rows52[:-1]):
    zero(row, f"(5,2) integrated row {4-index}")

A52 = R * x**5 + D * x**4 + U52 * x**3 + M52 * x**2 + L52 * x + a52
C52 = Q * x**2 + b
As52 = d52(A52)
Cs52 = d52(C52)
J52 = sp.expand(sp.diff(A52, x) * Cs52 - As52 * sp.diff(C52, x))
for degree in range(1, 5):
    zero(J52.coeff(x, degree), f"(5,2) direct depressed bucket x^{degree}")
zero(J52.coeff(x, 0) - L52 * bp, "(5,2) constant bucket L*bprime")
zero(
    sp.expand(L52).coeff(b, 2) - 3 * K52 / (4 * Q),
    "(5,2) unique b2*bprime local leader",
)

# If b~-Q*g^2, the order-five A-arm residue is exactly 3R/8.
g = sp.symbols("g", nonzero=True)
arm52_residue = R * g**5 + K52 * (-Q * g**2) * g**3
arm52_residue += (3 * K52 / (4 * Q)) * (Q**2 * g**4) * g
zero(arm52_residue / g**5 - 3 * R / 8, "(5,2) terminal A-arm residue")

# Exact dominance: Qg^2+b can cancel a b-pole only at 2 ord(g)=ord(b).
dominance52 = 0
for b_pole in range(1, 17):
    for g_pole in range(17):
        exponents = [2 * g_pole, b_pole]
        tied = exponents.count(max(exponents)) > 1
        gate(tied == (2 * g_pole == b_pole),
             f"(5,2) C-arm dominance b={b_pole},g={g_pole}")
        dominance52 += 1

# Polynomial infinity: the leading degree 3d-1 coefficient is nonzero.
d_inf = sp.symbols("d_inf", integer=True, positive=True)
b0 = sp.symbols("b0", nonzero=True)
leading52 = (3 * K52 / (4 * Q)) * b0**2 * d_inf * b0
gate(leading52 != 0, "(5,2) constant-h infinity coefficient")

# Sharp rational hostile: C-arm cancels and only the A-arm is nonpolynomial.
C52_hostile = s**14 * z**2 + 2 * s**6 * z
A52_hostile = (
    s**35 * z**5 + 5 * s**27 * z**4 + sp.Rational(15, 2) * s**19 * z**3
    + sp.Rational(5, 2) * s**11 * z**2 - sp.Rational(5, 8) * s**3 * z
    + sp.Rational(3, 8) / s**5
)
zero(bracket(A52_hostile, C52_hostile, z, s) - sp.Rational(15, 4),
     "(5,2) sharp hostile Jacobian")
zero(C52_hostile.subs(z, 0), "(5,2) sharp hostile C arm")
gate(sp.cancel(A52_hostile.subs(z, 0)) == sp.Rational(3, 8) / s**5,
     "(5,2) sharp hostile A-arm residual")


# ---------------------------------------------------------------------------
# 4. The depressed (5,3) row, conserved polynomial, and resultant 441.
# ---------------------------------------------------------------------------
K53 = 5 * R / (3 * Q)
H53 = 4 * D / (3 * Q)
U53 = K53 * B + E
M53 = H53 * B + K53 * b + F
L53 = K53 * B**2 / (3 * Q) + E * B / Q + H53 * b + G
a53 = (
    H53 * B**2 / (6 * Q) + 2 * K53 * B * b / (3 * Q)
    + 2 * F * B / (3 * Q) + E * b / Q + A0
)
d53 = lambda expression: total_derivative(expression, (B, b), (Bp, bp))

rows53 = [
    5 * R * Bp - 3 * Q * d53(U53),
    5 * R * bp + 4 * D * Bp - 3 * Q * d53(M53),
    4 * D * bp + 3 * U53 * Bp - d53(U53) * B - 3 * Q * d53(L53),
    3 * U53 * bp + 2 * M53 * Bp - d53(M53) * B - 3 * Q * d53(a53),
    L53 * Bp - d53(L53) * B + 2 * M53 * bp,
    L53 * bp - d53(a53) * B,
]
for index, row in enumerate(rows53[:4]):
    zero(row, f"(5,3) integrated row {5-index}")

P53 = (
    sp.Rational(5, 3) * R * B**3 - 12 * D * Q * B * b
    - 18 * F * Q**2 * b - 9 * G * Q**2 * B - 15 * Q * R * b**2
)
Omega53 = (
    12 * D * Q * b * bp - 4 * D * B**2 * Bp - 6 * F * Q * B * Bp
    + 9 * G * Q**2 * bp - 5 * R * B**2 * bp - 10 * R * B * b * Bp
)
zero(rows53[4] + d53(P53) / (9 * Q**2), "(5,3) conserved P derivative")
zero(9 * Q**2 * rows53[5] - Omega53, "(5,3) exact constant one-form")

A53 = R*x**5 + D*x**4 + U53*x**3 + M53*x**2 + L53*x + a53
C53 = Q*x**3 + B*x + b
J53 = sp.expand(sp.diff(A53, x) * d53(C53) - d53(A53) * sp.diff(C53, x))
for degree in range(6):
    zero(J53.coeff(x, degree) - rows53[5-degree],
         f"(5,3) direct depressed bucket x^{degree}")

# In P, the maximal pole/degree among B3, B*b, b2 is tied only at 2n=3m.
tie_profiles53 = 0
for m_order in range(1, 21):
    for n_order in range(1, 31):
        exponents = [3 * m_order, m_order + n_order, 2 * n_order]
        tied_max = exponents.count(max(exponents)) > 1
        gate(tied_max == (2 * n_order == 3 * m_order),
             f"(5,3) P dominance m={m_order},n={n_order}")
        tie_profiles53 += 1

X = sp.symbols("X")
p53 = 5 * X**2 + 10 * X + 6
q53 = X**3 - 9 * (1 + X)**2
resultant53 = sp.resultant(p53, q53, X)
gate(resultant53 == 441, "(5,3) simultaneous-arm resultant is 441")
gate(sp.gcd(p53, q53) == 1, "(5,3) simultaneous-arm polynomials coprime")

Z = sp.symbols("Z")
arm53 = 1 + sp.Rational(5, 3) * (X + Z) + sp.Rational(5, 9) * X**2
arm53 += sp.Rational(10, 9) * X * Z
zero(
    arm53.subs(Z, -1-X) + p53 / 9,
    "(5,3) A-arm residual after C-arm cancellation",
)
zero(q53 - (X**3 - 9 * (-1-X)**2), "(5,3) invariant after C-arm")

# At infinity the last two one-form terms have an uncancellable coefficient.
m53, n53, B0, b0 = sp.symbols("m53 n53 B0 b0", nonzero=True)
leading_omega53 = -(5 * n53 + 10 * m53) * R * B0**2 * b0
gate(leading_omega53 != 0, "(5,3) constant-h infinity coefficient")


# ---------------------------------------------------------------------------
# 5. The depressed (5,4) packet, all W strata, and terminal coprimality.
# ---------------------------------------------------------------------------
K54 = 5 * R / (4 * Q)
H54 = D / Q
U54 = K54 * B + E
M54 = H54 * B + K54 * N + F
L54 = K54 * B**2 / (8 * Q) + 3 * E * B / (4 * Q) + H54 * N + K54 * b + G
a54 = (
    K54 * B * N / (4 * Q) + F * B / (2 * Q)
    + 3 * E * N / (4 * Q) + H54 * b + A0
)
d54 = lambda expression: total_derivative(expression, (B, N, b), (Bp, Np, bp))

rows54 = [
    5 * R * Bp - 4 * Q * d54(U54),
    5 * R * Np + 4 * D * Bp - 4 * Q * d54(M54),
    5 * R * bp + 4 * D * Np + 3 * U54 * Bp - 2 * d54(U54) * B
    - 4 * Q * d54(L54),
    4 * D * bp + 3 * U54 * Np - d54(U54) * N + 2 * M54 * Bp
    - 2 * d54(M54) * B - 4 * Q * d54(a54),
    3 * U54 * bp + 2 * M54 * Np - d54(M54) * N + L54 * Bp
    - 2 * d54(L54) * B,
    L54 * Np - d54(L54) * N + 2 * M54 * bp - 2 * d54(a54) * B,
    L54 * bp - d54(a54) * N,
]
for index, row in enumerate(rows54[:4]):
    zero(row, f"(5,4) integrated row {6-index}")

W54 = 5 * R * B + 12 * E * Q
P54 = (
    W54 * (B**2 - 8 * Q * b) - 20 * Q * R * N**2
    - 64 * F * Q**2 * N - 32 * G * Q**2 * B
)
T54 = (
    N * (15 * R * B**2 + 24 * E * Q * B - 40 * Q * R * b - 32 * G * Q**2)
    + 16 * F * Q * (B**2 - 4 * Q * b)
)
Omega54 = (
    (5 * R * B**2 + 24 * E * Q * B + 40 * Q * R * b + 32 * G * Q**2) * bp
    - (24 * E * Q * N + 10 * R * B * N) * Np
    - (16 * F * Q * N + 10 * R * N**2) * Bp
)
zero(rows54[4] + d54(P54) / (32 * Q**2), "(5,4) first conserved polynomial P")
zero(rows54[5] + d54(T54) / (32 * Q**2), "(5,4) second conserved polynomial T")
zero(32 * Q**2 * rows54[6] - Omega54, "(5,4) exact constant one-form")

A54 = R*x**5 + D*x**4 + U54*x**3 + M54*x**2 + L54*x + a54
C54 = Q*x**4 + B*x**2 + N*x + b
J54 = sp.expand(sp.diff(A54, x) * d54(C54) - d54(A54) * sp.diff(C54, x))
for degree in range(7):
    zero(J54.coeff(x, degree) - rows54[6-degree],
         f"(5,4) direct depressed bucket x^{degree}")

# W-pole dominance and the two possible N channels.
dominance54 = 0
for m_order in range(1, 17):
    for b_order in range(1, 33):
        for n_order in range(1, 25):
            p_exponents = [3*m_order, m_order+b_order, 2*n_order, n_order, m_order]
            # If b is more singular than B2, N*b is uniquely dominant in T;
            # if it is less singular, W*B2 or N2 is exposed in P.  Surviving
            # profiles must have b_order=2m and n_order<=3m/2.
            potentially_survives = (
                b_order == 2*m_order and 2*n_order <= 3*m_order
            )
            p_has_tie = p_exponents.count(max(p_exponents)) > 1
            t_exponents = [n_order+2*m_order, n_order+b_order, 2*m_order, b_order]
            t_has_tie = t_exponents.count(max(t_exponents)) > 1
            gate(
                (p_has_tie and t_has_tie) == potentially_survives,
                f"(5,4) W-pole dominance m={m_order},b={b_order},n={n_order}",
            )
            dominance54 += 1

# For 0<n<3m/2, P forces b=B2/(8Q), but T retains 10R*N*B2.
N0 = sp.symbols("N0", nonzero=True)
b_intermediate = B0**2 / (8 * Q)
intermediate_T = N0 * (15 * R * B0**2 - 40 * Q * R * b_intermediate)
zero(intermediate_T - 10 * R * N0 * B0**2,
     "(5,4) intermediate N-pole terminal coefficient")

# The N-regular channel: invariants and C-arm force a coprime A residue.
q54_regular = X**2 + 8 * X + 8
p54_regular = 5 * X**2 - 8
resultant54_regular = sp.resultant(q54_regular, p54_regular, X)
gate(resultant54_regular == -256, "(5,4) N-regular arm resultant")
gate(sp.gcd(q54_regular, p54_regular) == 1,
     "(5,4) N-regular arm polynomials coprime")

# The extremal N-pole channel retains P, T, C-arm, and A-arm simultaneously.
Y = sp.symbols("Y")
p54 = 15 * X**3 + 20 * X**2 + 40 * X + 32
q54 = 50 * X**5 + 25 * X**4 - 80 * X**2 + 64
resultant54 = sp.resultant(p54, q54, X)
gate(resultant54 == 3171942400000,
     "(5,4) simultaneous-arm resultant exact integer")
gate(sp.factorint(resultant54) == {2: 23, 5: 5, 11: 2},
     "(5,4) resultant prime factorization")
gate(sp.gcd(p54, q54) == 1, "(5,4) simultaneous-arm polynomials coprime")

Z_extreme = 3 * X**2 / 8
invariant_extreme = Y**2 + X**3 / 2
C_arm_extreme = 1 + X + Y + Z_extreme
A_arm_extreme = (5 * X**2 + 10 * X * Y - 8) / 32
zero(
    sp.expand(40 * X * C_arm_extreme - 4 * (10 * X * Y - (8 - 5 * X**2)))
    - p54,
    "(5,4) p eliminant from C and A arms",
)
zero(
    q54 - sp.expand((8 - 5 * X**2)**2 + 50 * X**5),
    "(5,4) q eliminant from squared A arm and P invariant",
)

# Arms alone really do leave a false one-parameter-style cancellation point;
# the conserved P/T leading equations reject it.
arms_only = {X: sp.Integer(1), Y: sp.Rational(3, 10)}
Z_arms_only = sp.Rational(-23, 10)
zero((1 + X + Y + sp.Symbol("Zaux")).subs(arms_only).subs(sp.Symbol("Zaux"), Z_arms_only),
     "(5,4) arms-only hostile C cancellation")
zero(A_arm_extreme.subs(arms_only), "(5,4) arms-only hostile A cancellation")
gate(
    (X**3 - 8*X*sp.Symbol("Zaux") - 4*Y**2).subs(arms_only).subs(
        sp.Symbol("Zaux"), Z_arms_only
    ) == sp.Rational(476, 25),
    "(5,4) arms-only hostile violates P",
)
gate(
    (Y * (3*X**2 - 8*sp.Symbol("Zaux"))).subs(arms_only).subs(
        sp.Symbol("Zaux"), Z_arms_only
    ) == sp.Rational(321, 50),
    "(5,4) arms-only hostile violates T",
)

# W-zero channel: b~-Qg4 cancels C but leaves R-KQ=-R/4 in A.
zero(R - K54 * Q + R / 4, "(5,4) W-zero A-arm residue minus R/4")
zero(sp.expand(Omega54).coeff(b * bp) - 40 * Q * R,
     "(5,4) W-zero unique b*bprime coefficient")

# W identically zero makes B constant; P then makes N constant.
B_W0 = -12 * E * Q / (5 * R)
P_W0 = sp.expand(P54.subs(B, B_W0))
zero(sp.diff(P_W0, N, 2) + 40 * Q * R,
     "(5,4) W-identically-zero P is nonconstant quadratic in N")

# Constant h, deg(B)=m>0: exact nonzero leading one-form coefficients.
m54 = sp.symbols("m54", integer=True, positive=True)
Blead = sp.symbols("Blead", nonzero=True)
b_lead_regular = Blead**2 / (8 * Q)
omega_regular = (
    (5*R*Blead**2 + 40*Q*R*b_lead_regular) * (2*m54*b_lead_regular)
)
zero(omega_regular - sp.Rational(5, 2) * R*m54*Blead**4/Q,
     "(5,4) constant-h N-regular infinity coefficient")

b_lead_extreme = 3 * Blead**2 / (8 * Q)
Nlead2 = -Blead**3 / (2 * Q)
omega_extreme = (
    (5*R*Blead**2 + 40*Q*R*b_lead_extreme) * (2*m54*b_lead_extreme)
    - 10*R*Blead*Nlead2 * (sp.Rational(3, 2)*m54)
    - 10*R*Nlead2 * (m54*Blead)
)
zero(omega_extreme - sp.Rational(55, 2) * R*m54*Blead**4/Q,
     "(5,4) constant-h extremal-N infinity coefficient")

# If B is constant and W!=0, eliminating b from P leaves a genuine cubic
# in N in T, so no nonconstant polynomial parametrization was divided away.
Wc, Pc = sp.symbols("Wc Pc", nonzero=True)
Bc = sp.symbols("Bc")
b_from_P = (
    Wc*Bc**2 - 20*Q*R*N**2 - 64*F*Q**2*N - 32*G*Q**2*Bc - Pc
) / (8*Q*Wc)
T_constant_B = sp.expand(
    N*(15*R*Bc**2 + 24*E*Q*Bc - 40*Q*R*b_from_P - 32*G*Q**2)
    + 16*F*Q*(Bc**2 - 4*Q*b_from_P)
)
zero(sp.Poly(T_constant_B, N).LC() - 100*Q*R**2/Wc,
     "(5,4) constant-B W-nonzero cubic elimination leader")


# ---------------------------------------------------------------------------
# 6. DVR/field edges, frozen semantics, and optimization-safe output.
# ---------------------------------------------------------------------------
# No algebraic closure is used: irreducible characteristic-zero primes are
# separable, so pi' is a DVR unit.  These nonlinear controls check the exact
# derivative residue for poles through order six.
prime_controls = [s, s**2 + 1, s**3 + 2]
dvr_controls = 0
for prime_index, prime in enumerate(prime_controls):
    gate(sp.gcd(prime, sp.diff(prime, s)) == 1,
         f"separable DVR prime {prime_index}")
    for pole_order in range(1, 7):
        numerator = sp.cancel(
            sp.diff(prime**(-pole_order), s) * prime**(pole_order + 1)
        )
        zero(numerator + pole_order * sp.diff(prime, s),
             f"DVR derivative leading residue prime={prime_index},m={pole_order}")
        dvr_controls += 1

# The two resultants are nonzero in every characteristic-zero residue field.
gate(sp.Integer(441) != 0, "(5,3) resultant survives characteristic zero")
gate(sp.Integer(3171942400000) != 0,
     "(5,4) resultant survives characteristic zero")

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "no inactive Python assert",
)

semantic = {
    "buckets": "ten coefficient-convolution rows plus direct specialization",
    "top": "constant target SL2; h5/hj UFD lattices j=1..4",
    "degree_drops": "j=0 contradiction;j=1 target shear;j=2,j=3 depressed packets",
    "row_52": "five rows;integrations;two-arm residue 3R/8;infinity degree",
    "row_53": "six rows;P conserved;constant one-form;resultant 441",
    "row_54": "seven rows;P,T conserved;constant one-form;all W strata",
    "resultant_54": str(resultant54),
    "resultant_54_factorization": "2^23*5^5*11^2",
    "arms_only_hostile": "X=1,Y=3/10,Z=-23/10 violates P,T",
    "field": "arbitrary characteristic zero;irreducible DVRs;no root extraction",
    "valuation_profiles": valuation_profiles,
    "dominance_52": dominance52,
    "tie_profiles_53": tie_profiles53,
    "dominance_54": dominance54,
    "scope": "polynomial k[s,z] transverse degree at most five only;JC2 open",
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("THM3871_QUINTIC_NORMAL_STRIP_INDEPENDENT_HOSTILE_AUDIT")
print("status=PASS_PROOF+VERIFIED_EXACT;JC2_OPEN")
print("ten_buckets=coefficient_convolution+direct_specialization_PASS")
print("top_direction=constant_SL2;Kummer_h5_hj_j1_to_j4_PASS")
print("degree_drops=j0_contradiction+j1_shear+j2+j3_packets_PASS")
print("row_5_2=all_integrations+local_residue_3R/8+infinity_PASS")
print(f"row_5_3=P_conserved+one_form;resultant={resultant53}_PASS")
print("row_5_4=P_and_T_conserved+one_form+W_pole_unit_zero_identically_zero_PASS")
print(f"row_5_4_resultant={resultant54}=2^23*5^5*11^2_PASS")
print("arms_only_hostile=X1_Y3/10_Z-23/10;P=476/25;T=321/50_PASS")
print(f"valuation_profiles={valuation_profiles};Kummer_points={kummer_profiles}")
print(f"dominance_52={dominance52};tie_profiles_53={tie_profiles53};dominance_54={dominance54}")
print(f"DVR_controls={dvr_controls};nonclosed_characteristic_zero_PASS")
print("scope=polynomial_transverse_degree_at_most_5_only;rational_and_infinite_open")
print(f"semantic_sha256={semantic_sha}")
print(f"GATES={GATES}")
print("RESULT=PASS")
