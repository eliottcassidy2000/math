#!/usr/bin/env python3
"""Exact symbolic audit of five-face renewal propagation.

For the fixed cubic inverse chart of THM-3495, suppose P has the two
renewal faces in the packet A(e,m), and put Q=L^e Norm(P).  This companion
checks the unique common endpoint, both hybrid inverse limits, the nonmonic
Vieta normalization, every exponent identity, and exact scalar propagation.

The checks are generic in e,m.  The pinned L->H, H->J, and J->G controls
include cases where the leading value is not root-independent.  The fixed
G->R5 and R5->R6 consequences are then evaluated exactly.

Scope: the fixed inverse chart, conditional on Q being polynomial.  This
does not prove polynomiality at another rung, a finite-sheet unit, image
status, degree 243, an all-level tower, or a general Jacobian claim.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression, label: str) -> None:
    require(sp.simplify(expression) == 0, f"{label}: {sp.factor(expression)}")


def polynomial_face(expression, vector: tuple[int, int, int], *, take_minimum: bool):
    polynomial = sp.Poly(expression, A, B, C, domain=sp.QQ)
    rows = [
        (sum(exponent * weight for exponent, weight in zip(monomial, vector)), monomial, coefficient)
        for monomial, coefficient in polynomial.terms()
    ]
    extreme = min(row[0] for row in rows) if take_minimum else max(row[0] for row in rows)
    face = sum(
        coefficient * A**monomial[0] * B**monomial[1] * C**monomial[2]
        for row_weight, monomial, coefficient in rows
        if row_weight == extreme
    )
    return extreme, sp.expand(face)


def require_face_ledger(actual, expected, label: str) -> None:
    require(actual.keys() == expected.keys(), f"{label}: invariant names changed")
    for name in actual:
        actual_weight, actual_face = actual[name]
        expected_weight, expected_face = expected[name]
        require(actual_weight == expected_weight, f"{label}: {name} weight changed")
        require_zero(actual_face - expected_face, f"{label}: {name} face changed")


def require_zero_mod_cubic(expression, cubic, label: str) -> None:
    numerator, _ = sp.together(expression).as_numer_denom()
    coefficient_field = sp.QQ.frac_field(A, B, C)
    remainder = sp.rem(
        sp.Poly(numerator, q, domain=coefficient_field),
        sp.Poly(cubic, q, domain=coefficient_field),
    )
    require(remainder.is_zero, f"{label}: residual remainder {remainder}")


# --- Generic packet arithmetic -----------------------------------------

e, m = sp.symbols("e m", integer=True, nonnegative=True)
third = sp.Rational(1, 3)

r = e - 2 * m * third
a = 2 * e - 4 * m * third
b = 2 * e - 2 * m * third
gamma0 = -8 * e + 2 * m
h = a - 3 * b
delta6 = gamma0 - b
delta8 = gamma0 - 3 * b

e1 = 7 * e - 2 * m
m1 = 3 * e - 2 * m
r1 = e1 - 2 * m1 * third
a1 = 2 * e1 - 4 * m1 * third
b1 = 2 * e1 - 2 * m1 * third
gamma1 = -8 * e1 + 2 * m1

# The endpoint x^a z^b lies on the input gamma face and is its unique
# largest-z term.  Since the complete largest-z face is a singleton, the
# inequalities gamma>=gamma0 and k<=b make it the unique minimizer of both
# gamma-k and gamma-3k.
require_zero(a - 5 * b - gamma0, "input endpoint gamma")
require_zero(delta6 - (-10 * e + 8 * m * third), "delta6 minimum")
require_zero(delta8 - (-14 * e + 4 * m), "delta8 minimum")
require_zero(h - (-4 * e + 2 * m * third), "endpoint q exponent")

# Exponent ledger for the c-top scaling.
require_zero(2 * e - 2 * h - a1, "c-top A exponent")
require_zero(2 * e + 3 * b - h - b1, "c-top C exponent")
require_zero(6 * e + 3 * (-a + 6 * b) - 3 * b1, "c-top scaling order")

# Exponent ledger for the target-gamma scaling.
require_zero(e + 3 * b - e1, "gamma C exponent")
require_zero(e - h - r1, "gamma D exponent")
require_zero(-8 * e + 3 * (a - 8 * b) - gamma1, "gamma scaling order")

# The endpoint coefficient of the output gamma face is its face scalar
# multiplied by 27^r1.  This is exactly the c-top face scalar below.
require_zero(r1 - (e - h), "output endpoint compatibility exponent")


# --- Exact inverse-chart initial forms ---------------------------------

A, B, C, q = sp.symbols("A B C q")

L = 27 * A**2 * C**2 - 18 * A * B * C + 16 * A + B**3 * C - B**2
T = 4 - 3 * B * C
S = 27 * A * C**2 - 9 * B * C + 8
K = 9 * A * C - B
Y0 = 81 * A * B * C**2 - 72 * A * C - 15 * B**2 * C + 16 * B
A1 = 27 * A * B * C**2 + 54 * A * C - 9 * B**2 * C + 2 * B
A2 = (
    27 * A * B**2 * C**2
    + 18 * A * B * C
    - 48 * A
    - 9 * B**3 * C
    + 10 * B**2
)
Z0 = (
    -2916 * A**3 * C**4
    + 2916 * A**2 * B * C**3
    - 4536 * A**2 * C**2
    + 621 * A * B**3 * C**3
    - 1026 * A * B**2 * C**2
    + 504 * A * B * C
    + 64 * A
    - 207 * B**4 * C**2
    + 454 * B**3 * C
    - 256 * B**2
)
invariants = {"L": L, "T": T, "S": S, "K": K, "Y0": Y0, "A1": A1, "A2": A2, "Z0": Z0}

# (a,b,c)=(A,B,C*s^3), w=q/s, s->infinity.
ctop_expected = {
    "L": (6, 27 * A**2 * C**2),
    "T": (3, -3 * B * C),
    "S": (6, 27 * A * C**2),
    "K": (3, 9 * A * C),
    "Y0": (6, 81 * A * B * C**2),
    "A1": (6, 27 * A * B * C**2),
    "A2": (6, 27 * A * B**2 * C**2),
    "Z0": (12, -2916 * A**3 * C**4),
}
ctop_faces = {
    name: polynomial_face(expression, (0, 0, 3), take_minimum=False)
    for name, expression in invariants.items()
}
require_face_ledger(ctop_faces, ctop_expected, "c-top")

D0 = 27 * A**2 * C
ctop_cubic = D0 * q**3 - 2
ctop_y = -3 * ctop_faces["K"][1] * ctop_faces["L"][1] * q**2 / (
    2 * ctop_faces["S"][1]
)
ctop_z = ctop_faces["Z0"][1] / (8 * ctop_faces["S"][1])
require_zero_mod_cubic(q * ctop_y + 1, ctop_cubic, "c-top y=-1/q")
require_zero_mod_cubic(q**3 * ctop_z + C, ctop_cubic, "c-top z=-C/q^3")

# (a,b,c)=(A*t,B/t,C/t^5), w=q*t, t->0.
D = 27 * A**2 * C + B**3
gamma_expected = {
    "L": (-8, C * D),
    "T": (-6, -3 * B * C),
    "S": (-9, 27 * A * C**2),
    "K": (-4, 9 * A * C),
    "Y0": (-10, 81 * A * B * C**2),
    "A1": (-10, 27 * A * B * C**2),
    "A2": (-11, 27 * A * B**2 * C**2),
    "Z0": (-17, -2916 * A**3 * C**4 + 621 * A * B**3 * C**3),
}
gamma_faces = {
    name: polynomial_face(expression, (1, -1, -5), take_minimum=True)
    for name, expression in invariants.items()
}
require_face_ledger(gamma_faces, gamma_expected, "target-gamma")

gamma_cubic = D * q**3 - 3 * B * q - 2
gamma_y = (
    gamma_faces["Y0"][1]
    - 3 * gamma_faces["K"][1] * gamma_faces["L"][1] * q**2
) / (2 * gamma_faces["S"][1])
gamma_z = (
    gamma_faces["Z0"][1]
    + 6 * gamma_faces["L"][1] * gamma_faces["A1"][1] * q
    - 9 * gamma_faces["L"][1] * gamma_faces["A2"][1] * q**2
) / (8 * gamma_faces["S"][1])
require_zero_mod_cubic(q * gamma_y + 1, gamma_cubic, "target-gamma y=-1/q")
require_zero_mod_cubic(q**3 * gamma_z + C, gamma_cubic, "target-gamma z=-C/q^3")

# A raw resultant with q is 2.  The function-field norm divides by the
# nonmonic leading coefficient once, giving product(q)=2/LC.  This remains
# valid for every integral h, including the negative, non-multiple-of-three
# exponents occurring for L and H.
require_zero(sp.resultant(ctop_cubic, q, q) - 2, "c-top raw resultant")
require_zero(sp.resultant(gamma_cubic, q, q) - 2, "gamma raw resultant")
require_zero(sp.discriminant(ctop_cubic, q) + 78732 * A**4 * C**2, "c-top discriminant")
require_zero(
    sp.discriminant(gamma_cubic, q) + 2916 * A**2 * C * D,
    "gamma discriminant",
)

gamma_q_matrix = sp.Matrix(
    [
        [0, 0, 2 / D],
        [1, 0, 3 * B / D],
        [0, 1, 0],
    ]
)
ctop_q_matrix = gamma_q_matrix.subs({B: 0, D: D0})
require_zero(gamma_q_matrix.det() - 2 / D, "gamma multiplication determinant")
require_zero(ctop_q_matrix.det() - 2 / D0, "c-top multiplication determinant")


# --- Integer hostiles and exact scalar calibration ----------------------

def packet_data(e_value: int, m_value: int) -> dict[str, int]:
    require(e_value >= 0, "negative e")
    require(0 <= m_value <= e_value, "admissible packet requires 0<=m<=e")
    require(m_value % 3 == 0, "packet requires 3|m")
    a_value = 2 * e_value - 4 * m_value // 3
    b_value = 2 * e_value - 2 * m_value // 3
    h_value = a_value - 3 * b_value
    e_next = 7 * e_value - 2 * m_value
    m_next = 3 * e_value - 2 * m_value
    return {
        "e": e_value,
        "m": m_value,
        "a": a_value,
        "b": b_value,
        "h": h_value,
        "e1": e_next,
        "m1": m_next,
        "a1": 2 * e_next - 4 * m_next // 3,
        "b1": 2 * e_next - 2 * m_next // 3,
        "r1": e_next - 2 * m_next // 3,
        "gamma1": -8 * e_next + 2 * m_next,
    }


def renewal_scalars(e_value: int, m_value: int, endpoint_coefficient: Fraction):
    row = packet_data(e_value, m_value)
    gamma_scalar = (
        Fraction((-1) ** row["b"])
        * endpoint_coefficient**3
        * Fraction(2) ** row["h"]
    )
    ctop_scalar = gamma_scalar * Fraction(27) ** row["r1"]
    require(ctop_scalar != 0 and gamma_scalar != 0, "renewal scalar vanished")
    return row, ctop_scalar, gamma_scalar


states = [(0, 0), (1, 0), (3, 3), (7, 3), (43, 15), (271, 99), (1699, 615)]
rows = [packet_data(*state) for state in states]
require(rows[2]["h"] == -10, "m=e boundary hostile changed")
require(rows[1]["h"] % 3 != 0 and rows[3]["h"] % 3 != 0, "root-label hostile vanished")

# P=L, Q=H/2^6.  The exponent h=-4 is not divisible by 3.
_, L_to_H_top, L_to_H_gamma = renewal_scalars(1, 0, Fraction(3**3))
require(L_to_H_gamma == Fraction(2**2 * 3**9, 2**6), "L->H gamma scalar")
require(L_to_H_top == Fraction(2**2 * 3**24, 2**6), "L->H c-top scalar")

# P=H, Q=J/2^35.  The exponent h=-26 is not divisible by 3.
_, H_to_J_top, H_to_J_gamma = renewal_scalars(7, 3, Fraction(2**2 * 3**24))
require(H_to_J_gamma == Fraction(2**15 * 3**72, 2**35), "H->J gamma scalar")
require(H_to_J_top == Fraction(2**15 * 3**171, 2**35), "H->J c-top scalar")

# P=J, Q=G.  This reproduces THM-3513 exactly.
_, G_top, G_gamma = renewal_scalars(43, 15, Fraction(2**15 * 3**171))
require(G_gamma == Fraction(3**513, 2**117), "J->G gamma scalar")
require(G_top == Fraction(3**1128, 2**117), "J->G c-top scalar")

# Fixed consequences once polynomiality is supplied by THM-3506/3521.
R5_row, R5_top, R5_gamma = renewal_scalars(271, 99, G_top)
require(R5_row["e1"] == 1699 and R5_row["m1"] == 615, "R5 packet pair")
require(R5_top == Fraction(3**7251, 2**1369), "R5 c-top scalar")
require(R5_gamma == Fraction(3**3384, 2**1369), "R5 gamma scalar")

R6_row, R6_top, R6_gamma = renewal_scalars(1699, 615, R5_top)
require(R6_row["e1"] == 10663 and R6_row["m1"] == 3867, "R6 packet pair")
require(R6_top == Fraction(3**46008, 2**10493), "R6 c-top scalar")
require(R6_gamma == Fraction(3**21753, 2**10493), "R6 gamma scalar")

semantic_payload = "\n".join(
    [
        "generic:a,b,h,e1,m1,a1,b1,r1,gamma1",
        *[
            ":".join(str(row[key]) for key in ("e", "m", "a", "b", "h", "e1", "m1", "a1", "b1", "r1", "gamma1"))
            for row in rows
        ],
        "G:top=3^1128/2^117:gamma=3^513/2^117",
        "R5:top=3^7251/2^1369:gamma=3^3384/2^1369",
        "R6:top=3^46008/2^10493:gamma=3^21753/2^10493",
    ]
).encode("ascii")
semantic_sha256 = hashlib.sha256(semantic_payload).hexdigest()
require(
    semantic_sha256 == "8b6a447c98e4e7f6bfc493818696d4a9193b4da47ab7b2f9e0368e9155940a91",
    "semantic ledger changed",
)


print("== fixed Keller five-face renewal propagation ==")
print("input endpoint: a=2e-4m/3, b=2e-2m/3, h=a-3b=-4e+2m/3")
print("unique hybrid minima: delta6=gamma-k=-10e+8m/3; delta8=gamma-3k=-14e+4m")
print("c-top residual: 27*A^2*C*q^3-2; Norm(q)=2/(27*A^2*C)")
print("target-gamma residual: D*q^3-3*B*q-2; D=27*A^2*C+B^3; Norm(q)=2/D")
print("raw resultant=2 in both charts; nonmonic leading coefficient division: PASS")
print("exact inverse-coordinate residual reductions and generic separability: PASS")
print("output pair: e'=7e-2m, m'=3e-2m")
print("output c-top: x^(2e'-4m'/3) z^(2e'-2m'/3)")
print("output gamma: z^e' (27x^2z+y^3)^(e'-2m'/3)")
print("endpoint compatibility: c_top_scalar=gamma_scalar*27^(e'-2m'/3)")
print("branchwise root-independence hostile: h(L)=-4 and h(H)=-26 are not divisible by 3")
print("L->H/2^6, H->J/2^35, and J->G exact renewal scalars: PASS")
print(
    "R5 full renewal faces: "
    "c-top=(3^7251/2^1369)*x^2578*z^2988; "
    "gamma=(3^3384/2^1369)*z^1699*(27x^2z+y^3)^1289"
)
print(
    "R6 full renewal faces: "
    "c-top=(3^46008/2^10493)*x^16170*z^18748; "
    "gamma=(3^21753/2^10493)*z^10663*(27x^2z+y^3)^8085"
)
print("polynomiality is the only global hypothesis; no new finite-sheet condition enters renewal")
print(f"semantic_sha256={semantic_sha256}")
print("scope: fixed inverse chart; no fifth image, degree243, future polynomiality, all-level, or general JC claim")
print("all exact checks passed")
