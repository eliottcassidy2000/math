#!/usr/bin/env python3
"""Exact companion for THM-3929's linear-color one-place obstruction.

Reproduction:
  python3 04-computation/jc2_root_regular_one_place_linear_color_cubic_thm3929.py
  python3 -O 04-computation/jc2_root_regular_one_place_linear_color_cubic_thm3929.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(sp.together(expression)) == 0, message)


def binary_discriminant(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(b**2 * c**2 - 4 * a * c**3 - 4 * b**3 * d - 27 * a**2 * d**2 + 18 * a * b * c * d)


A, C, T, U, V, s = sp.symbols("A C T U V s")
a, c, d = sp.symbols("a c d", nonzero=True)


# ---------------------------------------------------------------------------
# Universal repeated-root elimination.
# ---------------------------------------------------------------------------

f = a * T**3 + C * T**2 + c * T + d
df = sp.diff(f, T)
zero(2 * f - T * df - (-a * T**3 + c * T + 2 * d), "repeated-root elimination")
zero(-(2 * f - T * df) - (a * T**3 - c * T - 2 * d), "incidence sign")


# ---------------------------------------------------------------------------
# The explicit t=s coordinate division.
# ---------------------------------------------------------------------------

alpha, beta, gamma, delta = sp.symbols("alpha beta gamma delta", nonzero=True)
aa, cc, dd = sp.symbols("aa cc dd")

# Reduce aa*s^3-cc*s-2dd using A=alpha*s^3+beta*s^2+gamma*s+delta.
reduced = sp.expand(
    aa * (A - beta * s**2 - gamma * s - delta) / alpha - cc * s - 2 * dd
)
poly_s = sp.Poly(reduced, s)
zero(poly_s.coeff_monomial(s**2) + aa * beta / alpha, "s^2 bucket")
zero(poly_s.coeff_monomial(s) + aa * gamma / alpha + cc, "s bucket")
zero(poly_s.coeff_monomial(1) - aa * (A - delta) / alpha + 2 * dd, "constant bucket")

# Insert the forced beta,c,d values and verify both repeated-root equations.
r0, r1, r2 = sp.symbols("r0 r1 r2")
aA = r0 + r1 * (A - delta) + r2 * (A - delta) ** 2
A_s = alpha * s**3 + gamma * s + delta
cA = -gamma * aA / alpha
dA = aA * (A - delta) / (2 * alpha)
a_s = sp.expand(aA.subs(A, A_s))
c_s = sp.expand(cA.subs(A, A_s))
d_s = sp.expand(dA.subs(A, A_s))
C_s = sp.factor(-a_s * (3 * s**2 - gamma / alpha) / (2 * s))
f_s = a_s * s**3 + C_s * s**2 + c_s * s + d_s
df_s = 3 * a_s * s**2 + 2 * C_s * s + c_s
zero(f_s, "coordinate repeated-root equation")
zero(df_s, "coordinate derivative equation")

pole_numerator = sp.expand(a_s * (3 * s**2 - gamma / alpha))
zero(pole_numerator.subs(s, 0) + r0 * gamma / alpha, "s=0 pole residue")
gate(sp.Poly(pole_numerator, s).degree() >= 0, "pole numerator is polynomial")

# Unit/monogenic q=0 positive control.
A_pos = s**3
t_pos = s
C_pos = -3 * s
a_pos, c_pos, d_pos = sp.Integer(2), sp.Integer(0), A
zero((a_pos * T**3 - c_pos * T - 2 * d_pos).subs({A: A_pos, T: t_pos}),
     "monogenic positive incidence")
zero((3 * a_pos * T**2 + 2 * C * T + c_pos).subs({A: A_pos, C: C_pos, T: t_pos}),
     "monogenic positive derivative")
gate(a_pos == 2, "positive index form represents scalar unit")

# Root-regular but nonunit coefficient-ideal control.
A_bad = s**3
C_bad = -3 * s**4
a_bad, c_bad, d_bad = 2 * A, sp.Integer(0), A**2
zero((a_bad * T**3 - c_bad * T - 2 * d_bad).subs({A: A_bad, T: s}),
     "nonunit control incidence")
zero((3 * a_bad * T**2 + 2 * C * T + c_bad).subs({A: A_bad, C: C_bad, T: s}),
     "nonunit control derivative")
gate(all(sp.expand(expr.subs({A: 0, C: 0})) == 0 for expr in (a_bad, C, c_bad, d_bad)),
     "nonunit control has common coefficient address")

# A regular repeated root need not itself be the normalization coordinate.
A_nc = s**3
t_nc = s**2
C_nc = -3 * s**2
a_nc, c_nc, d_nc = sp.Integer(2), sp.Integer(0), A**2
zero((a_nc * T**3 - c_nc * T - 2 * d_nc).subs({A: A_nc, T: t_nc}),
     "noncoordinate regular incidence")
zero((3 * a_nc * T**2 + 2 * C * T + c_nc).subs({A: A_nc, C: C_nc, T: t_nc}),
     "noncoordinate regular derivative")
zero(A_nc / t_nc - s, "A and t recover the normalization coordinate")


# ---------------------------------------------------------------------------
# Centered-Mobius first finite-pole seam.
# ---------------------------------------------------------------------------

u, p, q, r, x, z = sp.symbols("u p q r x z")
w = z + q / (3 * x)
raw_w_equation = x * w**3 - q * w**2 - p * w - 1
centered = sp.factor(27 * x**2 * raw_w_equation)
centered_expected = (
    27 * x**3 * z**3
    - 9 * x * (3 * p * x + q**2) * z
    - (27 * x**2 + 9 * p * q * x + 2 * q**3)
)
zero(centered - centered_expected, "centered Mobius incidence")

x_u = u * (u**2 + p * u + q)
t_u = 1 / u - q / (3 * x_u)
a_m = 27 * x_u**3
c_m = 9 * x_u * (3 * p * x_u + q**2)
C_m = sp.factor(-(3 * a_m * t_u**2 + c_m) / (2 * t_u))
D = 3 * u**2 + 3 * p * u + 2 * q
E = 9 * u**4 + 21 * p * u**3 + (12 * p**2 + 12 * q) * u**2 + 15 * p * q * u + 5 * q**2
C_expected = -27 * u**2 * (u**2 + p * u + q) ** 2 * E / (2 * D)
zero(C_m - C_expected, "centered Mobius color pullback")

remainder = sp.rem(E, D, u)
zero(remainder - q * (p * u + q), "Mobius denominator remainder")
res_prefactor = sp.factor(sp.resultant(D, u * (u**2 + p * u + q), u))
zero(res_prefactor - 2 * q**3, "Mobius denominator/prefactor resultant")
res_E = sp.factor(sp.resultant(D, E, u))
zero(res_E + 27 * q**3 * (p**2 - 3 * q), "Mobius exceptional resultant")
gate(sp.Poly(q * (p * u + q), u).degree() < sp.Poly(D, u).degree(),
     "nonzero remainder is too small for full denominator cancellation")

# At p^2=3q one, but not both, denominator roots can cancel.
q_exception = p**2 / 3
D_exception = sp.factor(D.subs(q, q_exception))
E_exception = sp.factor(E.subs(q, q_exception))
gcd_exception = sp.gcd(sp.Poly(D_exception, u), sp.Poly(E_exception, u)).as_expr()
gate(sp.Poly(gcd_exception, u).degree() == 1, "exception cancels exactly one root")
gate(sp.Poly(D_exception, u).degree() == 2, "exception denominator remains quadratic")

# q=0 primitive content collapses to x*z^3-p*z-1 and a scalar endpoint.
q0_centered = sp.factor(centered_expected.subs(q, 0) / (27 * x**2))
zero(q0_centered - (x * z**3 - p * z - 1), "q=0 primitive incidence")
x_q0 = u**2 * (u + p)
C_q0 = sp.factor(-(3 * x_q0 * (1 / u) ** 2 + p) / (2 / u))
zero(C_q0 + u * (3 * u + 4 * p) / 2, "q=0 polynomial color")


# ---------------------------------------------------------------------------
# THM-3927 higher-degree finite-root-pole hostile and its one-place compression.
# ---------------------------------------------------------------------------

u = sp.symbols("u")
A_h = -(u - 1) ** 2 * (u + 2) / sp.Integer(108)
t_h = -4 * (u + 2) / ((u - 1) * (u + 1))
C_h = (u - 1) * (u + 1) * (u**4 + 6 * u**3 - 22 * u - 21) / (sp.Integer(72) * (u + 2))
a_h = A_h * (1 + 27 * A_h)
c_h = -1 - 24 * A_h
d_h = 8 * A_h
zero(a_h * t_h**3 - c_h * t_h - 2 * d_h, "higher-pole hostile incidence")
zero(a_h * t_h**3 + C_h * t_h**2 + c_h * t_h + d_h, "higher-pole hostile cubic")
zero(3 * a_h * t_h**2 + 2 * C_h * t_h + c_h, "higher-pole hostile derivative")

t_num, t_den = (sp.factor(part) for part in sp.fraction(sp.cancel(t_h)))
gate(sp.gcd(sp.Poly(t_num, u), sp.Poly(t_den, u)) == 1, "hostile root fraction is reduced")
gate(sp.factor(t_den) == (u - 1) * (u + 1), "hostile has two finite root poles")

B_h = sp.factor(A_h * C_h)
gate(not sp.denom(B_h).has(u), "compressed color is polynomial over Q")
gate(sp.Poly(A_h, u).degree() == 3, "compressed first coordinate degree three")
gate(sp.Poly(B_h, u).degree() == 8, "compressed color degree eight")

# Original coefficient ideal is primitive; scaled compression has a common address.
A0 = sp.symbols("A0")
gcd_original = sp.gcd(
    sp.Poly(A0 * (1 + 27 * A0), A0),
    sp.gcd(sp.Poly(-1 - 24 * A0, A0), sp.Poly(8 * A0, A0)),
)
gate(gcd_original.degree() == 0, "hostile original coefficient triple is primitive")
scaled_coefficients = (A0**2 * (1 + 27 * A0), C, A0 * (-1 - 24 * A0), 8 * A0**2)
gate(all(expr.subs({A0: 0, C: 0}) == 0 for expr in scaled_coefficients),
     "one-place scaling creates a common coefficient address")

Delta_h = binary_discriminant(a_h, C_h, c_h, d_h)
zero(Delta_h, "hostile parametrizes its discriminant")


summary = {
    "checks": CHECKS,
    "elimination": "2f-tf'=-a*t^3+c*t+2d",
    "root_regular_mechanism": "t integral over k[A] forces a|c,d; unit coefficient ideal forces a scalar",
    "coordinate_fold": "beta=gamma=0; c=0; d=a(A-delta)/(2alpha); C=-3as/2",
    "mobius_remainder": str(remainder),
    "mobius_resultant": str(res_E),
    "mobius_survivor": "q=0 gives x*t^3-p*t-1 and scalar endpoint d",
    "higher_pole_hostile": {
        "root_denominator": str(t_den),
        "compressed_degrees": [sp.Poly(A_h, u).degree(), sp.Poly(B_h, u).degree()],
    },
    "scope": "degree-at-least-two homogeneous root maps remain open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3929 exact root-regular and centered-Mobius companion")
print(f"CHECKS={CHECKS}")
print("ROOT_REGULAR=a divides c,d; unit coefficient ideal forces scalar a")
print(f"MOBIUS_REMAINDER={remainder}")
print(f"MOBIUS_RESULTANT={res_E}")
print("MOBIUS_SURVIVOR=q=0;primitive=x*t^3-p*t-1;scalar_endpoint=d")
print(f"HIGHER_POLE_DENOMINATOR={t_den}")
print(f"COMPRESSED_DEGREES=A:{sp.Poly(A_h, u).degree()},B:{sp.Poly(B_h, u).degree()}")
print("OPEN=homogeneous repeated-root maps of degree at least two")
print(f"SEMANTIC_SHA256={semantic}")
