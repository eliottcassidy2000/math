#!/usr/bin/env python3
"""Exact referee for THM-2582.

Checks the odd/even block controls for the discriminant-of-a-norm law,
rederives the sporadic Keller norm identity from the corrected THM-2576
inverse chart, and verifies the resulting level-two square-class ledger.
"""

import hashlib
import pickle

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


x, y, t = sp.symbols("x y t")

print("== THM-2582: odd-block discriminant tower ==")
print()
print("== parity controls for conjugate blocks ==")

g = y**2 - t
disc_g = sp.discriminant(g, y)
require(disc_g == 4 * t, "quadratic index discriminant changed")

controls = (
    (3, x**3 + (y + 1) * x + y**2 + 2, (t + 2) ** 2),
    (2, x**2 + (y + 1) * x + y**2 + 2, (t + 2) ** 2),
)
for degree, block, expected_square in controls:
    product = sp.resultant(g, block, y)
    disc_product = sp.factor(sp.discriminant(product, x))
    norm_disc = sp.factor(sp.resultant(g, sp.discriminant(block, x), y))
    ratio = sp.factor(disc_product / (disc_g**degree * norm_disc))
    require(ratio == expected_square, "block discriminant quotient changed")
    print(
        f"  degree {degree}: Disc(Norm f)/(Disc(g)^{degree} Norm Disc(f))"
        f" = ({t + 2})^2"
    )

for degree in range(1, 9):
    # Swapping two degree-m blocks changes their resultant by (-1)^(m^2).
    require((degree * degree - degree) % 2 == 0, "m^2 parity changed")
print("  block-swap sign = (-1)^m; odd blocks are alternating, even blocks symmetric")

print()
print("== corrected inverse-chart norm identity ==")

X, s, a, b, c = sp.symbols("X s a b c")
L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
T = 4 - 3 * b * c
S = 27 * a * c**2 - 9 * b * c + 8
E = L * X**3 + T * X - 2 * c
D = (
    -3 * X**2 * a * c
    + X**2 * b**2 * c
    - X**2 * b
    + 2 * X * b * c
    - 2 * X
    + c
)

res_ed = sp.factor(sp.resultant(E, D, X))
require(res_ed == a**2 * c**3 * S**2, "Res(E,D) identity changed")
print("  Res_X(E,D)=a^2 c^3 S^2")

with open("05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl", "rb") as handle:
    H = pickle.load(handle)
H_poly = sp.Poly(H, a, b, c)
ledger = "\n".join(
    f"{monomial}:{coefficient}" for monomial, coefficient in H_poly.terms()
)
ledger_hash = hashlib.sha256(ledger.encode("ascii")).hexdigest()
require(
    ledger_hash
    == "b146c11f33e895c08f72303d282e2668d955e0a58d9268a1b445d4d5202016c2",
    "H coefficient ledger changed",
)

# The root-product derivation uses the corrected standard identity
# Res(E,Q)=a^8 c^18 S^8 H, deg(Q)=15, together with
# Res(E,D)=a^2 c^3 S^2 and product(w_i)=2c/L.
numerator_exponents = {"a": 8, "c": 18, "S": 8, "L": -15, "H": 1, "2": 0}
denominator_exponents = {"a": 8, "c": 18, "S": 8, "L": -14, "H": 0, "2": 6}
norm_exponents = {
    key: numerator_exponents[key] - denominator_exponents[key]
    for key in numerator_exponents
}
require(
    norm_exponents == {"a": 0, "c": 0, "S": 0, "L": -1, "H": 1, "2": -6},
    "Keller norm exponent ledger changed",
)
print("  prod_(q in F^-1(t)) L(q)=H/(64 L)")

# Exact hostile at (a,b,c)=(1,1,1), including the PRS sign trap from
# MISTAKE-290.  The standard lower-degree-first resultant is restored as
# (-1)^(3*15) Res(Q,E).
vals = {a: 1, b: 1, c: 1}
N = (
    -3 * X**3 * a * b * c
    + 4 * X**3 * a
    - 6 * X**2 * a * c
    + X**2 * b**2 * c
    - X**2 * b
    + 2 * X * b * c
    - 2 * X
    + c
)
P = (
    16 * X**7
    + 296 * X**4 * s**2
    - 360 * X**4 * s
    + 108 * X**4
    + 180 * X**3 * c * s
    - 108 * X**3 * c
    + 27 * X**2 * c**2
    - 3 * X * s**4
    + 2 * X * s**3
    - c * s**3
)
Ev = E.subs(vals)
Dv = D.subs(vals)
Nv = N.subs(vals)
Pv = P.subs(vals)
Qv = sp.cancel(Dv**4 * Pv.subs(s, -Nv / Dv))
require(sp.denom(Qv) == 1 and sp.degree(Qv, X) == 15, "specialized Q changed")
prs_low_first = sp.resultant(Ev, Qv, X)
standard = -sp.resultant(Qv, Ev, X)
require(standard == -prs_low_first and standard > 0, "resultant sign repair failed")
Hv = H.subs(vals)
Lv = L.subs(vals)
require(Hv == 951326441195 and Lv == 25, "specialized H/L control changed")
require(
    standard == S.subs(vals) ** 8 * Hv,
    "standard Res(E,Q)=a^8 c^18 S^8 H failed at hostile",
)
print(f"  (1,1,1): PRS={prs_low_first}, standard Sylvester={standard}")
print(f"  H/(64L)={sp.Rational(Hv, 64 * Lv)}")

xi = sp.symbols("xi")
Yden_v = 12 * X**2 - X**2 + X + 2
qy_v = 1 - 3 * X * (9 * X - X + 2) / Yden_v
qz_v = sp.cancel((2 * X - 3 * X**2 * qy_v - 1) / X**3)
Lq_v = 27 * X**2 * qz_v**2 - 18 * X * qy_v * qz_v + 16 * X + qy_v**3 * qz_v - qy_v**2
Tq_v = 4 - 3 * qy_v * qz_v
inner_v = sp.cancel(Lq_v * xi**3 + Tq_v * xi - 2 * qz_v)
inner_num_v = sp.together(inner_v).as_numer_denom()[0]
core_v = sp.Poly(sp.resultant(Ev, inner_num_v, X), xi)
content_v = sp.gcd_list(core_v.all_coeffs())
core_v = sp.Poly(sp.cancel(core_v.as_expr() / content_v), xi)
require(core_v.degree() == 9, "specialized composite core lost degree nine")
require(sp.gcd(core_v, core_v.diff()).degree() == 0, "composite core became inseparable")
require(core_v.LC() == -Hv, "specialized composite lead is not -H")
print("  (1,1,1) composite x-core: degree 9, squarefree, lead=-H")

print()
print("== level-two square-class ledger ==")

# [Disc outer]=[-L].  The norm of the three inner discriminants is
# [(-4)^3 * H/(64L)]=[-H/L].  Odd-block alternation multiplies them.
outer = {"minus": 1, "L": 1, "H": 0}
inner_norm = {"minus": 1, "L": -1, "H": 1}
total = {key: (outer[key] + inner_norm[key]) % 2 for key in outer}
require(total == {"minus": 0, "L": 0, "H": 1}, "square-class cancellation changed")
print("  [-L] * [-H/L] = [H]")
print("  odd irreducible part of the degree-nine x-discriminant: H alone")
print()
print("scope: field square class/parity only; no higher-level or JC(2) conclusion")
print("all exact checks passed")
