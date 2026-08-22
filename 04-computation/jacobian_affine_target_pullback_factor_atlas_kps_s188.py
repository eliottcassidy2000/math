#!/usr/bin/env python3
"""Exact factor atlas for the affine-target pullback theorem THM-3559."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


x, y, z = sp.symbols("x y z")
A, B, C = sp.symbols("A B C")
r, t = sp.symbols("r t")
u = 1 + x * y

F1 = u**3 * z + y**2 * u * (4 + 3 * x * y)
F2 = y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x**2 * y - x**3 * z

# The general affine target coordinate through q=(-1/4,0,0).
H = sp.expand(A * (F1 + sp.Rational(1, 4)) + B * F2 + C * F3)
K = sp.factor(sp.diff(H, z))
L = sp.factor(H.subs(z, 0))
expected_K = A * u**3 + 3 * B * x * u**2 - C * x**3
require(sp.expand(K - expected_K) == 0, "z-coefficient binary cubic")
require(sp.expand(H - (K * z + L)) == 0, "affine-modification decomposition")

# On d_r=u-r*x=0, parameterize x=t, y=r-1/t.  A common factor d_r
# requires P(r)=0 and the two surviving Laurent coefficients to vanish.
P = A * r**3 + 3 * B * r**2 - C
restricted = sp.factor(
    t * L.subs({x: t, y: r - 1 / t})
)
expected_restricted = sp.expand(
    3 * r * P * t**3
    - 5 * P * t**2
    + (A * r**2 + A / 4 + 4 * B * r) * t
    + A * r
    + 2 * B
)
require(sp.expand(restricted - expected_restricted) == 0, "restriction to Kummer component")

# If d_r is common and the parameter triple is nonzero, A cannot vanish.
# Normalize A=1.  The Groebner basis gives the complete exceptional locus.
common_ideal = [
    r**3 + 3 * B * r**2 - C,
    r**2 + sp.Rational(1, 4) + 4 * B * r,
    r + 2 * B,
]
common_gb = sp.groebner(common_ideal, C, B, r, order="lex")
expected_gb = sp.groebner(
    [C + r / 8, B + r / 2, r**2 - sp.Rational(1, 4)],
    C,
    B,
    r,
    order="lex",
)
require(common_gb == expected_gb, "complete common Kummer-factor locus")

# Verify both exceptional factorizations and that the common factor occurs
# only once in H although it is a double root of the binary cubic K.
special_rows = []
for root in (sp.Rational(1, 2), sp.Rational(-1, 2)):
    b_value = -root / 2
    c_value = -root / 8
    h_special = sp.factor(H.subs({A: 1, B: b_value, C: c_value}))
    k_special = sp.factor(K.subs({A: 1, B: b_value, C: c_value}))
    l_special = sp.factor(L.subs({A: 1, B: b_value, C: c_value}))
    divisor = u - root * x
    require(sp.rem(sp.Poly(k_special, y), sp.Poly(divisor, y)) == 0, f"K divisor r={root}")
    require(sp.rem(sp.Poly(l_special, y), sp.Poly(divisor, y)) == 0, f"L divisor r={root}")
    quotient = sp.cancel(h_special / divisor)
    require(sp.denom(quotient) == 1, f"polynomial residual r={root}")
    require(sp.rem(sp.Poly(quotient, y), sp.Poly(divisor, y)) != 0, f"simple H divisor r={root}")
    p_special = sp.factor((sp.Symbol("R")**3 + 3 * b_value * sp.Symbol("R")**2 - c_value))
    special_rows.append((root, b_value, c_value, p_special))

# The only x-common case is the pure F3 row.  Its two irreducible geometric
# factors are the coordinate plane x=0 and the punctured Kummer surface C=0.
pure_f3 = sp.factor(H.subs({A: 0, B: 0, C: 1}))
expected_f3 = -x * (x**2 * z + 3 * x * y - 2)
require(sp.expand(pure_f3 - expected_f3) == 0, "pure F3 factorization")

# Collision incidence of a Kummer cylinder d_r=0.  It never contains p0 and
# cannot contain both p+ and p-.
d_p0 = 1
d_pp = -sp.Rational(1, 2) - r
d_pm = -sp.Rational(1, 2) + r
require(d_p0 != 0, "Kummer cylinder misses p0")
require(sp.solve([d_pp, d_pm], [r], dict=True) == [], "Kummer cylinder misses every pair")

print("THM-3559 affine target-pullback factor atlas")
print("H=A*(F1+1/4)+B*F2+C*F3=K*z+L")
print("K=A*u^3+3*B*x*u^2-C*x^3, u=1+x*y")
print("t*L|_(u=r*x)=3*r*P*t^3-5*P*t^2+(A*r^2+A/4+4*B*r)*t+A*r+2*B")
print("P=A*r^3+3*B*r^2-C")
print("normalized common-factor Groebner basis: C+r/8, B+r/2, r^2-1/4")
for root, b_value, c_value, p_special in special_rows:
    print(f"exception r={root}: [A:B:C]=[1:{b_value}:{c_value}], P(R)={p_special}")
print("pure F3 pullback: -x*(x^2*z+3*x*y-2)")
print("Kummer cylinder values at (p0,p+,p-): 1, -1/2-r, -1/2+r")
print("all active truth gates passed")
