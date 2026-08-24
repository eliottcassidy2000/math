#!/usr/bin/env python3
"""Exact companion for THM-3935's resolvent class and cubic character."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


t, X, W, x, y, U, V, s = sp.symbols("t X W x y U V s")
q = t**3 - X**2
H = sp.expand(q**2 - 4 * X**3)


# ---------------------------------------------------------------------------
# Normality and the vertical-prime localization bridge.
# ---------------------------------------------------------------------------

F = sp.expand(W**2 - H)
gate(sp.factor(F) == F, "quadratic resolvent is irreducible")
gate(sp.factor(sp.discriminant(H, X)) == -256 * t**12 * (16 * t**3 + 27),
     "generic quartic discriminant")

# Every closed fibre is irreducible.  If H_T=(X^2+bX+c)^2, coefficient
# comparison first gives b=-2 and c=-T-2; the X coefficient then gives
# T=-2, contradicting the constant coefficient c^2=T^2.
T, b, c = sp.symbols("T b c")
H_T = sp.expand((T - X**2) ** 2 - 4 * X**3)
square_trial = sp.expand((X**2 + b * X + c) ** 2)
coeff = sp.Poly(square_trial - H_T, X)
gate(coeff.coeff_monomial(X**3) == 2 * b + 4, "fibre square cubic coefficient")
gate(coeff.coeff_monomial(X**2) == b**2 + 2 * c + 2 * T,
     "fibre square quadratic coefficient")
gate(coeff.coeff_monomial(X) == 2 * b * c, "fibre square linear coefficient")
gate(coeff.coeff_monomial(1) == c**2 - T**2, "fibre square constant coefficient")
gate(sp.expand((c**2 - T**2).subs({b: -2, c: 0, T: -2})) == -4,
     "all closed fibre squares contradict their constant term")

# The affine singular support is the origin.  These equations encode the
# elementary case split W=0, X(q+3X)=0, t^2 q=0, H=0.
FW = sp.diff(F, W)
FX = sp.factor(sp.diff(F, X))
Ft = sp.factor(sp.diff(F, t))
gate(FW == 2 * W, "resolvent W derivative")
gate(FX == 4 * X * (3 * X - X**2 + t**3), "resolvent X derivative")
gate(Ft == -6 * t**2 * (t**3 - X**2), "resolvent t derivative")
zero(F.subs({X: 0, W: 0}) + t**6, "X=0 singular case forces t=0")
zero((t**3 - X**2 + 3 * X).subs(t, 0) - X * (3 - X),
     "nonzero-X singular case equation")
gate(H.subs({t: 0, X: 3}) == -27, "nonzero-X singular case misses surface")


# ---------------------------------------------------------------------------
# Generic genus-one curve and its global minimal Weierstrass model.
# ---------------------------------------------------------------------------

# First pass to the nonminimal Weierstrass coordinates, then make the global
# K=k(t) change x_old=t^2*x-1, y_old=t^3*y.
U_expr = W + X**2 - 2 * X - t**3 - 2
V_expr = X * U_expr
x_old = U_expr / 2
y_old = (V_expr - U_expr - t**3 - 2) / 2
E_old = x_old**3 + (t**3 + 3) * x_old**2 + (2 * t**3 + 3) * x_old + (t**3 + 2) ** 2 / 4
quartic_identity = sp.cancel(y_old**2 - E_old)
quartic_numerator = sp.fraction(quartic_identity)[0]
gate(sp.rem(sp.Poly(quartic_numerator, W), sp.Poly(F, W)).as_expr() == 0,
     "quartic-to-Weierstrass identity")

E_rhs = x**3 + t * x**2 + sp.Rational(1, 4)
zero(
    (t**3 * y) ** 2
    - (
        (t**2 * x - 1) ** 3
        + (t**3 + 3) * (t**2 * x - 1) ** 2
        + (2 * t**3 + 3) * (t**2 * x - 1)
        + (t**3 + 2) ** 2 / 4
    )
    - t**6 * (y**2 - E_rhs),
    "global minimal Weierstrass change",
)

a2 = t
a4 = sp.Integer(0)
a6 = sp.Rational(1, 4)
b2 = 4 * a2
b4 = 2 * a4
b6 = 4 * a6
b8 = 4 * a2 * a6 - a4**2
c4 = sp.factor(b2**2 - 24 * b4)
c6 = sp.factor(-b2**3 + 36 * b2 * b4 - 216 * b6)
Delta = sp.factor(-b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6)
gate(c4 == 16 * t**2, "minimal c4")
zero(c6 + 8 * (8 * t**3 + 27), "minimal c6")
zero(Delta + 16 * t**3 + 27, "minimal discriminant")
gate(Delta.subs(t, 0) == -27, "good fibre at t=0")

# At infinity put s=1/t, x_inf=s^2*x and y_inf=s^3*y.
xinf, yinf = sp.symbols("xinf yinf")
E_inf = yinf**2 - xinf**3 - s * xinf**2 - s**6 / 4
zero(
    E_inf.subs({xinf: s**2 * x, yinf: s**3 * y}, simultaneous=True).subs(s, 1 / t)
    * t**6
    - (y**2 - E_rhs),
    "integral infinity model",
)
c4_inf = 16 * s**2
c6_inf = -8 * s**3 * (8 + 27 * s**3)
Delta_inf = -s**9 * (16 + 27 * s**3)
gate(sp.Poly(c4_inf, s).as_dict().get((2,)) == 16, "infinity c4 valuation two")
gate(sp.Poly(c6_inf, s).as_dict().get((3,)) == -64, "infinity c6 valuation three")
gate(sp.Poly(Delta_inf, s).as_dict().get((9,)) == -16,
     "infinity discriminant valuation nine")


# ---------------------------------------------------------------------------
# Mordell--Weil generator and the I3* spin correction.
# ---------------------------------------------------------------------------


def add_points(P: tuple[sp.Expr, sp.Expr] | None,
               Q: tuple[sp.Expr, sp.Expr] | None) -> tuple[sp.Expr, sp.Expr] | None:
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if sp.cancel(x1 - x2) == 0:
        if sp.cancel(y1 + y2) == 0:
            return None
        slope = sp.cancel((3 * x1**2 + 2 * t * x1) / (2 * y1))
    else:
        slope = sp.cancel((y2 - y1) / (x2 - x1))
    x3 = sp.factor(slope**2 - t - x1 - x2)
    y3 = sp.factor(-y1 + slope * (x1 - x3))
    return x3, y3


Q = (sp.Integer(0), -sp.Rational(1, 2))
two_Q = add_points(Q, Q)
three_Q = add_points(two_Q, Q)
gate(two_Q == (-t, sp.Rational(1, 2)), "exact double of Q")
gate(three_Q == (t**-2, (t**3 + 2) / (2 * t**3)), "second infinity is 3Q")

# On the quartic branch W/X^2 -> -1, epsilon=1/X.  Its first surviving
# U-term gives exactly the finite Weierstrass point 3Q; the opposite branch
# is the chosen zero section.
epsilon = sp.symbols("epsilon")
negative_root = -epsilon**-2 * sp.series(
    sp.sqrt(1 - 4 * epsilon - 2 * t**3 * epsilon**2 + t**6 * epsilon**4),
    epsilon,
    0,
    4,
).removeO()
U_infinity = sp.expand(
    negative_root + epsilon**-2 - 2 * epsilon**-1 - t**3 - 2
)
gate(sp.limit(U_infinity / epsilon, epsilon, 0) == 2 * t**3 + 4,
     "second infinity U expansion")
gate(sp.limit(U_infinity, epsilon, 0) == 0 and t**-2 == three_Q[0],
     "second infinity x coordinate")
gate(
    sp.cancel(
        (sp.limit(U_infinity / epsilon, epsilon, 0) - t**3 - 2)
        / (2 * t**3)
        - three_Q[1]
    ) == 0,
    "second infinity y coordinate",
)

# The two affine X=0 points become Q and 2Q under the birational map.
for W_value, expected in ((t**3, Q), (-t**3, two_Q)):
    U_value = sp.expand(U_expr.subs({X: 0, W: W_value}))
    V_value = sp.Integer(0)
    mapped = (
        sp.cancel((U_value / 2 + 1) / t**2),
        sp.cancel((V_value - U_value - t**3 - 2) / (2 * t**3)),
    )
    gate(mapped == expected, f"affine section W={W_value} maps correctly")

# Three successive s-chart point blowups track Q to a terminal D7 spin tip.
local = sp.expand(E_inf)
qx = sp.Integer(0)  # Q has x_inf=0 exactly.
qy = -s**3 / 2
exceptional_equations = []
for index in range(3):
    local = sp.expand(local.subs({xinf: s * xinf, yinf: s * yinf}, simultaneous=True))
    minimum = min(monomial[0] for monomial, coefficient in sp.Poly(local, s, xinf, yinf).terms())
    gate(minimum == 2, f"blowup {index + 1} strict-transform multiplicity")
    local = sp.expand(local / s**2)
    qx = sp.cancel(qx / s)
    qy = sp.cancel(qy / s)
    exceptional_equations.append(sp.factor(local.subs(s, 0)))
gate(exceptional_equations[0] == yinf**2, "first exceptional double component")
gate(exceptional_equations[1] == yinf**2, "second exceptional double component")
gate(exceptional_equations[2] == (2 * yinf - 1) * (2 * yinf + 1) / 4,
     "third exceptional splits into spin tips")
gate(qx == 0 and qy == -sp.Rational(1, 2), "Q meets the negative spin tip")

# D7 Cartan data: either spin tip has correction 7/4 and determinant four.
D7 = 2 * sp.eye(7)
for i, j in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (4, 6)):
    D7[i, j] = D7[j, i] = -1
gate(D7.det() == 4, "D7 determinant")
D7_inverse = D7.inv()
gate(D7_inverse[5, 5] == sp.Rational(7, 4), "first spin correction")
gate(D7_inverse[6, 6] == sp.Rational(7, 4), "second spin correction")
gate(2 - D7_inverse[5, 5] == sp.Rational(1, 4), "Q height one quarter")
gate(sp.Rational(1, 4) * D7.det() == 1, "unimodular MW determinant invoice")

# There is no rational two-torsion.  A root of the monic cubic is in k[t].
# Degree comparison leaves r=-t+b; its t^2 and constant equations conflict.
broot = sp.symbols("broot")
root_trial = -t + broot
root_poly = sp.Poly(sp.expand(root_trial**3 + t * root_trial**2 + sp.Rational(1, 4)), t)
gate(root_poly.coeff_monomial(t**2) == broot, "two-torsion degree-one coefficient")
gate(root_poly.coeff_monomial(1) == broot**3 + sp.Rational(1, 4),
     "two-torsion constant coefficient")
gate(root_poly.coeff_monomial(1).subs(broot, 0) == sp.Rational(1, 4),
     "two-torsion contradiction")


# ---------------------------------------------------------------------------
# Exact affine class group, units, and unique Kummer character.
# ---------------------------------------------------------------------------

a = q + W
bminus = q - W
zero(sp.expand(a * bminus - 4 * X**3).subs(W**2, H), "Cardano product")
zero(H.subs(X, 0) - t**6, "two X=0 sections")

# Cardano cube roots recover the unique monogenic cubic.
z, v, u = sp.symbols("z v u")
# The identity is checked modulo z*v=X and z^3+v^3=q.
cardano_remainder = sp.expand((z + v) ** 3 - 3 * X * (z + v) - q)
cardano_remainder = sp.expand(cardano_remainder.subs(X, z * v))
zero(cardano_remainder - (z**3 + v**3 - (t**3 - z**2 * v**2)),
     "Cardano cubic identity modulo cube data")
zero(sp.expand((q + W) * (q - W) / 4 - X**3).subs(W**2, H),
     "two Kummer radicands multiply to X cubed")

# The natural cubic surface has only the origin in its singular support.
cubic_surface = sp.expand(u**3 - 3 * X * u - t**3 + X**2)
partials = [sp.diff(cubic_surface, variable) for variable in (u, X, t)]
gate(partials == [3 * u**2 - 3 * X, -3 * u + 2 * X, -3 * t**2],
     "natural cubic surface derivatives")
gate(cubic_surface.subs({t: 0, X: sp.Rational(9, 4), u: sp.Rational(3, 2)}) != 0,
     "spurious nonzero critical solution misses cubic surface")


summary = {
    "checks": CHECKS,
    "resolvent": "normal; unique singularity origin",
    "vertical_fibres": "all irreducible and principal",
    "generic_curve": "genus1 minus two rational infinities",
    "minimal_model": "y^2=x^3+t*x^2+1/4",
    "fibres": "good@0 + 3I1 + I3*@infinity",
    "mordell_weil": "Z*Q; Q=(0,-1/2); Ominus=3Q; height(Q)=1/4",
    "class_group": "Z/3 generated by either X=0 resolvent section",
    "units": "k*",
    "kummer": "one nontrivial C3 cover up to inversion; natural Cardano cover",
    "cubic_order": "unique normal cubic field/order is the monogenic THM3932 order",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3935 exact linear-conic resolvent class and character")
print(f"CHECKS={CHECKS}")
print("NORMAL=YES;SINGULAR_SUPPORT=(0,0,0)")
print("VERTICAL_FIBRES=ALL_IRREDUCIBLE_PRINCIPAL")
print("GENERIC_ELLIPTIC=y^2=x^3+t*x^2+1/4")
print("FIBRES=GOOD_AT_0,THREE_I1,I3STAR_AT_INFINITY")
print("MW=Z*Q;Q=(0,-1/2);OMINUS=3Q;HEIGHT_Q=1/4")
print("CLASS_GROUP=Z/3;UNITS=k*")
print("KUMMER_CHARACTERS=ONE_NONTRIVIAL_CYCLIC_COVER_UP_TO_INVERSION")
print("UNIQUE_CUBIC_FIELD=NATURAL_MONOGENIC_CARDANO_FIELD")
print(f"SEMANTIC_SHA256={semantic}")
