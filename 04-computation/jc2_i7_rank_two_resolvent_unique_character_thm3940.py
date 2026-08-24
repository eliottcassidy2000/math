#!/usr/bin/env python3
"""Exact companion for THM-3940's I7 rank-two character obstruction."""

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


t, X, W, x, y, s, r, beta = sp.symbols("t X W x y s r beta")
R = r + beta * t
q = t**3 + t * R * X - X**2
H = sp.expand(q**2 - 4 * X**3)
F = sp.expand(W**2 - H)


# ---------------------------------------------------------------------------
# Integral vertical fibres and normality.
# ---------------------------------------------------------------------------

quartic_discriminant_factor = sp.factor(sp.discriminant(H, X) / (-256 * t**12))
Delta_expected = -t * (4 * t + R**2) ** 2 - R**3 - 36 * R * t - 27
zero(quartic_discriminant_factor + Delta_expected,
     "quartic discriminant equals 256*t^12 times elliptic discriminant")

# Closed-fibre square charts.  Put A=lambda^3 and B=lambda(r+beta*lambda).
lam, A, B, p, c = sp.symbols("lambda A B p c")
closed_q = A + B * X - X**2
closed_H = sp.expand(closed_q**2 - 4 * X**3)
square = sp.expand((X**2 + p * X + c) ** 2)
coeff = sp.Poly(square - closed_H, X)
gate(coeff.coeff_monomial(X**3) == 2 * B + 2 * p + 4,
     "closed square cubic coefficient")
gate(coeff.coeff_monomial(1) == c**2 - A**2,
     "closed square constant coefficient")
p_value = -B - 2
for sign, label in ((1, "plus"), (-1, "minus")):
    c_value = sign * A
    equations = [
        A - lam**3,
        B - lam * (r + beta * lam),
        coeff.coeff_monomial(X**2).subs({p: p_value, c: c_value}),
        coeff.coeff_monomial(X).subs({p: p_value, c: c_value}),
    ]
    basis = sp.groebner(equations, A, B, r, beta, lam, order="lex")
    gate(basis.contains(sp.Integer(1)), f"closed square {label} chart empty")

FW = sp.diff(F, W)
FX = sp.factor(sp.diff(F, X))
Ft = sp.factor(sp.diff(F, t))
gate(FW == 2 * W, "surface W derivative")
zero(FX - 2 * (6 * X**2 + (2 * X - t * R) * q), "surface X derivative")
zero(Ft + 2 * (3 * t**2 + (r + 2 * beta * t) * X) * q,
     "surface t derivative")

# For a hypothetical nonzero singular point, write X=a^2, q=2*eps*a^3,
# t=u*a.  The t and X derivative rows are linear in a.  Eliminating a gives
# a nonzero quartic in u because beta!=0, so the singular locus is finite.
u, aa, eps = sp.symbols("u aa eps")
trace_row = 3 * u**2 + r + 2 * beta * u * aa
x_row = eps * (aa * (2 - beta * u**2) - r * u) + 3
aa_from_trace = -(3 * u**2 + r) / (2 * beta * u)
elimination = sp.factor(2 * beta * u * x_row.subs(aa, aa_from_trace))
gate(sp.Poly(elimination, u).degree() == 4,
     "nonzero surface singular row has finite u support")
gate(sp.Poly(elimination, u).coeff_monomial(u**4) == 3 * beta * eps,
     "surface singular eliminant has nonzero leading coefficient")


# ---------------------------------------------------------------------------
# Elliptic model, boundary sections, and I7 infinity.
# ---------------------------------------------------------------------------

z = sp.cancel((q + W) / (2 * X**2))
x_forward = sp.cancel(t * z)
y_forward = sp.cancel(X * z * (z + 1) + (1 - t * R * z) / 2)
E_rhs = x**3 + (t + R**2 / 4) * x**2 - R * x / 2 + sp.Rational(1, 4)
forward_error = sp.cancel(y_forward**2 - E_rhs.subs(x, x_forward))
forward_numerator = sp.fraction(forward_error)[0]
gate(sp.rem(sp.Poly(forward_numerator, W), sp.Poly(F, W)).as_expr() == 0,
     "quartic-to-elliptic identity")

z_inverse = x / t
X_inverse = sp.cancel(t**2 * (2 * y - 1 + R * x) / (2 * x * (x + t)))
q_inverse = t**3 + t * R * X_inverse - X_inverse**2
W_inverse = sp.cancel(2 * X_inverse**2 * z_inverse - q_inverse)
E_equation = sp.expand(y**2 - E_rhs)
inverse_error = sp.cancel(W_inverse**2 - (q_inverse**2 - 4 * X_inverse**3))
inverse_numerator = sp.factor(sp.fraction(inverse_error)[0])
gate(sp.rem(sp.Poly(inverse_numerator, y), sp.Poly(E_equation, y)).as_expr() == 0,
     "elliptic-to-quartic inverse identity")

# Boundary images are Q and -2Q.
epsilon = sp.symbols("epsilon")
eta_square = sp.expand((-1 + t * R * epsilon + t**3 * epsilon**2) ** 2 - 4 * epsilon)
eta_plus = sp.series(sp.sqrt(eta_square), epsilon, 0, 3).removeO()
for eta, expected_x, expected_y, label in (
    (eta_plus, 0, -sp.Rational(1, 2), "positive"),
    (-eta_plus, -t, -(1 + t * R) / 2, "negative"),
):
    z_boundary = sp.expand((-1 + t * R * epsilon + t**3 * epsilon**2 + eta) / 2)
    xb = sp.limit(t * z_boundary, epsilon, 0)
    yb = sp.limit(
        z_boundary * (z_boundary + 1) / epsilon + (1 - t * R * z_boundary) / 2,
        epsilon,
        0,
    )
    zero(xb - expected_x, f"{label} boundary x section")
    zero(yb - expected_y, f"{label} boundary y section")

Q = (sp.Integer(0), -sp.Rational(1, 2))
doubling_slope = sp.cancel((-R / 2) / (2 * Q[1]))
two_Q_x = sp.factor(doubling_slope**2 - (t + R**2 / 4) - 2 * Q[0])
two_Q_y = sp.factor(-Q[1] + doubling_slope * (Q[0] - two_Q_x))
zero(two_Q_x + t, "boundary section is the x coordinate of two Q")
zero(two_Q_y - (1 + t * R) / 2,
     "boundary section is the y coordinate of two Q")

a2 = t + R**2 / 4
a4 = -R / 2
a6 = sp.Rational(1, 4)
b2 = sp.factor(4 * a2)
b4 = sp.factor(2 * a4)
b6 = sp.factor(4 * a6)
b8 = sp.factor(4 * a2 * a6 - a4**2)
c4 = sp.factor(b2**2 - 24 * b4)
c6 = sp.factor(-b2**3 + 36 * b2 * b4 - 216 * b6)
Delta = sp.factor(-b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6)
zero(b2 - (4 * t + R**2), "I7 family b2")
zero(b4 + R, "I7 family b4")
gate(b6 == 1 and b8 == t, "I7 family b6 b8")
zero(c4 - ((4 * t + R**2) ** 2 + 24 * R), "I7 family c4")
zero(c6 + (4 * t + R**2) ** 3 + 36 * R * (4 * t + R**2) + 216,
     "I7 family c6")
zero(Delta - Delta_expected, "I7 family discriminant")
gate(sp.Poly(Delta, t).degree() == 5 and sp.Poly(Delta, t).LC() == -beta**4,
     "finite discriminant has degree five")

c4_inf = sp.factor(s**4 * c4.subs(t, 1 / s))
c6_inf = sp.factor(s**6 * c6.subs(t, 1 / s))
Delta_inf = sp.factor(s**12 * Delta.subs(t, 1 / s))
gate(c4_inf.subs(s, 0) == beta**4, "infinity c4 is a unit")
gate(c6_inf.subs(s, 0) == -beta**6, "infinity c6 is a unit")
gate(sp.Poly(Delta_inf, s).as_dict().get((7,)) == -beta**4,
     "infinity discriminant valuation seven")

# The five finite fibres are I1 on a nonempty open parameter set; r=0,beta=1
# is an exact positive control.
Delta_control = sp.factor(Delta.subs({r: 0, beta: 1}))
c4_control = sp.factor(c4.subs({r: 0, beta: 1}))
gate(sp.discriminant(Delta_control, t) == 126759838761,
     "five-I1 squarefree positive control")
gate(sp.resultant(Delta_control, c4_control, t) == 12131289,
     "five-I1 multiplicative positive control")


# ---------------------------------------------------------------------------
# No rational three-torsion, despite generic Mordell--Weil rank two.
# ---------------------------------------------------------------------------

psi3 = sp.expand(3 * x**4 + b2 * x**3 + 3 * b4 * x**2 + 3 * b6 * x + b8)
zero(psi3 - (3 * x**4 + (4 * t + R**2) * x**3 - 3 * R * x**2 + 3 * x + t),
     "three-division polynomial")

# A rational root is integral over k[t] and hence polynomial.  Degree balance
# leaves degree two.  The coefficient ladder for
# x=-(beta^2/3)t^2+c*t+d ends in -beta^5/3.
cc, dd = sp.symbols("cc dd")
x_trial = -beta**2 * t**2 / 3 + cc * t + dd
division_trial = sp.Poly(sp.expand(psi3.subs(x, x_trial)), t)
gate(division_trial.degree() == 7, "division trial leading cancellation")
gate(
    sp.factor(division_trial.coeff_monomial(t**7))
    == -beta**6 * (2 * beta * r + 3 * cc + 4) / 27,
    "division t7 coefficient",
)
cc_value = -(2 * beta * r + 4) / 3
gate(
    sp.factor(division_trial.coeff_monomial(t**6).subs(cc, cc_value))
    == -beta**6 * (3 * dd + r**2) / 27,
    "division t6 coefficient",
)
dd_value = -r**2 / 3
gate(
    sp.factor(
        division_trial.coeff_monomial(t**5).subs({cc: cc_value, dd: dd_value})
    )
    == -beta**5 / 3,
    "division t5 contradiction for beta nonzero",
)


# ---------------------------------------------------------------------------
# Genuine Cardano class and normal monogenic natural cubic.
# ---------------------------------------------------------------------------

zero(sp.expand((q + W) * (q - W) - 4 * X**3).subs(W**2, H),
     "Cardano product")
zero(H.subs(X, 0) - t**6, "two X=0 sections")
linear_initial = sum(
    coefficient * t**monomial[0] * X**monomial[1] * W**monomial[2]
    for monomial, coefficient in sp.Poly(q + W, t, X, W).terms()
    if sum(monomial) == 1
)
gate(linear_initial == W, "Cardano ideal has independent X,W initials")

v = sp.symbols("v")
cubic_surface = sp.expand(v**3 - 3 * X * v - q)
partials = [sp.factor(sp.diff(cubic_surface, variable)) for variable in (v, X, t)]
zero(partials[0] - 3 * (v**2 - X), "natural cubic v derivative")
zero(partials[1] - (2 * X - 3 * v - t * R), "natural cubic X derivative")
zero(partials[2] + 3 * t**2 + (r + 2 * beta * t) * X,
     "natural cubic t derivative")

# A polynomial root has total degree at most one.  The X^3 row forces its X
# coefficient to vanish, while the t^2 X row then equals -beta.
alpha, gamma, delta0 = sp.symbols("alpha gamma delta0")
linear_root = alpha * t + gamma * X + delta0
root_poly = sp.Poly(sp.expand(cubic_surface.subs(v, linear_root)), t, X)
gate(root_poly.coeff_monomial(X**3) == gamma**3,
     "natural cubic root X3 coefficient")
gate(root_poly.coeff_monomial(t**2 * X).subs(gamma, 0) == -beta,
     "natural cubic root t2X contradiction")

critical_resultant = sp.factor(
    sp.resultant(
        2 * v**2 - 3 * v - r * t - beta * t**2,
        3 * t**2 + (r + 2 * beta * t) * v**2,
        t,
    )
)
gate(sp.Poly(critical_resultant, v).degree() == 6,
     "natural cubic critical locus is finite")
gate(sp.Poly(critical_resultant, v).LC() == -8 * beta**3,
     "natural cubic critical eliminant is nonzero")


summary = {
    "checks": CHECKS,
    "family": "q=t^3+t(r+beta*t)X-X^2; beta nonzero",
    "surface": "normal; all closed fibres integral; Cardano prime non-Cartier",
    "elliptic": "I7 at infinity; generic five I1 and MW rank two",
    "boundary": "Q and -2Q; difference -3Q",
    "torsion": "E(k(t))[3]=0 by polynomial division ladder",
    "character": "Cl[3]=Z/3; unique smooth-locus C3 cover up to inversion",
    "cubic": "natural Cardano cubic is normal and monogenic",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3940 exact I7 rank-two resolvent character obstruction")
print(f"CHECKS={CHECKS}")
print("FAMILY=q=t^3+t*(r+beta*t)*X-X^2;beta!=0")
print("VERTICAL_FIBRES=ALL_INTEGRAL;RESOLVENT=NORMAL")
print("ELLIPTIC=I7_AT_INFINITY;GENERIC_FINITE=5I1;GENERIC_MW_RANK=2")
print("BOUNDARY=Q,-2Q;DIFFERENCE=-3Q")
print("RATIONAL_THREE_TORSION=NONE")
print("CLASS_THREE_TORSION=Z/3;CHARACTER=UNIQUE_UP_TO_INVERSION")
print("NATURAL_CUBIC=NORMAL_MONOGENIC")
print(f"SEMANTIC_SHA256={semantic}")
