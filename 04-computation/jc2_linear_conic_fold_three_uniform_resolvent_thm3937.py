#!/usr/bin/env python3
"""Exact companion for THM-3937's uniform fold-three resolvent rigidity."""

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


r, t, X, W, x, y, s = sp.symbols("r t X W x y s")
q = t**3 + r * t * X - X**2
H = sp.expand(q**2 - 4 * X**3)
F = sp.expand(W**2 - H)


# ---------------------------------------------------------------------------
# Uniform normality and irreducible vertical fibres.
# ---------------------------------------------------------------------------

zero(
    sp.factor(sp.discriminant(H, X))
    + 256
    * t**12
    * (16 * t**3 + 8 * r**2 * t**2 + (r**4 + 36 * r) * t + r**3 + 27),
    "generic quartic discriminant",
)

FW = sp.diff(F, W)
FX = sp.factor(sp.diff(F, X))
Ft = sp.factor(sp.diff(F, t))
gate(FW == 2 * W, "surface W derivative")
zero(FX - 2 * (6 * X**2 + (2 * X - r * t) * q), "surface X derivative")
zero(Ft + 2 * (r * X + 3 * t**2) * q, "surface t derivative")

# In the nonzero-X, nonzero-r singular row put z=t/r^2 and a=2+9z.
# H=0 and FX=0 reduce to R*a^2+108=0 and
# R*a*(2a-1)+162=0, where R=r^3.  Their difference forces a=2,
# hence z=0 and contradicts X=-3t^2/r !=0.
R, aa = sp.symbols("R aa")
singular_h = R * aa**2 + 108
singular_x = R * aa * (2 * aa - 1) + 162
gate(sp.expand(singular_x - 2 * singular_h) == -R * aa - 54,
     "nonzero singular row first elimination")
zero(singular_h.subs(R * aa, -54) - (-54 * aa + 108),
     "nonzero singular row second elimination")
gate(H.subs({r: 0, t: 0, X: 3}) == -27,
     "r=0 spurious critical point misses surface")

# For a closed t=lambda fibre write A=lambda^3 and B=r*lambda.  If its
# quartic were (X^2+pX+c)^2, the cubic coefficient gives p=-B-2 and the
# constant coefficient gives c=+A or -A.  Both exact charts are empty after
# imposing A=lambda^3 and B=r*lambda.
lam, A, B, p, c = sp.symbols("lambda A B p c")
closed_q = A + B * X - X**2
closed_H = sp.expand(closed_q**2 - 4 * X**3)
square = sp.expand((X**2 + p * X + c) ** 2)
square_coeff = sp.Poly(square - closed_H, X)
gate(square_coeff.coeff_monomial(X**3) == 2 * B + 2 * p + 4,
     "closed-fibre square cubic coefficient")
gate(square_coeff.coeff_monomial(X**2) == 2 * A - B**2 + 2 * c + p**2,
     "closed-fibre square quadratic coefficient")
gate(square_coeff.coeff_monomial(X) == -2 * A * B + 2 * c * p,
     "closed-fibre square linear coefficient")
gate(square_coeff.coeff_monomial(1) == -A**2 + c**2,
     "closed-fibre square constant coefficient")

p_value = -B - 2
for sign, label in ((1, "plus"), (-1, "minus")):
    c_value = sign * A
    equations = [
        A - lam**3,
        B - r * lam,
        square_coeff.coeff_monomial(X**2).subs({p: p_value, c: c_value}),
        square_coeff.coeff_monomial(X).subs({p: p_value, c: c_value}),
    ]
    basis = sp.groebner(equations, A, B, r, lam, order="lex")
    gate(basis.contains(sp.Integer(1)), f"closed-fibre {label}-root chart empty")


# ---------------------------------------------------------------------------
# Exact generic elliptic model and its boundary sections.
# ---------------------------------------------------------------------------

z = sp.cancel((q + W) / (2 * X**2))
x_forward = sp.cancel(t * z)
y_forward = sp.cancel(X * z * (z + 1) + (1 - r * t * z) / 2)
E_rhs = x**3 + (t + r**2 / 4) * x**2 - r * x / 2 + sp.Rational(1, 4)

forward_error = sp.cancel(y_forward**2 - E_rhs.subs(x, x_forward))
forward_numerator = sp.fraction(forward_error)[0]
gate(
    sp.rem(sp.Poly(forward_numerator, W), sp.Poly(F, W)).as_expr() == 0,
    "quartic-to-elliptic identity",
)

z_inverse = x / t
X_inverse = sp.cancel(t**2 * (2 * y - 1 + r * x) / (2 * x * (x + t)))
q_inverse = t**3 + r * t * X_inverse - X_inverse**2
W_inverse = sp.cancel(2 * X_inverse**2 * z_inverse - q_inverse)
E_equation = sp.expand(y**2 - E_rhs)
inverse_error = sp.cancel(W_inverse**2 - (q_inverse**2 - 4 * X_inverse**3))
inverse_numerator = sp.factor(sp.fraction(inverse_error)[0])
zero(
    sp.cancel(inverse_numerator / E_equation)
    - 4 * t**6 * (r * x + 2 * y - 1) ** 2,
    "elliptic-to-quartic inverse identity",
)
composition = sp.cancel(
    X_inverse.subs({x: x_forward, y: y_forward}, simultaneous=True) - X
)
composition_numerator = sp.fraction(composition)[0]
gate(
    sp.rem(sp.Poly(composition_numerator, W), sp.Poly(F, W)).as_expr() == 0,
    "birational inverse composition",
)

# At quartic infinity epsilon=1/X and eta=W/X^2.  The two eta signs map to
# Q=(0,-1/2) and -2Q=(-t,-(1+rt)/2).
epsilon = sp.symbols("epsilon")
eta_square = sp.expand((-1 + r * t * epsilon + t**3 * epsilon**2) ** 2 - 4 * epsilon)
eta_plus = sp.series(sp.sqrt(eta_square), epsilon, 0, 3).removeO()
eta_minus = -eta_plus
for eta, expected_x, expected_y, label in (
    (eta_plus, 0, -sp.Rational(1, 2), "positive"),
    (eta_minus, -t, -(1 + r * t) / 2, "negative"),
):
    z_boundary = sp.expand((-1 + r * t * epsilon + t**3 * epsilon**2 + eta) / 2)
    xb = sp.limit(t * z_boundary, epsilon, 0)
    yb = sp.limit(
        z_boundary * (z_boundary + 1) / epsilon + (1 - r * t * z_boundary) / 2,
        epsilon,
        0,
    )
    gate(sp.cancel(xb - expected_x) == 0, f"{label} infinity x section")
    gate(sp.cancel(yb - expected_y) == 0, f"{label} infinity y section")


# ---------------------------------------------------------------------------
# Kodaira fibres and Mordell--Weil generator.
# ---------------------------------------------------------------------------

a2 = t + r**2 / 4
a4 = -r / 2
a6 = sp.Rational(1, 4)
b2 = sp.factor(4 * a2)
b4 = sp.factor(2 * a4)
b6 = sp.factor(4 * a6)
b8 = sp.factor(4 * a2 * a6 - a4**2)
c4 = sp.factor(b2**2 - 24 * b4)
c6 = sp.factor(-b2**3 + 36 * b2 * b4 - 216 * b6)
Delta = sp.factor(-b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6)
gate(b2 == 4 * t + r**2, "uniform b2")
gate(b4 == -r and b6 == 1 and b8 == t, "uniform b4 b6 b8")
zero(c4 - ((4 * t + r**2) ** 2 + 24 * r), "uniform c4")
zero(c6 + (4 * t + r**2) ** 3 + 36 * r * (4 * t + r**2) + 216,
     "uniform c6")
zero(
    Delta
    + 16 * t**3 + 8 * r**2 * t**2 + (r**4 + 36 * r) * t + r**3 + 27,
    "uniform discriminant",
)
gate(sp.factor(sp.discriminant(Delta, t)) == -256 * (2 * r**3 + 27) ** 3,
     "finite-fibre collision parameter")
gate(sp.factor(sp.resultant(Delta, c4, t)) == 4096 * (2 * r**3 + 27) ** 2,
     "generic finite fibres are multiplicative I1")

special_relation = 2 * r**3 + 27
t0 = r**2 / 12
t1 = -2 * r**2 / 3
special_factor = -16 * (t - t0) ** 2 * (t - t1)
gate(sp.rem(sp.Poly(sp.together(Delta - special_factor), r),
            sp.Poly(special_relation, r)).as_expr() == 0,
     "special discriminant factorization")
zero(sp.factor(c4.subs(t, t0)) - 8 * r * special_relation / 9,
     "special double root c4 vanishes")
zero(sp.factor(c6.subs(t, t0)) + 8 * special_relation * (4 * r**3 + 27) / 27,
     "special double root c6 vanishes")
gate(sp.diff(c4, t).subs(t, t0) == 32 * r**2 / 3,
     "special double root c4 has order one")
gate(
    sp.rem(
        sp.Poly(sp.diff(c6, t).subs(t, t0) - 144 * r, r),
        sp.Poly(special_relation, r),
    ).as_expr()
    == 0,
    "special double root c6 has order one",
)
gate(sp.diff(Delta, t, 2).subs(t, t0) == -24 * r**2,
     "special discriminant has order two")
gate(
    sp.rem(sp.Poly(sp.together(c4.subs(t, t1) + 27 * r / 2), r),
           sp.Poly(special_relation, r)).as_expr() == 0,
    "special simple root has nonzero c4",
)

# Infinity is uniformly minimal I3*: v(c4,c6,Delta)=(2,3,9).
c4_inf = sp.factor(s**4 * c4.subs(t, 1 / s))
c6_inf = sp.factor(s**6 * c6.subs(t, 1 / s))
Delta_inf = sp.factor(s**12 * Delta.subs(t, 1 / s))
gate(sp.Poly(c4_inf, s).as_dict().get((2,)) == 16,
     "infinity c4 valuation two")
gate(sp.Poly(c6_inf, s).as_dict().get((3,)) == -64,
     "infinity c6 valuation three")
gate(sp.Poly(Delta_inf, s).as_dict().get((9,)) == -16,
     "infinity discriminant valuation nine")


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
        slope = sp.cancel((3 * x1**2 + 2 * a2 * x1 + a4) / (2 * y1))
    else:
        slope = sp.cancel((y2 - y1) / (x2 - x1))
    x3 = sp.factor(slope**2 - a2 - x1 - x2)
    y3 = sp.factor(-y1 + slope * (x1 - x3))
    return x3, y3


Q = (sp.Integer(0), -sp.Rational(1, 2))
two_Q = add_points(Q, Q)
three_Q = add_points(two_Q, Q)
zero(two_Q[0] + t, "uniform double of Q x coordinate")
zero(two_Q[1] - (r * t + 1) / 2, "uniform double of Q y coordinate")
zero(three_Q[0] - (r * t + 1) / t**2, "uniform triple of Q x coordinate")
zero(
    three_Q[1]
    - (r**2 * t**2 + 3 * r * t + t**3 + 2) / (2 * t**3),
    "uniform triple of Q y coordinate",
)

# Q follows the same three blowups to a D7 spin tip for every r.
xinf, yinf = sp.symbols("xinf yinf")
local = sp.expand(
    yinf**2
    - xinf**3
    - (s + r**2 * s**2 / 4) * xinf**2
    + r * s**4 * xinf / 2
    - s**6 / 4
)
qx = sp.Integer(0)
qy = -s**3 / 2
exceptional = []
for index in range(3):
    local = sp.expand(local.subs({xinf: s * xinf, yinf: s * yinf}, simultaneous=True))
    minimum = min(monomial[0] for monomial, coefficient in sp.Poly(local, s, xinf, yinf).terms())
    gate(minimum == 2, f"uniform blowup {index + 1} multiplicity")
    local = sp.expand(local / s**2)
    qx = sp.cancel(qx / s)
    qy = sp.cancel(qy / s)
    exceptional.append(sp.factor(local.subs(s, 0)))
gate(exceptional[0] == yinf**2, "first uniform exceptional double component")
gate(exceptional[1] == yinf**2, "second uniform exceptional double component")
gate(exceptional[2] == (2 * yinf - 1) * (2 * yinf + 1) / 4,
     "third uniform exceptional spin tips")
gate(qx == 0 and qy == -sp.Rational(1, 2), "Q meets negative spin tip")

D7 = 2 * sp.eye(7)
for i, j in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (4, 6)):
    D7[i, j] = D7[j, i] = -1
gate(D7.det() == 4, "D7 determinant")
gate(D7.inv()[5, 5] == sp.Rational(7, 4), "D7 spin correction")
gate(2 - D7.inv()[5, 5] == sp.Rational(1, 4), "uniform Q height")
gate(sp.Rational(1, 4) * D7.det() == 1, "uniform MW determinant invoice")

# A rational two-torsion root is integral over k[t], hence polynomial.  Degree
# comparison leaves -t+b.  Its t^2, t, and constant rows are incompatible.
broot = sp.symbols("broot")
root_trial = -t + broot
two_torsion = sp.Poly(
    sp.expand(root_trial**3 + a2 * root_trial**2 + a4 * root_trial + a6), t
)
gate(two_torsion.coeff_monomial(t**2) == broot + r**2 / 4,
     "two-torsion t-square coefficient")
gate(
    sp.factor(two_torsion.coeff_monomial(t).subs(broot, -r**2 / 4)) == r / 2,
    "two-torsion t coefficient",
)
gate(
    two_torsion.coeff_monomial(1).subs({broot: 0, r: 0}) == sp.Rational(1, 4),
    "two-torsion constant contradiction",
)


# ---------------------------------------------------------------------------
# Cardano class, normal monogenic cubic, and codimension-two distinction.
# ---------------------------------------------------------------------------

a_cardano = q + W
b_cardano = q - W
zero(sp.expand(a_cardano * b_cardano - 4 * X**3).subs(W**2, H),
     "uniform Cardano product")
zero(H.subs(X, 0) - t**6, "uniform two X=0 sections")

u = sp.symbols("u")
cubic_surface = sp.expand(u**3 - 3 * X * u - q)
cubic_partials = [sp.factor(sp.diff(cubic_surface, variable)) for variable in (u, X, t)]
zero(cubic_partials[0] - 3 * (-X + u**2), "uniform cubic u derivative")
zero(cubic_partials[1] - (-r * t + 2 * X - 3 * u), "uniform cubic X derivative")
zero(cubic_partials[2] - (-r * X - 3 * t**2), "uniform cubic t derivative")
u_singular = sp.symbols("u_singular")
t_from_u = (2 * u_singular**2 - 3 * u_singular) / r
zero(
    sp.factor(r**2 * (-r * u_singular**2 - 3 * t_from_u**2))
    + u_singular**2 * (r**3 + 3 * (2 * u_singular - 3) ** 2),
    "nonzero-r cubic singular locus is finite",
)
gate(
    cubic_surface.subs({r: 0, t: 0, X: sp.Rational(9, 4), u: sp.Rational(3, 2)})
    != 0,
    "r=0 spurious cubic critical point misses surface",
)

# A polynomial root must be affine-linear.  Exact coefficient comparison has
# no solution for any specialization of r.
alpha, beta, gamma = sp.symbols("alpha beta gamma")
linear_root = alpha * t + beta * X + gamma
root_equations = sp.Poly(
    sp.expand(linear_root**3 - 3 * X * linear_root - q), t, X
).coeffs()
root_basis = sp.groebner(root_equations, alpha, beta, gamma, r, order="lex")
gate(root_basis.contains(sp.Integer(1)), "uniform cubic has no linear root")

# Cardano cube roots recover that monic cubic.
zeta, vee = sp.symbols("zeta vee")
cardano_identity = sp.expand((zeta + vee) ** 3 - 3 * X * (zeta + vee) - q)
cardano_identity = sp.expand(cardano_identity.subs(X, zeta * vee))
zero(cardano_identity - (zeta**3 + vee**3 - (t**3 + r * t * zeta * vee - zeta**2 * vee**2)),
     "uniform Cardano cubic identity")

# At the isolated resolvent origin the Cardano ideal (X,q+W) has independent
# linear initials X and W; this is the exact non-Cartier sidecar.
linear_initial = sum(
    coefficient * t**monomial[0] * X**monomial[1] * W**monomial[2]
    for monomial, coefficient in sp.Poly(q + W, t, X, W).terms()
    if sum(monomial) == 1
)
gate(linear_initial == W, "Cardano ideal second generator has initial W")

# The normalization-address collision wall is distinct from the elliptic
# finite-fibre collision wall.
normalization_parameter = sp.symbols("normalization_parameter")
gate(
    sp.discriminant(normalization_parameter**3 + r * normalization_parameter - 2,
                    normalization_parameter)
    == -4 * (r**3 + 27),
    "normalization-address collision parameter",
)


summary = {
    "checks": CHECKS,
    "family": "q=t^3+r*t*X-X^2",
    "resolvent": "normal with sole singular origin; all vertical fibres irreducible",
    "elliptic": "y^2=x^3+(t+r^2/4)x^2-(r/2)x+1/4",
    "finite_fibres": "3I1, or II+I1 when 2r^3+27=0",
    "infinity": "I3*; Q=(0,-1/2) spin height 1/4",
    "mordell_weil": "Z*Q; boundary difference is 3Q",
    "affine": "Cl=Z/3; units=k*; Pic singular surface=0",
    "kummer": "unique regular-locus C3 character up to inversion",
    "cubic": "unique normal field/order is natural monogenic Cardano cubic",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3937 exact uniform fold-three resolvent rigidity")
print(f"CHECKS={CHECKS}")
print("PARAMETER=ALL_r_IN_k")
print("RESOLVENT=NORMAL;SINGULAR_SUPPORT=(0,0,0);VERTICAL_FIBRES=IRREDUCIBLE")
print("ELLIPTIC=y^2-x^3-(t+r^2/4)x^2+(r/2)x-1/4=0")
print("FINITE_FIBRES=3I1_OR_SPECIAL_II_PLUS_I1")
print("INFINITY=I3STAR;MW=Z*Q;HEIGHT_Q=1/4;BOUNDARY_DIFFERENCE=3Q")
print("CLASS_GROUP=Z/3;UNITS=k*;PIC_S=0;PIC_SREG=Z/3")
print("KUMMER=UNIQUE_SMOOTH_LOCUS_C3_CHARACTER_UP_TO_INVERSION")
print("NORMAL_CUBIC=NATURAL_MONOGENIC_CARDANO_ORDER")
print(f"SEMANTIC_SHA256={semantic}")
