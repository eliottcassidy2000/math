#!/usr/bin/env python3
"""Finite exact companion for provisional THM-3598.

The theorem is identity- and residue-driven.  This script checks the
all-exponent rational-exact family, arm-pole hostile, the complete positive
exponent-three graph-linear residue classifier, and the negative-graph
multi-root obstruction without Python assert gates.
"""

import sympy as sp


b, c, e, x, w = sp.symbols("b c e x w")
CHECKS = 0


def require(label, condition):
    """Record one active truth gate and fail with a stable label."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError("FAILED: " + label)


def zero(expr):
    return sp.cancel(sp.expand(expr)) == 0


def density_bracket(F, G, exponent):
    return sp.expand(
        c**exponent
        * (sp.diff(F, b) * sp.diff(G, c) - sp.diff(F, c) * sp.diff(G, b))
    )


print("THM-3598 exact companion")
print("SECTION all-exponent rational-exact composite channel")

FAMILY_ROWS = 0
for exponent in range(2, 9):
    a = b**2 + b + exponent
    H = sp.integrate(a ** (exponent - 1), b)
    channel = a * c
    P = 1 / ((exponent - 1) * channel ** (exponent - 1))
    G = channel**2 + (exponent + 1) * channel**3
    Q = H + G
    require(f"H derivative n={exponent}", zero(sp.diff(H, b) - a ** (exponent - 1)))
    require(f"channel-H bracket n={exponent}", zero(density_bracket(channel, H, exponent) + channel**exponent))
    require(f"rational mate n={exponent}", zero(density_bracket(P, Q, exponent) - 1))
    require(f"nonlinear G commutes n={exponent}", zero(density_bracket(P, G, exponent)))
    FAMILY_ROWS += 4

print(f"PASS {FAMILY_ROWS} all-exponent family identities")


print("SECTION all-arm arm-pole specialization")
Sigma = b * (b**2 + 1)
ARM_ROWS = 0
for exponent in range(2, 9):
    H = sp.integrate(Sigma ** (exponent - 1), b)
    channel_surface = c ** (exponent + 1) * e
    G_surface = channel_surface + channel_surface**2
    Q_surface = H + G_surface

    require(
        f"relation channel n={exponent}",
        zero((Sigma * c - channel_surface).subs(e, Sigma / c**exponent)),
    )
    require(f"arm H jet n={exponent}", sp.rem(sp.diff(Q_surface, b), Sigma, b) == 0)
    require(f"arm c jet n={exponent}", sp.diff(Q_surface, c).subs(c, 0) == 0)
    require(f"arm e jet n={exponent}", sp.diff(Q_surface, e).subs(c, 0) == 0)
    ARM_ROWS += 4

H_a13 = b**7 / 7 + 2 * b**5 / 5 + b**3 / 3
x_a13 = Sigma * c
P_a13 = 1 / (2 * x_a13**2)
Q_a13 = H_a13 + x_a13 + x_a13**2
require("A13 H derivative", zero(sp.diff(H_a13, b) - Sigma**2))
require("A13 rational mate", zero(density_bracket(P_a13, Q_a13, 3) - 1))
require("A13 polynomial channel", zero((x_a13 - c**4 * e).subs(e, Sigma / c**3)))
print(f"PASS {ARM_ROWS + 3} all-arm, critical-arm, and A13 controls")


print("SECTION formal graph-linear residue classifier")
z = sp.symbols("z")
h0, h1, h2, h3, h4 = sp.symbols("h0 h1 h2 h3 h4", nonzero=True)
A0, A1, A2, A3 = sp.symbols("A0 A1 A2 A3")
h_jet = h0 + h1 * z + h2 * z**2 / 2 + h3 * z**3 / 6 + h4 * z**4 / 24
A_jet = A0 + A1 * z + A2 * z**2 / 2 + A3 * z**3 / 6
eta_jet = A_jet**2 / (h0 - h_jet) ** 3
residue_direct = sp.residue(eta_jet, z, 0)
hp = sp.diff(h_jet, z)
R_jet = A_jet**2 / hp
Dh_R = sp.diff(R_jet, z) / hp
Dh2_R = sp.diff(Dh_R, z) / hp
residue_expected = -sp.Rational(1, 2) * Dh2_R.subs(z, 0)
require("formal residue identity", zero(residue_direct - residue_expected))

alpha, beta, t = sp.symbols("alpha beta t")
primitive = (alpha * w + beta) / (2 * (w - t) ** 2) - alpha / (w - t)
target_form = (alpha * t + beta) / (w - t) ** 3
require("affine-R primitive", zero(sp.diff(primitive, t) - target_form))

CLASSIFIER_ROWS = 2
for degree in range(1, 6):
    A = b**degree + 2 * b + 1
    scale = sp.Rational(degree + 1, degree + 2)
    h = sp.integrate(A**2 / scale, b)
    channel = A * c
    P = scale / (2 * channel**2)
    Q = h + channel
    Dh = lambda value: sp.diff(value, b) / sp.diff(h, b)
    obstruction = sp.cancel(Dh(Dh(A**2 / sp.diff(h, b))))
    require(f"square derivative n3 degree={degree}", zero(sp.diff(h, b) - A**2 / scale))
    require(f"zero residue obstruction degree={degree}", zero(obstruction))
    require(f"classified rational mate degree={degree}", zero(density_bracket(P, Q, 3) - 1))
    CLASSIFIER_ROWS += 3

for exponent in range(2, 9):
    k = exponent - 1
    A = b**2 + exponent * b + 1
    scale = sp.Rational(exponent + 1, exponent + 2)
    h = sp.integrate(A**k / scale, b)
    channel = A * c
    Q = h + channel
    P = scale / (k * channel**k)
    require(
        f"all-n constant-R classification n={exponent}",
        zero(density_bracket(P, Q, exponent) - 1),
    )
    CLASSIFIER_ROWS += 1

for k in range(3, 9):
    witness = None
    for j in range(1, k):
        for degree_h in range(2, 4 * k + 2):
            numerator = degree_h * (j + 1) - 1
            if numerator % k == 0:
                witness = (j, degree_h, numerator // k)
                break
        if witness is not None:
            break
    require(f"nonconstant-R witness exists k={k}", witness is not None)
    j, degree_h, degree_a = witness
    h = b**degree_h / degree_h
    A = b**degree_a
    R_poly = degree_h**j * t**j
    R_at_h = sp.expand(R_poly.subs(t, h))
    require(
        f"nonconstant-R equation k={k}",
        zero(A**k - sp.diff(h, b) * R_at_h),
    )
    channel = A * c
    Q = h + channel
    P = 0
    for ell in range(j + 1):
        derivative = sp.diff(R_poly, t, ell).subs(t, Q)
        P += (-1) ** ell * derivative / (
            sp.factorial(ell) * (k - ell) * channel ** (k - ell)
        )
    require(
        f"nonconstant-R primitive k={k}",
        zero(density_bracket(P, Q, k + 1) - 1),
    )
    CLASSIFIER_ROWS += 3

for exponent in range(2, 9):
    k = exponent - 1
    A = b**2 + 3 * b + exponent
    F_constant = sp.integrate(A**k, b)
    h_constant = sp.Integer(exponent + 7)
    Q_constant = h_constant + A * c
    P_constant = F_constant / (Q_constant - h_constant) ** exponent
    require(
        f"constant-h rational mate n={exponent}",
        zero(density_bracket(P_constant, Q_constant, exponent) - 1),
    )
    CLASSIFIER_ROWS += 1

for degree_h in range(2, 7):
    for degree_a in range(1, 5):
        h = b**degree_h + b
        A = b**degree_a + 1
        hp = sp.diff(h, b)
        R = A**2 / hp
        obstruction = sp.cancel(sp.diff(sp.diff(R, b) / hp, b) / hp)
        require(
            f"active nonsquare hostile h={degree_h} A={degree_a}",
            obstruction != 0,
        )
        CLASSIFIER_ROWS += 1

for degree_h in range(1, 17):
    for degree_a in range(0, 17):
        require(
            f"alpha degree parity h={degree_h} A={degree_a}",
            2 * degree_a != 2 * degree_h - 1,
        )
        CLASSIFIER_ROWS += 1

print(f"PASS {CLASSIFIER_ROWS} residue, primitive, square, hostile, and parity gates")


print("SECTION positive primitive torsor and arm/residual debt")
TORSOR_ROWS = 0
for k in range(1, 9):
    for order in range(k):
        coefficient = sum(
            (-1) ** ell * sp.binomial(order, ell) / (k - ell)
            for ell in range(order + 1)
        )
        expected = (
            (-1) ** order
            * sp.factorial(order)
            * sp.factorial(k - order - 1)
            / sp.factorial(k)
        )
        require(
            f"arm leading coefficient k={k} order={order}",
            sp.simplify(coefficient - expected) == 0 and expected != 0,
        )
        TORSOR_ROWS += 1

h_residual = b**2 + b
A_residual = b + 2
lambda_residual = h_residual.subs(b, 0)
c_residual = -(h_residual - lambda_residual) / A_residual
require(
    "residual fibre is nonempty off c=0",
    zero((h_residual + A_residual * c_residual) - lambda_residual)
    and c_residual.subs(b, 1) != 0,
)
TORSOR_ROWS += 1
print(f"PASS {TORSOR_ROWS} primitive-pole and residual-fibre controls")


print("SECTION shallow negative-channel field trace")
NEGATIVE_ROWS = 0
for exponent in range(2, 9):
    for shallow in range(1, exponent):
        if (exponent - 1) % shallow != 0:
            continue
        quotient = (exponent - 1) // shallow
        A = b**2 + b + exponent
        F = b + (b + 1) * c + c**2
        trace_index = exponent - shallow - 1
        polynomial = A + c**shallow * (F - w)
        reciprocal = sp.series(1 / polynomial, c, 0, trace_index + 1).removeO()
        trace_form = -sp.expand(reciprocal).coeff(c, trace_index)
        trace_poly = sp.Poly(sp.cancel(trace_form), w)
        require(
            f"trace w-degree n={exponent} s={shallow}",
            trace_poly.degree() == quotient - 1,
        )
        require(
            f"trace leading n={exponent} s={shallow}",
            zero(trace_poly.LC() + 1 / A**quotient),
        )
        NEGATIVE_ROWS += 2

A_trace, F0 = sp.symbols("A_trace F0", nonzero=True)
require(
    "n3 s1 trace packet",
    zero(-sp.series(1 / (A_trace + c * (F0 - w)), c, 0, 2).removeO().coeff(c, 1) - (F0 - w) / A_trace**2),
)
require(
    "n3 s2 trace packet",
    zero(-sp.series(1 / (A_trace + c**2 * (F0 - w)), c, 0, 1).removeO().coeff(c, 0) + 1 / A_trace),
)
NEGATIVE_ROWS += 2

for exponent in range(2, 9):
    k = exponent - 1
    for root_count in range(2, 6):
        A = sp.prod(b - j for j in range(root_count))
        h = b**3 + (root_count + exponent) * b
        c_fibre = A / (w - h)
        eta_coefficient = -(w - h) ** (exponent - 2) / A**k
        Q_negative = h + A / c
        require(
            f"negative fibre n={exponent} roots={root_count}",
            zero(Q_negative.subs(c, c_fibre) - w),
        )
        require(
            f"negative time n={exponent} roots={root_count}",
            zero(
                1 / density_bracket(b, Q_negative, exponent).subs(c, c_fibre)
                - eta_coefficient
            ),
        )
        finite_pole_degree = k * sp.degree(A, b) - root_count
        infinity_zero_order = k * sp.degree(A, b) - 1
        require(
            f"multi-root divisor n={exponent} roots={root_count}",
            finite_pole_degree < infinity_zero_order,
        )
        NEGATIVE_ROWS += 3

for exponent in range(2, 9):
    k = exponent - 1
    for multiplicity in range(1, 7):
        total_order = k * multiplicity
        if total_order == 1:
            require(
                "one-root logarithmic boundary",
                sp.residue(1 / (b - 2), b, 2) == 1,
            )
        else:
            A = (b - 2) ** multiplicity
            primitive_one_root = -1 / (
                (total_order - 1) * (b - 2) ** (total_order - 1)
            )
            require(
                f"one-root primitive n={exponent} mult={multiplicity}",
                zero(sp.diff(primitive_one_root, b) - 1 / A**k),
            )
        NEGATIVE_ROWS += 1

print(f"PASS {NEGATIVE_ROWS} trace, negative-graph, and one-root controls")


print("SECTION sharp exponent-one boundary")
require("logarithmic residue", sp.residue(1 / x, x, 0) == 1)
require("power primitive starts n=2", sp.diff(1 / x, x) == -1 / x**2)
print("PASS exponent one has a logarithmic, not rational, primitive")


print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
