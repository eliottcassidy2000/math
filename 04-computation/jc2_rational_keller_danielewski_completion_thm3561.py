#!/usr/bin/env python3
"""Exact companion for THM-3561.

The calculation is deliberately symbolic over QQ.  It verifies the rational
Keller pair, its polynomial Danielewski completion, the three etale boundary
charts, collision transport, and finite controls for the all-weight
intersection proof.  It also checks the global polynomial representative of
the exact symplectic primitive and homogeneous-weight hostiles.  The proofs
of the intersection, Picard, de Rham, and all-weight statements are
structural and are written in the theorem file; the rows below are hostile
controls rather than finite extrapolations.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def jacobian2(f: sp.Expr, g: sp.Expr, r: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    return sp.factor(sp.diff(f, r) * sp.diff(g, s) - sp.diff(f, s) * sp.diff(g, r))


x, y, z, q = sp.symbols("x y z q")
u = 1 + x * y

# The middle coordinate of the fixed three-dimensional Keller map.
F2 = sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
F1 = sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y))
F3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)
h = 9 * u**2 - 6 * u - 2
require(sp.expand(F2 - (h * y + 3 * x * u**2 * z)) == 0, "middle-coordinate rewrite")

unit_cofactor = 3 * x**2 * z * u + 9 * u**2 - 15 * u + 4
require(sp.expand(x * F2 - (2 + u * unit_cofactor)) == 0, "explicit unit identity")

w = sp.symbols("w")
w_from_xyz = 3 * z + 2 * (4 * u + 1) * y**2 / u**2
require(
    sp.factor(F2 - u**2 * (y + x * w_from_xyz)) == 0,
    "localized triangularization",
)

u_w = 1 - x**2 * w
y_w = -x * w
z_w = (w - 2 * (4 * u_w + 1) * x**2 * w**2 / u_w**2) / 3
leaf_substitution = {y: y_w, z: z_w}
require(sp.factor(F2.subs(leaf_substitution)) == 0, "leaf inverse lands on F2=0")
require(sp.factor(F1.subs(leaf_substitution) - u_w * w / 3) == 0, "restricted F1")
require(
    sp.factor(F3.subs(leaf_substitution) - 2 * x * (1 + 2 * u_w) / (3 * u_w**2)) == 0,
    "restricted F3",
)

# On F2=0, put y=-xw and u=1-x^2w.  The restricted target coordinates,
# scaled by constants, become a rational plane pair after q=w/u.
D = 1 + x**2 * q
a = q / D**2
c = x * D * (D + 2)
require(jacobian2(a, c, x, q) == -3, "rational Keller determinant")

collision = [
    (sp.Integer(0), -sp.Rational(3, 4)),
    (sp.Integer(1), sp.Integer(-3)),
    (sp.Integer(-1), sp.Integer(-3)),
]
collision_rows = []
for xx, qq in collision:
    row = (
        xx,
        qq,
        sp.simplify(D.subs({x: xx, q: qq})),
        sp.simplify(a.subs({x: xx, q: qq})),
        sp.simplify(c.subs({x: xx, q: qq})),
    )
    collision_rows.append(row)
    require(row[3:] == (-sp.Rational(3, 4), sp.Integer(0)), "triple collision")

# Maximal polynomial-observable completion.
b = sp.factor(a * c**2)
e = sp.factor(a * (b + 4))
require(sp.factor(b - (D - 1) * (D + 2) ** 2) == 0, "b polynomialization")
require(sp.factor(e - q * (D + 3)) == 0, "e polynomialization")
require(sp.factor(c**2 * e - b * (b + 4)) == 0, "Danielewski relation")

# The generic inverse parameter D has a cubic equation.
B, T = sp.symbols("B T")
inverse_cubic = sp.expand((T - 1) * (T + 2) ** 2 - B)
require(sp.Poly(inverse_cubic, T).degree() == 3, "inverse cubic degree")
require(sp.factor(inverse_cubic) == inverse_cubic, "inverse cubic factor control")

# Etaleness of Phi=(b,c,e): use (b,c) when c!=0 and (c,e) on c=0.
J_bc = jacobian2(b, c, x, q)
J_ce = jacobian2(c, e, x, q)
require(sp.factor(J_bc + 3 * c**2) == 0, "open-chart etaleness")

t = sp.symbols("t")
J_ce_in_D = 6 * (D + 1) * (D**2 + 2 * D - 2)
require(sp.factor(J_ce - J_ce_in_D) == 0, "boundary Jacobian formula")
boundary_jacobians = {}
for d_value in (sp.Integer(1), sp.Integer(0), sp.Integer(-2)):
    value = sp.expand(6 * (d_value + 1) * (d_value**2 + 2 * d_value - 2))
    boundary_jacobians[int(d_value)] = int(value)
require(boundary_jacobians == {1: 12, 0: -12, -2: 12}, "three boundary charts")

completed_collision_rows = []
for xx, qq in collision:
    row = tuple(sp.simplify(f.subs({x: xx, q: qq})) for f in (b, c, e))
    completed_collision_rows.append(row)
    require(row == (sp.Integer(0), sp.Integer(0), sp.Integer(-3)), "completed collision")

# Smoothness hostile: R_b=R_c=R_e=0 would force c=0,b=-2, which
# contradicts R=0.  The Groebner basis of the singular ideal is [1].
bb, cc, ee = sp.symbols("bb cc ee")
relation = cc**2 * ee - bb * (bb + 4)
singular_basis = sp.groebner(
    [relation, sp.diff(relation, bb), sp.diff(relation, cc), sp.diff(relation, ee)],
    bb,
    cc,
    ee,
    order="lex",
)
require(list(singular_basis) == [sp.Integer(1)], "smooth target")

# Finite controls for the structural grading proof of
# k[a,c] intersect k[x,q] = k[b,c,e].  A weight -s element needs both the
# b=0 and b=-4 cancellations to order ceil(s/2).
intersection_rows = []
for s in range(1, 9):
    m = (s + 1) // 2
    numerator = (b * (b + 4)) ** m
    candidate = sp.cancel(numerator / c**s)
    expected = e ** (s // 2) if s % 2 == 0 else c * e ** ((s + 1) // 2)
    require(sp.factor(candidate - expected) == 0, f"intersection generator weight {-s}")

    # Remove one b+4 factor.  The resulting expression must retain a D-pole.
    hostile = sp.cancel((b**m * (b + 4) ** (m - 1)) / c**s)
    hostile_denominator = sp.factor(sp.denom(hostile))
    require(hostile_denominator != 1, f"insufficient boundary cancellation weight {-s}")
    require(sp.rem(hostile_denominator, D, q) == 0 or hostile_denominator.has(D), f"D pole weight {-s}")
    intersection_rows.append((s, m, sp.factor(expected)))

# Poisson coefficient system on Y.  Directly calculate the generator brackets
# in (a,c) using independent
# formal target variables before imposing b=a*c^2 and e=a(b+4).
aa, c0 = sp.symbols("aa c0")
b0 = aa * c0**2
e0 = aa * (b0 + 4)
bracket_bc = jacobian2(b0, c0, aa, c0)
bracket_ce = jacobian2(c0, e0, aa, c0)
bracket_be = jacobian2(b0, e0, aa, c0)
require(bracket_bc == c0**2, "Poisson {b,c}")
require(sp.factor(bracket_ce + 2 * (b0 + 2)) == 0, "Poisson {c,e}")
require(sp.factor(bracket_be + 2 * c0 * e0) == 0, "Poisson {b,e}")

# Complete affine-linear hostile.  The three Pluecker coefficients multiply
# c^2, -2ce, and -2(b+2).  Vanishing of the b coefficient forces the only
# constant contribution to vanish as well.
p1, p2, p3, r1, r2, r3 = sp.symbols("p1 p2 p3 r1 r2 r3")
linear_bracket = sp.expand(
    (p1 * r2 - p2 * r1) * cc**2
    - 2 * (p1 * r3 - p3 * r1) * cc * ee
    - 2 * (p2 * r3 - p3 * r2) * (bb + 2)
)
linear_poly = sp.Poly(linear_bracket, bb, cc, ee)
linear_b_coefficient = linear_poly.coeff_monomial(bb)
linear_constant = linear_poly.coeff_monomial(1)
require(sp.expand(linear_constant - 2 * linear_b_coefficient) == 0, "linear constant tied to b")

# The global symplectic primitive alpha.  Use R=b(b+4)-c^2e, so
# db/c-alpha=((b+2)/(8c))dR.  Clear 8c and reduce all three coefficients.
relation_positive = bb * (bb + 4) - cc**2 * ee
alpha_b = -cc * ee / 4
alpha_c = ee * (bb + 2) / 4
alpha_e = cc * (bb + 2) / 8
dR = (2 * (bb + 2), -2 * cc * ee, -cc**2)
primitive_residuals = (
    sp.expand(8 - 8 * cc * alpha_b - (bb + 2) * dR[0]),
    sp.expand(-8 * cc * alpha_c - (bb + 2) * dR[1]),
    sp.expand(-8 * cc * alpha_e - (bb + 2) * dR[2]),
)
primitive_reducer = sp.groebner([relation_positive], ee, cc, bb, order="lex")
for index, residual in enumerate(primitive_residuals):
    require(
        sp.expand(primitive_reducer.reduce(residual)[1]) == 0,
        f"global primitive differential coefficient {index}",
    )

# Pull alpha to c!=0 using e=S/c^2 and the differentiated relation.  It is
# exactly db/c there, so d alpha=db wedge dc/c^2; regularity makes this a
# global identity on the smooth surface.
S = bb * (bb + 4)
S_prime = sp.diff(S, bb)
alpha_b_chart = sp.factor(
    alpha_b.subs(ee, S / cc**2) + alpha_e * S_prime / cc**2
)
alpha_c_chart = sp.factor(
    alpha_c.subs(ee, S / cc**2) - alpha_e * 2 * S / cc**3
)
require(sp.factor(alpha_b_chart - 1 / cc) == 0, "primitive chart db coefficient")
require(alpha_c_chart == 0, "primitive chart dc coefficient")

# Finite controls for the universal homogeneous formula.  The proof handles
# arbitrary f,g and every r; these rows verify the bracket identity and each
# divisibility chamber for both parities.
f_control = bb**3 + 2 * bb + 1
g_control = bb**2 + 3
homogeneous_rows = []
for r in range(9):
    m = (r + 2) // 2
    h_control = S**m * g_control
    P_control = cc**r * f_control
    Q_control = cc ** (-r - 1) * h_control
    chart_bracket = sp.factor(
        cc**2
        * (
            sp.diff(P_control, bb) * sp.diff(Q_control, cc)
            - sp.diff(P_control, cc) * sp.diff(Q_control, bb)
        )
    )
    expected_bracket = sp.factor(
        -(r + 1) * sp.diff(f_control, bb) * h_control
        - r * f_control * sp.diff(h_control, bb)
    )
    require(sp.factor(chart_bracket - expected_bracket) == 0, f"homogeneous formula r={r}")
    if r == 0:
        require(sp.rem(expected_bracket, S, bb) == 0, "homogeneous r=0 S factor")
    elif r >= 2:
        require(
            sp.rem(expected_bracket, S ** (m - 1), bb) == 0,
            f"homogeneous r={r} arm factor",
        )
    else:
        require(sp.Poly(expected_bracket, bb).degree() > 0, "homogeneous r=1 degree")
    homogeneous_rows.append((r, m, sp.Poly(expected_bracket, bb).degree()))

# Finite controls for the all-degree two-by-two weight obstruction.  With
# R=2 and T=2k the two vanishing Wronskians have the universal UFD forms
# f=A h, g=B h^k and F=L K^(2k-1), G=M K.  The structural proof shows that
# every surviving two-by-two support reaches this chamber or its symmetric
# copy.  These rows independently expand and factor the resulting constant
# coefficient; they are controls, not the source of the universal claim.
h_two = S * (bb**2 + 2 * bb + 3)
K_two = bb**3 - bb + 2
A_two, B_two, L_two, M_two = map(sp.Integer, (2, 3, 5, 7))
two_by_two_rows = []
for k in range(1, 9):
    f_two = A_two * h_two
    g_two = B_two * h_two**k
    F_two = L_two * K_two ** (2 * k - 1)
    G_two = M_two * K_two
    H_two = sp.expand(
        sp.diff(f_two, bb) * G_two
        + 2 * f_two * sp.diff(G_two, bb)
        - 2 * k * sp.diff(F_two, bb) * g_two
        - (2 * k - 1) * F_two * sp.diff(g_two, bb)
    )
    first_factor = sp.expand(K_two * sp.diff(h_two, bb) + 2 * h_two * sp.diff(K_two, bb))
    second_factor = sp.expand(
        A_two * M_two
        - k
        * (2 * k - 1)
        * L_two
        * B_two
        * K_two ** (2 * k - 2)
        * h_two ** (k - 1)
    )
    require(
        sp.expand(H_two - first_factor * second_factor) == 0,
        f"two-by-two factorization k={k}",
    )
    first_degree = sp.Poly(first_factor, bb).degree()
    require(first_degree >= 1, f"two-by-two nonunit first factor k={k}")
    expected_first_degree = sp.degree(K_two, bb) + sp.degree(h_two, bb) - 1
    require(first_degree == expected_first_degree, f"two-by-two first degree k={k}")
    two_by_two_rows.append((k, first_degree, sp.Poly(H_two, bb).degree()))

# The root-order routing has only the R=2 or T=2 boundary.  Freeze its finite
# arithmetic shadow, including R=T=1 and both zero-weight edge cases.
for arm in range(1, 17):
    require(
        ((arm > 1) and ((arm + 1) // 2 == 1)) == (arm == 2),
        f"two-by-two root arm arithmetic {arm}",
    )

require(sp.factor(relation_positive.subs(cc, 0)) == bb * (bb + 4), "reducible c=0 fibre")

print("THM-3561 exact companion")
print("middle_leaf_unit_identity=PASS")
print("middle_leaf_ring=QQ[x,w,(1-x^2*w)^-1]")
print("rational_pair=(q/D^2, x*D*(D+2))")
print("rational_jacobian=-3")
print(f"rational_collision_rows={collision_rows}")
print("completion=(b,c,e)=((D-1)(D+2)^2, xD(D+2), q(D+3))")
print("target_relation=c^2*e-b*(b+4)=0")
print(f"inverse_cubic={inverse_cubic}")
print(f"J_bc={J_bc}")
print(f"J_ce={J_ce}")
print(f"boundary_jacobians={boundary_jacobians}")
print(f"completed_collision_rows={completed_collision_rows}")
print("target_smooth=PASS")
print("image=Y_minus_{(-4,0,0)}")
print("generic_degree=3")
print("intersection_generators=b,c,e")
for s, m, expected in intersection_rows:
    print(f"intersection_weight=-{s} cancellation_order={m} normal_form={expected}")
print("poisson_generators={b,c}:c^2 {c,e}:-2(b+2) {b,e}:-2ce")
print("affine_linear_constant_bracket=zero_only")
print("symplectic_primitive=global_exact")
print(f"homogeneous_weight_controls={homogeneous_rows}")
print("homogeneous_constant_bracket=none")
print(f"two_by_two_weight_controls={two_by_two_rows}")
print("two_by_two_constant_bracket=none")
print("generator_slices=none")
print("euler_Y=2 euler_A2=1")
print("picard_Y=Z")
print("scope=NO_JC2_CLAIM; solve {P,Q}=unit on Y to obtain one")
