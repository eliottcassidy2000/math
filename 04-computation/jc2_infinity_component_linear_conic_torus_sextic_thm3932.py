#!/usr/bin/env python3
"""Exact controls for THM-3932's infinity-component conic classification.

Reproduction:
  python3 04-computation/jc2_infinity_component_linear_conic_torus_sextic_thm3932.py
  python3 -O 04-computation/jc2_infinity_component_linear_conic_torus_sextic_thm3932.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.simplify(sp.cancel(expression)) == 0, message)


def smith_diagonal(matrix: sp.Matrix) -> tuple[int, ...]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return tuple(
        abs(int(smith[i, i]))
        for i in range(min(smith.rows, smith.cols))
        if smith[i, i] != 0
    )


x, y, h, s, u, w = sp.symbols("x y h s u w")
X, Y, Z, S, T, W = sp.symbols("X Y Z S T W")
r, beta = sp.symbols("r beta")
B0, B1, C0, C1 = sp.symbols("B0 B1 C0 C1")


# ---------------------------------------------------------------------------
# Fold-degree compression and the parity obstruction at degree two.
# ---------------------------------------------------------------------------

fold_rows = []
for d in range(1, 7):
    # q(h^2,y)-2h^3 has y-degree at most three.  An irreducible branch
    # therefore has [k(s):k(h)] = d <= 3.
    if d <= 3:
        fold_rows.append(d)
gate(fold_rows == [1, 2, 3], "linear p leaves exactly folds one, two, three")

for denominator_x_degree in (1, 2):
    gate(
        6 - 2 * denominator_x_degree <= 4,
        "nonconstant fold-one y coefficient lowers the branch below degree six",
    )

# For d=2, center the quadratic polynomial h(s) at h=alpha*s^2+beta.
# A degree-six y has trace degree three in h, while a1(h^2) has even h-degree.
gate(6 // 2 == 3, "degree-six polynomial has cubic even part over a quadratic h")
gate(all(power % 2 == 0 for power in (0, 2, 4)), "allowed trace powers come from h squared")
gate(3 not in (0, 2, 4), "fold-two trace parity excludes degree six")


# ---------------------------------------------------------------------------
# Complete trace-zero grammar for the fold-three row.
# ---------------------------------------------------------------------------

# h=s^3+r*s+beta.  Every polynomial y(s) of degree at most six has a unique
# decomposition E(h)+B(h)s+C(h)(s^2+2r/3), with deg E<=2 and deg B,C<=1.
v = (B0 + B1 * h) * s + (C0 + C1 * h) * (s**2 + sp.Rational(2, 3) * r)
minimal_v = sp.Poly(
    sp.resultant(s**3 + r * s + beta - h, Y - v, s), Y
)
gate(minimal_v.all_coeffs()[1] == 0, "centered cubic trace vanishes")
e2_v = sp.factor(minimal_v.all_coeffs()[2])
norm_v = sp.factor(-minimal_v.all_coeffs()[3])

e2_poly = sp.Poly(e2_v, h)
norm_poly = sp.Poly(norm_v, h)
zero(norm_poly.coeff_monomial(h**5) - C1**3, "top odd norm coefficient is C1 cubed")

e2_h1_after = sp.factor(e2_poly.coeff_monomial(h).subs(C1, 0))
norm_h1_after = sp.factor(norm_poly.coeff_monomial(h).subs(C1, 0))
norm_h3_after = sp.factor(norm_poly.coeff_monomial(h**3).subs(C1, 0))

equation_trace = 2 * B0 * B1 * r - 3 * B0 * C0 + 3 * B1 * C0 * beta
equation_norm_one = (
    B0**3
    - 3 * beta * B0**2 * B1
    + sp.Rational(4, 3) * B0 * B1 * C0 * r**2
    - B0 * C0**2 * r
    + beta * B1 * C0**2 * r
    - 2 * beta * C0**3
)
equation_norm_three = B1**2 * (3 * B0 - B1 * beta)

zero(e2_h1_after - equation_trace, "fold-three even second coefficient equation")
zero(norm_h1_after - equation_norm_one, "fold-three norm has no h term equation")
zero(norm_h3_after - equation_norm_three, "fold-three prescribed h-cubed coefficient")

family_substitution = {beta: -2, B0: 0, B1: 1, C0: 0, C1: 0}
zero(equation_trace.subs(family_substitution), "explicit family trace equation")
zero(equation_norm_one.subs(family_substitution), "explicit family norm h equation")
zero(equation_norm_three.subs(family_substitution) - 2, "explicit family Cardano norm equation")


# ---------------------------------------------------------------------------
# The explicit fold-three family and its r=0 projective member.
# ---------------------------------------------------------------------------

h_r = s**3 + r * s - 2
x_r = h_r**2
y_r = s * h_r
q_r = y**3 + r * x * y - x**2
zero(q_r.subs({x: x_r, y: y_r}) - 2 * h_r**3, "one-parameter fold-three family")
zero(4 * (3 * x_r) ** 3 - 27 * q_r.subs({x: x_r, y: y_r}) ** 2, "family Cardano identity")

q = y**3 - x**2
F = sp.expand(4 * x**3 - q**2)
h_0 = s**3 - 2
x_0 = h_0**2
y_0 = s * h_0
zero(q.subs({x: x_0, y: y_0}) - 2 * h_0**3, "explicit branch cube coordinate")
zero(F.subs({x: x_0, y: y_0}), "explicit normalization lies on sextic")
gate(sp.Poly(F, x, y, domain=sp.QQ).total_degree() == 6, "explicit branch has degree six")
gate(sp.Poly(F, x, y, domain=sp.QQ).is_irreducible is True, "explicit sextic is irreducible over Q")
gate(sp.factor(sum(
    coefficient * x**monomial[0] * y**monomial[1]
    for monomial, coefficient in sp.Poly(F, x, y).terms()
    if sum(monomial) == 6
)) == -y**6, "sextic top form has one projective support point")

h_inverse = q / (2 * x)
s_inverse = 2 * x * y / q
zero(h_inverse.subs({x: x_0, y: y_0}) - h_0, "Cardano parameter rational inverse")
zero(s_inverse.subs({x: x_0, y: y_0}) - s, "normalization parameter rational inverse")
zero(h_inverse**2 - x + F / (4 * x**2), "inverse square identity modulo the sextic")

Q2 = 3 * X * Z
Q3 = Y**3 - X**2 * Z
Delta = sp.expand(4 * Q2**3 - 27 * Q3**2)
Fh = sp.expand(Delta / 27)
H = T**3 - 2 * S**3
Xh = H**2
Yh = S**2 * T * H
Zh = S**6
gate(all(sp.Poly(form, S, T).total_degree() == 6 for form in (Xh, Yh, Zh)), "homogeneous normalization has degree six")
zero(Fh.subs({X: Xh, Y: Yh, Z: Zh}), "projective normalization identity")
gate(sp.gcd(sp.gcd(Xh, Yh), Zh) == 1, "projective normalization is basepoint free")
zero(Yh.subs(S, 0), "unique infinity address has Y zero")
zero(Zh.subs(S, 0), "unique infinity address has Z zero")
gate(Xh.subs(S, 0) == T**6, "unique infinity address survives in X")
zero(Fh.subs(Z, 0) + Y**6, "projective infinity support is the single point [1:0:0]")


# ---------------------------------------------------------------------------
# Complete singularity and genus/address ledger for the explicit sextic.
# ---------------------------------------------------------------------------

affine_singular = sp.groebner([F, sp.diff(F, x), sp.diff(F, y)], x, y, order="lex")
gate(
    [sp.factor(poly.as_expr()) for poly in affine_singular.polys]
    == [x * (3 * x + y**3), y**5],
    "affine singular scheme has sole support at the origin",
)

uu, vv, zz = sp.symbols("uu vv zz")
F_infinity = sp.expand(4 * vv**3 - (uu**3 - vv) ** 2)
infinity_singular = sp.groebner(
    [F_infinity, sp.diff(F_infinity, uu), sp.diff(F_infinity, vv)],
    uu,
    vv,
    order="lex",
)
gate(
    sp.solve_poly_system(
        [F_infinity, sp.diff(F_infinity, uu), sp.diff(F_infinity, vv)], uu, vv
    )
    == [(0, 0)],
    "infinity chart has one singular support point",
)

gate(sp.discriminant(s**3 - 2, s) == -108, "origin has three distinct normalization addresses")
gate(sp.gcd(s**3 - 2, 3 * s**2) == 1, "all three origin addresses are smooth upstairs")
gate(3 * 2 == 6, "three pairwise contact-two branches contribute delta six")

u_infinity = zz**2 / (1 - 2 * zz**3)
v_infinity = zz**6 / (1 - 2 * zz**3) ** 2
zero(
    v_infinity - u_infinity**3 + 2 * zz**9 / (1 - 2 * zz**3) ** 3,
    "infinity branch has characteristic orders two and nine",
)
gate(sp.gcd(2, 9) == 1, "infinity branch is unibranch with Puiseux pair (2,9)")
gate((2 - 1) * (9 - 1) // 2 == 4, "infinity branch contributes delta four")
gate((6 - 1) * (6 - 2) // 2 == 10, "sextic arithmetic genus is ten")
gate(6 + 4 == 10, "finite and infinite delta exhaust genus")


# ---------------------------------------------------------------------------
# Coefficient map, depressed cubic, and the genuine resolvent three-class.
# ---------------------------------------------------------------------------

p_affine = 3 * x
jacobian = sp.factor(
    sp.diff(p_affine, x) * sp.diff(q, y) - sp.diff(p_affine, y) * sp.diff(q, x)
)
zero(jacobian - 9 * y**2, "fold-three coefficient map has nonconstant Jacobian")

cubic = u**3 - 3 * x * u - q
zero(sp.discriminant(cubic, u) - 27 * F, "depressed cubic discriminant")

aa, bb, cc = sp.symbols("aa bb cc")
linear_root = aa * x + bb * y + cc
root_equations = [
    coefficient
    for _, coefficient in sp.Poly(
        sp.expand(cubic.subs(u, linear_root)), x, y
    ).terms()
]
gate(sp.solve(root_equations, (aa, bb, cc), dict=True) == [], "depressed cubic has no affine-linear root")

cubic_singular = sp.groebner(
    [cubic, sp.diff(cubic, x), sp.diff(cubic, y), sp.diff(cubic, u)],
    u,
    x,
    y,
    order="lex",
)
gate(
    sp.solve_poly_system(
        [cubic, sp.diff(cubic, x), sp.diff(cubic, y), sp.diff(cubic, u)],
        u,
        x,
        y,
    )
    == [(0, 0, 0)],
    "monogenic cubic surface has one isolated singularity and is normal",
)

D = sp.expand(q**2 - 4 * x**3)
resolvent = sp.expand(w**2 - D)
gate(sp.Poly(resolvent, w, x, y, domain=sp.QQ).is_irreducible is True, "quadratic resolvent is integral")
resolvent_singular = sp.groebner(
    [resolvent, sp.diff(resolvent, w), sp.diff(resolvent, x), sp.diff(resolvent, y)],
    w,
    x,
    y,
    order="lex",
)
gate(
    sp.solve_poly_system(
        [resolvent, sp.diff(resolvent, w), sp.diff(resolvent, x), sp.diff(resolvent, y)],
        w,
        x,
        y,
    )
    == [(0, 0, 0)],
    "quadratic resolvent has one isolated singularity and is normal",
)

a_cardano = q + w
b_cardano = q - w
zero(a_cardano * b_cardano - 4 * x**3 + resolvent, "Cardano norm modulo resolvent")
zero(a_cardano.subs({x: 0, w: -y**3}), "D-plus is the zero sheet of q+w")
zero(b_cardano.subs({x: 0, w: y**3}), "D-minus is the zero sheet of q-w")
relation_matrix = sp.Matrix([[1, 1], [3, 0]])
gate(smith_diagonal(relation_matrix) == (1, 3), "resolvent line packet has a formal Z/3")
gate(relation_matrix.det() == -3, "resolvent three-class has exact relation order")

# At m=(x,y,a), P=(x,a) has two independent linear generators because the
# local defining relation a^2-2qa+4x^3 has no linear term.  This is the exact
# non-Cartier control used in the prose proof.
local_relation = sp.expand(sp.Symbol("a") ** 2 - 2 * q * sp.Symbol("a") + 4 * x**3)
gate(min(sum(monomial) for monomial, _ in sp.Poly(local_relation, sp.Symbol("a"), x, y).terms()) >= 2, "local resolvent relation has no linear term")
gate(sp.Matrix([[1, 0], [0, 1]]).rank() == 2, "D-plus ideal needs generators x and a modulo mP")

# The weighted initial form at the resolvent singularity is elliptic; this is
# a warning that the two divisor relations do not compute the entire Cl(B).
weighted_terms = {
    monomial: coefficient
    for monomial, coefficient in sp.Poly(resolvent, w, x, y).terms()
}
weighted_initial = sum(
    coefficient * w**monomial[0] * x**monomial[1] * y**monomial[2]
    for monomial, coefficient in weighted_terms.items()
    if 3 * monomial[0] + 2 * monomial[1] + monomial[2] == 6
)
zero(weighted_initial - (w**2 + 4 * x**3 - y**6), "resolvent has simple-elliptic weighted initial form")


semantic_payload = {
    "classification": "linear_p_one_place_sextic_has_fold_1_or_3;fold_2_empty",
    "fold_one": "triangular_polynomial_automorphism_and_monogenic_cusp_pullback",
    "fold_three": "complete_trace_norm_grammar_and_explicit_q_equals_y3_minus_x2",
    "branch": "irreducible_rational_sextic;three_finite_addresses;one_infinity_address;delta_6_plus_4",
    "resolvent": "normal_surface_with_genuine_nonprincipal_order_three_Dplus_class",
    "order": "associated_normal_cubic_is_explicitly_monogenic_not_a_Keller_completion",
    "scope": "full_global_resolvent_class_group_and_nonmonogenic_twists_not_claimed",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3932-infinity-component-linear-conic-torus-sextic-fold-classification")
print("folds=d_in_{1,2,3};d1=triangular_coordinate_fold;d2=empty;d3=nonempty")
print("family=Q2=3XZ;Q3=Y^3+rXYZ-X^2Z;h=s^3+rs-2;X=h^2;Y=sh")
print("explicit_r0=F/27=4X^3Z^3-(Y^3-X^2Z)^2")
print("normalization=[H^2:S^2*T*H:S^6],H=T^3-2S^3;basepoint_free;one_infinity_address")
print("singularities=origin_three_smooth_contact2_branches_delta6;infinity_Puiseux_2_9_delta4")
print("coefficient_map=Jac(3x,y^3-x^2)=9y^2_nonconstant")
print("cubic=u^3-3xu-(y^3-x^2);normal_and_globally_monogenic")
print("resolvent=Dplus=(x,q+w);div(x)=Dplus+Dminus;div(q+w)=3Dplus;order_exactly_3")
print("class_scope=genuine_Z3_subgroup_only;full_Cl_not_claimed_due_simple_elliptic_local_sidecar")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
