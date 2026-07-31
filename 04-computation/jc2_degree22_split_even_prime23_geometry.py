#!/usr/bin/env python3
"""Exact geometry scout for the split even-Faber prime-23 curve.

This script builds on the frozen lambda scout and studies

    f2(t,v,zeta)=0,       zeta*f1(t,v,zeta)^4=eta*t^23.

It separates three claims.

* GENERIC-EXACT: the infinity resultant over Q(c,d,e,w,eta) is an
  irreducible degree-23 polynomial with constant leading coefficient.
* SPECIALIZED-EXACT: at c=d=e=w=eta=1 the full (t,v) eliminant is
  Q-irreducible, the degree-six weighted surface is smooth, and all 23
  infinity points are transverse.
* THEORETICAL GENERIC CONSEQUENCE: weighted adjunction plus Bertini and the
  five universal (4,23) cusps give normalization genus 89 on a nonempty
  Zariski-open parameter set.  This is not a uniform statement for every
  coefficient specialization and not a closure of the split branch.

The q-plane birational model is also reconstructed.  Its Newton polygon has
158 interior points, and an exact subresultant/Hessian audit finds precisely
69 ordinary all-one torus nodes and no other torus singular value.  Thus this
second model independently gives normalization genus 89.
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
from pathlib import Path

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "jc2_degree22_split_lambda_prime23_scout.py"
BASE_SHA256 = "1bbadb900e27112f57f600b0b5b73a3b85fc5b02f23e3f8f687dac4fa1c41fc3"
require(
    hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256,
    "audited split-lambda scout changed",
)

spec = importlib.util.spec_from_file_location("split_lambda_prime23", BASE_PATH)
require(spec is not None and spec.loader is not None, "cannot load lambda scout")
base = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(base)


def positive_primitive(poly: s.Poly) -> s.Poly:
    """Return the primitive integral associate with positive leading term."""
    _, primitive = s.primitive(poly.as_expr(), *poly.gens)
    result = s.Poly(primitive, *poly.gens, domain=s.QQ)
    if result.LC() < 0:
        result = -result
    return result


def polygon_hull(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Monotone-chain convex hull, without repeating the initial vertex."""
    points = sorted(set(points))

    def cross(o: tuple[int, int], a: tuple[int, int], b: tuple[int, int]) -> int:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[tuple[int, int]] = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def face_polynomial(
    poly: s.Poly, start: tuple[int, int], end: tuple[int, int]
) -> s.Poly:
    """Univariate polynomial of a Newton edge after stripping its monomial."""
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = int(s.gcd(abs(dx), abs(dy)))
    step = (dx // length, dy // length)
    coefficients = dict(poly.terms())
    X = s.symbols("X")
    expression = sum(
        coefficients.get((start[0] + k * step[0], start[1] + k * step[1]), 0)
        * X**k
        for k in range(length + 1)
    )
    return s.Poly(expression, X, domain=s.QQ)


t, v, zeta = base.tau, base.v, base.zeta
c, d, e, w = base.c, base.dd, base.e, base.w
f1, f2 = base.f1, base.f2
eta = s.symbols("eta")

# ---------------------------------------------------------------------------
# 1. Weighted complete-intersection closure and a smooth surface witness.
# ---------------------------------------------------------------------------

h = s.symbols("h")


def weighted_homogenize(expression: s.Expr, degree: int) -> s.Expr:
    result = 0
    for (it, iv, iz), coefficient in s.Poly(expression, t, v, zeta).terms():
        weight = it + 2 * iv + 3 * iz
        require(weight <= degree, "weighted degree bound changed")
        result += coefficient * t**it * v**iv * zeta**iz * h ** (degree - weight)
    return s.expand(result)


F1 = weighted_homogenize(f1, 5)
F2 = weighted_homogenize(f2, 6)
H = s.expand(zeta * F1**4 - eta * t**23)
require(
    all(
        ih + it + 2 * iv + 3 * iz == degree
        for polynomial, degree in ((F1, 5), (F2, 6), (H, 23))
        for (ih, it, iv, iz), _ in s.Poly(polynomial, h, t, v, zeta).terms()
    ),
    "weighted homogenization is not homogeneous",
)

ones = {c: 1, d: 1, e: 1, w: 1}
f2_one = s.expand(f2.subs(ones))
f2_one_factors = s.factor_list(f2_one)[1]
require(
    len(f2_one_factors) == 1 and f2_one_factors[0][1] == 1,
    "all-one degree-six surface became reducible",
)
affine_surface_jacobian = [
    f2_one,
    s.diff(f2_one, t),
    s.diff(f2_one, v),
    s.diff(f2_one, zeta),
]
affine_surface_gb = s.groebner(
    affine_surface_jacobian, zeta, v, t, order="grevlex"
)
require(affine_surface_gb.contains(s.Integer(1)), "all-one affine surface is singular")

F2_one = s.expand(F2.subs(ones))
F2_infinity = s.expand(F2_one.subs(h, 0))
infinity_partials = [s.diff(F2_one, x).subs(h, 0) for x in (h, t, v, zeta)]
for patch_variable in (t, v, zeta):
    patch_gb = s.groebner(
        [F2_infinity, *infinity_partials, patch_variable - 1],
        zeta,
        v,
        t,
        order="grevlex",
    )
    require(patch_gb.contains(s.Integer(1)), "all-one surface singular at infinity")

# The only singular points of P(1,1,2,3) are the v- and zeta-coordinate
# points.  The fixed v^3 and zeta^2 coefficients keep the surface away from
# both, uniformly in the four parameters.
require(F2.subs({h: 0, t: 0, v: 1, zeta: 0}) != 0, "surface met weight-2 point")
require(F2.subs({h: 0, t: 0, v: 0, zeta: 1}) != 0, "surface met weight-3 point")

# The pencil generated on the surface by zeta*F1^4 and t^23 has no base
# point on the chart h=t=0.  (The affine h=1 base scheme is exactly the
# already-audited G3 plus L5 section.)  Checking the two weighted patches is
# enough because v=zeta=0 is not a projective point.
F1_corner = s.expand(F1.subs({h: 0, t: 0}))
F2_corner = s.expand(F2.subs({h: 0, t: 0}))
for patch_variable in (v, zeta):
    corner_gb = s.groebner(
        [F2_corner, zeta * F1_corner**4, patch_variable - 1],
        zeta,
        v,
        order="grevlex",
    )
    require(corner_gb.contains(s.Integer(1)), "pencil gained a corner base point")

# ---------------------------------------------------------------------------
# 2. Exact specialized full eliminant and its arithmetic/geometric scope.
# ---------------------------------------------------------------------------

f1_one = s.expand(f1.subs(ones))
H_affine_one = s.expand(zeta * f1_one**4 - t**23)
R_tv_raw = s.Poly(s.resultant(f2_one, H_affine_one, zeta), t, v, domain=s.QQ)
R_tv = positive_primitive(R_tv_raw)
require((R_tv.degree(t), R_tv.degree(v)) == (46, 23), "global bidegree changed")
require(len(R_tv.terms()) == 553, "global eliminant support changed")
require(
    not s.Poly(R_tv.as_expr(), v).LC().free_symbols,
    "global eliminant lost its constant v-leading coefficient",
)
require(
    all(it + 2 * iv <= 46 for (it, iv), _ in R_tv.terms()),
    "global Newton weight escaped triangle",
)
require(
    polygon_hull([monomial for monomial, _ in R_tv.terms()])
    == [(0, 0), (46, 0), (0, 23)],
    "global Newton triangle changed",
)
R_tv_factors = s.factor_list(R_tv.as_expr())[1]
require(
    len(R_tv_factors) == 1
    and R_tv_factors[0][1] == 1
    and s.Poly(R_tv_factors[0][0], t, v).degree(t) == 46
    and s.Poly(R_tv_factors[0][0], t, v).degree(v) == 23,
    "all-one full eliminant became reducible over Q",
)

# The resultant occurs to exponent one and the penultimate subresultant is
# linear in zeta.  Since its two coefficients have v-degree below 23, they
# cannot both vanish on the irreducible R_tv divisor.  Thus f2 and H have
# generic gcd degree one over Q(R_tv), proving that projection to (t,v) is
# birational rather than a hidden two-sheet quotient.
zeta_subresultants = s.subresultants(f2_one, H_affine_one, zeta)
require(
    [s.Poly(polynomial, zeta).degree() for polynomial in zeta_subresultants]
    == [5, 2, 1, 0],
    "zeta subresultant degree ladder changed",
)
linear_subresultant = s.Poly(zeta_subresultants[-2], zeta)
require(
    max(
        s.Poly(coefficient, t, v).degree(v)
        for coefficient in linear_subresultant.all_coeffs()
    )
    == 11,
    "linear zeta subresultant could vanish generically on R_tv",
)

fixed_tv = positive_primitive(s.Poly(R_tv.as_expr().subs(t, 0), v, domain=s.QQ))
require(
    fixed_tv.as_expr() == s.expand(base.G3 * base.L5**4),
    "all-one fixed fibre changed",
)

# ---------------------------------------------------------------------------
# 3. Generic infinity face and exact all-one transversality.
# ---------------------------------------------------------------------------

V, S = s.symbols("V S")
A_infinity = (
    9370240 * c * V
    - 14992384 * e
    + 819896 * S
    - 1449459 * V * S
)
B_infinity = (
    15944049 * S**2
    - 206145280 * c * S
    - 1978994688 * d * V
    - 1319329792 * w
    + 1443016960 * V**2
    - 1190488992 * V**3
)
H_infinity = s.expand(S * A_infinity**4 - eta)
R_infinity_raw = s.resultant(B_infinity, H_infinity, S)
R_infinity_content, R_infinity_primitive_expr = s.primitive(
    R_infinity_raw, V, c, d, e, w, eta
)
R_infinity_generic = s.Poly(
    R_infinity_primitive_expr, V, c, d, e, w, eta, domain=s.QQ
)
if R_infinity_generic.LC() < 0:
    R_infinity_generic = -R_infinity_generic
require(R_infinity_content == 11**30, "generic infinity content changed")
require(R_infinity_generic.degree(V) == 23, "generic infinity degree changed")
require(len(R_infinity_generic.terms()) == 1416, "generic infinity support changed")
require(
    s.Poly(R_infinity_generic.as_expr(), V).LC()
    == 2**25 * 3**21 * 7**5 * 11**40,
    "generic infinity leading coefficient changed",
)
generic_infinity_factors = s.factor_list(R_infinity_generic.as_expr())[1]
require(
    len(generic_infinity_factors) == 1
    and generic_infinity_factors[0][1] == 1
    and s.Poly(generic_infinity_factors[0][0], V).degree() == 23,
    "generic infinity resultant became reducible",
)

infinity_one_substitution = {c: 1, d: 1, e: 1, w: 1, eta: 1}
P_infinity = positive_primitive(
    s.Poly(R_infinity_raw.subs(infinity_one_substitution), V, domain=s.QQ)
)
require(P_infinity.degree() == 23, "all-one infinity degree changed")
mod17_factors = s.factor_list(P_infinity.as_expr(), modulus=17)[1]
require(
    len(mod17_factors) == 1
    and mod17_factors[0][1] == 1
    and s.Poly(mod17_factors[0][0], V, modulus=17).degree() == 23,
    "all-one infinity polynomial is not irreducible mod 17",
)

A_one = A_infinity.subs(infinity_one_substitution)
B_one = B_infinity.subs(infinity_one_substitution)
H_infinity_one = H_infinity.subs(infinity_one_substitution)
infinity_jacobian = s.det(
    s.Matrix(
        [
            [s.diff(B_one, V), s.diff(B_one, S)],
            [s.diff(H_infinity_one, V), s.diff(H_infinity_one, S)],
        ]
    )
)
infinity_jacobian_resultant = s.Poly(
    s.resultant(B_one, infinity_jacobian, S), V, domain=s.QQ
)
require(
    s.gcd(P_infinity, infinity_jacobian_resultant).degree() == 0,
    "an all-one infinity point is not transverse",
)

# ---------------------------------------------------------------------------
# 4. Birational q-plane model and the projection-genus loss.
# ---------------------------------------------------------------------------

q = s.symbols("q")
q_substitution = {**ones, zeta: q**4 * t**3}
G1_q = s.expand(q * f1.subs(q_substitution) + t**5)
G2_q = s.expand(f2.subs(q_substitution))
R_tq = positive_primitive(
    s.Poly(s.resultant(G1_q, G2_q, v), t, q, domain=s.QQ)
)
require((R_tq.degree(t), R_tq.degree(q)) == (15, 23), "q-plane bidegree changed")
require(len(R_tq.terms()) == 75, "q-plane support changed")
q_hull = polygon_hull([monomial for monomial, _ in R_tq.terms()])
require(q_hull == [(0, 3), (15, 0), (15, 23)], "q-plane Newton hull changed")
R_tq_factors = s.factor_list(R_tq.as_expr())[1]
require(
    len(R_tq_factors) == 1
    and R_tq_factors[0][1] == 1
    and s.Poly(R_tq_factors[0][0], t, q).degree(q) == 23,
    "q-plane eliminant became reducible over Q",
)
v_subresultants = s.subresultants(G1_q, G2_q, v)
require(
    [s.Poly(polynomial, v).degree() for polynomial in v_subresultants]
    == [3, 2, 1, 0],
    "q-plane v-subresultant degree ladder changed",
)
v_linear_subresultant = s.Poly(v_subresultants[-2], v)
v_linear_content = s.factor(s.gcd(*v_linear_subresultant.all_coeffs()))
require(
    v_linear_content != 0
    and q in v_linear_content.free_symbols
    and not s.factor(v_linear_content / q).free_symbols,
    "q-plane linear subresultant lost its sole q content",
)

# The vanishing of the normalized linear-subresultant coefficients is the
# double-v collision locus.  Eliminating t splits it into a contracted q=0
# boundary factor, one extraneous quartic boundary factor, and a squarefree
# degree-69 torus factor.  The latter is the exact node polynomial.
v_coefficient, constant_coefficient = v_linear_subresultant.all_coeffs()
collision_a = s.cancel(v_coefficient / v_linear_content)
collision_b = s.cancel(constant_coefficient / v_linear_content)
collision_resultant = s.factor_list(s.resultant(collision_a, collision_b, t), q)
collision_factor_degrees = sorted(
    (s.Poly(factor, q).degree(), multiplicity)
    for factor, multiplicity in collision_resultant[1]
)
require(
    collision_factor_degrees == [(1, 5), (4, 1), (69, 1)],
    "q-plane collision factor degrees changed",
)
quartic_collision_factor = next(
    s.Poly(factor, q, domain=s.QQ)
    for factor, multiplicity in collision_resultant[1]
    if s.Poly(factor, q).degree() == 4 and multiplicity == 1
)
require(
    s.factor(quartic_collision_factor.as_expr() - (99 * q**4 - 640)) == 0
    or s.factor(quartic_collision_factor.as_expr() + (99 * q**4 - 640)) == 0,
    "extraneous collision quartic changed",
)
Q69 = positive_primitive(
    next(
        s.Poly(factor, q, domain=s.QQ)
        for factor, multiplicity in collision_resultant[1]
        if s.Poly(factor, q).degree() == 69 and multiplicity == 1
    )
)
require(s.gcd(Q69, Q69.diff()).degree() == 0, "Q69 is not squarefree")
Q69_SHA256 = hashlib.sha256(s.srepr(Q69.as_expr()).encode()).hexdigest()
require(
    Q69_SHA256 == "239dcceebc110b863b8008071e66835264406bf39346d61654acbd409e21d9bf",
    "Q69 canonical digest changed",
)

# At every Q69 collision: the two v-preimages are finite and distinct; the
# collision point in (t,q) is transverse; and the plane Hessian is nonzero.
# The four coprime resultants certify those assertions simultaneously over
# the irreducible Q69 field, without enumerating 69 algebraic roots.
G1_q_poly = s.Poly(G1_q, v)
G1_q_discriminant = s.discriminant(G1_q_poly.as_expr(), v)
collision_jacobian = s.expand(
    s.diff(collision_a, t) * s.diff(collision_b, q)
    - s.diff(collision_a, q) * s.diff(collision_b, t)
)
q_plane_hessian = s.expand(
    s.diff(R_tq.as_expr(), t, 2) * s.diff(R_tq.as_expr(), q, 2)
    - s.diff(R_tq.as_expr(), t, q) ** 2
)
node_controls = (
    G1_q_discriminant,
    collision_jacobian,
    G1_q_poly.LC(),
    q_plane_hessian,
)
node_control_degrees = []
for control in node_controls:
    control_resultant = s.Poly(s.resultant(collision_a, control, t), q, domain=s.QQ)
    node_control_degrees.append(control_resultant.degree())
    require(
        s.gcd(Q69, control_resultant).degree() == 0,
        "a Q69 collision failed an ordinary-node control",
    )
require(node_control_degrees == [66, 131, 6, 292], "node-control degrees changed")

# Completeness of the affine torus singularity census.  A common singular
# q-value must divide both eliminants below.  Their exact gcd is q^30*Q69^2;
# q=0 is the contracted boundary already handled by the edge audit, leaving
# precisely the 69 ordinary Q69 nodes in t*q != 0.
singular_t_eliminant = positive_primitive(
    s.Poly(s.resultant(R_tq.as_expr(), s.diff(R_tq.as_expr(), t), t), q)
)
singular_q_eliminant = positive_primitive(
    s.Poly(s.resultant(R_tq.as_expr(), s.diff(R_tq.as_expr(), q), t), q)
)
singular_value_gcd = positive_primitive(s.gcd(singular_t_eliminant, singular_q_eliminant))
singular_value_factors = s.factor_list(singular_value_gcd.as_expr(), q)[1]
require(
    [(s.Poly(factor, q).degree(), multiplicity) for factor, multiplicity in singular_value_factors]
    == [(1, 30), (69, 2)],
    "q-plane singular-value census changed",
)
require(
    s.Poly(singular_value_factors[0][0], q).monic().as_expr() == q
    and positive_primitive(s.Poly(singular_value_factors[1][0], q)).as_expr()
    == Q69.as_expr(),
    "q-plane singular factor identity changed",
)
for start, end in zip(q_hull, q_hull[1:] + q_hull[:1]):
    face = face_polynomial(R_tq, start, end)
    require(face.degree() == int(s.gcd(abs(end[0] - start[0]), abs(end[1] - start[1]))),
            "q-plane edge degree changed")
    require(s.gcd(face, s.Poly(s.diff(face.as_expr()), *face.gens)).degree() == 0,
            "q-plane boundary face is not squarefree")

# Weighted adjunction for a (6,23) complete intersection in P(1,1,2,3).
# The curve avoids the two quotient singular points, so its dualizing degree is
# (6+23-7)*(6*23/6)=506 and p_a=254.  The same number follows from the
# weighted Hilbert-series count h_0(22)-h_0(16).
def weighted_monomial_count(degree: int) -> int:
    return sum(
        degree - 2 * iv - 3 * iz + 1
        for iz in range(degree // 3 + 1)
        for iv in range((degree - 3 * iz) // 2 + 1)
    )


arithmetic_genus = weighted_monomial_count(22) - weighted_monomial_count(16)
require(arithmetic_genus == 254, "weighted complete-intersection genus changed")
cusp_delta = (4 - 1) * (23 - 1) // 2
require(cusp_delta == 33, "(4,23) cusp delta changed")
generic_normalization_genus = arithmetic_genus - 5 * cusp_delta
require(generic_normalization_genus == 89, "generic normalization genus changed")

# Pick's theorem for conv{(0,3),(15,0),(15,23)}.
twice_q_area = 15 * 23
q_boundary_lattice_length = 3 + 23 + 5
q_newton_interior = (twice_q_area - q_boundary_lattice_length + 2) // 2
require(q_newton_interior == 158, "q-plane Newton interior count changed")
q_projection_delta = q_newton_interior - generic_normalization_genus
require(q_projection_delta == 69, "q-plane projection delta changed")

cover_degree = 23
total_ramification = 2 * generic_normalization_genus - 2 + 2 * cover_degree
fixed_ramification = 5 * (4 - 1)
infinity_ramification = 0
remaining_ramification = total_ramification - fixed_ramification
require(
    (total_ramification, fixed_ramification, remaining_ramification)
    == (222, 15, 207),
    "Riemann--Hurwitz accounting changed",
)

print("degree-22 split even-Faber prime-23 intrinsic geometry")
print(f"base_split_lambda_sha256={BASE_SHA256}")
print("scope=ALL_ONE_EXACT_GENUS89_PLUS_GENERIC_OPEN_NOT_UNIFORM_NOT_FULL_SPLIT")
print("weighted_closure=P(1,1,2,3)_complete_intersection_degrees_6_23")
print("all_one_degree6_surface=SMOOTH_EXACT_AFFINE_AND_INFINITY")
print("all_one_degree6_surface_irreducible=True")
print("ambient_weight_singular_points_avoided=True")
print("generic_complete_intersection=lci_Gorenstein_on_smooth_surface")
print("pencil_base_support=G3_plus_L5:no_h=t=0_base_point")
print("all_one_tv_eliminant=irreducible_over_Q:bidegree_46_23:terms_553")
print("tv_zeta_subresultant_degrees=5,2,1,0:generic_gcd_degree1:birational")
print("all_one_tv_newton_triangle=(0,0),(46,0),(0,23):interior_484")
print("all_one_geometric_irreducibility=prime23_plus_nontrivial_4cycle_inertia")
print("generic_infinity_resultant=irreducible_degree23:terms_1416:content_11^30")
print("generic_infinity_leading_coefficient=2^25*3^21*7^5*11^40:constant")
print("all_one_infinity=irreducible_mod17:23_transverse_points:ord_t=-1")
print("generic_infinity_ramification=0")
print("universal_t0_local_types=five_x^4=t^23_cusps_plus_three_smooth_points")
print("weighted_arithmetic_genus=254")
print("each_old_cusp_delta=(4-1)(23-1)/2=33")
print("bertini_generic_other_singularities=none")
print("generic_normalization_genus=254-5*33=89")
print("generic_RH=degree23:total_ramification222:t0_15:infinity_0:elsewhere_207")
print("fixed_plus_infinity_alone_do_not_force_positive_genus:genus0_budget_elsewhere29")
print("q_plane=birational_Q_irreducible:bidegree_15_23:terms_75")
print("q_plane_v_subresultant_degrees=3,2,1,0:linear_with_sole_q_content")
print("q_plane_newton_hull=(0,3),(15,0),(15,23):interior_158")
print("q_plane_all_three_boundary_faces=squarefree")
print(f"q_plane_Q69_sha256={Q69_SHA256}:squarefree=True")
print("q_plane_Q69_controls=two_finite_distinct_v_preimages:transverse_collision:nonzero_Hessian")
print("q_plane_singular_value_gcd=q^30*Q69^2:no_other_torus_singular_values")
print("all_one_q_plane_nodes=69_ordinary:normalization_genus_158-69=89")
print("q_plane_projection_delta_on_generic_genus89_locus=158-89=69")
print("Kummer_reformulation=globally_trivial_from_base_scout")
print("ALL CHECKS PASSED")
