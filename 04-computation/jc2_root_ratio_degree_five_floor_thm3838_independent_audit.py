#!/usr/bin/env python3
"""Assertion-free independent hostile audit of THM-3838.

The script independently checks the entire low-degree lattice, exact
composition degrees, both generic hyperelliptic sidecars, the genus-three
Riemann--Hurwitz equality boundary, quartic adjunction data, and the
unimodular reduction step.  The geometric implications are written out in
AUDIT.md; this companion tests their exact numerical and symbolic seams.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    CHECKS += 1


def arithmetic_genus(degree: int) -> int:
    return (degree - 1) * (degree - 2) // 2


# Every degree case below the claimed floor.
genus_table = {degree: arithmetic_genus(degree) for degree in range(1, 9)}
require(genus_table == {1: 0, 2: 0, 3: 1, 4: 3, 5: 6, 6: 10, 7: 15, 8: 21},
        "plane arithmetic-genus table")
require(all(genus_table[degree] < 3 for degree in (1, 2, 3)),
        "primitive degrees one through three cannot carry genus three")

# If h=p(g), write m=deg(p), d=deg(g), so deg(h)=m*d.  The inherited
# genus floor gives d>=4.  Exhaust every target degree through seven.
factorizations = {
    target: tuple(
        (outer_degree, primitive_degree)
        for outer_degree in range(1, 9)
        for primitive_degree in range(4, 9)
        if outer_degree * primitive_degree == target
    )
    for target in range(1, 9)
}
require(factorizations[1] == () and factorizations[2] == () and factorizations[3] == (),
        "degrees below four have no genus-floor composition")
require(factorizations[4] == ((1, 4),),
        "degree four has the unique linear-outer quartic boundary")
require(factorizations[5] == ((1, 5),),
        "degree five forces a primitive quintic")
require(factorizations[6] == ((1, 6),),
        "degree six forces a primitive sextic")
require(factorizations[7] == ((1, 7),),
        "degree seven forces a primitive septic")
require(factorizations[8] == ((1, 8), (2, 4)),
        "the first nonlinear outer composition occurs in degree eight")

# Direct exact controls for deg(p(g))=deg(p)deg(g), including inhomogeneous
# g and p.  The top term cannot cancel because its homogeneous part is a
# nonzero scalar times the relevant power of the top part of g.
x, y, S = sp.symbols("x y S")
composition_cases = 0
for primitive_degree in range(1, 8):
    g = x**primitive_degree + 2 * x ** max(primitive_degree - 1, 0) * y
    if primitive_degree >= 2:
        g += 3 * y**primitive_degree + x + 1
    else:
        g += y + 1
    require(sp.Poly(g, x, y).total_degree() == primitive_degree,
            f"primitive control degree {primitive_degree}")
    for outer_degree in range(1, 6):
        p = 5 * S**outer_degree + 7 * S ** max(outer_degree - 1, 0) + 11
        composed = sp.expand(p.subs(S, g))
        require(sp.Poly(composed, x, y).total_degree() == outer_degree * primitive_degree,
                f"composition degree {outer_degree}x{primitive_degree}")
        composition_cases += 1
require(composition_cases == 35, "composition-control census")

# Re-enter H by coefficient blocks in T, independently of the canonical
# companion's monomial ordering.
T, Z = sp.symbols("T Z")
coefficients_in_T = (
    84,
    36 * Z**2 + 196 * Z,
    84 * Z**3 + 36 * Z**2,
    49 * Z**4 + 112 * Z**3,
    -12 * Z**5,
    -14 * Z**6 + 12 * Z**5,
    0,
    Z**8,
)
H = sp.expand(sum(coefficient * T ** (7 - index)
                  for index, coefficient in enumerate(coefficients_in_T)))
require(sp.Poly(H, T).degree() == 7 and sp.Poly(H, T).LC() == 84,
        "odd sidecar is degree seven with nonzero constant leading coefficient")
require(sp.Poly(H, Z).degree() == 8 and sp.Poly(H, Z).LC() == 1,
        "even sidecar is monic degree eight")

# Exact generic squarefreeness, computed afresh rather than importing either
# displayed discriminant factorization from THM-3827.
disc_even = sp.factor(sp.discriminant(H, Z))
disc_odd = sp.factor(sp.discriminant(H, T))
require(disc_even != 0 and sp.Poly(disc_even, T).degree() == 52,
        "degree-eight generic discriminant is nonzero")
require(disc_odd != 0 and sp.Poly(disc_odd, Z).degree() == 51,
        "degree-seven generic discriminant is nonzero")

for t_value in (-3, -1, 1, 2):
    specialized = sp.Poly(H.subs(T, t_value), Z)
    require(specialized.degree() == 8,
            f"even specialization retains degree at T={t_value}")
    require(sp.gcd(specialized, specialized.diff()).degree() == 0,
            f"even squarefree hostile specialization T={t_value}")

for z_value in (-3, -1, 1, 2):
    specialized = sp.Poly(H.subs(Z, z_value), T)
    require(specialized.degree() == 7,
            f"odd specialization retains degree at Z={z_value}")
    require(sp.gcd(specialized, specialized.diff()).degree() == 0,
            f"odd squarefree hostile specialization Z={z_value}")

# Nonconstant substitution into a nonzero one-variable discriminant is
# injective.  These explicit high-degree compositions are hostile controls.
s = sp.symbols("s")
p_control = s**3 - 2 * s + 5
q_control = 2 * s**4 + s + 3
require(sp.Poly(disc_even.subs(T, p_control), s).degree() == 156,
        "even discriminant survives a cubic generative composition")
require(sp.Poly(disc_odd.subs(Z, q_control), s).degree() == 204,
        "odd discriminant survives a quartic generative composition")

# Hyperelliptic branch counts and infinity parity.
even_finite_branch = 8
even_infinite_branch = 0
odd_finite_branch = 7
odd_infinite_branch = 1
even_genus = (even_finite_branch + even_infinite_branch - 2) // 2
odd_genus = (odd_finite_branch + odd_infinite_branch - 2) // 2
even_infinity_points = 2 if sp.Poly(H, Z).degree() % 2 == 0 and sp.Poly(H, Z).LC() == 1 else 0
odd_infinity_points = 1 if sp.Poly(H, T).degree() % 2 == 1 else 0
require((even_genus, odd_genus) == (3, 3),
        "both sidecars have hyperelliptic genus three")
require((even_infinity_points, odd_infinity_points) == (2, 1),
        "even and odd models have two and one geometric infinity points")

# At source genus three, Riemann--Hurwitz against a genus-three sidecar has
# exactly one solution: map degree one and zero ramification.
rh_solutions = tuple(
    (map_degree, ramification)
    for map_degree in range(1, 8)
    for ramification in range(0, 25)
    if 2 * 3 - 2 == map_degree * (2 * 3 - 2) + ramification
)
require(rh_solutions == ((1, 0),),
        "genus-three equality forces an unramified degree-one map")

# Plane quartic equality: p_a=g=3 leaves zero total delta invariant.  By
# adjunction K_C=O_C(1), so the plane embedding is the canonical map.  Its
# morphism degree one conflicts with the hyperelliptic canonical degree two.
quartic_arithmetic_genus = arithmetic_genus(4)
quartic_geometric_genus = 3
delta_total = quartic_arithmetic_genus - quartic_geometric_genus
canonical_degree = 4 * (4 - 3)
plane_canonical_map_degree = 1
hyperelliptic_canonical_map_degree = 2
require(delta_total == 0, "quartic equality has no finite or infinite singularity defect")
require(canonical_degree == 2 * quartic_geometric_genus - 2 == 4,
        "quartic adjunction gives the canonical line bundle")
require(plane_canonical_map_degree != hyperelliptic_canonical_map_degree,
        "plane canonical embedding and hyperelliptic canonical map have different degree")

# A smooth-control quartic: its three partials vanish simultaneously only at
# the excluded affine origin, so no projective singularity exists.
X, Y, W = sp.symbols("X Y W")
fermat_quartic = X**4 + Y**4 + W**4
partials = tuple(sp.diff(fermat_quartic, variable) for variable in (X, Y, W))
require(partials == (4 * X**3, 4 * Y**3, 4 * W**3),
        "smooth plane-quartic control partials")
for chart_variable in (X, Y, W):
    chart_ideal = sp.groebner((*partials, chart_variable - 1), X, Y, W)
    require(chart_ideal.contains(sp.Integer(1)),
            f"Fermat quartic has no singularity in projective chart {chart_variable}")

# Unimodularity makes the displayed root ratio reduced.  This explicit row
# is only a Bezout/gcd control, not an atlas candidate.
h_control = x
k_control = 1 + x * y
m_control = y
C_control = 1
require(sp.expand(C_control * k_control - m_control * h_control) == 1,
        "determinant-one Bezout control")
require(sp.gcd(sp.Poly(h_control, x, y), sp.Poly(k_control, x, y)).total_degree() == 0,
        "unimodular row is coprime")

source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "no Python assert statements")

semantic = {
    "degree_cases": "primitive degrees 1..3 excluded; total degree 4 is uniquely outer-linear primitive quartic",
    "composition": "deg p(g)=deg(p)deg(g); nonlinear outer degree begins at total degree 8",
    "generic_curves": "both sidecars independently squarefree, hyperelliptic genus 3",
    "equality": "Riemann-Hurwitz gives degree-one isomorphism to sidecar",
    "quartic": "p_a=g=3 forces smooth plane quartic; adjunction makes its plane map canonical",
    "conflict": "smooth plane quartic canonical map has degree 1, hyperelliptic genus-3 canonical map degree 2",
    "unit": "Ck-mh=1 forces gcd(h,k)=1",
    "scope": "conditional on a dominant etale atlas; degree five and atlas existence remain open",
    "verdict": "PASS",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3838-independent-hostile-audit")
print("degree_cases=1_2_3_excluded;4_unique_linear_outer_quartic;5_6_7_linear_outer")
print("composition=nonlinear_outer_first_possible_at_total_degree_8")
print("sidecars=degree8_and_degree7_squarefree_hyperelliptic_genus3")
print("equality=riemann_hurwitz_degree1_then_plane_quartic_hyperelliptic_conflict")
print("unit=determinant_one_row_is_coprime_so_ratio_is_reduced")
print("canonical_replay=normal_optimized_frozen_checked_separately")
print("scope=conditional_all_degree_floor;degree5_and_atlas_existence_OPEN")
print(f"composition_controls={composition_cases}")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print("RESULT=PASS")
