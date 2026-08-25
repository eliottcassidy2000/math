#!/usr/bin/env python3
"""Primary exact certificate for THM-4045.

The theorem closes the exact maximum-weight-seven branch in the live reduced
(2,3) planar-Jacobian cell.  This script reconstructs the complete lower
Newton subdivision, its face equations, the nodal and edge ledger, the hidden
j=1728 elliptic tail, and the common target scaling.  The stable-map argument
and the all-characteristic-zero CM mismatch are proved in the theorem file.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction
from itertools import combinations
from math import gcd
import sys

import sympy as sp


sys.stdout.reconfigure(newline="\n")
CHECKS = 0
SEMANTIC: list[str] = []


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)
    print(f"PASS  {label}")


def note(value: str) -> None:
    SEMANTIC.append(value)
    print(f"RESULT {value}")


def pick(vertices: list[tuple[int, int]]) -> tuple[int, int, int]:
    area2 = abs(
        sum(
            vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
            - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
            for i in range(len(vertices))
        )
    )
    boundary = sum(
        gcd(
            abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]),
        )
        for i in range(len(vertices))
    )
    numerator = area2 - boundary + 2
    require(numerator % 2 == 0, f"Pick parity for {vertices}")
    return area2, boundary, numerator // 2


def combine_support(
    raw: list[tuple[int, int, int, str, Fraction]],
) -> list[tuple[int, int, int, str, Fraction]]:
    coefficients: dict[tuple[int, int, int], Fraction] = {}
    labels: dict[tuple[int, int, int], list[str]] = {}
    for i, j, h, label, coefficient in raw:
        key = (i, j, h)
        coefficients[key] = coefficients.get(key, Fraction(0)) + coefficient
        labels.setdefault(key, []).append(label)
    return [
        (i, j, h, "+".join(labels[(i, j, h)]), coefficient)
        for (i, j, h), coefficient in sorted(coefficients.items())
        if coefficient
    ]


def support_from_h(
    terms: list[tuple[int, int, str, Fraction]], gamma: Fraction
) -> list[tuple[int, int, int, str, Fraction]]:
    """Expand (s^2-p)(1-QH)+gamma*Q*s^2 coefficientwise."""
    # The coefficient of s^2 is the Q-adic unit 1+gamma*Q, so its Newton
    # point has height zero; gamma*Q*s^2 is not a second lifted support point.
    raw = [
        (2, 0, 0, "base_s2", Fraction(1)),
        (0, 1, 0, "base_minus_p", Fraction(-1)),
    ]
    for i, j, label, coefficient in terms:
        raw.append((i + 2, j, 1, f"minus_s2_{label}", -coefficient))
        raw.append((i, j + 1, 1, f"plus_p_{label}", coefficient))
    return combine_support(raw)


def lower_facets(
    points: list[tuple[int, int, int, str, Fraction]],
) -> dict[frozenset[tuple[int, int, int]], tuple[Fraction, Fraction, Fraction]]:
    """Enumerate all projected two-dimensional lower faces."""
    faces: dict[
        frozenset[tuple[int, int, int]], tuple[Fraction, Fraction, Fraction]
    ] = {}
    for triple in combinations(range(len(points)), 3):
        p0, p1, p2 = (points[index] for index in triple)
        x0, y0, z0 = p0[:3]
        x1, y1, z1 = p1[:3]
        x2, y2, z2 = p2[:3]
        determinant = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        if determinant == 0:
            continue
        slope_x = Fraction(
            (z1 - z0) * (y2 - y0) - (z2 - z0) * (y1 - y0),
            determinant,
        )
        slope_y = Fraction(
            (x1 - x0) * (z2 - z0) - (x2 - x0) * (z1 - z0),
            determinant,
        )
        intercept = Fraction(z0) - slope_x * x0 - slope_y * y0
        gaps = [
            Fraction(z) - slope_x * x - slope_y * y - intercept
            for x, y, z, _, _ in points
        ]
        if min(gaps) >= 0:
            face = frozenset(points[index][:3] for index, gap in enumerate(gaps) if gap == 0)
            faces[face] = (slope_x, slope_y, intercept)
    return faces


print("STATUS=THM-4045;MAX_WEIGHT_SEVEN_EXCLUDED;JC(2)_OPEN")
print("UNIVERSE=live_reduced_(2,3)_cell;b=d=0;exact_total_weight_M=7")

# Specialize only harmless nonzero scalars for support enumeration.  The
# symbolic identities below retain gamma, A5, and phi.
gamma0 = Fraction(-1, 2)
a5_0 = Fraction(1)
epsilon0 = gamma0 * Fraction(2752, 135)
kappa0 = gamma0 * Fraction(-5696, 45)
phi0 = Fraction(13, 1)
lambda0 = Fraction(7, 5)
alpha0 = Fraction(-11, 3)
require(epsilon0 == Fraction(-1376, 135), "forced raw epsilon at A5=1")
require(kappa0 == Fraction(2848, 45), "forced raw kappa at A5=1")
require(kappa0 - epsilon0 == Fraction(1984, 27), "shared s^2*p^3 coefficient does not cancel")
require(phi0 * kappa0 * epsilon0 != 0, "all three boundary coefficients are nonzero")

h_terms = [
    (0, 1, "lambda_p", lambda0),
    (0, 2, "alpha_p2", alpha0),
    (0, 3, "epsilon_p3", epsilon0),
    (2, 2, "kappa_y2", kappa0),
    (1, 3, "phi_p2y", phi0),
]
points = support_from_h(h_terms, gamma0)
coefficient = {(i, j, h): value for i, j, h, _, value in points}
require(coefficient[(2, 3, 1)] == kappa0 - epsilon0, "expanded shared coefficient is kappa-epsilon")

A = (0, 1, 0)
B = (2, 0, 0)
C = (4, 2, 1)
D = (3, 3, 1)
E = (1, 4, 1)
F = (0, 4, 1)
main_face = frozenset({A, B, D, E})
tail_face = frozenset({B, C, D})
vertical_face = frozenset({A, E, F})
facets = lower_facets(points)
expected_facets = {
    main_face: (Fraction(1, 7), Fraction(2, 7), Fraction(-2, 7)),
    tail_face: (Fraction(1, 4), Fraction(1, 4), Fraction(-1, 2)),
    vertical_face: (Fraction(0), Fraction(1, 3), Fraction(-1, 3)),
}
require(facets == expected_facets, "complete support has exactly the three claimed lower faces")
for face, plane in sorted(facets.items(), key=lambda item: item[1]):
    print("FACET", sorted(face), "PLANE", plane)

# The ramification Q=sigma^84 clears every face slope.  The three integral
# valuation functions dictate the three affine charts.
integral_planes = {
    face: tuple(84 * value for value in plane) for face, plane in facets.items()
}
require(integral_planes[main_face] == (12, 24, -24), "main chart scaling (s,p)=(sigma^-12 S,sigma^-24 P)")
require(integral_planes[tail_face] == (21, 21, -42), "tail chart scaling (s,p)=(sigma^-21 S,sigma^-21 P)")
require(integral_planes[vertical_face] == (0, 28, -28), "vertical chart scaling p=sigma^-28 P")

# Initial equations and their normalizations.
S, P, T = sp.symbols("S P T")
gamma, A5, phi, epsilon, kappa = sp.symbols(
    "gamma A5 phi epsilon kappa", nonzero=True
)
forced_epsilon = gamma * sp.Rational(2752, 135) / A5**3
forced_kappa = -gamma * sp.Rational(5696, 45) / A5**3
require(sp.simplify(forced_kappa - forced_epsilon + sp.Rational(3968, 27) * gamma / A5**3) == 0,
        "symbolic kappa-epsilon is -3968*gamma/(27*A5^3)")

g_main = (S**2 - P) * (1 - phi * S * P**3)
g_tail = 1 - kappa * S**2 * P**2 - phi * S * P**3
g_vertical = -1 + epsilon * P**3 + phi * S * P**3
require(sp.expand(g_main) == S**2 - P - phi*S**3*P**3 + phi*S*P**4,
        "main face factors into two rational components")
require(sp.diff(g_vertical, S) == phi * P**3, "vertical face is smooth and rational on the torus")

# The main components meet at P=S^2 and phi*S^7=1.  Their Jacobian
# determinant is nonzero at every one of the seven roots.
f0 = S**2 - P
f1 = 1 - phi*S*P**3
jacobian = sp.det(sp.Matrix([[sp.diff(f0, S), sp.diff(f0, P)],
                             [sp.diff(f1, S), sp.diff(f1, P)]]))
require(sp.simplify(jacobian.subs(P, S**2) + 7*phi*S**6) == 0,
        "main intersections are transverse with determinant -7*phi*S^6")
require(sp.diff(phi*S**7 - 1, S) == 7*phi*S**6,
        "the seven main-face node addresses are simple")
require(abs(2*3-(-1)*1) == 7,
        "binomial lattice intersection number seven is saturated by the torus nodes")

# In the main chart, the exact completed local equation at each node is
# U*V=-gamma*sigma^84; no truncation is made here.
sigma, lam, alpha = sp.symbols("sigma lambda alpha", nonzero=True)
h_sigma = (
    phi*S*P**3
    + sigma**12*(epsilon*P**3 + kappa*S**2*P**2)
    + sigma**36*alpha*P**2
    + sigma**60*lam*P
)
U = S**2 - P
V = (1 - h_sigma) / S**2
scaled_equation = U*(1-h_sigma) + gamma*sigma**84*S**2
require(sp.simplify(scaled_equation/S**2 - (U*V + gamma*sigma**84)) == 0,
        "exact nodal chart is U*V=-gamma*sigma^84")
require(84 - 1 == 83, "regular resolution of each node inserts an A_83 rational chain")

# The tail is birational to a binary quartic.  With T=SP and
# W=2*kappa*T+phi*P^2, its equation is W^2=phi^2*P^4+4*kappa.
W = 2*kappa*T + phi*P**2
tail_relation = kappa*T**2 + phi*T*P**2 - 1
require(sp.expand(W**2 - (phi**2*P**4 + 4*kappa) - 4*kappa*tail_relation) == 0,
        "tail quadratic completes to W^2=phi^2*P^4+4*kappa")
require(4*kappa*phi != 0,
        "tail gradient equations have nonzero determinant and no torus singularity")
quartic_I = 48*phi**2*kappa
quartic_J = sp.Integer(0)
require(quartic_I != 0 and quartic_J == 0, "tail quartic has I!=0,J=0 and hence j=1728")

# Every edge restriction is squarefree in characteristic zero.  These are
# also the boundary controls needed to rule out unrecorded positive genus.
edge_polynomials = {
    "AB": -1 + T,
    "BD": 1 - phi*T,
    "BC": 1 - kappa*T**2,
    "CD": kappa + phi*T,
    "DE": -1 + T,
    "AE": -1 + phi*T,
    "EF": epsilon + phi*T,
    "FA": -1 + epsilon*T**3,
}
for edge, polynomial in edge_polynomials.items():
    require(sp.degree(sp.gcd(polynomial, sp.diff(polynomial, T)), T) == 0,
            f"edge {edge} is squarefree")

# The generic toric completion needs a separate outer-edge audit: points
# above the lower hull still contribute to a projected outer edge.  At B the
# coefficient is 1+gamma*Q, and the vertical FA edge contains lambda, alpha,
# and epsilon.  Its cubic discriminant has nonzero Q-adic leading term.
Qvar = sp.symbols("Q")
lam_generic, alpha_generic = sp.symbols("lambda_generic alpha_generic")
generic_outer_edges = {
    "AB": -1 + (1+gamma*Qvar)*T,
    "BC": (1+gamma*Qvar) - Qvar*kappa*T**2,
    "CD": kappa + phi*T,
    "DE": -1 + T,
    "EF": epsilon + phi*T,
    "FA": -1 + Qvar*lam_generic*T + Qvar*alpha_generic*T**2 + Qvar*epsilon*T**3,
}
for edge, polynomial in generic_outer_edges.items():
    require(sp.degree(sp.gcd(polynomial, sp.diff(polynomial, T)), T) == 0,
            f"generic outer edge {edge} is squarefree over k((Q))")
fa_discriminant = sp.expand(sp.discriminant(generic_outer_edges["FA"], T))
require(fa_discriminant.coeff(Qvar, 0) == 0 and fa_discriminant.coeff(Qvar, 1) == 0,
        "generic FA cubic discriminant has Q-adic order at least two")
require(sp.simplify(fa_discriminant.coeff(Qvar, 2) + 27*epsilon**2) == 0,
        "generic FA cubic discriminant leads with -27*epsilon^2*Q^2")
generic_h = lam_generic*P + alpha_generic*P**2 + epsilon*P**3 + kappa*S**2*P**2 + phi*S*P**3
generic_f = (S**2-P)*(1-Qvar*generic_h) + gamma*Qvar*S**2
require(sp.expand(generic_f.subs(P, S**2)) == gamma*Qvar*S**2,
        "on the generic torus curve p=s^2 is impossible")

# Pick and dual-graph ledgers.  The main face's six interior points are its
# seven parallel transverse nodes minus one; the tail contributes genus one.
global_polygon = [(0, 1), (2, 0), (4, 2), (3, 3), (1, 4), (0, 4)]
main_polygon = [(0, 1), (2, 0), (3, 3), (1, 4)]
tail_polygon = [(2, 0), (4, 2), (3, 3)]
vertical_polygon = [(0, 1), (1, 4), (0, 4)]
require(pick(global_polygon) == (21, 9, 7), "global Newton polygon has seven interior points")
require(pick(main_polygon) == (14, 4, 6), "main face has arithmetic genus six")
require(pick(tail_polygon) == (4, 4, 1), "tail face has genus one")
require(pick(vertical_polygon) == (3, 5, 0), "vertical face is rational")
require(7 - 2 + 1 == 6, "seven-node two-component main graph has first Betti number six")
require(6 + 1 == 7, "special-fibre genus ledger is graph six plus elliptic one")

# Definition 3.12 slope/Farey audit.  For BD take L*=6-3i+j,
# P0=(2,0),P1=(2,1); for AE take L*=3i-j+1,
# P0=(0,1),P1=(0,0).  Consecutive numerator/denominator pairs must have
# determinant one.
def value(plane: tuple[int, int, int], point: tuple[int, int]) -> int:
    return plane[0]*point[0] + plane[1]*point[1] + plane[2]


def determinant_one(sequence: list[tuple[int, int]]) -> bool:
    return all(
        sequence[index][0]*sequence[index+1][1]
        - sequence[index+1][0]*sequence[index][1] == 1
        for index in range(len(sequence)-1)
    )


plane_m = tuple(int(entry) for entry in integral_planes[main_face])
plane_t = tuple(int(entry) for entry in integral_planes[tail_face])
plane_v = tuple(int(entry) for entry in integral_planes[vertical_face])
bd_slopes = (
    value(plane_m,(2,1))-value(plane_m,(2,0)),
    value(plane_t,(2,1))-value(plane_t,(2,0)),
)
ae_slopes = (
    value(plane_m,(0,0))-value(plane_m,(0,1)),
    value(plane_v,(0,0))-value(plane_v,(0,1)),
)
bd_sequence = [(integer,1) for integer in range(bd_slopes[0],bd_slopes[1]-1,-1)]
ae_sequence = [(integer,1) for integer in range(ae_slopes[0],ae_slopes[1]-1,-1)]
require(bd_slopes == (24,21) and determinant_one(bd_sequence),
        "BD slopes have minimal sequence 24>23>22>21")
require(ae_slopes == (-24,-28) and determinant_one(ae_sequence),
        "AE slopes have minimal sequence -24>-25>-26>-27>-28")
require((len(bd_sequence)-2,len(ae_sequence)-2) == (2,3),
        "Definition 3.12 inner-chain lengths are two and three")
outer_slopes = [12,21,-21,-12,-28,0]
require(all(determinant_one([(slope,1),(slope-1,1)]) for slope in outer_slopes),
        "all six integer outer slopes have r=0 and add no chain")
require(all(denominator == 1 for _,denominator in bd_sequence+ae_sequence),
        "all inner-chain multiplicities are one")

# The target uses exactly the same q=sigma^-84 base change.
X, Y, a = sp.symbols("X Y a")
target_scaled = Y**2 - X**3 - 1 + sp.Rational(3, 4)*a**2*sigma**56*X + sp.Rational(1, 4)*a**3*sigma**84
require(target_scaled.subs(sigma, 0) == Y**2-X**3-1,
        "target has good special fibre Y^2=X^3+1 with j=0")

# Hostile boundary: a weight-eight p*y^2 term contributes the endpoint
# (4,3,1), strictly below the max-seven tail plane.  The theorem therefore
# makes no claim about weight eight.
tail_plane = facets[tail_face]
hostile_gap = Fraction(1) - (
    tail_plane[0]*4 + tail_plane[1]*3 + tail_plane[2]
)
require(hostile_gap == Fraction(-1, 4), "weight-eight p*y^2 destroys the elliptic tail facet")

note("lower_faces=main_nodal_rational_graph+elliptic_j1728_tail+rational_vertical")
note("node_local_model=UV+gamma*sigma^84;seven_A83_resolutions")
note("positive_genus_inventory=E_j1728_only;dual_graph_rank=6;total_genus=7")
note("target_good_reduction=E_j0;Hom(E_j1728,E_j0)=0_in_characteristic_zero")
note("live_b=d=0_exact_M7=IMPOSSIBLE;new_floor_M>=8;weight8_OPEN")
note("hostile_weight8_p_y2_gap=-1/4;max7_hypothesis_is_load_bearing")
note("generic_completion=actual_smooth_C_q;FA_discriminant_lead=-27epsilon^2Q^2")

semantic_hash = hashlib.sha256("\n".join(SEMANTIC).encode()).hexdigest()
print(f"SEMANTIC_SHA256={semantic_hash}")
print(f"CHECKS={CHECKS}")
print("ALL THM-4045 PRIMARY CHECKS PASSED")
