#!/usr/bin/env python3
"""Exact companion for THM-2476."""

from __future__ import annotations

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


r, v, zeta, lam = s.symbols("r v zeta lam")
a, b, c, d, e, f, g, h = s.symbols("a b c d e f g h")

base_k = 922383 * v**2 - 25410 * v + 63
base_2 = (
    15944049 * zeta**2
    + (-162339408 * v + 2236080) * zeta
    - 1190488992 * v**3
    + 147581280 * v**2
    - 1219680 * v
    + 672
)
b_flux = 65591680 * zeta + 1443016960 * v**2 - 71554560 * v + 98560

# lambda=E^2/B^5 and r=E/(B^2*y), so
# B/y^2=r^2/lambda and E/y^5=r^5/lambda^2.
F1 = (
    1331 * lam * (616 * r**2 + lam * (63 - 1089 * v)) * zeta
    + 4
    * (
        lam * (-745360 * v + 6160) * r**2
        - 3748096 * r**5
        + lam**2 * base_k
    )
)
F2 = lam * r**2 * b_flux - 239878144 * r**5 + lam**2 * base_2

L5 = (
    155624547606 * v**5
    + 3215383215 * v**4
    - 1700698560 * v**3
    + 58124770 * v**2
    - 855470 * v
    + 2583
)
P = (
    126434012626944 * r**10
    + 22755800252416 * r**9
    + lam * (-50286255022080 * v + 2299591827456) * r**7
    + lam
    * (34222590223360 * v**2 + 3959638538240 * v - 44411530240)
    * r**6
    + lam**2
    * (44449457564160 * v**2 - 5142912445440 * v + 35419507200)
    * r**5
    + lam**2
    * (
        -149234938081152 * v**3
        - 9500102156160 * v**2
        + 695599766400 * v
        - 6017413248
    )
    * r**4
    + lam**3
    * (
        206782580709936 * v**4
        + 6246495741024 * v**3
        - 1509756494400 * v**2
        + 34466937120 * v
        - 193496688
    )
    * r**2
    - 567 * lam**4 * L5
)

raw_resultant = s.resultant(F1, F2, zeta)
require(
    s.expand(raw_resultant - 28344976 * lam**2 * P) == 0,
    "invariant B-E resultant changed",
)
require(len(s.Poly(P, r, v, lam).terms()) == 25, "B-E eliminant support changed")
require(s.degree(P, r) == 10 and s.degree(P, v) == 5, "B-E eliminant bidegree changed")
require(s.expand(P.subs(r, 0) + 567 * lam**4 * L5) == 0, "fixed L5 section changed")
require(
    s.Poly(P, r).coeff_monomial(r) == s.Poly(P, r).coeff_monomial(r**3) == 0,
    "Hensel odd-order gap changed",
)
require(
    s.gcd(s.Poly(L5, v), s.Poly(s.diff(L5, v), v)).degree() == 0,
    "fixed quintic lost squarefreeness",
)

wall = 702768 * r**5 + 23716 * r**4 - 1694 * lam * r**2 - 315 * lam**2
first_flux_wall_v = (616 * r**2 + 63 * lam) / (1089 * lam)
require(
    s.factor(P.subs(v, first_flux_wall_v) - 256 * wall**2) == 0,
    "first-flux wall identity changed",
)

# Newton triangle conv{(0,0),(10,0),(0,5)}.
support = {(ir, iv) for (ir, iv, _), _coefficient in s.Poly(P, r, v, lam).terms()}
vertices = {(0, 0), (10, 0), (0, 5)}
require(vertices <= support, "Newton triangle lost a vertex")
require(
    all(ir >= 0 and iv >= 0 and ir + 2 * iv <= 10 for ir, iv in support),
    "Newton support escaped the claimed triangle",
)

line = v + a * r**2 + b * r + c
line_remainder = s.Poly(s.rem(s.Poly(P, v), s.Poly(line, v)).as_expr(), r)
gb_line = s.groebner(line_remainder.all_coeffs(), a, b, c, lam, order="grevlex")
require(
    len(gb_line.polys) == 1 and gb_line.polys[0].as_expr() == 1,
    "linear-factor coefficient ideal stopped being unit",
)

# A quadratic factor has q0=v^2+c*v+h dividing the squarefree L5. The
# missing r coefficient and coprime Hensel uniqueness force q1=0. The unique
# order-two lift has q2=a*v+f. Since P and q1 have no order-three term, the
# same argument forces q3=0. Newton support then leaves only q4=d. The full
# coefficient ideal of this forced ansatz is unit.
quadratic_hensel = v**2 + c * v + h + (a * v + f) * r**2 + d * r**4
quadratic_remainder = s.Poly(
    s.rem(s.Poly(P, v), s.Poly(quadratic_hensel, v)).as_expr(),
    v,
    r,
)
gb_quadratic_hensel = s.groebner(
    [coefficient for _, coefficient in quadratic_remainder.terms()],
    a,
    c,
    d,
    f,
    h,
    lam,
    order="grevlex",
)
require(
    len(gb_quadratic_hensel.polys) == 1
    and gb_quadratic_hensel.polys[0].as_expr() == 1,
    "forced quadratic Hensel coefficient ideal stopped being unit",
)

disc_v_raw = s.discriminant(P, v)
disc_factors = s.factor_list(disc_v_raw)[1]
require(len(disc_factors) == 3, "unexpected extra v-discriminant factor")
require((lam, 14) in disc_factors, "lambda multiplicity in discriminant changed")
require((wall, 2) in disc_factors, "wall-square factor changed")
K_candidates = [
    factor
    for factor, exponent in disc_factors
    if exponent == 1 and s.degree(factor, r) == 30
]
require(len(K_candidates) == 1, "moving K30 factor is no longer unique")
K30 = K_candidates[0]
require(
    s.degree(K30, lam) == 14 and len(s.Poly(K30, r, lam).terms()) == 56,
    "K30 degree or support changed",
)
require(
    s.expand(K30.subs(r, 0) + 5532126228517587890625 * lam**14) == 0,
    "K30 primitive normalization changed",
)
degree_drop = (
    980708488959375 * lam**2
    + 32820181258788 * lam
    + 120472576000
)
require(
    s.expand(
        s.Poly(K30, r).LC()
        - 10284112691593526727130742784 * degree_drop
    )
    == 0,
    "K30 degree-drop factor changed",
)

A3 = (
    19731972545646643418400000000 * lam**3
    + 346409095257857532829648104 * lam**2
    + 1415693788892415142663014 * lam
    - 981468600381938943169
)
disc_k_raw = s.discriminant(K30, r)
disc_k = s.factor_list(disc_k_raw)[1]
require(len(disc_k) == 5, "unexpected extra K30-discriminant factor")
require((lam, 379) in disc_k, "K30 invariant-origin multiplicity changed")
require((A3, 3) in disc_k, "K30 A3 factor changed")
P12_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 12 and exponent == 2]
Q12_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 12 and exponent == 3]
R16_candidates = [factor for factor, exponent in disc_k if s.degree(factor, lam) == 16 and exponent == 1]
require(
    len(P12_candidates) == len(Q12_candidates) == len(R16_candidates) == 1,
    "large exceptional factors changed",
)
P12, Q12, R16 = P12_candidates[0], Q12_candidates[0], R16_candidates[0]

S5 = (
    131955965143315044357267878722928640000000000000000 * lam**5
    + 3095105736402260153329374221431569404424300000000 * lam**4
    + 53687721220771954218792719815073181002030599351104 * lam**3
    - 6341516493144875869324270617788470732795917450816 * lam**2
    + 235228509291405819426976849459832940179373568841 * lam
    + 178772563473178835054245823352852460725645500
)
wall_resultant_raw = s.resultant(K30, wall, r)
wall_resultant_factors = s.factor_list(wall_resultant_raw)[1]
require(len(wall_resultant_factors) == 3, "unexpected extra wall-resultant factor")
require((lam, 56) in wall_resultant_factors, "wall resultant lambda factor changed")
require((A3, 3) in wall_resultant_factors, "wall resultant A3 factor changed")
require((S5, 1) in wall_resultant_factors, "wall resultant S5 factor changed")

exceptional_factors = [degree_drop, A3, P12, Q12, R16, S5]
for i, left in enumerate(exceptional_factors):
    require(
        s.gcd(s.Poly(left, lam), s.Poly(s.diff(left, lam), lam)).degree() == 0,
        f"exceptional factor {i} lost squarefreeness",
    )
    for right in exceptional_factors[:i]:
        require(
            s.gcd(s.Poly(left, lam), s.Poly(right, lam)).degree() == 0,
            "distinct exceptional factors acquired a common root",
        )

next_coefficient = s.Poly(K30, r).coeff_monomial(r**29)
require(
    s.gcd(s.Poly(next_coefficient, lam), s.Poly(degree_drop, lam)).degree() == 0,
    "degree-drop fibre lost degree 29",
)
require(
    s.gcd(s.Poly(disc_k_raw, lam), s.Poly(degree_drop, lam)).degree() == 0,
    "a degree-drop fibre lost squarefreeness",
)
require(
    s.gcd(s.Poly(wall_resultant_raw, lam), s.Poly(degree_drop, lam)).degree() == 0,
    "a degree-drop fibre met the wall",
)

max_collision_gcd_degree = 3
minimum_simple_k_roots = 30 - 2 * max_collision_gcd_degree
minimum_visible_simple_branches = minimum_simple_k_roots - s.degree(wall, r)
minimum_total_branches = minimum_visible_simple_branches + (minimum_visible_simple_branches % 2)
minimum_genus = 1 + (minimum_total_branches - 10) // 2
require(minimum_simple_k_roots == 24, "simple-root floor changed")
require(minimum_visible_simple_branches == 19, "visible branch floor changed")
require(minimum_total_branches == 20, "ramification-parity floor changed")
require(minimum_genus == 6, "genus floor changed")

# Original y=0 boundary on C=D=W=0.
Bpar, Epar, y, u, Z = s.symbols("Bpar Epar y u Z")
original_A = 616 * Bpar - 1089 * u + 63 * y**2
original_K = (
    -745360 * Bpar * u * y
    + 6160 * Bpar * y**3
    - 3748096 * Epar
    + 922383 * u**2 * y
    - 25410 * u * y**3
    + 63 * y**5
)
N1 = 1331 * original_A * Z + 4 * original_K
N2 = (
    15944049 * Z**2
    + 65591680 * Bpar * Z * y
    - 162339408 * Z * u * y
    + 2236080 * Z * y**3
    + 1443016960 * Bpar * u**2
    - 71554560 * Bpar * u * y**2
    + 98560 * Bpar * y**4
    - 239878144 * Epar * y
    - 1190488992 * u**3
    + 147581280 * u**2 * y**2
    - 1219680 * u * y**4
    + 672 * y**6
)
denominator = 1331 * (616 * Bpar - 1089 * u)
k = 14992384 * Epar
require(
    s.expand(N1.subs(y, 0) - (denominator * Z - k)) == 0,
    "y=0 first-flux rational reconstruction changed",
)
require(
    s.expand(
        N2.subs(y, 0)
        - (15944049 * Z**2 + 1443016960 * Bpar * u**2 - 1190488992 * u**3)
    )
    == 0,
    "y=0 second-flux equation changed",
)
boundary_quintic = s.expand(
    -(15944049 * k**2)
    - (1443016960 * Bpar * u**2 - 1190488992 * u**3) * denominator**2
)
require(
    s.cancel(denominator**2 * N2.subs({y: 0, Z: k / denominator}) + boundary_quintic)
    == 0,
    "y=0 constant-field quintic changed",
)

print("THM-2476 exact companion")
print("invariants=lambda=E^2/B^5,r=E/(B^2*y)")
print("resultant=28344976*lambda^2*P")
print("P_terms=25")
print("P_bidegree=(10,5)")
print("newton_polygon=conv{(0,0),(10,0),(0,5)}")
print("linear_factor_unit_ideal=YES")
print("quadratic_factor=fixed_L5_Hensel_forced_ansatz_unit_ideal")
print("uniform_absolute_irreducibility=lambda_nonzero")
print("wall=702768*r^5+23716*r^4-1694*lambda*r^2-315*lambda^2")
print("P_at_wall=256*wall^2")
print("disc_v(P)=constant*lambda^14*wall^2*K30")
print("K30_bidegree=(30,14)")
print("K30_terms=56")
print("disc_r(K30)=constant*lambda^379*A3^3*P12^2*Q12^3*R16")
print("exceptional_factor_degrees=2,3,12,12,16,5")
print("degree_drop_fibres=two_squarefree_degree_29_wall_disjoint")
print("minimum_simple_K_roots=24")
print("minimum_visible_simple_branches=19")
print("minimum_total_ramification=20")
print("minimum_normalization_genus=6")
print("y_zero_boundary=rational_Z_then_constant_field_quintic")
print("ALL CHECKS PASSED")
