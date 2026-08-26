#!/usr/bin/env python3
"""Exact companion for THM-4146, extending THM-4139.

The certificate has four independent layers:

* the complete rational preperiodic graph of x^2-29/16;
* the universal central-sign lift attached to a quadratic three-cycle;
* the exact 3:4:5 and Delta-V horizontal-section intersections; and
* the exact-period-six dynatomic and Mersenne/Zsigmondy controls.

No Python ``assert`` is used, so optimized mode retains every truth gate.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def matrix_mod(matrix: sp.Matrix, modulus: int) -> sp.Matrix:
    return matrix.applyfunc(lambda entry: int(entry) % modulus)


def poly_pow_mod(
    base: sp.Poly, exponent: int, modulus_poly: sp.Poly, prime: int
) -> sp.Poly:
    result = sp.Poly(1, base.gens[0], modulus=prime)
    power = base.rem(modulus_poly)
    value = exponent
    while value:
        if value & 1:
            result = (result * power).rem(modulus_poly)
        power = (power * power).rem(modulus_poly)
        value >>= 1
    return result


# ---------------------------------------------------------------------------
# The exact rational dynamics and complete preperiodic graph.
# ---------------------------------------------------------------------------

x, y = sp.symbols("x y")
f = x**2 - sp.Rational(29, 16)
orbit_x = (
    -sp.Rational(7, 4),
    sp.Rational(5, 4),
    -sp.Rational(1, 4),
)
require(tuple(sp.factor(f.subs(x, value)) for value in orbit_x) == orbit_x[1:] + orbit_x[:1], "wrong rational three-cycle")

g = (y**2 - 29) / 4
candidate_numerators = tuple(range(-7, 8, 2))
candidate_graph = {
    value: int(g.subs(y, value)) for value in candidate_numerators
}
expected_graph = {
    -7: 5,
    -5: -1,
    -3: -5,
    -1: -7,
    1: -7,
    3: -5,
    5: -1,
    7: 5,
}
require(candidate_graph == expected_graph, "wrong bounded odd-numerator graph")
require(sp.factor((y**2 - 29) - 4 * y).subs(y, 9) > 0, "positive escape boundary failed")
require(sp.factor((y**2 - 29) - 4 * (-y)).subs(y, -9) > 0, "negative escape boundary failed")

# At an odd numerator, y^2-29 is 4 mod 8, so the next numerator is odd.
for residue in (1, 3, 5, 7):
    require(((residue * residue - 29) // 4) % 2 == 1, "odd numerator parity failed")


# ---------------------------------------------------------------------------
# Universal trace parameter and the SL_2 central-sign lift.
# ---------------------------------------------------------------------------

s, X = sp.symbols("s X")
c_s = -(s**2 + s + 2)
p_s = -(s**2 + 2 * s + 3)
r_s = -(s**3 + 2 * s**2 + 3 * s + 1)
cycle_cubic = X**3 - s * X**2 + p_s * X - r_s
mobius_numerator = (s + 1) * X - (s**2 + s + 1)
mobius_denominator = X - s
require(
    sp.factor((X**2 + c_s) - mobius_numerator / mobius_denominator)
    == sp.factor(cycle_cubic / mobius_denominator),
    "universal quadratic/Mobius interpolation failed",
)
L_s = sp.Matrix(
    [[s + 1, -(s**2 + s + 1)], [1, -s]]
)
I2 = sp.eye(2)
require(sp.factor(L_s.det()) == 1, "universal lift lost determinant one")
require(sp.factor(sp.trace(L_s)) == 1, "universal lift lost trace one")
require(sp.simplify(L_s**2 - L_s + I2) == sp.zeros(2), "Phi_6 identity failed")
require(sp.simplify(L_s**3 + I2) == sp.zeros(2), "central sign failed")
require(
    sp.factor(sp.discriminant(cycle_cubic, X)) == (4 * s**2 + 6 * s + 9) ** 2,
    "cycle-cubic discriminant changed",
)

# Specific y=4x model.
A = sp.Matrix([[1, -13], [1, 3]])
B = A / 4
require(A.det() == 16 and sp.trace(A) == 4, "integer lift invariants changed")
require(A**3 == -64 * I2, "integer lift cube changed")
require(B.det() == 1 and B**3 == -I2 and B**6 == I2, "order-six lift failed")

mobius_y = (y - 13) / (y + 3)
require(
    sp.factor(
        (g - mobius_y)
        - (y + 7) * (y - 5) * (y + 1) / (4 * (y + 3))
    )
    == 0,
    "specific interpolation divisor failed",
)

C = sp.Matrix([[0, -1], [1, 1]])
P = sp.Matrix([[4, 1], [0, 1]])
require(P * C * P.inv() == B, "standard Phi_6 conjugacy failed")
require(C**2 - C + I2 == sp.zeros(2) and C**3 == -I2, "standard lift failed")

z = sp.symbols("z")
h = z**2 + z / 2 - 2
projective_rotation = -1 / (z + 1)
require(
    sp.factor(
        (h - projective_rotation)
        - (z + 2) * (z - 1) * (z + sp.Rational(1, 2)) / (z + 1)
    )
    == 0,
    "normalized cycle interpolation failed",
)

# The six integral lifts form the complete q=48 hexagon.
lift_orbit = (
    sp.Matrix([-7, 1]),
    sp.Matrix([-5, -1]),
    sp.Matrix([2, -2]),
    sp.Matrix([7, -1]),
    sp.Matrix([5, 1]),
    sp.Matrix([-2, 2]),
)
for index, vector in enumerate(lift_orbit):
    require(B * vector == lift_orbit[(index + 1) % 6], "wrong integral six-lift")

Qmat = sp.Matrix([[1, 1], [1, 13]])
require(B.T * Qmat * B == Qmat, "invariant quadratic form failed")
require(all((vector.T * Qmat * vector)[0] == 48 for vector in lift_orbit), "hexagon norm changed")
integral_q48 = []
for xx in range(-16, 17):
    for yy in range(-8, 9):
        if xx * xx + 2 * xx * yy + 13 * yy * yy == 48:
            integral_q48.append((xx, yy))
require(set(integral_q48) == {tuple(map(int, vector)) for vector in lift_orbit}, "q=48 has an unlisted integral point")

H_X, H_Y = sp.symbols("H_X H_Y")
line_section = (H_X + 7 * H_Y) * (H_X - 5 * H_Y) * (H_X + H_Y)
BX, BY = B * sp.Matrix([H_X, H_Y])
require(sp.factor(line_section.subs({H_X: BX, H_Y: BY}, simultaneous=True)) == -line_section, "defining-section character is not minus one")


# ---------------------------------------------------------------------------
# Arithmetic-progression and Pythagorean uniqueness.
# ---------------------------------------------------------------------------

m, rr, cc = sp.symbols("m rr cc", nonzero=True)
ap_points = (m - rr, m + rr, m)
ap_equations = (
    sp.expand((ap_points[0] ** 2 + cc) - ap_points[1]),
    sp.expand((ap_points[1] ** 2 + cc) - ap_points[2]),
    sp.expand((ap_points[2] ** 2 + cc) - ap_points[0]),
)
ap_solution = {m: -sp.Rational(1, 4), rr: sp.Rational(3, 2), cc: -sp.Rational(29, 16)}
require(all(eq.subs(ap_solution) == 0 for eq in ap_equations), "AP cycle solution failed")
# The first two difference equations give rr*(-4m-1)=0 and rr*(2m+rr-1)=0.
require(sp.factor(ap_equations[0] - ap_equations[1]) == -rr * (4 * m + 1), "first AP difference changed")
require(sp.factor(ap_equations[1] - ap_equations[2]) == rr * (2 * m + rr - 1), "second AP difference changed")

a, b, hyp, D = sp.symbols("a b hyp D", positive=True)
D_tail = a**2 + 2 * b**2 - a * b
first_arrow_D = (a + b) ** 2 - b * hyp
require(
    sp.factor(first_arrow_D - D_tail) == b * (3 * a - b - hyp),
    "Pythagorean cycle compatibility changed",
)
require(
    sp.factor((3 * a - b) ** 2 - (a**2 + b**2)) == 2 * a * (4 * a - 3 * b),
    "3:4:5 forcing identity changed",
)
k = sp.symbols("k", positive=True)
pythagorean_solution = {a: 3 * k, b: 4 * k, hyp: 5 * k, D: 29 * k**2}
pythagorean_cycle = (-(a + b), hyp, -(b - a))
pythagorean_map = lambda value: sp.factor((value**2 - D) / b)
require(
    tuple(pythagorean_map(value).subs(pythagorean_solution) for value in pythagorean_cycle)
    == tuple(value.subs(pythagorean_solution) for value in pythagorean_cycle[1:] + pythagorean_cycle[:1]),
    "3:4:5 signed cycle failed",
)
require(D_tail.subs(pythagorean_solution) == 29 * k**2, "Pythagorean D changed")


# ---------------------------------------------------------------------------
# The exact THM-4134 horizontal-section intersection.
# ---------------------------------------------------------------------------

a0, rho, t = sp.symbols("a0 rho t", nonzero=True)
section_U = a0 / 2 + 16 * rho**2 / (9 * a0**2)
section_V = -rho - 64 * rho**3 / (27 * a0**3)
section_q = a0**3 / 2 + rho**2
elliptic_response = section_V**2 - section_U**3 + sp.Rational(3, 4) * a0**2 * section_U + a0**3 / 4
require(sp.factor(elliptic_response - section_q) == 0, "THM-4134 polynomial section failed")

rho_from_t = 3 * a0 * t / 4
U_t = sp.factor(section_U.subs(rho, rho_from_t))
V_t = sp.factor(section_V.subs(rho, rho_from_t))
q_t = sp.factor(section_q.subs(rho, rho_from_t))
require(sp.factor(U_t - (a0 / 2 + t**2)) == 0, "horizontal U normalization failed")
require(sp.factor(V_t - (-t**3 - 3 * a0 * t / 4)) == 0, "horizontal V normalization failed")
require(sp.factor(V_t.subs(a0, -48)) == -t * (t - 6) * (t + 6), "a=-48 zero divisor failed")

# Uniformly, the two outer roots of the section cubic collide in (U,V,q),
# while the central root lies on a different q-fibre for every a!=0.
r2 = -3 * a0 / 4
q_center = sp.factor(q_t.subs(t, 0))
q_outer = sp.factor(q_t.subs(t**2, r2))
U_center = sp.factor(U_t.subs(t, 0))
U_outer = sp.factor(U_t.subs(t**2, r2))
require(q_center == a0**3 / 2 and q_outer == 5 * a0**3 / 64, "uniform horizontal q split changed")
require(U_center == a0 / 2 and U_outer == -a0 / 4, "uniform horizontal U split changed")
require(sp.factor(q_center - q_outer) == 27 * a0**3 / 64, "horizontal fibres unexpectedly merged")

horizontal_roots = (-6, 6, 0)
require(tuple(root - 1 for root in horizontal_roots) == (-7, 5, -1), "horizontal/dynamics translation failed")
horizontal_images = tuple(
    (
        sp.factor(U_t.subs({a0: -48, t: root})),
        sp.factor(V_t.subs({a0: -48, t: root})),
        sp.factor(q_t.subs({a0: -48, t: root})),
    )
    for root in horizontal_roots
)
require(horizontal_images == ((12, 0, -8640), (12, 0, -8640), (-24, 0, -55296)), "horizontal image/fibre ledger changed")

T0, Z0 = sp.symbols("T0 Z0")
B_t = sp.Matrix([[sp.Rational(1, 2), -3], [sp.Rational(1, 4), sp.Rational(1, 2)]])
T1, Z1 = B_t * sp.Matrix([T0, Z0])
horizontal_line = T0 * (T0 - 6 * Z0) * (T0 + 6 * Z0)
require(B_t**3 == -I2 and B_t**6 == I2, "horizontal lift order changed")
require(B_t.T * sp.diag(1, 12) * B_t == sp.diag(1, 12), "horizontal norm failed")
require(
    sp.factor(horizontal_line.subs({T0: T1, Z0: Z1}, simultaneous=True)) == -horizontal_line,
    "horizontal defining section lost its sign",
)


# ---------------------------------------------------------------------------
# Exact period six over R and Q.
# ---------------------------------------------------------------------------

iterates = {0: y}
current = y
for depth in range(1, 7):
    current = sp.cancel(current.subs(y, g))
    iterates[depth] = current

phi6 = sp.cancel(
    (iterates[6] - y) * (iterates[1] - y)
    / ((iterates[3] - y) * (iterates[2] - y))
)
phi6_numerator = sp.Poly(
    sp.together(phi6).as_numer_denom()[0], y, domain=sp.QQ
).clear_denoms(convert=True)[1].primitive()[1]
require(phi6_numerator.degree() == 54, "exact-period-six degree changed")
for proper_depth in (1, 2, 3):
    lower = sp.Poly(sp.together(iterates[proper_depth] - y).as_numer_denom()[0], y, domain=sp.QQ)
    require(sp.gcd(phi6_numerator, lower).degree() == 0, "period-six factor meets a lower period")

# Explicit Rabin irreducibility certificate modulo 11.
prime = 11
phi_mod = sp.Poly(phi6_numerator.as_expr(), y, modulus=prime).monic()
frobenius = sp.Poly(y, y, modulus=prime)
frobenius_rows: dict[int, sp.Poly] = {}
for step in range(1, 55):
    frobenius = poly_pow_mod(frobenius, prime, phi_mod, prime)
    if step in (18, 27, 54):
        frobenius_rows[step] = frobenius
require(sp.gcd(phi_mod, frobenius_rows[18] - sp.Poly(y, y, modulus=prime)).degree() == 0, "mod-11 degree-18 Rabin gate failed")
require(sp.gcd(phi_mod, frobenius_rows[27] - sp.Poly(y, y, modulus=prime)).degree() == 0, "mod-11 degree-27 Rabin gate failed")
require(frobenius_rows[54] == sp.Poly(y, y, modulus=prime), "mod-11 Frobenius closure failed")
require(phi_mod.is_irreducible, "SymPy and explicit mod-11 irreducibility disagree")

real_period_six_points = phi6_numerator.count_roots(-sp.oo, sp.oo)
require(real_period_six_points == 18, "real period-six root count changed")
coefficient_payload = ",".join(str(int(value)) for value in phi6_numerator.all_coeffs())
phi6_sha256 = hashlib.sha256((coefficient_payload + "\n").encode()).hexdigest()

# The Mersenne exponent 63 is also the exact degree of the sixth-iterate
# derivative.  It is not the number of exact period-six points (that is 54).
g6_numerator = sp.together(iterates[6]).as_numer_denom()[0]
require(sp.degree(sp.diff(g6_numerator, y), y) == 63, "sixth-iterate derivative degree changed")
derivative_chain = sp.Rational(1, 2**6)
for depth in range(6):
    derivative_chain *= iterates[depth]
require(
    sp.factor(sp.diff(iterates[6], y) - derivative_chain) == 0,
    "sixth-iterate chain-rule product failed",
)


# ---------------------------------------------------------------------------
# The exact 63 / no-new-prime-factor sidecar and the F_2^2 palette.
# ---------------------------------------------------------------------------

require(2**6 - 1 == 63 and sp.factorint(63) == {3: 2, 7: 1}, "Mersenne factorization changed")
require(sp.n_order(2, 3) == 2 and sp.n_order(2, 7) == 3, "63 acquired a primitive order-six prime")
require(matrix_mod(A**3, 63) == matrix_mod(-I2, 63), "A^3 mod 63 changed")
require(matrix_mod(A**6, 63) == I2, "A does not have order dividing six mod 63")
for exponent in (1, 2, 3):
    require(matrix_mod(A**exponent, 63) != I2, "A has premature order mod 63")

require(sp.n_order(2, 63) == 6 and sp.n_order(-2, 63) == 6, "signed two lost order six mod 63")
inv4_mod63 = pow(4, -1, 63)
map_mod63 = lambda value: ((value * value - 29) * inv4_mod63) % 63
seen_mod63: set[int] = set()
cycles_mod63: list[tuple[int, ...]] = []
for seed in range(63):
    path: list[int] = []
    position: dict[int, int] = {}
    value = seed
    while value not in position and value not in seen_mod63:
        position[value] = len(path)
        path.append(value)
        value = map_mod63(value)
    if value in position:
        cycle_values = path[position[value] :]
        pivot = min(range(len(cycle_values)), key=lambda index: cycle_values[index])
        cycle_values = cycle_values[pivot:] + cycle_values[:pivot]
        cycles_mod63.append(tuple(cycle_values))
    seen_mod63.update(path)
require(
    sorted(cycles_mod63) == [(5, 62, 56), (14, 26, 20), (35, 47, 41)],
    "scalar mod-63 cycle census changed",
)

pair_sums_x = tuple(orbit_x[index] + orbit_x[(index + 1) % 3] for index in range(3))
pair_sums_y = tuple(4 * value for value in pair_sums_x)
require(pair_sums_x == (-sp.Rational(1, 2), 1, -2), "normalized pair sums changed")
require(pair_sums_y == (-2, 4, -8) and sp.prod(pair_sums_y) == 64, "integer pair-sum product changed")
cycle_polynomial_y = (y + 7) * (y - 5) * (y + 1)
require(sp.expand(cycle_polynomial_y - ((y + 3) * (y**2 - 33) + 64)) == 0, "cycle polynomial scalar-64 identity failed")

C2 = matrix_mod(C, 2)
palette = (sp.Matrix([1, 0]), sp.Matrix([0, 1]), sp.Matrix([1, 1]))
for index, color in enumerate(palette):
    require(matrix_mod(C2 * color, 2) == palette[(index + 1) % 3], "F2^2 palette rotation failed")

semantic_lines = (
    "finite_affine_rational_preperiodic_numerators=-7,-5,-3,-1,1,3,5,7;denominator=4",
    "unique_finite_affine_rational_cycle=-7/4->5/4->-1/4;no_affine_rational_period6",
    "universal_three_cycle_lift=trace1;det1;B^3=-I;B^6=I",
    "specific_A=1,-13;1,3;A^3=-64I;invariant_q=X^2+2XY+13Y^2",
    "integral_q48_hexagon=6_points;projective_cycle=-7,5,-1",
    "line_divisor_preserved;defining_cubic_character=-1",
    "AP_three_cycle_unique=c=-29/16",
    "Pythagorean_template_unique=3:4:5;D=29_up_to_scale",
    "horizontal_section=a:-48;t:-6,6,0;target_q:-8640,-8640,-55296",
    "period6_dynatomic=degree54;irreducible_mod11;real_points18;real_cycles3",
    "Mersenne63=derivative_degree;3^2*7;orders_of_2_at_primes=2,3;scalar_mod63_no_period6",
    "THM4138_DeltaV_degree15_16_wall=CLOSED;JC2=OPEN;no_new_JC2_consequence",
)
semantic_sha256 = hashlib.sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()

print("THM4146_PRIMARY_AUDIT")
print("finite_affine_rational_cycle=-7/4->5/4->-1/4->-7/4")
print("finite_affine_rational_preperiodic_numerators=" + ",".join(str(value) for value in candidate_numerators))
print("bounded_odd_graph=" + ",".join(f"{key}:{candidate_graph[key]}" for key in candidate_numerators))
print("universal_cycle_parameter=c=-(s^2+s+2);disc=(4s^2+6s+9)^2")
print("universal_lift=trace1;det1;L^2-L+I=0;L^3=-I;L^6=I")
print("specific_integer_lift=A[[1,-13],[1,3]];det=16;trace=4;A^3=-64I")
print("specific_invariant=q=X^2+2XY+13Y^2=(X+Y)^2+3*(2Y)^2;level=48")
print("integral_q48_hexagon=" + str(tuple(tuple(map(int, vector)) for vector in lift_orbit)))
print("model_defining_section_character=L(Bz)/L(z)=-1;not_BC_normal_character")
print("AP_uniqueness=(m,r,c)=(-1/4,3/2,-29/16)")
print("Pythagorean_uniqueness=(a,b,h,D)=k*(3,4,5,29k)")
print("horizontal_normalization=t=4rho/(3a);U=a/2+t^2;V=-t^3-3at/4")
print("horizontal_a_minus48_roots=-6,6,0;translated=-7,5,-1")
print("horizontal_images=(U,V,q)=" + str(horizontal_images))
print("horizontal_firewall=three_parameter_points_do_not_lie_in_one_q_fibre")
print("period6_dynatomic_degree=54;irreducible_mod11=PASS")
print("period6_real_points=18;period6_real_cycles=3;period6_affine_rational_points=0")
print("period6_coefficient_sha256=" + phi6_sha256)
print("mersenne63=degree((g^6)')=63=3^2*7;ord_3(2)=2;ord_7(2)=3;A_mod63_order=6")
print("pair_sums_x=-1/2,1,-2;pair_sums_y=-2,4,-8;integer_product=64")
print("scalar_mod63_cycles=" + str(tuple(sorted(cycles_mod63))) + ";no_length6")
print("palette=F2^2_nonzero_C3_matches_normalized_cycle_-2,1,-1/2")
print("semantic_sha256=" + semantic_sha256)
print("verdict=EXACT_ORDER6_LINEAR_LIFT_AND_HORIZONTAL_FIBRE_HOSTILE;NO_BC_ACTION;JC2_OPEN")
