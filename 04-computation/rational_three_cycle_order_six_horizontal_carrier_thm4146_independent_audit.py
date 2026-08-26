#!/usr/bin/env python3
"""Independent exact audit for THM-4146, extending THM-4139.

This file imports no primary certificate code.  It uses direct rational
iteration, a separate polynomial-division construction of Phi_6, elementary
matrix multiplication, and a uniform symbolic horizontal-fibre calculation.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib

import sympy as sp


def gate(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def matmul(left: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]], right: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def matvec(matrix: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]], vector: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


# Direct rational dynamics.
def f(value: Fraction) -> Fraction:
    return value * value - Fraction(29, 16)


cycle = (Fraction(-7, 4), Fraction(5, 4), Fraction(-1, 4))
gate(tuple(f(value) for value in cycle) == cycle[1:] + cycle[:1], "cycle")

odd_graph: dict[int, int] = {}
for numerator in range(-7, 8, 2):
    next_numerator = (numerator * numerator - 29) // 4
    gate((numerator * numerator - 29) % 4 == 0, "integral odd graph")
    gate(next_numerator % 2 != 0, "odd graph parity")
    odd_graph[numerator] = next_numerator
gate(
    odd_graph
    == {-7: 5, -5: -1, -3: -5, -1: -7, 1: -7, 3: -5, 5: -1, 7: 5},
    "bounded graph",
)
for boundary in range(9, 100, 2):
    gate((boundary * boundary - 29) // 4 > boundary, "positive escape")
    gate((boundary * boundary - 29) // 4 > abs(-boundary), "negative escape")

# The delicate two-adic boundary is a complete unit-class calculation:
# for x=u/4 and every 2-adic unit u, u^2-29 has valuation exactly two.
for odd_unit_class in (1, 3, 5, 7):
    gate((odd_unit_class * odd_unit_class - 29) % 8 == 4, "two-adic unit class")


# Matrix lift, its orbit, invariant form, and central sign.
B = (
    (Fraction(1, 4), Fraction(-13, 4)),
    (Fraction(1, 4), Fraction(3, 4)),
)
identity = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)))
B2 = matmul(B, B)
B3 = matmul(B2, B)
B6 = matmul(B3, B3)
gate(B3 == ((Fraction(-1), Fraction(0)), (Fraction(0), Fraction(-1))), "B cube")
gate(B6 == identity, "B sixth power")

vectors = ((-7, 1), (-5, -1), (2, -2), (7, -1), (5, 1), (-2, 2))
for index, vector in enumerate(vectors):
    image = matvec(B, (Fraction(vector[0]), Fraction(vector[1])))
    gate(image == tuple(Fraction(value) for value in vectors[(index + 1) % 6]), "six-vector orbit")

def quadratic(vector: tuple[int, int]) -> int:
    xx, yy = vector
    return xx * xx + 2 * xx * yy + 13 * yy * yy


gate({quadratic(vector) for vector in vectors} == {48}, "invariant level")
enumerated = {
    (xx, yy)
    for xx in range(-20, 21)
    for yy in range(-10, 11)
    if quadratic((xx, yy)) == 48
}
gate(enumerated == set(vectors), "complete integral hexagon")

# Rebuild the defining-cubic character rather than inheriting its semantic row.
xx, yy = sp.symbols("xx yy")
line_section = (xx + 7 * yy) * (xx - 5 * yy) * (xx + yy)
bx = (xx - 13 * yy) / 4
by = (xx + 3 * yy) / 4
pulled_line_section = sp.expand(line_section.subs({xx: bx, yy: by}, simultaneous=True))
gate(sp.expand(pulled_line_section + line_section) == 0, "defining-section character")
gate(sp.expand((bx + 7 * by) - 2 * (xx + yy)) == 0, "first factor multiplier")
gate(sp.expand((bx - 5 * by) + (xx + 7 * yy)) == 0, "second factor multiplier")
gate(sp.expand((bx + by) - (xx - 5 * yy) / 2) == 0, "third factor multiplier")


# Independent symbolic derivation of the trace-parameter identities.
u, trace = sp.symbols("u trace")
cs = -(trace**2 + trace + 2)
poly = (
    u**3
    - trace * u**2
    - (trace**2 + 2 * trace + 3) * u
    + trace**3 + 2 * trace**2 + 3 * trace + 1
)
numerator = (trace + 1) * u - (trace**2 + trace + 1)
denominator = u - trace
gate(sp.cancel((u**2 + cs) - numerator / denominator - poly / denominator) == 0, "universal interpolant")
gate(sp.factor(sp.discriminant(poly, u)) == (4 * trace**2 + 6 * trace + 9) ** 2, "universal discriminant")


# The AP equations are solved from symbolic differences, not by importing the
# proposed values from the primary.
m_symbol, r_symbol, c_symbol = sp.symbols("m_symbol r_symbol c_symbol")
e1 = (m_symbol - r_symbol) ** 2 + c_symbol - (m_symbol + r_symbol)
e2 = (m_symbol + r_symbol) ** 2 + c_symbol - m_symbol
e3 = m_symbol**2 + c_symbol - (m_symbol - r_symbol)
gate(sp.factor(e1 - e2) == -r_symbol * (4 * m_symbol + 1), "first AP difference")
gate(sp.factor(e2 - e3) == r_symbol * (2 * m_symbol + r_symbol - 1), "second AP difference")
m = Fraction(-1, 4)
r = Fraction(3, 2)
c = Fraction(-29, 16)
ap = (m - r, m + r, m)
gate(tuple(value * value + c for value in ap) == ap[1:] + ap[:1], "AP solution")
gate(-4 * m == 1 and 2 * m + r == 1 and m * m + c == m - r, "AP forcing solution")


# Pythagorean template: c_h=3a-b and Pythagoras force 4a=3b.
a, b = sp.symbols("a b", positive=True)
hyp = 3 * a - b
gate(sp.factor(hyp**2 - a**2 - b**2) == 2 * a * (4 * a - 3 * b), "Pythagorean forcing")
scale = sp.symbols("scale", positive=True)
D = a**2 + 2 * b**2 - a * b
gate(sp.factor(D.subs({a: 3 * scale, b: 4 * scale})) == 29 * scale**2, "Pythagorean 29")


# Horizontal section and the uniform q-fibre obstruction.
a0, t = sp.symbols("a0 t", nonzero=True)
rho = 3 * a0 * t / 4
U = a0 / 2 + 16 * rho**2 / (9 * a0**2)
V = -rho - 64 * rho**3 / (27 * a0**3)
q = a0**3 / 2 + rho**2
E = V**2 - U**3 + sp.Rational(3, 4) * a0**2 * U + a0**3 / 4
gate(sp.factor(E - q) == 0, "horizontal section")
gate(sp.factor(U - (a0 / 2 + t**2)) == 0, "horizontal U")
gate(sp.factor(V - (-t**3 - 3 * a0 * t / 4)) == 0, "horizontal V")
gate(sp.factor(V.subs(a0, -48)) == -t * (t - 6) * (t + 6), "horizontal roots")

r_squared = -3 * a0 / 4
q_center = sp.factor(q.subs(t, 0))
q_outer = sp.factor(q.subs(t**2, r_squared))
gate(q_center == a0**3 / 2 and q_outer == 5 * a0**3 / 64, "uniform fibre split")
gate(sp.factor(q_center - q_outer) == 27 * a0**3 / 64, "uniform noncollision")


# Build Phi_6 by exact quotient in QQ[y], then certify its two load-bearing
# properties independently.
y = sp.symbols("y")
gy = (y**2 - 29) / 4
iterate = y
F: dict[int, sp.Poly] = {}
for depth in range(1, 7):
    iterate = sp.cancel(iterate.subs(y, gy))
    numerator_depth = sp.together(iterate - y).as_numer_denom()[0]
    F[depth] = sp.Poly(numerator_depth, y, domain=sp.QQ)

product = F[6] * F[1]
quotient_3, remainder_3 = sp.div(product, F[3], domain=sp.QQ)
gate(remainder_3.is_zero, "Phi6/F3 quotient")
quotient_2, remainder_2 = sp.div(quotient_3, F[2], domain=sp.QQ)
gate(remainder_2.is_zero, "Phi6/F2 quotient")
phi6 = quotient_2.clear_denoms(convert=True)[1].primitive()[1]
gate(phi6.degree() == 54, "Phi6 degree")
gate(sp.Poly(phi6.as_expr(), y, modulus=11).is_irreducible, "Phi6 mod11")
gate(phi6.count_roots(-sp.oo, sp.oo) == 18, "Phi6 real roots")
payload = ",".join(str(int(value)) for value in phi6.all_coeffs())
phi6_sha256 = hashlib.sha256((payload + "\n").encode()).hexdigest()
gate(phi6_sha256 == "e48d6769fc26e83b238513c8dbebbb16dcdf1533a8680ea0bc6bc08e97070c8a", "Phi6 coefficient digest")


# Mersenne and scalar modulo-63 hostile.
gate(2**6 - 1 == 63 and sp.factorint(63) == {3: 2, 7: 1}, "Mersenne factorization")
gate(sp.n_order(2, 3) == 2 and sp.n_order(2, 7) == 3, "primitive-prime audit")
inv4 = pow(4, -1, 63)
step = lambda value: ((value * value - 29) * inv4) % 63
seen: set[int] = set()
cycles: list[tuple[int, ...]] = []
for seed in range(63):
    path: list[int] = []
    positions: dict[int, int] = {}
    value = seed
    while value not in positions and value not in seen:
        positions[value] = len(path)
        path.append(value)
        value = step(value)
    if value in positions:
        row = path[positions[value] :]
        pivot = min(range(len(row)), key=lambda index: row[index])
        cycles.append(tuple(row[pivot:] + row[:pivot]))
    seen.update(path)
gate(sorted(cycles) == [(5, 62, 56), (14, 26, 20), (35, 47, 41)], "mod63 scalar cycles")

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
gate(semantic_sha256 == "014f11897d8ca216faa897f117db118a9760b791d393d664289b621d99456d3c", "semantic digest")

print("THM4146_INDEPENDENT_AUDIT")
print("rational_graph=PASS;escape=PASS;two_adic_unit_classes=PASS")
print("matrix_central_sign=PASS;integral_q48_hexagon=COMPLETE;model_defining_section_character=PASS")
print("universal_trace_parameter=PASS;AP_difference_derivation=PASS;Pythagorean_uniqueness=PASS")
print("horizontal_section=PASS;uniform_q_fibre_split=27*a^3/64_nonzero")
print("period6_dynatomic=degree54;irreducible_mod11;real_points18;real_cycles3")
print("period6_coefficient_sha256=" + phi6_sha256)
print("scalar_mod63_cycles=" + str(tuple(sorted(cycles))) + ";all_length3")
print("semantic_sha256=" + semantic_sha256)
print("verdict=ACCEPT;SECTION_POLYNOMIAL_COLLISION_NOT_ONE_FIBRE;JC2_OPEN")
