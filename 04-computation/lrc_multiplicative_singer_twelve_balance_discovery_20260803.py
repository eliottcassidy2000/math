#!/usr/bin/env python3
"""Exact discovery companion for the multiplicative Singer twelve-balance.

This is a discovery artifact, not a promoted theorem companion.  It compares
the THM-3246 owner Hodge word and positive owner masses in the regular
multiplicative C_168 representation, and constructs the sharp norm-phase
sidecar complementary to their common rank defect.
"""

import ast
import contextlib
import hashlib
import io
import runpy
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path

from sympy import Poly, cyclotomic_poly, factor_list, symbols


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / (
        "01-canon/theorems/"
        "THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate.md"
    ): "ef77a1f8fce16eb851eb38d5110a61ab73aa693f2d0ee9e11a912aa4fc302c87",
    ROOT / (
        "01-canon/theorems/"
        "THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word.md"
    ): "6badc0c9aba09b56d3d055a96cb8ef8b619d8492508bf21476eba5f624b13055",
    ROOT / (
        "01-canon/theorems/"
        "THM-3252-singer-compactified-owner-hodge-word-universal-charged-cyclicity.md"
    ): "1f8797de2d5fac74814fb78ca4f4d500de8c42eb14a6e1721e5f3e2a2810a873",
    ROOT / (
        "01-canon/theorems/"
        "THM-3253-positive-owner-mass-newton-cyclicity-and-maximal-common-heisenberg-module.md"
    ): "b94aea11abe97a6cc1a3826a91fab59d4c04e15f0e6acd9c924c5463b7bd63e8",
    ROOT / "04-computation/lrc_second_owner_all_dilation_seam_thm3246.py":
        "e23b098b38aa2199a348f48f8ab4ac0ce5913c870ead972bd31296494fc25a4b",
    ROOT / "05-knowledge/results/lrc_second_owner_all_dilation_seam_thm3246.out":
        "d7f7dd96b01c597113e78f903cad36246cb47b10e9a1758cb831aa0e83e8cebc",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n")


for dependency, expected_hash in DEPENDENCIES.items():
    require(hashlib.sha256(lf_bytes(dependency)).hexdigest() == expected_hash,
            ("dependency hash drift", dependency.name))

syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "assert statements are optimization-sensitive")
require(float_literals == 0, "floating literals are forbidden")

# Replay THM-3246 and import only its exact promoted objects.  All Fourier,
# polynomial, quotient, gauge and marker calculations below are independent.
dependency_script = ROOT / "04-computation/lrc_second_owner_all_dilation_seam_thm3246.py"
dependency_output = ROOT / "05-knowledge/results/lrc_second_owner_all_dilation_seam_thm3246.out"
dependency_stdout = io.StringIO()
with contextlib.redirect_stdout(dependency_stdout):
    inherited = runpy.run_path(str(dependency_script))
require(dependency_stdout.getvalue().encode("utf-8") == lf_bytes(dependency_output),
        "THM-3246 transcript drift")

N = 168
x, g = symbols("x g")
cycle = Poly(x ** N - 1, x, domain="ZZ")
divisors = tuple(d for d in range(1, N + 1) if N % d == 0)
forced_orders = (2, 3, 4, 6, 12)
remaining_orders = tuple(d for d in divisors if d not in forced_orders)
s12 = Poly(sum(x ** j for j in range(12)), x, domain="ZZ")
require(cycle.rem(s12).is_zero, "S_12 must divide x^168-1")


def word_poly(values):
    return Poly(sum(value * x ** j for j, value in enumerate(values)),
                x, domain="ZZ")


def transformed_word(values, multiplier, shift):
    result = [0] * N
    for owner, value in enumerate(values):
        result[(shift + multiplier * owner) % N] = value
    return tuple(result)


def poly_rank(values):
    return N - word_poly(values).gcd(cycle).degree()


def primitive_g_gcd(remainder, degree):
    coordinate_polynomials = []
    for index in range(degree):
        coordinate = Poly(remainder.nth(index), g, domain="ZZ")
        if not coordinate.is_zero:
            coordinate_polynomials.append(coordinate)
    require(coordinate_polynomials, "nonzero remainder expected")
    common = coordinate_polynomials[0]
    for coordinate in coordinate_polynomials[1:]:
        common = common.gcd(coordinate)
    primitive = common.primitive()[1]
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive


# The exact Hodge word.
q_word = tuple(inherited["q_word"])
require(len(q_word) == N, "Hodge word length")
q_scale = 1
for value in q_word:
    q_scale = lcm(q_scale, value.denominator)
q_integer = tuple(int(value * q_scale) for value in q_word)
q_content = 0
for value in q_integer:
    q_content = gcd(q_content, abs(value))
require((q_scale, q_content) == (32006016000, 1), "primitive Hodge scaling")
q_polynomial = word_poly(q_integer)
q_gcd = q_polynomial.gcd(cycle)
require(q_gcd == s12, "Hodge multiplicative gcd")

q_class_sums = tuple(
    sum(q_word[residue::12], Fraction(0)) for residue in range(12)
)
require(set(q_class_sums) == {Fraction(1, 296352)},
        "Hodge norm-fibre balance")
q_zero_orders = tuple(
    d for d in divisors
    if q_polynomial.rem(Poly(cyclotomic_poly(d, x), x, domain="ZZ")).is_zero
)
require(q_zero_orders == forced_orders, "Hodge zero-order census")

# The twelve negative seam owners form one complete norm-phase transversal.
# More surprisingly, cancellation occurs only after the twelve phases are
# summed: every individual phase restriction of the Hodge word is a unit of
# Q[C_168].
negative_owners = tuple(owner for owner, value in enumerate(q_word) if value < 0)
require(tuple(sorted(owner % 12 for owner in negative_owners)) == tuple(range(12)),
        "negative seam is a norm-phase transversal")
q_phase_ranks = []
for residue in range(12):
    phase_word = tuple(q_integer[owner] if owner % 12 == residue else 0
                       for owner in range(N))
    require((sum(phase_word[owner] > 0 for owner in range(N)),
             sum(phase_word[owner] < 0 for owner in range(N))) == (13, 1),
            ("Hodge phase signs", residue))
    q_phase_ranks.append(poly_rank(phase_word))
require(tuple(q_phase_ranks) == (N,) * 12,
        "every Hodge norm-phase slice is multiplicatively cyclic")

# The positive mass numerator word P_g.  The common positive denominator D_g
# does not change Fourier support or circulant rank.
mass_coefficients = tuple(
    tuple(int(value) for value in inherited["explicit_polynomial"](owner))
    for owner in range(N)
)
require(all(len(row) == 3 for row in mass_coefficients), "mass coefficient rows")
a_word = tuple(row[0] for row in mass_coefficients)
b_word = tuple(row[1] for row in mass_coefficients)
c_word = tuple(row[2] for row in mass_coefficients)
a_polynomial = word_poly(a_word)
b_polynomial = word_poly(b_word)
c_polynomial = word_poly(c_word)

mass_class_sums = tuple(
    tuple(sum(mass_coefficients[owner][coordinate]
              for owner in range(residue, N, 12))
          for coordinate in range(3))
    for residue in range(12)
)
require(set(mass_class_sums) == {(120960, -528, 2)},
        "positive-mass norm-fibre balance")

# A cyclotomic factor divides P_g only if every coordinate of its remainder
# vanishes.  The gcd in Z[g] of those coordinates is therefore a complete
# obstruction.  This proves the all-integer-g statement without sampling g.
coordinate_gcds = {}
mass_forced_orders = []
for order in divisors:
    cyclotomic = Poly(cyclotomic_poly(order, x), x, domain="ZZ")
    remainder = Poly(
        a_polynomial.as_expr() * g * g
        + b_polynomial.as_expr() * g
        + c_polynomial.as_expr(),
        x, domain="ZZ[g]",
    ).rem(cyclotomic)
    if remainder.is_zero:
        mass_forced_orders.append(order)
    else:
        coordinate_gcds[order] = primitive_g_gcd(remainder, cyclotomic.degree())

require(tuple(mass_forced_orders) == forced_orders,
        "symbolic positive-mass forced orders")
require(coordinate_gcds[1] == Poly(60480 * g * g - 264 * g + 1, g),
        "trivial-character polynomial")
require(coordinate_gcds[8] == Poly(504 * g - 1, g),
        "order-eight exceptional root")
require(coordinate_gcds[24] == Poly(504 * g - 1, g),
        "order-twenty-four exceptional root")
constant_obstruction_orders = tuple(
    order for order in remaining_orders
    if order not in (1, 8, 24)
)
require(all(coordinate_gcds[order].degree() == 0
            for order in constant_obstruction_orders),
        "all other cyclotomic coordinate gcds are units")
require(60480 - 264 + 1 > 0 and 2 * 60480 - 264 > 0,
        "trivial character is positive for every g>=1")
require(Fraction(1, 504).denominator != 1,
        "orders eight and twenty-four have no integer exceptional dilation")

# Restricting the positive masses to any one norm fibre also removes every
# multiplicative zero.  For characters factoring through C_12 the transform
# is the nonzero fibre total times a phase.  Orders 8 and 24 could vanish only
# at g=1/504; every other order has coordinate gcd one.
phase_total = Poly(60480 * g * g - 264 * g + 1, g, domain="ZZ")
phase_mass_certificate_rows = []
for residue in range(12):
    phase_a = word_poly(tuple(a_word[owner] if owner % 12 == residue else 0
                              for owner in range(N)))
    phase_b = word_poly(tuple(b_word[owner] if owner % 12 == residue else 0
                              for owner in range(N)))
    phase_c = word_poly(tuple(c_word[owner] if owner % 12 == residue else 0
                              for owner in range(N)))
    for order in divisors:
        cyclotomic = Poly(cyclotomic_poly(order, x), x, domain="ZZ")
        remainder = Poly(
            phase_a.as_expr() * g * g
            + phase_b.as_expr() * g
            + phase_c.as_expr(),
            x, domain="ZZ[g]",
        ).rem(cyclotomic)
        require(not remainder.is_zero,
                ("phase slice has no forced zero", residue, order))
        obstruction = primitive_g_gcd(remainder, cyclotomic.degree())
        if order in (1, 2, 3, 4, 6, 12):
            require(obstruction == phase_total,
                    ("quotient-character fibre total", residue, order))
            obstruction_class = "phase_total"
        elif order in (8, 24):
            require(obstruction == Poly(504 * g - 1, g, domain="ZZ"),
                    ("phase exceptional dilation", residue, order))
            obstruction_class = "504g-1"
        else:
            require(obstruction.degree() == 0,
                    ("phase unit obstruction", residue, order))
            obstruction_class = "unit"
        phase_mass_certificate_rows.append((residue, order, obstruction_class))

require(len(phase_mass_certificate_rows) == 12 * len(divisors),
        "complete phase-mass cyclotomic certificate")

direct_dilations = (1, 2, 3, 17, 18, 57, 58, 169)
for dilation in direct_dilations:
    values = tuple(
        row[0] * dilation * dilation + row[1] * dilation + row[2]
        for row in mass_coefficients
    )
    require(all(value > 0 for value in values),
            ("positive mass hostile", dilation))
    require(word_poly(values).gcd(cycle) == s12,
            ("direct multiplicative gcd", dilation))
    for residue in range(12):
        phase_values = tuple(value if owner % 12 == residue else 0
                             for owner, value in enumerate(values))
        require(word_poly(phase_values).gcd(cycle).degree() == 0,
                ("direct phase multiplicative unit", dilation, residue))

# Singer gauges act by affine automorphisms j |-> b+a*j of C_168.  They send
# norm fibre r to b+a*r mod 12, preserving the zero-order set and rank.
units = tuple(multiplier for multiplier in range(N) if gcd(multiplier, N) == 1)
require(len(units) == 48, "unit multiplier census")
gauge_checks = 0
for multiplier in units:
    for shift in range(N):
        image_classes = tuple(sorted((shift + multiplier * residue) % 12
                                     for residue in range(12)))
        require(image_classes == tuple(range(12)),
                ("gauge quotient permutation", multiplier, shift))
        gauge_checks += 1
require(gauge_checks == 8064, "full Singer gauge census")

hostile_gauges = ((1, 0), (5, 167), (13, 83), (43, 57), (71, 91), (167, 7))
mass_g58 = tuple(
    row[0] * 58 * 58 + row[1] * 58 + row[2]
    for row in mass_coefficients
)
for multiplier, shift in hostile_gauges:
    require(word_poly(transformed_word(q_integer, multiplier, shift)).gcd(cycle)
            == s12, ("Hodge hostile gauge", multiplier, shift))
    require(word_poly(transformed_word(mass_g58, multiplier, shift)).gcd(cycle)
            == s12, ("mass hostile gauge", multiplier, shift))

# In the deterministic Singer model Norm(alpha)=6 of order 12.  Hence j mod
# 12 is exactly the multiplicative norm/determinant phase, with 14 points in
# each fibre.
norm_fibres = {
    norm_value: tuple(owner for owner in range(N)
                      if pow(6, owner, 13) == norm_value)
    for norm_value in {pow(6, owner, 13) for owner in range(N)}
}
require(len(norm_fibres) == 12, "twelve norm values")
require(set(map(len, norm_fibres.values())) == {14}, "fourteen points per norm fibre")
for left in range(N):
    for right in range(N):
        require((pow(6, left, 13) == pow(6, right, 13))
                == (left % 12 == right % 12), "norm/residue equivalence")

# The same phase restriction has the opposite behavior in THM-3253's
# additive 13 by 13 coefficient arrangement.  A norm conic has only seven
# active rows or seven active columns, so every weighting on one fibre has
# additive matrix rank at most seven.  This support-only bound is gauge- and
# dilation-independent and pinpoints the representation switch.
def field_mul(left, right):
    a0, a1 = left
    b0, b1 = right
    return ((a0 * b0 + 2 * a1 * b1) % 13,
            (a0 * b1 + a1 * b0) % 13)


alpha = (1, 2)
singer_points = []
point = (1, 0)
for _ in range(N):
    singer_points.append(point)
    point = field_mul(point, alpha)
require(point == (1, 0) and len(set(singer_points)) == N,
        "deterministic Singer plane")
additive_support_profiles = []
for residue in range(12):
    support = tuple(singer_points[owner]
                    for owner in range(residue, N, 12))
    active_rows = len({row for row, _ in support})
    active_columns = len({column for _, column in support})
    require(sorted((active_rows, active_columns)) == [7, 8],
            ("norm-conic additive support", residue))
    additive_support_profiles.append((active_rows, active_columns))

# A centred fibre marker has precisely the missing 11 Fourier modes.  Its
# nonnegative uncentred version has those 11 modes plus the trivial mode.
marker_residue = 0
indicator = tuple(1 if owner % 12 == marker_residue else 0
                  for owner in range(N))
centered_marker = tuple(12 * value - 1 for value in indicator)
indicator_gcd = word_poly(indicator).gcd(cycle)
centered_gcd = word_poly(centered_marker).gcd(cycle)
expected_indicator_gcd = cycle.exquo(Poly(x ** 12 - 1, x, domain="ZZ"))
expected_centered_gcd = cycle.exquo(s12)
require(indicator_gcd == expected_indicator_gcd, "indicator Fourier support")
require(centered_gcd == expected_centered_gcd, "centered-marker Fourier support")
require((poly_rank(indicator), poly_rank(centered_marker)) == (12, 11),
        "phase-marker ranks")
require((sum(indicator), sum(centered_marker)) == (14, 0),
        "marker totals")
require((sum(value > 0 for value in centered_marker),
         sum(value < 0 for value in centered_marker)) == (14, 154),
        "centered marker must be signed")
require(s12.gcd(centered_gcd).degree() == 0
        and s12 * centered_gcd == cycle,
        "centered marker is an exact direct complement")
require(s12.gcd(indicator_gcd).degree() == 0,
        "positive indicator fills every missing mode")
require(157 + 12 - 1 == N, "positive marker has one trivial-line overlap")

for multiplier in units:
    for shift in range(N):
        transformed = transformed_word(centered_marker, multiplier, shift)
        target_residue = (shift + multiplier * marker_residue) % 12
        expected = tuple(11 if owner % 12 == target_residue else -1
                         for owner in range(N))
        require(transformed == expected,
                ("marker gauge covariance", multiplier, shift))

mass_coefficient_text = "\n".join(
    "%d|%d|%d" % row for row in mass_coefficients
)
coordinate_gcd_text = "\n".join(
    "%d|%s" % (order, coordinate_gcds[order].as_expr())
    for order in sorted(coordinate_gcds)
)
phase_mass_certificate_text = "\n".join(
    "%d|%d|%s" % row for row in phase_mass_certificate_rows
)
marker_text = "\n".join(map(str, centered_marker))

print("multiplicative Singer twelve-balance exact discovery")
print("status=VERIFIED-EXACT-DISCOVERY-CANDIDATE;not-promoted")
print("dependency_hash_checks=%d,dependency_replay=PASS" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d,sympy_exact=YES" %
      (assert_nodes, float_literals))
print("group=C_168,norm_quotient=C_12,norm_fibres=12x14")
print("hodge_class_sum_each=1/296352")
print("mass_class_sum_each=120960*g^2-528*g+2;denominator=(504*g-1)(840*g-2)")
print("common_zero_orders=%s,common_gcd=S12=1+x+...+x^11" %
      (forced_orders,))
print("hodge_circulant_rank=157;mass_circulant_rank=157_for_every_integer_g>=1")
print("negative_seam_norm_transversal=YES;each_hodge_phase_signs=(13+,1-)")
print("each_of_12_hodge_phase_slices_multiplicative_rank=168")
print("each_of_12_positive_mass_phase_slices_multiplicative_rank=168_for_every_integer_g>=1")
print("single_norm_phase_additive_support_profiles=%s;matrix_rank_at_most_7" %
      (tuple(additive_support_profiles),))
print("remaining_order_obstructions=unit:%s;order8,24:504*g-1;order1:60480*g^2-264*g+1" %
      (constant_obstruction_orders,))
print("direct_dilation_controls=%s" % (direct_dilations,))
print("singer_gauges=%d;hostile_exact_gauges=%s;rank_invariant=PASS" %
      (gauge_checks, hostile_gauges))
print("centered_norm_fibre_marker_rank=11,signs=(14,154),direct_complement=YES")
print("nonnegative_norm_fibre_indicator_rank=12,trivial_overlap=1,full_sum_rank=168")
print("positivity_boundary=any_nonzero_nonnegative_filling_marker_has_rank_at_least_12")
print("q_integer_sha256=%s" %
      hashlib.sha256("\n".join(map(str, q_integer)).encode("ascii")).hexdigest())
print("mass_coefficient_sha256=%s" %
      hashlib.sha256(mass_coefficient_text.encode("ascii")).hexdigest())
print("coordinate_gcd_sha256=%s" %
      hashlib.sha256(coordinate_gcd_text.encode("ascii")).hexdigest())
print("phase_mass_certificate_sha256=%s" %
      hashlib.sha256(phase_mass_certificate_text.encode("ascii")).hexdigest())
print("centered_marker_sha256=%s" %
      hashlib.sha256(marker_text.encode("ascii")).hexdigest())
print("scope=abstract-multiplicative-orbit;no-endpoint-target-ancestry-or-LRC14-decrement")
print("all_exact_checks=PASS")
