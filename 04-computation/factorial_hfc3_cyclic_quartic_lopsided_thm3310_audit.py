"""Exact phase-boundary and lopsided-axis certificates for cyclic quartics.

This studies only the cubic moment of the five-dimensional C3-eigen quartic
cell from THM-3304.  It does not use the finite-modular higher moments and
does not claim projective emptiness in characteristic zero.
"""

from __future__ import annotations

import hashlib
import math
from fractions import Fraction

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def weak_compositions(total, length, prefix=()):
    if length == 1:
        yield prefix + (total,)
        return
    for first in range(total + 1):
        yield from weak_compositions(total - first, length - 1,
                                     prefix + (first,))


def multinomial(counts):
    result = math.factorial(sum(counts))
    for count in counts:
        result //= math.factorial(count)
    return result


def kernel_coefficient(r, s):
    """[u^r v^s] (1-u^3-v^3-3uv)^(-1)."""
    result = 0
    for ell in range(min(r, s) + 1):
        if (r - ell) % 3 or (s - ell) % 3:
            continue
        result += (multinomial(((r - ell) // 3,
                                (s - ell) // 3, ell)) * 3**ell)
    return result


def dirichlet_moment(r, s):
    coefficient = kernel_coefficient(r, s)
    return Fraction(2 * math.factorial(r) * math.factorial(s) * coefficient,
                    math.factorial(r + s + 2))


c = sp.symbols("c0:5")
rho = sp.symbols("rho")

# On e1=1 the basis is (b, a^2, a*b^2, b^4, a^3*b).  Each pair records
# its (a,b)-exponents.  Expanding the third moment therefore needs one exact
# Fourier--Dirichlet kernel coefficient per parameter monomial.
EXPONENT_PAIRS = ((0, 1), (2, 0), (1, 2), (0, 4), (3, 1))
raw_expression = 0
raw_coefficients = {}
for counts in weak_compositions(3, 5):
    r = sum(count * pair[0]
            for count, pair in zip(counts, EXPONENT_PAIRS))
    s = sum(count * pair[1]
            for count, pair in zip(counts, EXPONENT_PAIRS))
    coefficient = multinomial(counts) * dirichlet_moment(r, s)
    require(coefficient > 0, f"positive Fourier coefficient at {counts}")
    raw_coefficients[counts] = coefficient
    monomial = sp.Rational(coefficient.numerator, coefficient.denominator)
    for variable, exponent in zip(c, counts):
        monomial *= variable**exponent
    raw_expression += monomial

moment = sp.Poly(raw_expression, *c, domain=sp.QQ)
common_denominator, cleared = moment.clear_denoms()
content, primitive = sp.Poly(cleared, *c, domain=sp.ZZ).primitive()
require(int(common_denominator) == 2802800 and int(content) == 1,
        "primitive cubic normalization")
require(len(primitive.terms()) == 35,
        "the cubic moment has complete degree-three support")
require(all(int(coefficient) > 0 for coefficient in primitive.coeffs()),
        "strict coefficient positivity")

# Closed phase-boundary repair of THM-3304's open-half-plane argument.
# Rotate a closed pi/3 coefficient-phase arc to [0,pi/3].  Every cubic
# monomial then lies in the closed upper half-plane.  Vanishing would force
# every term onto its boundary.  Pure cubes force every live coefficient
# phase to be 0 or pi/3; complete mixed support supplies c_i^2*c_j between
# any two occupied endpoints, giving positive imaginary part.  If only one
# endpoint is occupied, every term is aligned and their positive magnitudes
# cannot cancel.
for i in range(5):
    require(primitive.coeff_monomial(c[i]**3) > 0,
            f"pure cube {i}")
    for j in range(5):
        if i != j:
            require(primitive.coeff_monomial(c[i]**2*c[j]) > 0,
                    f"boundary mixed term {(i, j)}")

# A magnitude-sensitive complement.  If |c_j| <= q*|c_i| for j != i,
# the sum of all non-pure terms is at most
#   |c_i|^3 (S1_i*q + S2_i*q^2 + S3_i*q^3).
# The following five positive margins make the pure cube strictly dominant,
# independently of all phases.
pure_weights = []
group_sums = []
dominance_polynomials = []
for i in range(5):
    weight = int(primitive.coeff_monomial(c[i]**3))
    sums = [0, 0, 0, 0]
    for powers, coefficient in primitive.terms():
        if powers[i] < 3:
            sums[3 - powers[i]] += int(coefficient)
    require(all(value > 0 for value in sums[1:]),
            f"positive dominance groups {i}")
    pure_weights.append(weight)
    group_sums.append(tuple(sums[1:]))
    dominance_polynomials.append(
        weight - sums[1]*rho - sums[2]*rho**2 - sums[3]*rho**3)

axis_radii = (sp.Rational(1, 10), sp.Rational(1, 16),
              sp.Rational(1, 18), sp.Rational(1, 23),
              sp.Rational(1, 22))
margins = tuple(sp.factor(poly.subs(rho, radius))
                for poly, radius in zip(dominance_polynomials, axis_radii))
require(all(margin > 0 for margin in margins),
        "strict lopsided-axis margins")
require(min(axis_radii) == sp.Rational(1, 23),
        "uniform second-largest-coordinate consequence")

# These exact sign brackets locate the unique positive zero of each
# decreasing dominance polynomial.  Thus they also show how close the clean
# reciprocal-integer boxes are to the sharp thresholds for this particular
# coefficientwise triangle-inequality certificate.
root_brackets = ((sp.Rational(27, 250), sp.Rational(1081, 10000)),
                 (sp.Rational(637, 10000), sp.Rational(319, 5000)),
                 (sp.Rational(293, 5000), sp.Rational(587, 10000)),
                 (sp.Rational(441, 10000), sp.Rational(221, 5000)),
                 (sp.Rational(237, 5000), sp.Rational(19, 400)))
for index, (poly, bracket) in enumerate(zip(dominance_polynomials,
                                             root_brackets)):
    lower, upper = bracket
    require(poly.subs(rho, lower) > 0 > poly.subs(rho, upper),
            f"root bracket {index}")
    require(sp.Poly(sp.diff(poly, rho), rho).all_coeffs()[0] < 0,
            f"negative leading derivative coefficient {index}")
    require(all(coefficient < 0
                for coefficient in sp.Poly(sp.diff(poly, rho), rho).all_coeffs()),
            f"strict decrease on positive axis {index}")


def coefficient_hash(poly):
    payload = ";".join(
        ",".join(map(str, powers)) + "=" + str(int(coefficient))
        for powers, coefficient in poly.terms())
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def rational_tuple(values):
    return "(" + ",".join(str(value) for value in values) + ")"


payload = "\n".join([
    "scope=the cubic moment I3 in the five-dimensional cyclic quartic C3-eigencell",
    "basis_on_e1=1=(b,a^2,a*b^2,b^4,a^3*b)",
    "I3=(" + str(primitive.as_expr()) + ")/2802800",
    "term_count=35; all coefficients strictly positive",
    "coefficient_sha256=" + coefficient_hash(primitive),
    "phase_boundary=closed coefficient-phase arc of width pi/3 is zero-free",
    "phase_consequence=the nonzero coefficient phases of every I3-zero have circular diameter strictly greater than pi/3",
    "pure_weights=" + rational_tuple(pure_weights),
    "axis_group_sums=" + repr(tuple(group_sums)),
    "axis_radii=" + rational_tuple(axis_radii),
    "axis_margins=" + rational_tuple(margins),
    "axis_rule=if max_{j!=i}|c_j|<=axis_radii[i]*|c_i| for some i, then I3!=0 for every phase choice",
    "global_modulus_consequence=at every projective I3-zero, the second-largest |c_i| is strictly greater than 1/23 of the largest",
    "sharp_triangle_threshold_brackets=" + repr(root_brackets),
    "interpretation=five full coamoeba fibres around the projective coordinate axes are excluded",
    "not_claimed=no characteristic-zero projective elimination or higher-moment common-zero theorem",
    "controls=complete support; pure and ordered mixed terms; exact positive margins; exact root sign brackets",
]) + "\n"
print(payload, end="")
print("payload_sha256=" + hashlib.sha256(payload.encode("ascii")).hexdigest())
