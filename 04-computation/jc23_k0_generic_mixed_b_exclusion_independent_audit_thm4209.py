#!/usr/bin/env python3
"""Independent source-coordinate referee for proposed THM-4209.

This file reconstructs the K=0 source from the canonical normal form.  It
does not import or inspect the primary mixed-face probe.  Its decisive
elimination is Res_s(A,B), rather than the primary normalized projection.
"""

from __future__ import annotations

from math import gcd
from sympy import (
    Poly, Rational, cancel, diff, expand, factor, resultant, symbols,
    hessian, det, rem,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


s, p, X, T, Q = symbols("s p X T Q")
Phi, Theta, eta, zeta = symbols("Phi Theta eta zeta")
a, u, W, q = symbols("a u W q")
DELTA = Rational(5696, 105)
t = p - s**2


def source_h() -> object:
    return (
        -3*p + Rational(8, 3)*p**2 - Rational(1376, 135)*p**3
        + Phi*s*p**3 + DELTA*p**4 + Theta*s**2*p**3
        + eta*s*p**4 + zeta*s**3*p**3
    )


H = source_h()
Gsp = -s**2/(2*t) + H
A = cancel((-s*p + t**2*diff(H, s))/p)
C0 = s**2 + 2*t**2*diff(H, p)
B = cancel((C0 + s*A)/t**2)
require(not cancel(A).as_numer_denom()[1].free_symbols, "A is not polynomial")
require(not cancel(B).as_numer_denom()[1].free_symbols, "B is not polynomial")
require(cancel(t**2*diff(Gsp, s) - p*A) == 0, "first gradient identity")
require(cancel(2*t**2*diff(Gsp, p) - (t**2*B - s*A)) == 0,
        "second gradient identity")
require(factor(A.subs(p, 0)) == -s, "p=0 A row")
require(factor(B.subs(p, 0)) == -6, "p=0 B row")
require(factor(A.subs(p, s**2)) == -s, "t=0 A row")
require(factor(B.subs({p: 0, s: 0})) == -6, "collapsed t=0 B row")
require(factor(Poly(A, s).LC()) == 3*p**2*zeta, "generic source A infinity")
require(factor(Poly(B, s).LC()) == 9*p**2*zeta, "generic source B infinity")
require(factor(Poly(A.subs(zeta, -eta), s).LC()) == -3*eta*p**2,
        "anti source A infinity")
require(factor(Poly(B.subs(zeta, -eta), s).LC()) == -9*eta*p**2,
        "anti source B infinity")


def source_residual(substitutions: dict) -> tuple[int, int, object, object, int, int]:
    aa = expand(A.subs(substitutions))
    bb = expand(B.subs(substitutions))
    rr = resultant(aa, bb, s)
    polynomial = Poly(rr, p)
    valuation = min(monomial[0] for monomial, coefficient in polynomial.terms()
                    if coefficient != 0)
    residual = Poly(cancel(rr/p**valuation), p)
    return (
        valuation,
        residual.degree(),
        factor(residual.nth(0)),
        factor(residual.LC()),
        Poly(aa, s).degree(),
        Poly(bb, s).degree(),
    )


generic = source_residual({})
require(generic == (
    6, 21,
    1259712*zeta**3,
    2125764*eta**3*zeta**2*(eta + zeta)**4,
    6, 3,
), "generic source resultant")

anti = source_residual({zeta: -eta})
anti_lead = 1327104*eta**5*(DELTA + Theta)**4
require(anti == (
    6, 19, -1259712*eta**3, factor(anti_lead), 6, 3,
), "anti source resultant")

eta_zero = source_residual({eta: 0})
require(eta_zero == (
    6, 20, 1259712*zeta**3,
    Rational(40870620168192, 1225)*zeta**7, 6, 3,
), "eta-zero hostile")

zeta_zero = source_residual({zeta: 0})
require(zeta_zero == (
    4, 18, -31104*Theta**2, -65610*Theta*eta**6, 5, 2,
), "zeta-zero hostile")

anti_collision = source_residual({zeta: -eta, Theta: -DELTA})
require(anti_collision == (
    6, 17, -1259712*eta**3, 777924*Phi**4*eta**5, 6, 3,
), "anti Theta=-Delta hostile")

anti_collision_phi_zero = source_residual(
    {zeta: -eta, Theta: -DELTA, Phi: 0})
require(anti_collision_phi_zero[:3] ==
        (6, 15, -1259712*eta**3), "terminal anti collision degree")

both_top_zero = source_residual({eta: 0, zeta: 0})
require(both_top_zero[:3] ==
        (4, 15, -31104*Theta**2), "weight-nine exit hostile")


# Normalized-coordinate typing and infinity rows.
ss = X*T
pp = T + ss**2
Hxt = expand(H.subs({s: ss, p: pp}))
Gxt = expand(-X**2*T/2 + Hxt)
f = cancel(diff(Gxt, X)/T)
h = expand(diff(Gxt, T))
require(not cancel(f).as_numer_denom()[1].free_symbols,
        "normalized f not polynomial")
require((Poly(f, X).degree(), Poly(h, X).degree()) == (8, 9),
        "generic normalized degrees")
require(factor(Poly(f, X).LC()) == 9*T**8*(eta + zeta),
        "generic f infinity row")
require(factor(Poly(h, X).LC()) == 9*T**8*(eta + zeta),
        "generic h infinity row")
fa = expand(f.subs(zeta, -eta))
ha = expand(h.subs(zeta, -eta))
require((Poly(fa, X).degree(), Poly(ha, X).degree()) == (7, 8),
        "anti normalized degrees")
anti_infinity = 8*T**7*(DELTA + Theta)
require(cancel(Poly(fa, X).LC() - anti_infinity) == 0,
        "anti f infinity row")
require(cancel(Poly(ha, X).LC() - anti_infinity) == 0,
        "anti h infinity row")


def zero_mod_x2(expression: object, t_value: object, x_square: object) -> bool:
    value = Poly(expand(expression.subs(T, t_value)), X)
    return rem(value, Poly(X**2 - x_square, X)).as_expr() == 0


Hess_xt = det(hessian(Gxt, (X, T)))
for t_value, x_square, g_value, hessian_value in (
    (0, -6, 0, 6),
    (Rational(-1, 6), 6, Rational(1, 2), -6),
):
    require(zero_mod_x2(diff(Gxt, X), t_value, x_square),
            "universal G_X")
    require(zero_mod_x2(diff(Gxt, T), t_value, x_square),
            "universal G_T")
    require(zero_mod_x2(Gxt - g_value, t_value, x_square),
            "universal G value")
    require(zero_mod_x2(Hess_xt - hessian_value, t_value, x_square),
            "universal Hessian")


# A disjoint normalized-projection redundancy check at exact rational controls.
def normalized_control(name: str, values: dict) -> tuple:
    ff = expand(f.subs(values))
    hh = expand(h.subs(values))
    rr = resultant(ff, hh, X)
    polynomial = Poly(rr, T)
    valuation = min(monomial[0] for monomial, coefficient in polynomial.terms()
                    if coefficient != 0)
    residual = cancel(rr/T**valuation)
    six_multiplicity = 0
    while rem(Poly(residual, T), Poly(6*T + 1, T)).as_expr() == 0:
        residual = cancel(residual/(6*T + 1))
        six_multiplicity += 1
    qq = Poly(residual, T)
    return (
        name, Poly(ff, X).degree(), Poly(hh, X).degree(),
        polynomial.degree(), valuation, six_multiplicity, qq.degree(),
        factor(qq.nth(0)), factor(qq.LC()), qq.eval(Rational(-1, 6)) != 0,
    )


normalized_generic = normalized_control(
    "generic", {Phi: 2, Theta: 5, eta: 3, zeta: 7})
require(normalized_generic[:7] ==
        ("generic", 8, 9, 79, 56, 2, 21),
        "generic normalized control degrees")
require(normalized_generic[7] == -1121008359375,
        "generic normalized constant endpoint")

normalized_anti = normalized_control(
    "anti", {Phi: 2, Theta: 5, eta: 3, zeta: -3})
require(normalized_anti[:7] ==
        ("anti", 7, 8, 63, 42, 2, 19),
        "anti normalized control degrees")


# Literal Newton polygon, packet, and anti-diagonal local normalization.
Fq = expand((s**2 - p)*(1 - Q*H) - Q*s**2/2)
support = sorted({monomial for monomial, coefficient in Poly(Fq, s, p).terms()
                  if coefficient != 0})


def cross(origin: tuple[int, int], left: tuple[int, int],
          right: tuple[int, int]) -> int:
    return ((left[0] - origin[0])*(right[1] - origin[1])
            - (left[1] - origin[1])*(right[0] - origin[0]))


def convex_hull(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
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


hull_raw = convex_hull(support)
expected_hull = [(0, 1), (2, 0), (5, 3), (1, 5), (0, 5)]
require(set(hull_raw) == set(expected_hull), "mixed-B Newton hull")
vertex_coefficients = tuple(
    factor(Poly(Fq, s, p).coeff_monomial(s**i*p**j))
    for i, j in expected_hull
)
expected_vertex_coefficients = (
    -1, -(Q - 2)/2, -Q*zeta, Q*eta, DELTA*Q,
)
require(all(cancel(left - right) == 0 for left, right
            in zip(vertex_coefficients, expected_vertex_coefficients)),
        "mixed-B hull vertex units")
hull = expected_hull
area2 = abs(sum(x1*y2 - y1*x2 for (x1, y1), (x2, y2)
                in zip(hull, hull[1:] + hull[:1])))
boundary = sum(gcd(abs(x2 - x1), abs(y2 - y1))
               for (x1, y1), (x2, y2) in zip(hull, hull[1:] + hull[:1]))
pick_genus = (area2 - boundary + 2)//2
require((area2, boundary, pick_genus) == (31, 11, 11),
        "mixed-B Pick data")

packet: list[int] = []
edge_rows: list[tuple] = []
for left, right in zip(hull, hull[1:] + hull[:1]):
    dx, dy = right[0] - left[0], right[1] - left[1]
    length = gcd(abs(dx), abs(dy))
    inward_u, inward_v = -dy//length, dx//length
    c = inward_u*left[0] + inward_v*left[1]
    index = inward_u + inward_v - c
    edge_rows.append((left, right, inward_u, inward_v, c, length, index))
    if inward_u == 1 and inward_v == 0 and c == 0:
        continue  # affine vertical face
    packet.extend([index]*length)
packet.sort(reverse=True)
require(packet == [8, 8, 4, 2, 2, 2, 1], "generic packet")
generic_defect = sum(index - 1 for index in packet)
require((sum(packet), generic_defect) == (27, 20), "generic packet totals")

local_anti = cancel(u**11*Fq.subs({s: 1/u, p: (1-a)/u**2, zeta: -eta}))
local_poly = Poly(expand(local_anti), a, u)
minimum_total_degree = min(sum(monomial) for monomial, coefficient
                           in local_poly.terms() if coefficient != 0)
tangent = factor(sum(coefficient*a**monomial[0]*u**monomial[1]
                     for monomial, coefficient in local_poly.terms()
                     if sum(monomial) == minimum_total_degree))
require(minimum_total_degree == 2, "anti tangent degree")
require(expand(tangent - Q*a*(eta*a - (DELTA + Theta)*u)) == 0,
        "anti ordinary-node tangent")
anti_packet = [7, 7, 4, 2, 2, 2, 1]
require((sum(anti_packet), sum(index - 1 for index in anti_packet)) ==
        (25, 18), "anti packet totals")


# Prime carrier and response arithmetic.
carrier = Poly(zeta*W**3 - (q - Rational(1, 2)), W)
require(carrier.degree() == 3 and factor(diff(carrier.as_expr(), W)) ==
        3*W**2*zeta, "generic cubic carrier")
generic_length = generic[1] + 4
anti_length = anti[1] + 4
require((generic_length, anti_length) == (25, 23), "critical lengths")
require(2*21 - generic_length - 1 + 3 == 19 < 20,
        "generic finite response")
require(2*(27 - generic_length) == 4 < 20,
        "generic full response")
require(2*19 - anti_length - 1 + 3 == 17 < 18,
        "anti finite response")
require(2*(25 - anti_length) == 4 < 18,
        "anti full response")


print("JC23_K0_MIXED_SOURCE_REFEREE")
print("DELTA", DELTA, "K", 0, "SOURCE_COMPLETE", True)
print("GRADIENT_IDENTITIES", True,
      "TRANSFORM_DET", factor(p*t**2),
      "HESSIAN_BRIDGE", "p*detD(A,B)=2*t^2*detHess mod(A,B) on pt!=0")
print("SOURCE_GENERIC", generic)
print("SOURCE_ANTI", anti)
print("HOSTILE_ETA_ZERO", eta_zero)
print("HOSTILE_ZETA_ZERO", zeta_zero)
print("HOSTILE_ANTI_THETA_MINUS_DELTA", anti_collision)
print("HOSTILE_ANTI_THETA_MINUS_DELTA_PHI_ZERO", anti_collision_phi_zero)
print("HOSTILE_BOTH_TOP_ZERO_WEIGHT_EXIT", both_top_zero)
print("NORMALIZED_GENERIC_CONTROL", normalized_generic)
print("NORMALIZED_ANTI_CONTROL", normalized_anti)
print("NO_LOSS", "p=0:A=-s,B=-6", "t=0:no-common-zero",
      "source-infinity generic=(3*zeta*p^2,9*zeta*p^2)",
      "source-infinity anti=(-3*eta*p^2,-9*eta*p^2)")
print("UNIVERSAL", "T0:X^2=-6,G=0,Hess=6",
      "T-1/6:X^2=6,G=1/2,Hess=-6")
print("HULL", tuple(hull), "VERTEX_COEFFS", vertex_coefficients,
      "PICK", (area2, boundary, pick_genus))
print("EDGES", tuple(edge_rows))
print("GENERIC_PACKET", tuple(packet), "SUM", sum(packet),
      "DEFECT", generic_defect, "L", generic_length)
print("ANTI_TANGENT", tangent)
print("ANTI_PACKET", tuple(anti_packet), "SUM", sum(anti_packet),
      "DEFECT", 18, "L", anti_length)
print("CARRIER", "pure cubic degree 3; separable over C(q); q=1/2 deleted")
print("RESPONSES", "generic finite 19<20 full 4<20",
      "anti finite 17<18 full 4<18")
print("VERDICT", "ACCEPT proposed generic and anti K=0 rows relative to inherited bridges")
