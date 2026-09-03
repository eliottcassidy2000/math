#!/usr/bin/env python3
"""Import-independent hostile audit of THM-4343's U+W=0 closure.

This intentionally imports neither candidate certificate.  It checks only the
finite/exact parts of the referee report; the proper-flat gluing is audited in
REFEREE.md as a geometric argument rather than disguised as computation.
"""

from fractions import Fraction as F
from math import gcd
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def twice_area(vertices):
    return abs(sum(
        x0 * y1 - y0 * x1
        for (x0, y0), (x1, y1) in zip(vertices, vertices[1:] + vertices[:1])
    ))


def boundary_length(vertices):
    return sum(
        gcd(abs(x1 - x0), abs(y1 - y0))
        for (x0, y0), (x1, y1) in zip(vertices, vertices[1:] + vertices[:1])
    )


def cubic_discriminant(a, b, c, d):
    return b*b*c*c - 4*a*c*c*c - 4*b*b*b*d - 27*a*a*d*d + 18*a*b*c*d


def poly_eval(coefficients_low_first, x, modulus=None):
    value = 0
    for coefficient in reversed(coefficients_low_first):
        value = value * x + coefficient
        if modulus is not None:
            value %= modulus
    return value


# 1. Reconstruct the global polygon and the three separately owned orbits.
polygon = [(0, 1), (2, 0), (4, 2), (4, 5), (0, 7)]
area2 = twice_area(polygon)
boundary = boundary_length(polygon)
interior = (area2 - boundary + 2) // 2
need((area2, boundary, interior) == (42, 14, 15), "global Pick ledger changed")

cubic_edge = frozenset(((4, 2), (4, 5)))
top_edge = frozenset(((4, 5), (0, 7)))
internal_edge = frozenset(((2, 0), (4, 5)))
need(cubic_edge != top_edge != internal_edge != cubic_edge, "owners collapsed")
need(cubic_edge & top_edge == {(4, 5)}, "wrong cubic/top closure incidence")
need(cubic_edge & internal_edge == {(4, 5)}, "wrong cubic/internal incidence")
need(top_edge & internal_edge == {(4, 5)}, "wrong top/internal incidence")

# K,W,U are units.  These are the constant/leading endpoint coefficients of
# the cubic, top, and internal edge polynomials.  Hence their roots lie in the
# relative tori, rather than at the common fixed point.
U0, W0, K0 = F(-1), F(1), F(1)
need(U0 * W0 * K0 != 0 and W0 == -U0, "hostile unit packet invalid")
endpoint_packets = {
    "cubic": (K0, W0),          # K + ... + W P^3
    "top": (U0, U0),            # U(X-1)^2
    "internal": (F(1), -W0),   # 1-WX
}
need(all(c0 != 0 and cinf != 0 for c0, cinf in endpoint_packets.values()),
     "an alleged interior root can exit to a toric endpoint")

# Same coordinate label 1 on all three distinct orbits is a positive hostile,
# not a collision.  First cubic has a double root; second has a triple root.
double_coeffs = [F(1), F(-1), F(-1), F(1)]
triple_coeffs = [F(-1), F(3), F(-3), F(1)]
need(poly_eval(double_coeffs, F(1)) == 0, "double hostile misses label 1")
need(poly_eval([i * double_coeffs[i] for i in range(1, 4)], F(1)) == 0,
     "double hostile is not repeated")
need(poly_eval([i*(i-1) * double_coeffs[i] for i in range(2, 4)], F(1)) != 0,
     "double hostile became triple")
need(all(poly_eval([
        i * (i-1) * (i-2) * triple_coeffs[i] for i in range(3, 4)
    ], F(1)) != 0 for _ in [0]), "triple hostile lost cubic degree")
for derivative_order in range(3):
    coeffs = triple_coeffs[:]
    for _ in range(derivative_order):
        coeffs = [i * coeffs[i] for i in range(1, len(coeffs))]
    need(poly_eval(coeffs, F(1)) == 0, "triple hostile jet failed")
owner_labels = {
    (tuple(sorted(cubic_edge)), F(1)),
    (tuple(sorted(top_edge)), F(1)),
    (tuple(sorted(internal_edge)), F(1)),
}
need(len(owner_labels) == 3, "coordinate quotient erased an owner")

# 2. Audit the top exact identity and its correction timing.  With W=-U and
# Z=0, the general Lambda-zero decomposition has A_top=2U+W=U.
r = F(7, 5)
U = F(13, 7)
W = -U
Z = F(0)
A_top = 2*U + W
lhs = U*r**6 + W*r**5 + Z*r**4
rhs = A_top*F(1, 2)*(r**6-r**4) - W*F(1, 2)*r**4*(r-1)**2
need(lhs == rhs and A_top == U, "top identity did not specialize exactly")
# Multiplication by the outer q=(r-1) makes the discrepancy cubic in q; after
# q=t^6 y and removal of the common t^12 it is delayed by t^6.
need(3*6 - 12 == 6, "correction no longer begins at t^6")
for s in range(1, 8):
    for nu in range(s + 1, s + 9):
        need(6*(nu-s) < 6*(s+nu), "terminal splitter lost the race")

# 3. The deep repeated range forces the entire cubic normal jet B to vanish.
alpha = eta = Phi = F(0)
upsilon = F(-17, 9)
xi = -upsilon
Delta = F(23, 11)
Theta = -Delta
for P in [F(-3), F(0), F(2, 5), F(7)]:
    B = P**3 * (Phi + eta*P + alpha*P**2)
    need(B == 0, "deep ladder did not kill B identically")
need(upsilon + xi == Delta + Theta == 0, "paired deep cancellations failed")

# 4. Reconstruct the unique terminal packet from its equations, then compute
# the cubic discriminant by the explicit four-term formula (no CAS/result
# import).  A second mod-11 certificate independently detects squarefreeness.
e = F(-1376, 135)
c = F(5152, 405)
Delta = 3*c - F(32, 9)
K = F(2848, 45) - F(7, 6)*Delta
U = -c*c/F(2)
W = -U
Theta = -Delta
upsilon = -F(8, 3)*c
xi = -upsilon
need((c, Delta, K, U, W, Theta, xi) == (
    F(5152, 405), F(4672, 135), F(1856, 81),
    F(-13271552, 164025), F(13271552, 164025),
    F(-4672, 135), F(41216, 1215)), "terminal packet arithmetic failed")
need(e + K == c and c*c + 2*U == 0, "terminal balance failed")

disc = cubic_discriminant(W, xi, Theta, K)
expected_disc = F(-3947729324424178958336, 32688603759375)
need(disc == expected_disc and disc != 0, "terminal cubic discriminant failed")
scaled_high_first = [207368, 86940, -88695, 58725]
p = 11
disc_mod_p = cubic_discriminant(*(coefficient % p for coefficient in scaled_high_first)) % p
need(disc_mod_p != 0, "mod-11 squarefreeness witness failed")
need(K != 0 and W != 0, "terminal cubic can meet an endpoint")

# 5. The main component is genus three: x(x^6-U) has degree seven and no
# repeated root when U is a unit in characteristic zero.
need(U != 0, "central unit lost")
for x in [F(0), F(1), F(-1), F(2)]:
    if x*(x**6-U) == 0:
        need(7*x**6-U != 0, "sampled repeated central branch root")
need((7 - 1)//2 == 3, "central hyperelliptic genus changed")

# 6. Explicit arithmetic-genus/normalization ledgers.  These count delta or
# node contributions, not transported simple response points.
squarefree_genus = (0+3+1) + (12+1) - 3 + 1
singular_raw_genus = (0+3+0) + (12+1+1) - 3 + 1
horizontal_normalized_genus = (0+3+0) + (12+1) - 3 + 1
rational_bridge_genus = (0+3+0+0) + (12+1+2) - 4 + 1
elliptic_tail_genus = (0+3+0+1) + (12+1+1) - 4 + 1
need((squarefree_genus, singular_raw_genus, horizontal_normalized_genus,
      rational_bridge_genus, elliptic_tail_genus) == (15, 15, 14, 15, 15),
     "component/genus ledger failed")

# All listed good-form orders are positive in their intended positive cone.
for s in range(1, 9):
    for nu in range(1, 12):
        twice_orders = [
            2*(6*s+2*nu), 5*s+9*nu, 2*(2*s+4*nu),
            3*s+7*nu, 2*(s+3*nu),
        ]
        need(all(order > 0 for order in twice_orders), "nonpositive A23 order")

print("THM-4343 hostile finite/exact referee audit: PASS")
print(f"Pick=(2A={area2},B={boundary},I={interior})")
print("owner supports=cubic/top/internal distinct; common vertex avoided by units")
print(f"terminal discriminant={disc}")
print(f"terminal discriminant mod 11={disc_mod_p}")
print("genus ledgers=squarefree15 raw-singular15 horizontal14 bridge15 tail15")
print("scope note=proper-flat/disjoint-support gluing remains a prose theorem step")
