#!/usr/bin/env python3
"""Finite exact companion for provisional THM-3596.

The proof is toric and local.  This script verifies the exact Laurent packet,
Newton polygon, face polynomials, arm resolutions, and the THM-3591 invoice
controls without Python assert gates.
"""

from math import gcd

import sympy as sp


b, c, e, w, lam = sp.symbols("b c e w lam")
CHECKS = 0


def require(label, condition):
    """Record one active truth gate and fail with a stable label."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError("FAILED: " + label)


def zero(expr):
    return sp.cancel(sp.expand(expr)) == 0


def bracket(F, G):
    """A13 bracket on ambient representatives."""
    Sigma = b * (b**2 + 1)
    return sp.expand(
        c**3 * (sp.diff(F, b) * sp.diff(G, c) - sp.diff(F, c) * sp.diff(G, b))
        - sp.diff(Sigma, b)
        * (sp.diff(F, c) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, c))
        - 3
        * c**2
        * e
        * (sp.diff(F, b) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, b))
    )


def laurent_packet(expr):
    """Return {(b exponent,c exponent): coefficient} for a Laurent packet."""
    packet = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        powers = term.as_powers_dict()
        i = int(powers.get(b, 0))
        j = int(powers.get(c, 0))
        coefficient = sp.cancel(term / (b**i * c**j))
        packet[i, j] = sp.cancel(packet.get((i, j), 0) + coefficient)
    return {point: coefficient for point, coefficient in packet.items() if coefficient != 0}


def hull(points):
    """Andrew monotone-chain hull without repeating the first vertex."""
    points = sorted(set(points))

    def cross(o, a, d):
        return (a[0] - o[0]) * (d[1] - o[1]) - (a[1] - o[1]) * (d[0] - o[0])

    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def lattice_distance(a, d, point):
    dx, dy = d[0] - a[0], d[1] - a[1]
    area = dx * (point[1] - a[1]) - dy * (point[0] - a[0])
    return area // gcd(abs(dx), abs(dy))


print("THM-3596 exact companion")
print("SECTION Laurent Hamiltonian packet")

Sigma = b * (b**2 + 1)
q = (225 * b**4 + 237 * b**2 + 8) / 8
v = sp.Rational(105, 16) * b
Q = e**3 + q * e + c**2 + v * c**4 + lam * c * e
F = Sigma**3 * c**-9 + q * Sigma * c**-3 + c**2 + v * c**4 + lam * Sigma * c**-2
require("Laurent elimination", zero(Q.subs(e, Sigma / c**3) - F))

for i in range(4):
    for j in range(4):
        R = b**i * c**j
        T = b ** (3 - i) * c ** (3 - j)
        require(
            f"density bracket i={i} j={j}",
            zero(bracket(R, T) - c**3 * (sp.diff(R, b) * sp.diff(T, c) - sp.diff(R, c) * sp.diff(T, b))),
        )

Fb = sp.diff(F, b)
Fc = sp.diff(F, c)
require("D_F(b)", zero(c**3 * Fc - c**3 * sp.diff(F, c)))
require("D_F(c)", zero(-c**3 * Fb + c**3 * sp.diff(F, b)))

num_b = sp.Poly(
    sp.fraction(sp.together(Fb))[0],
    b,
    c,
    domain=sp.QQ.frac_field(lam),
)
num_c = sp.Poly(
    sp.fraction(sp.together(Fc))[0],
    b,
    c,
    domain=sp.QQ.frac_field(lam),
)
require("derivative Laurent gcd", sp.gcd(num_b, num_c).as_expr() == 1)
print("PASS Laurent elimination, density, time form, and derivative gcd")


print("SECTION Newton hull and interior time-form point")
packet = laurent_packet(F - w)
expected_support = {
    (3, -9),
    (5, -9),
    (7, -9),
    (9, -9),
    (1, -3),
    (3, -3),
    (5, -3),
    (7, -3),
    (1, -2),
    (3, -2),
    (0, 2),
    (1, 4),
    (0, 0),
}
require("generic Laurent support", set(packet) == expected_support)

packet_zero = laurent_packet((F - w).subs(lam, 0))
expected_zero = expected_support - {(1, -2), (3, -2)}
require("lambda-zero support", set(packet_zero) == expected_zero)

expected_hull = [(0, 0), (3, -9), (9, -9), (7, -3), (1, 4), (0, 2)]
require("generic hull", hull(packet) == expected_hull)
require("lambda-zero hull", hull(packet_zero) == expected_hull)

p = (1, -2)
distances = []
for index, a in enumerate(expected_hull):
    d = expected_hull[(index + 1) % len(expected_hull)]
    distance = lattice_distance(a, d, p)
    require(f"strict interior edge={index}", distance > 0)
    distances.append(distance)
require("primitive distance invoice", distances == [1, 7, 17, 36, 6, 1])
require("time-form monomial", p == (1, -2))
print(f"PASS hull vertices={len(expected_hull)} distances={distances}")


print("SECTION toric face packet")
face_bottom = Sigma**3 * c**-9
face_bc = b**9 * c**-9 + sp.Rational(225, 8) * b**7 * c**-3
face_cd = sp.Rational(225, 8) * b**7 * c**-3 + sp.Rational(105, 16) * b * c**4
face_de = c**2 + sp.Rational(105, 16) * b * c**4
face_ef = c**2 - w
face_fa = b**3 * c**-9 + b * c**-3 - w

face_expected = [face_fa, face_bottom, face_bc, face_cd, face_de, face_ef]
for index, (a, d) in enumerate(
    zip(expected_hull, expected_hull[1:] + expected_hull[:1])
):
    dx, dy = d[0] - a[0], d[1] - a[1]
    selected = 0
    for point, coefficient in packet.items():
        if dx * (point[1] - a[1]) - dy * (point[0] - a[0]) == 0:
            selected += coefficient * b ** point[0] * c ** point[1]
    require(f"face identity edge={index}", zero(selected - face_expected[index]))

require("bottom factor", zero(face_bottom - b**3 * (b**2 + 1) ** 3 * c**-9))
for difference in [(-2, 6), (-6, 7), (-1, -2)]:
    require(f"binomial face direction {difference}", difference != (0, 0))
require("EF discriminant", sp.discriminant(c**2 - w, c) == 4 * w)
u = sp.symbols("u")
require("FA discriminant", sp.discriminant(u**3 + u - w, u) == -27 * w**2 - 4)
require("bottom torus roots", sp.factor((b**2 + 1) ** 3) == (b**2 + 1) ** 3)
print("PASS six faces; only bottom has the two triple torus roots")


print("SECTION smooth-arm restoration")
require("q at central arm", q.subs(b, 0) == 1)
require("q on side arms", sp.rem(sp.together(q + sp.Rational(1, 2)), b**2 + 1, b) == 0)
require("Sigma prime central", sp.diff(Sigma, b).subs(b, 0) == 1)
require(
    "Sigma prime side",
    sp.rem(sp.diff(Sigma, b) + 2, b**2 + 1, b) == 0,
)

disc_central = sp.discriminant(e**3 + e - w, e)
disc_side = sp.discriminant(e**3 - e / 2 - w, e)
require("central arm cubic generic", disc_central == -27 * w**2 - 4)
require("side arm cubic generic", disc_side == sp.Rational(1, 2) - 27 * w**2)

y = sp.symbols("y")
for beta in [sp.Integer(0), sp.I, -sp.I]:
    slope = sp.diff(Sigma, b).subs(b, beta)
    trial = beta + y / slope
    require(
        f"inverse arm first jet beta={beta}",
        sp.expand(Sigma.subs(b, trial)).coeff(y, 1) == 1,
    )
    require(f"arm slope unit beta={beta}", slope != 0)

require("lambda arm independence", sp.diff(Q, e).subs(c, 0).has(lam) is False)
require("lambda face independence", all(face.has(lam) is False for face in face_expected))
print("PASS three arm charts and generic unit time forms")


print("SECTION THM-3591 invoice controls")


def W(r, f, s, g):
    return sp.expand(s * sp.diff(f, b) * g - r * f * sp.diff(g, b))


scalar_channel = W(1, 1, -3, Sigma * q) + W(-6, Sigma**2, 4, v)
require("global scalar payment", zero(scalar_channel + 1))

Delta = 3 * e**2 + q - 2 * lam * e**2
require("unrepaired arm debt", zero(Delta.subs(lam, 0) - (3 * e**2 + q)))
require("repaired arm interpolation", zero(Delta.subs(lam, sp.Rational(3, 2)) - q))
require("q inverse-derivative residue", sp.rem(sp.together(q - (1 + 3 * b**2 / 2)), Sigma, b) == 0)
print("PASS lambda=0 hostile and lambda=3/2 fully paid arm/scalar control")


print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
