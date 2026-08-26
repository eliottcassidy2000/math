#!/usr/bin/env python3
"""Primary exact certificate for THM-4183, the P-only Delta-zero wall.

The calculation stays in normalized (X,T) coordinates for the critical
scheme.  It also rebuilds both Newton polygons directly in (s,p), including
the Theta=0 blowdown where the index-three and index-five edges are replaced
by one primitive index-seven edge.  No residual discriminant or selected
projection-fibre separation is assumed.
"""

from math import gcd

import sympy as sp


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def valuation(poly, variable):
    terms = sp.Poly(poly, variable).terms()
    need(bool(terms), "zero polynomial has no valuation")
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


def exact_quotient(numerator, denominator, variable, message):
    quotient = sp.cancel(numerator / denominator)
    need(sp.denom(quotient) == 1, message + " quotient is not polynomial")
    return sp.Poly(quotient, variable)


def convex_hull(points):
    points = sorted(set(points))

    def cross(origin, left, right):
        return ((left[0] - origin[0]) * (right[1] - origin[1])
                - (left[1] - origin[1]) * (right[0] - origin[0]))

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
    return tuple(lower[:-1] + upper[:-1])


def polygon_ledger(points):
    vertices = convex_hull(points)
    area2 = abs(sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
        for index in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(index + 1) % len(vertices)][0] - vertices[index][0]),
            abs(vertices[(index + 1) % len(vertices)][1] - vertices[index][1]))
        for index in range(len(vertices))
    )
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return vertices, area2, boundary, (area2 - boundary + 2) // 2


def selected_face(poly, x, y, u, v, level):
    answer = 0
    for (i, j), coefficient in sp.Poly(poly, x, y).terms():
        if u * i + v * j == level:
            answer += coefficient * x**i * y**j
    return sp.factor(answer)


def reduced_mod(expr, modulus, variable):
    numerator, denominator = sp.together(sp.cancel(expr)).as_numer_denom()
    need(sp.resultant(denominator, modulus, variable) != 0,
         "denominator meets reduction modulus")
    return sp.rem(sp.Poly(numerator, variable),
                  sp.Poly(modulus, variable)).as_expr()


X, T = sp.symbols("X T")
Phi, Theta, eta = sp.symbols("Phi Theta eta")
P = T + X**2 * T**2
Y = X * T * P
K0 = sp.Rational(2848, 45)
epsilon = -sp.Rational(1376, 135)

# Complete P-only, Delta-zero exact-M=9 source: zeta=Delta=0, eta!=0.
G = sp.expand(
    -X**2 * T / 2 - 3 * P + sp.Rational(8, 3) * P**2
    + epsilon * P**3 + K0 * Y**2 + Phi * P**2 * Y
    + Theta * P * Y**2 + eta * P**3 * Y
)
f = sp.expand(sp.cancel(sp.diff(G, X) / T))
h = sp.expand(sp.diff(G, T))
Hessian = sp.factor(sp.det(sp.hessian(G, (X, T))))
critical_jacobian = sp.factor(sp.det(sp.Matrix((
    (sp.diff(f, X), sp.diff(f, T)),
    (sp.diff(h, X), sp.diff(h, T)),
))))

need(sp.factor(T * critical_jacobian - Hessian
               - f * sp.diff(G, X, T)) == 0,
     "Morse-resultant bridge changed")
need((sp.degree(f, X), sp.degree(h, X)) == (8, 9),
     "normalized critical degrees changed")
need(sp.factor(sp.Poly(f, X).LC() - 9 * eta * T**8) == 0
     and sp.factor(sp.Poly(h, X).LC() - 9 * eta * T**8) == 0,
     "normalized X-infinity sidecar changed")
need(sp.factor(f.subs(T, 0) + X) == 0,
     "universal T-zero first equation changed")
need(sp.factor(h.subs(T, 0) + (X**2 + 6) / 2) == 0,
     "universal T-zero second equation changed")
need(reduced_mod(Hessian.subs(T, 0) - 6, X**2 + 6, X) == 0,
     "T-zero universal pair is not Morse")

for expression, expected, label in (
        (f, 0, "f"), (h, 0, "h"), (G, sp.Rational(1, 2), "value")):
    remainder = reduced_mod(
        expression.subs(T, -sp.Rational(1, 6)) - expected,
        X**2 - 6,
        X,
    )
    need(sp.factor(remainder) == 0, "T=-1/6 " + label + " changed")
need(reduced_mod(Hessian.subs(T, -sp.Rational(1, 6)) + 6,
                 X**2 - 6, X) == 0,
     "T=-1/6 universal pair is not Morse")

# The exact symbolic resultant works simultaneously for arbitrary Phi and
# both exhaustive Theta rows.
resultant = sp.resultant(f, h, X)
need(valuation(resultant, T) == 56, "normalized T-artifact changed")
Q20 = exact_quotient(resultant, T**56 * (6 * T + 1)**2,
                     T, "generic normalized resultant")
need(Q20.degree() == 20, "Theta-nonzero residual degree changed")
need(sp.factor(Q20.LC() - 72900 * K0**4 * eta**5 * Theta**4) == 0,
     "Theta-nonzero top endpoint changed")
need(sp.factor(Q20.TC() + sp.Rational(3**15, 2**7) * eta**7) == 0,
     "Theta-nonzero bottom endpoint changed")

f0 = sp.expand(f.subs(Theta, 0))
h0 = sp.expand(h.subs(Theta, 0))
resultant0 = sp.resultant(f0, h0, X)
need(sp.factor(resultant0 - resultant.subs(Theta, 0)) == 0,
     "Theta-zero direct specialization changed")
Q19 = exact_quotient(resultant0, T**56 * (6 * T + 1)**2,
                     T, "Theta-zero normalized resultant")
need(Q19.degree() == 19, "Theta-zero residual degree changed")
need(sp.factor(Q19.LC() - 944784 * K0**6 * eta**7) == 0,
     "Theta-zero top endpoint changed")
need(sp.factor(Q19.TC() + sp.Rational(3**15, 2**7) * eta**7) == 0,
     "Theta-zero bottom endpoint changed")

# First normal coefficient at the Theta wall.  In reciprocal u=1/T, the
# nominal degree-twenty endpoint has order Theta^4, while the u coefficient
# is a unit on Theta=0.  Exactly one projected root escapes.
need(sp.factor(Q20.nth(19).subs(Theta, 0)
               - 944784 * K0**6 * eta**7) == 0,
     "first reciprocal normal coefficient changed")
u_escape = -sp.cancel(
    (72900 * K0**4 * eta**5 * Theta**4)
    / (944784 * K0**6 * eta**7)
)
need(sp.factor(u_escape + sp.Rational(25, 324)
               * Theta**4 / (K0**2 * eta**2)) == 0,
     "reciprocal escape scale changed")

# Direct Delta-zero boundary reconstruction in (s,p).
s, p, Q = sp.symbols("s p Q")
H = sp.expand(
    -3 * p + sp.Rational(8, 3) * p**2 + epsilon * p**3
    + K0 * s**2 * p**2 + Phi * s * p**3
    + Theta * s**2 * p**3 + eta * s * p**4
)
F = sp.expand((s**2 - p) * (1 - Q * H) - Q * s**2 / 2)

generic_vertices = ((0, 1), (2, 0), (4, 2), (4, 3),
                    (3, 4), (1, 5), (0, 4))
wall_vertices = ((0, 1), (2, 0), (4, 2),
                 (3, 4), (1, 5), (0, 4))
for specialization, expected in (
        ({Phi: 0, Theta: 1, eta: 1},
         (generic_vertices, 28, 10, 10)),
        ({Phi: 0, Theta: 0, eta: 1},
         (wall_vertices, 27, 9, 10))):
    support = tuple(
        powers for powers, coefficient
        in sp.Poly(F.subs(specialization), s, p).terms()
        if coefficient != 0
    )
    need(polygon_ledger(support) == expected,
         "Delta-zero P-only polygon changed")

low_face = s**2 * (1 - Q / 2) - p
carrier_face = s**2 * ((1 - Q / 2) - K0 * Q * (s * p)**2)
theta_face = -Q * p**2 * s**4 * (K0 + Theta * p)
mixed_face = -Q * p**3 * s**3 * (eta * p + Theta * s)
top_face = Q * eta * p**4 * s * (p - s**2)
diagonal_face = Q * p**4 * (epsilon + eta * s * p)
vertical_face = p * (-1 + Q * (-3 * p + sp.Rational(8, 3) * p**2
                               + epsilon * p**3))
for normal, expected, label in (
        ((1, 2, 2), low_face, "low"),
        ((-1, 1, -2), carrier_face, "quadratic carrier"),
        ((-1, 0, -4), theta_face, "Theta edge"),
        ((-1, -1, -7), mixed_face, "mixed edge"),
        ((-1, -2, -11), top_face, "top edge"),
        ((1, -1, -4), diagonal_face, "diagonal edge"),
        ((1, 0, 0), vertical_face, "vertical affine edge")):
    need(sp.factor(selected_face(F, s, p, *normal) - expected) == 0,
         label + " changed")

F0 = sp.expand(F.subs(Theta, 0))
merged_face = -Q * p**2 * s**3 * (K0 * s + eta * p**2)
need(sp.factor(selected_face(F0, s, p, -2, -1, -10)
               - merged_face) == 0,
     "Theta-zero primitive merged edge changed")
need(sp.Poly(K0 * s + eta * p**2, s).degree() == 1,
     "Theta-zero merged edge is not torus-smooth")

# Packets, Pick genus, response degrees, and exact strict contradictions.
generic_packet = (8, 5, 4, 3, 2, 2, 1)
wall_packet = (8, 7, 4, 2, 2, 1)
for packet, genus, full_n, finite_n, length in (
        (generic_packet, 10, 25, 21, 24),
        (wall_packet, 10, 24, 20, 23)):
    defect = sum(index - 1 for index in packet)
    need(sum(packet) == full_n and defect == 2 * genus - 2,
         "packet/Pick ledger changed")
    need(full_n - 4 == finite_n,
         "quadratic-carrier response degree changed")
    need(2 * (full_n - length) < defect,
         "full response contradiction changed")
    need(2 * finite_n - length - 1 + 2 < finite_n - 1,
         "finite response contradiction changed")

print("P_ONLY_DELTA_ZERO_PRIMARY_EXACT_CERTIFICATE")
print("wall=zeta:0,Delta:0,eta!=0;Phi_arbitrary")
print("rows=Theta_nonzero_or_Theta_zero")
print("bridge=T*detD(f,h)-detHess(G)-f*G_XT=0")
print("critical_degrees=Theta_nonzero:20;Theta_zero:19")
print("critical_lengths=Theta_nonzero:24;Theta_zero:23")
print("normal_escape=u=1/T~-25*Theta^4/(324*K0^2*eta^2)")
print("packets=Theta_nonzero:(8,5,4,3,2,2,1);Theta_zero:(8,7,4,2,2,1)")
print("carrier=K0*(s*p)^2=q-1/2;degree=2;beta=2")
print("responses=full_and_finite_bounds_strict_in_both_rows")
print(f"checks={CHECKS}")
print("verdict=P_ONLY_DELTA_ZERO_CLOSES")
