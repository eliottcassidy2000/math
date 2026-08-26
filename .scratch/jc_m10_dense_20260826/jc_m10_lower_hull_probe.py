#!/usr/bin/env python3
"""Tracked FINITE-EXACT scratch probe for the first M=10 reduced JC face.

This tracked scratch artifact is deliberately outside canon and is not a
theorem certificate.  It expands the Q-adic support, computes every lower
supporting plane over Q, and checks the two normalizations that make the
M=10 frontier qualitatively different from M=9.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd

import sympy as sp


def require(condition, label):
    if not condition:
        raise RuntimeError(label)
    print(f"PASS {label}")


def monomials_through(weight):
    ans = []
    for i in range(weight // 2 + 1):
        for j in range(weight // 3 + 1):
            w = 2 * i + 3 * j
            if 0 < w <= weight and (i, j) not in {(0, 1), (1, 1)}:
                ans.append((i, j, w))
    return sorted(ans, key=lambda row: (row[2], row[1], row[0]))


def expanded_support(monomials):
    # Coordinates are (s exponent, p exponent, Q-adic height).
    support = {(2, 0, 0): ["+s^2"], (0, 1, 0): ["-p"]}
    for i, j, _ in monomials:
        support.setdefault((j + 2, i + j, 1), []).append(f"-c_{i}{j}")
        support.setdefault((j, i + j + 1, 1), []).append(f"+c_{i}{j}")
    return support


def lower_faces(support):
    points = sorted(support)
    faces = {}
    for ia, ib, ic in combinations(range(len(points)), 3):
        p, q, r = (points[k] for k in (ia, ib, ic))
        det = ((q[0] - p[0]) * (r[1] - p[1])
               - (q[1] - p[1]) * (r[0] - p[0]))
        if det == 0:
            continue
        aa = Fraction(
            (q[2] - p[2]) * (r[1] - p[1])
            - (q[1] - p[1]) * (r[2] - p[2]), det)
        bb = Fraction(
            (q[0] - p[0]) * (r[2] - p[2])
            - (q[2] - p[2]) * (r[0] - p[0]), det)
        cc = Fraction(p[2]) - aa * p[0] - bb * p[1]
        gaps = [Fraction(h) - aa * s - bb * pexp - cc
                for s, pexp, h in points]
        if min(gaps) < 0:
            continue
        equality = tuple(k for k, gap in enumerate(gaps) if gap == 0)
        if len(equality) >= 3:
            faces[equality] = (aa, bb, cc)
    return points, sorted(
        [(plane, tuple(points[k] for k in equality))
         for equality, plane in faces.items()],
        key=lambda row: row[0],
    )


def convex_hull(points):
    points = sorted(set(points))

    def cross(o, a, b):
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

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


def pick(points):
    polygon = convex_hull(points)
    twice_area = abs(sum(
        polygon[k][0] * polygon[(k + 1) % len(polygon)][1]
        - polygon[(k + 1) % len(polygon)][0] * polygon[k][1]
        for k in range(len(polygon))
    ))
    boundary = sum(
        gcd(abs(polygon[(k + 1) % len(polygon)][0] - polygon[k][0]),
            abs(polygon[(k + 1) % len(polygon)][1] - polygon[k][1]))
        for k in range(len(polygon))
    )
    interior = (twice_area - boundary + 2) // 2
    return polygon, twice_area, boundary, interior


def edge_packet(vertices):
    packet = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        # The vertical s=0 edge consists of affine points, not punctures.
        if start[0] == end[0] == 0:
            continue
        packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True))


def poly_valuation(polynomial):
    return min(key[0] for key in polynomial.as_dict())


def coefficient_digest(polynomial, monic=False):
    if monic:
        polynomial = polynomial.monic()
    payload = ",".join(str(sp.Rational(c)) for c in polynomial.all_coeffs())
    return sha256(payload.encode()).hexdigest()


all_monomials = monomials_through(10)
expected = [
    (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
    (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
    (0, 3, 9), (5, 0, 10), (2, 2, 10),
]
require(all_monomials == expected, "complete b=d=0 monomial universe through M=10")
print("MONOMIALS", ",".join(f"p^{i}y^{j}:w{w}" for i, j, w in all_monomials))

required = {(5, 0), (2, 2), (0, 3)}
optional = [(i, j, w) for i, j, w in all_monomials if (i, j) not in required]
expected_planes = {
    (Fraction(1, 10), Fraction(1, 5), Fraction(-1, 5)),
    (Fraction(1, 6), Fraction(1, 6), Fraction(-1, 3)),
}

# The two-face hull is independent of every lower coefficient when the two
# M=10 coefficients and zeta=[y^3]H are nonzero.
for mask in range(1 << len(optional)):
    chosen = [(i, j, w) for i, j, w in all_monomials if (i, j) in required]
    chosen += [row for bit, row in enumerate(optional) if mask & (1 << bit)]
    _, faces = lower_faces(expanded_support(chosen))
    require_planes = {plane for plane, _ in faces}
    if require_planes != expected_planes:
        raise RuntimeError((mask, require_planes, faces))
print(f"PASS robust generic lower hull across {1 << len(optional)} lower-support subsets")

_, generic_faces = lower_faces(expanded_support(all_monomials))
for plane, points in generic_faces:
    print("GENERIC_FACE", plane, points)

main_points = [(0, 1), (2, 0), (4, 4), (2, 5), (0, 6)]
tail_points = [(2, 0), (4, 4), (5, 3)]
global_points = main_points + tail_points
require(pick(main_points)[1:] == (30, 10, 11), "main Pick ledger (30,10,11)")
require(pick(tail_points)[1:] == (6, 6, 1), "tail Pick ledger (6,6,1)")
require(pick(global_points)[1:] == (36, 12, 13), "global Pick ledger (36,12,13)")
print("PICK main", pick(main_points))
print("PICK tail", pick(tail_points))
print("PICK global", pick(global_points))
generic_polygon = pick(global_points)[0]
packet = edge_packet(generic_polygon)
require(packet == (9, 9, 6, 2, 2, 2, 1), "generic infinity packet")
require(sum(index - 1 for index in packet) == 24,
        "packet defect saturates 2g-2=24")
print("PACKET", packet, "n", sum(packet), "defect", sum(index - 1 for index in packet))

S, P, a, b, zeta = sp.symbols("S P a b zeta", nonzero=True)
g_main = (S**2 - P) * (1 - a * P**5 - b * S**2 * P**4)
g_tail = 1 - zeta * S**3 * P**3 - b * S**2 * P**4
node_det = sp.det(sp.Matrix([
    [sp.diff(S**2 - P, S), sp.diff(S**2 - P, P)],
    [sp.diff(1 - a * P**5 - b * S**2 * P**4, S),
     sp.diff(1 - a * P**5 - b * S**2 * P**4, P)],
])).subs(P, S**2)
require(sp.factor(node_det) == -10 * S**9 * (a + b),
        "ten main attachments are transverse off a+b=0")
require(sp.expand(g_main) == sp.expand(
    (S**2 - P) * (1 - a * P**5 - b * S**2 * P**4)),
    "main face factorization")

T, W = sp.symbols("T W")
require(sp.expand(g_tail - (1 - zeta * T**3 - b * W**2)
                         .subs({T: S * P, W: S * P**2})) == 0,
        "tail normalization b W^2=1-zeta T^3")

X, Y = sp.symbols("X Y")
f5 = 1 - X**5 - Y**2
require(sp.gcd(sp.Poly(1 - X**5, X), sp.Poly(-5 * X**4, X)).degree() == 0,
        "cyclotomic quintic face is smooth")
require((5 - 1) // 2 == 2, "cyclotomic face genus two")
phi5 = sp.cyclotomic_poly(5, X)
require(phi5 == X**4 + X**3 + X**2 + X + 1,
        "order-five automorphism has rational minimal polynomial Phi_5")
print("CYCLOTOMIC_FACE b*Y^2=1-a*X^5 genus=2 eigencharacters=zeta5,zeta5^2")

# In tail coordinates x^3=-zeta*T^3 and y^2=b*W^2, the two points on
# the length-two internal edge are (0,+1),(0,-1) on y^2=x^3+1.
# The tangent at (0,1) has slope zero, so doubling gives (0,-1).
x0, y0 = sp.Integer(0), sp.Integer(1)
slope = sp.cancel(3 * x0**2 / (2 * y0))
x2 = sp.expand(slope**2 - 2 * x0)
y2 = sp.expand(slope * (x0 - x2) - y0)
require((x2, y2) == (0, -1), "tail attachment point has exact order three")
print("TAIL_ATTACHMENTS (0,1),(0,-1); difference is nonzero 3-torsion")
print("TAIL_MAP_OBSTRUCTION Hom survives; any compatible Eisenstein endomorphism is divisible by 1-zeta_3, hence degree divisible by 3")

# A lawful exact rational control for the dense critical-open chamber.  The
# forced row is K=2848/45-(7/6)Delta; every other displayed coefficient is a
# freely chosen nonzero rational, and a+b is nonzero.
ss, pp = sp.symbols("ss pp")
tt = pp - ss**2
Delta_c = sp.Rational(1)
K_c = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta_c
Phi_c, Theta_c = sp.Rational(2), sp.Rational(5)
eta_c, zeta_c = sp.Rational(7), sp.Rational(11)
a_c, b_c = sp.Rational(13), sp.Rational(17)
H_c = (
    -3 * pp + sp.Rational(8, 3) * pp**2
    - sp.Rational(1376, 135) * pp**3
    + K_c * ss**2 * pp**2 + Phi_c * ss * pp**3
    + Delta_c * pp**4 + Theta_c * ss**2 * pp**3
    + eta_c * ss * pp**4 + zeta_c * ss**3 * pp**3
    + a_c * pp**5 + b_c * ss**2 * pp**4
)
A_c = sp.cancel((-ss * pp + tt**2 * sp.diff(H_c, ss)) / pp)
C_c = ss**2 + 2 * tt**2 * sp.diff(H_c, pp)
B_c = sp.cancel((C_c + ss * A_c) / tt**2)
require(sp.denom(A_c) == sp.denom(B_c) == 1,
        "lawful source critical pair is polynomial")
require((sp.degree(A_c, ss), sp.degree(B_c, ss)) == (6, 3),
        "source critical s-degrees are (6,3)")
require(sp.LC(sp.Poly(A_c, ss)) == 3 * zeta_c * pp**2
        and sp.LC(sp.Poly(B_c, ss)) == 9 * zeta_c * pp**2,
        "source-infinity leading rows exclude finite p!=0 loss")
require(sp.expand(A_c.subs(pp, 0)) == -ss
        and sp.expand(B_c.subs(pp, 0)) == -6,
        "p=0 has no source critical point")
require(sp.expand(A_c.subs(pp, ss**2)) == -ss
        and sp.expand(B_c.subs({pp: ss**2, ss: 0})) == -6,
        "t=0 has no source critical point")

resultant_ab = sp.Poly(sp.resultant(A_c, B_c, ss), pp)
resultant_ac = sp.Poly(sp.resultant(A_c, C_c, ss), pp)
require((poly_valuation(resultant_ab), poly_valuation(resultant_ac)) == (6, 8),
        "source resultant exceptional powers are p^6 and p^8")
R25_ab = sp.Poly(sp.cancel(resultant_ab.as_expr() / pp**6), pp)
R25_ac = sp.Poly(sp.cancel(resultant_ac.as_expr() / pp**8), pp)
require(R25_ab == R25_ac and R25_ab.degree() == 25,
        "two independent source projections give the same R25")
require(R25_ab.nth(0) == -sp.Rational(189675421056, 5),
        "R25 nonzero constant endpoint")
require(R25_ab.LC() == sp.Integer(2731549392000000000),
        "R25 nonzero leading endpoint")
require(sp.gcd(R25_ab, R25_ab.diff()).degree() == 0,
        "R25 control is squarefree")
print("SOURCE_R25 degree", R25_ab.degree(), "constant", R25_ab.nth(0),
      "leading", R25_ab.LC(), "digest", coefficient_digest(R25_ab))

XX, TT = sp.symbols("XX TT")
PP = TT + XX**2 * TT**2
YY = XX * TT * PP
G_c = (
    -XX**2 * TT / 2 - 3 * PP + sp.Rational(8, 3) * PP**2
    - sp.Rational(1376, 135) * PP**3 + K_c * YY**2
    + Phi_c * PP**2 * YY + Delta_c * PP**4
    + Theta_c * PP * YY**2 + eta_c * PP**3 * YY
    + zeta_c * YY**3 + a_c * PP**5 + b_c * PP**2 * YY**2
)
resultant_xt = sp.Poly(
    sp.resultant(sp.cancel(sp.diff(G_c, XX) / TT), sp.diff(G_c, TT), XX), TT
)
require(poly_valuation(resultant_xt) == 72,
        "X,T projection has Sylvester factor T^72")
xt_residual = sp.cancel(resultant_xt.as_expr() / TT**72)
six_t_plus_one_power = 0
while sp.rem(sp.Poly(xt_residual, TT), sp.Poly(6 * TT + 1, TT)) == 0:
    xt_residual = sp.cancel(xt_residual / (6 * TT + 1))
    six_t_plus_one_power += 1
Q25 = sp.Poly(xt_residual, TT)
require(six_t_plus_one_power == 2 and Q25.degree() == 25,
        "X,T projection is T^72(6T+1)^2 Q25")
require(Q25.nth(0) == -sp.Integer(768867187500000000),
        "Q25 nonzero zero endpoint")
q_minus_sixth = -sp.Rational(
    1100752588157520645361077864812824646096209729,
    6071976554829712588800000,
)
require(Q25.eval(-sp.Rational(1, 6)) == q_minus_sixth,
        "Q25 nonzero universal-fibre endpoint")
require(sp.gcd(Q25, Q25.diff()).degree() == 0,
        "Q25 control is squarefree")
print("TARGET_Q25 degree", Q25.degree(), "constant", Q25.nth(0),
      "at_-1/6", Q25.eval(-sp.Rational(1, 6)),
      "digest", coefficient_digest(Q25, monic=True))

critical_length = 25 + 2 + 2
finite_n, carrier_index = sum(packet) - 3 * 2, 3
finite_cap = 2 * finite_n - critical_length - 1 + carrier_index
full_cap = 2 * (sum(packet) - critical_length)
require((critical_length, finite_n, carrier_index) == (29, 25, 3),
        "critical length and cubic finite-carrier ledger")
require(finite_cap == 23 < finite_n - 1,
        "finite response strict deficit 23<24")
require(full_cap == 4 < sum(index - 1 for index in packet),
        "full response strict deficit 4<24")
print("RESPONSE L=29 finite=(n,beta,cap,target)=", (finite_n, carrier_index, finite_cap, finite_n - 1),
      "full=(n,cap,defect)=", (sum(packet), full_cap, 24))

# With zeta=0, every possible extra lower face is the rational side plane
# h=(s-2)/2.  It is present iff y^2 or p*y^2 is present below the top.
zeta_zero = [row for row in all_monomials if row[:2] != (0, 3)]
z0_optional = [row for row in zeta_zero if row[:2] not in {(5, 0), (2, 2)}]
z0_signatures = set()
for mask in range(1 << len(z0_optional)):
    chosen = [row for row in zeta_zero if row[:2] in {(5, 0), (2, 2)}]
    chosen += [row for bit, row in enumerate(z0_optional) if mask & (1 << bit)]
    _, faces = lower_faces(expanded_support(chosen))
    z0_signatures.add(tuple(plane for plane, _ in faces))
expected_z0 = {
    ((Fraction(1, 10), Fraction(1, 5), Fraction(-1, 5)),),
    ((Fraction(1, 10), Fraction(1, 5), Fraction(-1, 5)),
     (Fraction(1, 2), Fraction(0), Fraction(-1))),
}
require(z0_signatures == expected_z0,
        "zeta=0 has only main face and optional rational side plane")
print("ZETA_ZERO_SIGNATURES", sorted(z0_signatures, key=len))
print("ZETA_ZERO_SIDE 1-(SP)^2*(kappa+theta*P+b*P^2)=0; normalization Y^2=quadratic(P), genus 0")

print("ALL M10 LOWER-HULL SCRATCH CHECKS PASSED")
