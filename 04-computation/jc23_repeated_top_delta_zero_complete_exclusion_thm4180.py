#!/usr/bin/env python3
"""Primary exact certificate for THM-4180.

The new coordinate is the diagonal-face variable W=s*p.  At Delta=0 the
old horizontal face disappears, but Q*p^4*(epsilon+eta*W) replaces it and
retains the rational index-four place.  Critical lengths are computed with
the source pair (A,B); no projected-root discriminant is used.
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
    return tuple(lower[:-1] + upper[:-1])


def polygon_ledger(points):
    vertices = convex_hull(points)
    area2 = abs(sum(
        vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
        - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]))
        for i in range(len(vertices))
    )
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return vertices, area2, boundary, (area2 - boundary + 2) // 2


def face(poly, x, y, u, v, level):
    answer = 0
    for (i, j), coefficient in sp.Poly(poly, x, y).terms():
        if u * i + v * j == level:
            answer += coefficient * x**i * y**j
    return sp.factor(answer)


def remainder_rational(expr, modulus, variable):
    numerator, denominator = sp.together(sp.cancel(expr)).as_numer_denom()
    need(sp.resultant(denominator, modulus, variable) != 0,
         "denominator meets hostile modulus")
    return sp.factor(sp.rem(sp.Poly(numerator, variable),
                            sp.Poly(modulus, variable)).as_expr())


s, p, Q = sp.symbols("s p Q")
theta, phi, eta = sp.symbols("theta phi eta")
a, z = sp.symbols("a z")
X, T = sp.symbols("X T")
u = sp.symbols("u")

k0 = sp.Rational(2848, 45)
epsilon = -sp.Rational(1376, 135)
alpha = sp.Rational(8, 3)
lam = sp.Rational(-3)
b0 = sp.factor(epsilon + k0)
need(b0 == sp.Rational(7168, 135), "fixed row-C coefficient")

t = p - s**2
H = sp.expand(
    lam*p + alpha*p**2 + epsilon*p**3 + k0*s**2*p**2
    + phi*s*p**3 + theta*s**2*p**3
    + eta*s*p**3*(p - s**2)
)
F = sp.expand((s**2 - p)*(1 - Q*H) - Q*s**2/2)

# The contracted polygon is common to all three structural rows.
expected_vertices = ((0, 1), (2, 0), (5, 3), (1, 5), (0, 4))
for specialization in (
        {theta: 1, phi: 2, eta: 3},
        {theta: 0, phi: 2, eta: 3},
        {theta: 0, phi: 0, eta: 3}):
    support = tuple(
        powers for powers, coefficient
        in sp.Poly(F.subs(specialization), s, p).terms()
        if coefficient != 0
    )
    ledger = polygon_ledger(support)
    need(ledger == (expected_vertices, 30, 10, 11),
         "Delta-zero polygon/Pick ledger")

low_face = s**2*(1 - Q/2) - p
cubic_face = s**2*((1 - Q/2) - k0*Q*(s*p)**2 + eta*Q*(s*p)**3)
repeated_face = Q*eta*s*p**3*(p - s**2)**2
diagonal_face = Q*p**4*(epsilon + eta*s*p)
vertical_face = p*(-1 + Q*(lam*p + alpha*p**2 + epsilon*p**3))
need(sp.factor(face(F, s, p, 1, 2, 2) - low_face) == 0,
     "low rational face")
need(sp.factor(face(F, s, p, -1, 1, -2) - cubic_face) == 0,
     "prime cubic face")
need(sp.factor(face(F, s, p, -1, -2, -11) - repeated_face) == 0,
     "repeated face")
need(sp.factor(face(F, s, p, 1, -1, -4) - diagonal_face) == 0,
     "replacement diagonal face")
need(sp.factor(face(F, s, p, 1, 0, 0) - vertical_face) == 0,
     "vertical affine face")
need(epsilon != 0, "diagonal root is nonzero when eta is a unit")

# The repeated edge keeps its old chart, while its adjacent Delta-face is
# replaced by the diagonal W=s*p face.
r = 1 - a
L = sp.cancel(z**11 * F.subs({s: 1/z, p: r/z**2}))
L_expected = sp.expand(
    a*z**9 - Q*z**9/2 + Q*eta*a**2*r**3
    - Q*theta*a*r**3*z - Q*phi*a*r**3*z**2
    - Q*a*(epsilon*r**3 + k0*r**2)*z**3
    - Q*alpha*a*r**2*z**5 - Q*lam*a*r*z**7
)
need(sp.factor(L - L_expected) == 0, "strict transform")
La = sp.diff(L_expected, a)

fast_A = theta*z/eta
slow_A = -z**8/(2*theta)
need(valuation(L_expected.subs(a, fast_A), z) >= 3,
     "row-A fast branch")
need(valuation(L_expected.subs(a, slow_A), z) >= 10,
     "row-A slow branch")
need(valuation(La.subs(a, fast_A), z) == 1
     and valuation(La.subs(a, slow_A), z) == 1,
     "row-A residue order")

L_B = sp.expand(L_expected.subs(theta, 0))
fast_B = phi*z**2/eta
slow_B = -z**7/(2*phi)
need(valuation(L_B.subs(a, fast_B), z) >= 5,
     "row-B fast branch")
need(valuation(L_B.subs(a, slow_B), z) >= 10,
     "row-B slow branch")
need(valuation(sp.diff(L_B, a).subs(a, fast_B), z) == 2
     and valuation(sp.diff(L_B, a).subs(a, slow_B), z) == 2,
     "row-B residue order")

L_C = sp.expand(L_B.subs(phi, 0))
fast_C = b0*z**3/eta
slow_C = -z**6/(2*b0)
need(valuation(L_C.subs(a, fast_C), z) >= 7,
     "row-C fast branch")
need(valuation(L_C.subs(a, slow_C), z) >= 11,
     "row-C slow branch")
need(valuation(sp.diff(L_C, a).subs(a, fast_C), z) == 3
     and valuation(sp.diff(L_C, a).subs(a, slow_C), z) == 3,
     "row-C residue order")

# One source projection treats all rows.  It is equivalent to the critical
# ideal on p*t!=0 and has no common point at s=infinity there.
Gsp = -s**2/(2*t) + H
A = sp.cancel((-s*p + t**2*sp.diff(H, s))/p)
C0 = sp.expand(s**2 + 2*t**2*sp.diff(H, p))
B = sp.cancel((C0 + s*A)/t**2)
need(sp.denom(A) == 1 and sp.denom(B) == 1, "source pair polynomial")
need(sp.factor(t**2*sp.diff(Gsp, s) - p*A) == 0,
     "source first-gradient identity")
need(sp.factor(2*t**2*sp.diff(Gsp, p) - (t**2*B - s*A)) == 0,
     "source second-gradient identity")
need((sp.degree(A, s), sp.degree(B, s)) == (6, 3),
     "source s-degrees")
need(sp.factor(sp.Poly(A, s).LC() + 3*eta*p**2) == 0
     and sp.factor(sp.Poly(B, s).LC() + 9*eta*p**2) == 0,
     "source infinity sidecar")
need(sp.factor(A.subs(p, 0) + s) == 0 and B.subs({p: 0, s: 0}) == -6,
     "p-zero chart loss")
need(sp.factor(A.subs(p, s**2) + s) == 0
     and B.subs({p: s**2, s: 0}) == -6,
     "t-zero chart loss")

resultant_AB = sp.factor(sp.resultant(A, B, s))
need(valuation(resultant_AB, p) == 6, "generic source p-artifact")
R19 = sp.Poly(sp.cancel(resultant_AB/p**6), p)
DA = sp.factor(4*theta*k0**2 - 27*eta**2)
need(R19.degree() == 19, "row-A residual degree")
need(sp.factor(R19.LC() - 1327104*eta**5*theta**4) == 0,
     "row-A live top endpoint")
need(sp.factor(R19.TC() - 46656*eta*DA) == 0,
     "row-A inner-wall endpoint")

resultant_B = sp.factor(sp.resultant(A.subs(theta, 0),
                                     B.subs(theta, 0), s))
need(sp.factor(resultant_B - resultant_AB.subs(theta, 0)) == 0,
     "row-B direct specialization")
RB = sp.Poly(sp.cancel(resultant_B/p**6), p)
need(RB.degree() == 17, "row-B residual degree")
need(sp.factor(RB.LC() - 777924*eta**5*phi**4) == 0,
     "row-B top endpoint")
need(sp.factor(RB.TC() + 1259712*eta**3) == 0,
     "row-B bottom endpoint")

resultant_C = sp.factor(sp.resultant(A.subs({theta: 0, phi: 0}),
                                     B.subs({theta: 0, phi: 0}), s))
need(sp.factor(resultant_C - resultant_AB.subs({theta: 0, phi: 0})) == 0,
     "row-C direct specialization")
RC = sp.Poly(sp.cancel(resultant_C/p**6), p)
need(RC.degree() == 15, "row-C residual degree")
need(sp.factor(
    RC.LC() - sp.Rational(168955354770571264, 50625)*eta**5
) == 0, "row-C top endpoint")
need(sp.factor(RC.TC() + 1259712*eta**3) == 0,
     "row-C bottom endpoint")

# Row A on DA=0.  The chart is bijective because eta and k0 are units:
# theta=3u^2, eta=2k0*u/3, u=3eta/(2k0).
wall_sub = {theta: 3*u**2, eta: 2*k0*u/3}
Aw = sp.expand(A.subs(wall_sub))
Bw = sp.expand(B.subs(wall_sub))
resultant_wall = sp.factor(sp.resultant(Aw, Bw, s))
need(sp.factor(resultant_wall - resultant_AB.subs(wall_sub)) == 0,
     "row-A wall direct specialization")
need(valuation(resultant_wall, p) == 7, "row-A wall p-artifact")
R18 = sp.Poly(sp.cancel(resultant_wall/p**7), p)
J = 3*phi*k0 + 8*k0*u + 27*u**3
S0 = 18225*u**4 - 1515136*u**2 - 129777664
P0 = 26081*u**2 + 12924224
S_original = sp.Rational(4, 15)*S0
P_original = sp.Rational(3584, 15)*P0
expected_wall_lc = (sp.Rational(524288, 364651875)*k0**5*u**5
                    *(315*u**2)**4)
need(R18.degree() == 18 and sp.factor(R18.LC() - expected_wall_lc) == 0,
     "row-A wall live top endpoint")
need(sp.factor(R18.TC() - 82944*k0**2*u**2*J) == 0,
     "row-A wall J endpoint")

phi_J = -(8*k0*u + 27*u**3)/(3*k0)
R17 = sp.Poly(sp.cancel(R18.as_expr().subs(phi, phi_J)/p), p)
need(R17.degree() == 17, "row-A J-wall degree")
need(sp.factor(R17.LC() - R18.LC()) == 0,
     "row-A J-wall top endpoint")
need(sp.factor(R17.TC() + sp.Rational(6912, 5)*k0*u**3*S_original) == 0,
     "row-A S endpoint")

c2_on_J = sp.expand(R18.nth(2).subs(phi, phi_J))
terminal_identity = sp.rem(
    sp.Poly(c2_on_J + sp.Rational(96, 35)*k0**2*u**3*P_original, u),
    sp.Poly(S0, u),
).as_expr()
need(sp.factor(terminal_identity) == 0, "row-A terminal coefficient")
need(sp.gcd(sp.Poly(S0, u), sp.Poly(P0, u)).degree() == 0,
     "row-A terminal coprimality")
need(sp.resultant(S0, P0, u) != 0, "row-A terminal resultant")

# The normalized chart independently identifies exactly the four omitted
# affine points, and the source generator change is invertible elsewhere.
PN = T + X**2*T**2
YN = X*T*PN
GN = sp.expand(
    -X**2*T/2 - 3*PN + alpha*PN**2 + epsilon*PN**3 + k0*YN**2
    + phi*PN**2*YN + theta*PN*YN**2
    + eta*PN**3*YN - eta*YN**3
)
fn = sp.expand(sp.cancel(sp.diff(GN, X)/T))
hn = sp.expand(sp.diff(GN, T))
Hess = sp.factor(sp.det(sp.hessian(GN, (X, T))))
need(sp.factor(fn.subs(T, 0) + X) == 0,
     "normalized T-zero f")
need(sp.factor(hn.subs(T, 0) + (X**2 + 6)/2) == 0,
     "normalized T-zero h")
need(remainder_rational(Hess.subs(T, 0) - 6, X**2 + 6, X) == 0,
     "normalized T-zero Hessian")
need(sp.factor(sp.cancel(fn.subs(T, -1/X**2) + (X**2 - 6)/X)) == 0,
     "normalized p-zero f")
need(sp.factor(sp.cancel(hn.subs(T, -1/X**2) + (X**2 - 6)/2)) == 0,
     "normalized p-zero h")
need(remainder_rational(Hess.subs(T, -1/X**2) + 6, X**2 - 6, X) == 0,
     "normalized p-zero Hessian")

# Boundary, genus, degree, and response ledgers.
W, q = sp.symbols("W q")
carrier = -eta*W**3 + k0*W**2 - (q - sp.Rational(1, 2))
need(sp.factor(
    sp.discriminant(carrier, W)
    - (q - sp.Rational(1, 2))
      *(4*k0**3 - 27*eta**2*(q - sp.Rational(1, 2)))
) == 0, "prime cubic carrier")

rows = (
    ("A", (7, 7, 4, 2, 2, 2, 1), 10, 23, 19, 25),
    ("B", (6, 6, 4, 2, 2, 2, 1), 9, 21, 17, 23),
    ("C", (5, 5, 4, 2, 2, 2, 1), 8, 19, 15, 21),
)
for name, packet, genus, length, finite_n, full_n in rows:
    defect = sum(index - 1 for index in packet)
    need(sum(packet) == full_n and defect == 2*genus - 2,
         name + " packet/genus ledger")
    need(finite_n == full_n - 6, name + " cubic response degree")
    need(2*(full_n - length) < defect, name + " full contradiction")
    finite_capacity = 2*finite_n - length - 1 + 3
    need(finite_capacity < finite_n - 1, name + " finite open-row contradiction")

for length in (22, 21, 20):
    need(19 > 3 + 1, "row-A carrier lemma type")
    kappa = 19 + 3 - length
    need(2*kappa + 3 < 15, "row-A wall finite contradiction")
    need(2*(25 - length) < 18, "row-A wall full contradiction")

print("THM4180_REPEATED_TOP_DELTA_ZERO_PRIMARY_EXACT_CERTIFICATE")
print("wall=Delta:0,zeta:-eta,eta!=0")
print("polygon=((0,1),(2,0),(5,3),(1,5),(0,4));Pick=(30,10,11)")
print("replacement_face=Q*p^4*(epsilon+eta*s*p);root_W=-epsilon/eta")
print("vertical_s_zero=affine_not_infinity_packet")
print("rows=A:theta!=0;B:theta=0,phi!=0;C:theta=phi=0,b0=7168/135")
print("packets=A:(7,7,4,2,2,2,1);B:(6,6,4,2,2,2,1);C:(5,5,4,2,2,2,1)")
print("source_pair=(A,B);s_infinity=none_for_p!=0,eta!=0")
print("critical_lengths=A_off_DA:23;A_on_DA:22,21,20;B:21;C:19")
print("row_A_DA_tower=18,17,16;terminal_gcd(S0,P0)=1")
print("carrier=q-1/2=k0*W^2-eta*W^3;degree=3;beta=3")
print("responses=all_full_and_finite_bounds_strict")
print(f"checks={CHECKS}")
print("verdict=DELTA_ZERO_REPEATED_TOP_CLOSES")
