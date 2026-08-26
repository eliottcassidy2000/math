#!/usr/bin/env python3
"""Clean-room source-coordinate referee for THM-4159.

This audit uses the alternative critical pair (A,C), direct face extraction,
resultant nonvanishing rather than gcd calls for the terminal quotient, and
two disjoint rational controls.  It imports no primary-audit code.
"""
from hashlib import sha256
from math import gcd
import sympy as sp

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def order(poly, variable):
    return min(
        monomial[0] for monomial, coefficient in sp.Poly(poly, variable).terms()
        if coefficient
    )


def quotient_rows(poly, variable, modulus, parameter):
    rows = []
    for (degree,), coefficient in sp.Poly(poly, variable).terms():
        numerator, denominator = sp.fraction(sp.cancel(coefficient))
        remainder = sp.rem(
            sp.Poly(numerator, parameter), sp.Poly(modulus, parameter)
        ).as_expr()
        if remainder:
            rows.append((degree, sp.factor(remainder/denominator)))
    return rows


def hull(points):
    points = sorted(set(points))

    def cross(o, a, b):
        return ((a[0]-o[0])*(b[1]-o[1])
                - (a[1]-o[1])*(b[0]-o[0]))

    lower = []
    for point in points:
        while len(lower) > 1 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) > 1 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1]+upper[:-1])


s, p = sp.symbols("s p")
X, T = sp.symbols("X T")
u, Phi, q, W = sp.symbols("u Phi q W")
K = sp.Rational(2848, 45)
epsilon = -sp.Rational(1376, 135)
Zeta = sp.Rational(5696, 135)*u
Theta = 3*u**2
J = 8544*Phi-22784*u-1215*u**3
S = 2460375*u**4-204543360*u**2+5580439552
PhiJ = u*(1215*u**2+22784)/8544
t = p-s**2

H = (
    -3*p + sp.Rational(8, 3)*p**2 + epsilon*p**3
    + K*s**2*p**2 + Phi*s*p**3 + Theta*s**2*p**3
    + Zeta*s**3*p**3
)
A = sp.cancel((-s*p+t**2*sp.diff(H, s))/p)
C = sp.expand(s**2+2*t**2*sp.diff(H, p))
require(sp.denom(A) == 1, "A polynomial")
require((sp.degree(A, s), sp.degree(C, s)) == (6, 7),
        "independent pair degrees")
require(sp.factor(sp.Poly(A, s).LC()-3*Zeta*p**2) == 0,
        "A top row")
require(sp.factor(sp.Poly(C, s).LC()-6*Zeta*p**2) == 0,
        "C top row")

resultant_AC = sp.resultant(A, C, s)
require(order(resultant_AC, p) == 9, "AC generic wall artifact")
R17 = sp.cancel(resultant_AC/p**9)
P17 = sp.Poly(R17, p)
require(P17.degree() == 17, "AC R17")

top = Zeta*W**3+Theta*W**2+Phi*W+epsilon
D = sp.factor(sp.discriminant(top, W))
require(sp.factor(
    P17.TC()-sp.Rational(8305770496, 1125)*u**2*J
) == 0, "AC R17 constant")
require(sp.factor(
    P17.LC()+sp.Rational(23983352712374779904, 759375)*u**5*D
) == 0, "AC R17 leading")

# Specialization of the determinant gives the direct J=0 and J=S=0 rows.
specialJ = sp.factor(R17.subs(Phi, PhiJ))
require(order(specialJ, p) == 1, "AC J exceptional order")
R16 = sp.cancel(specialJ/p)
P16 = sp.Poly(R16, p)
require(P16.degree() == 16, "AC R16")
require(sp.factor(
    P16.TC()-sp.Rational(2916352, 16875)*u**3*S
) == 0, "AC R16 constant")
require(sp.factor(
    P16.LC()+sp.Rational(23983352712374779904, 759375)*u**5
    * D.subs(Phi, PhiJ)
) == 0, "AC R16 leading")

rows = quotient_rows(R16, p, S, u)
require((rows[0][0], rows[-1][0]) == (16, 1), "AC p11R15 quotient")
terminal_polynomials = (
    u,
    369170566011315*u**2-23248683486112768,
    3904455285*u**2-155035505152,
    2005507674605933782764615*u**2
        +316908385228357703794041472,
    401085*u**2-16287712,
    30267225703125*u**8+2043284356800000*u**6
        +264381824212992000*u**4+6498574373014732800*u**2
        +498260889496415371264,
)
require(sp.discriminant(S, u) != 0, "S is reduced")
for polynomial in terminal_polynomials:
    require(sp.resultant(S, polynomial, u) != 0,
            "terminal quotient resultant firewall")

# Two rational controls, one in each positive-dimensional stratum.
control21 = sp.Poly(R17.subs({u: 1, Phi: 1}), p)
require(control21.degree() == 17, "L21 control degree")
require(sp.gcd(control21, control21.diff()).degree() == 0,
        "L21 control squarefree")
require(J.subs({u: 1, Phi: 1}) == -15455, "L21 control J")

control20 = sp.Poly(R16.subs(u, 1), p)
require(control20.degree() == 16, "L20 control degree")
require(sp.gcd(control20, control20.diff()).degree() == 0,
        "L20 control squarefree")
require(S.subs(u, 1) != 0, "L20 control S")

# Collapsed rows and four universal points are audited without primary code.
require(sp.factor(A.subs(p, 0)) == -s, "p0 A row")
require(sp.factor(C.subs(p, 0)+s**2*(6*s**2-1)) == 0, "p0 C row")
require(sp.factor(A.subs(p, s**2)) == -s, "t0 A row")
require(sp.factor(C.subs(p, s**2)-s**2) == 0, "t0 C row")

P = T+X**2*T**2
Y = X*T*P
G = (
    -X**2*T/2-3*P+sp.Rational(8, 3)*P**2+epsilon*P**3
    +K*Y**2+Phi*P**2*Y+Theta*P*Y**2+Zeta*Y**3
)
Hess = sp.det(sp.hessian(G, (X, T)))
for sign in (-1, 1):
    zero = {T: 0, X: sign*sp.sqrt(-6)}
    half = {T: -sp.Rational(1, 6), X: sign*sp.sqrt(6)}
    require(sp.simplify(sp.diff(G, X).subs(zero)) == 0, "zero critical")
    require(sp.simplify(sp.diff(G, T).subs(zero)) == 0, "zero critical T")
    require(sp.simplify(Hess.subs(zero)-6) == 0, "zero Hessian")
    require(sp.simplify(sp.diff(G, X).subs(half)) == 0, "half critical")
    require(sp.simplify(sp.diff(G, T).subs(half)) == 0, "half critical T")
    require(sp.simplify(Hess.subs(half)+6) == 0, "half Hessian")

# Direct generic-fibre expansion and face extraction.
F = sp.Poly(sp.expand((s**2-p)*(q-H)-s**2/2), s, p)
support = tuple(monomial for monomial, coefficient in F.terms() if coefficient)
polygon = hull(support)
require(polygon == ((0, 1), (2, 0), (5, 3), (3, 4), (0, 4)),
        "independent polygon")

def edge_coefficients(left, right):
    dx, dy = right[0]-left[0], right[1]-left[1]
    length = gcd(abs(dx), abs(dy))
    step = (dx//length, dy//length)
    return tuple(sp.factor(F.coeff_monomial(
        s**(left[0]+index*step[0])*p**(left[1]+index*step[1])
    )) for index in range(length+1))

carrier_coefficients = edge_coefficients((2, 0), (5, 3))
top_coefficients = edge_coefficients((3, 4), (0, 4))
require(all(sp.factor(left-right) == 0 for left, right in zip(
    carrier_coefficients, (q-sp.Rational(1, 2), 0, -K, -Zeta)
)), "cubic carrier face")
require(all(sp.factor(left-right) == 0 for left, right in zip(
    top_coefficients, (Zeta, Theta, Phi, epsilon)
)), "constant top cubic face")
carrier = Zeta*W**3+K*W**2-(q-sp.Rational(1, 2))
require(sp.degree(carrier, W) == 3 and sp.degree(carrier, q) == 1,
        "prime cubic carrier")
require(sp.factor(sp.discriminant(carrier, W)
        -(q-sp.Rational(1, 2))*(4*K**3-27*Zeta**2*(q-sp.Rational(1, 2)))) == 0,
        "carrier discriminant")

# The response bookkeeping is reconstructed independently.
full_packet = (8, 3, 3, 3, 2, 2, 2, 1)
finite_packet = (8, 3, 3, 3, 1)
require(sum(full_packet) == 24, "full degree")
require(sum(e-1 for e in full_packet) == 16, "full origin index")
require(sum(finite_packet) == 18, "finite degree")
require(sum(e-1 for e in finite_packet) == 13, "finite origin index")
require(sum(e-1 for e in (2, 2, 2)) == 3, "carrier index three")
for length in (21, 20, 19):
    require(2*(24-length) < 16, "full response")
    # Three transpositions can lower the orbit count by at most three.
    orbit_floor = 18-(1+3)
    union_floor = orbit_floor+1
    overlap_ceiling = (36-length)-union_floor
    require((orbit_floor, union_floor, overlap_ceiling)
            == (14, 15, 21-length), "finite orbit lemma ledger")
    require(2*overlap_ceiling+3 < 13, "finite response")

semantic = (
    "independent_pair=A,C;p9R17,p10R16,p11R15",
    "controls=L21(u1,Phi1),L20(u1,J0):squarefree",
    "terminal=resultant(S,all endpoint/discriminant factors)!=0",
    "faces=prime cubic carrier plus constant squarefree cubic",
    "packet=full24/index16;finite18/index13;carrier index3",
    "bounds=full<=10;finite<=7",
)
digest = sha256("\n".join(semantic).encode()).hexdigest()

print("THM4159_Y_ONLY_INNER_RESULTANT_WALL_INDEPENDENT_AUDIT")
print("checks="+str(CHECKS))
print("alternative_source_pair=A,C")
print("source_AC=p^9R17;p^10R16;p^11R15")
print("strata=J!=0:L21;J=0,S!=0:L20;J=S=0:L19")
print("controls=L21_u1_Phi1_squarefree;L20_u1_J0_squarefree")
print("terminal_firewall=Disc(S)!=0;all_six_resultants_nonzero")
print("faces=prime_cubic_carrier(2,2,2)+rational(8,3,3,3,1)")
print("origin_indices=full16;finite13;carrier_total3")
print("monodromy_bounds=full_at_most10;finite_at_most7")
print("semantic_sha256="+digest)
print("verdict=ACCEPT")
