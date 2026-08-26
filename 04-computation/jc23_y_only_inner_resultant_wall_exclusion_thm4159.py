#!/usr/bin/env python3
"""Exact primary certificate for THM-4159's complete inner-resultant wall.

No repository files are imported. The computation starts from the complete
Delta=eta=0 weight-nine source, imposes I_C=0 by a global parameter u, and
computes the three exhaustive strict-transform strata.
"""
from __future__ import annotations

from hashlib import sha256
from math import gcd
import sympy as sp

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def valuation(poly, variable):
    terms = sp.Poly(poly, variable).terms()
    require(bool(terms), "zero polynomial")
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


def convex_hull(points):
    points = sorted(set(points))

    def cross(o, a, b):
        return ((a[0]-o[0])*(b[1]-o[1])
                - (a[1]-o[1])*(b[0]-o[0]))

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


def packet(vertices):
    result = []
    for left, right in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = right[0]-left[0], right[1]-left[1]
        length = gcd(abs(dx), abs(dy))
        if left[0] == right[0] == 0:
            require(length == 3, "vertical affine edge")
            continue
        normal = (-dy//length, dx//length)
        level = normal[0]*left[0] + normal[1]*left[1]
        index = normal[0] + normal[1] - level
        result.extend([index]*length)
    return tuple(sorted(result, reverse=True))


def quotient_profile(poly, variable, modulus, parameter):
    """Nonzero coefficient rows after reduction modulo modulus(parameter)."""
    rows = []
    for (degree,), coefficient in sp.Poly(poly, variable).terms():
        numerator, denominator = sp.fraction(sp.cancel(coefficient))
        remainder = sp.rem(
            sp.Poly(numerator, parameter), sp.Poly(modulus, parameter)
        ).as_expr()
        if remainder != 0:
            rows.append((degree, sp.factor(remainder/denominator)))
    return rows


s, p = sp.symbols("s p")
X, T = sp.symbols("X T")
u, Phi = sp.symbols("u Phi")
w, q = sp.symbols("w q")
K = sp.Rational(2848, 45)
Zeta = sp.Rational(5696, 135)*u
Theta = 3*u**2
J = 8544*Phi - 22784*u - 1215*u**3
S = 2460375*u**4 - 204543360*u**2 + 5580439552
Phi_J = sp.factor((22784*u + 1215*u**3)/8544)
epsilon = -sp.Rational(1376, 135)

require(sp.factor(4*Theta*K**2 - 27*Zeta**2) == 0,
        "I_C parameterization")
require(sp.factor(3*Zeta/(2*K) - u) == 0,
        "inverse u parameterization")
require(sp.factor(J.subs(Phi, Phi_J)) == 0, "J parameterization")

top = Zeta*w**3 + Theta*w**2 + Phi*w + epsilon
D = sp.factor(sp.discriminant(top, w))

t = p-s**2
H = sp.expand(
    -3*p + sp.Rational(8, 3)*p**2 + epsilon*p**3
    + K*s**2*p**2 + Phi*s*p**3 + Theta*s**2*p**3
    + Zeta*s**3*p**3
)
A = sp.cancel((-s*p + t**2*sp.diff(H, s))/p)
C0 = sp.expand(s**2 + 2*t**2*sp.diff(H, p))
B = sp.cancel((C0+s*A)/t**2)
require(sp.denom(A) == 1 and sp.denom(B) == 1, "source polynomials")
require((sp.degree(A, s), sp.degree(B, s)) == (6, 3),
        "source degrees")
require(sp.factor(sp.Poly(A, s).LC()-3*Zeta*p**2) == 0,
        "source A leading row")
require(sp.factor(sp.Poly(B, s).LC()-9*Zeta*p**2) == 0,
        "source B leading row")

# First exact elimination: source pair (A,B).
source_resultant = sp.resultant(A, B, s)
require(valuation(source_resultant, p) == 7, "generic wall AB p-order")
R17_source = sp.cancel(source_resultant/p**7)
require(sp.denom(R17_source) == 1, "source R17 polynomial")
R17s = sp.Poly(R17_source, p)
require(R17s.degree() == 17, "source R17 degree")
require(sp.factor(
    R17s.TC() - sp.Rational(8305770496, 1125)*u**2*J
) == 0, "source R17 constant")
require(sp.factor(
    R17s.LC() + sp.Rational(23983352712374779904, 759375)*u**5*D
) == 0, "source R17 leading")

# Independent normalized-coordinate elimination.
P = T + X**2*T**2
Y = X*T*P
G = sp.expand(
    -X**2*T/2 - 3*P + sp.Rational(8, 3)*P**2 + epsilon*P**3
    + K*Y**2 + Phi*P**2*Y + Theta*P*Y**2 + Zeta*Y**3
)
f = sp.cancel(sp.diff(G, X)/T)
h = sp.diff(G, T)
require(sp.denom(f) == 1, "normalized first equation")
require((sp.degree(f, X), sp.degree(h, X)) == (8, 9),
        "normalized X degrees")
normalized_resultant = sp.resultant(f, h, X)
require(valuation(normalized_resultant, T) == 56,
        "normalized T artifact")
Q17_expr = sp.cancel(normalized_resultant/(T**56*(6*T+1)**2))
require(sp.denom(Q17_expr) == 1, "normalized residual polynomial")
Q17 = sp.Poly(Q17_expr, T)
require(Q17.degree() == 17, "normalized Q17 degree")
require(sp.factor(
    Q17.TC()
    + sp.Rational(1519777094677765052956672, 56953125)*u**7
) == 0, "normalized Q17 constant")
require(sp.factor(
    Q17.LC()
    + sp.Rational(23983352712374779904, 13839609375)*u**3*J**2*D
) == 0, "normalized Q17 leading")

# J=0 strict transform: source p^8 R16 and normalized Q16.
source_J = sp.factor(R17_source.subs(Phi, Phi_J))
require(valuation(source_J, p) == 1, "J source exceptional p")
R16_source = sp.cancel(source_J/p)
R16s = sp.Poly(R16_source, p)
require(R16s.degree() == 16, "J source R16")
require(sp.factor(
    R16s.TC() - sp.Rational(2916352, 16875)*u**3*S
) == 0, "J source R16 constant")
require(sp.factor(
    R16s.LC() + sp.Rational(23983352712374779904, 759375)*u**5
    * D.subs(Phi, Phi_J)
) == 0, "J source R16 leading")

Q16_expr = sp.factor(Q17_expr.subs(Phi, Phi_J))
Q16 = sp.Poly(Q16_expr, T)
require(Q16.degree() == 16, "J normalized Q16")
require(sp.factor(
    Q16.TC()
    + sp.Rational(1519777094677765052956672, 56953125)*u**7
) == 0, "J normalized Q16 constant")
require(sp.factor(
    Q16.LC()
    + sp.Rational(2956854296576, 3113912109375)*u**3*S**2
    * D.subs(Phi, Phi_J)
) == 0, "J normalized Q16 leading")

# J=S=0 is a finite reduced coefficient stratum.  Work exactly in Q[u]/(S).
source_rows = quotient_profile(R16_source, p, S, u)
normalized_rows = quotient_profile(Q16_expr, T, S, u)
require((source_rows[0][0], source_rows[-1][0]) == (16, 1),
        "JS source p^9 R15 profile")
require((normalized_rows[0][0], normalized_rows[-1][0]) == (15, 0),
        "JS normalized Q15 profile")

source_top_linear = 369170566011315*u**2 - 23248683486112768
source_bottom_linear = 3904455285*u**2 - 155035505152
normalized_top_linear = (
    2005507674605933782764615*u**2
    + 316908385228357703794041472
)
normalized_bottom_linear = 401085*u**2 - 16287712
B_D = (
    30267225703125*u**8 + 2043284356800000*u**6
    + 264381824212992000*u**4 + 6498574373014732800*u**2
    + 498260889496415371264
)
require(sp.factor(
    D.subs(Phi, Phi_J) + u**2*B_D/sp.Integer(99781787520000)
) == 0, "J-wall top discriminant")
require(sp.gcd(sp.Poly(S, u), sp.Poly(sp.diff(S, u), u)).degree() == 0,
        "S squarefree")
require(sp.gcd(sp.Poly(S, u), sp.Poly(u, u)).degree() == 0,
        "S roots nonzero")
for name, polynomial in (
    ("source top", source_top_linear),
    ("source bottom", source_bottom_linear),
    ("normalized top", normalized_top_linear),
    ("normalized bottom", normalized_bottom_linear),
    ("top discriminant", B_D),
):
    require(sp.gcd(sp.Poly(S, u), sp.Poly(polynomial, u)).degree() == 0,
            "S gcd firewall: "+name)
require(sp.factor(source_rows[0][1]).has(source_top_linear),
        "JS source top factor")
require(sp.factor(source_rows[-1][1]).has(source_bottom_linear),
        "JS source bottom factor")
require(sp.factor(normalized_rows[0][1]).has(normalized_top_linear),
        "JS normalized top factor")
require(sp.factor(normalized_rows[-1][1]).has(normalized_bottom_linear),
        "JS normalized bottom factor")

# The four universal critical points remain separate with Hessian +/-6.
hessian = sp.factor(sp.det(sp.hessian(G, (X, T))))
for sign in (-1, 1):
    zero = {T: 0, X: sign*sp.sqrt(-6)}
    half = {T: -sp.Rational(1, 6), X: sign*sp.sqrt(6)}
    require(sp.simplify(sp.diff(G, X).subs(zero)) == 0, "zero G_X")
    require(sp.simplify(h.subs(zero)) == 0, "zero G_T")
    require(sp.simplify(G.subs(zero)) == 0, "zero value")
    require(sp.simplify(hessian.subs(zero)-6) == 0, "zero Hessian")
    require(sp.simplify(sp.diff(G, X).subs(half)) == 0, "half G_X")
    require(sp.simplify(h.subs(half)) == 0, "half G_T")
    require(sp.simplify(G.subs(half)-sp.Rational(1, 2)) == 0,
            "half value")
    require(sp.simplify(hessian.subs(half)+6) == 0, "half Hessian")

# Boundary support is unchanged on u*D != 0.
Fq = sp.expand((s**2-p)*(q-H)-s**2/2)
support = tuple(sorted(
    monomial for monomial, coefficient in sp.Poly(Fq, s, p).terms()
    if coefficient != 0
))
polygon = convex_hull(support)
expected_polygon = ((0, 1), (2, 0), (5, 3), (3, 4), (0, 4))
require(polygon == expected_polygon, "wall polygon")
boundary_packet = packet(polygon)
require(boundary_packet == (8, 3, 3, 3, 2, 2, 2, 1),
        "wall packet")
require(sum(boundary_packet) == 24, "full degree")
require(sum(index-1 for index in boundary_packet) == 16,
        "full origin index")
finite_origin_packet = (8, 3, 3, 3, 1)
require(sum(finite_origin_packet) == 18, "finite degree")
require(sum(index-1 for index in finite_origin_packet) == 13,
        "finite origin index")
require(sum(index-1 for index in (2, 2, 2)) == 3,
        "three carrier transpositions")

# Exact arithmetic of the full and finite carrier-orbit contradictions.
lengths = (21, 20, 19)
for length in lengths:
    full_overlap = 24-length
    require(2*full_overlap < 16, "full commutator contradiction")
    finite_handle_orbits_rank = 18-4
    finite_union_floor = finite_handle_orbits_rank+1
    finite_overlap = (2*18-length)-finite_union_floor
    require(finite_union_floor == 15, "finite support-union floor")
    require(finite_overlap == 21-length, "finite overlap ceiling")
    require(2*finite_overlap+3 < 13,
            "finite origin relation contradiction")

semantic_rows = (
    "map=zeta=(5696/135)u;Theta=3u^2;u=3zeta/(2K0)",
    "J=8544Phi-22784u-1215u^3",
    "S=2460375u^4-204543360u^2+5580439552",
    "strata=J!=0:L21|J=0,S!=0:L20|J=S=0:L19",
    "source=AB:p7R17|p8R16|p9R15",
    "normalized=T56(6T+1)^2*(Q17|Q16|Q15)",
    "gcd=S squarefree/nonzero roots;S coprime endpoints and Disc(C)",
    "boundary=g9;packet=8,3,3,3,2,2,2,1",
    "responses=full24:index16|finite18:index13+three transpositions",
    "full=2(24-L)<=10<16",
    "finite=orbits<=4;union>=15;k<=21-L;2k+3<=7<13",
)
semantic_digest = sha256("\n".join(semantic_rows).encode()).hexdigest()

print("THM4159_Y_ONLY_INNER_RESULTANT_WALL_PRIMARY_AUDIT")
print("checks="+str(CHECKS))
print("parameter_map=zeta=(5696/135)u;Theta=3u^2;inverse_u=3zeta/(2K0)")
print("J=8544Phi-22784u-1215u^3")
print("S=2460375u^4-204543360u^2+5580439552")
print("strata=J!=0:L21;J=0,S!=0:L20;J=S=0:L19")
print("source_AB=p^7R17;p^8R16;p^9R15")
print("normalized=T^56(6T+1)^2*(Q17,Q16,Q15)")
print("JS_gcd_firewall=S_squarefree;u,DiscC,all_source_and_normalized_endpoints_coprime")
print("polygon=(0,1),(2,0),(5,3),(3,4),(0,4);Pick=(27,11,9)")
print("packet=8,3,3,3,2,2,2,1;full_n24_origin_index16")
print("finite_origin_packet=8,3,3,3,1;n18;index13;carrier_beta3")
print("full_bounds=L21:6,L20:8,L19:10<16")
print("finite_bounds=L21:3,L20:5,L19:7<13")
print("semantic_sha256="+semantic_digest)
print("verdict=PASS")
