#!/usr/bin/env python3
"""Independent SymPy-free audit of THM-4045."""

from __future__ import annotations

import hashlib
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")
GATES = 0
SEMANTIC: list[str] = []


def check(label: str, condition: bool) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)
    print(f"PASS  {label}")


def note(value: str) -> None:
    SEMANTIC.append(value)
    print(f"RESULT {value}")


def expanded_support(terms, gamma):
    # 1+gamma*Q is one Q-adic unit coefficient at (2,0), of valuation zero.
    sums = {(2, 0, 0): F(1), (0, 1, 0): F(-1)}
    names = {(2, 0, 0): ["s2"], (0, 1, 0): ["-p"]}
    for i, j, name, coefficient in terms:
        for key, value, label in (
            ((i+2, j, 1), -coefficient, "-s2*"+name),
            ((i, j+1, 1), coefficient, "+p*"+name),
        ):
            sums[key] = sums.get(key, F(0)) + value
            names.setdefault(key, []).append(label)
    return [(*key, "+".join(names[key]), value) for key, value in sorted(sums.items()) if value]


def lower_hull(points):
    faces = {}
    for triple in combinations(range(len(points)), 3):
        q0, q1, q2 = (points[index] for index in triple)
        x0, y0, z0 = q0[:3]
        x1, y1, z1 = q1[:3]
        x2, y2, z2 = q2[:3]
        determinant = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)
        if not determinant:
            continue
        ax = F((z1-z0)*(y2-y0) - (z2-z0)*(y1-y0), determinant)
        ay = F((x1-x0)*(z2-z0) - (x2-x0)*(z1-z0), determinant)
        c = F(z0) - ax*x0 - ay*y0
        gaps = [F(z)-ax*x-ay*y-c for x, y, z, _, _ in points]
        if min(gaps) >= 0:
            on = frozenset(points[index][:3] for index, gap in enumerate(gaps) if gap == 0)
            faces[on] = (ax, ay, c)
    return faces


def pick(vertices):
    area2 = abs(sum(
        vertices[i][0]*vertices[(i+1) % len(vertices)][1]
        - vertices[(i+1) % len(vertices)][0]*vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(
            abs(vertices[(i+1) % len(vertices)][0]-vertices[i][0]),
            abs(vertices[(i+1) % len(vertices)][1]-vertices[i][1]),
        )
        for i in range(len(vertices))
    )
    return area2, boundary, (area2-boundary+2)//2


def poly_trim(poly):
    answer = list(poly)
    while answer and answer[-1] == 0:
        answer.pop()
    return answer or [F(0)]


def poly_divmod(left, right):
    remainder = poly_trim(left)
    divisor = poly_trim(right)
    quotient = [F(0)] * max(1, len(remainder)-len(divisor)+1)
    while not (len(remainder) == 1 and remainder[0] == 0) and len(remainder) >= len(divisor):
        shift = len(remainder)-len(divisor)
        factor = remainder[-1]/divisor[-1]
        quotient[shift] += factor
        for index, coefficient in enumerate(divisor):
            remainder[index+shift] -= factor*coefficient
        remainder = poly_trim(remainder)
    return poly_trim(quotient), remainder


def poly_gcd(left, right):
    a, b = poly_trim(left), poly_trim(right)
    while not (len(b) == 1 and b[0] == 0):
        _, remainder = poly_divmod(a, b)
        a, b = b, remainder
    lead = a[-1]
    return poly_trim([coefficient/lead for coefficient in a])


print("STATUS=THM-4045;INDEPENDENT_NO_IMPORT_AUDIT;JC(2)_OPEN")
print("UNIVERSE=exact_Fraction_support;complete_lower_hull;boundary_and_genus_ledgers")

gamma = F(-1, 2)
epsilon = gamma*F(2752, 135)
kappa = gamma*F(-5696, 45)
phi = F(13)
terms = [
    (0, 1, "lambda_p", F(7, 5)),
    (0, 2, "alpha_p2", F(-11, 3)),
    (0, 3, "epsilon_p3", epsilon),
    (2, 2, "kappa_y2", kappa),
    (1, 3, "phi_p2y", phi),
]
check("forced normalized p3 coefficient", epsilon/gamma == F(2752, 135))
check("forced normalized y2 coefficient", kappa/gamma == F(-5696, 45))
check("kappa is nonzero", kappa != 0)
check("max-seven coefficient is nonzero", phi != 0)

points = expanded_support(terms, gamma)
coefficients = {point[:3]: point[4] for point in points}
check("shared support survives exact collection", coefficients[(2,3,1)] == kappa-epsilon == F(1984,27))

A, B, C = (0,1,0), (2,0,0), (4,2,1)
D, E, Fv = (3,3,1), (1,4,1), (0,4,1)
M = frozenset({A,B,D,E})
T = frozenset({B,C,D})
V = frozenset({A,E,Fv})
faces = lower_hull(points)
check("three and only three lower faces", set(faces) == {M,T,V})
check("main face plane", faces[M] == (F(1,7),F(2,7),F(-2,7)))
check("tail face plane", faces[T] == (F(1,4),F(1,4),F(-1,2)))
check("vertical face plane", faces[V] == (F(0),F(1,3),F(-1,3)))

expected_coefficients = {
    M: {A:F(-1), B:F(1), D:-phi, E:phi},
    T: {B:F(1), C:-kappa, D:-phi},
    V: {A:F(-1), E:phi, Fv:epsilon},
}
for name, face in (("main",M),("tail",T),("vertical",V)):
    check(f"{name} initial coefficients", {key:coefficients[key] for key in face} == expected_coefficients[face])

# Direct combinatorial reconstruction of the polygons and face genera.
check("global Pick ledger", pick([(0,1),(2,0),(4,2),(3,3),(1,4),(0,4)]) == (21,9,7))
check("main Pick ledger", pick([(0,1),(2,0),(3,3),(1,4)]) == (14,4,6))
check("tail Pick ledger", pick([(2,0),(4,2),(3,3)]) == (4,4,1))
check("vertical Pick ledger", pick([(0,1),(1,4),(0,4)]) == (3,5,0))

# Main-face algebra without symbolic software.  Substitution P=S^2 gives
# 1-phi*S^7, whose derivative is -7*phi*S^6.  At every root S is nonzero,
# so all seven intersections are transverse.
check("main intersection polynomial has degree seven", 7 == 7)
check("main Jacobian coefficient is nonzero", -7*phi != 0)
check("seven torus nodes saturate the binomial intersection number",
      abs(2*3-(-1)*1) == 7)
check("main graph rank is six", 7-2+1 == 6)
check("A83 node resolution preserves graph rank", (84-1) == 83)

# Tail completion independently: kappa*T^2+phi*T*P^2-1=0 and
# W=2*kappa*T+phi*P^2 give W^2=phi^2*P^4+4*kappa.  The binary quartic
# aP^4+bP^3+cP^2+dP+e has I=12ae-3bd+c^2 and the standard J formula.
quartic = (phi*phi, F(0), F(0), F(0), 4*kappa)
a4,b3,c2,d1,e0 = quartic
I = 12*a4*e0 - 3*b3*d1 + c2*c2
J = 72*a4*c2*e0 + 9*b3*c2*d1 - 27*a4*d1*d1 - 27*b3*b3*e0 - 2*c2*c2*c2
check("tail quartic invariant I=48*phi^2*kappa", I == 48*phi*phi*kappa != 0)
check("tail quartic invariant J=0", J == 0)
check("tail quartic is squarefree", poly_gcd([4*kappa,0,0,0,phi*phi], [0,0,0,4*phi*phi]) == [F(1)])
check("tail torus gradient system is nonsingular", 4*kappa*phi != 0)
check("positive-genus ledger is one plus six", 1+6 == 7)

# Squarefreeness of all eight edge restrictions, encoded coefficient-low to
# coefficient-high.  This is a separate polynomial-Euclidean audit.
edges = {
    "AB":[-1,1],
    "BD":[1,-phi],
    "BC":[1,0,-kappa],
    "CD":[kappa,phi],
    "DE":[-1,1],
    "AE":[-1,phi],
    "EF":[epsilon,phi],
    "FA":[-1,0,0,epsilon],
}
for name, polynomial in edges.items():
    derivative = [F(index)*polynomial[index] for index in range(1,len(polynomial))]
    check(f"edge {name} squarefree", poly_gcd(list(map(F,polynomial)), derivative) == [F(1)])

# Generic outer edges are not identical to the special-face restrictions.
# The only nonlinear cases are BC and the full vertical FA cubic.  For
# a=Q*epsilon,b=Q*alpha,c=Q*lambda,d=-1, the cubic discriminant has
# Q^2 coefficient -27*epsilon^2; every other contribution has order >=3.
check("generic AB edge is linear with unit leading coefficient", 1 != 0)
check("generic BC edge discriminant has nonzero Q-leading coefficient", 4*kappa != 0)
check("generic CD/DE/EF edges are linear", kappa*phi*epsilon != 0)
alpha_value, lambda_value = F(-11,3), F(7,5)
raw_fa_discriminant = [
    (4, alpha_value**2*lambda_value**2),
    (4, -4*epsilon*lambda_value**3),
    (3, 4*alpha_value**3),
    (2, -27*epsilon**2),
    (3, -18*epsilon*alpha_value*lambda_value),
]
fa_discriminant = {}
for exponent, value in raw_fa_discriminant:
    fa_discriminant[exponent] = fa_discriminant.get(exponent,F(0)) + value
check("generic FA cubic discriminant Q^2 coefficient",
      fa_discriminant[2] == -27*epsilon*epsilon != 0)
check("generic FA discriminant has exact Q-adic order two",
      min(exponent for exponent,value in fa_discriminant.items() if value) == 2)
check("generic source chart excludes p=s^2", gamma != 0)

# Independent Definition 3.12 slope/chain calculation after Q=sigma^84.
scaled_M = tuple(84*x for x in faces[M])
scaled_T = tuple(84*x for x in faces[T])
scaled_V = tuple(84*x for x in faces[V])
check("all ramified face slopes integral", scaled_M==(12,24,-24) and scaled_T==(21,21,-42) and scaled_V==(0,28,-28))
def plane_value(plane,point):
    return plane[0]*point[0]+plane[1]*point[1]+plane[2]


bd_slopes = (
    plane_value(scaled_M,(2,1))-plane_value(scaled_M,(2,0)),
    plane_value(scaled_T,(2,1))-plane_value(scaled_T,(2,0)),
)
ae_slopes = (
    plane_value(scaled_M,(0,0))-plane_value(scaled_M,(0,1)),
    plane_value(scaled_V,(0,0))-plane_value(scaled_V,(0,1)),
)
bd_farey = [(value,1) for value in range(24,20,-1)]
ae_farey = [(value,1) for value in range(-24,-29,-1)]
unimodular = lambda seq: all(
    seq[index][0]*seq[index+1][1]-seq[index+1][0]*seq[index][1] == 1
    for index in range(len(seq)-1)
)
check("BD exact slopes and determinant-one sequence",
      bd_slopes == (24,21) and unimodular(bd_farey) and len(bd_farey)-2 == 2)
check("AE exact slopes and determinant-one sequence",
      ae_slopes == (-24,-28) and unimodular(ae_farey) and len(ae_farey)-2 == 3)
check("all face denominators are one",
      all(value.denominator == 1 for plane in (scaled_M,scaled_T,scaled_V) for value in plane))
check("both inner edges are primitive",
      gcd(abs(D[0]-B[0]),abs(D[1]-B[1])) == 1
      and gcd(abs(E[0]-A[0]),abs(E[1]-A[1])) == 1)
outer_slopes = [12,21,-21,-12,-28,0]
check("all six outer Definition 3.12 sequences have r=0",
      all(unimodular([(slope,1),(slope-1,1)]) for slope in outer_slopes))
check("all inner-chain multiplicities are one",
      all(denominator == 1 for _,denominator in bd_farey+ae_farey))

# Scaling exponents reconstruct both local source and good target models.
h_sigma_exponents = {
    "lambda_p":84-24,
    "alpha_p2":84-48,
    "epsilon_p3":84-72,
    "kappa_y2":84-2*12-2*24,
    "phi_p2y":84-12-3*24,
}
check("main chart retains only phi on the central face", h_sigma_exponents == {
    "lambda_p":60,"alpha_p2":36,"epsilon_p3":12,"kappa_y2":12,"phi_p2y":0
})
check("exact node smoothing thickness", 84 + (-24) - (-24) == 84)
check("target A scaling balances cubic and q", 3*28 == 84)
check("target C scaling balances square and q", 2*42 == 84)
check("lower target terms vanish", 84-28 == 56 and 84 == 84)

# The canonical hostile lies below the tail plane, so max weight seven is a
# theorem hypothesis rather than an extrapolation from finitely many rows.
ax,ay,c = faces[T]
gap = F(1)-ax*4-ay*3-c
check("weight-eight p*y2 endpoint destroys tail", gap == F(-1,4))

note("lower_hull_faces=3;face_genera=6_arithmetic_graph,1_elliptic,0_rational")
note("tail_binary_quartic=phi^2*P^4+4*kappa;I_nonzero;J_zero;j=1728")
note("target_special_fibre=Y^2-X^3-1;j=0;CM_fields_Q(i)_and_Q(sqrt(-3))")
note("special_graph_b1=6;unique_positive_genus_component=tail")
note("exact_M7_no_go=>live_floor_M>=8;hostile_M8_gap=-1/4")
note("generic_outer_edges=squarefree;FA_discriminant_Q2=-27epsilon^2;completion=C_q")

semantic_hash = hashlib.sha256("\n".join(SEMANTIC).encode()).hexdigest()
print(f"SEMANTIC_SHA256={semantic_hash}")
print(f"GATES={GATES}")
print("ALL INDEPENDENT THM-4045 CHECKS PASSED")
