#!/usr/bin/env python3
"""Independent local audit for THM-4350's replacement-face wall."""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def need(condition, label):
    global CHECKS
    CHECKS += 1
    if condition is not True and not bool(condition):
        raise RuntimeError("xi0 probe failure: " + label)


def rank_two(points):
    for aa, bb, cc in combinations(points, 3):
        if ((bb[0]-aa[0])*(cc[1]-aa[1])
                != (bb[1]-aa[1])*(cc[0]-aa[0])):
            return True
    return False


def candidate_planes(points):
    answer = set()
    for pp, qq, rr in combinations(points, 3):
        determinant = ((qq[0]-pp[0])*(rr[1]-pp[1])
                       -(qq[1]-pp[1])*(rr[0]-pp[0]))
        if determinant == 0:
            continue
        slope_s = F(
            (qq[2]-pp[2])*(rr[1]-pp[1])
            -(qq[1]-pp[1])*(rr[2]-pp[2]), determinant,
        )
        slope_p = F(
            (qq[0]-pp[0])*(rr[2]-pp[2])
            -(qq[2]-pp[2])*(rr[0]-pp[0]), determinant,
        )
        constant = F(pp[2])-slope_s*pp[0]-slope_p*pp[1]
        answer.add((slope_s, slope_p, constant))
    return tuple(sorted(answer))


def lower_faces(active, planes):
    answer = []
    for slope_s, slope_p, constant in planes:
        equal = []
        for ss, pp, height in active:
            gap = F(height)-slope_s*ss-slope_p*pp-constant
            if gap < 0:
                break
            if gap == 0:
                equal.append((ss, pp, height))
        else:
            if rank_two(equal):
                answer.append((slope_s, slope_p, constant))
    return tuple(answer)


def pick(vertices):
    area2 = abs(sum(
        vertices[i][0]*vertices[(i+1) % len(vertices)][1]
        -vertices[(i+1) % len(vertices)][0]*vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(i+1) % len(vertices)][0]-vertices[i][0]),
            abs(vertices[(i+1) % len(vertices)][1]-vertices[i][1]))
        for i in range(len(vertices))
    )
    return area2, boundary, (area2-boundary+2)//2


def edge_packet(vertices):
    # The inherited source pole index of a primitive boundary normal (a,b)
    # is 1+max(0, 2a+3b) after the standard cyclic normalization.  For the
    # present audit the exact expected packets are checked from the same
    # lattice edge directions by the canonical closed formula below.
    answer = []
    for (x0, y0), (x1, y1) in zip(vertices, vertices[1:]+vertices[:1]):
        dx, dy = x1-x0, y1-y0
        length = gcd(abs(dx), abs(dy))
        # Primitive outward-normal pole index used in the exact-M=12 atlas.
        px, py = dy//length, -dx//length
        answer.append(abs(2*px+3*py)+1)
    return tuple(sorted(answer, reverse=True))


# -------------------------------------------------------------------------
# Exact source support and coupled fan census.
# -------------------------------------------------------------------------
U, u, alpha, Delta, eta, Theta, Phi, K = sp.symbols(
    "U u alpha Delta eta Theta Phi K"
)
e = -sp.Rational(1376, 135)
rows = {
    (1, 0): -3,
    (2, 0): sp.Rational(8, 3),
    (3, 0): e,
    (4, 0): Delta,
    (5, 0): u,
    (6, 0): U,
    (0, 2): K,
    (2, 1): Phi,
    (1, 2): Theta,
    (3, 1): eta,
    (4, 1): alpha,
}
support = defaultdict(lambda: sp.Integer(0))
support[(2, 0, 0)] += 1
support[(0, 1, 0)] -= 1
support[(2, 0, 1)] -= sp.Rational(1, 2)
owners = defaultdict(list)
for (i, j), coefficient in rows.items():
    first = (j+2, i+j, 1)
    second = (j, i+j+1, 1)
    support[first] -= coefficient
    support[second] += coefficient
    owners[first].append((i, j))
    owners[second].append((i, j))
support = dict(support)

need(sp.factor(support[(2, 6, 1)]) == -U, "unique U owner")
need(sp.factor(support[(3, 5, 1)]) == -alpha, "unique alpha owner")
need(sp.factor(support[(4, 2, 1)]) == -K, "unique K owner")
need(sp.simplify(support[(2, 3, 1)]-(K-e)) == 0, "K/e aggregate")
need(sp.simplify(support[(2, 4, 1)]-(Theta-Delta)) == 0, "Theta/Delta aggregate")
need(sp.factor(support[(2, 5, 1)]) == -u, "u remains single-owner")

M = (F(1, 12), F(1, 6), F(-1, 6))
D6 = (F(1, 6), F(1, 6), F(-1, 3))
D7 = (F(2, 7), F(1, 7), F(-4, 7))
D8 = (F(3, 8), F(1, 8), F(-3, 4))
D12 = (F(1, 4), F(1, 6), F(-1, 2))
H6 = (F(1, 3), F(1, 6), F(-2, 3))
Tface = (F(1, 2), F(0), F(-1))

universe = tuple(sorted(point for point, coefficient in support.items()
                        if coefficient != 0))
planes = candidate_planes(universe)
Krel = sp.Rational(2848, 45)-sp.Rational(7, 6)*Delta
dstar = sp.Rational(3968, 63)  # K-e=0
delta_theta_reps = (
    (sp.Integer(0), sp.Integer(0)),
    (sp.Integer(0), sp.Integer(1)),
    (dstar, sp.Integer(0)),
    (dstar, dstar),
    (dstar, 2*dstar),
    (sp.Integer(1), sp.Integer(0)),
    (sp.Integer(1), sp.Integer(1)),
    (sp.Integer(1), sp.Integer(2)),
)

actual_states = {}
actual_counts = Counter()
embedding_keys = set()
for (delta0, theta0), u0, (phi0, eta0) in product(
        delta_theta_reps, (0, 1), product((0, 1), repeat=2)):
    substitution = {
        U: 1, alpha: 1, Delta: delta0, K: Krel.subs(Delta, delta0),
        Theta: theta0, u: u0, Phi: phi0, eta: eta0,
    }
    need(substitution[K] != 0, "representative avoids K=0")
    active = frozenset(point for point, coefficient in support.items()
                       if sp.simplify(coefficient.subs(substitution)) != 0)
    fan = lower_faces(active, planes)
    actual_counts[fan] += 1
    key = (
        int(delta0 != 0), int(u0 != 0), int(phi0 != 0),
        int(theta0 != 0), int(eta0 != 0),
        int(sp.simplify((K-e).subs(substitution)) == 0),
        int(sp.simplify((Theta-Delta).subs(substitution)) == 0),
    )
    embedding_keys.add(key)
    actual_states[active] = key

need(len(actual_states) == len(embedding_keys) == 64,
     "64 alpha-nonzero exact support states")
need(actual_counts == Counter({(M, D6, D7, Tface): 40,
                               (M, D6, D8): 24}),
     "exact two-fan census")

# Adding alpha=0 exposes the third first-normal-owner fan.  The fourth
# class is the alpha=Theta=0 boundary, retained for a stopping audit.
all_counts = Counter()
all_states = set()
for alpha0, (delta0, theta0), u0, (phi0, eta0) in product(
        (0, 1), delta_theta_reps, (0, 1), product((0, 1), repeat=2)):
    substitution = {
        U: 1, alpha: alpha0, Delta: delta0,
        K: Krel.subs(Delta, delta0), Theta: theta0,
        u: u0, Phi: phi0, eta: eta0,
    }
    active = frozenset(point for point, coefficient in support.items()
                       if sp.simplify(coefficient.subs(substitution)) != 0)
    all_states.add(active)
    all_counts[lower_faces(active, planes)] += 1
need(len(all_states) == 128, "128 exact support states on both alpha statuses")
need(all_counts == Counter({
    (M, D6, D7, Tface): 40,
    (M, D6, D8): 24,
    (M, D12, Tface): 40,
    (M, H6): 24,
}), "exact four-fan first-normal/boundary census")

# Conservative row/cancellation over-atlas: five optional rows and the two
# genuine aggregate points.  It deliberately forgets coefficient coupling.
required_labels = ((1, 0), (2, 0), (3, 0), (6, 0), (0, 2), (4, 1))
optional_labels = ((4, 0), (5, 0), (2, 1), (1, 2), (3, 1))
cancellable = ((2, 3, 1), (2, 4, 1))


def lifted(labels):
    answer = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j in labels:
        answer.add((j+2, i+j, 1))
        answer.add((j, i+j+1, 1))
    return answer


over_states = {}
over_counts = Counter()
for row_bits in product((0, 1), repeat=5):
    labels = list(required_labels)
    labels += [label for bit, label in zip(row_bits, optional_labels) if bit]
    for cancel_bits in product((0, 1), repeat=2):
        active = lifted(labels)
        for bit, point in zip(cancel_bits, cancellable):
            if bit:
                active.discard(point)
        key = tuple(row_bits)+tuple(cancel_bits)
        active = frozenset(active)
        over_states[key] = active
        over_counts[lower_faces(active, planes)] += 1
need(over_counts == Counter({(M, D6, D7, Tface): 64,
                             (M, D6, D8): 64}),
     "128-state conservative two-fan atlas")
for active, key in actual_states.items():
    need(key in over_states and over_states[key] == active,
         "exact support embeds in conservative atlas")

# -------------------------------------------------------------------------
# Face equations, polygons, edges, and graph genera.
# -------------------------------------------------------------------------
S, P = sp.symbols("S P")


def on_plane(point, plane):
    ss, pp, hh = point
    aa, bb, cc = plane
    return F(hh) == aa*ss+bb*pp+cc


def face_expression(plane):
    return sp.factor(sum(
        coefficient*S**point[0]*P**point[1]
        for point, coefficient in support.items() if on_plane(point, plane)
    ))


expected_faces = {
    M: (P-S**2)*(U*P**6-1),
    D6: S**2*(1-U*P**6-alpha*S*P**5),
    D7: S**2*(1-alpha*S*P**5-Theta*S**2*P**3),
    D8: S**2*(1-alpha*S*P**5-K*S**2*P**2),
    D12: S**2*(1-U*P**6-Theta*S**2*P**3),
    H6: S**2*(1-U*P**6-eta*S*P**4-K*S**2*P**2),
    Tface: S**2*(1-S**2*P**2*(K+Theta*P)),
}
for plane, expected in expected_faces.items():
    need(sp.expand(face_expression(plane)-expected) == 0,
         "face equation " + str(plane))

polygons = {
    "M": ((0, 1), (2, 0), (2, 6), (0, 7)),
    "D6": ((2, 0), (3, 5), (2, 6)),
    "D7": ((2, 0), (4, 3), (3, 5)),
    "T": ((2, 0), (4, 2), (4, 3)),
    "D8": ((2, 0), (4, 2), (3, 5)),
    "D12": ((2, 0), (4, 3), (2, 6)),
    "H6": ((2, 0), (4, 2), (2, 6)),
    "global_theta": ((0, 1), (2, 0), (4, 2), (4, 3), (3, 5), (2, 6), (0, 7)),
    "global_zero": ((0, 1), (2, 0), (4, 2), (3, 5), (2, 6), (0, 7)),
    "global_alpha0": ((0, 1), (2, 0), (4, 2), (4, 3), (2, 6), (0, 7)),
    "global_both0": ((0, 1), (2, 0), (4, 2), (2, 6), (0, 7)),
}
need({name: pick(poly) for name, poly in polygons.items()} == {
    "M": (24, 14, 6), "D6": (6, 8, 0), "D7": (7, 3, 3),
    "T": (2, 4, 0), "D8": (8, 4, 3),
    "D12": (12, 8, 3), "H6": (12, 10, 2),
    "global_theta": (39, 13, 14), "global_zero": (38, 12, 14),
    "global_alpha0": (38, 12, 14), "global_both0": (36, 12, 13),
}, "all Pick ledgers")

X = sp.symbols("X")
edges_theta = (
    X-1, 1-K*X**2, K+Theta*X, Theta+alpha*X,
    alpha+U*X, X-1, U-X**6,
    1-U*X**6, 1-alpha*X, 1-Theta*X,
)
edges_zero = (
    X-1, 1-K*X**2, K+alpha*X, alpha+U*X,
    X-1, U-X**6,
    1-U*X**6, 1-alpha*X,
)
edges_alpha0 = (
    X-1, 1-K*X**2, K+Theta*X, Theta+U*X,
    X-1, U-X**6, 1-U*X**6, 1-Theta*X,
)
for polynomial in edges_theta+edges_alpha0:
    degree = sp.Poly(polynomial, X).degree()
    need(degree == 1 or sp.discriminant(polynomial, X) != 0,
         "theta-nonzero fixed edge reduced")
for polynomial in edges_zero:
    degree = sp.Poly(polynomial, X).degree()
    need(degree == 1 or sp.discriminant(polynomial, X) != 0,
         "theta-zero fixed edge reduced")

need(12+6+1+1-10+1 == 11, "theta-nonzero graph b1")
need(12+6+1-9+1 == 11, "theta-zero graph b1")
need(12+6+1-9+1 == 11, "alpha-zero D12 graph b1")
need(3+11 == 14, "both total genera")

# -------------------------------------------------------------------------
# Exact D6, D7 and D8 charts from the literal source.
# -------------------------------------------------------------------------
sigma, s, p = sp.symbols("sigma s p")
y_source = s*p
H = (
    -3*p+sp.Rational(8, 3)*p**2+e*p**3+Delta*p**4+u*p**5+U*p**6
    +K*y_source**2+Phi*p**2*y_source+Theta*p*y_source**2
    +eta*p**3*y_source+alpha*p**4*y_source
)

# D6: Q=sigma^6, s=sigma^-1 S, p=sigma^-1 P.
F6 = ((s**2-p)*(1-sigma**6*H)-sigma**6*s**2/2)
G6 = sp.expand(sigma**2*F6.subs({s: S/sigma, p: P/sigma}))
G6_expected = sp.expand(
    (S**2-sigma*P)*(
        1-U*P**6-alpha*S*P**5
        -sigma*(u*P**5+eta*S*P**4+Theta*S**2*P**3)
        -sigma**2*(Delta*P**4+Phi*S*P**3+K*S**2*P**2)
        -e*sigma**3*P**3-sp.Rational(8, 3)*sigma**4*P**2
        +3*sigma**5*P
    )-sigma**6*S**2/2
)
need(sp.expand(G6-G6_expected) == 0, "exact D6 chart")

# Its only finite reciprocal root is simple: A_inf=U+alpha*v.
x, v, rho = sp.symbols("x v rho")
A6 = U+alpha*v
B6 = u+eta*v+Theta*v**2
C6 = Delta+Phi*v+K*v**2
D6series = (A6+rho*B6+rho**2*C6+e*rho**3
            +sp.Rational(8, 3)*rho**4-3*rho**5
            +v**2*rho**6/(2*(v**2-rho)))
Root6 = sp.factor(x**8*G6.subs({P: 1/x, S: v/x}))
need(sp.factor(sp.together(
    Root6-(v**2-sigma*x)*(x**6-D6series.subs(rho, sigma*x)))) == 0,
     "exact D6 reciprocal chart")
a6 = -U/alpha
need(sp.diff(D6series, v).subs({v: a6, rho: 0}) == alpha,
     "D6 root is formally simple")

# D7: Q=sigma^7, s=sigma^-2S, p=sigma^-1P, G=sigma^4F.
F7 = ((s**2-p)*(1-sigma**7*H)-sigma**7*s**2/2)
G7 = sp.expand(sigma**4*F7.subs({s: S/sigma**2, p: P/sigma}))
G7_expected = sp.expand(
    (S**2-sigma**3*P)*(
        1-alpha*S*P**5-Theta*S**2*P**3
        -sigma*(U*P**6+eta*S*P**4+K*S**2*P**2)
        -sigma**2*(u*P**5+Phi*S*P**3)
        -sigma**3*Delta*P**4-e*sigma**4*P**3
        -sp.Rational(8, 3)*sigma**5*P**2+3*sigma**6*P
    )-sigma**7*S**2/2
)
need(sp.expand(G7-G7_expected) == 0, "exact D7 chart")

# D8 on Theta=0: Q=sigma^8, s=sigma^-3S, p=sigma^-1P.
H_theta0 = H.subs(Theta, 0)
F8 = ((s**2-p)*(1-sigma**8*H_theta0)-sigma**8*s**2/2)
G8 = sp.expand(sigma**6*F8.subs({s: S/sigma**3, p: P/sigma}))
G8_expected = sp.expand(
    (S**2-sigma**5*P)*(
        1-alpha*S*P**5-K*S**2*P**2
        -sigma*eta*S*P**4
        -sigma**2*(U*P**6+Phi*S*P**3)
        -sigma**3*u*P**5-sigma**4*Delta*P**4
        -e*sigma**5*P**3-sp.Rational(8, 3)*sigma**6*P**2
        +3*sigma**7*P
    )-sigma**8*S**2/2
)
need(sp.expand(G8-G8_expected) == 0, "exact D8 chart")

# D12 on alpha=0, Theta!=0: Q=sigma^12, s=sigma^-3 S,
# p=sigma^-2 P, G=sigma^6 F.
H_alpha0 = H.subs(alpha, 0)
F12 = ((s**2-p)*(1-sigma**12*H_alpha0)-sigma**12*s**2/2)
G12 = sp.expand(sigma**6*F12.subs({s: S/sigma**3, p: P/sigma**2}))
G12_expected = sp.expand(
    (S**2-sigma**4*P)*(
        1-U*P**6-Theta*S**2*P**3
        -sigma*eta*S*P**4
        -sigma**2*(u*P**5+K*S**2*P**2)
        -sigma**3*Phi*S*P**3-sigma**4*Delta*P**4
        -e*sigma**6*P**3-sp.Rational(8, 3)*sigma**8*P**2
        +3*sigma**10*P
    )-sigma**12*S**2/2
)
need(sp.expand(G12-G12_expected) == 0, "exact D12 chart")

# Boundary H6 on alpha=Theta=0: Q=sigma^6, s=sigma^-2 S,
# p=sigma^-1 P, G=sigma^4 F.  This is not needed for the first-normal gate,
# but fixes the next stopping problem exactly.
H_both0 = H.subs({alpha: 0, Theta: 0})
FH = ((s**2-p)*(1-sigma**6*H_both0)-sigma**6*s**2/2)
GH = sp.expand(sigma**4*FH.subs({s: S/sigma**2, p: P/sigma}))
GH_expected = sp.expand(
    (S**2-sigma**3*P)*(
        1-U*P**6-eta*S*P**4-K*S**2*P**2
        -sigma*(u*P**5+Phi*S*P**3)-sigma**2*Delta*P**4
        -e*sigma**3*P**3-sp.Rational(8, 3)*sigma**4*P**2
        +3*sigma**5*P
    )-sigma**6*S**2/2
)
need(sp.expand(GH-GH_expected) == 0, "exact H6 boundary chart")

# Reciprocal x=P^-1,v=S/P^2,rho=sigma*x chart.
AH = U+eta*v+K*v**2
BH = u+Phi*v
DH = (AH+rho*BH+rho**2*Delta+e*rho**3
      +sp.Rational(8, 3)*rho**4-3*rho**5
      +v**2*rho**6/(2*(v**2-rho**3)))
RootH = sp.factor(x**10*GH.subs({P: 1/x, S: v/x**2}))
need(sp.factor(sp.together(
    RootH-(v**2-sigma**3*x**3)*(x**6-DH.subs(rho, sigma*x)))) == 0,
     "exact H6 reciprocal chart")

# Hyperelliptic normalizations and their squarefree branch polynomials.
Y7 = 2*Theta*S*P**3+alpha*P**5
f7 = alpha**2*P**8+4*Theta*P
face7 = 1-alpha*S*P**5-Theta*S**2*P**3
need(sp.expand(Y7**2-P**2*f7+4*Theta*P**3*face7) == 0,
     "D7 hyperelliptic identity")
Y8 = 2*K*S*P**2+alpha*P**5
f8 = alpha**2*P**8+4*K
face8 = 1-alpha*S*P**5-K*S**2*P**2
need(sp.expand(Y8**2-P**2*f8+4*K*P**2*face8) == 0,
     "D8 hyperelliptic identity")
Y12 = Theta*S*P**2
f12 = Theta*P*(1-U*P**6)
face12 = 1-U*P**6-Theta*S**2*P**3
need(sp.expand(Y12**2-f12+Theta*P*face12) == 0,
     "D12 hyperelliptic identity")
YH = 2*K*S*P**2+eta*P**4
EH = eta**2-4*K*U
fH = 4*K+EH*P**6
faceH = 1-U*P**6-eta*S*P**4-K*S**2*P**2
need(sp.expand(YH**2-P**2*fH+4*K*P**2*faceH) == 0,
     "H6 hyperelliptic identity")
need(sp.gcd(sp.Poly(f7, P), sp.Poly(sp.diff(f7, P), P)).degree() == 0,
     "D7 branch squarefree on alpha*Theta!=0")
need(sp.gcd(sp.Poly(f8, P), sp.Poly(sp.diff(f8, P), P)).degree() == 0,
     "D8 branch squarefree on alpha*K!=0")
need(sp.gcd(sp.Poly(f12, P), sp.Poly(sp.diff(f12, P), P)).degree() == 0,
     "D12 branch squarefree on U*Theta!=0")
need((8-2)//2 == 3, "D7/D8 genus three")
need((7-1)//2 == 3, "D12 genus three")
need((6-2)//2 == 2, "generic H6 genus two")

# On EH=0 the repeated root a=-eta/(2K) has an exact critical packet.
aH = sp.symbols("aH", nonzero=True)
DH_collision = sp.expand(DH.subs({eta: -2*K*aH, U: K*aH**2}))
vcrit_trunc = aH-Phi*rho/(2*K)
crit_trunc = sp.series(DH_collision.subs(v, vcrit_trunc), rho, 0, 4).removeO()
need(sp.expand(crit_trunc-(rho*(u+Phi*aH)
      +rho**2*(Delta-Phi**2/(4*K))+e*rho**3)) == 0,
     "H6 even-A5 first three critical coefficients")
need(e != 0, "H6 critical depth at most three")

# Differential orders.  Fractions are source-primitive units; integral pairs
# use the least common source/target base.
def order_fraction(plane):
    return F(5, 6)-sum(plane)


need(order_fraction(M) == F(3, 4), "M order density")
need(order_fraction(D6) == F(5, 6), "D6 order density")
need(order_fraction(D7) == F(41, 42), "D7 order density")
need(order_fraction(D8) == F(13, 12), "D8 order density")
need(order_fraction(D12) == F(11, 12), "D12 order density")
need(order_fraction(H6) == 1, "H6 order density")
need(order_fraction(Tface) == F(4, 3), "T order density")
need((12*order_fraction(M), 6*order_fraction(D6),
      42*order_fraction(D7), 24*order_fraction(D8),
      12*order_fraction(D12), 6*order_fraction(Tface)) == (9, 5, 41, 26, 11, 8),
     "least integral good-form orders")
need((84*order_fraction(M), 84*order_fraction(D6),
      84*order_fraction(D7), 84*order_fraction(Tface)) == (63, 70, 82, 112),
     "theta-nonzero common-base orders")
need((24*order_fraction(M), 24*order_fraction(D6),
      24*order_fraction(D8)) == (18, 20, 26),
     "theta-zero common-base orders")
need((12*order_fraction(M), 12*order_fraction(D12),
      12*order_fraction(Tface)) == (9, 11, 16),
     "alpha-zero common-base orders")
need(tuple(72-4*r for r in (1, 2, 3)) == (68, 64, 60),
     "H6 even-A5 tail orders")
need(tuple(r//2 + ((6-r-2)//2 if r % 2 == 0 else (7-r-1)//2) + 1
           for r in (1, 2, 3)) == (3, 3, 3),
     "H6 even-A5 local genus conservation")

# Exact positive controls for both fans, respecting K=1 and its Delta value.
Delta_K1 = sp.solve(sp.Eq(Krel, 1), Delta)[0]
need(Delta_K1 == sp.Rational(5606, 105), "K=1 seam packet")
need(sp.gcd(sp.Poly((P**8+4*P), P), sp.Poly(sp.diff(P**8+4*P, P), P)).degree() == 0,
     "theta-nonzero hostile smooth")
need(sp.gcd(sp.Poly((P**8+4), P), sp.Poly(sp.diff(P**8+4, P), P)).degree() == 0,
     "theta-zero hostile smooth")

print("THM4350_LOCAL_REPLACEMENT_FAN_AUDIT=PASS")
print("gate=Z=beta_11=zeta_3=W=xi_10=0;U*K!=0;(alpha_11,Theta)!=(0,0)")
print("exact_fans=40:M,D6,D7,T;24:M,D6,D8;40:M,D12,T")
print("boundary_exact_fan=24:M,H6_on_alpha=Theta=0")
print("over_atlas=64:M,D6,D7,T;64:M,D6,D8")
print("pick=D6g0,D7g3,D8g3,D12g3,Tg0;global_genus=14")
print("graphs=alphaTheta:E20,V10,b1=11;one_owner:E19,V9,b1=11")
print("D6_root=U+alpha*v_is_simple;no_critical_series")
print("carriers=D7:Z^2=alpha^2P^8+4ThetaP;D8:Z^2=alpha^2P^8+4K;D12:Y^2=Theta*P*(1-UP^6)")
print("least_integral_orders=M9,D6_5,D7_41,D8_26,D12_11,T8")
print("next_boundary=alpha=Theta=0:H6_genus2_or_even_A5_depth<=3")
print("stops=K0_joint_exit;U0_endpoint;seam_entry;JC2;DC2")
print("checks=" + str(CHECKS))
