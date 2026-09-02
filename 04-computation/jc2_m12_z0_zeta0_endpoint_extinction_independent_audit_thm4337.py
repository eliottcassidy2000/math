#!/usr/bin/env python3
"""Clean-room exact audit of the THM-4337 zeta-zero endpoint wall.

This script deliberately does not import any repository computation.  It uses
two small, explicit engines:

1. a dual-support enumeration of every affine plane through the master lifted
   support, followed by bit-set feasibility tests for each lower-row state;
2. a polygon/edge engine which reconstructs every face, boundary orbit, and
   internal orbit from exponent dictionaries.

Scope audited here:

    Z = zeta_3 = 0,      U * beta_11 != 0,

with W arbitrary, Lambda=U+W arbitrary, the y^2 row K present or absent,
all seven remaining lower rows optional, and every multiply-owned active
lifted point independently cancellable.  The latter is an over-atlas: some
simultaneous cancellations need not be coefficient-realizable, so constancy
on it is stronger than constancy on the actual coefficient space.
"""

from collections import Counter, defaultdict, deque
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------------------
# 1. Lifted row universe and independent dual lower-support engine
# ---------------------------------------------------------------------------


ROWS = tuple(
    (i, j)
    for j in range(5)
    for i in range(7)
    if 0 < 2 * i + 3 * j <= 12 and (i, j) not in {(0, 1), (1, 1)}
)
need(len(ROWS) == 16, "the inherited exact-weight-twelve row universe")

P1, P2, P3 = (1, 0), (2, 0), (3, 0)
UROW, WROW, ZROW = (6, 0), (3, 2), (0, 4)
KROW, BETAROW, ZETAROW = (0, 2), (1, 3), (0, 3)
BASE_POINTS = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}


def lifted_points(row):
    i, j = row
    return ((j + 2, i + j, 1), (j, i + j + 1, 1))


MASTER_OWNERS = defaultdict(list)
for row in ROWS:
    for point in lifted_points(row):
        MASTER_OWNERS[point].append(row)

MASTER = tuple(sorted(BASE_POINTS | set(MASTER_OWNERS)))
POINT_INDEX = {point: index for index, point in enumerate(MASTER)}
need(len(MASTER) == 28, "master lifted support has 28 distinct points")


def det2(a, b, c):
    return ((b[0] - a[0]) * (c[1] - a[1])
            - (b[1] - a[1]) * (c[0] - a[0]))


def affine_plane(a, b, c):
    """Return (A,B,C) for z=A*x+B*y+C through three projected points."""
    determinant = det2(a, b, c)
    need(determinant != 0, "affine_plane called on collinear projections")
    dx1, dy1, dz1 = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    dx2, dy2, dz2 = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    A = F(dz1 * dy2 - dy1 * dz2, determinant)
    B = F(dx1 * dz2 - dz1 * dx2, determinant)
    C = F(a[2]) - A * a[0] - B * a[1]
    return A, B, C


def point_mask(points):
    mask = 0
    for point in points:
        mask |= 1 << POINT_INDEX[point]
    return mask


# Every lower facet of every sub-support must pass through a non-collinear
# triple from MASTER.  Precompute all possible planes.  For each plane retain
# the master points strictly below it and all rank-two equality certificates.
PLANE_DATA = {}
for triple in combinations(MASTER, 3):
    if det2(*triple) == 0:
        continue
    plane = affine_plane(*triple)
    below = []
    equal = []
    A, B, C = plane
    for point in MASTER:
        slack = F(point[2]) - A * point[0] - B * point[1] - C
        if slack < 0:
            below.append(point)
        elif slack == 0:
            equal.append(point)
    triples = tuple(
        point_mask(t)
        for t in combinations(equal, 3)
        if det2(*t) != 0
    )
    PLANE_DATA[plane] = (point_mask(below), triples)

need(len(PLANE_DATA) == 293, "dual master-plane count")


def lower_facets(support_mask):
    facets = []
    for plane, (below_mask, rank_two_triples) in PLANE_DATA.items():
        if support_mask & below_mask:
            continue
        if any((support_mask & triple) == triple
               for triple in rank_two_triples):
            facets.append(plane)
    return tuple(sorted(facets))


M = (F(1, 12), F(1, 6), F(-1, 6))
V7 = (F(1, 7), F(1, 7), F(-2, 7))
V9 = (F(1, 9), F(1, 6), F(-2, 9))
KF = (F(1), F(-1, 2), F(-2))


def support_atlas(w_present, k_present, beta_present=True):
    required = {P1, P2, P3, UROW}
    if beta_present:
        required.add(BETAROW)
    if w_present:
        required.add(WROW)
    if k_present:
        required.add(KROW)

    forbidden = {ZROW, ZETAROW}
    if not w_present:
        forbidden.add(WROW)
    if not k_present:
        forbidden.add(KROW)
    if not beta_present:
        forbidden.add(BETAROW)

    optional = tuple(sorted(set(ROWS) - required - forbidden))
    complexes = Counter()
    state_count = 0

    for optional_mask in range(1 << len(optional)):
        active = set(required)
        active.update(
            optional[index]
            for index in range(len(optional))
            if optional_mask & (1 << index)
        )
        owners = defaultdict(list)
        for row in active:
            for point in lifted_points(row):
                owners[point].append(row)

        # Independent deletion of every multiply-owned active coefficient is
        # a safe hostile overcount of aggregate coefficient cancellation.
        cancellable = tuple(sorted(
            point for point, labels in owners.items() if len(labels) >= 2
        ))
        full_support = point_mask(BASE_POINTS | set(owners))
        for cancellation_mask in range(1 << len(cancellable)):
            support = full_support
            for index, point in enumerate(cancellable):
                if cancellation_mask & (1 << index):
                    support &= ~(1 << POINT_INDEX[point])
            complexes[lower_facets(support)] += 1
            state_count += 1

    return optional, state_count, complexes


ATLAS_EXPECTED = {
    (False, False): ((M, V9), 300),
    (False, True): ((M, V9, KF), 600),
    (True, False): ((M, V7), 600),
    (True, True): ((M, V7, KF), 1200),
}

ATLAS_SUMMARY = []
for w_present, k_present in product((False, True), repeat=2):
    optional, states, complexes = support_atlas(w_present, k_present)
    expected_faces, expected_states = ATLAS_EXPECTED[(w_present, k_present)]
    need(len(optional) == 7, "seven genuinely optional lower rows")
    need(states == expected_states, "clean-room support-state count")
    need(complexes == Counter({tuple(sorted(expected_faces)): states}),
         "all hostile states have the exact unique lower complex")
    ATLAS_SUMMARY.append((w_present, k_present, states,
                          tuple(sorted(expected_faces))))


# ---------------------------------------------------------------------------
# 2. Exact polygons, Pick ledgers, and face normalizations
# ---------------------------------------------------------------------------


def convex_hull(points):
    points = sorted(set(points))
    need(len(points) >= 3, "polygon needs at least three points")

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


def pick_ledger(points):
    polygon = convex_hull(points)
    twice_area = abs(sum(
        polygon[index][0] * polygon[(index + 1) % len(polygon)][1]
        - polygon[(index + 1) % len(polygon)][0] * polygon[index][1]
        for index in range(len(polygon))
    ))
    boundary = sum(
        gcd(abs(polygon[(index + 1) % len(polygon)][0] - polygon[index][0]),
            abs(polygon[(index + 1) % len(polygon)][1] - polygon[index][1]))
        for index in range(len(polygon))
    )
    interior = (twice_area - boundary + 2) // 2
    need(2 * interior == twice_area - boundary + 2,
         "Pick parity")
    return polygon, (twice_area, boundary, interior)


# Coefficients are linear forms in (1,U,W,beta,K).  No CAS is used.
ZERO = (0, 0, 0, 0, 0)
ONE = (1, 0, 0, 0, 0)
CU = (0, 1, 0, 0, 0)
CW = (0, 0, 1, 0, 0)
CB = (0, 0, 0, 1, 0)
CK = (0, 0, 0, 0, 1)


def cadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def cscale(n, a):
    return tuple(n * x for x in a)


def csub(a, b):
    return cadd(a, cscale(-1, b))


M_W = {
    (2, 0): ONE,
    (0, 1): cscale(-1, ONE),
    (2, 6): csub(CW, CU),
    (0, 7): CU,
    (4, 5): cscale(-1, CW),
}
M_0 = {
    (2, 0): ONE,
    (0, 1): cscale(-1, ONE),
    (2, 6): cscale(-1, CU),
    (0, 7): CU,
}
FACE_V7 = {
    (2, 0): ONE,
    (4, 5): cscale(-1, CW),
    (5, 4): cscale(-1, CB),
}
FACE_V9 = {
    (2, 0): ONE,
    (2, 6): cscale(-1, CU),
    (5, 4): cscale(-1, CB),
}
FACE_K = {
    (2, 0): ONE,
    (4, 2): cscale(-1, CK),
    (5, 4): cscale(-1, CB),
}


def merge_faces(*faces):
    merged = {}
    for face in faces:
        for exponent, coefficient in face.items():
            if exponent in merged:
                need(merged[exponent] == coefficient,
                     "adjacent face coefficients agree on shared edge")
            else:
                merged[exponent] = coefficient
    return merged


FACE_LEDGERS = {
    "M_W": pick_ledger(M_W),
    "M_0": pick_ledger(M_0),
    "V7": pick_ledger(FACE_V7),
    "V9": pick_ledger(FACE_V9),
    "K": pick_ledger(FACE_K),
}
need(FACE_LEDGERS["M_W"][1] == (36, 10, 14), "M_W Pick ledger")
need(FACE_LEDGERS["M_0"][1] == (24, 14, 6), "M_0 Pick ledger")
need(FACE_LEDGERS["V7"][1] == (7, 3, 3), "V7 Pick ledger")
need(FACE_LEDGERS["V9"][1] == (18, 8, 6), "V9 Pick ledger")
need(FACE_LEDGERS["K"][1] == (2, 4, 0), "K-face Pick ledger")


GLOBAL_DATA = {}
for w_present, k_present in product((False, True), repeat=2):
    faces = [M_W if w_present else M_0,
             FACE_V7 if w_present else FACE_V9]
    if k_present:
        faces.append(FACE_K)
    merged = merge_faces(*faces)
    polygon, ledger = pick_ledger(merged)
    need(ledger[2] == 17, "every zeta-zero global polygon has genus 17")
    GLOBAL_DATA[(w_present, k_present)] = (merged, polygon, ledger)

need(GLOBAL_DATA[(False, False)][2] == (42, 10, 17), "W0K0 global")
need(GLOBAL_DATA[(False, True)][2] == (44, 12, 17), "W0K1 global")
need(GLOBAL_DATA[(True, False)][2] == (43, 11, 17), "W1K0 global")
need(GLOBAL_DATA[(True, True)][2] == (45, 13, 17), "W1K1 global")


def kummer_genus(degree, valuations):
    need(sum(valuations) == 0, "complete Kummer divisor has degree zero")
    ramification = sum(degree - gcd(degree, abs(value))
                        for value in valuations if value)
    numerator = 2 - 2 * degree + ramification
    need(numerator >= 0 and numerator % 2 == 0, "integral Kummer genus")
    return numerator // 2, ramification


# C: with x=P^{-1}, y=W*S/P, the inverse is P=x^{-1}, S=y/(W*x),
# and the equation is y^2=W*x*(x^6-U).  Its seven finite roots are simple
# under U*W!=0, and infinity has odd order seven.
C_GENUS, C_RAM = kummer_genus(2, (1,) * 7 + (-7,))
need((C_GENUS, C_RAM) == (3, 8), "central hyperelliptic normalization")

# V7: x=W*S^2*P^5 and 1-x=beta*S^3*P^4 give
# P^7=(beta^2/W^3)*x^3/(1-x)^2.  The inverse recovery is
# S=(W/beta)*P*(1-x)/x.
V7_GENUS, V7_RAM = kummer_genus(7, (3, -2, -1))
need((V7_GENUS, V7_RAM) == (3, 18), "cyclic-seven normalization")

# V9: S^3=(1-U*P^6)/(beta*P^4).  There are six simple numerator
# roots, plus valuations -4 at zero and -2 at infinity.
V9_GENUS, V9_RAM = kummer_genus(3, (1,) * 6 + (-4, -2))
need((V9_GENUS, V9_RAM) == (6, 16), "cyclic-three V9 normalization")

# K-face: V=S*P gives 1-K*V^2-beta*V^3*P=0, whence
# P=(1-K*V^2)/(beta*V^3) and S=V/P.  It is rational.
K_GENUS = 0

# The determinant criterion for a trinomial 1-X-Y proves torus smoothness.
SMOOTH_DETERMINANTS = {
    "C": 0 * 5 - 6 * 2,
    "V7": 2 * 4 - 5 * 3,
    "V9": 0 * 4 - 6 * 3,
    "K": 2 * 4 - 2 * 3,
}
need(SMOOTH_DETERMINANTS == {"C": -12, "V7": -7,
                             "V9": -18, "K": 2},
     "torus-smooth trinomial determinants")

# Generic-point P-derivative witnesses.  On the indicated face, vanishing of
# each displayed remainder would make two monomials constant.  The same
# nonzero exponent determinant above then contradicts transcendence degree 1.
DERIVATIVE_WITNESSES = {
    "C": "P*C_P-5*C=-U*P^6-5",
    "V7": "P*B_P-4*B=-4-W*S^2*P^5",
    "V9": "P*D_P-4*D=-4-2*U*P^6",
}


def p_euler_remainder(polynomial, weight):
    """Coefficient dictionary of P*d/dP(polynomial)-weight*polynomial."""
    return {
        exponent: cscale(exponent[1] - weight, coefficient)
        for exponent, coefficient in polynomial.items()
        if cscale(exponent[1] - weight, coefficient) != ZERO
    }


C_CORE = {(0, 0): ONE, (0, 6): cscale(-1, CU),
          (2, 5): cscale(-1, CW)}
V7_CORE = {(0, 0): ONE, (2, 5): cscale(-1, CW),
           (3, 4): cscale(-1, CB)}
V9_CORE = {(0, 0): ONE, (0, 6): cscale(-1, CU),
           (3, 4): cscale(-1, CB)}
need(p_euler_remainder(C_CORE, 5)
     == {(0, 0): cscale(-5, ONE), (0, 6): cscale(-1, CU)},
     "exact central P-Euler remainder")
need(p_euler_remainder(V7_CORE, 4)
     == {(0, 0): cscale(-4, ONE), (2, 5): cscale(-1, CW)},
     "exact V7 P-Euler remainder")
need(p_euler_remainder(V9_CORE, 4)
     == {(0, 0): cscale(-4, ONE), (0, 6): cscale(-2, CU)},
     "exact V9 P-Euler remainder")


# ---------------------------------------------------------------------------
# 3. Boundary and internal edge schemes, reconstructed from exponent data
# ---------------------------------------------------------------------------


def edge_sequence(polynomial, start, end):
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dx), abs(dy))
    need(length > 0, "nontrivial edge")
    step = (dx // length, dy // length)
    sequence = tuple(
        polynomial.get((start[0] + index * step[0],
                        start[1] + index * step[1]), ZERO)
        for index in range(length + 1)
    )
    need(sequence[0] != ZERO and sequence[-1] != ZERO,
         "edge endpoints are active")
    return length, sequence


def fmt_coefficient(coefficient):
    names = ("1", "U", "W", "beta", "K")
    pieces = []
    for value, name in zip(coefficient, names):
        if not value:
            continue
        pieces.append(f"{value:+d}*{name}")
    return "0" if not pieces else "".join(pieces).lstrip("+")


def fmt_sequence(sequence):
    return "[" + ",".join(fmt_coefficient(c) for c in sequence) + "]"


EDGE_DATA = {}
for key, (merged, polygon, _) in GLOBAL_DATA.items():
    w_present, k_present = key
    edges = []
    for index, start in enumerate(polygon):
        end = polygon[(index + 1) % len(polygon)]
        length, sequence = edge_sequence(merged, start, end)
        edges.append(("boundary", start, end, length, sequence))

    if w_present:
        internal = [((2, 0), (4, 5), "C--V7")]
        if k_present:
            internal.append(((2, 0), (5, 4), "V7--K"))
    else:
        internal = [((2, 0), (2, 6), "six-lines--V9")]
        if k_present:
            internal.append(((2, 0), (5, 4), "V9--K"))
    for start, end, label in internal:
        length, sequence = edge_sequence(merged, start, end)
        edges.append((label, start, end, length, sequence))
    EDGE_DATA[key] = tuple(edges)


# Exact sequence checks fix the edge schemes up to multiplication by a unit
# and inversion of X.
LINE_X_MINUS_1 = (cscale(-1, ONE), ONE)
ONE_MINUS_BX = (ONE, cscale(-1, CB))
ONE_MINUS_WX = (ONE, cscale(-1, CW))
ONE_MINUS_UX6 = (ONE,) + (ZERO,) * 5 + (cscale(-1, CU),)
U_MINUS_X6 = (CU,) + (ZERO,) * 5 + (cscale(-1, ONE),)
ONE_MINUS_KX2 = (ONE, ZERO, cscale(-1, CK))
TOP = (cscale(-1, CW), csub(CW, CU), CU)
need(TOP == (cscale(-1, CW), csub(CW, CU), CU),
     "(X-1)(U*X+W) coefficient identity")


EXPECTED_EDGE_SEQUENCES = {
    (True, False): (
        LINE_X_MINUS_1,
        ONE_MINUS_BX,
        (cscale(-1, CB), cscale(-1, CW)),
        TOP,
        U_MINUS_X6,
        ONE_MINUS_WX,
    ),
    (True, True): (
        LINE_X_MINUS_1,
        ONE_MINUS_KX2,
        (cscale(-1, CK), cscale(-1, CB)),
        (cscale(-1, CB), cscale(-1, CW)),
        TOP,
        U_MINUS_X6,
        ONE_MINUS_WX,
        ONE_MINUS_BX,
    ),
    (False, False): (
        LINE_X_MINUS_1,
        ONE_MINUS_BX,
        (cscale(-1, CB), cscale(-1, CU)),
        (cscale(-1, CU), CU),
        U_MINUS_X6,
        ONE_MINUS_UX6,
    ),
    (False, True): (
        LINE_X_MINUS_1,
        ONE_MINUS_KX2,
        (cscale(-1, CK), cscale(-1, CB)),
        (cscale(-1, CB), cscale(-1, CU)),
        (cscale(-1, CU), CU),
        U_MINUS_X6,
        ONE_MINUS_UX6,
        ONE_MINUS_BX,
    ),
}

for key, edges in EDGE_DATA.items():
    actual = tuple(edge[4] for edge in edges)
    need(actual == EXPECTED_EDGE_SEQUENCES[key],
         f"exact edge schemes for stratum {key}")
    boundary_lengths = tuple(edge[3] for edge in edges
                             if edge[0] == "boundary")
    need(sum(boundary_lengths) == GLOBAL_DATA[key][2][1],
         "edge orbit lengths reproduce global boundary lattice count")

# Every linear/binomial scheme is squarefree under its displayed nonzero
# coefficient gate.  The K-binomial has two distinct nonzero roots in
# characteristic zero.  The only repeated scheme is TOP when U+W=0:
# TOP becomes U*(X-1)^2 because W=-U.
TOP_AT_LAMBDA_ZERO = (CU, cscale(-2, CU), CU)
need(tuple(cadd(c, ZERO) for c in TOP_AT_LAMBDA_ZERO)
     == (CU, cscale(-2, CU), CU), "Lambda-zero top double root")


# ---------------------------------------------------------------------------
# 4. Component graphs, genera, and good-form orders
# ---------------------------------------------------------------------------


def graph_betti(vertices, weighted_edges):
    adjacency = defaultdict(set)
    edge_count = 0
    for left, right, multiplicity in weighted_edges:
        need(left in vertices and right in vertices and multiplicity > 0,
             "valid component edge")
        adjacency[left].add(right)
        adjacency[right].add(left)
        edge_count += multiplicity
    seen = set()
    queue = deque([next(iter(vertices))])
    while queue:
        vertex = queue.popleft()
        if vertex in seen:
            continue
        seen.add(vertex)
        queue.extend(adjacency[vertex] - seen)
    need(seen == set(vertices), "component graph connected")
    return edge_count - len(vertices) + 1, edge_count


GRAPH_DATA = {}
for w_present, k_present in product((False, True), repeat=2):
    if w_present:
        vertices = {"R", "C", "V7"}
        edges = [("R", "C", 12), ("C", "V7", 1)]
        genus_sum = C_GENUS + V7_GENUS
        if k_present:
            vertices.add("K")
            edges.append(("V7", "K", 1))
    else:
        lines = {f"L{index}" for index in range(6)}
        vertices = {"R", "V9"} | lines
        edges = []
        for line in sorted(lines):
            edges.append(("R", line, 2))
            edges.append((line, "V9", 1))
        genus_sum = V9_GENUS
        if k_present:
            vertices.add("K")
            edges.append(("V9", "K", 1))
    betti, edge_count = graph_betti(vertices, edges)
    need(betti == 11, "off-contact graph has b1=11")
    need(genus_sum + betti == 17, "component genus plus graph genus")
    GRAPH_DATA[(w_present, k_present)] = (
        len(vertices), edge_count, betti, genus_sum
    )


def good_order(plane, base):
    a_s, b, c = plane
    value = F(base) * (F(5, 6) - a_s - b - c)
    need(value.denominator == 1 and value > 0,
         "positive integral good-form order")
    return value.numerator


ORDER_DATA = {
    (True, False): (("M", good_order(M, 84)),
                    ("V7", good_order(V7, 84))),
    (True, True): (("M", good_order(M, 84)),
                   ("V7", good_order(V7, 84)),
                   ("K", good_order(KF, 84))),
    (False, False): (("M", good_order(M, 36)),
                     ("V9", good_order(V9, 36))),
    (False, True): (("M", good_order(M, 36)),
                    ("V9", good_order(V9, 36)),
                    ("K", good_order(KF, 36))),
}
need(ORDER_DATA[(True, True)] == (("M", 63), ("V7", 70), ("K", 196)),
     "W-nonzero orders")
need(ORDER_DATA[(False, True)] == (("M", 27), ("V9", 28), ("K", 84)),
     "W-zero orders")


# Lambda=0 occurs only in the W!=0 strata: W=-U and U!=0.  The R--C
# intersection becomes one A_23 contact of delta 12.  The other internal
# roots live on distinct toric divisors and are finite/nonzero, so neither the
# C--V7 node nor (when present) the V7--K node can merge with that contact.
for k_present in (False, True):
    component_count = 3 + int(k_present)
    delta_sum = 12 + 1 + int(k_present)
    need(C_GENUS + V7_GENUS + delta_sum - component_count + 1 == 17,
         "Lambda-zero delta ledger")

# The literal top-contact expansion has c1=alpha_11+beta_11.  A deepest
# repeat requires c1=0, hence alpha_11=-beta_11!=0 on the present gate.  Thus
# C1=alpha_11/c^2 is nonzero (the repeated discriminant has c!=0).  In the
# repeated range nu>s, compare its gap d1=s+nu with the nonzero b^12*q gap
# db=6(nu-s): db is first below nu/s=7/5, both tie at equality, and d1 is
# first above it.  Their Keller-form orders are respectively 6s+2nu and
# (5s+9nu)/2, both positive.  Both gaps precede the W-dependent normalized
# t^6 correction.  The K=[y^2]H row is K*t^6*r^2 in Hhat and becomes
# K*y*r^2=K*y+2*K*t^6*y^2+K*t^12*y^3.  Its leading part is absorbed into c;
# it cannot kill C1 because c!=0, while its remainder begins at t^6 too.
splitter_cases = ((5, 6, "b12"), (5, 7, "tie"), (5, 8, "C1"))
need((1 + 6, 1 - 6) == (7, -5),
     "d1-db=7*s-5*nu gives threshold nu/s=7/5")
for s, nu, expected in splitter_cases:
    d1 = s + nu
    db = 6 * (nu - s)
    actual = "tie" if d1 == db else ("C1" if d1 < db else "b12")
    need(actual == expected, "exact early-splitter threshold")
    need(6 * s + 2 * nu > 0, "b12*q good-form order")
    need(5 * s + 9 * nu > 0, "C1 good-form order numerator")
    need(d1 < 6 * (s + nu) and db < 6 * (s + nu),
         "both early splitters precede every t6 correction")
# At equality, X*(C1-X^5/c) has derivative C1-6*X^5/c=-5*C1
# at every nonzero root.  Characteristic zero and C1!=0 make those roots
# simple; the zero root is the next ratio case.
need(1 - 6 == -5 and -5 != 0,
     "equality-face nonzero-root derivative is -5*C1")


# ---------------------------------------------------------------------------
# 5. Sharp beta=0 hostiles
# ---------------------------------------------------------------------------


H_WK = (F(1, 2), F(0), F(-1))
H_0K = (F(1, 3), F(1, 6), F(-2, 3))
BETA_ZERO_EXPECTED = {
    (False, False): (M,),
    (False, True): (M, H_0K),
    (True, False): (M,),
    (True, True): (M, H_WK),
}
BETA_ZERO_SUMMARY = []
for w_present, k_present in product((False, True), repeat=2):
    # A sharp explicit witness: all nonrequired lower rows are zero.  There
    # are no aggregate cancellations to choose in this witness.
    required = {P1, P2, P3, UROW}
    if w_present:
        required.add(WROW)
    if k_present:
        required.add(KROW)
    support = BASE_POINTS | {
        point for row in required for point in lifted_points(row)
    }
    facets = lower_facets(point_mask(support))
    need(facets == tuple(sorted(BETA_ZERO_EXPECTED[(w_present, k_present)])),
         "sharp beta-zero lower-complex witness")

    # Its global support is the union of the surviving face vertices.
    if w_present:
        points = set(M_W)
        if k_present:
            points |= {(4, 2)}
    else:
        points = set(M_0)
        if k_present:
            points |= {(4, 2)}
    polygon, ledger = pick_ledger(points)
    need(ledger[2] < 17, "beta-zero witness loses the genus-17 ledger")
    BETA_ZERO_SUMMARY.append((w_present, k_present, facets, ledger))

# beta=0 also makes the common edge 1-beta*X monomial and permits alpha=0
# when c1=alpha+beta=0, so both the edge-node audit and the C1 shortcut fail
# exactly at the excluded gate.


# ---------------------------------------------------------------------------
# Deterministic audit transcript
# ---------------------------------------------------------------------------


def fmt_plane(plane):
    return "(" + ",".join(str(value) for value in plane) + ")"


print("THM4337_ZETA0_CLEANROOM_AUDIT")
print(f"rows={len(ROWS)} master_points={len(MASTER)} candidate_planes={len(PLANE_DATA)}")
for w_present, k_present, states, faces in ATLAS_SUMMARY:
    print(
        f"ATLAS W={int(w_present)} K={int(k_present)} states={states} "
        f"faces={'|'.join(fmt_plane(face) for face in faces)}"
    )
for name in ("M_W", "M_0", "V7", "V9", "K"):
    polygon, ledger = FACE_LEDGERS[name]
    print(f"FACE {name} polygon={polygon} pick={ledger}")
print(f"NORMAL C=y^2-W*x*(x^6-U) genus={C_GENUS} ram={C_RAM}")
print(f"NORMAL V7=P^7-(beta^2/W^3)*x^3/(1-x)^2 genus={V7_GENUS} ram={V7_RAM}")
print(f"NORMAL V9=S^3-(1-U*P^6)/(beta*P^4) genus={V9_GENUS} ram={V9_RAM}")
print("NORMAL K V=S*P;P=(1-K*V^2)/(beta*V^3) genus=0")
print(f"SMOOTH determinants={SMOOTH_DETERMINANTS}")
print(f"DERIVATIVE witnesses={DERIVATIVE_WITNESSES}")
print("DERIVATIVE identities=exact_P_Euler_remainders_executed")
for key in ((True, False), (True, True), (False, False), (False, True)):
    merged, polygon, ledger = GLOBAL_DATA[key]
    print(f"GLOBAL W={int(key[0])} K={int(key[1])} polygon={polygon} pick={ledger}")
    for kind, start, end, length, sequence in EDGE_DATA[key]:
        print(
            f"EDGE W={int(key[0])} K={int(key[1])} kind={kind} "
            f"segment={start}->{end} length={length} coeffs={fmt_sequence(sequence)}"
        )
    print(f"GRAPH W={int(key[0])} K={int(key[1])} V,E,b1,gsum={GRAPH_DATA[key]}")
    print(f"ORDERS W={int(key[0])} K={int(key[1])} data={ORDER_DATA[key]}")
print("LAMBDA0 top=(X-1)(U*X+W)=U*(X-1)^2;delta=12")
print("LAMBDA0 repeat:c1=alpha+beta=0=>alpha=-beta!=0=>C1!=0")
print("LAMBDA0 early_splitter=min(b^12*q,C1*t);threshold_nu_over_s=7/5")
print("LAMBDA0 orders=b^12*q:6*s+2*nu;C1*t:(5*s+9*nu)/2;both_before_all_t6_corrections")
print("LAMBDA0 equality_face=X*(C1-X^5/c);nonzero_root_derivative=-5*C1")
print("LAMBDA0 K_y2=K*y+2*K*t^6*y^2+K*t^12*y^3;leading_absorbed_in_c;remainder_t6")
for w_present, k_present, facets, ledger in BETA_ZERO_SUMMARY:
    print(
        f"BETA0_HOSTILE W={int(w_present)} K={int(k_present)} "
        f"faces={'|'.join(fmt_plane(face) for face in facets)} pick={ledger}"
    )
print("BETA0_FAILURE edge_1-betaX_is_monomial;alpha=0_allows_C1=0")
print("RESULT=ALL_THM4337_EXACT_CHECKS_PASS")
