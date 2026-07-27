#!/usr/bin/env python3
"""Exact finite audit for THM-2620.

The calculation is dependency-free.  It enumerates every endpoint matrix over
F_13 and independently checks the translation quotient, determinant fibres,
projective transvections, chronological gluing, homogeneity, and the intrinsic
row/column V4 action.
"""

from collections import Counter, defaultdict
from itertools import product


P = 13
V = tuple(product(range(P), repeat=2))
ZERO = (0, 0)
I = (1, 0, 0, 1)
MINUS_I = (P - 1, 0, 0, P - 1)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def add(x, y):
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def neg(x):
    return ((-x[0]) % P, (-x[1]) % P)


def sub(x, y):
    return add(x, neg(y))


def scale(c, x):
    return ((c * x[0]) % P, (c * x[1]) % P)


def swap(x):
    return (x[1], x[0])


def det(x, y):
    return (x[0] * y[1] - x[1] * y[0]) % P


def dot(x, y):
    return (x[0] * y[0] + x[1] * y[1]) % P


def inv(a):
    require(a % P != 0, "attempted to invert zero")
    return pow(a, P - 2, P)


def mat(a, b, c, d):
    """Row-major 2 by 2 matrix over F_13."""
    return (a % P, b % P, c % P, d % P)


def mmul(a, b):
    aa, ab, ac, ad = a
    ba, bb, bc, bd = b
    return mat(aa * ba + ab * bc, aa * bb + ab * bd,
               ac * ba + ad * bc, ac * bb + ad * bd)


def mdet(a):
    return (a[0] * a[3] - a[1] * a[2]) % P


def minv(a):
    d = inv(mdet(a))
    return mat(d * a[3], -d * a[1], -d * a[2], d * a[0])


def mpow(a, n):
    require(n >= 0, "negative matrix power")
    ans = I
    base = a
    while n:
        if n & 1:
            ans = mmul(ans, base)
        base = mmul(base, base)
        n //= 2
    return ans


def mneg(a):
    return tuple((-x) % P for x in a)


def psl_key(a):
    require(mdet(a) == 1, "PSL representative is not in SL2")
    return min(a, mneg(a))


def columns(left, right):
    """Matrix whose first and second columns are left and right."""
    return mat(left[0], right[0], left[1], right[1])


def mat_key(a):
    """Serialize a linear map by its two image columns."""
    return ((a[0], a[2]), (a[1], a[3]))


def unipotent(a=1):
    return mat(1, a, 0, 1)


def owner(t):
    return mat(0, 1, -1, t)


def conjugate(a, b):
    return mmul(mmul(a, b), minv(a))


def mproduct(items):
    ans = I
    for item in items:
        ans = mmul(ans, item)
    return ans


def point(x):
    """Canonical representative of [x] in P^1(F_13)."""
    require(x != ZERO, "zero vector has no projective point")
    if x[0] != 0:
        z = inv(x[0])
    else:
        z = inv(x[1])
    return scale(z, x)


PROJECTIVE = tuple(sorted({point(x) for x in V if x != ZERO}))


def phi(q, delta, x):
    return det(q, x) * inv(delta) % P


def transvection(q, delta, x, power=1):
    """(I + power*q*phi)(x), with phi(x)=det(q,x)/delta."""
    return add(x, scale(power * phi(q, delta, x) % P, q))


def projective_transvection(q, delta, x, power=1):
    return point(transvection(q, delta, x, power))


def endpoint_fibre(q, delta):
    ans = []
    for right in V:
        left = add(right, q)
        if det(left, right) == delta:
            ans.append((left, right))
    return tuple(ans)


def matrix_key(q, delta):
    """The full linear action of T_(q,delta), serialized on the basis."""
    return (transvection(q, delta, (1, 0)),
            transvection(q, delta, (0, 1)))


checks = 0


def checked(condition, message):
    global checks
    require(condition, message)
    checks += 1


# The fourteen-point projective line.
checked(len(PROJECTIVE) == P + 1, "projective-line cardinality changed")

# Complete endpoint-matrix census and the row-difference/determinant identity.
matrix_census = Counter()
for left in V:
    for right in V:
        q = sub(left, right)
        delta = det(left, right)
        checked(delta == det(q, right), "det(L,R)=det(q,R) failed")
        if q == ZERO:
            checked(delta == 0, "equal endpoints had nonzero determinant")
            matrix_census["q_zero"] += 1
        elif delta == 0:
            matrix_census["q_nonzero_delta_zero"] += 1
        else:
            matrix_census["q_nonzero_delta_nonzero"] += 1

expected_census = Counter({
    "q_zero": P * P,
    "q_nonzero_delta_zero": (P * P - 1) * P,
    "q_nonzero_delta_nonzero": (P * P - 1) * P * (P - 1),
})
checked(matrix_census == expected_census, "endpoint-matrix census changed")
checked(sum(matrix_census.values()) == P ** 4, "matrix census is incomplete")

# Common row translation fixes q and makes Delta exactly uniform for q!=0.
translation_histograms = set()
for q in V:
    if q == ZERO:
        continue
    right = (0, 0)
    left = q
    delta = det(left, right)
    hist = Counter()
    for t in V:
        moved_left = add(left, t)
        moved_right = add(right, t)
        checked(sub(moved_left, moved_right) == q,
                "common translation changed row difference")
        moved_delta = det(moved_left, moved_right)
        checked(moved_delta == (delta + det(q, t)) % P,
                "translation determinant law failed")
        hist[moved_delta] += 1
    checked(hist == Counter({a: P for a in range(P)}),
            "translation did not make determinant uniform")
    translation_histograms.add(tuple(sorted(hist.items())))
checked(len(translation_histograms) == 1,
        "translation histogram depends on the nonzero target vector")

# The 169 difference characters are only the diagonal slice
# (alpha,-alpha) of the 13^4 joint endpoint character bank.  A two-point
# hostile has the same difference and every diagonal character, but different
# determinant and an off-diagonal character witness.
hostile_q = (1, 0)
hostile_t = (0, 1)
hostile_zero = (hostile_q, ZERO)
hostile_one = (add(hostile_q, hostile_t), hostile_t)
checked(sub(*hostile_zero) == hostile_q == sub(*hostile_one),
        "joint-endpoint hostile changed target difference")
checked(det(*hostile_zero) == 0 and det(*hostile_one) == 1,
        "joint-endpoint hostile did not change determinant")
diagonal_character_matches = 0
for alpha in V:
    exponent_zero = (dot(alpha, hostile_zero[0])
                     - dot(alpha, hostile_zero[1])) % P
    exponent_one = (dot(alpha, hostile_one[0])
                    - dot(alpha, hostile_one[1])) % P
    checked(exponent_zero == exponent_one,
            "difference character detected a common endpoint translation")
    diagonal_character_matches += 1
offdiag_left = (0, 1)
offdiag_right = ZERO
offdiag_zero = (dot(offdiag_left, hostile_zero[0])
                + dot(offdiag_right, hostile_zero[1])) % P
offdiag_one = (dot(offdiag_left, hostile_one[0])
               + dot(offdiag_right, hostile_one[1])) % P
checked(offdiag_zero != offdiag_one,
        "off-diagonal joint character missed the endpoint hostile")
checked(diagonal_character_matches == P * P,
        "diagonal character census changed")

# Every fixed (q,Delta) fibre has thirteen points.  Nonzero Delta is exactly
# one parabolic graph on P^1 minus [q]; Delta=0 is the collapsed boundary.
nondegenerate_edges = 0
for q in V:
    if q == ZERO:
        continue
    q_point = point(q)
    affine_points = tuple(x for x in PROJECTIVE if x != q_point)
    checked(len(affine_points) == P, "projective complement is not affine")
    for delta in range(P):
        fibre = endpoint_fibre(q, delta)
        checked(len(fibre) == P, "fixed target/determinant fibre changed size")
        if delta == 0:
            for left, right in fibre:
                checked(right in tuple(scale(c, q) for c in range(P)),
                        "degenerate right endpoint left the q-line")
                checked(left in tuple(scale(c, q) for c in range(P)),
                        "degenerate left endpoint left the q-line")
            continue

        nondegenerate_edges += len(fibre)
        by_line = {}
        edge_graph = set()
        for left, right in fibre:
            checked(left != ZERO and right != ZERO,
                    "nonzero determinant fibre contains a zero endpoint")
            checked(point(left) != q_point and point(right) != q_point,
                    "nonzero determinant fibre hit the fixed boundary point")
            checked(det(q, right) == delta and det(q, left) == delta,
                    "endpoint normalization changed along one edge")
            right_line = point(right)
            checked(right_line not in by_line,
                    "two normalized representatives occupy one line")
            by_line[right_line] = right
            edge_graph.add((right_line, point(left)))
            checked(transvection(q, delta, right) == left,
                    "transvection did not carry right endpoint to left")

        checked(set(by_line) == set(affine_points),
                "fixed determinant missed an affine projective point")
        expected_graph = {
            (x, projective_transvection(q, delta, x)) for x in affine_points
        }
        checked(edge_graph == expected_graph,
                "endpoint fibre is not the transvection graph")

        # Nilpotence gives T^m=I+mN.  Projectively, [q] is the unique fixed
        # point and the other thirteen points form one 13-cycle.
        checked(all(transvection(q, delta, x, P) == x for x in V),
                "transvection did not have order dividing thirteen")
        for m in range(1, P):
            checked(any(transvection(q, delta, x, m) != x for x in V),
                    "nontrivial transvection power became the identity")
        fixed = {
            x for x in PROJECTIVE
            if projective_transvection(q, delta, x) == x
        }
        checked(fixed == {q_point}, "transvection fixed-point set changed")
        orbit = []
        x = affine_points[0]
        for _ in range(P):
            checked(x not in orbit, "affine transvection orbit closed early")
            orbit.append(x)
            x = projective_transvection(q, delta, x)
        checked(x == orbit[0] and set(orbit) == set(affine_points),
                "affine transvection orbit is not one thirteen-cycle")

        # Fixed q and Delta make projective gluing equal exact vector gluing:
        # the determinant normalization chooses one representative per line.
        for left, right in fibre:
            next_right = by_line[point(left)]
            checked(next_right == left,
                    "projective adjacency failed exact chronological gluing")

checked(nondegenerate_edges == (P * P - 1) * (P - 1) * P,
        "nondegenerate edge total changed")

# Homogeneity and the exact count of nontrivial transvections.  Simultaneous
# scaling (q,Delta)->(c q,c^2 Delta) leaves the linear map unchanged.
transvection_parameters = Counter()
transvection_representative = {}
transvection_square_classes = defaultdict(set)
for q in V:
    if q == ZERO:
        continue
    for delta in range(1, P):
        key = matrix_key(q, delta)
        transvection_parameters[key] += 1
        transvection_representative.setdefault(key, (q, delta))
        transvection_square_classes[key].add(pow(delta, (P - 1) // 2, P))
        for c in range(1, P):
            checked(matrix_key(scale(c, q), c * c * delta % P) == key,
                    "projective endpoint scaling changed the transvection")
checked(len(transvection_parameters) == P * P - 1,
        "nontrivial transvection count changed")
checked(set(transvection_parameters.values()) == {P - 1},
        "transvection parameter fibres are not scalar orbits")
checked(all(len(classes) == 1 for classes in transvection_square_classes.values()),
        "one transvection crossed determinant square classes")
square_class_hist = Counter(next(iter(classes))
                            for classes in transvection_square_classes.values())
checked(square_class_hist == Counter({1: 84, P - 1: 84}),
        "PSL transvection classes did not split 84+84")

# The endpoint fibre is the column-matrix conjugacy B U B^-1.  Modulo common
# scalar, a nondegenerate B is a PGL2 element, and retaining [R] makes the map
# to a transvection plus a nonfixed marked point bijective.  The determinant-
# square half is PSL2 and selects one of the two 84-element transvection classes.
pointed_pgl = Counter()
pointed_psl = Counter()
marks_by_transvection = defaultdict(set)
for q in V:
    if q == ZERO:
        continue
    for right in V:
        delta = det(q, right)
        if delta == 0:
            continue
        b = columns(q, right)
        t_matrix = mmul(mmul(b, unipotent()), minv(b))
        key = matrix_key(q, delta)
        checked(mat_key(t_matrix) == key, "B U B^-1 disagrees with T_(q,Delta)")
        marked = point(right)
        checked(marked != point(q), "pointed transvection used its fixed point")
        pointed_pgl[(key, marked)] += 1
        marks_by_transvection[key].add(marked)
        if pow(delta, (P - 1) // 2, P) == 1:
            pointed_psl[(key, marked)] += 1

checked(len(pointed_pgl) == 2184 and set(pointed_pgl.values()) == {P - 1},
        "PGL pointed-transvection bijection failed")
checked(len(pointed_psl) == 1092 and set(pointed_psl.values()) == {P - 1},
        "PSL pointed-transvection bijection failed")
checked(all(len(marks) == P for marks in marks_by_transvection.values()),
        "a transvection did not have thirteen marked nonfixed points")

# On PSL2, the order-seven owner acts on the left and the target deck on the
# right.  These actions commute as actions and their joint C7 x C13 action is
# free, giving twelve 91-point torsor orbits without an internal order-91
# element.
sl2 = {
    mat(a, b, c, d)
    for a, b, c, d in product(range(P), repeat=4)
    if (a * d - b * c) % P == 1
}
psl2 = {psl_key(a) for a in sl2}
checked(len(sl2) == 2184 and len(psl2) == 1092,
        "SL2/PSL2 group census changed")
g_owner = mat(7, 5, 10, 11)
u_deck = unipotent()
checked(mpow(g_owner, 7) == MINUS_I and mpow(u_deck, P) == I,
        "owner/deck orders changed")
remaining = set(psl2)
owner_deck_orbits = []
while remaining:
    g = next(iter(remaining))
    orbit = {
        psl_key(mmul(mmul(mpow(g_owner, i), g), mpow(u_deck, j)))
        for i in range(7) for j in range(P)
    }
    checked(len(orbit) == 91, "left-owner/right-deck action was not free")
    owner_deck_orbits.append(orbit)
    remaining -= orbit
checked(len(owner_deck_orbits) == 12 and not remaining,
        "PSL2 did not split into twelve C7 x C13 torsor orbits")

# The row/column swaps commute and form an intrinsic V4.  Row reversal inverts
# the transvection; target-column swap conjugates it by the coordinate swap.
for left in V:
    for right in V:
        q = sub(left, right)
        delta = det(left, right)
        row_left, row_right = right, left
        col_left, col_right = swap(left), swap(right)
        both_left, both_right = swap(right), swap(left)
        checked(sub(row_left, row_right) == neg(q)
                and det(row_left, row_right) == (-delta) % P,
                "row reversal sign law failed")
        checked(sub(col_left, col_right) == swap(q)
                and det(col_left, col_right) == (-delta) % P,
                "column swap sign law failed")
        checked(sub(both_left, both_right) == neg(swap(q))
                and det(both_left, both_right) == delta,
                "double-swap sign law failed")
        checked((both_left, both_right)
                == (swap(row_left), swap(row_right)),
                "row and column swaps failed to commute")

        if q != ZERO and delta != 0:
            for x in ((1, 0), (0, 1), add(left, right)):
                checked(transvection(neg(q), (-delta) % P, x)
                        == transvection(q, delta, x, -1),
                        "endpoint reversal did not invert the transvection")
                checked(transvection(swap(q), (-delta) % P, swap(x))
                        == swap(transvection(q, delta, x)),
                        "target swap did not conjugate the transvection")


def orbit_histogram(universe, operations):
    universe = set(universe)
    remaining = set(universe)
    histogram = Counter()
    while remaining:
        x = next(iter(remaining))
        orbit = {operation(x) for operation in operations}
        checked(orbit <= universe, "V4 operation left its declared universe")
        histogram[len(orbit)] += 1
        remaining -= orbit
    fixed = tuple(sum(operation(x) == x for x in universe)
                  for operation in operations)
    return fixed, histogram


endpoint_operations = (
    lambda pair: pair,
    lambda pair: (pair[1], pair[0]),
    lambda pair: (swap(pair[0]), swap(pair[1])),
    lambda pair: (swap(pair[1]), swap(pair[0])),
)
all_endpoint_pairs = {(left, right) for left in V for right in V}
nondegenerate_pairs = {
    pair for pair in all_endpoint_pairs
    if sub(*pair) != ZERO and det(*pair) != 0
}
all_v4_fixed, all_v4_orbits = orbit_histogram(
    all_endpoint_pairs, endpoint_operations)
nondeg_v4_fixed, nondeg_v4_orbits = orbit_histogram(
    nondegenerate_pairs, endpoint_operations)
checked(all_v4_fixed == (28561, 169, 169, 169),
        "all-matrix V4 Burnside counts changed")
checked(all_v4_orbits == Counter({1: 13, 2: 234, 4: 7020}),
        "all-matrix V4 orbit sizes changed")
checked(nondeg_v4_fixed == (26208, 0, 0, 144),
        "nondegenerate V4 Burnside counts changed")
checked(nondeg_v4_orbits == Counter({2: 72, 4: 6516}),
        "nondegenerate V4 orbit sizes changed")

parameter_operations = (
    lambda item: item,
    lambda item: (neg(item[0]), (-item[1]) % P),
    lambda item: (swap(item[0]), (-item[1]) % P),
    lambda item: (neg(swap(item[0])), item[1]),
)
parameter_pairs = {
    (q, delta) for q in V if q != ZERO for delta in range(1, P)
}
parameter_v4_fixed, parameter_v4_orbits = orbit_histogram(
    parameter_pairs, parameter_operations)
checked(parameter_v4_fixed == (2016, 0, 0, 144),
        "parameter-pair V4 Burnside counts changed")
checked(parameter_v4_orbits == Counter({2: 72, 4: 468}),
        "parameter-pair V4 orbit sizes changed")


def key_operation(which, key):
    q, delta = transvection_representative[key]
    return matrix_key(*parameter_operations[which]((q, delta)))


transvection_operations = tuple(
    (lambda which: lambda key: key_operation(which, key))(which)
    for which in range(4)
)
transvection_keys = set(transvection_parameters)
transvection_v4_fixed, transvection_v4_orbits = orbit_histogram(
    transvection_keys, transvection_operations)
checked(transvection_v4_fixed == (168, 0, 0, 24),
        "transvection V4 Burnside counts changed")
checked(transvection_v4_orbits == Counter({2: 12, 4: 36}),
        "transvection V4 orbit sizes changed")

# Signed SL2 refinement.  Projective order-seven closure may be +I or -I;
# the trace of the moving matrix determines the central sign.  This is the
# exact C2 scale invoice lost on passage to PSL2.
signed_expected = {
    (3, "F"): ({10, 11}, {6, 8, 9}),
    (3, "R"): ({2, 3}, {4, 5, 7}),
    (5, "F"): ({2, 12}, {8, 10, 11}),
    (5, "R"): ({1, 11}, {2, 3, 5}),
    (6, "F"): ({1, 3}, {9, 11, 12}),
    (6, "R"): ({10, 12}, {1, 2, 4}),
}
signed_closures = {}
for t in (3, 5, 6):
    g = owner(t)
    checked(mpow(g, 7) == MINUS_I, f"g_{t}^7 changed sign")
    for orientation in ("F", "R"):
        plus, minus = set(), set()
        for a in range(1, P):
            factors = [conjugate(mpow(g, k), unipotent(a)) for k in range(7)]
            norm = mproduct(factors if orientation == "F"
                            else tuple(reversed(factors)))
            if norm == I:
                plus.add(a)
            elif norm == MINUS_I:
                minus.add(a)
            else:
                continue
            moving_trace = (t - a) % P if orientation == "F" else (t + a) % P
            checked(moving_trace in {3, 5, 6, 7, 8, 10},
                    "central norm had a non-order-seven moving trace")
        checked((plus, minus) == signed_expected[t, orientation],
                "signed six-state closure atlas changed")
        signed_closures[t, orientation] = (plus, minus)

canonical_factors = [
    conjugate(mpow(g_owner, k), unipotent()) for k in range(7)
]
checked(mproduct(tuple(reversed(canonical_factors))) == MINUS_I,
        "canonical trace-two endpoint lift lost its -I defect")
checked(mproduct(tuple(reversed([mneg(a) for a in canonical_factors]))) == I,
        "actual negative-unipotent THM-2603 lift did not close at +I")

# A fixed nondegenerate chronological edge iterates by q.  Seven edges do not
# close, while thirteen do: this is an obstruction, not an LRC trivialization.
q_control = (1, 0)
delta_control = 1
right_control = endpoint_fibre(q_control, delta_control)[0][1]
walk = [right_control]
for _ in range(P):
    walk.append(transvection(q_control, delta_control, walk[-1]))
checked(walk[7] == add(right_control, scale(7, q_control)),
        "seven-edge displacement changed")
checked(point(walk[7]) != point(right_control),
        "seven-edge parabolic walk closed unexpectedly")
checked(walk[P] == right_control, "thirteen-edge parabolic walk did not close")

print("== THM-2620 endpoint-pair parabolic transvection ==")
print(f"field: F_{P}")
print(f"projective line: {len(PROJECTIVE)} = {P} affine points + 1 fixed boundary")
print("endpoint matrix census:", dict(matrix_census))
print("common-translation determinant histogram:",
      dict(Counter({a: P for a in range(P)})))
print(f"difference DFT: {P * P} diagonal characters inside "
      f"{P ** 4} joint endpoint characters")
print("same-q hostile: all 169 diagonal characters agree; "
      "one off-diagonal character separates Delta=0 from Delta=1")
print(f"fixed (q,Delta) fibre size: {P}")
print(f"nondegenerate endpoint edges: {nondegenerate_edges}")
print(f"distinct nontrivial transvections: {len(transvection_parameters)}"
      f" = {P + 1} fixed points x {P - 1} powers")
print(f"parameter multiplicity per transvection: {P - 1}")
print("pointed bundles: PGL2=2184=168*13; PSL2 square class=1092=84*13")
print("PSL2 left-C7/right-C13 carrier: 12 free orbits of size 91")
print("V4 fixed counts all/nondegenerate/parameters/transvections:",
      all_v4_fixed, nondeg_v4_fixed, parameter_v4_fixed,
      transvection_v4_fixed)
print("V4 orbit histograms:", dict(all_v4_orbits),
      dict(nondeg_v4_orbits), dict(parameter_v4_orbits),
      dict(transvection_v4_orbits))
for t in (3, 5, 6):
    for orientation in ("F", "R"):
        plus, minus = signed_closures[t, orientation]
        print(f"signed closure {t}{orientation}: +I={sorted(plus)} "
              f"-I={sorted(minus)}")
print("central C2 invoice: canonical +U lift=-I; inherited -U lift=+I")
print("chronology: fixed q,Delta projective gluing equals vector gluing")
print("parabolic hostile: seven steps nonclosing; thirteen steps closing")
print(f"exact assertions passed: {checks}")
