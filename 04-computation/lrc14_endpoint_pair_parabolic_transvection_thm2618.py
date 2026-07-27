#!/usr/bin/env python3
"""Exact finite audit for THM-2618.

The calculation is dependency-free.  It enumerates every endpoint matrix over
F_13 and independently checks the translation quotient, determinant fibres,
projective transvections, chronological gluing, homogeneity, and the intrinsic
row/column V4 action.
"""

from collections import Counter
from itertools import product


P = 13
V = tuple(product(range(P), repeat=2))
ZERO = (0, 0)


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
for q in V:
    if q == ZERO:
        continue
    for delta in range(1, P):
        key = matrix_key(q, delta)
        transvection_parameters[key] += 1
        for c in range(1, P):
            checked(matrix_key(scale(c, q), c * c * delta % P) == key,
                    "projective endpoint scaling changed the transvection")
checked(len(transvection_parameters) == P * P - 1,
        "nontrivial transvection count changed")
checked(set(transvection_parameters.values()) == {P - 1},
        "transvection parameter fibres are not scalar orbits")

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

print("== THM-2618 endpoint-pair parabolic transvection ==")
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
print("V4 endpoint symmetries: row reversal, target swap, double swap PASS")
print("chronology: fixed q,Delta projective gluing equals vector gluing")
print("parabolic hostile: seven steps nonclosing; thirteen steps closing")
print(f"exact assertions passed: {checks}")
