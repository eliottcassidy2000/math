#!/usr/bin/env python3
"""Exact referee for THM-2626.

The script uses only integer arithmetic in F_13.  It builds PSL_2(F_13),
the 84-state P_13-coset frame bundle, all ordered seven-edge norms, and the
affine tangent-frame lifts.  Every truth-bearing check uses ``require`` so
that ``python -O`` performs the same verification.
"""

from itertools import combinations


P = 13
INF = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def raw(entries):
    return tuple(int(x) % P for x in entries)


I_RAW = raw((1, 0, 0, 1))
MINUS_I_RAW = raw((-1, 0, 0, -1))


def raw_mul(x, y):
    a, b, c, d = x
    e, f, g, h = y
    return raw((a * e + b * g, a * f + b * h,
                c * e + d * g, c * f + d * h))


def raw_det(x):
    a, b, c, d = x
    return (a * d - b * c) % P


def raw_pow(x, n):
    answer = I_RAW
    base = x
    while n:
        if n & 1:
            answer = raw_mul(answer, base)
        base = raw_mul(base, base)
        n //= 2
    return answer


def raw_inv(x):
    a, b, c, d = x
    require((a * d - b * c) % P == 1, "inverse requires determinant one")
    return raw((d, -b, -c, a))


def is_scalar(x):
    a, b, c, d = x
    return b == 0 and c == 0 and a == d


def canonical(x):
    x = raw(x)
    minus_x = raw((-x[0], -x[1], -x[2], -x[3]))
    return min(x, minus_x)


def psl_mul(x, y):
    return canonical(raw_mul(x, y))


def psl_inv(x):
    return canonical(raw_inv(x))


I = canonical(I_RAW)


def U(a):
    return raw((1, a, 0, 1))


def g_matrix(t):
    return raw((0, 1, -1, t))


def projective_action(matrix, x):
    a, b, c, d = matrix
    if x == INF:
        if c == 0:
            return INF
        return a * pow(c, -1, P) % P
    denominator = (c * x + d) % P
    if denominator == 0:
        return INF
    return (a * x + b) * pow(denominator, -1, P) % P


def projective_cycles(matrix):
    unseen = set(range(P)) | {INF}
    cycles = []
    while unseen:
        start = min(unseen)
        cycle = []
        x = start
        while x not in cycle:
            cycle.append(x)
            unseen.discard(x)
            x = projective_action(matrix, x)
        require(x == start, "projective orbit did not close at its start")
        cycles.append(tuple(cycle))
    return cycles


def rotate_to(cycle, start):
    j = cycle.index(start)
    return tuple(cycle[j:] + cycle[:j])


def derivative_factor(matrix, x):
    require(x != INF, "affine derivative requested at infinity")
    c, d = matrix[2], matrix[3]
    denominator = (c * x + d) % P
    require(denominator != 0, "affine derivative requested at a pole")
    return pow(denominator * denominator % P, -1, P)


# Build SL_2(F_13) and its projective quotient.
sl2 = []
psl2 = set()
for a in range(P):
    for b in range(P):
        for c in range(P):
            for d in range(P):
                if (a * d - b * c) % P == 1:
                    matrix = raw((a, b, c, d))
                    sl2.append(matrix)
                    psl2.add(canonical(matrix))

require(len(sl2) == 2184, "wrong SL2(F13) order")
require(len(psl2) == 1092, "wrong PSL2(F13) order")

S_RAW = g_matrix(0)
R_RAW = g_matrix(1)
S = canonical(S_RAW)
R = canonical(R_RAW)

require(canonical(raw_pow(S_RAW, 2)) == I, "S does not have projective order two")
require(canonical(raw_pow(R_RAW, 3)) == I, "R does not have projective order three")
require(canonical(raw_mul(raw_inv(R_RAW), S_RAW)) == canonical(U(1)),
        "U != R^-1 S projectively")

for t in (3, 5, 6):
    g = g_matrix(t)
    sr_power = raw_pow(raw_mul(S_RAW, R_RAW), t)
    require(canonical(g) == canonical(raw_mul(S_RAW, sr_power)),
            "g_t != S(SR)^t projectively")
    require(raw_pow(g, 7) == MINUS_I_RAW, "g_t^7 != -I")


# P_13, its Borel normalizer, and the 84-state frame bundle.
p13 = {canonical(U(a)) for a in range(P)}
borel = set()
for r in range(1, P):
    r_inv = pow(r, -1, P)
    for b in range(P):
        borel.add(canonical(raw((r, b, 0, r_inv))))

require(len(p13) == 13, "wrong P13 order")
require(len(borel) == 78, "wrong Borel order")


def left_cosets(group, subgroup):
    remaining = set(group)
    answer = []
    while remaining:
        representative = min(remaining)
        coset = frozenset(psl_mul(representative, h) for h in subgroup)
        require(len(coset) == len(subgroup), "bad left coset size")
        answer.append(coset)
        remaining.difference_update(coset)
    return answer


p13_cosets = left_cosets(psl2, p13)
borel_cosets = left_cosets(psl2, borel)
require(len(p13_cosets) == 84, "wrong number of P13 cosets")
require(len(borel_cosets) == 14, "wrong number of Borel cosets")

infinity_frame_cosets = {
    coset for coset in p13_cosets if coset.issubset(borel)
}
require(len(infinity_frame_cosets) == 6,
        "wrong number of frames over infinity")
infinity_start = next(iter(infinity_frame_cosets))
infinity_borel_orbit = {
    frozenset(psl_mul(element, member) for member in infinity_start)
    for element in borel
}
require(infinity_borel_orbit == infinity_frame_cosets,
        "Borel is not transitive on the infinity frame fibre")

fibre_sizes = []
for large_coset in borel_cosets:
    fibre_sizes.append(sum(1 for small_coset in p13_cosets
                           if small_coset.issubset(large_coset)))
require(fibre_sizes == [6] * 14, "frame bundle does not have fibre six")

qr13 = {pow(x, 2, P) for x in range(1, P)}
require(qr13 == {1, 3, 4, 9, 10, 12}, "wrong quadratic-residue fibre")

# The affine part of G/P13 is not merely a 13-by-6 set.  The upper Borel
# acts on root/frame coordinates by the Frobenius law
#
#     (c,kappa).(x,eta)=(kappa*x+c,kappa*eta),
#     (c,kappa)(d,lambda)=(c+kappa*d,kappa*lambda).
#
# Here kappa is a quadratic residue.  This action is regular on the 78
# affine states, while the fibre over infinity is one six-state Borel orbit.
def borel_parameter(matrix):
    a, b, c, d = matrix
    require(c == 0 and a * d % P == 1, "non-Borel matrix in parameter map")
    return (a * b % P, a * a % P)


borel_parameters = {borel_parameter(matrix) for matrix in borel}
affine_states = {(x, eta) for x in range(P) for eta in qr13}
require(borel_parameters == {(c, kappa) for c in range(P) for kappa in qr13},
        "Borel parameters do not form C13 semidirect QR13")
orbit_of_origin_frame = {
    ((kappa * 0 + c) % P, kappa * 1 % P)
    for c, kappa in borel_parameters
}
require(orbit_of_origin_frame == affine_states,
        "Borel action is not regular on affine root/frame states")
for left in borel:
    c, kappa = borel_parameter(left)
    for right in borel:
        d, lam = borel_parameter(right)
        product = psl_mul(left, right)
        require(borel_parameter(product)
                == ((c + kappa * d) % P, kappa * lam % P),
                "Borel semidirect-product law failed")

# The same regular B-set is the arc set of the Paley graph P(13).
def paley_adjacent(x, y):
    return x != y and (y - x) % P in qr13


paley_arcs = {(x, (x + eta) % P) for x in range(P) for eta in qr13}
paley_edges = {frozenset((x, y)) for x, y in paley_arcs}
require(len(paley_arcs) == 78 and len(paley_edges) == 39,
        "wrong Paley edge/arc count")
require(all(sum(paley_adjacent(x, y) for y in range(P)) == 6
            for x in range(P)), "Paley graph is not six-regular")
for x in range(P):
    for y in range(x + 1, P):
        common = sum(paley_adjacent(x, z) and paley_adjacent(y, z)
                     for z in range(P))
        require(common == (2 if paley_adjacent(x, y) else 3),
                "Paley strongly-regular parameters failed")
local_cycle_edges = {
    frozenset(edge) for edge in ((1, 4), (4, 3), (3, 12),
                                 (12, 9), (9, 10), (10, 1))
}
actual_local_edges = {
    frozenset((x, y)) for x in qr13 for y in qr13
    if x < y and paley_adjacent(x, y)
}
require(actual_local_edges == local_cycle_edges,
        "Paley local frame graph is not the frozen C6")

# Full affine slopes are regular on every non-diagonal ordered pair; a
# nonsquare slope swaps the Paley graph with its complement.
all_affine_pair_images = {
    (shift, (shift + slope) % P)
    for shift in range(P) for slope in range(1, P)
}
require(len(all_affine_pair_images) == P * (P - 1),
        "full affine group is not regular on ordered distinct pairs")
require(all(not paley_adjacent(2 * x % P, 2 * y % P)
            for edge in paley_edges for x, y in [tuple(edge)]),
        "nonsquare multiplier did not swap Paley edge classes")

chi = {1: 1, 4: 3, 3: 2, 12: 6, 9: 4, 10: 5}
require(set(chi) == qr13 and set(chi.values()) == set(range(1, 7)),
        "owner-frame identification is not bijective")
for x in qr13:
    for y in qr13:
        require(chi[x * y % P] == chi[x] * chi[y] % 7,
                "owner-frame identification is not multiplicative")

# THM-2605's fixed phase v=k*r+q realizes the normal C13 translation
# layer.  A Borel root scaling r'=s*r+t would need the displayed affine
# action on q.  It is a uniform physical shift of the whole phase graph only
# when s=1; the QR scaling is genuinely extra data.
phase_action_cells = 0
for k in range(1, P):
    for v in range(P):
        for slope in qr13:
            for shift in range(P):
                deltas = set()
                for root in range(P):
                    q = (v - k * root) % P
                    root_next = (slope * root + shift) % P
                    q_next = (v - k * root_next) % P
                    require(q_next == (slope * q + (1 - slope) * v
                                       - k * shift) % P,
                            "Borel fixed-phase action failed")
                    deltas.add((root_next - root) % P)
                    phase_action_cells += 1
                require((len(deltas) == 1) == (slope == 1),
                        "nontranslation Borel slope became a uniform C13 shift")
require(phase_action_cells == 158184, "wrong Borel phase-action census")

# The positive two-time Paley kernel inherited from arbitrary root-word
# accessibility is constant in the six frame directions.  Verify exactly in
# Z[zeta_6], represented by pairs a+b*zeta with zeta^2=zeta-1, that all five
# nontrivial C6 Fourier modes vanish.
def zeta6_mul(left, right):
    a, b = left
    c, d = right
    return (a * c - b * d, a * d + b * c + b * d)


def zeta6_pow(exponent):
    answer = (1, 0)
    base = (0, 1)
    for _ in range(exponent % 6):
        answer = zeta6_mul(answer, base)
    return answer


for colour in range(1, 6):
    total = [0, 0]
    for frame_exponent in range(6):
        value = zeta6_pow(colour * frame_exponent)
        total[0] += value[0]
        total[1] += value[1]
    require(tuple(total) == (0, 0), "flat Paley kernel has charged C6 mass")


# G is perfect: the normal closure of [S,R] is all of G.
commutator = psl_mul(psl_mul(psl_mul(S, R), psl_inv(S)), psl_inv(R))
conjugates = {
    psl_mul(psl_mul(g, commutator), psl_inv(g))
    for g in psl2
}
derived = {I}
frontier = [I]
generators = sorted(conjugates)
while frontier:
    h = frontier.pop()
    for generator in generators:
        product = psl_mul(h, generator)
        if product not in derived:
            derived.add(product)
            frontier.append(product)
require(len(derived) == 1092, "normal closure of [S,R] is not all of G")


# Trace recurrence and the three projective order-seven trace classes.
def f7(x):
    return (x**6 - 5 * x**4 + 6 * x**2 - 1) % P


trace_roots = [x for x in range(P) if f7(x) == 0]
require(trace_roots == [3, 5, 6, 7, 8, 10], "wrong order-seven trace set")
for x in range(P):
    product = 1
    for root in trace_roots:
        product = product * (x - root) % P
    require(product == f7(x), "trace polynomial factorization failed")

trace_square_cycle = []
y = 9
for _ in range(3):
    trace_square_cycle.append(y)
    y = (y - 2) ** 2 % P
require(trace_square_cycle == [9, 10, 12] and y == 9,
        "arithmetic C3 trace-square cycle failed")


# Ordered norm calculation.
expected_forward = {
    3: [6, 8, 9, 10, 11],
    5: [2, 8, 10, 11, 12],
    6: [1, 3, 9, 11, 12],
}
expected_reverse = {
    3: [2, 3, 4, 5, 7],
    5: [1, 2, 3, 5, 11],
    6: [1, 2, 4, 10, 12],
}

forward_sets = {}
reverse_sets = {}
closing_records = []

for t in (3, 5, 6):
    g = g_matrix(t)
    g_inv = raw_inv(g)
    forward = []
    reverse = []
    for exponent in range(P):
        unipotent = U(exponent)
        factors = []
        for k in range(7):
            factor = raw_mul(raw_mul(raw_pow(g, k), unipotent),
                             raw_pow(g_inv, k))
            factors.append(factor)

        norm_forward = I_RAW
        for factor in factors:
            norm_forward = raw_mul(norm_forward, factor)
        norm_reverse = I_RAW
        for factor in reversed(factors):
            norm_reverse = raw_mul(norm_reverse, factor)

        telescope_forward = raw_mul(raw_pow(raw_mul(unipotent, g), 7),
                                    raw_pow(g_inv, 7))
        telescope_reverse = raw_mul(raw_pow(g, 7),
                                    raw_pow(raw_mul(g_inv, unipotent), 7))
        require(norm_forward == telescope_forward, "forward telescope failed")
        require(norm_reverse == telescope_reverse, "reverse telescope failed")
        require((raw_mul(unipotent, g)[0] + raw_mul(unipotent, g)[3]) % P
                == (t - exponent) % P, "forward trace failed")
        require((raw_mul(g_inv, unipotent)[0] + raw_mul(g_inv, unipotent)[3]) % P
                == (t + exponent) % P, "reverse trace failed")

        if is_scalar(norm_forward):
            if exponent:
                forward.append(exponent)
                closing_records.append((t, "+", exponent,
                                        raw_mul(unipotent, g)))
        if is_scalar(norm_reverse):
            if exponent:
                reverse.append(exponent)
                closing_records.append((t, "-", exponent,
                                        raw_mul(g_inv, unipotent)))

    forward_sets[t] = forward
    reverse_sets[t] = reverse
    require(forward == expected_forward[t], "wrong forward norm set")
    require(reverse == expected_reverse[t], "wrong reverse norm set")

require(len(closing_records) == 30, "wrong number of labelled closing pairs")


# Exact t=3 off-diagonal factors, checked as polynomial functions on F_13.
g3 = g_matrix(3)
g3_inv = raw_inv(g3)
for exponent in range(P):
    factors = [raw_mul(raw_mul(raw_pow(g3, k), U(exponent)),
                       raw_pow(g3_inv, k)) for k in range(7)]
    norm = I_RAW
    for factor in factors:
        norm = raw_mul(norm, factor)
    c_value = exponent
    for root in (6, 8, 9, 10, 11):
        c_value = c_value * (exponent - root) % P
    b_value = 10 * exponent
    for root in (4, 6, 8, 9, 10, 11):
        b_value = b_value * (exponent - root) % P
    require(norm[2] == c_value % P, "t=3 c(a) factor failed")
    require(norm[1] == b_value % P, "t=3 b(a) factor failed")


# Sharp minimal oriented chart cover.
chart_sets = {
    (3, "+"): set(forward_sets[3]),
    (3, "-"): set(reverse_sets[3]),
    (5, "+"): set(forward_sets[5]),
    (5, "-"): set(reverse_sets[5]),
    (6, "+"): set(forward_sets[6]),
    (6, "-"): set(reverse_sets[6]),
}
universe = set(range(1, P))
covering_families = {}
keys = tuple(chart_sets)
for size in range(0, 7):
    covering_families[size] = []
    for family in combinations(keys, size):
        covered = set()
        for key in family:
            covered.update(chart_sets[key])
        if covered == universe:
            covering_families[size].append(family)

require(not covering_families[1] and not covering_families[2],
        "oriented norm atlas unexpectedly uses fewer than three charts")
expected_minimal_covers = [
    ((3, "+"), (3, "-"), (6, "+")),
    ((3, "+"), (3, "-"), (6, "-")),
]
require(covering_families[3] == expected_minimal_covers,
        "wrong minimal oriented norm covers")


# Exact projective cycles of the three Hurwitz matrices.
expected_g_cycles = {
    3: ((INF, 0, 9, 2, 1, 7, 3), (4, 12, 10, 11, 8, 5, 6)),
    5: ((INF, 0, 8, 4, 1, 10, 5), (2, 9, 3, 7, 6, 12, 11)),
    6: ((INF, 0, 11, 5, 1, 8, 6), (2, 10, 3, 9, 4, 7, 12)),
}
for t in (3, 5, 6):
    cycles = projective_cycles(g_matrix(t))
    infinity_cycle = next(c for c in cycles if INF in c)
    finite_cycle = next(c for c in cycles if INF not in c)
    infinity_cycle = rotate_to(list(infinity_cycle), INF)
    finite_cycle = rotate_to(list(finite_cycle), min(finite_cycle))
    require((infinity_cycle, finite_cycle) == expected_g_cycles[t],
            "wrong Hurwitz projective cycles")


# Index each P13 coset and record its projective base point.
p13_coset_index = {coset: j for j, coset in enumerate(p13_cosets)}
p13_coset_base = []
for coset in p13_cosets:
    representative = min(coset)
    p13_coset_base.append(projective_action(representative, INF))
require({x: p13_coset_base.count(x) for x in range(P + 1)}
        == {x: 6 for x in range(P + 1)},
        "affine/projective frame fibre count failed")


labelled_lifted_finite_cycles = 0
example_cycle = None
example_qr_factors = None

for t, orientation, exponent, step in closing_records:
    require(not is_scalar(step), "closing step is unexpectedly scalar")
    require(is_scalar(raw_pow(step, 7)), "closing step lacks projective order seven")
    cycles = projective_cycles(step)
    require(sorted(len(cycle) for cycle in cycles) == [7, 7],
            "closing step does not have two seven-cycles")
    finite_cycles = [cycle for cycle in cycles if INF not in cycle]
    require(len(finite_cycles) == 1, "closing step lacks a unique affine cycle")
    finite_cycle = finite_cycles[0]

    # The fixed-point discriminant is nonsquare.
    trace = (step[0] + step[3]) % P
    discriminant = (trace * trace - 4) % P
    require(discriminant not in qr13 and discriminant != 0,
            "closing order-seven step has a projective fixed point")

    # Direct tangent lift.
    qr_factors = []
    product = 1
    x = finite_cycle[0]
    for _ in range(7):
        factor = derivative_factor(step, x)
        require(factor in qr13, "tangent multiplier left the six-frame fibre")
        qr_factors.append(factor)
        product = product * factor % P
        x = projective_action(step, x)
    require(x == finite_cycle[0] and product == 1,
            "affine tangent frame failed to return")
    for eta in qr13:
        transported = eta
        for factor in qr_factors:
            transported = transported * factor % P
        require(transported == eta, "owner frame failed to close")

    # Independent lift by exact G/P13 coset action.
    step_psl = canonical(step)
    permutation = []
    for coset in p13_cosets:
        image = frozenset(psl_mul(step_psl, element) for element in coset)
        require(image in p13_coset_index, "coset action left G/P13")
        permutation.append(p13_coset_index[image])
    seen = set()
    frame_orbits = []
    for start in range(84):
        if start in seen:
            continue
        orbit = []
        point = start
        while point not in orbit:
            orbit.append(point)
            seen.add(point)
            point = permutation[point]
        require(point == start and len(orbit) == 7,
                "frame action does not have a seven-cycle")
        frame_orbits.append(tuple(orbit))
    require(len(frame_orbits) == 12, "wrong number of frame cycles")
    finite_set = set(finite_cycle)
    finite_frame_orbits = [
        orbit for orbit in frame_orbits
        if all(p13_coset_base[index] in finite_set for index in orbit)
    ]
    require(len(finite_frame_orbits) == 6,
            "wrong number of lifted all-affine frame cycles")
    labelled_lifted_finite_cycles += len(finite_frame_orbits)

    if (t, orientation, exponent) == (3, "+", 6):
        cycle_at_zero = rotate_to(list(finite_cycle), 0)
        factors_at_zero = []
        x = 0
        for _ in range(7):
            factors_at_zero.append(derivative_factor(step, x))
            x = projective_action(step, x)
        example_cycle = cycle_at_zero
        example_qr_factors = tuple(factors_at_zero)

require(labelled_lifted_finite_cycles == 180,
        "wrong labelled chart/exponent frame-cycle count")
require(example_cycle == (0, 2, 7, 9, 8, 11, 1),
        "wrong frozen all-affine example cycle")
require(example_qr_factors == (3, 1, 9, 4, 12, 12, 10),
        "wrong frozen tangent multiplier sequence")
example_owner_factors = tuple(chi[x] for x in example_qr_factors)
require(example_owner_factors == (2, 1, 4, 3, 6, 6, 5),
        "wrong frozen owner multiplier sequence")


# Orientation-sensitive matrices.
def ordered_norms(t, exponent):
    g = g_matrix(t)
    g_inv = raw_inv(g)
    factors = [raw_mul(raw_mul(raw_pow(g, k), U(exponent)),
                       raw_pow(g_inv, k)) for k in range(7)]
    forward = I_RAW
    reverse = I_RAW
    for factor in factors:
        forward = raw_mul(forward, factor)
    for factor in reversed(factors):
        reverse = raw_mul(reverse, factor)
    return forward, reverse


forward_6, reverse_6 = ordered_norms(3, 6)
forward_7, reverse_7 = ordered_norms(3, 7)
require(forward_6 == MINUS_I_RAW and reverse_6 == raw((3, 4, 1, 6)),
        "wrong a=6 orientation hostile")
require(forward_7 == raw((6, 9, 12, 3)) and reverse_7 == MINUS_I_RAW,
        "wrong a=7 orientation hostile")


# Same root word, different tangent holonomy.
root_cycle = (0, 2, 7, 9, 8, 11, 1, 0)
translation_edges = [U(root_cycle[j + 1] - root_cycle[j]) for j in range(7)]
twisted_edges = list(translation_edges)
twisted_edges[0] = raw((2, 1, 0, 7))


def follow_with_derivative(edges):
    point = root_cycle[0]
    tangent = 1
    visited = [point]
    for edge in edges:
        tangent = tangent * derivative_factor(edge, point) % P
        point = projective_action(edge, point)
        visited.append(point)
    return tuple(visited), tangent


plain_visited, plain_tangent = follow_with_derivative(translation_edges)
twisted_visited, twisted_tangent = follow_with_derivative(twisted_edges)
require(plain_visited == root_cycle and twisted_visited == root_cycle,
        "tangent hostile changed the designated root word")
require(plain_tangent == 1 and twisted_tangent == 4,
        "tangent hostile did not change frame holonomy")
require(chi[plain_tangent] == 1 and chi[twisted_tangent] == 3,
        "tangent hostile owner colours are wrong")

# The Paley second endpoint follows the tangent only for affine Borel maps.
# A genuinely Mobius closing step separates the two transports.
paley_hostile = raw((7, 6, 12, 3))
hostile_base = projective_action(paley_hostile, 0)
hostile_natural_endpoint = projective_action(paley_hostile, 1)
hostile_tangent_endpoint = (
    hostile_base + derivative_factor(paley_hostile, 0)
) % P
require((hostile_base, hostile_natural_endpoint, hostile_tangent_endpoint)
        == (2, 0, 5), "Paley affine-only endpoint hostile failed")


# Oriented endpoint Gram-frame lift.  The determinant is the orientation
# sidecar that fails in characteristic two.
def det_columns(endpoint, prime):
    lx, ly, rx, ry = endpoint
    return (lx * ry - ly * rx) % prime


def gram_key(endpoint, prime, retain_det=True):
    lx, ly, rx, ry = endpoint
    vx, vy = (lx - rx) % prime, (ly - ry) % prime
    aa = (lx * lx + ly * ly) % prime
    bb = (lx * rx + ly * ry) % prime
    cc = (rx * rx + ry * ry) % prime
    key = (vx, vy, aa, bb, cc)
    if retain_det:
        key += (det_columns(endpoint, prime),)
    return key


def reconstruct_odd_gram(key, prime):
    vx, vy, aa, bb, cc, delta = key
    norm = (aa - 2 * bb + cc) % prime
    mu = (bb - cc) % prime
    if norm:
        norm_inv = pow(norm, -1, prime)
        rx = (vx * mu - vy * delta) * norm_inv % prime
        ry = (vy * mu + vx * delta) * norm_inv % prime
    else:
        require(vx != 0, "nonzero isotropic vector lost its first coordinate")
        vx_inv = pow(vx, -1, prime)
        hx, hy = 0, vx_inv
        kappa = vy * vx_inv % prime
        require(kappa != 0 and mu == delta * kappa % prime,
                "isotropic Gram compatibility failed")
        denominator = 2 * delta * kappa % prime
        t = (cc - delta * delta * (hx * hx + hy * hy)) * pow(
            denominator, -1, prime) % prime
        rx = (delta * hx + t * vx) % prime
        ry = (delta * hy + t * vy) % prime
    return ((rx + vx) % prime, (ry + vy) % prime, rx, ry)


def endpoint_gram_census(prime):
    endpoints = [
        (lx, ly, rx, ry)
        for lx in range(prime) for ly in range(prime)
        for rx in range(prime) for ry in range(prime)
        if det_columns((lx, ly, rx, ry), prime)
    ]
    full = {}
    no_det = {}
    isotropic_vectors = set()
    for endpoint in endpoints:
        key = gram_key(endpoint, prime, True)
        full.setdefault(key, []).append(endpoint)
        no_det.setdefault(key[:-1], []).append(endpoint)
        vx, vy = key[:2]
        if (vx * vx + vy * vy) % prime == 0:
            isotropic_vectors.add((vx, vy))
        require(reconstruct_odd_gram(key, prime) == endpoint,
                "oriented Gram reconstruction failed")
        require((key[2] - 2 * key[3] + key[4]) % prime
                == (vx * vx + vy * vy) % prime,
                "Gram difference identity failed")
        require((key[2] * key[4] - key[3] * key[3]) % prime
                == key[5] * key[5] % prime,
                "Gram determinant identity failed")
    require(all(len(fibre) == 1 for fibre in full.values()),
            "oriented Gram key stopped being injective")
    histogram = {
        size: sum(1 for fibre in no_det.values() if len(fibre) == size)
        for size in {len(fibre) for fibre in no_det.values()}
    }
    expected_iso = 2 * (prime - 1) if pow(prime - 1, (prime - 1) // 2, prime) == 1 else 0
    expected_single = expected_iso * prime * (prime - 1)
    expected_double = ((prime * prime - 1 - expected_iso)
                       * prime * (prime - 1)) // 2
    require(len(isotropic_vectors) == expected_iso,
            "isotropic-vector count changed")
    require(histogram.get(1, 0) == expected_single
            and histogram.get(2, 0) == expected_double,
            "determinant-forgetting fibre census changed")
    return len(endpoints), len(full), expected_iso, expected_single, expected_double


gram_censuses = {prime: endpoint_gram_census(prime)
                 for prime in (3, 5, 7, 13)}
require(gram_censuses[13] == (26208, 26208, 24, 3744, 11232),
        "F13 oriented Gram census changed")

# Normalized SL2 frame and exact 84-frame endpoint invoice.
frame_of = {
    element: index
    for index, coset in enumerate(p13_cosets)
    for element in coset
}
frame_counts = [0] * len(p13_cosets)
for lx in range(P):
    for ly in range(P):
        for rx in range(P):
            for ry in range(P):
                endpoint = (lx, ly, rx, ry)
                delta = det_columns(endpoint, P)
                if not delta:
                    continue
                delta_inv = pow(delta, -1, P)
                frame_matrix = raw(((lx - rx) % P, rx * delta_inv,
                                    (ly - ry) % P, ry * delta_inv))
                require(raw_det(frame_matrix) == 1,
                        "normalized endpoint frame left SL2")
                frame_counts[frame_of[canonical(frame_matrix)]] += 1
require(set(frame_counts) == {312},
        "endpoint pairs stopped distributing uniformly over 84 frames")

# Right-coordinate covariance under modular generators.
def transpose_raw(matrix):
    a, b, c, d = matrix
    return raw((a, c, b, d))


for endpoint in ((1, 0, 0, 1), (1, 1, 2, 3), (4, 7, 9, 2)):
    lx, ly, rx, ry = endpoint
    endpoint_matrix = raw((lx, rx, ly, ry))
    if raw_det(endpoint_matrix) == 0:
        continue
    gram = raw_mul(transpose_raw(endpoint_matrix), endpoint_matrix)
    for change in (S, R, U(1), raw((2, 0, 0, 7))):
        moved = raw_mul(endpoint_matrix, change)
        moved_gram = raw_mul(transpose_raw(moved), moved)
        predicted = raw_mul(raw_mul(transpose_raw(change), gram), change)
        require(moved_gram == predicted,
                "right-coordinate Gram covariance failed")
        require(raw_det(moved) == raw_det(endpoint_matrix) * raw_det(change) % P,
                "right-coordinate determinant covariance failed")

# Characteristic-two orientation-loss hostile.
f2_endpoints = [
    endpoint
    for endpoint in ((lx, ly, rx, ry)
                     for lx in range(2) for ly in range(2)
                     for rx in range(2) for ry in range(2))
    if det_columns(endpoint, 2)
]
f2_full = {}
for endpoint in f2_endpoints:
    f2_full.setdefault(gram_key(endpoint, 2, True), []).append(endpoint)
f2_collisions = [fibre for fibre in f2_full.values() if len(fibre) > 1]
require(len(f2_endpoints) == 6 and len(f2_full) == 5
        and f2_collisions == [[(0, 1, 1, 0), (1, 0, 0, 1)]],
        "characteristic-two identity/swap hostile changed")


def fmt(values):
    return "{" + ",".join(str(x) for x in values) + "}"


def fmt_cycle(cycle):
    return "(" + ",".join("inf" if x == INF else str(x) for x in cycle) + ")"


print("THM-2626 PSL2(F13) FRAME/NORM AUDIT")
print("group SL2=2184 PSL2=1092 P13=13 Borel=78 frames=84 points=14 fibre=6")
print("affine Borel_torsor=78 semidirect=C13xC6 twisted_law=1 infinity_orbit=6")
print("Paley13 srg=(13,6,2,3) edges=39 arcs=78 local_C6=1 mobius_endpoint_hostile=2:0,5")
print("physical THM2605_layer=C13 phase_cells=158184 Borel_scalings_extra=5 flat_Paley_modes=5")
print("modular S_order=2 R_order=3 U=R^-1S hurwitz_t=3,5,6")
print("perfect normal_closure_[S,R]=1092")
print("trace7 roots=" + fmt(trace_roots) + " trace_square_C3=(9,10,12)")
for t in (3, 5, 6):
    print("norm t={} forward={} reverse={}".format(
        t, fmt(forward_sets[t]), fmt(reverse_sets[t])))
print("minimum_cardinality oriented_charts=3 optimal_covering_triples=2")
for t in (3, 5, 6):
    print("projective t={} cycles={} {}".format(
        t, fmt_cycle(expected_g_cycles[t][0]), fmt_cycle(expected_g_cycles[t][1])))
print("closing labelled_pairs=30 frame_orbits_per_pair=12 affine_lifts_per_pair=6")
print("closing labelled_affine_frame_cycles=180")
print("example t=3 a=6 forward cycle={} qr={} owner={}".format(
    fmt_cycle(example_cycle), fmt(example_qr_factors), fmt(example_owner_factors)))
print("orientation a=6 forward=-I reverse=(3,4;1,6) a=7 roles_reversed=1")
print("root_word_hostile same_endpoints=1 tangent_holonomy=1,4 owner_holonomy=1,3")
for prime in (3, 5, 7, 13):
    total, full, iso, singles, doubles = gram_censuses[prime]
    print("Gram p={} total={} full={} isotropic={} no_delta={}x1+{}x2".format(
        prime, total, full, iso, singles, doubles))
print("Gram p=2 total=6 full=5 hostile=identity_swap")
print("Gram F13 frames=84 endpoints_per_frame=312 covariance=1 metric_intertwiner=0")
print("scope physical_carrier=0 relation_current=0 lrc_rows_removed=0")
print("THM-2626 PASS")
