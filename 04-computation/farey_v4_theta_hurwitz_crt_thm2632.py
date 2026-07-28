#!/usr/bin/env python3
"""Exact referee for THM-2632.

The script is dependency-free, deterministic, and keeps every truth-bearing
check active under ``python -O``.
"""

from collections import Counter, deque
from fractions import Fraction
from itertools import permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mmul(a, b, modulus):
    return tuple(
        tuple(
            sum(a[i][k] * b[k][j] for k in range(len(b))) % modulus
            for j in range(len(b[0]))
        )
        for i in range(len(a))
    )


def meye(n):
    return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))


def mpow(a, exponent, modulus):
    require(exponent >= 0, "matrix exponent must be nonnegative")
    out = meye(len(a))
    base = a
    e = exponent
    while e:
        if e & 1:
            out = mmul(out, base, modulus)
        base = mmul(base, base, modulus)
        e //= 2
    return out


def minv2(a, modulus):
    det = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) % modulus
    require(det != 0, "singular 2 by 2 matrix")
    inv_det = pow(det, -1, modulus)
    return (
        ((a[1][1] * inv_det) % modulus, (-a[0][1] * inv_det) % modulus),
        ((-a[1][0] * inv_det) % modulus, (a[0][0] * inv_det) % modulus),
    )


def mreduce(a, modulus):
    return tuple(tuple(x % modulus for x in row) for row in a)


def mneg(a, modulus):
    return tuple(tuple((-x) % modulus for x in row) for row in a)


def pcanon(a, modulus):
    a = mreduce(a, modulus)
    negative = mneg(a, modulus)
    return min(a, negative)


def pmul(a, b, modulus):
    return pcanon(mmul(a, b, modulus), modulus)


def pinv_matrix(a, modulus):
    return pcanon(minv2(a, modulus), modulus)


def projective_order(a, modulus, limit=100):
    identity = pcanon(meye(2), modulus)
    x = identity
    for order in range(1, limit + 1):
        x = pmul(x, a, modulus)
        if x == identity:
            return order
    raise RuntimeError("projective order exceeded audit limit")


def sl2(modulus):
    return sorted(
        (
            ((a, b), (c, d))
            for a, b, c, d in product(range(modulus), repeat=4)
            if (a * d - b * c) % modulus == 1
        )
    )


def mvec(a, v, modulus):
    return tuple(sum(a[i][j] * v[j] for j in range(len(v))) % modulus for i in range(len(a)))


def vadd(a, b, modulus=2):
    return tuple((x + y) % modulus for x, y in zip(a, b))


def det2(a, modulus):
    return (a[0][0] * a[1][1] - a[0][1] * a[1][0]) % modulus


def gl2(modulus):
    return sorted(
        (
            ((a, b), (c, d))
            for a, b, c, d in product(range(modulus), repeat=4)
            if det2(((a, b), (c, d)), modulus) != 0
        )
    )


V4 = ((0, 0), (1, 0), (0, 1), (1, 1))
# This reference order makes pi(I)=(0,1,2).  It is a gauge; the theta cuts are not.
DIRS = ((1, 0), (1, 1), (0, 1))
DIR_INDEX = {v: i for i, v in enumerate(DIRS)}
PAIR_COORDS = ((0, 1), (0, 2), (1, 2))
L2 = ((1, 1), (0, 1))
R2 = ((1, 0), (1, 1))


def flank_order(u):
    left = (u[0][0], u[1][0])
    right = (u[0][1], u[1][1])
    middle = vadd(left, right)
    return tuple(DIR_INDEX[x] for x in (left, middle, right))


def inversion_bits(order):
    pos = {value: i for i, value in enumerate(order)}
    return tuple(int(pos[a] > pos[b]) for a, b in PAIR_COORDS)


def bfs_distance(adjacency, source):
    distance = {source: 0}
    queue = deque([source])
    while queue:
        x = queue.popleft()
        for y in adjacency[x]:
            if y not in distance:
                distance[y] = distance[x] + 1
                queue.append(y)
    return distance


states = gl2(2)
require(len(states) == 6, "GL2(F2) must have six states")
state_index = {u: i for i, u in enumerate(states)}
orders = [flank_order(u) for u in states]
require(set(orders) == set(permutations(range(3))), "flank orders do not exhaust S3")
bits = [inversion_bits(order) for order in orders]
require(len(set(bits)) == 6, "inversion embedding is not injective")
require(set(product(range(2), repeat=3)) - set(bits) == {(0, 1, 0), (1, 0, 1)}, "wrong omitted cube vertices")

edges = {}
adjacency = {i: set() for i in range(6)}
for i, u in enumerate(states):
    for name, generator in (("L", L2), ("R", R2)):
        child = mmul(u, generator, 2)
        j = state_index[child]
        key = tuple(sorted((i, j)))
        order = orders[i]
        theta = order[0] if name == "L" else order[2]
        flipped_pair = tuple(sorted((order[1], order[2]))) if name == "L" else tuple(sorted((order[0], order[1])))
        require(vadd(DIRS[flipped_pair[0]], DIRS[flipped_pair[1]]) == DIRS[theta], "theta is not the omitted direction")
        coordinate = PAIR_COORDS.index(flipped_pair)
        changed = [k for k in range(3) if bits[i][k] != bits[j][k]]
        require(changed == [coordinate], "Farey edge does not flip its unique inversion coordinate")
        if key in edges:
            require(edges[key] == (name, theta, coordinate), "edge label is orientation-dependent")
        else:
            edges[key] = (name, theta, coordinate)
        adjacency[i].add(j)

require(len(edges) == 6, "Farey residue graph must have six edges")
require(all(len(adjacency[i]) == 2 for i in adjacency), "Farey residue graph must be 2-regular")
require(len(bfs_distance(adjacency, 0)) == 6, "Farey residue graph must be connected")
for i in range(6):
    graph_distance = bfs_distance(adjacency, i)
    for j in range(6):
        hamming = sum(x != y for x, y in zip(bits[i], bits[j]))
        require(graph_distance[j] == hamming, "inversion embedding is not isometric")

edge_label_census = Counter((name, theta) for name, theta, _ in edges.values())
require(set(edge_label_census.values()) == {1} and len(edge_label_census) == 6, "edges are not C2 times theta")
theta_census = Counter(coordinate for _, _, coordinate in edges.values())
require(theta_census == Counter({0: 2, 1: 2, 2: 2}), "wrong theta-class sizes")

# The two omitted orientations are exactly the cyclic tournaments.
def tournament_outdegrees(bit_word):
    out = [0, 0, 0]
    for bit, (a, b) in zip(bit_word, PAIR_COORDS):
        winner, loser = (b, a) if bit else (a, b)
        out[winner] += 1
    return tuple(sorted(out))


require(all(tournament_outdegrees(x) == (0, 1, 2) for x in bits), "a residue state is not transitive")
require(all(tournament_outdegrees(x) == (1, 1, 1) for x in ((0, 1, 0), (1, 0, 1))), "omitted vertices are not cyclic")

# Left GL2(F2) action preserves the L/R edge colour and permutes theta labels.
left_orbits = []
unseen = set(edges)
while unseen:
    seed = min(unseen)
    orbit = set()
    for h in states:
        image = tuple(sorted((state_index[mmul(h, states[seed[0]], 2)], state_index[mmul(h, states[seed[1]], 2)])))
        orbit.add(image)
    left_orbits.append(orbit)
    unseen -= orbit
require(sorted(len(orbit) for orbit in left_orbits) == [3, 3], "left S3 edge orbits must have sizes 3,3")
require(all(len({edges[e][0] for e in orbit}) == 1 for orbit in left_orbits), "left action mixes L/R colours")

# THM-2056 defect parity.
parity_rows = []
for w in V4:
    lam = vadd((1, 1), w)
    if lam == (0, 0):
        for order in orders:
            values = tuple((lam[0] * DIRS[d][0] + lam[1] * DIRS[d][1]) % 2 for d in order)
            require(values == (0, 0, 0), "zero channel failed")
        parity_rows.append((w, "none", (0, 0, 0)))
        continue
    kernel = [d for d in DIRS if (lam[0] * d[0] + lam[1] * d[1]) % 2 == 0]
    require(len(kernel) == 1, "nonzero covector has wrong kernel")
    for order in orders:
        values = tuple((lam[0] * DIRS[d][0] + lam[1] * DIRS[d][1]) % 2 for d in order)
        require(sorted(values) == [0, 1, 1], "Farey triangle parity is not 0,1,1")
    parity_rows.append((w, DIR_INDEX[kernel[0]], (0, 1, 1)))

# The ternary parameter branches reduce to swap, swap, stay; the triple action is trivial.
ternary_2 = (
    ((0, 1), (-1, 2)),
    ((0, 1), (1, 2)),
    ((1, 0), (2, 1)),
)
swap2 = ((0, 1), (1, 0))
identity2 = meye(2)
ternary_shadows = tuple(tuple(tuple(x % 2 for x in row) for row in a) for a in ternary_2)
require(ternary_shadows == (swap2, swap2, identity2), "wrong ternary mod-two shadows")
for i, u in enumerate(states):
    antipode = state_index[mmul(u, swap2, 2)]
    require(bfs_distance(adjacency, i)[antipode] == 3, "ternary swap is not the hexagon antipode")
    require(sum(x != y for x, y in zip(bits[i], bits[antipode])) == 3, "ternary swap does not complement comparison bits")
    require(state_index[mmul(u, identity2, 2)] == i, "ternary stay branch moved the state")
berggren = (
    ((1, -2, 2), (2, -1, 2), (2, -2, 3)),
    ((1, 2, 2), (2, 1, 2), (2, 2, 3)),
    ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3)),
)
require(all(tuple(tuple(x % 2 for x in row) for row in a) == meye(3) for a in berggren), "Berggren shadows are not trivial")

# Ordinary graceful boundary: C6 fails, every one-edge deletion is P6 and succeeds.
cycle_edges = tuple(sorted(edges))


def graceful_witness(vertices, graph_edges):
    q = len(graph_edges)
    for labels in permutations(range(q + 1), len(vertices)):
        if len(set(labels)) != len(labels):
            continue
        differences = {abs(labels[u] - labels[v]) for u, v in graph_edges}
        if differences == set(range(1, q + 1)):
            return labels
    return None


require(graceful_witness(tuple(range(6)), cycle_edges) is None, "C6 unexpectedly graceful")
path_witnesses = {}
for deleted in cycle_edges:
    remaining = tuple(edge for edge in cycle_edges if edge != deleted)
    witness = graceful_witness(tuple(range(6)), remaining)
    require(witness is not None, "one-edge deletion is not graceful")
    path_witnesses[deleted] = witness

# Hurwitz norm parity and the two mod-26 lifts.
P13 = 13


def u_matrix(a, modulus):
    return ((1, a % modulus), (0, 1))


def g_matrix(t, modulus):
    return ((0, 1), ((-1) % modulus, t % modulus))


def ordered_norm(t, orientation, a, modulus):
    g = g_matrix(t, modulus)
    g_inv = minv2(g, modulus)
    u = u_matrix(a, modulus)
    factors = []
    for k in range(7):
        factors.append(mmul(mmul(mpow(g, k, modulus), u, modulus), mpow(g_inv, k, modulus), modulus))
    if orientation == "R":
        factors.reverse()
    out = meye(2)
    for factor in factors:
        out = mmul(out, factor, modulus)
    if orientation == "F":
        telescope = mmul(mpow(mmul(u, g, modulus), 7, modulus), mpow(g_inv, 7, modulus), modulus)
    else:
        telescope = mmul(mpow(g, 7, modulus), mpow(mmul(g_inv, u, modulus), 7, modulus), modulus)
    require(out == telescope, "ordered norm does not telescope")
    return out


def projectively_identity(a, modulus):
    return a[0][1] == 0 and a[1][0] == 0 and a[0][0] == a[1][1] and a[0][0] != 0


success = []
success_table = {}
for t in (3, 5, 6):
    for orientation in ("F", "R"):
        row = []
        for a in range(1, 13):
            norm13 = ordered_norm(t, orientation, a, 13)
            norm2 = ordered_norm(t, orientation, a, 2)
            require(norm2 == u_matrix(a, 2), "seven-edge norm lost its parity shadow")
            require(ordered_norm(t, orientation, a + 13, 13) == norm13, "mod-13 lift changed the norm")
            require(ordered_norm(t, orientation, a + 13, 2) == u_matrix(a + 13, 2), "second CRT lift has wrong shadow")
            require(u_matrix(a, 2) != u_matrix(a + 13, 2), "the two CRT lifts have equal parity")
            require(mpow(u_matrix(a, 2), 7, 2) == norm2, "restored holonomy has different parity shadow")
            if projectively_identity(norm13, 13):
                row.append(a)
                success.append((t, orientation, a))
        success_table[(t, orientation)] = tuple(row)

expected_success = {
    (3, "F"): (6, 8, 9, 10, 11),
    (3, "R"): (2, 3, 4, 5, 7),
    (5, "F"): (2, 8, 10, 11, 12),
    (5, "R"): (1, 2, 3, 5, 11),
    (6, "F"): (1, 3, 9, 11, 12),
    (6, "R"): (1, 2, 4, 10, 12),
}
require(success_table == expected_success, "Hurwitz success atlas changed")
require(len(success) == 30, "wrong number of successful norms")


def fixed_vectors(a):
    return sum(mvec(a, x, 2) == x for x in V4)


fixed_census = Counter()
for _, _, a in success:
    fixed_census[(a % 2, fixed_vectors(u_matrix(a, 2)))] += 1
require(set(fixed_census) == {(0, 4), (1, 2)}, "wrong pure-linear fixed-section boundary")

# The full simultaneous 2/3 shadow is PSL2(Z/6)=S3 x A4, not merely parity.
classes6 = sorted({pcanon(a, 6) for a in sl2(6)})
classes3 = sorted({pcanon(a, 3) for a in sl2(3)})
require(len(sl2(6)) == 144 and len(classes6) == 72, "wrong PSL2(Z/6) order")
require(len(classes3) == 12, "wrong PSL2(F3) order")
crt_images = {
    (pcanon(mreduce(a, 2), 2), pcanon(mreduce(a, 3), 3))
    for a in classes6
}
require(len(crt_images) == 72, "CRT map to S3 x A4 is not bijective")
require(all(projective_order(a, 6) in (1, 2, 3, 6) for a in classes6), "mod-six quotient has exponent above six")

identity6 = pcanon(meye(2), 6)


def pcommutator(a, b, modulus):
    return pmul(pmul(pmul(a, b, modulus), pinv_matrix(a, modulus), modulus), pinv_matrix(b, modulus), modulus)


commutator_generators = {pcommutator(a, b, 6) for a in classes6 for b in classes6}
derived = {identity6}
queue = deque([identity6])
while queue:
    x = queue.popleft()
    for generator in commutator_generators:
        y = pmul(x, generator, 6)
        if y not in derived:
            derived.add(y)
            queue.append(y)
require(len(derived) == 12, "wrong derived subgroup order")
require(all(pmul(a, b, 6) == pmul(b, a, 6) for a in derived for b in derived), "derived subgroup is not C3 x V4")
derived_orders = Counter(projective_order(a, 6) for a in derived)
require(derived_orders == Counter({6: 6, 2: 3, 3: 2, 1: 1}), "wrong C3 x V4 order census")

u6 = pcanon(u_matrix(1, 6), 6)
u6_powers = tuple(pcanon(mpow(u6, k, 6), 6) for k in range(6))
require(len(set(u6_powers)) == 6 and projective_order(u6, 6) == 6, "standard parabolic is not C6")
cosets = []
for power in u6_powers:
    coset = frozenset(pmul(power, n, 6) for n in derived)
    if coset not in cosets:
        cosets.append(coset)
require(len(cosets) == 6 and len(set().union(*cosets)) == 72, "U does not generate the C6 abelianization")
require(set(u6_powers).intersection(derived) == {identity6}, "C6 complement meets the derived group")

# Identify the C3 and V4 factors and the simultaneous order-two/order-three action.
c3_factor = {
    a for a in derived
    if pcanon(mreduce(a, 3), 3) == pcanon(meye(2), 3)
}
v4_factor = {
    a for a in derived
    if pcanon(mreduce(a, 2), 2) == pcanon(meye(2), 2)
}
require(len(c3_factor) == 3 and len(v4_factor) == 4, "derived factors have wrong sizes")
conjugated_c3 = {pmul(pmul(u6, a, 6), pinv_matrix(u6, 6), 6) for a in c3_factor}
conjugated_v4 = {pmul(pmul(u6, a, 6), pinv_matrix(u6, 6), 6) for a in v4_factor}
require(conjugated_c3 == c3_factor and conjugated_v4 == v4_factor, "U does not normalize the two derived factors")
nontrivial_c3 = c3_factor - {identity6}
require(all(pmul(pmul(u6, a, 6), pinv_matrix(u6, 6), 6) == pinv_matrix(a, 6) for a in nontrivial_c3), "U does not invert C3")
v4_nontrivial = v4_factor - {identity6}
v4_orbit = []
x = min(v4_nontrivial)
for _ in range(3):
    v4_orbit.append(x)
    x = pmul(pmul(u6, x, 6), pinv_matrix(u6, 6), 6)
require(set(v4_orbit) == v4_nontrivial and x == v4_orbit[0], "U does not cycle V4 nonzero elements")

# The norm law is universal in every exponent-six group.
norm6_checks = 0
for g in classes6:
    g_inverse = pinv_matrix(g, 6)
    for a in range(6):
        ua = pcanon(u_matrix(a, 6), 6)
        forward = pcanon(mmul(mpow(mmul(ua, g, 6), 7, 6), mpow(g_inverse, 7, 6), 6), 6)
        reverse = pcanon(mmul(mpow(g, 7, 6), mpow(mmul(g_inverse, ua, 6), 7, 6), 6), 6)
        require(forward == ua and reverse == ua, "exponent-six norm did not return its root step")
        norm6_checks += 2
require(norm6_checks == 72 * 6 * 2, "wrong mod-six norm audit size")

for t, orientation, a in success:
    shadows = []
    norm13 = ordered_norm(t, orientation, a, 13)
    for k in range(6):
        lift = a + 13 * k
        require(ordered_norm(t, orientation, lift, 13) == norm13, "CRT lift changed normalized mod-13 closure")
        shadow = pcanon(ordered_norm(t, orientation, lift, 6), 6)
        expected = pcanon(u_matrix(lift, 6), 6)
        require(shadow == expected, "normalized CRT lift has wrong C6 shadow")
        shadows.append(shadow)
    require(set(shadows) == set(u6_powers), "six exponent lifts do not exhaust the C6 sidecar")


def p1_action(a, x, modulus=3):
    if x is None:
        numerator, denominator = a[0][0] % modulus, a[1][0] % modulus
    else:
        numerator = (a[0][0] * x + a[0][1]) % modulus
        denominator = (a[1][0] * x + a[1][1]) % modulus
    if denominator == 0:
        return None
    return numerator * pow(denominator, -1, modulus) % modulus


p1f3 = (None, 0, 1, 2)
fixed_profile6 = {}
for a in range(6):
    ua2 = u_matrix(a, 2)
    ua3 = u_matrix(a, 3)
    fixed_profile6[a] = (
        fixed_vectors(ua2),
        sum(p1_action(ua3, x) == x for x in p1f3),
    )
require(fixed_profile6 == {0: (4, 4), 1: (2, 1), 2: (4, 1), 3: (2, 4), 4: (4, 1), 5: (2, 1)}, "wrong C2/C3 fixed-count profile")

# The original Hurwitz pair has a different integral CRT shadow; its mod-13
# normalization matrix is singular mod 2 and mod 6, so the shadow cannot move.
original_a = ((-1, 2), (3, -7))
original_c = ((5, -17), (-17, 58))
require(det2(original_a, 13) == 1 and det2(original_c, 13) == 1, "original pair is not determinant one")
for modulus in (2, 6, 13):
    a0 = mreduce(original_a, modulus)
    c0 = mreduce(original_c, modulus)
    descending = mmul(mpow(a0, 7, modulus), mpow(mmul(minv2(a0, modulus), c0, modulus), 7, modulus), modulus)
    if modulus in (2, 6):
        require(pcanon(descending, modulus) == pcanon(c0, modulus), "original exponent-six norm forgot its own root step")
    else:
        require(projectively_identity(descending, 13), "original Hurwitz norm did not close mod 13")
require(projective_order(mreduce(original_c, 2), 2) == 3, "original mod-two root shadow is not order three")
require(projective_order(u_matrix(1, 2), 2) == 2, "normalized root shadow is not order two")
normalizing_q = ((4, 5), (0, 10))
require(det2(normalizing_q, 13) == 1, "normalization is not invertible mod 13")
require(det2(normalizing_q, 2) == 0 and det2(normalizing_q, 6) not in (1, 5), "normalization unexpectedly transports the CRT shadow")

# AGL(2,2)=S4 and the seventh-power affine detector.
def pcompose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def ppow(p, exponent):
    out = tuple(range(len(p)))
    base = p
    e = exponent
    while e:
        if e & 1:
            out = pcompose(out, base)
        base = pcompose(base, base)
        e //= 2
    return out


def pinverse(p):
    out = [0] * len(p)
    for i, image in enumerate(p):
        out[image] = i
    return tuple(out)


def cycle_type(p):
    seen = set()
    lengths = []
    for start in range(len(p)):
        if start in seen:
            continue
        x = start
        length = 0
        while x not in seen:
            seen.add(x)
            length += 1
            x = p[x]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def affine_perm(a, c):
    return tuple(V4.index(vadd(mvec(a, x, 2), c)) for x in V4)


def translation_perm(c):
    return tuple(V4.index(vadd(x, c)) for x in V4)


affine_maps = [(a, c, affine_perm(a, c)) for a in states for c in V4]
require(len({p for _, _, p in affine_maps}) == 24, "AGL2(F2) is not faithful of order 24")
cycle_census = Counter(cycle_type(p) for _, _, p in affine_maps)
require(cycle_census == Counter({(3, 1): 8, (2, 1, 1): 6, (4,): 6, (2, 2): 3, (1, 1, 1, 1): 1}), "wrong S4 cycle census")
delta_census = Counter()
for a, c, h in affine_maps:
    delta = pcompose(ppow(h, 7), pinverse(h))
    require(delta == ppow(h, 6), "delta7 is not h^6")
    require(ppow(h, 13) == h, "h^13 is not h")
    ctype = cycle_type(h)
    delta_census[(ctype, delta != tuple(range(4)))] += 1
    if ctype == (4,):
        fixed_directions = [v for v in DIRS if mvec(a, v, 2) == v]
        require(len(fixed_directions) == 1, "four-cycle linear shadow has wrong fixed channel")
        omitted = fixed_directions[0]
        require(delta == translation_perm(omitted), "delta7 misses the omitted V4 direction")
    else:
        require(delta == tuple(range(4)), "delta7 fires off the four-cycle class")
require(sum(count for (ctype, nonzero), count in delta_census.items() if nonzero) == 6, "delta7 must fire six times")
for a in states:
    if a == meye(2) or mpow(a, 2, 2) != meye(2):
        continue
    fibre_types = Counter(cycle_type(affine_perm(a, c)) for c in V4)
    require(fibre_types == Counter({(2, 1, 1): 2, (4,): 2}), "same linear reflection has wrong affine split")
identity_fibre = [(c, affine_perm(meye(2), c)) for c in V4]
require(cycle_type(identity_fibre[0][1]) == (1, 1, 1, 1), "identity affine map changed class")
require(all(ppow(h, 6) == tuple(range(4)) for _, h in identity_fibre), "delta7 sees a pure V4 translation")

# Exact closed/a.e. two-comb radius and the Farey-mediant hostile.
def closed_radius(a, b):
    intervals = []
    for speed in (a, b):
        for k in range(0, speed + 2):
            intervals.append((Fraction(14 * k - 1, 14 * speed), Fraction(14 * k + 1, 14 * speed)))
    intervals.sort()
    right = Fraction(0)
    started = False
    for left, endpoint in intervals:
        if left <= 0 <= endpoint:
            right = max(right, endpoint)
            started = True
        elif started and left <= right:
            right = max(right, endpoint)
        elif started and left > right:
            break
    require(started, "central component was not found")
    return right


r_1_13 = closed_radius(1, 13)
r_1_14 = closed_radius(1, 14)
r_2_27 = closed_radius(2, 27)
r_open_1_13 = Fraction(1, 14)
require(r_1_13 == Fraction(15, 182), "wrong (1,13) radius")
require(r_1_14 == Fraction(15, 196), "wrong (1,14) radius")
require(r_2_27 == Fraction(5, 126), "wrong mediant radius")
require(abs(Fraction(13, 14) - 1) == Fraction(1, 14), "open-radius seam hostile failed")
require(1 * 14 - 1 * 13 == 1 and (1 + 1, 13 + 14) == (2, 27), "not a Farey flank and mediant")
require(r_2_27 < min(r_1_13, r_1_14), "radius did not violate Farey quasi-concavity")


print("THM-2632 FAREY/V4/CRT EXACT REFEREE")
print(f"GL2_F2_states={len(states)} residue_edges={len(edges)} theta_sizes={dict(sorted(theta_census.items()))}")
for i, (u, order, bit_word) in enumerate(zip(states, orders, bits)):
    print(f"state {i}: U={u} order={order} inversions={bit_word}")
print(f"omitted_cube_vertices={sorted(set(product(range(2), repeat=3)) - set(bits))}")
print(f"left_S3_edge_orbits={sorted(len(x) for x in left_orbits)} labels={sorted(sorted({edges[e][0] for e in x}) for x in left_orbits)}")
print("edge_addresses=" + repr(sorted((name, theta, coordinate) for name, theta, coordinate in edges.values())))
print("owner_parity=" + repr(parity_rows))
print(f"ternary_parameter_shadows={ternary_shadows} berggren_shadows=identity")
print("graceful_C6=0 graceful_one_edge_deletions=6")
for deleted in cycle_edges:
    print(f"path delete={edges[deleted][:2]} edge={deleted} witness={path_witnesses[deleted]}")
for key in sorted(success_table):
    print(f"norm_success {key}={success_table[key]}")
print(f"norm_success_total={len(success)} pure_linear_fixed_census={dict(sorted(fixed_census.items()))}")
print(f"Q6_PSL2_Z6={len(classes6)} derived={len(derived)} derived_orders={dict(sorted(derived_orders.items()))} abelianization=C6")
print(f"Q6_norm_law_checks={norm6_checks} C2xC3_fixed_profile={fixed_profile6}")
print("CRT_lifts=a+13k(k=0..5) same_mod13=YES exhaust_C6=YES")
print("original_Hurwitz_shadow=root_C_order3 normalized_shadow=U_order2 normalization_CRT_singular=YES")
print(f"AGL2_F2_cycle_census={dict(sorted(cycle_census.items()))}")
print(f"delta7_nonzero=6 delta7_class=(4,) h13_equals_h=YES")
print(f"closed_radii (1,13)={r_1_13} (1,14)={r_1_14} mediant(2,27)={r_2_27}")
print(f"literal_open_radius (1,13)={r_open_1_13}")
print("Farey_radius_quasiconcavity=FALSE")
print("ALL CHECKS PASSED")
