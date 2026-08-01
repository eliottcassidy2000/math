#!/usr/bin/env python3
"""Exact finite referee for THM-3035's S3 orbit diamond."""

from collections import Counter, deque


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


INF = "inf"


def inv(a, p):
    return pow(a % p, -1, p)


def s_map(r, p):
    if r == INF:
        return 0
    if r == 0:
        return INF
    return inv(r, p)


def c_map(r, p):
    if r == INF:
        return 0
    if (r + 1) % p == 0:
        return INF
    return (-inv(r + 1, p)) % p


def alpha_map(r, p):
    if r == INF:
        return 0
    if r == 0:
        return INF
    return (-inv(r, p)) % p


def compose(left, right):
    return tuple(left[i] for i in right)


def anharmonic_group(p):
    points = list(range(p)) + [INF]
    index = {point: i for i, point in enumerate(points)}
    s_perm = tuple(index[s_map(point, p)] for point in points)
    c_perm = tuple(index[c_map(point, p)] for point in points)
    identity = tuple(range(p + 1))
    group = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for generator in (s_perm, c_perm):
            nxt = compose(generator, current)
            if nxt not in group:
                group.add(nxt)
                queue.append(nxt)
    require(len(group) == 6, f"anharmonic group at p={p} is not S3")
    return points, index, tuple(sorted(group)), s_perm, c_perm


def orbit(point, points, index, group):
    i = index[point]
    return frozenset(points[g[i]] for g in group)


def stabilizer_size(point, index, group):
    i = index[point]
    return sum(g[i] == i for g in group)


def orbit_partition(p):
    points, index, group, s_perm, c_perm = anharmonic_group(p)
    unseen = set(points)
    orbits = []
    while unseen:
        point = min(unseen, key=lambda q: p if q == INF else q)
        orb = orbit(point, points, index, group)
        orbits.append((orb, stabilizer_size(point, index, group)))
        unseen -= set(orb)
    orbits.sort(
        key=lambda item: (
            len(item[0]),
            sorted(p if x == INF else x for x in item[0]),
        )
    )
    return points, index, group, s_perm, c_perm, orbits


def expected_named_orbits(p):
    boundary = frozenset((INF, 0, p - 1))
    harmonic = (
        frozenset((1 % p, (-2) % p, (-inv(2, p)) % p))
        if p != 2
        else boundary
    )
    equianharmonic = frozenset(
        r for r in range(p) if (r * r + r + 1) % p == 0
    )
    return boundary, harmonic, equianharmonic


def fixed_character(orb, points, index, permutation):
    return sum(points[permutation[index[x]]] == x for x in orb)


prime_rows = []
for p in (2, 3, 5, 7, 11, 13, 17, 19):
    points, index, group, s_perm, c_perm, orbits = orbit_partition(p)
    boundary, harmonic, equianharmonic = expected_named_orbits(p)
    require(boundary in [item[0] for item in orbits], f"boundary missing at p={p}")
    if p >= 5:
        require(harmonic in [item[0] for item in orbits], f"harmonic missing at p={p}")
        require(
            stabilizer_size(next(iter(harmonic)), index, group) == 2,
            f"harmonic stabilizer failed at p={p}",
        )
        expected_equi = 2 if p % 3 == 1 else 0
        require(
            len(equianharmonic) == expected_equi,
            f"equianharmonic count failed at p={p}",
        )
        if expected_equi:
            require(
                equianharmonic in [item[0] for item in orbits],
                f"equianharmonic orbit missing at p={p}",
            )
            require(
                stabilizer_size(next(iter(equianharmonic)), index, group) == 3,
                f"equianharmonic stabilizer failed at p={p}",
            )
        named = set(boundary) | set(harmonic) | set(equianharmonic)
        generic = [item for item in orbits if set(item[0]).isdisjoint(named)]
        expected_generic = (p - 7) // 6 if p % 3 == 1 else (p - 5) // 6
        require(len(generic) == expected_generic, f"generic count failed at p={p}")
        require(
            all(len(o) == 6 and stab == 1 for o, stab in generic),
            f"generic stabilizer failed at p={p}",
        )
    prime_rows.append((p, tuple(sorted((len(o), stab) for o, stab in orbits))))

# The two bad-prime degenerations.
_, index2, group2, _, _, orbits2 = orbit_partition(2)
require(len(orbits2) == 1 and len(orbits2[0][0]) == 3, "p2 exhaustion")
points3, index3, group3, s3, c3, orbits3 = orbit_partition(3)
require(
    sorted((len(o), stab) for o, stab in orbits3) == [(1, 6), (3, 2)],
    "p3 collapse",
)
require(s_map(1, 3) == 1 and c_map(1, 3) == 1, "p3 S3-fixed singleton")
require(alpha_map(1, 3) == 2, "p3 PSL/PGL hostile")

# Exact p=13 orbit and character atlas.
points13, index13, group13, s13, c13, orbits13 = orbit_partition(13)
boundary13, harmonic13, equi13 = expected_named_orbits(13)
generic13 = next(o for o, stab in orbits13 if len(o) == 6)
require(boundary13 == frozenset((INF, 0, 12)), "p13 boundary")
require(harmonic13 == frozenset((1, 6, 11)), "p13 harmonic")
require(equi13 == frozenset((3, 9)), "p13 equianharmonic")
require(generic13 == frozenset((2, 4, 5, 7, 8, 10)), "p13 generic")

characters = {
    "regular": (
        len(generic13),
        fixed_character(generic13, points13, index13, s13),
        fixed_character(generic13, points13, index13, c13),
    ),
    "natural3": (
        len(harmonic13),
        fixed_character(harmonic13, points13, index13, s13),
        fixed_character(harmonic13, points13, index13, c13),
    ),
    "parity2": (
        len(equi13),
        fixed_character(equi13, points13, index13, s13),
        fixed_character(equi13, points13, index13, c13),
    ),
    "trivial1": (1, 1, 1),
}
require(
    characters
    == {
        "regular": (6, 0, 0),
        "natural3": (3, 1, 0),
        "parity2": (2, 0, 2),
        "trivial1": (1, 1, 1),
    },
    "permutation characters",
)

# C_r={(g,g.r)} retains the regular six-state frame.  Its second projection
# forgets precisely |Stab(r)| copies.
cooccurrence_fibres = {}
for label, point in (("generic", 2), ("harmonic", 1), ("equianharmonic", 3)):
    i = index13[point]
    pairs = [(g, points13[g[i]]) for g in group13]
    require(len({g for g, _ in pairs}) == 6, f"first projection {label}")
    counts = Counter(image for _, image in pairs)
    require(len(set(counts.values())) == 1, f"second fibres {label}")
    cooccurrence_fibres[label] = next(iter(counts.values()))
pairs_p3 = [(g, points3[g[index3[1]]]) for g in group3]
require(Counter(image for _, image in pairs_p3) == Counter({1: 6}), "p3 fibre6")
cooccurrence_fibres["collapsed_p3"] = 6
require(
    cooccurrence_fibres
    == {"generic": 1, "harmonic": 2, "equianharmonic": 3, "collapsed_p3": 6},
    "cooccurrence multiplicities",
)

# det(s)=-1, det(c)=1.  The projective determinant square class is invariant.
psl_intersections = []
for p in (3, 5, 7, 11, 13, 17, 19):
    minus_one_square = pow(p - 1, (p - 1) // 2, p) == 1
    intersection_size = 6 if minus_one_square else 3
    require(intersection_size == (6 if p % 4 == 1 else 3), f"PSL split p={p}")
    psl_intersections.append((p, intersection_size))
require(
    all(alpha_map(r, 2) == s_map(r, 2) for r in (0, 1, INF)),
    "alpha and s must alias only at p2 control",
)

# The six orders of the three nonzero V4 directions are the regular S3 frame.
vectors = ((1, 0), (0, 1), (1, 1))


def add2(a, b):
    return (a[0] ^ b[0], a[1] ^ b[1])


bases = [(u0, v0) for u0 in vectors for v0 in vectors if u0 != v0]
frames = {(u0, add2(u0, v0), v0) for u0, v0 in bases}
require(len(bases) == 6 and len(frames) == 6, "regular V4 frame")
seen = {((1, 0), (0, 1))}
queue = deque(seen)
while queue:
    u0, v0 = queue.popleft()
    for nxt in ((u0, add2(u0, v0)), (add2(u0, v0), v0)):
        if nxt not in seen:
            seen.add(nxt)
            queue.append(nxt)
require(len(seen) == 6, "Farey L/R frame orbit")

# Exact THM-2056 hostile: same mod-two frame and endpoint gate bits, opposite
# child gate.  The determinant-one and acute conditions are checked literally.
def farey_defect(w, d):
    return d[0] ** 2 + d[1] ** 2 - 91 * (w[0] * d[0] + w[1] * d[1])


w = (1, 0)
u0 = (1, 0)
v0 = (88, 1)
v1 = (90, 1)
require(
    (v0[0] % 2, v0[1] % 2) == (v1[0] % 2, v1[1] % 2) == (0, 1),
    "hostile frame mismatch",
)
require(u0[0] * v0[1] - u0[1] * v0[0] == 1, "hostile determinant0")
require(u0[0] * v1[1] - u0[1] * v1[0] == 1, "hostile determinant1")
require(u0[0] * v0[0] + u0[1] * v0[1] > 0, "hostile acute0")
require(u0[0] * v1[0] + u0[1] * v1[1] > 0, "hostile acute1")
require(
    farey_defect(w, u0) < 0
    and farey_defect(w, v0) < 0
    and farey_defect(w, v1) < 0,
    "hostile endpoint bits",
)
child0 = (u0[0] + v0[0], u0[1] + v0[1])
child1 = (u0[0] + v1[0], u0[1] + v1[1])
farey_values = (
    farey_defect(w, u0),
    farey_defect(w, v0),
    farey_defect(w, v1),
    farey_defect(w, child0),
    farey_defect(w, child1),
)
require(farey_values == (-90, -263, -89, -177, 1), "Farey hostile values")

# Stronger hostile: even the endpoint signs and the exact Gram product do not
# determine the child.  What is missing is the full integral defect lift.
ua, va = (1, -8), (1, -7)
ub, vb = (1, 0), (57, -1)
require(
    tuple(x % 2 for x in ua) == tuple(x % 2 for x in ub) == (1, 0)
    and tuple(x % 2 for x in va) == tuple(x % 2 for x in vb) == (1, 1),
    "fixed-Gram hostile frame",
)
det_a = ua[0] * va[1] - ua[1] * va[0]
det_b = ub[0] * vb[1] - ub[1] * vb[0]
dot_a = ua[0] * va[0] + ua[1] * va[1]
dot_b = ub[0] * vb[0] + ub[1] * vb[1]
require((det_a, det_b, dot_a, dot_b) == (1, -1, 57, 57), "fixed-Gram geometry")
fixed_gram_values = (
    farey_defect(w, ua),
    farey_defect(w, va),
    farey_defect(w, ub),
    farey_defect(w, vb),
    farey_defect(w, (ua[0] + va[0], ua[1] + va[1])),
    farey_defect(w, (ub[0] + vb[0], ub[1] + vb[1])),
)
require(
    fixed_gram_values == (-26, -41, -90, -1937, 47, -1913),
    "fixed-Gram hostile values",
)


print("modular/Farey/anharmonic/quartic orbit-diamond referee")
for p, signature in prime_rows:
    print(f"p={p}:orbit_(size,stabilizer)={signature}")
print("p13_boundary=[0,12,inf];harmonic=[1,6,11];generic=[2,4,5,7,8,10];equianharmonic=[3,9]")
print(f"permutation_characters_id_transposition_3cycle={characters}")
print(f"cooccurrence_second_projection_fibres={cooccurrence_fibres}")
print(f"H_p_intersection_PSL_sizes={psl_intersections}")
print("PSL_involution_alpha=-1/r;PGL_reflection_s=1/r;alias_only_at_p2=PASS")
print("Farey_mod2_frames=6=regular_S3")
print(f"Farey_same_frame_hostile_Fu_Fv0_Fv1_Fchild0_Fchild1={farey_values}")
print(f"Farey_fixed_Gram57_hostile_endpoints_and_children={fixed_gram_values}")
print("all_exact_checks=PASS")
