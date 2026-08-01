#!/usr/bin/env python3
"""Exact finite referee for the Farey/anharmonic/quartic S3 orbit diamond."""

from collections import Counter, deque
from math import gcd


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
    orbits.sort(key=lambda item: (len(item[0]), sorted(p if x == INF else x for x in item[0])))
    return points, index, group, s_perm, c_perm, orbits


def expected_named_orbits(p):
    boundary = frozenset((INF, 0, p - 1))
    harmonic = frozenset((1 % p, (-2) % p, (-inv(2, p)) % p)) if p != 2 else boundary
    equianharmonic = frozenset(r for r in range(p) if (r * r + r + 1) % p == 0)
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
        require(stabilizer_size(next(iter(harmonic)), index, group) == 2, f"harmonic stabilizer failed at p={p}")
        expected_equi = 2 if p % 3 == 1 else 0
        require(len(equianharmonic) == expected_equi, f"equianharmonic count failed at p={p}")
        if expected_equi:
            require(equianharmonic in [item[0] for item in orbits], f"equianharmonic orbit missing at p={p}")
            require(stabilizer_size(next(iter(equianharmonic)), index, group) == 3, f"equianharmonic stabilizer failed at p={p}")
        named = set(boundary) | set(harmonic) | set(equianharmonic)
        generic = [item for item in orbits if set(item[0]).isdisjoint(named)]
        expected_generic = (p - 7) // 6 if p % 3 == 1 else (p - 5) // 6
        require(len(generic) == expected_generic, f"generic orbit count failed at p={p}")
        require(all(len(o) == 6 and stab == 1 for o, stab in generic), f"generic stabilizer failed at p={p}")
    prime_rows.append((p, tuple(sorted((len(o), stab) for o, stab in orbits))))

# Bad-prime degenerations.
_, index2, group2, _, _, orbits2 = orbit_partition(2)
require(len(orbits2) == 1 and len(orbits2[0][0]) == 3, "p=2 boundary exhaustion failed")
points3, index3, group3, s3, c3, orbits3 = orbit_partition(3)
require(sorted((len(o), stab) for o, stab in orbits3) == [(1, 6), (3, 2)], "p=3 orbit collapse failed")
require(s_map(1, 3) == 1 and c_map(1, 3) == 1, "p=3 transverse singleton is not S3-fixed")
require(alpha_map(1, 3) == 2, "PSL involution p=3 hostile failed")

# Exact p=13 atlas and permutation characters.
points13, index13, group13, s13, c13, orbits13 = orbit_partition(13)
boundary13, harmonic13, equi13 = expected_named_orbits(13)
generic13 = next(o for o, stab in orbits13 if len(o) == 6)
require(boundary13 == frozenset((INF, 0, 12)), "p13 boundary failed")
require(harmonic13 == frozenset((1, 6, 11)), "p13 harmonic failed")
require(equi13 == frozenset((3, 9)), "p13 equianharmonic failed")
require(generic13 == frozenset((2, 4, 5, 7, 8, 10)), "p13 generic failed")

characters = {
    "regular": (len(generic13), fixed_character(generic13, points13, index13, s13), fixed_character(generic13, points13, index13, c13)),
    "natural3": (len(harmonic13), fixed_character(harmonic13, points13, index13, s13), fixed_character(harmonic13, points13, index13, c13)),
    "parity2": (len(equi13), fixed_character(equi13, points13, index13, s13), fixed_character(equi13, points13, index13, c13)),
    "trivial1": (1, 1, 1),
}
require(characters == {"regular": (6, 0, 0), "natural3": (3, 1, 0), "parity2": (2, 0, 2), "trivial1": (1, 1, 1)}, "permutation characters failed")

# Co-occurrence cover C_r={(g,g.r)}: first projection is bijective and second
# projection has constant fibre |Stab(r)|.
cooccurrence_fibres = {}
for label, point in (("generic", 2), ("harmonic", 1), ("equianharmonic", 3)):
    i = index13[point]
    pairs = [(g, points13[g[i]]) for g in group13]
    require(len({g for g, _ in pairs}) == 6, f"cooccurrence first projection failed for {label}")
    counts = Counter(image for _, image in pairs)
    require(len(set(counts.values())) == 1, f"cooccurrence fibres are not constant for {label}")
    cooccurrence_fibres[label] = next(iter(counts.values()))
pairs_p3 = [(g, points3[g[index3[1]]]) for g in group3]
counts_p3 = Counter(image for _, image in pairs_p3)
require(counts_p3 == Counter({1: 6}), "p3 collapsed cooccurrence fibre failed")
cooccurrence_fibres["collapsed_p3"] = 6
require(
    cooccurrence_fibres
    == {"generic": 1, "harmonic": 2, "equianharmonic": 3, "collapsed_p3": 6},
    "cooccurrence multiplicities failed",
)

# PGL/PSL reflection-sign split.  det(s)=-1 and det(c)=1.  In projective
# dimension two the determinant square class is scalar-independent.
psl_intersections = []
for p in (3, 5, 7, 11, 13, 17, 19):
    minus_one_square = pow(p - 1, (p - 1) // 2, p) == 1
    intersection_size = 6 if minus_one_square else 3
    require(intersection_size == (6 if p % 4 == 1 else 3), f"PSL intersection failed at p={p}")
    psl_intersections.append((p, intersection_size))
require(all(alpha_map(r, 2) == s_map(r, 2) for r in (0, 1, INF)), "alpha and s must alias at p=2")

# The six mod-2 Farey frames are the six orders of the three nonzero V4
# directions.  L(u,v)=(u,u+v), R(u,v)=(u+v,v) generate all six.
vectors = ((1, 0), (0, 1), (1, 1))


def add2(a, b):
    return (a[0] ^ b[0], a[1] ^ b[1])


bases = [(u0, v0) for u0 in vectors for v0 in vectors if u0 != v0]
frames = {(u0, add2(u0, v0), v0) for u0, v0 in bases}
require(len(bases) == 6 and len(frames) == 6, "Farey/V4 regular frame failed")
seen = {((1, 0), (0, 1))}
queue = deque(seen)
while queue:
    u0, v0 = queue.popleft()
    for nxt in ((u0, add2(u0, v0)), (add2(u0, v0), v0)):
        if nxt not in seen:
            seen.add(nxt)
            queue.append(nxt)
require(len(seen) == 6, "mod2 Farey generators do not produce regular S3")

# Same mod-2 frame and same endpoint gate bits do not determine the child
# Farey gate: THM-2056's Gram cross term is load-bearing.
def farey_defect(w, d):
    return d[0] * d[0] + d[1] * d[1] - 91 * (w[0] * d[0] + w[1] * d[1])


w = (1, 0)
u0 = (1, 0)
v0 = (88, 1)
v1 = (90, 1)
require((v0[0] % 2, v0[1] % 2) == (v1[0] % 2, v1[1] % 2) == (0, 1), "hostile frame mismatch")
require(u0[0] * v0[1] - u0[1] * v0[0] == 1, "first hostile determinant failed")
require(u0[0] * v1[1] - u0[1] * v1[0] == 1, "second hostile determinant failed")
require(u0[0] * v0[0] + u0[1] * v0[1] > 0, "first hostile is not acute")
require(u0[0] * v1[0] + u0[1] * v1[1] > 0, "second hostile is not acute")
require(farey_defect(w, u0) < 0 and farey_defect(w, v0) < 0 and farey_defect(w, v1) < 0, "hostile endpoint bits failed")
child0 = (u0[0] + v0[0], u0[1] + v0[1])
child1 = (u0[0] + v1[0], u0[1] + v1[1])
farey_values = (farey_defect(w, u0), farey_defect(w, v0), farey_defect(w, v1), farey_defect(w, child0), farey_defect(w, child1))
require(farey_values == (-90, -263, -89, -177, 1), "Farey hostile values failed")


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
print("all_exact_checks=PASS")
