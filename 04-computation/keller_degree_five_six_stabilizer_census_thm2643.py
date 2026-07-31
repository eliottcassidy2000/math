#!/usr/bin/env python3
"""Exact transitive-group referee for THM-2643.

The generators are the standard 5T/6T permutation representatives.  The
script is dependency-free, deterministic, and optimized-mode safe.
"""

from collections import Counter, deque
from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(p, q):
    """p after q."""
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    out = [0] * len(p)
    for i, image in enumerate(p):
        out[image] = i
    return tuple(out)


def from_cycles(n, *cycles):
    out = list(range(n))
    for cycle in cycles:
        c = tuple(x - 1 for x in cycle)
        for i, x in enumerate(c):
            out[x] = c[(i + 1) % len(c)]
    return tuple(out)


def subgroup_closure(n, generators):
    identity = tuple(range(n))
    steps = tuple(dict.fromkeys(tuple(generators) + tuple(inverse(g) for g in generators)))
    group = {identity}
    queue = deque([identity])
    while queue:
        x = queue.popleft()
        for g in steps:
            y = compose(x, g)
            if y not in group:
                group.add(y)
                queue.append(y)
    return frozenset(group)


def conjugate(g, h):
    return compose(compose(g, h), inverse(g))


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


def normalize_partition(blocks):
    return tuple(sorted(tuple(sorted(block)) for block in blocks))


def two_triples():
    out = []
    for other in combinations(range(1, 6), 2):
        first = frozenset((0,) + other)
        second = frozenset(set(range(6)) - first)
        out.append(normalize_partition((first, second)))
    return tuple(sorted(set(out)))


def perfect_matchings(points):
    points = tuple(sorted(points))
    if not points:
        return ((),)
    first = points[0]
    out = []
    for partner in points[1:]:
        rest = tuple(x for x in points if x not in (first, partner))
        for tail in perfect_matchings(rest):
            out.append(((first, partner),) + tail)
    return tuple(out)


def induced_action(group, partition):
    blocks = tuple(frozenset(block) for block in partition)
    lookup = {block: i for i, block in enumerate(blocks)}
    induced = set()
    for g in group:
        images = []
        for block in blocks:
            image = frozenset(g[x] for x in block)
            if image not in lookup:
                return None
            images.append(lookup[image])
        induced.add(tuple(images))
    return frozenset(induced)


def normal_orbit_partition(n, normal):
    unseen = set(range(n))
    blocks = []
    while unseen:
        x = min(unseen)
        block = frozenset(g[x] for g in normal)
        blocks.append(block)
        unseen -= block
    return normalize_partition(blocks)


def analyse(label, name, n, cycle_generators, expected):
    generators = tuple(from_cycles(n, *spec) for spec in cycle_generators)
    group = subgroup_closure(n, generators)
    require({g[0] for g in group} == set(range(n)), f"{label} is not transitive")
    stabilizer = frozenset(g for g in group if g[0] == 0)
    normal_generators = {conjugate(g, h) for g in group for h in stabilizer}
    normal = subgroup_closure(n, tuple(normal_generators))
    require(normal <= group and stabilizer <= normal, f"{label} bad normal closure")
    triple = (len(group), len(stabilizer), len(normal))
    require(triple == expected, f"{label} size triple changed: {triple}")

    outside = group - normal
    require(all(g[i] != i for g in outside for i in range(n)), f"{label} quotient coset has a fixed sheet")
    outside_cycles = Counter(cycle_type(g) for g in outside)

    partitions = [normalize_partition(tuple((i,) for i in range(n)))]
    if n == 6:
        partitions.extend(two_triples())
        partitions.extend(normalize_partition(x) for x in perfect_matchings(tuple(range(6))))
        require(len(two_triples()) == 10 and len(set(perfect_matchings(tuple(range(6))))) == 15, "wrong block universe")
    invariant = []
    for partition in sorted(set(partitions)):
        induced = induced_action(group, partition)
        if induced is not None:
            invariant.append((len(partition), len(induced)))
    block_census = Counter(invariant)
    has_regular_quotient = any(blocks > 1 and order == blocks for blocks, order in invariant)
    require(has_regular_quotient == (normal != group), f"{label} regular quotient criterion failed")

    normal_partition = normal_orbit_partition(n, normal)
    normal_induced = induced_action(group, normal_partition)
    require(normal_induced is not None, f"{label} normal orbits are not invariant")
    require(len(normal_induced) == len(group) // len(normal), f"{label} quotient action is not faithful")
    if normal != group:
        require(len(normal_induced) == len(normal_partition) > 1, f"{label} normal quotient is not regular")
    return {
        "label": label,
        "name": name,
        "n": n,
        "group": group,
        "H": stabilizer,
        "N": normal,
        "quotient_order": len(group) // len(normal),
        "outside_cycles": outside_cycles,
        "block_census": block_census,
        "normal_blocks": len(normal_partition),
    }


GROUPS5 = (
    ("5T1", "C5", (((1, 2, 3, 4, 5),),), (5, 1, 1)),
    ("5T2", "D10", (((1, 2, 3, 4, 5),), ((1, 4), (2, 3))), (10, 2, 10)),
    ("5T3", "F20", (((1, 2, 3, 4, 5),), ((1, 2, 4, 3),)), (20, 4, 20)),
    ("5T4", "A5", (((1, 2, 3),), ((3, 4, 5),)), (60, 12, 60)),
    ("5T5", "S5", (((1, 2),), ((1, 2, 3, 4, 5),)), (120, 24, 120)),
)

GROUPS6 = (
    ("6T1", "C6", (((1, 2, 3, 4, 5, 6),),), (6, 1, 1)),
    ("6T2", "S3_regular", (((1, 3, 5), (2, 4, 6)), ((1, 4), (2, 3), (5, 6))), (6, 1, 1)),
    ("6T3", "D12_hexagon", (((1, 4), (2, 3), (5, 6)), ((1, 2, 3, 4, 5, 6),)), (12, 2, 6)),
    ("6T4", "A4_C2_cosets", (((1, 4), (2, 5)), ((1, 3, 5), (2, 4, 6))), (12, 2, 4)),
    ("6T5", "S3xC3", (((1, 4), (2, 5), (3, 6)), ((2, 4, 6),)), (18, 3, 9)),
    ("6T6", "A4xC2", (((1, 3, 5), (2, 4, 6)), ((3, 6),)), (24, 4, 8)),
    ("6T7", "S4_edges", (((1, 4), (2, 5)), ((1, 5), (2, 4)), ((1, 3, 5), (2, 4, 6))), (24, 4, 24)),
    ("6T8", "S4_C4_cosets", (((1, 5), (2, 4), (3, 6)), ((1, 4), (2, 5)), ((1, 3, 5), (2, 4, 6))), (24, 4, 24)),
    ("6T9", "S3xS3", (((1, 4), (2, 5), (3, 6)), ((2, 4, 6),), ((1, 5), (2, 4))), (36, 6, 18)),
    ("6T10", "C3sq_C4", (((1, 4, 5, 2), (3, 6)), ((2, 4, 6),), ((1, 5), (2, 4))), (36, 6, 18)),
    ("6T11", "S4xC2", (((1, 5), (2, 4)), ((1, 3, 5), (2, 4, 6)), ((3, 6),)), (48, 8, 48)),
    ("6T12", "A5_degree6", (((1, 4), (5, 6)), ((1, 2, 3, 4, 6),)), (60, 10, 60)),
    ("6T13", "S3wrC2", (((2, 4),), ((1, 4), (2, 5), (3, 6)), ((2, 4, 6),)), (72, 12, 36)),
    ("6T14", "S5_exotic", (((1, 2), (3, 4), (5, 6)), ((1, 2, 3, 4, 6),)), (120, 20, 120)),
    ("6T15", "A6", (((1, 2, 3),), ((1, 2), (3, 4, 5, 6))), (360, 60, 360)),
    ("6T16", "S6", (((1, 2),), ((1, 2, 3, 4, 5, 6),)), (720, 120, 720)),
)

rows5 = [analyse(label, name, 5, generators, expected) for label, name, generators, expected in GROUPS5]
rows6 = [analyse(label, name, 6, generators, expected) for label, name, generators, expected in GROUPS6]
require([row["label"] for row in rows5 if row["N"] != row["group"]] == ["5T1"], "wrong quintic exclusions")
require(
    [row["label"] for row in rows6 if row["N"] != row["group"]]
    == ["6T1", "6T2", "6T3", "6T4", "6T5", "6T6", "6T9", "6T10", "6T13"],
    "wrong sextic exclusions",
)

# Exact action-sensitive quartic controls.
quartic_a4 = analyse("4T4", "A4_natural", 4, (((1, 2, 3),), ((2, 3, 4),)), (12, 3, 12))
quartic_s4 = analyse("4T5", "S4_natural", 4, (((1, 2),), ((1, 2, 3, 4),)), (24, 6, 24))
require(quartic_a4["quotient_order"] == quartic_s4["quotient_order"] == 1, "natural quartic controls failed")

print("THM-2643 KELLER STABILIZER CENSUS EXACT REFEREE")
for row in rows5:
    print(
        f"{row['label']} {row['name']} |G,H,N,Q|=({len(row['group'])},{len(row['H'])},{len(row['N'])},{row['quotient_order']}) "
        f"gate={'EXCLUDE' if row['N'] != row['group'] else 'PASS'} outside={dict(sorted(row['outside_cycles'].items()))}"
    )
for row in rows6:
    print(
        f"{row['label']} {row['name']} |G,H,N,Q|=({len(row['group'])},{len(row['H'])},{len(row['N'])},{row['quotient_order']}) "
        f"gate={'EXCLUDE' if row['N'] != row['group'] else 'PASS'} blocks={dict(sorted(row['block_census'].items()))} "
        f"outside={dict(sorted(row['outside_cycles'].items()))}"
    )
print("degree5_excluded=1/5 degree6_excluded=9/16 all_failure_cosets=derangements")
print("regular_quotient_criterion=PASS singleton_quotients_included=YES")
print("quartic_controls=A4_natural_PASS,S4_natural_PASS action_sensitive=YES")
print("ALL CHECKS PASSED")
