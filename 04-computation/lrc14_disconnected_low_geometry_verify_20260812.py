#!/usr/bin/env python3
"""Exact arithmetic and realizability checks for disconnected-low geometry."""
from fractions import Fraction as F
from hashlib import sha256
from math import comb, gcd, lcm


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


shape_counts = {1: 1, 2: 8, 3: 94, 4: 1295, 5: 19389}
EXPECTED_REALIZABILITY = "9ad5d29562c8aa99fd6b852c8f9786042495cc93cfc5c968df2a42947970e324"
LOW_UP = frozenset((F(4, 3), F(3, 2), F(2), F(5, 2), F(3), F(4), F(5), F(6)))
LOW_ORIENTED = tuple(sorted(LOW_UP | {1 / x for x in LOW_UP}))

# Coefficient of x^6 in product_n (1-x^n)^(-c_n), by bounded multiset DP.
dp = [0] * 7
dp[0] = 1
for size, colours in shape_counts.items():
    nxt = [0] * 7
    for used in range(7):
        for multiplicity in range((6 - used) // size + 1):
            ways = comb(colours + multiplicity - 1, multiplicity)
            nxt[used + size * multiplicity] += dp[used] * ways
    dp = nxt

partition_rows = (
    ((5, 1), 19389),
    ((4, 2), 10360),
    ((4, 1, 1), 1295),
    ((3, 3), 4465),
    ((3, 2, 1), 752),
    ((3, 1, 1, 1), 94),
    ((2, 2, 2), 120),
    ((2, 2, 1, 1), 36),
    ((2, 1, 1, 1, 1), 8),
    ((1, 1, 1, 1, 1, 1), 1),
)
require(dp[6] == 36520, ("profile coefficient", dp))
require(sum(value for _, value in partition_rows) == 36520, partition_rows)


def low_adjacent(x, y):
    a, b = sorted((F(x), F(y)))
    return b / a in LOW_UP


def primitive(row):
    denominator = lcm(*(x.denominator for x in row))
    values = tuple(int(x * denominator) for x in row)
    divisor = gcd(*values)
    return tuple(x // divisor for x in values)


def connected_shapes():
    """Enumerate every normalized connected low shape through size five."""
    layers = {1: ((F(1),),)}
    layer = {(F(1),)}
    for size in range(2, 6):
        children = set()
        for row in layer:
            have = set(row)
            for x in row:
                for ratio in LOW_ORIENTED:
                    y = x * ratio
                    if y >= 1 and y not in have:
                        children.add(tuple(sorted((*row, y))))
        layer = children
        layers[size] = tuple(sorted(layer))
    require({size: len(rows) for size, rows in layers.items()} == shape_counts,
            ("connected layers", {size: len(rows) for size, rows in layers.items()}))
    return layers


def low_components(row):
    unseen = set(range(len(row)))
    answer = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        component = {seed}
        frontier = [seed]
        while frontier:
            i = frontier.pop()
            hit = {j for j in unseen if low_adjacent(row[i], row[j])}
            unseen -= hit
            component |= hit
            frontier.extend(hit)
        answer.append(primitive(tuple(F(row[i]) for i in sorted(component))))
    return tuple(sorted(answer, key=lambda x: (len(x), x)))


def profile_realizability():
    """Realize every component multiset by recursive factor->6 separation.

    If every old coordinate is at most M and the next normalized component is
    multiplied by 6M+1, every cross ratio is strictly greater than six.
    Hence no cross pair is intrinsically low, while scaling preserves all
    internal ratios.  This proves that the generating-function coefficient is
    a count of realizable profiles, not merely a formal upper universe.
    """
    layers = connected_shapes()
    tokens = tuple((size, row) for size in range(1, 6) for row in layers[size])
    profiles = []

    def visit(start, remaining, chosen):
        if remaining == 0:
            profiles.append(tuple(chosen))
            return
        for index in range(start, len(tokens)):
            size, row = tokens[index]
            if size > remaining:
                break
            chosen.append(row)
            visit(index, remaining - size, chosen)
            chosen.pop()

    visit(0, 6, [])
    require(len(profiles) == 36520, ("profile enumeration", len(profiles)))
    semantic = sha256()
    maximum = 0
    for profile in profiles:
        realized = []
        old_maximum = 0
        signature = tuple(sorted((primitive(row) for row in profile), key=lambda x: (len(x), x)))
        for row in profile:
            component = primitive(row)
            scale = 1 if not realized else 6 * old_maximum + 1
            placed = tuple(scale * x for x in component)
            require(not realized or min(placed) > 6 * old_maximum,
                    ("separation", profile, realized, placed))
            realized.extend(placed)
            old_maximum = max(realized)
        realized = tuple(sorted(realized))
        require(len(realized) == len(set(realized)) == 6, ("collision", profile, realized))
        recovered = low_components(realized)
        require(recovered == signature, ("profile mismatch", signature, recovered, realized))
        maximum = max(maximum, realized[-1])
        semantic.update((repr((signature, realized)) + "\n").encode())
    return len(profiles), maximum, semantic.hexdigest()


def multipartite_tree_count(parts):
    vertices = sum(parts)
    return vertices ** (len(parts) - 2) * __import__("math").prod(
        (vertices - part) ** (part - 1) for part in parts
    )


tree_counts = tuple((parts, multipartite_tree_count(parts)) for parts, _ in partition_rows)
expected_tree_counts = (1, 32, 48, 81, 216, 324, 384, 576, 864, 1296)
require(tuple(value for _, value in tree_counts) == expected_tree_counts, tree_counts)

debt_max = F(186636088362, 11773143757375)
homogeneous_gap = F(1, 21) - debt_max
require(homogeneous_gap == F(1121969414539, 35319431272125), homogeneous_gap)

far_ratio_floor = F(1, 49) - F(6 * 168, 49 * (8 * 168 - 14))
require(far_ratio_floor == F(23, 4655), far_ratio_floor)
require(far_ratio_floor - F(1, 294) == F(43, 27930), far_ratio_floor)
far_ratio_tree_gap = F(5, 294) - debt_max
require(far_ratio_tree_gap == F(570672686921, 494472037809750), far_ratio_tree_gap)

realizability = profile_realizability()
require(realizability == (36520, 9331, EXPECTED_REALIZABILITY), realizability)

print("disconnected_profile_count", dp[6])
print("realizable_profile_count", realizability[0])
print("scale_separated_maximum_coordinate", realizability[1])
print("realizability_semantic_sha256", realizability[2])
print("partition_rows", partition_rows)
print("multipartite_tree_counts", tree_counts)
print("homogeneous_gap", homogeneous_gap)
print("q_ge_8p_physical_floor", far_ratio_floor)
print("q_ge_8p_gap_over_1_294", far_ratio_floor - F(1, 294))
print("five_q_ge_8p_edges_minus_debt", far_ratio_tree_gap)
print("all_exact_checks=PASS")
