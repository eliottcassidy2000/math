#!/usr/bin/env python3
"""Exact order-certificate audit for THM-3244's reset Rips boundary."""

import ast
import base64
import hashlib
import struct
import zlib
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / "04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py":
        "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa",
    ROOT / "05-knowledge/results/gmc_complete_physical_bank_unique_reset_thm3238.out":
        "77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b",
}


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


for dependency, expected in DEPENDENCIES.items():
    require(lf_sha256(dependency) == expected, ("dependency drift", dependency.name))

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assert_nodes == 0, "optimization-sensitive assert")
require(float_literals == 0, "floating literal")


# This compressed rank vector is aligned with the canonical nonempty
# submultiset order below.  Rank zero is the largest exact THM-3238
# coordinate.  It was independently extracted from the full 4,319-coordinate
# integer vector; the dependency script/output and raw certificate are pinned.
RANK_CERTIFICATE = (
    "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"
)
RANK_RAW_SHA256 = "ac7aeaf94958c04034137dbbffaf18f494485446421ad9006b6657a15b3487d7"

# For each of the 22 lawful response rows in the immutable THM-3238
# certificate, this row-major bit bank records the states admitting a strict
# one-pole, Q-monotone ascent in that row.  As with the H-rank certificate
# above, provenance is established by comparison with the complete exact
# response bank, not by this certificate's self-digest alone.
ROW_COVER_CERTIFICATE = (
    "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"
    "sA=="
)
ROW_COVER_RAW_SHA256 = "cf8a66dce650ea90fbdc6333fc571506c9bb9d626a9dc2f7ccdaedd2e7ab6cdb"
ROW_COVER_SIZES = (
    4127, 4215, 4094, 4113, 4208, 4154, 4137, 4039, 4211, 4014, 4211,
    4157, 4039, 4218, 3963, 4207, 4146, 4224, 4154, 4038, 4143, 4227,
)
ROW_COVERING_PAIRS = (
    (2, 7), (2, 10), (2, 13), (2, 17), (2, 19), (2, 21),
    (3, 9), (3, 11), (3, 16), (3, 22), (7, 14), (7, 18),
    (7, 22), (10, 11), (10, 16), (10, 22), (11, 13), (11, 17),
    (11, 21), (12, 13), (12, 19), (13, 14), (13, 18), (13, 22),
    (14, 19), (16, 17), (16, 21), (17, 22), (18, 19), (19, 22),
    (21, 22),
)

# Coordinate bits 0,...,7 record which Q-directed one-pole moves raise row 2
# or row 10.  The first 4,319 bytes are row 2 and the second 4,319 are row 10.
ROW_TWO_TEN_DIRECTION_CERTIFICATE = (
    "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"
    "hkqkmHzqY/F8m4KBAfc5D/8BT3KHKw=="
)
ROW_TWO_TEN_DIRECTION_RAW_SHA256 = (
    "2253c0f3b59b57966f9b2357f46c7d396fb8803b9f63bf285cbc9f105ef5388c"
)

# Exact (Delta row 2, Delta row 10) values on the Q-directed edges of two
# sharp states.  Together they cover every positive constant blend ratio.
BLEND_TRAP_A = (2, 2, 2, 3, 3, 4, 5, 5, 6, 7, 8)
BLEND_TRAP_B = (1, 3, 3, 4, 5)
BLEND_TRAP_A_DELTAS = (
    (-685954250014060056, -13986661762200),       # insert 1
    (647427527551915200, -82016379613632),         # delete 2
    (-26061685060581752400, -36145322647680),      # delete 5
)
BLEND_TRAP_B_DELTAS = (
    (-42501120955424475929920368, 6415016001245710560),  # insert 6
    (-61290274303868901672848544, 7303226904731459520),  # insert 7
    (-30235956848857099232025408, 3561194191531466880),  # insert 8
)

VALUES = tuple(range(1, 9))
CAPACITIES = (4, 3, 2, 2, 2, 1, 1, 1)
RESET = (1, 3, 3, 4, 5, 6, 7, 8)


def state_from_counts(counts):
    return tuple(
        value
        for value, count in zip(VALUES, counts)
        for _ in range(count)
    )


states = tuple(sorted(
    (
        state_from_counts(counts)
        for counts in product(*(range(capacity + 1) for capacity in CAPACITIES))
        if any(counts)
    ),
    key=lambda state: (len(state), state),
))
require(len(states) == 4319, "physical state census")
require(
    hashlib.sha256(repr(states).encode("ascii")).hexdigest()
    == "c060e22f900232b608ad3bce1f6b24cae51b6eb45c138b679f4c698fbad2c6a2",
    "canonical state-bank digest",
)
state_index = {state: index for index, state in enumerate(states)}
count_vectors = tuple(
    tuple(Counter(state)[value] for value in VALUES)
    for state in states
)
reset_counts = count_vectors[state_index[RESET]]

mask_bytes = (len(states) + 7) // 8
row_cover_raw = zlib.decompress(base64.b64decode("".join(ROW_COVER_CERTIFICATE)))
require(hashlib.sha256(row_cover_raw).hexdigest() == ROW_COVER_RAW_SHA256,
        "row-cover certificate digest")
require(len(row_cover_raw) == 22 * mask_bytes, "row-cover byte length")
row_cover_sets = tuple(
    frozenset(
        index for index in range(len(states))
        if row_cover_raw[(row - 1) * mask_bytes + index // 8]
        & (1 << (index % 8))
    )
    for row in range(1, 23)
)
require(tuple(map(len, row_cover_sets)) == ROW_COVER_SIZES,
        "lawful-row cover sizes")
nonreset_indices = frozenset(range(len(states))) - {state_index[RESET]}
require(all(cover <= nonreset_indices for cover in row_cover_sets),
        "row cover contains reset")
covering_pairs = tuple(
    (left + 1, right + 1)
    for left in range(22) for right in range(left + 1, 22)
    if row_cover_sets[left] | row_cover_sets[right] == nonreset_indices
)
require(covering_pairs == ROW_COVERING_PAIRS, "covering-pair census")
require(max(map(len, row_cover_sets)) == 4227 < len(nonreset_indices),
        "one-row selector boundary")
row_two, row_ten = row_cover_sets[1], row_cover_sets[9]
require(
    (len(row_two), len(row_ten), len(row_two & row_ten),
     len(row_two - row_ten), len(row_ten - row_two),
     len(nonreset_indices - (row_two | row_ten)))
    == (4215, 4014, 3911, 304, 103, 0),
    "distinguished two-row atlas",
)

direction_raw = zlib.decompress(base64.b64decode(
    "".join(ROW_TWO_TEN_DIRECTION_CERTIFICATE)
))
require(hashlib.sha256(direction_raw).hexdigest()
        == ROW_TWO_TEN_DIRECTION_RAW_SHA256, "direction certificate digest")
require(len(direction_raw) == 2 * len(states), "direction byte length")
row_direction_masks = (
    tuple(direction_raw[:len(states)]),
    tuple(direction_raw[len(states):]),
)
require(
    tuple(frozenset(index for index, mask in enumerate(masks) if mask)
          for masks in row_direction_masks)
    == (row_two, row_ten),
    "direction supports disagree with row covers",
)
common_direction_histogram = tuple(sorted(Counter(
    (row_direction_masks[0][index] & row_direction_masks[1][index]).bit_count()
    for index in row_two & row_ten
).items()))
require(common_direction_histogram
        == ((0, 3453), (1, 294), (2, 113), (3, 45), (4, 6)),
        "overlap common-direction histogram")


def positive_blend_trap_interval(deltas):
    lower = Fraction(0)
    upper = None
    for delta_two, delta_ten in deltas:
        if delta_two > 0:
            bound = Fraction(-delta_ten, delta_two)
            upper = bound if upper is None else min(upper, bound)
        elif delta_two < 0:
            lower = max(lower, Fraction(-delta_ten, delta_two))
        else:
            require(delta_ten <= 0, "empty constant-blend trap interval")
    require(upper is None or lower <= upper, "inconsistent blend interval")
    return lower, upper


blend_interval_a = positive_blend_trap_interval(BLEND_TRAP_A_DELTAS)
blend_interval_b = positive_blend_trap_interval(BLEND_TRAP_B_DELTAS)
require(blend_interval_a == (
    Fraction(0), Fraction(427168643821, 3372018372666225),
), "small-ratio blend trap")
require(blend_interval_b == (
    Fraction(44548722230872990, 295146673301558860624447), None,
), "large-ratio blend trap")
require(blend_interval_b[0] < blend_interval_a[1],
        "two-state blend intervals lost overlap")
require(state_index[BLEND_TRAP_A] in row_two - row_ten
        and state_index[BLEND_TRAP_B] in row_ten - row_two,
        "blend traps lost exclusive-chart typing")

rank_raw = zlib.decompress(base64.b64decode("".join(RANK_CERTIFICATE)))
require(hashlib.sha256(rank_raw).hexdigest() == RANK_RAW_SHA256, "rank digest")
require(len(rank_raw) == 2 * len(states), "rank byte length")
ranks = struct.unpack("<%dH" % len(states), rank_raw)
require(tuple(sorted(ranks)) == tuple(range(len(states))), "rank permutation")
require(states[ranks.index(0)] == RESET, "unique exact maximum")


def edit_distance(left, right):
    return sum(abs(a - b) for a, b in zip(left, right))


def reset_monotone(source, target):
    return all(
        abs(target_i - reset_i) <= abs(source_i - reset_i)
        for source_i, target_i, reset_i in zip(source, target, reset_counts)
    )


def physical_neighbors(counts):
    answer = []
    for coordinate, capacity in enumerate(CAPACITIES):
        if counts[coordinate] < capacity:
            changed = list(counts)
            changed[coordinate] += 1
            answer.append(tuple(changed))
        if counts[coordinate] and sum(counts) > 1:
            changed = list(counts)
            changed[coordinate] -= 1
            answer.append(tuple(changed))
    return tuple(answer)


def q_neighbor_index(index, coordinate):
    counts = list(count_vectors[index])
    if counts[coordinate] < reset_counts[coordinate]:
        counts[coordinate] += 1
    elif counts[coordinate] > reset_counts[coordinate]:
        counts[coordinate] -= 1
    else:
        return None
    if not any(counts):
        return None
    return state_index[state_from_counts(tuple(counts))]


# The direction certificate turns the two-chart cover into a finite directed
# atlas.  Dynamic programming along the strictly decreasing Q-distance grades
# computes the least possible number of chart switches on a complete route.
q_distances = tuple(edit_distance(counts, reset_counts) for counts in count_vectors)
infinity = len(states)
switch_cost_by_first_chart = [[infinity, infinity] for _ in states]
for index in sorted(range(len(states)), key=lambda item: q_distances[item]):
    if index == state_index[RESET]:
        continue
    for chart in range(2):
        candidates = []
        mask = row_direction_masks[chart][index]
        for coordinate in range(8):
            if not mask & (1 << coordinate):
                continue
            target_index = q_neighbor_index(index, coordinate)
            require(target_index is not None, "illegal certified direction")
            require(q_distances[target_index] + 1 == q_distances[index],
                    "certified direction not Q-monotone")
            if target_index == state_index[RESET]:
                candidates.append(0)
            else:
                candidates.append(min(
                    switch_cost_by_first_chart[target_index][next_chart]
                    + int(next_chart != chart)
                    for next_chart in range(2)
                ))
        if candidates:
            switch_cost_by_first_chart[index][chart] = min(candidates)

minimum_switches = tuple(
    min(costs) for index, costs in enumerate(switch_cost_by_first_chart)
    if index != state_index[RESET]
)
require(all(cost < infinity for cost in minimum_switches),
        "two-chart atlas lost a routed state")
switch_histogram = tuple(sorted(Counter(minimum_switches).items()))
require(switch_histogram == ((0, 716), (1, 3600), (2, 2)),
        "two-chart switch-depth histogram")
fixed_row_routes = tuple(
    sum(switch_cost_by_first_chart[index][chart] == 0
        for index in nonreset_indices)
    for chart in range(2)
)
require(fixed_row_routes == (534, 182), "fixed-row route census")
require(sum(min(costs) == 0 for index, costs in enumerate(
    switch_cost_by_first_chart) if index != state_index[RESET]) == 716,
        "zero-switch union census")
sharp_switch_states = tuple(sorted(
    states[index] for index in nonreset_indices
    if min(switch_cost_by_first_chart[index]) == 2
))
require(sharp_switch_states == (
    (1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8),
    (2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8),
), "sharp two-switch states")


def induced_region_stats(vertices):
    vertices = frozenset(vertices)
    edge_twice = 0
    unseen = set(vertices)
    component_sizes = []
    while unseen:
        start = unseen.pop()
        stack = [start]
        size = 0
        while stack:
            vertex = stack.pop()
            size += 1
            neighbors = tuple(
                state_index[state_from_counts(neighbor)]
                for neighbor in physical_neighbors(count_vectors[vertex])
            )
            edge_twice += sum(neighbor in vertices for neighbor in neighbors)
            new = set(neighbors) & unseen
            unseen.difference_update(new)
            stack.extend(new)
        component_sizes.append(size)
    edges = edge_twice // 2
    components = tuple(sorted(component_sizes, reverse=True))
    return (len(vertices), edges, components,
            edges - len(vertices) + len(components))


row_two_only_stats = induced_region_stats(row_two - row_ten)
row_ten_only_stats = induced_region_stats(row_ten - row_two)
row_overlap_stats = induced_region_stats(row_two & row_ten)
require(row_two_only_stats == (
    304, 603,
    (155, 56, 11, 7, 7, 7, 6, 5, 5, 5, 4, 4, 4, 4, 4, 3,
     2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    328,
), "row-two-exclusive transition geometry")
require(row_ten_only_stats == (
    103, 179, (37, 35, 12, 7, 7, 5), 82,
), "row-ten-exclusive transition geometry")
require(row_overlap_stats == (3911, 18554, (3911,), 14644),
        "row-overlap transition geometry")


# Strict one-pole local maxima of the exact supporting functional.
local_traps = []
for index, counts in enumerate(count_vectors):
    if states[index] == RESET:
        continue
    if all(
        ranks[state_index[state_from_counts(neighbor)]] > ranks[index]
        for neighbor in physical_neighbors(counts)
    ):
        local_traps.append(states[index])
require(len(local_traps) == 32, "one-pole local-trap census")


# Minimal reset-monotone edit radius to any strictly higher exact coordinate.
order = tuple(sorted(range(len(states)), key=lambda index: ranks[index]))
escape_radii = [None] * len(states)
for index, source in enumerate(count_vectors):
    if states[index] == RESET:
        continue
    best = None
    for target_index in order[:ranks[index]]:
        target = count_vectors[target_index]
        if reset_monotone(source, target):
            distance = edit_distance(source, target)
            best = distance if best is None else min(best, distance)
    require(best is not None, ("missing monotone escape", states[index]))
    escape_radii[index] = best

escape_histogram = tuple(sorted(Counter(
    radius for radius in escape_radii if radius is not None
).items()))
require(
    escape_histogram
    == ((1, 4231), (2, 45), (3, 23), (4, 10), (5, 1),
        (6, 2), (7, 2), (8, 2), (9, 1), (10, 1)),
    "monotone escape histogram",
)
sink_spectrum = tuple(
    1 + sum(radius > scale for radius in escape_radii if radius is not None)
    for scale in range(1, 11)
)
require(sink_spectrum == (88, 43, 20, 10, 9, 7, 5, 3, 2, 1),
        "Rips sink spectrum")
sharp_traps = tuple(
    states[index] for index, radius in enumerate(escape_radii) if radius == 10
)
require(
    sharp_traps == ((1, 1, 1, 1, 2, 2, 2, 3, 3, 4),),
    "sharp radius-ten trap",
)


# A minimal two-hub atlas gives a stronger two-jump routing.
hub_low = (3,)
hub_high = (1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8)
hub_indices = (state_index[hub_low], state_index[hub_high])
far_indices = tuple(
    index for index, counts in enumerate(count_vectors)
    if edit_distance(counts, reset_counts) > 10
)
require(len(far_indices) == 133, "far-state census")


def hub_covers(hub_index, source_index):
    source = count_vectors[source_index]
    hub = count_vectors[hub_index]
    return (
        ranks[hub_index] < ranks[source_index]
        and reset_monotone(source, hub)
        and edit_distance(source, hub) <= 10
    )


cover_sets = tuple(
    frozenset(index for index in far_indices if hub_covers(hub, index))
    for hub in hub_indices
)
require(
    (len(cover_sets[0]), len(cover_sets[1]), len(cover_sets[0] & cover_sets[1]))
    == (36, 127, 30),
    "two-hub cover counts",
)
require(cover_sets[0] | cover_sets[1] == frozenset(far_indices),
        "two-hub atlas lost a far state")
require(
    tuple(edit_distance(count_vectors[hub], reset_counts) for hub in hub_indices)
    == (7, 3),
    "hub-to-reset radii",
)

candidate_hubs = tuple(
    index for index, counts in enumerate(count_vectors)
    if edit_distance(counts, reset_counts) <= 10
)
single_cover_sizes = tuple(
    sum(hub_covers(hub, index) for index in far_indices)
    for hub in candidate_hubs
)
require(max(single_cover_sizes) == 127 < len(far_indices),
        "one-hub minimality boundary")
maximal_single_hubs = tuple(
    states[hub] for hub, size in zip(candidate_hubs, single_cover_sizes)
    if size == 127
)
require(
    maximal_single_hubs
    == (
        (1, 1, 2, 3, 4, 5, 6, 7, 8),
        (1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8),
    ),
    "maximal single-hub boundary",
)

# Near states jump directly to Q.  Far states use one of the two hubs and
# then Q.  Every step raises H, is reset-monotone, and has radius at most ten.
for index, counts in enumerate(count_vectors):
    if states[index] == RESET:
        continue
    if edit_distance(counts, reset_counts) <= 10:
        require(ranks[state_index[RESET]] < ranks[index], "direct reset ascent")
    else:
        require(any(hub_covers(hub, index) for hub in hub_indices),
                ("missing first hub", states[index]))
require(all(ranks[state_index[RESET]] < ranks[hub] for hub in hub_indices),
        "hub-to-reset ascent")


# Exact H_0 merge tree of the one-pole superlevel filtration.  Vertices enter
# in decreasing H order (increasing certificate rank); elder-rule union-find
# records every finite component bar.
adjacency = tuple(
    tuple(state_index[state_from_counts(neighbor)] for neighbor in physical_neighbors(counts))
    for counts in count_vectors
)
parent = list(range(len(states)))
active = [False] * len(states)
birth_rank = [None] * len(states)
birth_vertex = [None] * len(states)
merge_bars = []
component_profile = []
component_count = 0


def find_component(index):
    while parent[index] != index:
        parent[index] = parent[parent[index]]
        index = parent[index]
    return index


for threshold_rank, vertex in enumerate(order):
    active[vertex] = True
    parent[vertex] = vertex
    higher_roots = sorted({
        find_component(neighbor) for neighbor in adjacency[vertex] if active[neighbor]
    })
    if not higher_roots:
        birth_rank[vertex] = threshold_rank
        birth_vertex[vertex] = vertex
        component_count += 1
    else:
        survivor = min(higher_roots, key=lambda root: birth_rank[root])
        parent[vertex] = survivor
        for root in higher_roots:
            root = find_component(root)
            survivor = find_component(survivor)
            if root == survivor:
                continue
            if birth_rank[root] < birth_rank[survivor]:
                root, survivor = survivor, root
            merge_bars.append((
                threshold_rank - birth_rank[root],
                states[birth_vertex[root]],
                birth_rank[root],
                states[vertex],
                threshold_rank,
            ))
            parent[root] = survivor
            component_count -= 1
    component_profile.append(component_count)

merge_bars.sort(reverse=True)
require(component_count == 1, "one-pole graph disconnected")
require(len(merge_bars) == 32, "merge-tree finite-bar census")
require(max(component_profile) == 19, "merge-tree maximum component count")
require(
    hashlib.sha256(repr(tuple(merge_bars)).encode("ascii")).hexdigest()
    == "6f6aeb54bb02aa89cc063396f25c7883dc254c7c97339f1459d1c6125ba65ce2",
    "merge-tree digest",
)
require(
    hashlib.sha256(repr(tuple(component_profile)).encode("ascii")).hexdigest()
    == "07012e89bdf8799b5b9928d1a25200e5924256f19e740635fd95ab5203cd88bb",
    "component-profile digest",
)
bar_by_birth = {bar[1]: bar for bar in merge_bars}
require(
    bar_by_birth[(1,)]
    == (2096, (1,), 1, (1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 8), 2097),
    "largest persistent trap",
)
require(
    bar_by_birth[hub_low] == (581, (3,), 16, (1, 3), 597),
    "low-hub persistence bar",
)
require(
    bar_by_birth[sharp_traps[0]]
    == (484, sharp_traps[0], 12, (1, 1, 1, 2, 2, 2, 3, 3, 4), 496),
    "radius-ten trap persistence bar",
)
require(hub_high not in bar_by_birth and ranks[state_index[hub_high]] == 36,
        "high hub should lie in an older superlevel component")


print("THM-3244 exact reset Rips/non-Morse order certificate")
print("dependency_hash_checks=2")
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("state_bank=4319,rank_certificate_sha256=%s" % RANK_RAW_SHA256)
print("one_pole_nonreset_local_maxima=%d" % len(local_traps))
print("lawful_row_cover_sizes=%s" % repr(ROW_COVER_SIZES))
print("lawful_row2_row10_split=(4215,4014,3911,304,103,0)")
print("lawful_two_row_covering_pairs=%d" % len(covering_pairs))
print("best_single_lawful_row_cover=4227_of_4318")
print("row_cover_certificate_sha256=%s" % ROW_COVER_RAW_SHA256)
print("row2_row10_one_pole_selector_routes_in_exact_Q_distance=PASS")
print("row2_row10_overlap_common_direction_histogram=%s" % repr(
    common_direction_histogram
))
print("row2_row10_switch_histogram=%s" % repr(switch_histogram))
print("row2_row10_fixed_routes=%s,sharp_states=%s" % (
    repr(fixed_row_routes), repr(sharp_switch_states),
))
print("row2_only_graph=(vertices=304,components=29,edges=603,cycle_rank=328)")
print("row10_only_graph=(vertices=103,components=6,edges=179,cycle_rank=82)")
print("overlap_graph=(vertices=3911,components=1,edges=18554,cycle_rank=14644)")
print("direction_certificate_sha256=%s" % ROW_TWO_TEN_DIRECTION_RAW_SHA256)
print("constant_blend_trap_A=(0,%s],state=%s" % (
    str(blend_interval_a[1]), repr(BLEND_TRAP_A),
))
print("constant_blend_trap_B=[%s,infinity),state=%s" % (
    str(blend_interval_b[0]), repr(BLEND_TRAP_B),
))
print("no_constant_positive_row2_row10_blend=PASS,two_states_sharp=PASS")
print("reset_monotone_escape_histogram=%s" % repr(escape_histogram))
print("reset_monotone_sink_spectrum_radius1_to10=%s" % repr(sink_spectrum))
print("sharp_radius10_trap=%s" % repr(sharp_traps[0]))
print("far_states_radius_gt10=%d" % len(far_indices))
print("two_hubs=%s" % repr((hub_low, hub_high)))
print("two_hub_cover_low_high_overlap=%s" % repr((
    len(cover_sets[0]), len(cover_sets[1]), len(cover_sets[0] & cover_sets[1])
)))
print("max_single_hub_cover=127_of_133,two_hubs_minimal=PASS")
print("two_jump_radius10_monotone_routing=4319/4319")
print("one_pole_superlevel_H0=(births=33,finite_bars=32,max_components=19)")
print("largest_rank_persistence_bar=%s" % repr(bar_by_birth[(1,)]))
print("low_hub_and_radius10_trap_bars=%s" % repr((
    bar_by_birth[hub_low], bar_by_birth[sharp_traps[0]]
)))
print("merge_tree_sha256=6f6aeb54bb02aa89cc063396f25c7883dc254c7c97339f1459d1c6125ba65ce2")
print("scope=fixed_THM3238_bank_state_dependent_row_selector_not_one_dual_or_simplicial_contraction")
print("all_exact_checks=PASS")
