#!/usr/bin/env python3
"""Exact order-certificate audit for THM-3244's reset Rips boundary."""

import ast
import base64
import hashlib
import struct
import zlib
from collections import Counter
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
print("scope=fixed_THM3238_order_not_local_deletion_flow_or_simplicial_contraction")
print("all_exact_checks=PASS")
