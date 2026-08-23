#!/usr/bin/env python3
"""Independent exact classifier for the F12/F13 common reset graph.

The inherited response computation identifies the five common coordinate
moves.  This script forgets all response magnitudes and independently proves
that the resulting 239-state rewrite graph is terminating and confluent,
classifies its six normal-form basins, and records precisely what that raw
confluence fails to preserve.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections import Counter, defaultdict
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import factorial
from pathlib import Path


GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    if not condition:
        raise RuntimeError(label)
    GATES += 1


def polynomial_add(target: list[int], source: list[int]) -> None:
    while len(target) < len(source):
        target.append(0)
    for index, coefficient in enumerate(source):
        target[index] += coefficient


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation" / "gmc_complete_physical_bank_unique_reset_thm3238.py"
BASE_SHA256 = "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa"

base_bytes = BASE.read_bytes()
gate(hashlib.sha256(base_bytes).hexdigest() == BASE_SHA256,
     "THM-3238 coefficient source hash")
spec = spec_from_file_location("thm3238_fc_confluence_primary", BASE)
base = module_from_spec(spec)
if spec.loader is None:
    raise RuntimeError("THM-3238 coefficient-source loader missing")
spec.loader.exec_module(base)


FACE_PARAMETERS = {"F12": (1, 2), "F13": (1, 3)}


def count_vector(state: tuple[int, ...]) -> tuple[int, ...]:
    census = Counter(state)
    return tuple(census[root] for root in range(1, 9))


def profile(face: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    a_value, b_value = FACE_PARAMETERS[face]
    poles = tuple(sorted(base.reduced_poles(1, base.BANK, a_value, b_value)[0]))
    reset = tuple(sorted(base.residual_roots(
        1, base.dominant_row(1), a_value, b_value
    )))
    return poles, reset


def l1(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(abs(a - b) for a, b in zip(left, right))


profiles = {face: profile(face) for face in FACE_PARAMETERS}
expected_profiles = {
    "F12": ((1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5), (1, 2, 2, 3, 4, 5)),
    "F13": ((1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8),
            (1, 3, 3, 4, 5, 6, 7, 8)),
}
gate(profiles == expected_profiles, "exact inherited F12/F13 profiles")
face_capacities = {
    face: count_vector(poles) for face, (poles, _reset) in profiles.items()
}
face_resets = {
    face: count_vector(reset) for face, (_poles, reset) in profiles.items()
}
common_capacity = tuple(
    min(face_capacities[face][index] for face in FACE_PARAMETERS)
    for index in range(8)
)
gate(common_capacity == (4, 3, 2, 1, 1, 0, 0, 0),
     "derived common physical box")


def counts_to_state(counts: tuple[int, int, int, int, int]) -> tuple[int, ...]:
    return tuple(
        root
        for root, count in enumerate(counts, start=1)
        for _ in range(count)
    )


def successors(counts: tuple[int, int, int, int, int]) -> tuple[tuple[int, ...], ...]:
    c1, c2, c3, c4, c5 = counts
    answer: list[tuple[int, ...]] = []
    if c1 != 1:
        answer.append((c1 + (1 if c1 == 0 else -1), c2, c3, c4, c5))
    if c2 == 3:
        answer.append((c1, 2, c3, c4, c5))
    if c3 == 0:
        answer.append((c1, c2, 1, c4, c5))
    if c4 == 0:
        answer.append((c1, c2, c3, 1, c5))
    if c5 == 0:
        answer.append((c1, c2, c3, c4, 1))
    return tuple(answer)


def face_successors(face: str,
                    counts: tuple[int, int, int, int, int]) -> tuple[tuple[int, ...], ...]:
    census = counts + (0, 0, 0)
    capacity = face_capacities[face]
    reset = face_resets[face]
    old_distance = l1(census, reset)
    answer: list[tuple[int, ...]] = []
    for index in range(8):
        for step in (-1, 1):
            changed = list(census)
            changed[index] += step
            target = tuple(changed)
            if not 0 <= target[index] <= capacity[index] or not any(target):
                continue
            if l1(target, reset) == old_distance - 1:
                answer.append(target[:5])
    return tuple(answer)


def normal_form(counts: tuple[int, int, int, int, int]) -> tuple[int, ...]:
    _, c2, c3, _, _ = counts
    return (1, 2 if c2 == 3 else c2, 1 if c3 == 0 else c3, 1, 1)


def root_one_debt(c1: int) -> int:
    return 1 if c1 == 0 else max(c1 - 1, 0)


def debt(counts: tuple[int, int, int, int, int]) -> int:
    c1, c2, c3, c4, c5 = counts
    return (
        root_one_debt(c1)
        + (c2 == 3)
        + (c3 == 0)
        + (c4 == 0)
        + (c5 == 0)
    )


states = tuple(
    census
    for census in product(*(range(capacity + 1) for capacity in common_capacity[:5]))
    if any(census)
)
state_set = set(states)
gate(len(states) == 239, "shared physical state count")

edges: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
outdegree = Counter()
basins: dict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
debt_histogram = Counter()
path_count_histogram = Counter()

for state in states:
    next_states = successors(state)
    common_face_successors = tuple(
        target for target in face_successors("F12", state)
        if target in set(face_successors("F13", state))
    )
    gate(next_states == common_face_successors,
         f"coordinate rule equals derived face intersection at {state}")
    gate(all(target in state_set for target in next_states),
         f"successor remains physical at {state}")
    gate(all(debt(target) == debt(state) - 1 for target in next_states),
         f"strict debt descent at {state}")
    outdegree[len(next_states)] += 1
    edges.extend((state, target) for target in next_states)
    terminal = normal_form(state)
    gate(terminal in state_set, f"normal form exists at {state}")
    gate(not successors(terminal), f"normal form is terminal at {state}")
    gate(all(normal_form(target) == terminal for target in next_states),
         f"normal form invariant on every outgoing edge at {state}")
    for left, right in combinations(next_states, 2):
        common_next = set(successors(left)) & set(successors(right))
        gate(len(common_next) == 1,
             f"immediate coordinate fork has a unique diamond join at {state}")
        join = next(iter(common_next))
        gate(normal_form(join) == terminal,
             f"fork join retains the common normal form at {state}")
    basins[terminal].append(state)
    debt_histogram[debt(state)] += 1

    c1, c2, c3, c4, c5 = state
    a = root_one_debt(c1)
    one_shot = int(c2 == 3) + int(c3 == 0) + int(c4 == 0) + int(c5 == 0)
    predicted_paths = factorial(a + one_shot) // factorial(a)

    # Dynamic recursion is acyclic because every edge lowers debt.  Comparing
    # the closed form with the successor sum is a path-count recurrence check;
    # confluence itself is checked above by edgewise NF invariance and fork joins.
    successor_path_sum = sum(
        factorial(root_one_debt(target[0])
                  + int(target[1] == 3) + int(target[2] == 0)
                  + int(target[3] == 0) + int(target[4] == 0))
        // factorial(root_one_debt(target[0]))
        for target in next_states
    )
    if next_states:
        gate(predicted_paths == successor_path_sum,
             f"path recursion at {state}")
    else:
        gate(predicted_paths == 1, f"terminal path count at {state}")
    path_count_histogram[predicted_paths] += 1

expected_sinks = {
    (1, c2, c3, 1, 1)
    for c2 in (0, 1, 2)
    for c3 in (1, 2)
}
gate(set(basins) == expected_sinks, "six normal forms")

basin_sizes = {sink: len(members) for sink, members in basins.items()}
expected_basin_sizes = {
    (1, 0, 1, 1, 1): 39,
    (1, 0, 2, 1, 1): 20,
    (1, 1, 1, 1, 1): 40,
    (1, 1, 2, 1, 1): 20,
    (1, 2, 1, 1, 1): 80,
    (1, 2, 2, 1, 1): 40,
}
gate(basin_sizes == expected_basin_sizes, "normal-form basin sizes")
gate(sum(basin_sizes.values()) == 239, "basins partition state set")

expected_outdegree = {0: 6, 1: 41, 2: 85, 3: 75, 4: 28, 5: 4}
gate(dict(sorted(outdegree.items())) == expected_outdegree,
     "outdegree polynomial")
gate(len(edges) == 568, "common graph edge count")

expected_debt = {0: 6, 1: 29, 2: 57, 3: 64, 4: 48, 5: 26, 6: 8, 7: 1}
gate(dict(sorted(debt_histogram.items())) == expected_debt,
     "graded debt polynomial")
gate(max(debt_histogram) == 7, "sharp rewrite depth seven")

max_path_count = max(path_count_histogram)
max_path_states = tuple(state for state in states if (
    factorial(
        root_one_debt(state[0])
        + int(state[1] == 3) + int(state[2] == 0)
        + int(state[3] == 0) + int(state[4] == 0)
    ) // factorial(root_one_debt(state[0]))
) == max_path_count)
gate(max_path_count == 840, "sharp path multiplicity 840")
gate(max_path_states == ((4, 3, 0, 0, 0),),
     "unique maximal-path state")
expected_path_histogram = {
    1: 47, 2: 51, 3: 17, 4: 17, 6: 41, 12: 17,
    20: 17, 24: 14, 60: 7, 120: 9, 360: 1, 840: 1,
}
gate(dict(sorted(path_count_histogram.items())) == expected_path_histogram,
     "complete normalization-history distribution")
gate(sum(count for paths, count in path_count_histogram.items() if paths > 1) == 192,
     "192 states have nonunique reset chronology")

# Derive the smallest hostile response from the inherited THM-3238 exact
# coefficient source.  In particular, no response fibre is supplied as an
# input literal to the finite classifier.
n0 = (0, 0, 1, 1, 1)  # multiset (3,4,5)
n1 = (1, 0, 1, 1, 1)  # multiset (1,3,4,5)
gate(n1 in successors(n0), "smallest hostile common raw edge")


def row_values(a_value: int, b_value: int,
               state: tuple[int, ...]) -> tuple[object, ...]:
    vectors = base.coefficient_vectors(1, base.BANK, a_value, b_value, state)
    return tuple(
        sum(vectors[entry[0]][shape] for shape in upset)
        for entry, upset in zip(base.CERTIFICATE, base.UPSETS)
    )


deltas: dict[str, tuple[object, ...]] = {}
positive: dict[str, frozenset[int]] = {}
for face, (a_value, b_value) in FACE_PARAMETERS.items():
    before = row_values(a_value, b_value, counts_to_state(n0))
    after = row_values(a_value, b_value, counts_to_state(n1))
    delta = tuple(right - left for left, right in zip(before, after))
    deltas[face] = delta
    positive[face] = frozenset(
        index + 1 for index, value in enumerate(delta) if value > 0
    )

expected_positive = {
    "F12": frozenset((2, 5, 8, 9, 11, 12, 14, 16, 18, 22)),
    "F13": frozenset((3, 4, 6, 7, 10, 13, 17, 19, 20, 21)),
}
gate(positive == expected_positive,
     "derived exact face-specific positive response fibres")
gate(not (positive["F12"] & positive["F13"]),
     "common raw edge has disjoint strict-positive row fibres")
gate(positive["F12"] | positive["F13"] == set(range(1, 23)) - {1, 15},
     "twenty positive rows and two common negative rows")
gate(deltas["F12"][0] < 0 and deltas["F13"][0] < 0
     and deltas["F12"][14] < 0 and deltas["F13"][14] < 0,
     "rows 1 and 15 are negative in both faces")
equal_positive = tuple(
    (left + 1, right + 1, value_left)
    for left, value_left in enumerate(deltas["F12"]) if value_left > 0
    for right, value_right in enumerate(deltas["F13"])
    if value_right > 0 and value_left == value_right
)
gate(equal_positive == (), "derived raw-gauge equal-positive relation is empty")
gate(normal_form(n0) == n1 and not successors(n1),
     "hostile common edge lands at a normal form")

semantic = {
    "universe": "nonempty count box [0,4]x[0,3]x[0,2]x[0,1]x[0,1]",
    "moves": "c1 toward 1; c2 3->2; c3 0->1; c4,c5 0->1",
    "rewrite": "terminating graded DAG; every state has a unique normal form",
    "basins": [
        [list(sink), expected_basin_sizes[sink]]
        for sink in sorted(expected_basin_sizes)
    ],
    "debt_polynomial": expected_debt,
    "outdegree_polynomial": expected_outdegree,
    "max_paths": [840, [4, 3, 0, 0, 0]],
    "destroyed": "inherited named response row, raw response magnitude, and reset chronology",
    "hostile": "(3,4,5)->(1,3,4,5) is raw-confluent but has disjoint positive row fibres",
    "scope": "two named I2 faces only; no FC3/HFC3/SFC3 conclusion",
}
blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "no inactive Python assert")

print("probe=F12_F13_common_reset_confluent_normal_forms")
print("graph=states:239;edges:568;sinks:6;max_depth:7")
print("outdegree=0:6,1:41,2:85,3:75,4:28,5:4")
print("debt=0:6,1:29,2:57,3:64,4:48,5:26,6:8,7:1")
print("basins=39,20,40,20,80,40")
print("path_histogram=1:47,2:51,3:17,4:17,6:41,12:17,20:17,24:14,60:7,120:9,360:1,840:1")
print("chronology=192_states_nonunique;unique_max_(4,3,0,0,0):840_directed_reset_histories")
print("boundary=raw_state_confluence_supplies_no_common_inherited_named_row_or_raw_magnitude_certificate")
print("scope=two_named_I2_faces;FC3_HFC3_SFC3_OPEN")
print(f"semantic_sha256={hashlib.sha256(blob).hexdigest()}")
print(f"GATES={GATES}")
print("PASS")
