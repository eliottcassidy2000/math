#!/usr/bin/env python3
"""Independent hostile checker for the F12/F13 common raw reset graph.

Profiles and responses are reconstructed from the older THM-3238 coefficient
source.  The candidate is never imported.  Graph normal forms and path counts
are computed recursively rather than assumed from the candidate formulas.
"""

from __future__ import annotations

import ast
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
import json
from math import factorial
from pathlib import Path


GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    if condition is not True:
        raise RuntimeError(label)
    GATES += 1


ROOT = Path(__file__).resolve().parents[1]
if not (ROOT / "04-computation").is_dir():
    ROOT = ROOT.parent
CANDIDATE = (
    ROOT / "04-computation" / "fc_common_reset_confluent_normal_form_thm3732.py"
)
if not CANDIDATE.is_file():
    CANDIDATE = ROOT / ".scratch" / "fc_common_reset_confluent_normal_form_20260823.py"
BASE = ROOT / "04-computation" / "gmc_complete_physical_bank_unique_reset_thm3238.py"
CANDIDATE_SHA256 = "e24b49dce9d8a99b318f5344339df63f743fac0c5bfaec7f9b3365f8e47e7e19"
BASE_SHA256 = "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa"

candidate_bytes = CANDIDATE.read_bytes()
base_bytes = BASE.read_bytes()
gate(sha256(candidate_bytes).hexdigest() == CANDIDATE_SHA256, "candidate hash")
gate(sha256(base_bytes).hexdigest() == BASE_SHA256, "THM-3238 base hash")
gate(not any(isinstance(node, ast.Assert)
             for node in ast.walk(ast.parse(candidate_bytes.decode("utf-8")))),
     "candidate is assertion-free")

spec = spec_from_file_location("thm3238_fc_common_audit", BASE)
base = module_from_spec(spec)
if spec.loader is None:
    raise RuntimeError("THM-3238 loader missing")
spec.loader.exec_module(base)


@lru_cache(maxsize=None)
def partitions(total: int, maximum: int | None = None) -> tuple[tuple[int, ...], ...]:
    if total == 0:
        return ((),)
    bound = total if maximum is None else min(total, maximum)
    return tuple(
        (head,) + tail
        for head in range(bound, 0, -1)
        for tail in partitions(total - head, head)
    )


@lru_cache(maxsize=None)
def coarsens(fine: tuple[int, ...], coarse: tuple[int, ...]) -> bool:
    if sum(fine) != sum(coarse) or len(fine) < len(coarse):
        return False
    pieces = tuple(sorted(fine, reverse=True))
    bins = tuple(sorted(coarse, reverse=True))

    @lru_cache(maxsize=None)
    def place(index: int, capacities: tuple[int, ...]) -> bool:
        if index == len(pieces):
            return not any(capacities)
        piece = pieces[index]
        tried: set[int] = set()
        for slot, capacity in enumerate(capacities):
            if capacity < piece or capacity in tried:
                continue
            tried.add(capacity)
            changed = list(capacities)
            changed[slot] -= piece
            changed.sort(reverse=True)
            if place(index + 1, tuple(changed)):
                return True
        return False

    return place(0, bins)


UPSETS = tuple(
    frozenset(
        shape for shape in partitions(degree)
        if any(coarsens(generator, shape) for generator in generators)
    )
    for degree, _expected, generators, _factor in base.CERTIFICATE
)
gate(tuple(map(len, UPSETS)) == tuple(row[1] for row in base.CERTIFICATE),
     "fresh upset reconstruction")


def row_values(a_value: int, b_value: int, state: tuple[int, ...]) -> tuple[int, ...]:
    vectors = base.coefficient_vectors(1, base.BANK, a_value, b_value, state)
    return tuple(
        sum(vectors[entry[0]][shape] for shape in upset)
        for entry, upset in zip(base.CERTIFICATE, UPSETS)
    )


FACE_PARAMETERS = {"F12": (1, 2), "F13": (1, 3)}


def profile(face: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    a_value, b_value = FACE_PARAMETERS[face]
    poles = tuple(sorted(base.reduced_poles(1, base.BANK, a_value, b_value)[0]))
    reset = tuple(sorted(base.residual_roots(
        1, base.dominant_row(1), a_value, b_value
    )))
    return poles, reset


def counts(state: tuple[int, ...]) -> tuple[int, ...]:
    census = Counter(state)
    return tuple(census[root] for root in range(1, 9))


def state_from_counts(census: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        root for root, multiplicity in enumerate(census, 1)
        for _ in range(multiplicity)
    )


def l1(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(abs(a - b) for a, b in zip(left, right))


profiles = {face: profile(face) for face in FACE_PARAMETERS}
expected_profiles = {
    "F12": ((1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5), (1, 2, 2, 3, 4, 5)),
    "F13": ((1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8),
            (1, 3, 3, 4, 5, 6, 7, 8)),
}
gate(profiles == expected_profiles, "exact F12/F13 profiles")

face_capacities = {
    face: counts(poles) for face, (poles, _reset) in profiles.items()
}
face_resets = {
    face: counts(reset) for face, (_poles, reset) in profiles.items()
}
common_capacity = tuple(min(face_capacities[face][index] for face in FACE_PARAMETERS)
                        for index in range(8))
gate(common_capacity == (4, 3, 2, 1, 1, 0, 0, 0), "common multiplicity box")

states = tuple(
    census for census in product(*(range(capacity + 1) for capacity in common_capacity))
    if any(census)
)
state_set = set(states)
gate(len(states) == 239, "nonempty common state universe")


def face_successors(face: str, census: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
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
                answer.append(target)
    return tuple(answer)


def coordinate_rule(census: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    c1, c2, c3, c4, c5, _c6, _c7, _c8 = census
    answer: list[tuple[int, ...]] = []
    moves = (
        (0, 1 if c1 == 0 else -1 if c1 > 1 else 0),
        (1, -1 if c2 == 3 else 0),
        (2, 1 if c3 == 0 else 0),
        (3, 1 if c4 == 0 else 0),
        (4, 1 if c5 == 0 else 0),
    )
    for index, step in moves:
        if step:
            changed = list(census)
            changed[index] += step
            answer.append(tuple(changed))
    return tuple(answer)


graph: dict[tuple[int, ...], tuple[tuple[int, ...], ...]] = {}
edges: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
for census in states:
    right = set(face_successors("F13", census))
    common = tuple(target for target in face_successors("F12", census)
                   if target in right)
    gate(common == coordinate_rule(census), f"derived coordinate rule {census}")
    gate(all(target in state_set for target in common), f"physical successor {census}")
    graph[census] = common
    edges.extend((census, target) for target in common)


def c1_debt(value: int) -> int:
    return 1 if value == 0 else max(value - 1, 0)


def debt(census: tuple[int, ...]) -> int:
    c1, c2, c3, c4, c5, _c6, _c7, _c8 = census
    return (c1_debt(c1) + int(c2 == 3) + int(c3 == 0)
            + int(c4 == 0) + int(c5 == 0))


for source, target in edges:
    gate(debt(target) == debt(source) - 1, f"strict debt descent {source}->{target}")

outdegree_histogram = Counter(len(graph[state]) for state in states)
debt_histogram = Counter(debt(state) for state in states)
gate(len(edges) == 568, "edge count")
gate(dict(sorted(outdegree_histogram.items()))
     == {0: 6, 1: 41, 2: 85, 3: 75, 4: 28, 5: 4},
     "outdegree distribution")
gate(dict(sorted(debt_histogram.items()))
     == {0: 6, 1: 29, 2: 57, 3: 64, 4: 48, 5: 26, 6: 8, 7: 1},
     "debt distribution")


@lru_cache(maxsize=None)
def reachable_sinks(census: tuple[int, ...]) -> frozenset[tuple[int, ...]]:
    if not graph[census]:
        return frozenset((census,))
    return frozenset().union(*(reachable_sinks(target) for target in graph[census]))


@lru_cache(maxsize=None)
def reachable_states(census: tuple[int, ...]) -> frozenset[tuple[int, ...]]:
    return frozenset((census,)).union(
        *(reachable_states(target) for target in graph[census])
    )


@lru_cache(maxsize=None)
def path_count(census: tuple[int, ...]) -> int:
    if not graph[census]:
        return 1
    return sum(path_count(target) for target in graph[census])


@lru_cache(maxsize=None)
def max_remaining_length(census: tuple[int, ...]) -> int:
    if not graph[census]:
        return 0
    return 1 + max(max_remaining_length(target) for target in graph[census])


normal_forms: dict[tuple[int, ...], tuple[int, ...]] = {}
for census in states:
    sinks = reachable_sinks(census)
    gate(len(sinks) == 1, f"unique dynamically computed normal form {census}")
    normal = next(iter(sinks))
    normal_forms[census] = normal
    gate(debt(census) == max_remaining_length(census), f"sharp debt length {census}")
    for target in graph[census]:
        gate(normal_forms.get(target, next(iter(reachable_sinks(target)))) == normal,
             f"normal form invariant on edge {census}->{target}")
    for left, right in combinations(graph[census], 2):
        gate(normal in reachable_states(left) & reachable_states(right),
             f"local fork joins at normal form {census}")

sinks = {state for state in states if not graph[state]}
expected_sinks = {
    (1, c2, c3, 1, 1, 0, 0, 0)
    for c2 in (0, 1, 2) for c3 in (1, 2)
}
gate(sinks == expected_sinks, "six raw normal forms")

basins: dict[tuple[int, ...], int] = Counter(normal_forms.values())
expected_basins = {
    (1, 0, 1, 1, 1, 0, 0, 0): 39,
    (1, 0, 2, 1, 1, 0, 0, 0): 20,
    (1, 1, 1, 1, 1, 0, 0, 0): 40,
    (1, 1, 2, 1, 1, 0, 0, 0): 20,
    (1, 2, 1, 1, 1, 0, 0, 0): 80,
    (1, 2, 2, 1, 1, 0, 0, 0): 40,
}
gate(dict(basins) == expected_basins, "normal-form basin sizes")

path_histogram = Counter(path_count(state) for state in states)
expected_path_histogram = {
    1: 47, 2: 51, 3: 17, 4: 17, 6: 41, 12: 17,
    20: 17, 24: 14, 60: 7, 120: 9, 360: 1, 840: 1,
}
gate(dict(sorted(path_histogram.items())) == expected_path_histogram,
     "complete normalizing-path distribution")

for census in states:
    c1, c2, c3, c4, c5, _c6, _c7, _c8 = census
    a = c1_debt(c1)
    b = int(c2 == 3) + int(c3 == 0) + int(c4 == 0) + int(c5 == 0)
    closed = factorial(a + b) // factorial(a)
    gate(path_count(census) == closed, f"closed path formula {census}")

max_paths = max(path_histogram)
max_states = tuple(state for state in states if path_count(state) == max_paths)
max_state = (4, 3, 0, 0, 0, 0, 0, 0)
gate(max_paths == 840 and max_states == (max_state,), "unique 840-path state")
gate(sum(count for paths, count in path_histogram.items() if paths > 1) == 192,
     "192 states have nonunique chronology")


def first_two_paths(census: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
    if not graph[census]:
        return ((census,),)
    answer: list[tuple[tuple[int, ...], ...]] = []
    for target in graph[census]:
        for tail in first_two_paths(target):
            answer.append((census,) + tail)
            if len(answer) == 2:
                return tuple(answer)
    return tuple(answer)


two_histories = first_two_paths(max_state)
gate(len(two_histories) == 2 and two_histories[0] != two_histories[1],
     "distinct reset chronologies exist")
gate(two_histories[0][-1] == two_histories[1][-1] == normal_forms[max_state],
     "distinct chronologies join only after forgetting history")


# Reconstruct the inherited THM-3732 raw-edge response hostile from the older
# coefficient source, rather than accepting the candidate's hard-coded sets.
n0_state = (3, 4, 5)
n1_state = (1, 3, 4, 5)
n0 = counts(n0_state)
n1 = counts(n1_state)
gate(n1 in graph[n0] and not graph[n1], "hostile edge lands at raw normal form")

deltas: dict[str, tuple[int, ...]] = {}
positive: dict[str, frozenset[int]] = {}
for face, (a_value, b_value) in FACE_PARAMETERS.items():
    before = row_values(a_value, b_value, n0_state)
    after = row_values(a_value, b_value, n1_state)
    delta = tuple(right - left for left, right in zip(before, after))
    deltas[face] = delta
    positive[face] = frozenset(index + 1 for index, value in enumerate(delta) if value > 0)

expected_positive = {
    "F12": frozenset((2, 5, 8, 9, 11, 12, 14, 16, 18, 22)),
    "F13": frozenset((3, 4, 6, 7, 10, 13, 17, 19, 20, 21)),
}
gate(positive == expected_positive, "exact face-specific positive response fibres")
gate(not (positive["F12"] & positive["F13"]),
     "named-row strict-positive intersection is empty")
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
gate(equal_positive == (), "raw-gauge equal-positive relation is empty")

semantic = {
    "candidate_sha256": CANDIDATE_SHA256,
    "universe": "intersection of exact F12/F13 multiplicity boxes; 239 nonempty states",
    "graph": "568 common L1-reset-descending edges; debt-graded terminating DAG",
    "normal_forms": "six; dynamic unique sink for every state; basins 39,20,40,20,80,40",
    "histograms": {
        "outdegree": dict(sorted(outdegree_histogram.items())),
        "debt": dict(sorted(debt_histogram.items())),
        "paths": dict(sorted(path_histogram.items())),
    },
    "max_paths": "840 at unique count state (4,3,0,0,0)",
    "response_hostile": "(3,4,5)->(1,3,4,5); disjoint named-row positive fibres; no equal raw increment",
    "typing": "raw graph confluence only; response certificates and reset histories are forgotten",
    "scope": "two named bank-I2 faces; no FC3/HFC3/SFC3 conclusion",
}
blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert)
             for node in ast.walk(ast.parse(source))),
     "independent checker is assertion-free")

print("audit=F12_F13_common_raw_reset_confluence")
print("graph=states:239;edges:568;sinks:6;max_depth:7")
print("outdegree=0:6,1:41,2:85,3:75,4:28,5:4")
print("debt=0:6,1:29,2:57,3:64,4:48,5:26,6:8,7:1")
print("basins=39,20,40,20,80,40")
print("path_histogram=1:47,2:51,3:17,4:17,6:41,12:17,20:17,24:14,60:7,120:9,360:1,840:1")
print("chronology=192_states_nonunique;unique_max_(4,3,0,0,0):840_directed_reset_histories")
print("response=hostile_raw_normalizing_edge_has_no_common_inherited_named_row_or_raw_magnitude_certificate")
print("scope=raw_graph_only;certificate_chronology_FC3_HFC3_SFC3_not_implied")
print(f"semantic_sha256={sha256(blob).hexdigest()}")
print(f"GATES={GATES}")
print("RESULT=PASS_CANON_READY")
