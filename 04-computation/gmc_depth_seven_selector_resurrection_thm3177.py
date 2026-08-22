#!/usr/bin/env python3
"""Exact depth-seven selector resurrection certificate for THM-3177.

For support ``(1,3)``, bank ``I2``, this companion reconstructs the complete
1,820-state physical prefix bank through depth seven and verifies an explicit
denominator-100,000,000 probability law.  In every degree 5 through 13 it
uses a separately implemented partition generator, bin-packing coarsening
test, and integer-capacity Dinic cut to minimize over *all* nontrivial
coarsening upsets.  Every resulting minimum is strictly positive.

The floating column-generation search which discovered the law is not part of
the proof.  All promoted checks below use integers and ``Fraction`` only.
"""

import ast
import hashlib
from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_pole_prefix_hasse_current_scout.py"
UPSTREAM_BANK = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
UPSTREAM_LF_SHA256 = (
    "151edb9b8cee4807d3f08dc17af32e45420021ba30dfd116c04d9fcaf8bbd5b7")
UPSTREAM_BANK_LF_SHA256 = (
    "15b94691d53afbcdc6aefda89fc7cd5497534ca70fb780a686892dabb5646d6f")


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n")


def load_upstream_prefix(maximum_degree):
    actual = hashlib.sha256(lf_bytes(UPSTREAM)).hexdigest()
    require(actual == UPSTREAM_LF_SHA256,
            ("pole-prefix helper hash drift", actual))
    bank_actual = hashlib.sha256(lf_bytes(UPSTREAM_BANK)).hexdigest()
    require(bank_actual == UPSTREAM_BANK_LF_SHA256,
            ("signed-bank helper hash drift", bank_actual))
    tree = ast.parse(UPSTREAM.read_text(encoding="utf-8"))
    prefix = []
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "MAXIMUM_DEGREE"
                        for target in node.targets)):
            node.value = ast.Constant(maximum_degree)
        prefix.append(node)
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "UNIVERSE"
                        for target in node.targets)):
            break
    module = ast.fix_missing_locations(
        ast.Module(body=prefix, type_ignores=[]))
    namespace = {"__file__": str(UPSTREAM)}
    exec(compile(module, str(UPSTREAM), "exec"), namespace)
    return namespace


UP = load_upstream_prefix(13)
UP["residual_roots"] = lru_cache(maxsize=None)(UP["residual_roots"])
UP["all_monomial_values"] = lru_cache(maxsize=None)(
    UP["all_monomial_values"])
upstream_partitions = UP["partitions"]
upstream_coarsens = UP["coarsens"]
coefficient_vectors = UP["coefficient_vectors"]
reduced_poles = UP["reduced_poles"]
BANK = UP["BANKS"][1]


# ---------------------------------------------------------------------------
# 1. Independent partition/refinement implementation.
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def independent_partitions(total, maximum=None):
    """Nonincreasing integer partitions, independent of the bank helper."""

    if total == 0:
        return ((),)
    if maximum is None or maximum > total:
        maximum = total
    answer = []
    for first in range(maximum, 0, -1):
        for tail in independent_partitions(total - first, first):
            answer.append((first,) + tail)
    return tuple(answer)


@lru_cache(maxsize=None)
def independent_coarsens(fine, coarse):
    """Whether ``fine`` can be packed exactly into bins ``coarse``."""

    if sum(fine) != sum(coarse) or len(fine) < len(coarse):
        return False
    pieces = tuple(sorted(fine, reverse=True))
    bins = tuple(sorted(coarse, reverse=True))

    @lru_cache(maxsize=None)
    def place(index, capacities):
        if index == len(pieces):
            return not any(capacities)
        piece = pieces[index]
        tried = set()
        for slot, capacity in enumerate(capacities):
            if capacity < piece or capacity in tried:
                continue
            tried.add(capacity)
            updated = list(capacities)
            updated[slot] -= piece
            updated.sort(reverse=True)
            if place(index + 1, tuple(updated)):
                return True
        return False

    return place(0, bins)


PARTITIONS = {
    degree: independent_partitions(degree)
    for degree in range(5, 14)
}
PARTITION_CHECKS = 0
COARSENING_CHECKS = 0
for degree, shapes in PARTITIONS.items():
    require(set(shapes) == set(upstream_partitions(degree)),
            ("independent partition bank mismatch", degree))
    PARTITION_CHECKS += len(shapes)
    for fine in shapes:
        for coarse in shapes:
            require(independent_coarsens(fine, coarse)
                    == upstream_coarsens(fine, coarse),
                    ("independent coarsening mismatch", degree,
                     fine, coarse))
            COARSENING_CHECKS += 1


# ---------------------------------------------------------------------------
# 2. The complete depth-seven state bank and exact law.
# ---------------------------------------------------------------------------

POLES, _ = reduced_poles(1, BANK, 1, 3)
MULTIPLICITY = Counter(POLES)
VALUES = tuple(sorted(MULTIPLICITY))
BY_DEPTH = tuple(
    tuple(
        state
        for state in combinations_with_replacement(VALUES, depth)
        if all(Counter(state)[value] <= MULTIPLICITY[value]
               for value in set(state))
    )
    for depth in range(1, 8)
)
STATES = tuple(state for layer in BY_DEPTH for state in layer)
require(POLES == (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1),
        "support-(1,3) pole multiset drift")
require(tuple(map(len, BY_DEPTH)) == (8, 33, 93, 200, 348, 507, 631)
        and len(STATES) == 1820,
        "depth-seven physical-state census drift")

DENOMINATOR = 100_000_000
LAW_ITEMS = (
    ((1,), 329341),
    ((2,), 506132),
    ((1, 1, 2), 353515),
    ((1, 1, 1, 1), 589431),
    ((5, 5, 6, 7, 8), 13469020),
    ((1, 1, 1, 1, 2, 2), 448483),
    ((1, 1, 1, 1, 2, 2, 2), 510392),
    ((1, 1, 2, 2, 2, 3, 3), 446099),
    ((1, 2, 2, 2, 3, 3, 4), 313313),
    ((1, 4, 5, 5, 6, 7, 8), 12884975),
    ((2, 3, 3, 5, 6, 7, 8), 17513294),
    ((3, 3, 4, 4, 6, 7, 8), 2642438),
    ((4, 4, 5, 5, 6, 7, 8), 49993567),
)
LAW = dict(LAW_ITEMS)
require(all(state in STATES and count > 0 for state, count in LAW_ITEMS),
        "law contains an illegal or zero state")
require(sum(LAW.values()) == DENOMINATOR
        and reduce(gcd, LAW.values()) == 1,
        "law lost normalization or primitivity")

VECTORS = {
    state: coefficient_vectors(1, BANK, 1, 3, state)
    for state in STATES
}

BALANCE_CHECKS = 0
for state in STATES:
    for degree in range(5, 14):
        require(sum(VECTORS[state][degree].values()) == 0,
                ("state response lost zero mass", state, degree))
        BALANCE_CHECKS += 1


# ---------------------------------------------------------------------------
# 3. Exact integer Dinic minimization over every coarsening upset.
# ---------------------------------------------------------------------------

class Dinic:
    def __init__(self, count):
        self.graph = [[] for _ in range(count)]

    def add_edge(self, source, target, capacity):
        require(isinstance(capacity, int) and capacity >= 0,
                "Dinic capacity is not a nonnegative integer")
        forward = [target, capacity, len(self.graph[target])]
        reverse = [source, 0, len(self.graph[source])]
        self.graph[source].append(forward)
        self.graph[target].append(reverse)

    def max_flow(self, source, sink):
        total = 0
        count = len(self.graph)
        while True:
            level = [-1] * count
            level[source] = 0
            queue = deque([source])
            while queue:
                vertex = queue.popleft()
                for target, capacity, _ in self.graph[vertex]:
                    if capacity and level[target] < 0:
                        level[target] = level[vertex] + 1
                        queue.append(target)
            if level[sink] < 0:
                return total
            cursor = [0] * count

            def send(vertex, amount):
                if vertex == sink:
                    return amount
                while cursor[vertex] < len(self.graph[vertex]):
                    edge = self.graph[vertex][cursor[vertex]]
                    target, capacity, reverse_index = edge
                    if capacity and level[target] == level[vertex] + 1:
                        pushed = send(target, min(amount, capacity))
                        if pushed:
                            edge[1] -= pushed
                            self.graph[target][reverse_index][1] += pushed
                            return pushed
                    cursor[vertex] += 1
                return 0

            infinity = sum(edge[1] for edges in self.graph for edge in edges)
            while True:
                pushed = send(source, infinity)
                if not pushed:
                    break
                total += pushed

    def reachable(self, source):
        seen = {source}
        queue = deque([source])
        while queue:
            vertex = queue.popleft()
            for target, capacity, _ in self.graph[vertex]:
                if capacity and target not in seen:
                    seen.add(target)
                    queue.append(target)
        return seen


def minimal_generators(degree, upset):
    return tuple(
        shape for shape in PARTITIONS[degree]
        if shape in upset and not any(
            other != shape and other in upset
            and independent_coarsens(other, shape)
            for other in PARTITIONS[degree]
        )
    )


def minimum_nontrivial_upset(degree, rational_weights):
    shapes = PARTITIONS[degree]
    denominator = 1
    for value in rational_weights.values():
        denominator = lcm(denominator, value.denominator)
    weights = {
        shape: value.numerator * (denominator // value.denominator)
        for shape, value in rational_weights.items()
    }
    require(sum(weights.values()) == 0,
            ("weighted response lost zero mass", degree))
    index = {shape: position for position, shape in enumerate(shapes)}
    source, sink = len(shapes), len(shapes) + 1
    graph = Dinic(len(shapes) + 2)
    finite_total = sum(abs(value) for value in weights.values())
    infinity = finite_total + 1
    graph.add_edge(source, index[(degree,)], infinity)
    graph.add_edge(index[(1,) * degree], sink, infinity)
    negative_constant = 0
    for shape, value in weights.items():
        if value >= 0:
            graph.add_edge(index[shape], sink, value)
        else:
            graph.add_edge(source, index[shape], -value)
            negative_constant += -value
    for fine in shapes:
        for coarse in shapes:
            if fine != coarse and independent_coarsens(fine, coarse):
                graph.add_edge(index[fine], index[coarse], infinity)
    flow = graph.max_flow(source, sink)
    reachable = graph.reachable(source)
    upset = frozenset(shape for shape in shapes
                      if index[shape] in reachable)
    require((degree,) in upset and (1,) * degree not in upset,
            ("Dinic cut lost nontrivial anchors", degree))
    require(all(not independent_coarsens(fine, coarse) or coarse in upset
                for fine in upset for coarse in shapes),
            ("Dinic cut is not an upset", degree))
    integer_mass = sum(weights[shape] for shape in upset)
    require(flow == negative_constant + integer_mass,
            ("Dinic cut identity failed", degree))
    return Fraction(integer_mass, denominator), upset


EXPECTED_MINIMA = (
    (5, 109628640, 1, ((5,),)),
    (6, 72639032313522840, 1, ((6,),)),
    (7, 2990897801796238560, 14, ((2, 1, 1, 1, 1, 1),)),
    (8, 318185095107504, 21, ((2, 1, 1, 1, 1, 1, 1),)),
    (9, 9346991684459633328, 29,
     ((2, 1, 1, 1, 1, 1, 1, 1),)),
    (10, 3288235113012380640, 41,
     ((2, 1, 1, 1, 1, 1, 1, 1, 1),)),
    (11, 2177378045882508776256, 55,
     ((2, 1, 1, 1, 1, 1, 1, 1, 1, 1),)),
    (12, 20818394379371373887136, 75,
     ((3, 1, 1, 1, 1, 1, 1, 1, 1, 1),
      (2, 2, 1, 1, 1, 1, 1, 1, 1, 1))),
    (13, 1303441821531365341151616, 96,
     ((3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
      (2, 2, 2, 2, 2, 1, 1, 1))),
)

MINIMA = []
for expected_degree, expected_raw, expected_size, expected_generators in EXPECTED_MINIMA:
    raw_weights = {
        shape: sum(
            count * VECTORS[state][expected_degree][shape]
            for state, count in LAW_ITEMS
        )
        for shape in PARTITIONS[expected_degree]
    }
    require(all(value.denominator == 1 for value in raw_weights.values()),
            ("law response lost raw integrality", expected_degree))
    mass, upset = minimum_nontrivial_upset(expected_degree, raw_weights)
    generators = minimal_generators(expected_degree, upset)
    require(mass == expected_raw and len(upset) == expected_size
            and generators == expected_generators,
            ("minimum-upset certificate drift", expected_degree,
             mass, len(upset), generators))
    require(mass > 0, ("selector law lost strictness", expected_degree))
    MINIMA.append((expected_degree, mass, len(upset), generators))


# ---------------------------------------------------------------------------
# 4. THM-3163's abstract posterior-chain sidecar, with exact census.
# ---------------------------------------------------------------------------

LABELS_BY_VALUE = defaultdict(tuple)
for label, value in enumerate(POLES):
    LABELS_BY_VALUE[value] = LABELS_BY_VALUE[value] + (label,)
VALUE_BY_LABEL = dict(enumerate(POLES))


def labelled_lifts(state):
    counts = Counter(state)
    choices = [combinations(LABELS_BY_VALUE[value], counts[value])
               for value in sorted(counts)]
    return tuple(frozenset().union(*parts) for parts in product(*choices))


LIFTS = tuple(len(labelled_lifts(state)) for state, _ in LAW_ITEMS)
require(LIFTS == (4, 3, 18, 1, 1, 3, 1, 6, 8, 8, 6, 1, 1),
        "labelled lift census drift")
TERMINAL_MASS = {}
for state, count in LAW_ITEMS:
    lifts = labelled_lifts(state)
    for terminal in lifts:
        require(terminal not in TERMINAL_MASS,
                "distinct multiplicity states shared a labelled lift")
        TERMINAL_MASS[terminal] = Fraction(
            count, DENOMINATOR * len(lifts))
require(len(TERMINAL_MASS) == 61
        and sum(TERMINAL_MASS.values()) == 1,
        "labelled terminal law drift")

REACHABLE = {frozenset()}
for terminal in TERMINAL_MASS:
    ordered = tuple(sorted(terminal))
    for size in range(len(ordered) + 1):
        REACHABLE.update(frozenset(piece)
                         for piece in combinations(ordered, size))


def reach_weight(state):
    size = len(state)
    return sum(
        mass / comb(len(terminal), size)
        for terminal, mass in TERMINAL_MASS.items()
        if state <= terminal
    )


REACH_WEIGHT = {state: reach_weight(state) for state in REACHABLE}
TRANSITIONS = {}
STOP = {}
for state in REACHABLE:
    size = len(state)
    reach = REACH_WEIGHT[state]
    require(reach > 0, "reachable labelled state has zero posterior weight")
    STOP[state] = TERMINAL_MASS.get(state, Fraction(0)) / reach
    transition = {}
    for label in range(len(POLES)):
        if label in state:
            continue
        numerator = sum(
            mass / (comb(len(terminal), size) * (len(terminal) - size))
            for terminal, mass in TERMINAL_MASS.items()
            if state | {label} <= terminal
        )
        if numerator:
            transition[label] = numerator / reach
    require(STOP[state] + sum(transition.values()) == 1,
            ("posterior kernel lost stochasticity", state))
    TRANSITIONS[state] = transition

ARRIVAL = defaultdict(Fraction)
ARRIVAL[frozenset()] = Fraction(1)
RECOVERED = defaultdict(Fraction)
for state in sorted(REACHABLE, key=lambda item: (len(item), tuple(item))):
    require(ARRIVAL[state] == REACH_WEIGHT[state],
            ("posterior reach induction failed", state))
    RECOVERED[state] += ARRIVAL[state] * STOP[state]
    for label, probability in TRANSITIONS[state].items():
        ARRIVAL[state | {label}] += ARRIVAL[state] * probability
require(all(RECOVERED[state] == TERMINAL_MASS.get(state, 0)
            for state in REACHABLE)
        and sum(RECOVERED.values()) == 1,
        "posterior chain did not recover the terminal law")


def multiplicity_state(labels):
    return tuple(sorted(VALUE_BY_LABEL[label] for label in labels))


LUMPED_STATES = {multiplicity_state(state) for state in REACHABLE}
LUMPED_EDGES = set()
LUMP_SIGNATURES = {}
for state in REACHABLE:
    source_type = multiplicity_state(state)
    signature = defaultdict(Fraction)
    for label, probability in TRANSITIONS[state].items():
        target_type = multiplicity_state(state | {label})
        signature[VALUE_BY_LABEL[label]] += probability
        LUMPED_EDGES.add((source_type, target_type))
    signature = (STOP[state], tuple(sorted(signature.items())))
    if source_type in LUMP_SIGNATURES:
        require(signature == LUMP_SIGNATURES[source_type],
                ("exchangeable lift failed strong lumpability", source_type))
    else:
        LUMP_SIGNATURES[source_type] = signature

LABELLED_EDGE_COUNT = sum(len(edges) for edges in TRANSITIONS.values())
require((len(REACHABLE), LABELLED_EDGE_COUNT,
         len(LUMPED_STATES), len(LUMPED_EDGES))
        == (1620, 6692, 289, 908),
        "posterior-chain state/edge census drift")


print("THM-3177 depth-seven degree-thirteen strict selector resurrection")
print("pole_prefix_dependency_sha256=" + UPSTREAM_LF_SHA256)
print("signed_bank_dependency_sha256=" + UPSTREAM_BANK_LF_SHA256)
print("poles=" + repr(POLES))
print("state_counts_depth1_depth7_total="
      + repr(tuple(map(len, BY_DEPTH)) + (len(STATES),)))
print("independent_partition_checks=" + repr(PARTITION_CHECKS))
print("independent_coarsening_checks=" + repr(COARSENING_CHECKS))
print("state_degree_zero_mass_checks=" + repr(BALANCE_CHECKS))
print("law_denominator=" + repr(DENOMINATOR))
print("law_items=" + repr(LAW_ITEMS))
print("minimum_nontrivial_upsets=" + repr(tuple(MINIMA)))
print("all_nontrivial_upsets_strict=1")
print("labelled_lift_counts=" + repr(LIFTS))
print("posterior_chain_counts="
      + repr((len(TERMINAL_MASS), len(REACHABLE), LABELLED_EDGE_COUNT,
              len(LUMPED_STATES), len(LUMPED_EDGES))))
print("barcode=(C13_depth6_empty,C13_depth7_nonempty)")
print("scope=fixed_support_bank_averaged_virtual_prefix_currents")
print("abstract_markov_sidecar=response_compatibility_not_proved")
print("all_exact_checks=PASS")
