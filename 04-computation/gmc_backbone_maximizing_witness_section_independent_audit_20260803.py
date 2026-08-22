#!/usr/bin/env python3
"""Independent hostile audit of the frozen THM-3287/3288 scouts.

This script does not import or execute either scout.  It reconstructs the
eleven-by-twenty-two reset-link response bank from THM-3238, solves the source
blend inequalities directly in both orientations, and then rebuilds all path,
section, selector-cut, and lifted-walk claims.  Its convention is explicit:

    an ordered row edge u -> v carries lambda*f_u + f_v;
    a state s traps Q exactly when lambda*f_u(s) + f_v(s) >= 0.

Thus a large-ratio witness has f_u(s)>0>f_v(s), while a small-ratio witness
has f_v(s)>0>f_u(s).  The maximizing relation points from the former to the
latter.  These relation arrows are static witness matchings, not physical
one-pole transitions.
"""

from __future__ import annotations

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import product
from math import gcd, lcm
from pathlib import Path

import sympy as sp


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
THM3238 = ROOT / "04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py"
THM3254 = ROOT / "04-computation/gmc_first_shell_pair_clutch_thm3254.py"
THM3277 = ROOT / "04-computation/gmc_weighted_critical_phase_path_atlas_scout_20260803.py"
THM3278 = ROOT / "04-computation/fc3_selector_origin_bipartition_phase_bridge_thm3278.py"
CANDIDATE = ROOT / "04-computation/gmc_backbone_maximizing_witness_section_scout_20260803.py"
CANDIDATE_OUT = ROOT / "05-knowledge/results/gmc_backbone_maximizing_witness_section_scout_20260803.out"
RECURRENCE_SCOUT = ROOT / "04-computation/gmc_maximizing_witness_lifted_walk_rational_series_thm3288.py"
RECURRENCE_OUT = ROOT / "05-knowledge/results/gmc_maximizing_witness_lifted_walk_rational_series_thm3288.out"

DEPENDENCIES = {
    THM3238: "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa",
    ROOT / "05-knowledge/results/gmc_complete_physical_bank_unique_reset_thm3238.out":
        "77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b",
    THM3254: "05efd37eeedeca7e3be581977a894592a7873d94a966f06d9533482cc8498fee",
    ROOT / "05-knowledge/results/gmc_first_shell_pair_clutch_thm3254.out":
        "dd415c8ce6e2e196c115421d3508addabb724305a16843509264a8b3205beee9",
    THM3277: "0beca08f9214bd6befeafdc0ccc7648be33dd8c1c2d2f12560d2e56708f2cfcf",
    ROOT / "05-knowledge/results/gmc_weighted_critical_phase_path_atlas_scout_20260803.out":
        "a55cf2bb5c130d06aaae1c84dc70d2e78e570d82253c8fdb624c4f9d86ede2e0",
    THM3278: "07cf5cc1056bdd978a3a1d1146fee8139bef9a84e3b83aad08a0522b4df63a0f",
    ROOT / "05-knowledge/results/fc3_selector_origin_bipartition_phase_bridge_thm3278.out":
        "5e67eec14698f5258d17c5ba1ed295de3b9ff355aa279adb6d523ce3df683c7a",
    CANDIDATE: "aed89d67ec7acabfe5b4feae4a83f7c57b78053928be44a6cc319d81fa4a9cc6",
    CANDIDATE_OUT: "89200bc6cff7284dd33f636352f4f7f56294d90bcd2902fa96092d9f967f5fe3",
    RECURRENCE_SCOUT: "f442ff0dfd999fc314f420aacfc7102ba2ba17e033d3fd13e91407727af8a8da",
    RECURRENCE_OUT: "44a85f317bf9cb35c437e468633a88014de268cc3fb452836960655c4b764951",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected in DEPENDENCIES.items():
    require(hashlib.sha256(lf_bytes(dependency)).hexdigest() == expected,
            ("dependency hash drift", dependency.name))

source_syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(source_syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(source_syntax)
)
require(assert_nodes == 0 and float_literals == 0,
        "optimization-sensitive or floating audit source")


def literal_assignment(tree, name):
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == name
                        for target in node.targets)):
            return ast.literal_eval(node.value)
    raise RuntimeError(("missing literal assignment", name))


def unique_transcript_suffix(path, prefix):
    matches = [
        line[len(prefix):]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, ("transcript key drift", path.name, prefix))
    return matches[0]


# Load only THM-3238's exact pre-main definitions.  In particular, no trap,
# clutch, row relation, path lift, global section, or transfer matrix is loaded.
thm3238_tree = ast.parse(THM3238.read_text(encoding="utf-8"))
thm3238_prefix = []
for node in thm3238_tree.body:
    if isinstance(node, ast.FunctionDef) and node.name == "main":
        break
    thm3238_prefix.append(node)
thm3238_namespace = {
    "__file__": str(THM3238),
    "__name__": "thm3287_independent_response_prefix",
}
exec(
    compile(
        ast.fix_missing_locations(ast.Module(body=thm3238_prefix, type_ignores=[])),
        str(THM3238),
        "exec",
    ),
    thm3238_namespace,
)

reset = thm3238_namespace["RESET"]
values = thm3238_namespace["VALUES"]
capacities = thm3238_namespace["MULTIPLICITY"]
certificate = thm3238_namespace["CERTIFICATE"]
upsets = thm3238_namespace["UPSETS"]
bank = thm3238_namespace["BANK"]
coefficient_vectors = thm3238_namespace["coefficient_vectors"]
reset_counts = Counter(reset)
require(reset == (1, 3, 3, 4, 5, 6, 7, 8), "reset drift")


def physical(state):
    counts = Counter(state)
    return bool(state) and all(counts[value] <= capacities[value] for value in values)


# Build the link by literal admissible add/delete operations, rather than by
# filtering the upstream state enumeration as the frozen candidate does.
add_states = []
remove_states = []
for value in values:
    if reset_counts[value] < capacities[value]:
        add_states.append(tuple(sorted(reset + (value,))))
    if reset_counts[value] > 0:
        changed = list(reset)
        changed.remove(value)
        remove_states.append(tuple(sorted(changed)))
link_states = tuple(sorted(add_states + remove_states))
require(len(add_states) == 4 and len(remove_states) == 7 and len(link_states) == 11,
        "independent reset-link census")
require(len(set(link_states)) == 11 and all(physical(state) for state in link_states),
        "reset-link construction failed")


def count_l1(left, right):
    left_counts = Counter(left)
    right_counts = Counter(right)
    return sum(abs(left_counts[value] - right_counts[value]) for value in values)


require(all(count_l1(state, reset) == 1 for state in link_states),
        "link state is not one pole edit from Q")


def response_values(state):
    vector = coefficient_vectors(1, bank, 1, 3, state)
    exact = tuple(
        sum((vector[degree][shape] for shape in upset), Fraction(0))
        for (degree, _, _, _), upset in zip(certificate, upsets)
    )
    require(len(exact) == 22 and all(value.denominator == 1 for value in exact),
            ("nonintegral response row", state))
    return tuple(value.numerator for value in exact)


responses = {state: response_values(state) for state in link_states}
require(response_values(reset) == (0,) * 22, "Q response is not zero")
response_bank_digest = hashlib.sha256(repr(tuple(
    (state, responses[state]) for state in link_states
)).encode("ascii")).hexdigest()


def edit_name(state):
    counts = Counter(state)
    added = tuple(
        value for value in values
        for _ in range(max(0, counts[value] - reset_counts[value]))
    )
    removed = tuple(
        value for value in values
        for _ in range(max(0, reset_counts[value] - counts[value]))
    )
    if len(added) == 1 and not removed:
        return "Q+%d" % added[0]
    if len(removed) == 1 and not added:
        return "Q-%d" % removed[0]
    raise RuntimeError(("state is not in the reset link", state))


require(tuple(sorted(map(edit_name, link_states)))
        == ("Q+1", "Q+2", "Q+4", "Q+5", "Q-1", "Q-3", "Q-4",
            "Q-5", "Q-6", "Q-7", "Q-8"),
        "reset-link edit alphabet drift")


# Recover only the proved row graph and atlas literals from their companions.
thm3254_tree = ast.parse(THM3254.read_text(encoding="utf-8"))
covering_pairs = tuple(literal_assignment(thm3254_tree, "ROW_COVERING_PAIRS"))
delayed_pairs = set(literal_assignment(thm3254_tree, "DELAYED_WITNESSES"))
link_edges = tuple(edge for edge in covering_pairs if edge not in delayed_pairs)
core_edges = tuple(edge for edge in link_edges if 14 not in edge)
core_vertices = tuple(sorted({vertex for edge in core_edges for vertex in edge}))
require((len(covering_pairs), len(delayed_pairs), len(link_edges), len(core_edges),
         len(core_vertices)) == (31, 8, 23, 22, 12),
        "core graph census drift")
require(all(left < right for left, right in core_edges), "core orientation gauge drift")

thm3277_tree = ast.parse(THM3277.read_text(encoding="utf-8"))
selected_tree = tuple(literal_assignment(thm3277_tree, "expected_best_edges"))
phase = dict(literal_assignment(thm3277_tree, "expected_phase"))
proved_minimum_paths = dict(literal_assignment(thm3277_tree, "expected_paths"))
backbone = tuple(literal_assignment(thm3277_tree, "expected_atlas_edges"))
require(len(selected_tree) == 11 and len(backbone) == 7,
        "selected-tree/backbone census drift")
require(set(backbone) <= set(selected_tree) <= set(core_edges),
        "proved edge containment drift")
require(set(phase) == set(core_vertices), "phase domain drift")

cut_line = unique_transcript_suffix(
    ROOT / "05-knowledge/results/fc3_selector_origin_bipartition_phase_bridge_thm3278.out",
    "small_core=",
)
small_text, full_text = cut_line.split(";full_core=", 1)
selector_small = frozenset(ast.literal_eval(small_text))
selector_full = frozenset(ast.literal_eval(full_text))
require(selector_small.isdisjoint(selector_full)
        and selector_small | selector_full == set(core_vertices)
        and all((left in selector_small) != (right in selector_small)
                for left, right in core_edges),
        "THM-3278 cut is not the full core bipartition")


def solve_ordered_edge(first_row, second_row):
    """Solve lambda*f_first+f_second >= 0 directly on every link state."""
    large = {}
    small = {}
    whole_line = []
    for state in link_states:
        first = responses[state][first_row - 1]
        second = responses[state][second_row - 1]
        if first > 0 and second < 0:
            large[state] = Fraction(-second, first)
        elif first < 0 and second > 0:
            small[state] = Fraction(second, -first)
        elif first >= 0 and second >= 0:
            whole_line.append(state)
        else:
            # Both negative, or a negative entry paired with zero: no
            # positive ratio traps Q at this state.
            require(first <= 0 and second <= 0,
                    ("unclassified blend signs", first_row, second_row, state))
    require(not whole_line and large and small,
            ("degenerate endpoint bank", first_row, second_row, whole_line))
    require(all(bound > 0 for bound in (*large.values(), *small.values())),
            ("nonpositive endpoint", first_row, second_row))
    strength = max(
        upper / lower
        for upper in small.values() for lower in large.values()
    )
    relation = frozenset(
        (large_state, small_state)
        for large_state, lower in large.items()
        for small_state, upper in small.items()
        if upper / lower == strength
    )
    require(strength > 1 and relation,
            ("nonoverlapping maximizing relation", first_row, second_row))

    # Hostile sign-convention audit.  Endpoints vanish at their own bounds,
    # and a strict interior ratio traps Q at both states with negative
    # Q-minus-source delta.
    for large_state, small_state in relation:
        lower = large[large_state]
        upper = small[small_state]
        require(upper > lower and upper / lower == strength,
                ("bad maximizing overlap", first_row, second_row))
        large_first = responses[large_state][first_row - 1]
        large_second = responses[large_state][second_row - 1]
        small_first = responses[small_state][first_row - 1]
        small_second = responses[small_state][second_row - 1]
        require(lower * large_first + large_second == 0,
                ("large endpoint sign error", first_row, second_row))
        require(upper * small_first + small_second == 0,
                ("small endpoint sign error", first_row, second_row))
        interior = (lower + upper) / 2
        large_source = interior * large_first + large_second
        small_source = interior * small_first + small_second
        require(large_source > 0 and small_source > 0,
                ("overlap is not a two-state trap", first_row, second_row))
        require(-large_source < 0 and -small_source < 0,
                ("Q-directed delta sign reversed", first_row, second_row))
    return strength, relation, large, small


oriented_solutions = {}
for left, right in core_edges:
    forward = solve_ordered_edge(left, right)
    reverse = solve_ordered_edge(right, left)
    oriented_solutions[(left, right)] = forward
    oriented_solutions[(right, left)] = reverse
    forward_strength, forward_relation, forward_large, forward_small = forward
    reverse_strength, reverse_relation, reverse_large, reverse_small = reverse
    require(reverse_strength == forward_strength,
            ("strength changed under orientation", left, right))
    require(reverse_relation == frozenset(
        (small_state, large_state)
        for large_state, small_state in forward_relation
    ), ("direct reverse solve did not invert relation", left, right))
    for large_state, small_state in forward_relation:
        require(reverse_large[small_state] == 1 / forward_small[small_state]
                and reverse_small[large_state] == 1 / forward_large[large_state],
                ("reverse endpoint is not reciprocal", left, right))


def oriented_relation(first, second):
    return oriented_solutions[(first, second)][1]


strength_bank = {
    edge: oriented_solutions[edge][0] for edge in core_edges
}
relation_bank = {
    edge: oriented_solutions[edge][1] for edge in core_edges
}
relation_digest = hashlib.sha256(repr(tuple(
    (edge, strength_bank[edge], tuple(sorted(relation_bank[edge])))
    for edge in core_edges
)).encode("ascii")).hexdigest()
require(relation_digest == unique_transcript_suffix(
    CANDIDATE_OUT, "relation_bank_sha256="),
    "independent relation bank differs from frozen candidate")


def enumerate_simple_paths(edges):
    graph = {vertex: set() for vertex in core_vertices}
    for left, right in edges:
        graph[left].add(right)
        graph[right].add(left)
    answer = set()
    for start in core_vertices:
        pending = [(start, (start,))]
        while pending:
            vertex, path = pending.pop()
            for neighbor in sorted(graph[vertex]):
                if neighbor in path:
                    continue
                extension = path + (neighbor,)
                answer.add(extension)
                pending.append((neighbor, extension))
    return tuple(sorted(answer))


@lru_cache(maxsize=None)
def lifted_endpoint_pairs(path):
    require(len(path) >= 2, "lift requested for a vertex path")
    active = set(oriented_relation(path[0], path[1]))
    for first, second in zip(path[1:], path[2:]):
        following = oriented_relation(first, second)
        active = {
            (start_state, end_state)
            for start_state, middle_state in active
            for left_state, end_state in following
            if middle_state == left_state
        }
        if not active:
            break
    return frozenset(active)


def path_lifts(path):
    return bool(lifted_endpoint_pairs(tuple(path)))


full_paths = enumerate_simple_paths(core_edges)
tree_paths = enumerate_simple_paths(selected_tree)
backbone_paths = enumerate_simple_paths(backbone)
full_lifts = tuple(path for path in full_paths if path_lifts(path))
tree_lifts = tuple(path for path in tree_paths if path_lifts(path))
backbone_lifts = tuple(path for path in backbone_paths if path_lifts(path))
require((len(full_paths), len(full_lifts)) == (21226, 2442),
        "full-core lift census")
require((len(tree_paths), len(tree_lifts)) == (132, 68),
        "selected-tree lift census")
require((len(backbone_paths), len(backbone_lifts)) == (36, 36),
        "backbone lift census")


def length_histogram(paths):
    return tuple(sorted(Counter(len(path) - 1 for path in paths).items()))


require(length_histogram(backbone_lifts)
        == ((1, 14), (2, 10), (3, 6), (4, 4), (5, 2)),
        "backbone lift length histogram")
full_lift_digest = hashlib.sha256(repr(tuple(
    (path, path_lifts(path)) for path in full_paths
)).encode("ascii")).hexdigest()
require(full_lift_digest == unique_transcript_suffix(
    CANDIDATE_OUT, "full_core_lift_sha256="),
    "independent full-core lift flags differ from frozen candidate")


# Recompute the weighted phase minima using the independently solved strengths.
phase_records = {target: [] for target in range(12)}
for path in full_paths:
    cost = Fraction(1)
    for first, second in zip(path, path[1:]):
        cost *= strength_bank[tuple(sorted((first, second)))]
    target = (phase[path[-1]] - phase[path[0]]) % 12
    phase_records[target].append((cost, path))
minimum_paths = {}
for target, records in phase_records.items():
    minimum = min(cost for cost, _ in records)
    minimum_paths[target] = tuple(sorted(
        path for cost, path in records if cost == minimum
    ))
require(minimum_paths == proved_minimum_paths,
        "independent weights changed THM-3277 phase minima")
oriented_minimizers = tuple(
    path for target in range(12) for path in minimum_paths[target]
)
require(len(oriented_minimizers) == 14
        and all(path_lifts(path) for path in oriented_minimizers),
        "not all fourteen oriented phase minimizers lift")


# Enumerate global sections component-by-component.  The backbone is a forest,
# so each component is rooted and propagated through the oriented relations.
backbone_vertices = tuple(sorted({vertex for edge in backbone for vertex in edge}))
backbone_graph = {vertex: set() for vertex in backbone_vertices}
for left, right in backbone:
    backbone_graph[left].add(right)
    backbone_graph[right].add(left)

components = []
unseen = set(backbone_vertices)
while unseen:
    root = min(unseen)
    component = {root}
    queue = [root]
    for vertex in queue:
        for neighbor in sorted(backbone_graph[vertex]):
            if neighbor not in component:
                component.add(neighbor)
                queue.append(neighbor)
    unseen -= component
    components.append(tuple(sorted(component)))
components = tuple(sorted(components))
require(tuple(sorted(map(len, components))) == (3, 6),
        "backbone nonisolated component census")


def component_sections(component):
    root = min(component)
    parent = {root: None}
    order = [root]
    for vertex in order:
        for neighbor in sorted(backbone_graph[vertex]):
            if neighbor not in parent:
                parent[neighbor] = vertex
                order.append(neighbor)
    require(len(order) == len(component), "component traversal failed")
    require(sum(len(backbone_graph[v]) for v in component) // 2
            == len(component) - 1, "backbone component is not a tree")
    answer = []

    def extend(index, assignment):
        if index == len(order):
            answer.append(tuple(sorted(assignment.items())))
            return
        vertex = order[index]
        predecessor = parent[vertex]
        for state in link_states:
            if predecessor is not None and (
                    assignment[predecessor], state
                    ) not in oriented_relation(predecessor, vertex):
                continue
            assignment[vertex] = state
            extend(index + 1, assignment)
            del assignment[vertex]

    extend(0, {})
    return tuple(sorted(set(answer)))


component_section_banks = tuple(component_sections(component)
                                for component in components)
sections = []
for choices in product(*component_section_banks):
    assignment = dict(item for choice in choices for item in choice)
    sections.append(tuple(assignment[vertex] for vertex in backbone_vertices))
sections = tuple(sorted(set(sections)))
require(len(sections) == 4, "global section census")
section_words = tuple(
    tuple(edit_name(state) for state in section) for section in sections
)
expected_section_words = (
    ("Q+4", "Q-1", "Q-7", "Q-1", "Q+4", "Q+1", "Q-1", "Q-1", "Q+4"),
    ("Q+4", "Q-1", "Q-7", "Q-1", "Q+4", "Q+2", "Q-1", "Q-1", "Q+4"),
    ("Q+4", "Q-1", "Q-7", "Q-1", "Q+4", "Q+4", "Q-1", "Q-1", "Q+4"),
    ("Q+4", "Q-1", "Q-7", "Q-1", "Q+4", "Q+5", "Q-1", "Q-1", "Q+4"),
)
require(section_words == expected_section_words, "global section anatomy")
section_digest = hashlib.sha256(repr((backbone_vertices, sections)).encode(
    "ascii"
)).hexdigest()
require(section_digest == unique_transcript_suffix(
    CANDIDATE_OUT, "backbone_section_sha256="),
    "independent sections differ from frozen candidate")


# Q+4 is a maximizing endpoint on every backbone edge.  Audit its side in the
# relation, then reconstruct the primitive equal-margin blends independently.
common_state = tuple(sorted(reset + (4,)))
require(common_state in link_states and edit_name(common_state) == "Q+4",
        "common state is not Q+4")
common_roles = []
common_multiplicities = []
for edge in backbone:
    relation = relation_bank[edge]
    hits = tuple(pair for pair in sorted(relation) if common_state in pair)
    require(hits, ("Q+4 absent from backbone relation", edge))
    sides = {
        "source" if left == common_state else "target"
        for left, right in hits
    }
    require(len(sides) == 1 and all(left != right for left, right in hits),
            ("Q+4 has ambiguous relation typing", edge))
    common_roles.append((edge, next(iter(sides))))
    common_multiplicities.append((edge, len(hits)))
common_roles = tuple(common_roles)
common_multiplicities = tuple(common_multiplicities)


def primitive_equal_margin_weights(rows):
    row_values = {row: responses[common_state][row - 1] for row in rows}
    require(all(value != 0 for value in row_values.values()),
            "Q+4 has a zero selected row")
    common_multiple = lcm(*(abs(value) for value in row_values.values()))
    raw = {
        row: (2 if value > 0 else 1) * common_multiple // abs(value)
        for row, value in row_values.items()
    }
    divisor = reduce(gcd, raw.values())
    weights = {row: weight // divisor for row, weight in raw.items()}
    contributions = {
        row: weights[row] * row_values[row] for row in rows
    }
    require(reduce(gcd, weights.values()) == 1,
            "equal-margin weights are not primitive")
    return row_values, weights, contributions


backbone_row_values, backbone_weights, backbone_contributions = (
    primitive_equal_margin_weights(backbone_vertices)
)
backbone_edge_margins = tuple(
    (edge, backbone_contributions[edge[0]] + backbone_contributions[edge[1]])
    for edge in backbone
)
require(len({margin for _, margin in backbone_edge_margins}) == 1,
        "backbone edge margins differ")
backbone_margin = backbone_edge_margins[0][1]
backbone_full_margin = sum(backbone_contributions.values())
require(backbone_margin > 0 and backbone_full_margin == 3 * backbone_margin,
        "backbone equal-margin sign convention")

core_row_values, core_weights, core_contributions = (
    primitive_equal_margin_weights(core_vertices)
)
positive_rows = frozenset(row for row, value in core_row_values.items()
                          if value > 0)
negative_rows = frozenset(row for row, value in core_row_values.items()
                          if value < 0)
require(positive_rows == selector_small and negative_rows == selector_full,
        "Q+4 sign cut differs from THM-3278")
core_edge_margins = tuple(
    ((left, right), core_contributions[left] + core_contributions[right])
    for left, right in core_edges
)
require(len({margin for _, margin in core_edge_margins}) == 1,
        "full-core edge margins differ")
core_margin = core_edge_margins[0][1]
core_full_margin = sum(core_contributions.values())
require(core_margin > 0 and core_full_margin == 3 * core_margin,
        "full-core equal-margin sign convention")

# At every edge, the primitive weighted sum is the source blend (positive),
# while its Q-directed delta is the negative of that margin.  This explicitly
# rules out the common blend/delta sign reversal.
for (left, right), margin in core_edge_margins:
    ratio = Fraction(core_weights[left], core_weights[right])
    source_blend = (
        ratio * responses[common_state][left - 1]
        + responses[common_state][right - 1]
    )
    require(core_weights[right] * source_blend == margin > 0,
            ("weighted source blend sign error", left, right))
    require(-core_weights[right] * source_blend == -margin < 0,
            ("weighted Q-directed delta sign error", left, right))

candidate_rows = dict(ast.literal_eval(unique_transcript_suffix(
    CANDIDATE_OUT, "common_state_row_values="
)))
candidate_backbone_weights = dict(ast.literal_eval(unique_transcript_suffix(
    CANDIDATE_OUT, "primitive_integral_row_weights="
)))
candidate_core_weights = dict(ast.literal_eval(unique_transcript_suffix(
    CANDIDATE_OUT, "full_core_primitive_row_weights="
)))
require(candidate_rows == backbone_row_values
        and candidate_backbone_weights == backbone_weights
        and candidate_core_weights == core_weights,
        "independent Q+4 values or weights differ from frozen candidate")


# Strong chronological no-go: every relation arrow changes two count
# coordinates in L1, whereas any physical one-pole transition has L1 one.
canonical_relation_arrows = {
    pair for relation in relation_bank.values() for pair in relation
}
all_oriented_relation_arrows = {
    pair for solution in oriented_solutions.values() for pair in solution[1]
}
require(all(left != right for left, right in all_oriented_relation_arrows),
        "relation contains a diagonal state arrow")
relation_l1_histogram = tuple(sorted(Counter(
    count_l1(left, right) for left, right in all_oriented_relation_arrows
).items()))
require(relation_l1_histogram == ((2, len(all_oriented_relation_arrows)),),
        "a dominance relation arrow is a one-pole transition")
physical_q_arrows = {(state, reset) for state in link_states}
require(all(count_l1(left, right) == 1 for left, right in physical_q_arrows),
        "Q-directed physical link arrow is not one-pole")
require(not all_oriented_relation_arrows & physical_q_arrows
        and not all_oriented_relation_arrows
        & {(right, left) for left, right in physical_q_arrows},
        "static relation intersects a physical Q edge")


# Selected-tree hostile sections: the two demanded middle fibres do not meet.
root_hostile = (2, 17, 16)
first_middle_fibre = {
    right for _, right in oriented_relation(2, 17)
}
second_middle_fibre = {
    left for left, _ in oriented_relation(17, 16)
}
require(tuple(sorted(edit_name(state) for state in first_middle_fibre)) == ("Q-1",)
        and tuple(sorted(edit_name(state) for state in second_middle_fibre)) == ("Q-8",),
        "canonical-root hostile fibres drift")
require(not first_middle_fibre & second_middle_fibre
        and not path_lifts(root_hostile),
        "canonical-root hostile unexpectedly lifts")
shortest_tree_nonlifts = tuple(
    path for path in tree_paths if len(path) == 3 and not path_lifts(path)
)
require(shortest_tree_nonlifts == (
    (2, 17, 16), (11, 3, 16), (16, 3, 11), (16, 17, 2),
), "shortest selected-tree hostile census")


# Independent THM-3288 sanity check.  Build each symmetric decorated adjacency
# matrix directly from the reconstructed relations, then derive the shortest
# scalar recurrence by exact Hankel linear algebra (not Berlekamp--Massey).
recurrence_tree = ast.parse(RECURRENCE_SCOUT.read_text(encoding="utf-8"))
scout_expected = literal_assignment(recurrence_tree, "EXPECTED")


def lifted_adjacency(edges):
    arrows = set()
    nodes = set()
    for left_row, right_row in edges:
        for left_state, right_state in oriented_relation(left_row, right_row):
            left = (left_row, left_state)
            right = (right_row, right_state)
            nodes.update((left, right))
            arrows.update(((left, right), (right, left)))
    nodes = tuple(sorted(nodes))
    index = {node: position for position, node in enumerate(nodes)}
    adjacency = sp.zeros(len(nodes), len(nodes))
    for left, right in arrows:
        adjacency[index[right], index[left]] = 1
    require(adjacency == adjacency.T and all(adjacency[i, i] == 0
                                             for i in range(len(nodes))),
            "decorated adjacency is not simple symmetric")
    return nodes, tuple(sorted(arrows)), adjacency


def adjacency_sequence(adjacency, terms):
    vector = sp.ones(adjacency.rows, 1)
    answer = []
    for _ in range(terms):
        answer.append(int(sum(vector)))
        vector = adjacency * vector
    return tuple(answer)


def hankel_recurrence(sequence, maximum_order):
    for order in range(1, maximum_order + 1):
        system = sp.Matrix([
            [sequence[order + shift - lag]
             for lag in range(1, order + 1)]
            for shift in range(order)
        ])
        if system.det() == 0:
            continue
        target = sp.Matrix(sequence[order:2 * order])
        solution = system.LUsolve(target)
        coefficients = tuple(
            Fraction(int(value.p), int(value.q)) for value in solution
        )
        if all(
            sequence[position] == sum(
                coefficients[lag - 1] * sequence[position - lag]
                for lag in range(1, order + 1)
            )
            for position in range(order, len(sequence))
        ):
            return order, coefficients
    raise RuntimeError("no exact Hankel recurrence found")


def cyclic_residual(adjacency, coefficients):
    """Return q(A)1 for q(z)=z^r-sum_lag c_lag*z^(r-lag)."""
    order = len(coefficients)
    one = sp.ones(adjacency.rows, 1)
    residual = adjacency ** order * one
    for lag, coefficient in enumerate(coefficients, 1):
        residual -= int(coefficient) * adjacency ** (order - lag) * one
    return residual


recurrence_records = []
for family, edges in (
    ("backbone", backbone),
    ("selected_tree", selected_tree),
    ("full_core", core_edges),
):
    nodes, arrows, adjacency = lifted_adjacency(edges)
    sequence = adjacency_sequence(adjacency, 4 * len(nodes) + 40)
    order, coefficients = hankel_recurrence(sequence, len(nodes))
    expected = scout_expected[family]
    expected_coefficients = tuple(Fraction(value) for value in expected["coefficients"])
    require(len(nodes) == expected["vertices"] and len(arrows) == expected["arrows"],
            ("lifted automaton census differs", family))
    require(order == expected["order"] and coefficients == expected_coefficients,
            ("independent Hankel recurrence differs", family))
    require(sequence[:15] == tuple(expected["initial"]),
            ("independent adjacency sequence differs", family))
    leading_hankel = sp.Matrix(order, order,
                               lambda row, column: sequence[row + column])
    extended_hankel = sp.Matrix(order + 1, order + 1,
                                lambda row, column: sequence[row + column])
    require(leading_hankel.rank() == order
            and extended_hankel.rank() == order,
            ("Hankel rank did not stabilize", family))
    if family in {"backbone", "selected_tree"}:
        residual = cyclic_residual(adjacency, coefficients)
        require(residual == sp.zeros(len(nodes), 1),
                ("reported cyclic polynomial does not annihilate 1", family))
        annihilator_record = "q%d(A)1=0" % order
    else:
        require(order == 15 and coefficients[-1] == 0,
                "full-core prefix recurrence is not z times an even q14")
        tail_coefficients = coefficients[:-1]
        tail_residual = cyclic_residual(adjacency, tail_coefficients)
        prefix_residual = cyclic_residual(adjacency, coefficients)
        require(tail_residual != sp.zeros(len(nodes), 1),
                "full-core q14 unexpectedly annihilates the initial vector")
        require(int(sum(tail_residual)) == -1392,
                "full-core q14 scalar initial residual")
        require(adjacency * tail_residual == sp.zeros(len(nodes), 1)
                and prefix_residual == sp.zeros(len(nodes), 1),
                "full-core q14 residual does not die after one adjacency step")
        annihilator_record = (
            "1Tq14(A)1=-1392,Aq14(A)1=0,zq14(A)1=0"
        )
    recurrence_records.append((
        family, len(nodes), len(arrows), order,
        tuple(int(value) for value in coefficients), sequence[:15],
        annihilator_record,
    ))
recurrence_records = tuple(recurrence_records)


print("THM-3287 MAXIMIZING-WITNESS SECTION INDEPENDENT HOSTILE AUDIT")
print("dependency_hash_checks=%d;candidate_imports=0;recurrence_scout_imports=0" %
      len(DEPENDENCIES))
print("source_ast=(assert_nodes=%d,float_literals=%d)" %
      (assert_nodes, float_literals))
print("response_bank=(11_link_states,22_rows,all_integral);sha256=%s" %
      response_bank_digest)
print("link_edits=%s" % (tuple(sorted(map(edit_name, link_states))),))
print("orientation=ordered_u_to_v_means_lambda*f_u+f_v;trap_iff_source_blend_nonnegative")
print("endpoint_typing=(large:first_positive,small:second_positive);direct_reverse_solve=PASS;reciprocal_bounds=PASS")
print("relation_bank_sha256=%s;frozen_candidate_match=PASS" % relation_digest)
print("backbone_paths=(36/36_lift);histogram=%s" %
      (length_histogram(backbone_lifts),))
print("phase_minimizers=(14/14_lift;independent_weight_replay=PASS)")
print("global_sections=4;row_order=%s;words=%s" %
      (backbone_vertices, section_words))
print("global_section_sha256=%s;component_section_counts=%s" %
      (section_digest, tuple(map(len, component_section_banks))))
print("Q+4_maximizing_roles=%s;multiplicities=%s" %
      (common_roles, common_multiplicities))
print("Q+4_backbone_values=%s" %
      (tuple((row, backbone_row_values[row]) for row in backbone_vertices),))
print("Q+4_backbone_primitive_weights=%s" %
      (tuple((row, backbone_weights[row]) for row in backbone_vertices),))
print("Q+4_backbone_margins=(7_edges_each=%d,full_9=%d)" %
      (backbone_margin, backbone_full_margin))
print("THM3278_Q+4_sign_cut=(positive_small=%s,negative_full=%s)" %
      (tuple(sorted(positive_rows)), tuple(sorted(negative_rows))))
print("Q+4_full_core_primitive_weights=%s" %
      (tuple((row, core_weights[row]) for row in core_vertices),))
print("Q+4_full_core_margins=(22_edges_each=%d,full_12=%d);source_positive_Q_delta_negative=PASS" %
      (core_margin, core_full_margin))
print("selected_tree_paths=(68/132_lift);shortest_nonlifts=%s" %
      (shortest_tree_nonlifts,))
print("full_core_paths=(2442/21226_lift);sha256=%s" % full_lift_digest)
print("canonical_root_hostile=2-17-16;middle_fibres=(Q-1,Q-8),disjoint")
print("chronology_no_go=(canonical_unique_arrows=%d,all_oriented_unique_arrows=%d,L1_histogram=%s;physical_one_pole_L1=1;all_Q_directed_link_arrows_end_at_Q)" %
      (len(canonical_relation_arrows), len(all_oriented_relation_arrows),
       relation_l1_histogram))
for record in recurrence_records:
    family, vertices, arrows, order, coefficients, initial, annihilator = record
    print("independent_lifted_walk=%s;vertices=%d;arrows=%d;hankel_order=%d;coefficients=%s;initial_0_to_14=%s;cyclic=%s" %
          (family, vertices, arrows, order, coefficients, initial, annihilator))
print("scope=static_reset_link_dominance_relations_only;not_physical_transitions_not_response_composition_no_GMC_FC_JC_or_LRC_consequence")
print("all_exact_checks=PASS")
