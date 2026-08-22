#!/usr/bin/env python3
"""Exact scout for maximizing-witness sections on the THM-3277 backbone.

The clutch weight on an ordered row pair (u,v) comes from the overlap between
a large-ratio trap for lambda*f_u+f_v and a small-ratio trap.  This script
rebuilds those traps from THM-3238's coefficient formulas, regards the
maximizing pairs as a relation from the u-dominant witness to the v-dominant
witness, and audits gluing along simple row paths.  The relation is a static
dominance handoff, not a physical response-state transition.
"""

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from functools import reduce
from math import gcd, lcm
from pathlib import Path


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
THM3238_SCRIPT = (
    ROOT / "04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py"
)
THM3254_SCRIPT = ROOT / "04-computation/gmc_first_shell_pair_clutch_thm3254.py"
THM3277_SCRIPT = (
    ROOT / "04-computation/gmc_weighted_critical_phase_path_atlas_scout_20260803.py"
)
THM3278_SCRIPT = (
    ROOT / "04-computation/fc3_selector_origin_bipartition_phase_bridge_thm3278.py"
)
DEPENDENCIES = {
    ROOT / "01-canon/theorems/THM-3238-complete-physical-product-gamma-bank-unique-reset-stitch.md":
        "3c60ad9f7b74a10df5e6bba5b999e2dffbeb08a8b9d0886bc88a6a923dce1ca1",
    THM3238_SCRIPT:
        "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa",
    ROOT / "05-knowledge/results/gmc_complete_physical_bank_unique_reset_thm3238.out":
        "77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b",
    ROOT / "01-canon/theorems/THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go.md":
        "c7c5948d5181dc845da8f40683517c06e3eecfe837e17cb84d356b91eae1a081",
    THM3254_SCRIPT:
        "05efd37eeedeca7e3be581977a894592a7873d94a966f06d9533482cc8498fee",
    ROOT / "05-knowledge/results/gmc_first_shell_pair_clutch_thm3254.out":
        "dd415c8ce6e2e196c115421d3508addabb724305a16843509264a8b3205beee9",
    ROOT / "01-canon/theorems/THM-3269-scale-invariant-clutch-strength-and-canonical-weighted-bispanning-polarization.md":
        "55ca134eece22299ebc0c0e997f67a3247c6ee70281c34dd4f9f87cf631647d0",
    ROOT / "04-computation/gmc_scale_invariant_weighted_bispanning_thm3269.py":
        "41ae9aeb01fea1384f59f3a2687b1a0482954bf202e5be2b6fc928ef579b116a",
    ROOT / "05-knowledge/results/gmc_scale_invariant_weighted_bispanning_thm3269.out":
        "65e6e1bb04f3d42d64c2a8e5322c0f3b37c9d05ed6b053670a2cd46293742e3c",
    ROOT / "01-canon/theorems/THM-3277-weighted-critical-phase-geodesic-backbone-and-exchange-subatlas.md":
        "a10483aefd514292a5277ac6e5e426e961da55aa6f2ff2c1249bf82e679c42ec",
    THM3277_SCRIPT:
        "0beca08f9214bd6befeafdc0ccc7648be33dd8c1c2d2f12560d2e56708f2cfcf",
    ROOT / "05-knowledge/results/gmc_weighted_critical_phase_path_atlas_scout_20260803.out":
        "a55cf2bb5c130d06aaae1c84dc70d2e78e570d82253c8fdb624c4f9d86ede2e0",
    ROOT / "01-canon/theorems/THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary.md":
        "d6c0fe74d61b11951c9efa734fadbba0bec7f86a76621812b584ec1590eba89e",
    THM3278_SCRIPT:
        "07cf5cc1056bdd978a3a1d1146fee8139bef9a84e3b83aad08a0522b4df63a0f",
    ROOT / "05-knowledge/results/fc3_selector_origin_bipartition_phase_bridge_thm3278.out":
        "5e67eec14698f5258d17c5ba1ed295de3b9ff355aa279adb6d523ce3df683c7a",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected in DEPENDENCIES.items():
    require(hashlib.sha256(lf_bytes(dependency)).hexdigest() == expected,
            ("dependency hash drift", dependency.name))

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assert_nodes == 0, "optimization-sensitive assert")
require(float_literals == 0, "floating literal")


def literal_assignment(module, name):
    for node in module.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == name
                        for target in node.targets)):
            return ast.literal_eval(node.value)
    raise RuntimeError(("missing literal assignment", name))


# Load only THM-3238's pre-main exact definitions.  No trap interval, clutch
# strength, maximizing pair, path lift, or section is imported.
upstream_tree = ast.parse(THM3238_SCRIPT.read_text(encoding="utf-8"))
upstream_prefix = []
for upstream_node in upstream_tree.body:
    if isinstance(upstream_node, ast.FunctionDef) and upstream_node.name == "main":
        break
    upstream_prefix.append(upstream_node)
upstream_namespace = {
    "__file__": str(THM3238_SCRIPT),
    "__name__": "backbone_witness_exact_prefix",
}
upstream_module = ast.fix_missing_locations(
    ast.Module(body=upstream_prefix, type_ignores=[])
)
exec(compile(upstream_module, str(THM3238_SCRIPT), "exec"), upstream_namespace)

reset = upstream_namespace["RESET"]
values = upstream_namespace["VALUES"]
multiplicity = upstream_namespace["MULTIPLICITY"]
certificate = upstream_namespace["CERTIFICATE"]
upsets = upstream_namespace["UPSETS"]
bank = upstream_namespace["BANK"]
coefficient_vectors = upstream_namespace["coefficient_vectors"]
states = upstream_namespace["STATES"]
reset_counts = Counter(reset)
require(reset == (1, 3, 3, 4, 5, 6, 7, 8), "reset drift")


def physical(state):
    counts = Counter(state)
    return bool(state) and all(
        counts[value] <= multiplicity[value] for value in values
    )


def reset_distance(state):
    counts = Counter(state)
    return sum(abs(counts[value] - reset_counts[value]) for value in values)


def q_neighbors(state):
    answer = []
    counts = Counter(state)
    for value in values:
        if counts[value] == reset_counts[value]:
            continue
        changed = list(state)
        if counts[value] < reset_counts[value]:
            changed.append(value)
        else:
            changed.remove(value)
        if changed:
            target = tuple(sorted(changed))
            require(physical(target), "Q-neighbor left physical bank")
            require(reset_distance(target) + 1 == reset_distance(state),
                    "Q-neighbor lost strict reset descent")
            answer.append((value, target))
    return tuple(answer)


link_states = tuple(sorted(
    state for state in states if reset_distance(state) == 1
))
require(len(link_states) == 11, "independent reset-link census")
require(all(physical(state) for state in link_states), "nonphysical link state")
require(all(len(q_neighbors(state)) == 1
            and q_neighbors(state)[0][1] == reset
            for state in link_states),
        "reset-link edge is not the unique Q-directed transition")


def response_values(state):
    vector = coefficient_vectors(1, bank, 1, 3, state)
    exact = tuple(
        sum((vector[degree][shape] for shape in upset), Fraction(0))
        for (degree, _, _, _), upset in zip(certificate, upsets)
    )
    require(all(value.denominator == 1 for value in exact),
            ("nonintegral response", state))
    return tuple(value.numerator for value in exact)


responses = {state: response_values(state) for state in link_states}
require(response_values(reset) == (0,) * 22, "reset response not zero")

# THM-3254 supplies the proved covering-pair universe.  THM-3277 supplies the
# already-proved selected tree, phase chart, phase minima, and seven-edge
# backbone.  Every new witness relation and gluing statement is rebuilt here.
thm3254_tree = ast.parse(THM3254_SCRIPT.read_text(encoding="utf-8"))
covering_pairs = tuple(literal_assignment(thm3254_tree, "ROW_COVERING_PAIRS"))
delayed_pairs = set(literal_assignment(thm3254_tree, "DELAYED_WITNESSES"))
reset_edges = tuple(edge for edge in covering_pairs if edge not in delayed_pairs)
core_edges = tuple(edge for edge in reset_edges if 14 not in edge)
core_vertices = tuple(sorted({vertex for edge in core_edges for vertex in edge}))
require((len(covering_pairs), len(delayed_pairs), len(reset_edges),
         len(core_edges), len(core_vertices)) == (31, 8, 23, 22, 12),
        "core universe drift")
selector_small = frozenset({2, 11, 16, 18, 22})
selector_full = frozenset({3, 7, 10, 13, 17, 19, 21})
require(selector_small.isdisjoint(selector_full)
        and selector_small | selector_full == set(core_vertices)
        and all((left in selector_small) != (right in selector_small)
                for left, right in core_edges),
        "THM-3278 selector cut drift")

thm3277_tree = ast.parse(THM3277_SCRIPT.read_text(encoding="utf-8"))
selected_tree = tuple(literal_assignment(thm3277_tree, "expected_best_edges"))
phase = dict(literal_assignment(thm3277_tree, "expected_phase"))
expected_minimum_paths = dict(literal_assignment(thm3277_tree, "expected_paths"))
backbone = tuple(literal_assignment(thm3277_tree, "expected_atlas_edges"))
require(len(selected_tree) == 11 and len(backbone) == 7,
        "selected-tree/backbone drift")
require(set(backbone) <= set(selected_tree) <= set(core_edges),
        "backbone containment drift")
require(set(phase) == set(core_vertices), "phase domain drift")


def trap_interval_direct(state, first_row, second_row):
    """Ratios lambda>0 trapping lambda*f_first+f_second at state."""
    source = responses[state]
    delta_first = -source[first_row - 1]
    delta_second = -source[second_row - 1]
    if delta_first > 0:
        bound = Fraction(-delta_second, delta_first)
        return None if bound < 0 else (Fraction(0), bound)
    if delta_first < 0:
        bound = max(Fraction(0), Fraction(-delta_second, delta_first))
        return bound, None
    return None if delta_second > 0 else (Fraction(0), None)


def clutch_relation(ordered_edge):
    """Return strength, maximizing pairs, and first->second dominance relation."""
    first_row, second_row = ordered_edge
    intervals = tuple(
        (state, interval)
        for state in link_states
        if (interval := trap_interval_direct(
            state, first_row, second_row)) is not None
    )
    require(not any(lower == 0 and upper is None
                    for _, (lower, upper) in intervals),
            ("one-state whole-line trap", ordered_edge))
    small_ratio = tuple(
        (state, upper) for state, (lower, upper) in intervals
        if lower == 0 and upper is not None
    )
    large_ratio = tuple(
        (state, lower) for state, (lower, upper) in intervals
        if upper is None and lower > 0
    )
    require(small_ratio and large_ratio,
            ("missing endpoint trap fibre", ordered_edge))
    strength = max(
        upper / lower
        for _, upper in small_ratio
        for _, lower in large_ratio
    )
    maximizing_pairs = tuple(sorted(
        (small_state, large_state)
        for small_state, upper in small_ratio
        for large_state, lower in large_ratio
        if upper / lower == strength
    ))
    # lambda large makes the first row dominant; lambda small makes the
    # second row dominant.  Thus an oriented row edge first->second carries
    # the reversed state pair large-witness -> small-witness.
    dominance_relation = frozenset(
        (large_state, small_state)
        for small_state, large_state in maximizing_pairs
    )
    require(strength > 1 and dominance_relation,
            ("degenerate clutch relation", ordered_edge))
    return strength, maximizing_pairs, dominance_relation


relation_bank = {}
strength_bank = {}
maximizer_bank = {}
for edge in core_edges:
    strength, maximizing_pairs, relation = clutch_relation(edge)
    strength_bank[edge] = strength
    maximizer_bank[edge] = maximizing_pairs
    relation_bank[edge] = relation
    reversed_strength, _, reversed_relation = clutch_relation(edge[::-1])
    require(reversed_strength == strength,
            ("orientation changed strength", edge))
    require(reversed_relation == frozenset((right, left)
                                           for left, right in relation),
            ("orientation did not reverse relation", edge))


def oriented_relation(first, second):
    edge = (first, second) if first < second else (second, first)
    relation = relation_bank[edge]
    if first < second:
        return relation
    return frozenset((right, left) for left, right in relation)


relation_digest = hashlib.sha256(repr(tuple(
    (edge, strength_bank[edge], tuple(sorted(relation_bank[edge])))
    for edge in core_edges
)).encode("ascii")).hexdigest()


def enumerate_simple_paths(edges):
    graph = {vertex: [] for vertex in core_vertices}
    for left, right in edges:
        graph[left].append(right)
        graph[right].append(left)
    answer = []
    for start in core_vertices:
        stack = [(start, (start,))]
        while stack:
            vertex, path = stack.pop()
            if len(path) > 1:
                answer.append(path)
            for neighbor in graph[vertex]:
                if neighbor not in path:
                    stack.append((neighbor, path + (neighbor,)))
    return tuple(sorted(answer))


def path_lifts(path):
    possible = None
    for first, second in zip(path, path[1:]):
        relation = oriented_relation(first, second)
        if possible is None:
            possible = {right for _, right in relation}
        else:
            possible = {
                right for left, right in relation if left in possible
            }
        if not possible:
            return False
    return True


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
        "backbone path lift census")


def length_histogram(paths):
    counts = Counter(len(path) - 1 for path in paths)
    return tuple(sorted(counts.items()))


full_lift_digest = hashlib.sha256(repr(tuple(
    (path, path_lifts(path)) for path in full_paths
)).encode("ascii")).hexdigest()

# Independently replay THM-3277's exact path ranking, then test all fourteen
# minimizing oriented paths (two reversal ties at targets zero and six).
weight_lookup = {
    frozenset(edge): strength_bank[edge] for edge in core_edges
}
phase_records = {target: [] for target in range(12)}
for path in full_paths:
    weight = Fraction(1)
    for first, second in zip(path, path[1:]):
        weight *= weight_lookup[frozenset((first, second))]
    target = (phase[path[-1]] - phase[path[0]]) % 12
    phase_records[target].append((weight, path))
minimum_paths = {}
for target, records in phase_records.items():
    minimum = min(weight for weight, _ in records)
    minimum_paths[target] = tuple(
        path for weight, path in records if weight == minimum
    )
require(minimum_paths == expected_minimum_paths,
        "independent phase-minimum path replay")
oriented_minimizers = tuple(
    path for target in range(12) for path in minimum_paths[target]
)
require(len(oriented_minimizers) == 14,
        "oriented phase-minimizer tie census")
require(all(path_lifts(path) for path in oriented_minimizers),
        "phase minimizer lacks dominance lift")

# A global section assigns one link state to each of the nine nonisolated
# backbone rows, with every adjacent pair in the oriented clutch relation.
backbone_vertices = tuple(sorted({vertex for edge in backbone for vertex in edge}))
backbone_graph = {vertex: [] for vertex in backbone_vertices}
for left, right in backbone:
    backbone_graph[left].append(right)
    backbone_graph[right].append(left)
section_order = tuple(sorted(backbone_vertices,
                             key=lambda vertex: (-len(backbone_graph[vertex]),
                                                 vertex)))
sections = []


def extend_section(assignment):
    if len(assignment) == len(backbone_vertices):
        sections.append(tuple(assignment[vertex]
                              for vertex in backbone_vertices))
        return
    vertex = next(candidate for candidate in section_order
                  if candidate not in assignment)
    for state in link_states:
        if all(neighbor not in assignment
               or (state, assignment[neighbor])
               in oriented_relation(vertex, neighbor)
               for neighbor in backbone_graph[vertex]):
            assignment[vertex] = state
            extend_section(assignment)
            del assignment[vertex]


extend_section({})
sections = tuple(sorted(set(sections)))
require(len(sections) == 4, "global backbone section census")
section_digest = hashlib.sha256(repr((backbone_vertices, sections)).encode(
    "ascii")).hexdigest()


def edit_name(state):
    counts = Counter(state)
    added = tuple(value for value in values
                  for _ in range(max(0, counts[value] - reset_counts[value])))
    removed = tuple(value for value in values
                    for _ in range(max(0, reset_counts[value] - counts[value])))
    if len(added) == 1 and not removed:
        return "Q+%d" % added[0]
    if len(removed) == 1 and not added:
        return "Q-%d" % removed[0]
    raise RuntimeError(("non-link edit", state))


section_words = tuple(
    tuple(edit_name(state) for state in section) for section in sections
)
require(section_words == (
    ("Q+4", "Q-1", "Q-7", "Q-1", "Q+4", "Q+1", "Q-1", "Q-1", "Q+4"),
    ("Q+4", "Q-1", "Q-7", "Q-1", "Q+4", "Q+2", "Q-1", "Q-1", "Q+4"),
    ("Q+4", "Q-1", "Q-7", "Q-1", "Q+4", "Q+4", "Q-1", "Q-1", "Q+4"),
    ("Q+4", "Q-1", "Q-7", "Q-1", "Q+4", "Q+5", "Q-1", "Q-1", "Q+4"),
), "global section anatomy drift")

# Q+{4} is a common maximizing-pair witness on every backbone edge.  Its row
# signs also give one explicit primitive positive integral blend trapped at
# this same physical state on every edge and on all nine rows at once.
common_state = tuple(sorted(reset + (4,)))
require(common_state in link_states, "common witness left physical link")
common_roles = []
for edge in backbone:
    relation = relation_bank[edge]
    source_role = any(left == common_state for left, _ in relation)
    target_role = any(right == common_state for _, right in relation)
    require(source_role != target_role,
            ("common state role not unique", edge, source_role, target_role))
    common_roles.append((edge, "source" if source_role else "target"))
common_roles = tuple(common_roles)

backbone_row_values = {
    row: responses[common_state][row - 1] for row in backbone_vertices
}
require(all(value != 0 for value in backbone_row_values.values()),
        "zero common-state row value")
require(all(backbone_row_values[left] * backbone_row_values[right] < 0
            for left, right in backbone),
        "backbone edge does not cross common-state sign cut")
common_multiple = lcm(*(abs(value)
                        for value in backbone_row_values.values()))
raw_weights = {
    row: (2 if value > 0 else 1) * common_multiple // abs(value)
    for row, value in backbone_row_values.items()
}
weight_gcd = reduce(gcd, raw_weights.values())
integral_weights = {
    row: weight // weight_gcd for row, weight in raw_weights.items()
}
require(reduce(gcd, integral_weights.values()) == 1,
        "integral blend not primitive")
weighted_contributions = {
    row: integral_weights[row] * backbone_row_values[row]
    for row in backbone_vertices
}
edge_trap_margins = tuple(
    (edge, weighted_contributions[edge[0]] + weighted_contributions[edge[1]])
    for edge in backbone
)
full_trap_margin = sum(weighted_contributions.values())
require(len({margin for _, margin in edge_trap_margins}) == 1,
        "edge trap margins differ")
edge_margin = edge_trap_margins[0][1]
require(edge_margin > 0 and full_trap_margin == 3 * edge_margin,
        "common-state blend margin")

# The newer THM-3278 availability bipartition has an exact analytic
# representative: at Q+{4}, its small-side rows are precisely the positive
# responses and its full-side rows precisely the negative responses.  Hence
# the preceding equal-margin construction extends from the backbone to all
# 22 core edges and all twelve rows.
core_row_values = {
    row: responses[common_state][row - 1] for row in core_vertices
}
positive_rows = frozenset(row for row, value in core_row_values.items()
                          if value > 0)
negative_rows = frozenset(row for row, value in core_row_values.items()
                          if value < 0)
require(positive_rows == selector_small and negative_rows == selector_full,
        "common-state sign cut differs from selector-origin cut")
core_common_multiple = lcm(*(abs(value) for value in core_row_values.values()))
core_raw_weights = {
    row: (2 if value > 0 else 1) * core_common_multiple // abs(value)
    for row, value in core_row_values.items()
}
core_weight_gcd = reduce(gcd, core_raw_weights.values())
core_integral_weights = {
    row: weight // core_weight_gcd
    for row, weight in core_raw_weights.items()
}
core_weighted_contributions = {
    row: core_integral_weights[row] * core_row_values[row]
    for row in core_vertices
}
core_edge_margins = tuple(
    core_weighted_contributions[left] + core_weighted_contributions[right]
    for left, right in core_edges
)
core_full_margin = sum(core_weighted_contributions.values())
require(reduce(gcd, core_integral_weights.values()) == 1,
        "full-core integral blend not primitive")
require(len(set(core_edge_margins)) == 1 and core_edge_margins[0] > 0,
        "full-core equal trap margin")
core_edge_margin = core_edge_margins[0]
require(core_full_margin == 3 * core_edge_margin,
        "full twelve-row trap margin")

# Relation pairs join static maximizers only.  Every state occurring in them
# is one edit from Q and its unique Q-directed physical transition ends at Q;
# no relation arrow is one of those physical transitions.
relation_states = {
    state for relation in relation_bank.values()
    for pair in relation for state in pair
}
require(relation_states <= set(link_states),
        "relation escaped reset link")
physical_transitions = {(state, reset) for state in link_states}
relation_arrows = {
    pair for relation in relation_bank.values() for pair in relation
}
require(not (relation_arrows & physical_transitions),
        "dominance arrow mistaken for physical transition")

# The canonical-root two-edge path is a cheapest hostile extension.  Its
# required middle fibres are Q-{1} and Q-{8}, hence disjoint.  There is one
# other unoriented length-two obstruction in T_*.
root_hostile = (2, 17, 16)
first_middle_fibre = {
    right for _, right in oriented_relation(2, 17)
}
second_middle_fibre = {
    left for left, _ in oriented_relation(17, 16)
}
require(tuple(sorted(edit_name(state) for state in first_middle_fibre)) == ("Q-1",),
        "first hostile middle fibre")
require(tuple(sorted(edit_name(state) for state in second_middle_fibre)) == ("Q-8",),
        "second hostile middle fibre")
require(not (first_middle_fibre & second_middle_fibre)
        and not path_lifts(root_hostile),
        "canonical-root hostile unexpectedly lifts")
shortest_tree_nonlifts = tuple(
    path for path in tree_paths
    if len(path) == 3 and not path_lifts(path)
)
require(shortest_tree_nonlifts == (
    (2, 17, 16), (11, 3, 16), (16, 3, 11), (16, 17, 2),
), "shortest selected-tree nonlift census")


print("BACKBONE MAXIMIZING-WITNESS SECTION EXACT SCOUT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("independent_response_rebuild=(11_link_states,22_rows,all_integral)")
print("core=(V12,E22),selected_tree=11_edges,backbone=7_edges_on_9_rows")
print("dominance_relation=large_ratio_trap_to_small_ratio_trap;orientation_reversal=PASS")
print("relation_bank_sha256=%s" % relation_digest)
print("backbone_paths=(36/36_lift);length_histogram=%s" %
      (length_histogram(backbone_lifts),))
print("phase_minimizers=(14/14_lift;target_0_and_6_reversal_ties)")
print("backbone_global_sections=4;row_order=%s;words=%s" %
      (backbone_vertices, section_words))
print("backbone_section_sha256=%s" % section_digest)
print("common_maximizer_state=Q+4;edge_roles=%s" % (common_roles,))
print("common_state_row_values=%s" %
      (tuple((row, backbone_row_values[row]) for row in backbone_vertices),))
print("primitive_integral_row_weights=%s" %
      (tuple((row, integral_weights[row]) for row in backbone_vertices),))
print("common_state_trap_margins=(each_edge=%d,full_9_row=%d)" %
      (edge_margin, full_trap_margin))
print("THM3278_cut_at_Q+4=(positive_small=(2,11,16,18,22),negative_full=(3,7,10,13,17,19,21))")
print("full_core_primitive_row_weights=%s" %
      (tuple((row, core_integral_weights[row]) for row in core_vertices),))
print("full_core_Q+4_trap_margins=(22_edges_each=%d,full_12_row=%d)" %
      (core_edge_margin, core_full_margin))
print("selected_tree_paths=(68/132_lift);lift_histogram=%s;nonlift_histogram=%s" %
      (length_histogram(tree_lifts),
       length_histogram(path for path in tree_paths if not path_lifts(path))))
print("full_core_paths=(2442/21226_lift);lift_histogram=%s" %
      (length_histogram(full_lifts),))
print("full_core_lift_sha256=%s" % full_lift_digest)
print("canonical_root_hostile=2-17-16;middle_fibres=(Q-1,Q-8),disjoint")
print("shortest_Tstar_nonlifts=%s" % (shortest_tree_nonlifts,))
print("physical_Q_directed_transitions=11_all_end_at_Q;dominance_arrows_are_not_transitions")
print("map=oriented_row_path_to_static_maximizing_trap-state_section")
print("preserved=edge_clutch_maximum_and_shared_middle_witness;lost=ratio_endpoint_and_chronological_response_state")
print("scope=support_(1,3)_bank_I2_reset_link_dominance_Cech_handoff_only,no_response_composition_no_GMC_or_LRC_consequence")
print("all_exact_checks=PASS")
