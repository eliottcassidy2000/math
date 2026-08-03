#!/usr/bin/env python3
"""Independent exact audit of the multi-edge visible-kernel scout.

Rebuild the decorated graphs directly from the older pinned THM-3287 relation
source.  Do not import or execute the audited multi-edge primary.  Visibility
is recovered from the cyclic minimal polynomial of A^2, rather than from a
nullspace Gram projection.  The bordered update is checked with a Dyson
coefficient recursion, rather than the primary's local inverse recursion.

The scope is finite static symmetric relation-walk graphs only.  There is no
tournament, chronology, response-composition, FC/GMC, LRC, asymptotic, or
bit-complexity conclusion.
"""

from __future__ import annotations

import ast
import hashlib
import io
import runpy
from collections import Counter, defaultdict
from contextlib import redirect_stdout
from itertools import combinations
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
RELATION_SOURCE = ROOT / (
    "04-computation/gmc_backbone_maximizing_witness_section_scout_20260803.py"
)
RELATION_OUTPUT = ROOT / (
    "05-knowledge/results/gmc_backbone_maximizing_witness_section_scout_20260803.out"
)
HALF_SOURCE = ROOT / (
    "04-computation/gmc_selector_paired_half_transfer_partial_scout_20260803.py"
)
HALF_OUTPUT = ROOT / (
    "05-knowledge/results/gmc_selector_paired_half_transfer_partial_scout_20260803.out"
)
TARGET_SOURCE = ROOT / (
    "04-computation/"
    "gmc_multi_edge_visible_kernel_helly_bordered_resolvent_scout_20260803.py"
)
TARGET_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "gmc_multi_edge_visible_kernel_helly_bordered_resolvent_scout_20260803.out"
)
TARGET_REFLECTION = ROOT / (
    "07-reflections/"
    "multi-edge-visible-kernel-helly-and-bordered-resolvent-compiler-codex-20260803.md"
)
DEPENDENCIES = {
    RELATION_SOURCE:
        "aed89d67ec7acabfe5b4feae4a83f7c57b78053928be44a6cc319d81fa4a9cc6",
    RELATION_OUTPUT:
        "89200bc6cff7284dd33f636352f4f7f56294d90bcd2902fa96092d9f967f5fe3",
    HALF_SOURCE:
        "da47f9db80558343bd6442196c65150417241acbed65385b876ca09a66f0490f",
    HALF_OUTPUT:
        "ef276b6a99a3ec1ac59d5cc352a3aa9a50fa64440a6615695d684a68a6683208",
    TARGET_SOURCE:
        "c2e59a49e036f90492f353e05bf5304d473f491df9a6d79dcc124b0ec1273a1c",
    TARGET_OUTPUT:
        "d0662551b16e1cab7b64fe4a9b733fbc10ff9373264e7b560f0a07a6e2dffbe0",
    TARGET_REFLECTION:
        "09a382b20ef07c6440118f13e4f6f408790b5dcc42699f1dc7a024c129d1379d",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected_hash in DEPENDENCIES.items():
    require(
        hashlib.sha256(lf_bytes(dependency)).hexdigest() == expected_hash,
        ("dependency or audit-target hash drift", dependency.name),
    )

# Execute only the older relation source.  The half-transfer source and the
# audited primary are hash targets, never executable dependencies of this audit.
captured = io.StringIO()
with redirect_stdout(captured):
    relation = runpy.run_path(str(RELATION_SOURCE))
require(
    captured.getvalue().encode("utf-8") == lf_bytes(RELATION_OUTPUT),
    "older THM-3287 relation replay differs from its stored output",
)

selected_tree = tuple(relation["selected_tree"])
core_edges = tuple(relation["core_edges"])
selected_set = frozenset(selected_tree)
extra_edges = tuple(edge for edge in core_edges if edge not in selected_set)
require(
    len(selected_tree) == len(extra_edges) == 11 and len(core_edges) == 22,
    "selected-tree/full-core edge census drift",
)


def row_relation(row_edge):
    """Return active decorated nodes and undirected links for one row edge."""
    left_row, right_row = row_edge
    nodes = set()
    links = set()
    for left_state, right_state in relation["oriented_relation"](
        left_row, right_row
    ):
        left = (left_row, left_state)
        right = (right_row, right_state)
        nodes.update((left, right))
        links.add(tuple(sorted((left, right))))
    require(len(nodes) > 0 and len(links) > 0, ("empty row relation", row_edge))
    return frozenset(nodes), frozenset(links)


ROW_RELATIONS = {edge: row_relation(edge) for edge in core_edges}


def rebuilt_graph(row_edges):
    """Independent decorated adjacency reconstruction from row relations."""
    nodes = set()
    links = set()
    for row_edge in row_edges:
        edge_nodes, edge_links = ROW_RELATIONS[row_edge]
        nodes.update(edge_nodes)
        links.update(edge_links)
    nodes = tuple(sorted(nodes))
    index = {node: position for position, node in enumerate(nodes)}
    adjacency = sp.zeros(len(nodes), len(nodes))
    for left, right in links:
        adjacency[index[left], index[right]] = 1
        adjacency[index[right], index[left]] = 1
    require(
        adjacency == adjacency.T
        and sum(int(value) for value in adjacency) == 2 * len(links),
        "independent decorated adjacency reconstruction failed",
    )
    return nodes, tuple(sorted(links)), adjacency


def cyclic_zero_projection(adjacency):
    """Recover P_ker(A)1 through the minimal polynomial of A^2 on 1.

    If q is the first Krylov relation for T=A^2 and q(0) is nonzero, then
    1 lies in im(T), so its kernel projection vanishes.  If q=t*r, the
    diagonalizability of symmetric T gives P_ker(A)1=r(T)1/r(0).
    """
    transfer = adjacency * adjacency
    observer = sp.ones(adjacency.rows, 1)
    columns = []
    vector = observer
    for order in range(adjacency.rows + 1):
        if columns:
            basis = sp.Matrix.hstack(*columns)
            augmented = basis.row_join(vector)
            if DomainMatrix.from_Matrix(augmented).rank() == order:
                coefficients = basis.gauss_jordan_solve(vector)[0]
                polynomial = tuple(
                    [-coefficients[index] for index in range(order)]
                    + [sp.Integer(1)]
                )
                relation_vector = sum(
                    (
                        polynomial[index] * columns[index]
                        for index in range(order)
                    ),
                    vector,
                )
                require(
                    relation_vector == sp.zeros(adjacency.rows, 1),
                    "cyclic relation reconstruction failed",
                )
                if polynomial[0] != 0:
                    preimage = -sum(
                        (
                            polynomial[index] * columns[index - 1]
                            for index in range(1, order + 1)
                        ),
                        sp.zeros(adjacency.rows, 1),
                    ) / polynomial[0]
                    require(
                        transfer * preimage == observer,
                        "nonvisible observer lacks the certified T-preimage",
                    )
                    return sp.zeros(adjacency.rows, 1), sp.Integer(0), order
                require(polynomial[1] != 0, "cyclic zero root is not simple")
                projection = sum(
                    (
                        polynomial[index + 1] * columns[index]
                        for index in range(order)
                    ),
                    sp.zeros(adjacency.rows, 1),
                ) / polynomial[1]
                mass = sp.factor(sum(projection))
                require(
                    adjacency * projection == sp.zeros(adjacency.rows, 1)
                    and sp.factor((projection.T * projection)[0]) == mass,
                    "cyclic zero projection identities failed",
                )
                return projection, mass, order
        columns.append(vector)
        vector = transfer * vector
    raise RuntimeError("cyclic minimal polynomial failed to close")


def edge_set(mask):
    return tuple(
        extra_edges[index]
        for index in range(len(extra_edges))
        if mask & (1 << index)
    )


# Full Boolean lattice, with kernel dimension computed independently over QQ
# by DomainMatrix rank and zero mass computed by the cyclic A^2 certificate.
lattice_records = []
mass_by_mask = {}
for mask in range(1 << len(extra_edges)):
    chosen = edge_set(mask)
    nodes, links, adjacency = rebuilt_graph(selected_tree + chosen)
    _projection, mass, _cyclic_order = cyclic_zero_projection(adjacency)
    rank = DomainMatrix.from_Matrix(adjacency).rank()
    record = (
        mask, mask.bit_count(), len(nodes), 2 * len(links),
        len(nodes) - rank, mass,
    )
    lattice_records.append(record)
    mass_by_mask[mask] = mass

lattice_serialized = "".join(
    ";".join(str(value) for value in record) + "\n"
    for record in lattice_records
).encode("utf-8")
LATTICE_DIGEST = hashlib.sha256(lattice_serialized).hexdigest()
require(
    len(lattice_serialized) == 37894
    and LATTICE_DIGEST
        == "53f740b8fd3e8b1baa88a460e210662d51b6aa69ff8a157fcdfc9fd5d8510ad0",
    "independent Boolean-lattice digest differs from the audited primary",
)

records_by_size = defaultdict(list)
for record in lattice_records:
    records_by_size[record[1]].append(record)
size_summary = {}
for size, records in records_by_size.items():
    nonzero = tuple(record[-1] for record in records if record[-1] != 0)
    size_summary[size] = (
        len(records), len(records) - len(nonzero), len(nonzero),
        len(frozenset(nonzero)),
    )
EXPECTED_SIZE_SUMMARY = {
    0: (1, 1, 0, 0),
    1: (11, 10, 1, 1),
    2: (55, 45, 10, 1),
    3: (165, 115, 50, 3),
    4: (330, 188, 142, 10),
    5: (462, 212, 250, 20),
    6: (462, 169, 293, 25),
    7: (330, 91, 239, 23),
    8: (165, 29, 136, 18),
    9: (55, 4, 51, 13),
    10: (11, 0, 11, 7),
    11: (1, 0, 1, 1),
}
require(size_summary == EXPECTED_SIZE_SUMMARY, "lattice size summary drift")

pair_masks = tuple(
    mask for mask in range(1 << len(extra_edges)) if mask.bit_count() == 2
)
visible_pairs = tuple(mask for mask in pair_masks if mass_by_mask[mask] != 0)
special_bit = 1 << extra_edges.index((7, 18))
require(
    len(visible_pairs) == 10
    and all(mask & special_bit for mask in visible_pairs),
    "visible pairs are not exactly inherited from the special singleton",
)

minimal_visible_masks = []
for mask in range(1, 1 << len(extra_edges)):
    if mass_by_mask[mask] == 0:
        continue
    submask = (mask - 1) & mask
    visible_proper_subset = False
    while submask:
        if mass_by_mask[submask] != 0:
            visible_proper_subset = True
            break
        submask = (submask - 1) & mask
    if not visible_proper_subset:
        minimal_visible_masks.append(mask)

minimal_visible_sets = tuple(edge_set(mask) for mask in minimal_visible_masks)
EXPECTED_SUPPORTS = {
    ((7, 18),): (
        ((18, "Q+1"), -3), ((18, "Q+2"), 1),
        ((18, "Q+4"), 1), ((18, "Q+5"), 1),
        ((22, "Q+1"), 1), ((22, "Q+2"), 1),
        ((22, "Q+5"), 1),
    ),
    ((2, 10), (13, 18), (17, 22)): (
        ((2, "Q+4"), -12),
        ((11, "Q+1"), 4), ((11, "Q+2"), 4), ((11, "Q+5"), 4),
        ((18, "Q+1"), 3), ((18, "Q+2"), 3),
        ((18, "Q+4"), 3), ((18, "Q+5"), 3),
        ((22, "Q+1"), -4), ((22, "Q+2"), -4),
        ((22, "Q+4"), 12), ((22, "Q+5"), -4),
    ),
    ((2, 10), (13, 22), (17, 22)): (
        ((2, "Q+4"), -12),
        ((11, "Q+1"), 4), ((11, "Q+2"), 4), ((11, "Q+5"), 4),
        ((18, "Q+1"), 3), ((18, "Q+2"), 3),
        ((18, "Q+4"), 3), ((18, "Q+5"), 3),
        ((22, "Q+1"), -4), ((22, "Q+2"), -4),
        ((22, "Q+4"), 12), ((22, "Q+5"), -4),
    ),
    ((10, 22), (13, 18), (17, 22)): (
        ((2, "Q+4"), 12),
        ((11, "Q+1"), 4), ((11, "Q+2"), 4), ((11, "Q+5"), 4),
        ((18, "Q+1"), -3), ((18, "Q+2"), -3),
        ((18, "Q+4"), -3), ((18, "Q+5"), -3),
        ((22, "Q+1"), 4), ((22, "Q+2"), 4),
        ((22, "Q+4"), -12), ((22, "Q+5"), 4),
    ),
    ((10, 22), (13, 22), (17, 22)): (
        ((2, "Q+4"), 12),
        ((11, "Q+1"), 4), ((11, "Q+2"), 4), ((11, "Q+5"), 4),
        ((18, "Q+1"), -3), ((18, "Q+2"), -3),
        ((18, "Q+4"), -3), ((18, "Q+5"), -3),
        ((22, "Q+1"), 4), ((22, "Q+2"), 4),
        ((22, "Q+4"), -12), ((22, "Q+5"), 4),
    ),
    ((13, 22), (17, 22), (19, 22)): (
        ((2, "Q+4"), 3),
        ((22, "Q+1"), 1), ((22, "Q+2"), 1),
        ((22, "Q+4"), -3), ((22, "Q+5"), 1),
    ),
}
require(
    len(minimal_visible_sets) == 6
    and frozenset(minimal_visible_sets) == frozenset(EXPECTED_SUPPORTS),
    "minimal visible edge-set census drift",
)


def primitive(vector):
    denominators = tuple(sp.denom(value) for value in vector)
    scale = sp.ilcm(*denominators)
    integers = tuple(int(scale * value) for value in vector)
    divisor = sp.igcd(*tuple(abs(value) for value in integers if value))
    result = sp.Matrix(tuple(value // divisor for value in integers))
    if sum(result) < 0:
        result = -result
    return result


def labelled_support(nodes, vector):
    return tuple(
        ((row, relation["edit_name"](state)), int(value))
        for (row, state), value in zip(nodes, vector)
        if value != 0
    )


circuit_records = []
for chosen in EXPECTED_SUPPORTS:
    nodes, _links, adjacency = rebuilt_graph(selected_tree + chosen)
    projection, mass, _order = cyclic_zero_projection(adjacency)
    direction = primitive(projection)
    support = labelled_support(nodes, direction)
    require(support == EXPECTED_SUPPORTS[chosen], ("support drift", chosen))
    circuit_records.append((chosen, mass, support, direction, nodes))

EXPECTED_CIRCUIT_MASSES = {
    ((7, 18),): sp.Rational(3, 5),
    ((2, 10), (13, 18), (17, 22)): sp.Rational(12, 35),
    ((2, 10), (13, 22), (17, 22)): sp.Rational(12, 35),
    ((10, 22), (13, 18), (17, 22)): sp.Rational(12, 35),
    ((10, 22), (13, 22), (17, 22)): sp.Rational(12, 35),
    ((13, 22), (17, 22), (19, 22)): sp.Rational(3, 7),
}
require(
    {record[0]: record[1] for record in circuit_records}
        == EXPECTED_CIRCUIT_MASSES,
    "minimal circuit mass drift",
)

# Fixed 26-node universe for persistence and local-update checks.
full_nodes, _full_links, _full_adjacency_local = rebuilt_graph(core_edges)
full_index = {node: index for index, node in enumerate(full_nodes)}


def embedded_graph(row_edges):
    nodes, links, _adjacency = rebuilt_graph(row_edges)
    adjacency = sp.zeros(len(full_nodes), len(full_nodes))
    for left, right in links:
        adjacency[full_index[left], full_index[right]] = 1
        adjacency[full_index[right], full_index[left]] = 1
    return frozenset(nodes), adjacency


base_active, base_adjacency = embedded_graph(selected_tree)
edge_deltas = {}
for edge in extra_edges:
    _active, adjacency = embedded_graph(selected_tree + (edge,))
    edge_deltas[edge] = adjacency - base_adjacency

special_record = next(record for record in circuit_records if record[0] == ((7, 18),))
special_direction = sp.zeros(len(full_nodes), 1)
for node, value in zip(special_record[4], special_record[3]):
    special_direction[full_index[node]] = value
special_compatible = tuple(
    edge for edge in extra_edges
    if edge != (7, 18)
    and edge_deltas[edge] * special_direction
        == sp.zeros(len(full_nodes), 1)
)
require(
    len(special_compatible) == 10
    and int(sum(special_direction)) == 3
    and int((special_direction.T * special_direction)[0]) == 15,
    "special persistent direction failed",
)
special_interval_masses = tuple(
    mass_by_mask[mask]
    for mask in range(1 << len(extra_edges)) if mask & special_bit
)
require(
    len(special_interval_masses) == 1024
    and min(special_interval_masses) == sp.Rational(3, 5),
    "persistent 3/5 upper-interval lower bound failed",
)

fragile = ((2, 10), (13, 18), (17, 22))
fragile_record = next(record for record in circuit_records if record[0] == fragile)
fragile_direction = sp.zeros(len(full_nodes), 1)
for node, value in zip(fragile_record[4], fragile_record[3]):
    fragile_direction[full_index[node]] = value
fragile_compatible = tuple(
    edge for edge in extra_edges
    if edge not in fragile
    and edge_deltas[edge] * fragile_direction
        == sp.zeros(len(full_nodes), 1)
)
require(
    fragile_compatible == ((10, 16), (16, 21)),
    "fragile triple compatibility set drift",
)
fragile_mask = sum(1 << extra_edges.index(edge) for edge in fragile)
fragile_fourth_masses = {
    edge: mass_by_mask[fragile_mask | (1 << extra_edges.index(edge))]
    for edge in extra_edges if edge not in fragile
}
EXPECTED_FRAGILE_FOURTH_MASSES = {
    (3, 22): sp.Rational(12, 59),
    (7, 18): sp.Rational(14, 15),
    (10, 16): sp.Rational(12, 35),
    (10, 22): sp.Integer(0),
    (11, 13): sp.Rational(12, 37),
    (13, 22): sp.Integer(0),
    (16, 21): sp.Rational(12, 35),
    (19, 22): sp.Integer(0),
}
require(
    fragile_fourth_masses == EXPECTED_FRAGILE_FOURTH_MASSES,
    "fragile fourth-edge mass profile drift",
)


def sequence(transfer, observer, terms):
    values = []
    vector = observer
    for _degree in range(terms):
        values.append(int((observer.T * vector)[0]))
        vector = transfer * vector
    return tuple(values)


# Independent bordered audit.  Instead of expanding (C-xK)^(-1) directly,
# solve h=s+x*K*C*h coefficientwise and use correction=x*s^T*C*h.
COMPILER_TERMS = 24
base_indicator = sp.Matrix([
    1 if node in base_active else 0 for node in full_nodes
])
base_sequence = sequence(base_adjacency, base_indicator, COMPILER_TERMS)
update_shape_counts = Counter()
compiler_records = []
for first, second in combinations(extra_edges, 2):
    active, adjacency = embedded_graph(selected_tree + (first, second))
    delta = adjacency - base_adjacency
    updates = tuple(
        (row, column)
        for row in range(len(full_nodes))
        for column in range(row + 1, len(full_nodes))
        if delta[row, column] != 0
    )
    require(
        delta == delta.T
        and all(delta[row, column] == 1 for row, column in updates),
        ("pair update is not a symmetric pure addition", first, second),
    )
    local_dimension = 2 * len(updates)
    births = len(active - base_active)
    update_shape_counts[(len(updates), local_dimension, births)] += 1

    endpoints = sp.zeros(len(full_nodes), local_dimension)
    swap = sp.zeros(local_dimension, local_dimension)
    for update_index, (left, right) in enumerate(updates):
        left_column = 2 * update_index
        right_column = left_column + 1
        endpoints[left, left_column] = 1
        endpoints[right, right_column] = 1
        swap[left_column, right_column] = 1
        swap[right_column, left_column] = 1
    require(
        endpoints * swap * endpoints.T == delta
        and swap * swap == sp.eye(local_dimension),
        ("local border factorization failed", first, second),
    )

    observer = sp.Matrix([1 if node in active else 0 for node in full_nodes])
    s_coefficients = []
    k_coefficients = []
    observer_power = observer
    endpoint_power = endpoints
    for _degree in range(COMPILER_TERMS):
        s_coefficients.append(endpoints.T * observer_power)
        k_coefficients.append(endpoints.T * endpoint_power)
        observer_power = base_adjacency * observer_power
        endpoint_power = base_adjacency * endpoint_power

    h_coefficients = []
    for degree in range(COMPILER_TERMS):
        coefficient = s_coefficients[degree]
        for left_degree in range(degree):
            coefficient += (
                k_coefficients[left_degree]
                * swap
                * h_coefficients[degree - 1 - left_degree]
            )
        h_coefficients.append(coefficient)

    corrections = [sp.Integer(0)] * COMPILER_TERMS
    for degree in range(1, COMPILER_TERMS):
        corrections[degree] = sp.expand(sum(
            (
                s_coefficients[left_degree].T
                * swap
                * h_coefficients[degree - 1 - left_degree]
            )[0]
            for left_degree in range(degree)
        ))
    compiled = tuple(
        int(
            base_sequence[degree]
            + (births if degree == 0 else 0)
            + corrections[degree]
        )
        for degree in range(COMPILER_TERMS)
    )
    direct = sequence(adjacency, observer, COMPILER_TERMS)
    require(compiled == direct, ("bordered Dyson mismatch", first, second))
    compiler_records.append((
        first, second, len(updates), local_dimension, births, direct,
    ))

EXPECTED_UPDATE_SHAPES = Counter({
    (2, 4, 0): 21,
    (2, 4, 1): 21,
    (2, 4, 2): 3,
    (5, 10, 0): 7,
    (5, 10, 1): 3,
})
require(
    update_shape_counts == EXPECTED_UPDATE_SHAPES
    and len(compiler_records) == 55,
    "bordered update shape census drift",
)
compiler_serialized = "".join(
    f"{first};{second};{updates};{dimension};{births};"
    + ",".join(str(value) for value in direct) + "\n"
    for first, second, updates, dimension, births, direct in compiler_records
).encode("utf-8")
compiler_digest = hashlib.sha256(compiler_serialized).hexdigest()
require(
    len(compiler_serialized) == 11560
    and compiler_digest
        == "d969c893735eca95eb16fec7c574bfdb9e24a68160156d7a8f2eaf12a3050b8d",
    "24-term bordered sequence census digest drift",
)

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assert_nodes == 0 and float_literals == 0, "source AST gate")

print("MULTI-EDGE VISIBLE-KERNEL AND BORDERED-RESOLVENT INDEPENDENT AUDIT")
print("status=VERIFIED_FINITE_EXACT;no_theorem_id;static_symmetric_relation_walks_only")
print("dependency_and_target_hashes=BEGIN")
for dependency, expected_hash in DEPENDENCIES.items():
    print(f"{expected_hash}  {dependency.relative_to(ROOT)}")
print("dependency_and_target_hashes=END")
print("independence=executed_only_older_relation_source;target_primary_not_imported_or_executed;zero_mass_via_A2_cyclic_polynomial_not_nullspace_Gram")
print(
    f"edge_lattice=2048;serialized_bytes={len(lattice_serialized)};"
    f"sha256={LATTICE_DIGEST};exact_match=true"
)
print(f"size_summary={tuple(sorted(size_summary.items()))}")
print("pair_census=55;visible_pairs=10;new_pair_minimal_circuits=0")
print("minimal_visible_circuits=6;arity_profile=(one:1,two:0,three:5)")
for chosen, mass, support, _direction, _nodes in circuit_records:
    print(f"circuit={chosen};mass={mass};primitive_support={support}")
print(
    "persistent_circuit=((7,18),);delta_annihilation_edges=10;"
    "upper_interval_subsets=1024;minimum_zero_mass=3/5"
)
print(
    f"fragile_triple={fragile};compatible_next_edges={fragile_compatible};"
    f"fourth_edge_masses={fragile_fourth_masses}"
)
print(
    f"bordered_dyson_audit=pairs:55;terms:{COMPILER_TERMS};"
    f"local_dimension_max:10;update_shape_counts={tuple(sorted(update_shape_counts.items()))};"
    f"serialized_bytes:{len(compiler_serialized)};sha256:{compiler_digest};"
    "exact_match:true"
)
print("mechanism=cyclic_A2_zero_projector_and_local_bordered_Dyson_recursion")
print("preserved=exact_static_all_ones_walk_counts_active_births_kernel_mass_and_minimal_edge_support")
print("lost=walk_identity_chronology_tournament_orientation_response_state_and_asymptotic_bit_complexity")
print("scope=no_tournament_no_chronology_no_response_composition_no_GMC_FC_or_LRC_consequence")
print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")
print("ALL INDEPENDENT EXACT CHECKS PASSED")
