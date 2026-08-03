#!/usr/bin/env python3
"""Finite-exact selector-paired half-transfer scout for THM-3287/3288.

The active vertices are the static maximizing-witness pairs (row,state) of
THM-3288.  This companion independently rebuilds their three decorated graphs
from the pinned THM-3287 relation artifact, then studies the exact selector
split supplied by THM-3278/3287.  It proves only finite matrix identities for
static relation walks.  It does not construct chronological response dynamics.
"""

from __future__ import annotations

import ast
import hashlib
import io
import runpy
from collections import Counter, defaultdict
from contextlib import redirect_stdout
from pathlib import Path

import sympy as sp


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
SERIES_SOURCE = ROOT / (
    "04-computation/gmc_maximizing_witness_lifted_walk_rational_series_thm3288.py"
)
SERIES_OUTPUT = ROOT / (
    "05-knowledge/results/gmc_maximizing_witness_lifted_walk_rational_series_thm3288.out"
)
DEPENDENCIES = {
    RELATION_SOURCE:
        "aed89d67ec7acabfe5b4feae4a83f7c57b78053928be44a6cc319d81fa4a9cc6",
    RELATION_OUTPUT:
        "89200bc6cff7284dd33f636352f4f7f56294d90bcd2902fa96092d9f967f5fe3",
    SERIES_SOURCE:
        "f442ff0dfd999fc314f420aacfc7102ba2ba17e033d3fd13e91407727af8a8da",
    SERIES_OUTPUT:
        "44a85f317bf9cb35c437e468633a88014de268cc3fb452836960655c4b764951",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected_hash in DEPENDENCIES.items():
    require(
        hashlib.sha256(lf_bytes(dependency)).hexdigest() == expected_hash,
        ("dependency hash drift", dependency.name),
    )

captured = io.StringIO()
with redirect_stdout(captured):
    relation_source = runpy.run_path(str(RELATION_SOURCE))
require(
    captured.getvalue().encode("utf-8") == lf_bytes(RELATION_OUTPUT),
    "THM-3287 relation replay differs from stored output",
)

selector_small = frozenset(relation_source["selector_small"])
selector_full = frozenset(relation_source["selector_full"])
require(
    selector_small == frozenset({2, 11, 16, 18, 22})
    and selector_full == frozenset({3, 7, 10, 13, 17, 19, 21}),
    "selector partition drift",
)


EXPECTED = {
    "backbone": {
        "nodes": 18,
        "arrows": 32,
        "order": 10,
        "initial": (18, 32, 84, 160, 412, 798, 2060, 4012, 10426,
                    20362, 53286, 104272, 274628, 538258, 1425644),
        "half": (5, 5),
        "equitable_cells": 11,
    },
    "selected_tree": {
        "nodes": 23,
        "arrows": 40,
        "order": 14,
        "initial": (23, 40, 104, 212, 536, 1174, 2950, 6700, 16850,
                    39058, 98502, 231232, 584846, 1384502, 3510318),
        "half": (7, 7),
        "equitable_cells": 14,
    },
    "full_core": {
        "nodes": 26,
        "arrows": 68,
        "order": 15,
        "initial": (26, 68, 274, 1018, 4134, 16400, 67270, 270718,
                    1114054, 4498562, 18529428, 74885666, 308530344,
                    1247186562, 5138801640),
        "half": (8, 7),
        "equitable_cells": 16,
    },
}


def state_name(state):
    return relation_source["edit_name"](state)


def node_label(node):
    row, state = node
    return row, state_name(state)


def rebuilt_graph(edges):
    nodes = set()
    arrows = set()
    for left_row, right_row in edges:
        relation = relation_source["oriented_relation"](left_row, right_row)
        for left_state, right_state in relation:
            left = (left_row, left_state)
            right = (right_row, right_state)
            nodes.update((left, right))
            arrows.update(((left, right), (right, left)))
    nodes = tuple(sorted(nodes))
    index = {node: position for position, node in enumerate(nodes)}
    adjacency = sp.zeros(len(nodes), len(nodes))
    for left, right in arrows:
        adjacency[index[left], index[right]] = 1
    require(adjacency == adjacency.T, "decorated adjacency is not symmetric")
    require(
        sum(int(value) for value in adjacency) == len(arrows),
        "decorated arrow multiplicity drift",
    )
    return nodes, tuple(sorted(arrows)), adjacency


def sequence(transfer, start, output, terms):
    values = []
    vector = start
    for _ in range(terms):
        values.append(int((output * vector)[0]))
        vector = transfer * vector
    return tuple(values)


def krylov_relation(transfer, start):
    columns = []
    vector = start
    for order in range(transfer.rows + 2):
        if columns:
            basis = sp.Matrix.hstack(*columns)
            if basis.row_join(vector).rank() == basis.rank():
                coefficients = basis.gauss_jordan_solve(vector)[0]
                require(
                    basis * coefficients == vector,
                    "Krylov relation solve failed",
                )
                return order, tuple(coefficients), basis
        columns.append(vector)
        vector = transfer * vector
    raise RuntimeError("Krylov relation did not close")


def relation_polynomial(order, coefficients, symbol):
    return sp.Poly(
        sp.expand(
            symbol**order
            - sum(coefficients[index] * symbol**index
                  for index in range(order))
        ),
        symbol,
        domain=sp.QQ,
    )


def companion_data(order, coefficients, output, basis):
    companion = sp.zeros(order, order)
    for column in range(order - 1):
        companion[column + 1, column] = 1
    for row, coefficient in enumerate(coefficients):
        companion[row, order - 1] = coefficient
    boundary = output * basis
    require(boundary.cols == order, "observer boundary width drift")
    digest_payload = (
        tuple(tuple(int(value) for value in companion.row(row))
              for row in range(order)),
        tuple(int(value) for value in boundary),
    )
    digest = hashlib.sha256(repr(digest_payload).encode("ascii")).hexdigest()
    return companion, boundary, digest


def stable_equitable_refinement(adjacency, partition):
    partition = [tuple(sorted(cell)) for cell in partition]
    while True:
        owner = {
            vertex: cell_index
            for cell_index, cell in enumerate(partition)
            for vertex in cell
        }
        groups = defaultdict(list)
        for vertex in range(adjacency.rows):
            signature = (
                owner[vertex],
                tuple(
                    sum(int(adjacency[vertex, neighbor]) for neighbor in cell)
                    for cell in partition
                ),
            )
            groups[signature].append(vertex)
        refined = [tuple(vertices) for _, vertices in sorted(groups.items())]
        if refined == partition:
            return tuple(refined)
        partition = refined


def equitable_quotient(adjacency, partition):
    quotient = sp.zeros(len(partition), len(partition))
    for source_cell, vertices in enumerate(partition):
        for target_cell, targets in enumerate(partition):
            counts = tuple(
                sum(int(adjacency[vertex, target]) for target in targets)
                for vertex in vertices
            )
            require(len(set(counts)) == 1, "partition is not equitable")
            quotient[source_cell, target_cell] = counts[0]
    indicator = sp.zeros(adjacency.rows, len(partition))
    for cell_index, vertices in enumerate(partition):
        for vertex in vertices:
            indicator[vertex, cell_index] = 1
    require(
        adjacency * indicator == indicator * quotient,
        "equitable intertwining identity failed",
    )
    return quotient, indicator


def matrix_polynomial_on_vector(poly, matrix, vector):
    result = sp.zeros(matrix.rows, 1)
    for coefficient in poly.all_coeffs():
        result = matrix * result + coefficient * vector
    return result


z, t = sp.symbols("z t")
family_edges = {
    "backbone": relation_source["backbone"],
    "selected_tree": relation_source["selected_tree"],
    "full_core": relation_source["core_edges"],
}
family_data = {}

print("SELECTOR-PAIRED HALF-TRANSFER FINITE-EXACT PARTIAL SCOUT")
print("status=FINITE-EXACT_PARTIAL;no_theorem_id;static_relation_walks_only")
print("dependency_hashes=BEGIN")
for dependency, expected_hash in DEPENDENCIES.items():
    print(f"{expected_hash}  {dependency.relative_to(ROOT)}")
print("dependency_hashes=END")

for family_name, edges in family_edges.items():
    nodes, arrows, adjacency = rebuilt_graph(edges)
    family_data[family_name] = (nodes, arrows, adjacency)
    expected = EXPECTED[family_name]
    require(len(nodes) == expected["nodes"], (family_name, "node census"))
    require(len(arrows) == expected["arrows"], (family_name, "arrow census"))

    plus = tuple(
        index for index, (row, _) in enumerate(nodes) if row in selector_small
    )
    minus = tuple(
        index for index, (row, _) in enumerate(nodes) if row in selector_full
    )
    require(
        set(plus).isdisjoint(minus)
        and set(plus) | set(minus) == set(range(len(nodes))),
        (family_name, "selector split does not cover active nodes"),
    )
    plus_block = adjacency.extract(plus, plus)
    minus_block = adjacency.extract(minus, minus)
    x_block = adjacency.extract(plus, minus)
    require(
        plus_block == sp.zeros(len(plus), len(plus))
        and minus_block == sp.zeros(len(minus), len(minus)),
        (family_name, "selector split is not block-off-diagonal"),
    )
    require(
        all(
            (left[0] in selector_small) != (right[0] in selector_small)
            for left, right in arrows
        ),
        (family_name, "decorated arrow does not cross selector cut"),
    )

    degree = tuple(int(sum(adjacency.row(index))) for index in range(len(nodes)))
    plus_degree_histogram = tuple(sorted(Counter(degree[index] for index in plus).items()))
    minus_degree_histogram = tuple(sorted(Counter(degree[index] for index in minus).items()))
    require(
        len(plus_degree_histogram) > 1 and len(minus_degree_histogram) > 1,
        (family_name, "selector-only hostile unexpectedly regular"),
    )

    def degree_hostile(indices):
        labelled = tuple(sorted(
            (degree[index], node_label(nodes[index])) for index in indices
        ))
        return labelled[0], labelled[-1]

    selector_hostile = (degree_hostile(plus), degree_hostile(minus))

    refined = stable_equitable_refinement(adjacency, (plus, minus))
    require(
        len(refined) == expected["equitable_cells"],
        (family_name, "equitable cell census"),
    )
    quotient, indicator = equitable_quotient(adjacency, refined)
    cell_labels = tuple(
        tuple(node_label(nodes[index]) for index in cell) for cell in refined
    )
    quotient_rows = tuple(
        tuple(int(value) for value in quotient.row(row))
        for row in range(quotient.rows)
    )
    quotient_digest = hashlib.sha256(
        repr((cell_labels, quotient_rows)).encode("ascii")
    ).hexdigest()

    all_ones = sp.ones(len(nodes), 1)
    all_ones_row = all_ones.T
    walk_sequence = sequence(adjacency, all_ones, all_ones_row, 80)
    require(
        walk_sequence[:15] == expected["initial"],
        (family_name, "THM-3288 prefix mismatch"),
    )

    quotient_start = sp.ones(len(refined), 1)
    quotient_output = sp.Matrix([[len(cell) for cell in refined]])
    require(
        sequence(quotient, quotient_start, quotient_output, 80)
        == walk_sequence,
        (family_name, "equitable quotient lost the all-ones observer"),
    )

    order, coefficients, basis = krylov_relation(adjacency, all_ones)
    require(order == expected["order"], (family_name, "observer order"))
    characteristic = relation_polynomial(order, coefficients, z)
    companion, observer_boundary, observer_digest = companion_data(
        order, coefficients, all_ones_row, basis
    )
    require(
        sequence(
            companion,
            sp.eye(order)[:, 0],
            observer_boundary,
            80,
        ) == walk_sequence,
        (family_name, "observer companion lost walk sequence"),
    )
    hankel = sp.Matrix(
        order, order,
        lambda row, column: walk_sequence[row + column],
    )
    extended_hankel = sp.Matrix(
        order + 1, order + 1,
        lambda row, column: walk_sequence[row + column],
    )
    require(
        hankel.rank() == order and extended_hankel.rank() == order,
        (family_name, "observer minimality failure"),
    )

    quotient_order, quotient_coefficients, _ = krylov_relation(
        quotient, quotient_start
    )
    quotient_characteristic = relation_polynomial(
        quotient_order, quotient_coefficients, z
    )
    require(
        quotient_order == order and quotient_characteristic == characteristic,
        (family_name, "equitable and graph observer quotients disagree"),
    )

    half_transfer = adjacency * adjacency
    even_order, even_coefficients, even_basis = krylov_relation(
        half_transfer, all_ones
    )
    odd_start = adjacency * all_ones
    odd_order, odd_coefficients, odd_basis = krylov_relation(
        half_transfer, odd_start
    )
    require(
        (even_order, odd_order) == expected["half"],
        (family_name, "half-transfer observer orders"),
    )
    even_characteristic = relation_polynomial(even_order, even_coefficients, t)
    odd_characteristic = relation_polynomial(odd_order, odd_coefficients, t)
    even_sequence = sequence(half_transfer, all_ones, all_ones_row, 40)
    odd_sequence = sequence(half_transfer, odd_start, all_ones_row, 40)
    require(
        even_sequence == walk_sequence[::2]
        and odd_sequence == walk_sequence[1::2],
        (family_name, "half-transfer parity split"),
    )
    even_hankel = sp.Matrix(
        even_order, even_order,
        lambda row, column: even_sequence[row + column],
    )
    odd_hankel = sp.Matrix(
        odd_order, odd_order,
        lambda row, column: odd_sequence[row + column],
    )
    require(
        even_hankel.rank() == even_order and odd_hankel.rank() == odd_order,
        (family_name, "half-transfer minimality failure"),
    )

    if family_name != "full_core":
        require(
            even_characteristic == odd_characteristic,
            (family_name, "parity halves have different nonzero observer"),
        )
        expected_full = sp.Poly(
            odd_characteristic.as_expr().subs(t, z**2), z, domain=sp.QQ
        )
        require(
            characteristic == expected_full,
            (family_name, "z^2 characteristic reconstruction"),
        )
        head = matrix_polynomial_on_vector(
            odd_characteristic, half_transfer, all_ones
        )
        require(head == sp.zeros(len(nodes), 1), (family_name, "unexpected head"))
        head_summary = "zero"
    else:
        require(
            even_characteristic
            == sp.Poly(t * odd_characteristic.as_expr(), t, domain=sp.QQ),
            "full-core even half lacks the single zero mode",
        )
        expected_full = sp.Poly(
            z * odd_characteristic.as_expr().subs(t, z**2), z, domain=sp.QQ
        )
        require(
            characteristic == expected_full,
            "full-core z*p(z^2) characteristic reconstruction",
        )
        head = matrix_polynomial_on_vector(
            odd_characteristic, half_transfer, all_ones
        )
        require(
            adjacency * head == sp.zeros(len(nodes), 1),
            "full-core head is not a zero mode",
        )
        head_sum = int((all_ones_row * head)[0])
        head_norm = int((head.T * head)[0])
        head_constant = int(odd_characteristic.nth(0))
        require(
            (head_sum, head_norm, head_constant) == (-1392, 2516736, -1808),
            "full-core head invariant drift",
        )
        projection = head / head_constant
        require(
            adjacency * projection == sp.zeros(len(nodes), 1)
            and (projection.T * projection)[0]
                == (all_ones_row * projection)[0],
            "full-core cyclic zero projection identity",
        )
        for kernel_vector in adjacency.nullspace():
            require(
                ((all_ones - projection).T * kernel_vector)[0] == 0,
                "full-core zero projection is not orthogonal",
            )
        zero_mass = sp.factor((all_ones_row * projection)[0])
        require(zero_mass == sp.Rational(87, 113), "zero spectral mass drift")
        head_support = tuple(
            (node_label(node), int(value))
            for node, value in zip(nodes, head)
            if value != 0
        )
        require(
            len(head_support) == 12
            and all(row in selector_small for (row, _), _ in head_support),
            "head support left the selector-small side",
        )
        head_digest = hashlib.sha256(
            repr(head_support).encode("ascii")
        ).hexdigest()
        residuals = []
        lag_coefficients = tuple(
            odd_coefficients[odd_order - lag] for lag in range(1, odd_order + 1)
        )
        for position in range(odd_order, len(even_sequence)):
            predicted = sum(
                lag_coefficients[lag - 1] * even_sequence[position - lag]
                for lag in range(1, odd_order + 1)
            )
            residuals.append((position, int(even_sequence[position] - predicted)))
        require(
            tuple(item for item in residuals if item[1] != 0) == ((7, -1392),),
            "full-core even head residual profile",
        )
        head_summary = (
            f"sum={head_sum};norm={head_norm};p0={head_constant};"
            f"zero_mass={zero_mass};support={head_support};digest={head_digest}"
        )

    selector_sizes = (len(plus), len(minus))
    equitable_sizes = tuple(sorted(len(cell) for cell in refined))
    quotient_charpoly = sp.factor(quotient.charpoly(z).as_expr())
    print(
        f"family={family_name};active_nodes={len(nodes)};directed_arrows={len(arrows)};"
        f"selector_split={selector_sizes};X_shape={x_block.shape};"
        f"block_off_diagonal=PASS"
    )
    print(
        f"family={family_name};selector_only_lumping=FAIL;"
        f"plus_degree_histogram={plus_degree_histogram};"
        f"minus_degree_histogram={minus_degree_histogram};"
        f"hostile={selector_hostile}"
    )
    print(
        f"family={family_name};equitable_cells={len(refined)};"
        f"cell_sizes={equitable_sizes};observer_order={order};"
        f"observer_unseen_equitable_modes={len(refined)-order};"
        f"quotient_digest={quotient_digest}"
    )
    print(f"family={family_name};equitable_cell_labels={cell_labels}")
    print(f"family={family_name};equitable_quotient_rows={quotient_rows}")
    print(
        f"family={family_name};equitable_charpoly_factor={quotient_charpoly};"
        f"observer_characteristic={sp.factor(characteristic.as_expr())};"
        f"observer_companion_digest={observer_digest}"
    )
    print(
        f"family={family_name};half_transfer_orders=(even={even_order},odd={odd_order});"
        f"even_characteristic={sp.factor(even_characteristic.as_expr())};"
        f"odd_characteristic={sp.factor(odd_characteristic.as_expr())}"
    )
    print(f"family={family_name};half_transfer_head={head_summary}")
    print(
        f"family={family_name};sequence_prefix={walk_sequence[:15]};"
        f"hankel_rank={order};parity_hankel_ranks=({even_order},{odd_order})"
    )

backbone_nodes = set(family_data["backbone"][0])
tree_nodes = set(family_data["selected_tree"][0])
core_nodes = set(family_data["full_core"][0])
require(backbone_nodes < tree_nodes < core_nodes, "active-node nesting drift")
tree_new_nodes = tuple(sorted(node_label(node) for node in tree_nodes - backbone_nodes))
core_new_nodes = tuple(sorted(node_label(node) for node in core_nodes - tree_nodes))
require(len(tree_new_nodes) == 5 and len(core_new_nodes) == 3,
        "active-node creation census")

first_middle = frozenset(
    right for _, right in relation_source["oriented_relation"](2, 17)
)
second_middle = frozenset(
    left for left, _ in relation_source["oriented_relation"](17, 16)
)
require(
    tuple(sorted(state_name(state) for state in first_middle)) == ("Q-1",)
    and tuple(sorted(state_name(state) for state in second_middle)) == ("Q-8",)
    and not (first_middle & second_middle),
    "2-17-16 hostile drift",
)

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assert_nodes == 0 and float_literals == 0, "source AST gate")

print(f"active_node_creation=backbone_to_tree:{tree_new_nodes};tree_to_core:{core_new_nodes}")
print("projected_path_hostile=2-17-16;middle_fibres=(Q-1,Q-8);decorated_length2_lifts=0")
print("mechanism=selector_bipartition_gives_A_block_(0,X;X^T,0)_and_two_step_half_transfer")
print("preserved=all_static_walk_counts_via_equitable_and_Krylov_observers")
print("lost=individual_walk_identity_and_all_chronological_response_state")
print("scope=no_tournament_no_response_composition_no_GMC_FC_or_LRC_consequence_no_bit_complexity_claim")
print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")
print("ALL EXACT CHECKS PASSED")
