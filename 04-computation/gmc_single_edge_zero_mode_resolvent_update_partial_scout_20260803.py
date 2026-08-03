#!/usr/bin/env python3
"""Finite-exact single-edge update scout for the THM-3288 walk series.

Starting from THM-3277's selected eleven-edge row tree, add each of the
eleven remaining full-core row relations separately.  Rebuild the active
static maximizing-witness graph, compute its all-ones Krylov and half-transfer
orders, and measure the all-ones projection onto the adjacency kernel.

The unique one-edge addition creating a visible zero-mode head is then treated
as an exact rank-two resolvent update.  This is a static relation-walk result;
it is not a chronological response transition, tournament, FC/GMC statement,
or LRC mechanism.
"""

from __future__ import annotations

import ast
import hashlib
import io
import runpy
from contextlib import redirect_stdout
from pathlib import Path

import sympy as sp


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
THM3288 = ROOT / (
    "01-canon/theorems/"
    "THM-3288-maximizing-witness-lifted-walk-rational-series.md"
)
HALF_SOURCE = ROOT / (
    "04-computation/gmc_selector_paired_half_transfer_partial_scout_20260803.py"
)
HALF_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "gmc_selector_paired_half_transfer_partial_scout_20260803.out"
)
DEPENDENCIES = {
    THM3288:
        "7f76fb48e4ba1b2f2dce67d88c9be48aca2dd1152de10a929efa5c2a423daa58",
    HALF_SOURCE:
        "da47f9db80558343bd6442196c65150417241acbed65385b876ca09a66f0490f",
    HALF_OUTPUT:
        "ef276b6a99a3ec1ac59d5cc352a3aa9a50fa64440a6615695d684a68a6683208",
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
    half = runpy.run_path(str(HALF_SOURCE))
require(
    captured.getvalue().encode("utf-8") == lf_bytes(HALF_OUTPUT),
    "selector half-transfer replay differs from stored output",
)

source = half["relation_source"]
selected_tree = tuple(source["selected_tree"])
core_edges = tuple(source["core_edges"])
selected_set = frozenset(selected_tree)
extra_edges = tuple(edge for edge in core_edges if edge not in selected_set)
require(
    len(selected_tree) == 11
    and len(core_edges) == 22
    and len(extra_edges) == 11,
    "selected-tree/full-core edge census drift",
)


def zero_projection(adjacency):
    all_ones = sp.ones(adjacency.rows, 1)
    kernel = adjacency.nullspace()
    if not kernel:
        return sp.zeros(adjacency.rows, 1), sp.Integer(0), 0
    basis = sp.Matrix.hstack(*kernel)
    gram = basis.T * basis
    require(gram.det() != 0, "kernel Gram matrix is singular")
    projection = basis * gram.inv() * basis.T * all_ones
    require(
        adjacency * projection == sp.zeros(adjacency.rows, 1),
        "computed vector is not a kernel projection",
    )
    for kernel_vector in kernel:
        require(
            ((all_ones - projection).T * kernel_vector)[0] == 0,
            "computed kernel projection is not orthogonal",
        )
    mass = sp.factor((all_ones.T * projection)[0])
    require(
        mass == sp.factor((projection.T * projection)[0]),
        "kernel projection mass/norm identity failed",
    )
    return projection, mass, len(kernel)


def observer_data(edges):
    nodes, arrows, adjacency = half["rebuilt_graph"](edges)
    all_ones = sp.ones(adjacency.rows, 1)
    order, coefficients, basis = half["krylov_relation"](
        adjacency, all_ones
    )
    half_transfer = adjacency * adjacency
    even_order = half["krylov_relation"](
        half_transfer, all_ones
    )[0]
    odd_order = half["krylov_relation"](
        half_transfer, adjacency * all_ones
    )[0]
    projection, zero_mass, kernel_dimension = zero_projection(adjacency)
    return {
        "nodes": nodes,
        "arrows": arrows,
        "adjacency": adjacency,
        "all_ones": all_ones,
        "order": order,
        "coefficients": coefficients,
        "basis": basis,
        "even_order": even_order,
        "odd_order": odd_order,
        "projection": projection,
        "zero_mass": zero_mass,
        "kernel_dimension": kernel_dimension,
    }


def rational_observer(data, terms=90):
    adjacency = data["adjacency"]
    all_ones = data["all_ones"]
    order = data["order"]
    coefficients = data["coefficients"]
    sequence = half["sequence"](
        adjacency, all_ones, all_ones.T, terms
    )
    lag_coefficients = tuple(
        coefficients[order - lag] for lag in range(1, order + 1)
    )
    x = sp.symbols("x")
    denominator = sp.expand(
        1 - sum(
            lag_coefficients[lag - 1] * x**lag
            for lag in range(1, order + 1)
        )
    )
    numerator_coefficients = []
    for position in range(order):
        value = sp.Integer(sequence[position]) - sum(
            lag_coefficients[lag - 1] * sequence[position - lag]
            for lag in range(1, min(position, order) + 1)
        )
        numerator_coefficients.append(value)
    numerator = sp.expand(sum(
        value * x**position
        for position, value in enumerate(numerator_coefficients)
    ))
    denominator_coefficients = sp.Poly(
        denominator, x
    ).all_coeffs()[::-1]
    product = tuple(
        sp.expand(sum(
            denominator_coefficients[lag] * sequence[position - lag]
            for lag in range(
                min(position, sp.degree(denominator, x)) + 1
            )
        ))
        for position in range(len(sequence))
    )
    numerator_poly = sp.Poly(numerator, x)
    require(
        all(
            product[position] == numerator_poly.nth(position)
            for position in range(sp.degree(numerator, x) + 1)
        )
        and all(
            value == 0
            for value in product[sp.degree(numerator, x) + 1:]
        ),
        "rational observer reconstruction failed",
    )
    hankel = sp.Matrix(
        order, order,
        lambda row, column: sequence[row + column],
    )
    extended = sp.Matrix(
        order + 1, order + 1,
        lambda row, column: sequence[row + column],
    )
    require(
        hankel.rank() == order and extended.rank() == order,
        "observer Hankel minimality failed",
    )
    return sequence, denominator, numerator


EXPECTED = {
    (2, 10): (23, 42, 14, 7, 7, 9, sp.Integer(0)),
    (3, 22): (23, 42, 14, 7, 7, 9, sp.Integer(0)),
    (7, 18): (23, 42, 15, 8, 7, 9, sp.Rational(3, 5)),
    (10, 16): (24, 42, 14, 7, 7, 10, sp.Integer(0)),
    (10, 22): (23, 42, 14, 7, 7, 9, sp.Integer(0)),
    (11, 13): (24, 42, 14, 7, 7, 8, sp.Integer(0)),
    (13, 18): (23, 48, 14, 7, 7, 9, sp.Integer(0)),
    (13, 22): (23, 42, 14, 7, 7, 9, sp.Integer(0)),
    (16, 21): (24, 42, 14, 7, 7, 10, sp.Integer(0)),
    (17, 22): (23, 42, 14, 7, 7, 9, sp.Integer(0)),
    (19, 22): (23, 42, 12, 6, 6, 9, sp.Integer(0)),
}
require(tuple(EXPECTED) == extra_edges, "extra-edge order drift")

base = observer_data(selected_tree)
require(
    (
        len(base["nodes"]), len(base["arrows"]), base["order"],
        base["even_order"], base["odd_order"],
        base["kernel_dimension"], base["zero_mass"],
    ) == (23, 40, 14, 7, 7, 9, 0),
    "selected-tree control drift",
)
base_nodes = frozenset(base["nodes"])

records = []
for edge in extra_edges:
    data = observer_data(selected_tree + (edge,))
    observed = (
        len(data["nodes"]), len(data["arrows"]), data["order"],
        data["even_order"], data["odd_order"],
        data["kernel_dimension"], data["zero_mass"],
    )
    require(observed == EXPECTED[edge], ("one-edge signature drift", edge))
    relation = tuple(sorted(
        (half["state_name"](left), half["state_name"](right))
        for left, right in source["oriented_relation"](*edge)
    ))
    created_nodes = tuple(sorted(
        half["node_label"](node)
        for node in frozenset(data["nodes"]) - base_nodes
    ))
    records.append((edge, relation, created_nodes, data))

visible = tuple(
    edge for edge, _relation, _created, data in records
    if data["zero_mass"] != 0
)
require(visible == ((7, 18),), "visible one-edge zero mode is not unique")
require(
    tuple(
        (edge, created)
        for edge, _relation, created, _data in records if created
    ) == (
        ((10, 16), ((10, "Q-8"),)),
        ((11, 13), ((13, "Q-8"),)),
        ((16, 21), ((21, "Q-8"),)),
    ),
    "active-node creation controls drift",
)

special_edge, special_relation, special_created, special = next(
    record for record in records if record[0] == (7, 18)
)
require(
    special_relation == (("Q-7", "Q+1"),)
    and not special_created,
    "special rank-two edge endpoints drift",
)

base_index = {node: index for index, node in enumerate(base["nodes"])}
left_state, right_state = next(iter(source["oriented_relation"](*special_edge)))
left_node = (special_edge[0], left_state)
right_node = (special_edge[1], right_state)
left_index = base_index[left_node]
right_index = base_index[right_node]
delta = special["adjacency"] - base["adjacency"]
require(
    delta.rank() == 2
    and sum(1 for value in delta if value != 0) == 2
    and delta[left_index, right_index] == 1
    and delta[right_index, left_index] == 1,
    "special edge is not the expected rank-two adjacency update",
)

x, z, t = sp.symbols("x z t")
base_sequence, base_denominator, base_numerator = rational_observer(base)
special_sequence, special_denominator, special_numerator = rational_observer(
    special
)
special_characteristic = half["relation_polynomial"](
    special["order"], special["coefficients"], z
)
require(
    sp.factor(special_characteristic.as_expr())
    == z * (z**2 - 2) * (z**4 - 5*z**2 + 3)
       * (z**8 - 14*z**6 + 62*z**4 - 94*z**2 + 30),
    "special observer characteristic drift",
)

half_transfer = special["adjacency"] * special["adjacency"]
odd_start = special["adjacency"] * special["all_ones"]
odd_order, odd_coefficients, _odd_basis = half["krylov_relation"](
    half_transfer, odd_start
)
odd_polynomial = half["relation_polynomial"](
    odd_order, odd_coefficients, t
)
head = half["matrix_polynomial_on_vector"](
    odd_polynomial, half_transfer, special["all_ones"]
)
head_support = tuple(
    (half["node_label"](node), int(value))
    for node, value in zip(special["nodes"], head)
    if value != 0
)
require(
    special["adjacency"] * head == sp.zeros(len(special["nodes"]), 1)
    and int((special["all_ones"].T * head)[0]) == -108
    and int((head.T * head)[0]) == 19440
    and int(odd_polynomial.nth(0)) == -180
    and special["zero_mass"] == sp.Rational(3, 5),
    "special zero-mode head invariants drift",
)
require(
    head / odd_polynomial.nth(0) == special["projection"],
    "head is not the exact all-ones zero projection multiple",
)
require(
    head_support == (
        ((18, "Q+1"), 108),
        ((18, "Q+2"), -36),
        ((18, "Q+4"), -36),
        ((18, "Q+5"), -36),
        ((22, "Q+1"), -36),
        ((22, "Q+2"), -36),
        ((22, "Q+5"), -36),
    ),
    "special head support drift",
)

# Combinatorial mechanism for the kernel vector.  The new edge gives the
# central small-side vertex two neighbours.  Three same-row siblings see only
# one neighbour and three row-22 siblings see only the other, so their
# incidence vectors satisfy one signed 3-versus-6 relation.  Unlike ordinary
# twin differences, its coefficient sum is nonzero and is therefore visible
# to the all-ones observer.
central_label = (18, "Q+1")
left_neighbour_label = (7, "Q-7")
right_neighbour_label = (19, "Q-1")
right_only_labels = ((18, "Q+2"), (18, "Q+4"), (18, "Q+5"))
left_only_labels = ((22, "Q+1"), (22, "Q+2"), (22, "Q+5"))
special_index = {
    half["node_label"](node): index
    for index, node in enumerate(special["nodes"])
}


def neighbour_labels(label):
    index = special_index[label]
    return tuple(
        half["node_label"](special["nodes"][target])
        for target in range(len(special["nodes"]))
        if special["adjacency"][index, target] != 0
    )


require(
    neighbour_labels(central_label)
        == (left_neighbour_label, right_neighbour_label)
    and all(
        neighbour_labels(label) == (right_neighbour_label,)
        for label in right_only_labels
    )
    and all(
        neighbour_labels(label) == (left_neighbour_label,)
        for label in left_only_labels
    ),
    "special local neighbour profile drift",
)
central_basis = sp.eye(len(special["nodes"]))[:, special_index[central_label]]
leaf_sum = sum(
    (
        sp.eye(len(special["nodes"]))[:, special_index[label]]
        for label in right_only_labels + left_only_labels
    ),
    sp.zeros(len(special["nodes"]), 1),
)
primitive_head = 3 * central_basis - leaf_sum
require(
    special["adjacency"] * primitive_head
        == sp.zeros(len(special["nodes"]), 1)
    and head == 36 * primitive_head
    and sum(primitive_head) == -3,
    "special 3-versus-6 incidence mechanism failed",
)

even_sequence = special_sequence[::2]
odd_lag_coefficients = tuple(
    odd_coefficients[odd_order - lag]
    for lag in range(1, odd_order + 1)
)
even_residuals = tuple(
    (
        position,
        sp.expand(
            even_sequence[position]
            - sum(
                odd_lag_coefficients[lag - 1]
                * even_sequence[position - lag]
                for lag in range(1, odd_order + 1)
            )
        ),
    )
    for position in range(odd_order, len(even_sequence))
)
require(
    tuple(item for item in even_residuals if item[1] != 0) == ((7, -108),),
    "special one-step recurrence residual drift",
)

# Exact scalar Woodbury update.  For U=[e_left,e_right] and
# C=((0,1),(1,0)), Delta=U C U^T and C^{-1}=C.  With
# R=(I-xA_tree)^{-1}, s=U^T R 1 and K=U^T R U,
#
# G_new-G_tree = x s^T (C-xK)^{-1} s.
base_resolvent = (sp.eye(base["adjacency"].rows) - x * base["adjacency"]).inv()
u_matrix = sp.zeros(base["adjacency"].rows, 2)
u_matrix[left_index, 0] = 1
u_matrix[right_index, 1] = 1
c_matrix = sp.Matrix(((0, 1), (1, 0)))
s_vector = sp.simplify(u_matrix.T * base_resolvent * base["all_ones"])
k_matrix = sp.simplify(u_matrix.T * base_resolvent * u_matrix)
woodbury_correction = sp.factor(
    x * (s_vector.T * (c_matrix - x * k_matrix).inv() * s_vector)[0]
)
rational_difference = sp.factor(
    special_numerator / special_denominator
    - base_numerator / base_denominator
)
require(
    sp.cancel(woodbury_correction - rational_difference) == 0,
    "rank-two Woodbury observer update failed",
)

local_denominator = 24*x**8 - 77*x**6 + 54*x**4 - 13*x**2 + 1
new_local_denominator = 30*x**8 - 94*x**6 + 62*x**4 - 14*x**2 + 1
require(
    all(
        sp.denom(sp.factor(entry)) == local_denominator
        for entry in tuple(s_vector) + tuple(k_matrix)
    ),
    "local resolvent denominator drift",
)
correction_numerator, correction_denominator = sp.together(
    woodbury_correction
).as_numer_denom()
require(
    sp.factor(correction_denominator)
    == local_denominator * new_local_denominator,
    "Woodbury correction denominator drift",
)

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(
    assert_nodes == 0 and float_literals == 0,
    "source AST gate found assert or floating-point literal",
)

print("SINGLE-EDGE ZERO-MODE RESOLVENT UPDATE FINITE-EXACT PARTIAL SCOUT")
print("status=FINITE-EXACT_PARTIAL;no_theorem_id;static_relation_walks_only")
print("dependency_hashes=BEGIN")
for dependency, expected_hash in DEPENDENCIES.items():
    print(f"{expected_hash}  {dependency.relative_to(ROOT)}")
print("dependency_hashes=END")
print(
    "base=selected_tree;active_nodes=23;directed_arrows=40;"
    "observer_order=14;half_orders=(7,7);kernel_dimension=9;zero_mass=0"
)
for edge, relation, created_nodes, data in records:
    print(
        f"edge={edge};relation={relation};created_nodes={created_nodes};"
        f"active_nodes={len(data['nodes'])};directed_arrows={len(data['arrows'])};"
        f"observer_order={data['order']};"
        f"half_orders=({data['even_order']},{data['odd_order']});"
        f"kernel_dimension={data['kernel_dimension']};"
        f"zero_mass={data['zero_mass']}"
    )
print("unique_visible_zero_mode_edge=(7,18);relation=((Q-7,Q+1),);rank_delta=2")
print(f"special_prefix={special_sequence[:18]}")
print(
    f"special_characteristic={sp.factor(special_characteristic.as_expr())};"
    f"hankel_rank={special['order']}"
)
print(f"special_generating_denominator={sp.factor(special_denominator)}")
print(f"special_generating_numerator={sp.factor(special_numerator)}")
print(
    f"special_odd_half_polynomial={sp.factor(odd_polynomial.as_expr())};"
    "even_half_polynomial=t*special_odd_half_polynomial"
)
print(
    "special_head=sum:-108;norm:19440;p0:-180;zero_mass:3/5;"
    f"support={head_support};first_even_residual=(7,-108)"
)
print(
    "special_head_mechanism=3*N(18,Q+1)="
    "sum_N((18,Q+2),(18,Q+4),(18,Q+5),(22,Q+1),(22,Q+2),(22,Q+5));"
    "central_neighbours=((7,Q-7),(19,Q-1));primitive_coefficient_sum=-3"
)
print(f"woodbury_local_s={tuple(sp.factor(entry) for entry in s_vector)}")
print(
    "woodbury_local_K="
    f"{tuple(tuple(sp.factor(entry) for entry in k_matrix.row(row)) for row in range(2))}"
)
print(f"woodbury_correction_numerator={sp.factor(correction_numerator)}")
print(f"woodbury_correction_denominator={sp.factor(correction_denominator)}")
print("mechanism=one_static_relation_edge_is_a_rank_two_resolvent_update")
print("preserved=all_ones_static_walk_series_and_exact_active_node_boundary")
print("lost=individual_walk_identity_row_path_simplicity_and_physical_chronology")
print("scope=no_tournament_no_response_composition_no_GMC_FC_or_LRC_consequence")
print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")
print("ALL EXACT CHECKS PASSED")
