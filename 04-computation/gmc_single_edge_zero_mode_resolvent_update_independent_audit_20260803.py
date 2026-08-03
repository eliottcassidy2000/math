#!/usr/bin/env python3
"""Independent audit of the finite-exact single-edge resolvent scout.

This companion imports or executes no newer single-edge or half-transfer
scout.  It pins and replays only THM-3287's inherited relation artifact, then
independently rebuilds the selected-tree adjacency and all eleven one-edge
extensions.  Separate exact linear algebra reconstructs the Krylov orders,
kernel projections, local incidence certificate, rational observers, and the
rank-two Woodbury update.

All graphs here encode static maximizing-witness compatibility.  They are not
chronological response transitions, tournaments, or FC/GMC/LRC mechanisms.
"""

from __future__ import annotations

import ast
import hashlib
import io
import runpy
from collections import Counter
from contextlib import redirect_stdout
from pathlib import Path

import sympy as sp


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
THM3287 = ROOT / (
    "01-canon/theorems/"
    "THM-3287-weighted-backbone-dominance-witness-section-and-selector-cut.md"
)
RELATION_SOURCE = ROOT / (
    "04-computation/gmc_backbone_maximizing_witness_section_scout_20260803.py"
)
RELATION_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "gmc_backbone_maximizing_witness_section_scout_20260803.out"
)
DEPENDENCIES = {
    THM3287:
        "3cb8258eb93dd1b3bd59d664ce42c819ff31ca020cc34c57c98e4cf3946512c7",
    RELATION_SOURCE:
        "aed89d67ec7acabfe5b4feae4a83f7c57b78053928be44a6cc319d81fa4a9cc6",
    RELATION_OUTPUT:
        "89200bc6cff7284dd33f636352f4f7f56294d90bcd2902fa96092d9f967f5fe3",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected_hash in DEPENDENCIES.items():
    require(
        hashlib.sha256(lf_bytes(dependency)).hexdigest() == expected_hash,
        ("dependency hash drift", dependency.name),
    )

# Replay THM-3287's exact relation constructor and its own transitive gates.
# No graph, recurrence, projection, or resolvent helper is imported from a
# later artifact.
captured = io.StringIO()
with redirect_stdout(captured):
    inherited = runpy.run_path(str(RELATION_SOURCE))
require(
    captured.getvalue().encode("utf-8") == lf_bytes(RELATION_OUTPUT),
    "THM-3287 relation replay differs from frozen output",
)

reset = tuple(inherited["reset"])
values = tuple(inherited["values"])
relation_bank = dict(inherited["relation_bank"])
selected_tree = tuple(inherited["selected_tree"])
core_edges = tuple(inherited["core_edges"])
selected_set = frozenset(selected_tree)
extra_edges = tuple(edge for edge in core_edges if edge not in selected_set)
require(
    len(selected_tree) == 11
    and len(core_edges) == 22
    and len(extra_edges) == 11,
    "selected-tree/full-core census drift",
)

reset_counts = Counter(reset)


def state_name(state):
    counts = Counter(state)
    added = tuple(
        value
        for value in values
        for _ in range(max(0, counts[value] - reset_counts[value]))
    )
    removed = tuple(
        value
        for value in values
        for _ in range(max(0, reset_counts[value] - counts[value]))
    )
    if len(added) == 1 and not removed:
        return f"Q+{added[0]}"
    if len(removed) == 1 and not added:
        return f"Q-{removed[0]}"
    raise RuntimeError(("non-link state", state))


def node_label(node):
    row, state = node
    return row, state_name(state)


def oriented_relation(first, second):
    edge = (first, second) if first < second else (second, first)
    relation = relation_bank[edge]
    if first < second:
        return frozenset(relation)
    return frozenset((right, left) for left, right in relation)


def rebuild_graph(edges):
    nodes = set()
    arrows = set()
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
        adjacency[index[left], index[right]] = 1
    require(adjacency == adjacency.T, "rebuilt graph is not symmetric")
    require(
        sum(1 for value in adjacency if value != 0) == len(arrows),
        "rebuilt directed-arrow census drift",
    )
    return nodes, tuple(sorted(arrows)), adjacency


def krylov_closure(transfer, start):
    columns = []
    vector = start
    for order in range(transfer.rows + 2):
        if columns:
            basis = sp.Matrix.hstack(*columns)
            if basis.row_join(vector).rank() == len(columns):
                coefficients = basis.gauss_jordan_solve(vector)[0]
                require(
                    basis * coefficients == vector,
                    "independent Krylov solve failed",
                )
                return order, tuple(coefficients), basis
        columns.append(vector)
        vector = transfer * vector
    raise RuntimeError("independent Krylov closure failed")


def exact_sequence(transfer, start, output, terms):
    values = []
    vector = start
    for _ in range(terms):
        values.append(sp.expand((output * vector)[0]))
        vector = transfer * vector
    return tuple(values)


def zero_projection(adjacency):
    one = sp.ones(adjacency.rows, 1)
    kernel = adjacency.nullspace()
    if not kernel:
        return 0, sp.Integer(0), sp.zeros(adjacency.rows, 1)
    basis = sp.Matrix.hstack(*kernel)
    gram = basis.T * basis
    require(gram.det() != 0, "kernel Gram matrix is singular")
    projection = basis * gram.inv() * basis.T * one
    require(
        adjacency * projection == sp.zeros(adjacency.rows, 1),
        "kernel projection left the kernel",
    )
    require(
        all(
            ((one - projection).T * vector)[0] == 0
            for vector in kernel
        ),
        "kernel projection is not orthogonal",
    )
    mass = sp.factor((one.T * projection)[0])
    require(
        mass == sp.factor((projection.T * projection)[0]),
        "kernel projection mass/norm identity failed",
    )
    return len(kernel), mass, projection


def observer_data(edges):
    nodes, arrows, adjacency = rebuild_graph(edges)
    one = sp.ones(adjacency.rows, 1)
    order, coefficients, basis = krylov_closure(adjacency, one)
    square = adjacency * adjacency
    even_order = krylov_closure(square, one)[0]
    odd_order = krylov_closure(square, adjacency * one)[0]
    kernel_dimension, zero_mass, projection = zero_projection(adjacency)
    return {
        "nodes": nodes,
        "arrows": arrows,
        "adjacency": adjacency,
        "one": one,
        "order": order,
        "coefficients": coefficients,
        "basis": basis,
        "even_order": even_order,
        "odd_order": odd_order,
        "kernel_dimension": kernel_dimension,
        "zero_mass": zero_mass,
        "projection": projection,
    }


def signature(data):
    return (
        len(data["nodes"]),
        len(data["arrows"]),
        data["order"],
        data["even_order"],
        data["odd_order"],
        data["kernel_dimension"],
        data["zero_mass"],
    )


def relation_polynomial(order, coefficients, symbol):
    return sp.Poly(
        sp.expand(
            symbol**order
            - sum(
                coefficients[index] * symbol**index
                for index in range(order)
            )
        ),
        symbol,
        domain=sp.QQ,
    )


def polynomial_on_vector(polynomial, matrix, vector):
    answer = sp.zeros(matrix.rows, 1)
    power = vector
    for degree in range(polynomial.degree() + 1):
        answer += polynomial.nth(degree) * power
        power = matrix * power
    return answer


def rational_observer(data, terms=80):
    adjacency = data["adjacency"]
    one = data["one"]
    order = data["order"]
    coefficients = data["coefficients"]
    sequence = exact_sequence(adjacency, one, one.T, terms)
    x = sp.symbols("x")
    lag_coefficients = tuple(
        coefficients[order - lag]
        for lag in range(1, order + 1)
    )
    denominator = sp.expand(
        1
        - sum(
            lag_coefficients[lag - 1] * x**lag
            for lag in range(1, order + 1)
        )
    )
    numerator = sp.expand(sum(
        (
            sequence[position]
            - sum(
                lag_coefficients[lag - 1] * sequence[position - lag]
                for lag in range(1, position + 1)
            )
        ) * x**position
        for position in range(order)
    ))
    denominator_coefficients = tuple(
        sp.Poly(denominator, x).nth(degree)
        for degree in range(sp.degree(denominator, x) + 1)
    )
    product = tuple(
        sp.expand(sum(
            denominator_coefficients[lag] * sequence[position - lag]
            for lag in range(
                min(position, len(denominator_coefficients) - 1) + 1
            )
        ))
        for position in range(len(sequence))
    )
    numerator_poly = sp.Poly(numerator, x)
    require(
        all(
            product[position] == numerator_poly.nth(position)
            for position in range(numerator_poly.degree() + 1)
        )
        and all(
            value == 0
            for value in product[numerator_poly.degree() + 1:]
        ),
        "independent rational-series reconstruction failed",
    )
    leading_hankel = sp.Matrix(
        order,
        order,
        lambda row, column: sequence[row + column],
    )
    extended_hankel = sp.Matrix(
        order + 1,
        order + 1,
        lambda row, column: sequence[row + column],
    )
    require(
        leading_hankel.rank() == order
        and extended_hankel.rank() == order,
        "independent Hankel rank/minimality failed",
    )
    return sequence, numerator, denominator


EXPECTED = {
    (2, 10): (
        (("Q+4", "Q-1"),), (),
        (23, 42, 14, 7, 7, 9, sp.Integer(0)),
    ),
    (3, 22): (
        (("Q-1", "Q+4"),), (),
        (23, 42, 14, 7, 7, 9, sp.Integer(0)),
    ),
    (7, 18): (
        (("Q-7", "Q+1"),), (),
        (23, 42, 15, 8, 7, 9, sp.Rational(3, 5)),
    ),
    (10, 16): (
        (("Q-8", "Q+1"),), ((10, "Q-8"),),
        (24, 42, 14, 7, 7, 10, sp.Integer(0)),
    ),
    (10, 22): (
        (("Q-1", "Q+4"),), (),
        (23, 42, 14, 7, 7, 9, sp.Integer(0)),
    ),
    (11, 13): (
        (("Q+1", "Q-8"),), ((13, "Q-8"),),
        (24, 42, 14, 7, 7, 8, sp.Integer(0)),
    ),
    (13, 18): (
        (
            ("Q-1", "Q+1"),
            ("Q-1", "Q+2"),
            ("Q-1", "Q+4"),
            ("Q-1", "Q+5"),
        ),
        (),
        (23, 48, 14, 7, 7, 9, sp.Integer(0)),
    ),
    (13, 22): (
        (("Q-1", "Q+4"),), (),
        (23, 42, 14, 7, 7, 9, sp.Integer(0)),
    ),
    (16, 21): (
        (("Q+1", "Q-8"),), ((21, "Q-8"),),
        (24, 42, 14, 7, 7, 10, sp.Integer(0)),
    ),
    (17, 22): (
        (("Q-1", "Q+4"),), (),
        (23, 42, 14, 7, 7, 9, sp.Integer(0)),
    ),
    (19, 22): (
        (("Q-1", "Q+4"),), (),
        (23, 42, 12, 6, 6, 9, sp.Integer(0)),
    ),
}
require(tuple(EXPECTED) == extra_edges, "extra-edge order drift")

base = observer_data(selected_tree)
require(
    signature(base) == (23, 40, 14, 7, 7, 9, sp.Integer(0)),
    "selected-tree control drift",
)
base_nodes = frozenset(base["nodes"])

records = []
for edge in extra_edges:
    data = observer_data(selected_tree + (edge,))
    relation = tuple(sorted(
        (state_name(left), state_name(right))
        for left, right in oriented_relation(*edge)
    ))
    created_nodes = tuple(sorted(
        node_label(node)
        for node in frozenset(data["nodes"]) - base_nodes
    ))
    observed = (relation, created_nodes, signature(data))
    require(observed == EXPECTED[edge], ("one-edge audit drift", edge))
    records.append((edge, relation, created_nodes, data))

visible_edges = tuple(
    edge
    for edge, _relation, _created, data in records
    if data["zero_mass"] != 0
)
require(
    visible_edges == ((7, 18),),
    "visible zero-mode edge is not unique in the eleven-case universe",
)

special = next(data for edge, _relation, _created, data in records
               if edge == (7, 18))
require(
    special["nodes"] == base["nodes"],
    "special edge changed the active-node order",
)
base_index = {node: index for index, node in enumerate(base["nodes"])}
left_state, right_state = next(iter(oriented_relation(7, 18)))
left_node = (7, left_state)
right_node = (18, right_state)
left_index = base_index[left_node]
right_index = base_index[right_node]
delta = special["adjacency"] - base["adjacency"]
require(
    delta.rank() == 2
    and sum(1 for value in delta if value != 0) == 2
    and delta[left_index, right_index] == 1
    and delta[right_index, left_index] == 1,
    "special adjacency difference is not the claimed rank-two update",
)

z, t, x = sp.symbols("z t x")
special_characteristic = relation_polynomial(
    special["order"], special["coefficients"], z
)
require(
    sp.factor(special_characteristic.as_expr())
    == z * (z**2 - 2) * (z**4 - 5*z**2 + 3)
       * (z**8 - 14*z**6 + 62*z**4 - 94*z**2 + 30),
    "special observer characteristic drift",
)

half_transfer = special["adjacency"] * special["adjacency"]
odd_order, odd_coefficients, _odd_basis = krylov_closure(
    half_transfer, special["adjacency"] * special["one"]
)
odd_polynomial = relation_polynomial(odd_order, odd_coefficients, t)
head = polynomial_on_vector(
    odd_polynomial, half_transfer, special["one"]
)
head_support = tuple(
    (node_label(node), int(value))
    for node, value in zip(special["nodes"], head)
    if value != 0
)
require(
    sp.factor(odd_polynomial.as_expr())
    == (t - 2) * (t**2 - 5*t + 3)
       * (t**4 - 14*t**3 + 62*t**2 - 94*t + 30),
    "special odd-half polynomial drift",
)
require(
    special["adjacency"] * head
        == sp.zeros(len(special["nodes"]), 1)
    and int((special["one"].T * head)[0]) == -108
    and int((head.T * head)[0]) == 19440
    and int(odd_polynomial.nth(0)) == -180
    and head / odd_polynomial.nth(0) == special["projection"],
    "special zero-mode head/projection drift",
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

special_index = {
    node_label(node): index
    for index, node in enumerate(special["nodes"])
}
central = (18, "Q+1")
right_only = ((18, "Q+2"), (18, "Q+4"), (18, "Q+5"))
left_only = ((22, "Q+1"), (22, "Q+2"), (22, "Q+5"))
left_neighbour = (7, "Q-7")
right_neighbour = (19, "Q-1")


def neighbour_labels(label):
    row = special_index[label]
    return tuple(
        node_label(special["nodes"][column])
        for column in range(len(special["nodes"]))
        if special["adjacency"][row, column] != 0
    )


require(
    neighbour_labels(central) == (left_neighbour, right_neighbour)
    and all(neighbour_labels(label) == (right_neighbour,)
            for label in right_only)
    and all(neighbour_labels(label) == (left_neighbour,)
            for label in left_only),
    "three-versus-six local neighbour profile drift",
)
primitive_kernel = sp.zeros(len(special["nodes"]), 1)
primitive_kernel[special_index[central]] = 3
for label in right_only + left_only:
    primitive_kernel[special_index[label]] = -1
require(
    special["adjacency"] * primitive_kernel
        == sp.zeros(len(special["nodes"]), 1)
    and sum(primitive_kernel) == -3
    and (primitive_kernel.T * primitive_kernel)[0] == 15
    and head == 36 * primitive_kernel
    and special["projection"] == -primitive_kernel / 5
    and special["zero_mass"] == sp.Rational(3, 5),
    "three-versus-six kernel/mass certificate failed",
)

base_sequence, base_numerator, base_denominator = rational_observer(base)
special_sequence, special_numerator, special_denominator = rational_observer(
    special
)
require(
    special_sequence[:18] == (
        23, 42, 116, 252, 692, 1602, 4404, 10518, 28966,
        70344, 194026, 475722, 1313706, 3239490, 8953208,
        22155444, 61266026, 151938198,
    ),
    "special walk prefix drift",
)

# Direct symbolic rank-two resolvent audit.  If B is the selected-tree
# adjacency and U selects the two old endpoint vertices, then
# A=B+U*C*U^T with C the two-by-two exchange matrix.  Woodbury predicts
# G_A-G_B=x*s^T*(C-x*K)^(-1)*s.
base_resolvent = (
    sp.eye(base["adjacency"].rows) - x * base["adjacency"]
).inv()
u_matrix = sp.zeros(base["adjacency"].rows, 2)
u_matrix[left_index, 0] = 1
u_matrix[right_index, 1] = 1
c_matrix = sp.Matrix(((0, 1), (1, 0)))
s_vector = sp.simplify(u_matrix.T * base_resolvent * base["one"])
k_matrix = sp.simplify(u_matrix.T * base_resolvent * u_matrix)
woodbury_correction = sp.cancel(
    x * (s_vector.T * (c_matrix - x * k_matrix).inv() * s_vector)[0]
)
direct_difference = sp.cancel(
    special_numerator / special_denominator
    - base_numerator / base_denominator
)
require(
    sp.cancel(woodbury_correction - direct_difference) == 0,
    "independent Woodbury/direct rational observers disagree",
)

common_denominator = (1 - 2*x**2) * (1 - 5*x**2 + 3*x**4)
old_local_denominator = 1 - 13*x**2 + 54*x**4 - 77*x**6 + 24*x**8
new_local_denominator = 1 - 14*x**2 + 62*x**4 - 94*x**6 + 30*x**8
require(
    sp.expand(base_denominator)
        == sp.expand(common_denominator * old_local_denominator)
    and sp.expand(special_denominator)
        == sp.expand(common_denominator * new_local_denominator),
    "old/new observer factor localization drift",
)
require(
    all(
        sp.denom(sp.factor(entry)) == old_local_denominator
        for entry in tuple(s_vector) + tuple(k_matrix)
    ),
    "local resolvent entries do not share the old local denominator",
)
correction_numerator, correction_denominator = sp.fraction(
    woodbury_correction
)
require(
    sp.factor(correction_denominator)
        == old_local_denominator * new_local_denominator
    and sp.Poly(correction_numerator, x).gcd(
        sp.Poly(correction_denominator, x)
    ).degree() == 0,
    "reduced Woodbury correction denominator drift",
)

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(
    isinstance(node, ast.Assert) for node in ast.walk(syntax)
)
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(
    assert_nodes == 0 and float_literals == 0,
    "source AST gate found assert or floating-point literal",
)

print("SINGLE-EDGE ZERO-MODE RESOLVENT UPDATE INDEPENDENT EXACT AUDIT")
print("status=INDEPENDENT_AUDIT_OF_FINITE_EXACT_PARTIAL;no_theorem_id")
print("dependency_hashes=BEGIN")
for dependency, expected_hash in DEPENDENCIES.items():
    print(f"{expected_hash}  {dependency.relative_to(ROOT)}")
print("dependency_hashes=END")
print("imports_or_executes_new_single_edge_or_half_transfer_scout=NO")
print(
    "base=selected_tree;signature="
    "(nodes=23,arrows=40,order=14,half=(7,7),kernel=9,mass=0)"
)
for edge, relation, created_nodes, data in records:
    observed = signature(data)
    print(
        f"edge={edge};relation={relation};created_nodes={created_nodes};"
        f"signature=(nodes={observed[0]},arrows={observed[1]},"
        f"order={observed[2]},half=({observed[3]},{observed[4]}),"
        f"kernel={observed[5]},mass={observed[6]})"
    )
print("unique_visible_edge=(7,18);relation=((Q-7,Q+1),);delta_rank=2")
print(f"special_prefix={special_sequence[:18]}")
print(
    f"special_characteristic={sp.factor(special_characteristic.as_expr())};"
    "hankel_rank=15"
)
print(f"special_generating_denominator={sp.factor(special_denominator)}")
print(f"special_generating_numerator={sp.factor(special_numerator)}")
print(
    f"special_odd_half_polynomial={sp.factor(odd_polynomial.as_expr())};"
    "half_orders=(8,7)"
)
print(
    "three_vs_six_kernel=(sum=-3,norm=15);head_factor=36;"
    "head=(sum=-108,norm=19440,p0=-180);projection=-kernel/5;mass=3/5"
)
print(f"special_head_support={head_support}")
print(
    "three_vs_six_neighbours=central:((7,Q-7),(19,Q-1));"
    "row18_siblings_only:(19,Q-1);row22_siblings_only:(7,Q-7)"
)
print(f"woodbury_local_s={tuple(sp.factor(entry) for entry in s_vector)}")
print(
    "woodbury_local_K="
    f"{tuple(tuple(sp.factor(k_matrix[row, column]) for column in range(2)) for row in range(2))}"
)
print(f"woodbury_correction_numerator={sp.factor(correction_numerator)}")
print(f"woodbury_correction_denominator={sp.factor(correction_denominator)}")
print("woodbury_direct_rational_identity=PASS;correction_fraction_reduced=PASS")
print("preserved=all_ones_static_relation_walk_series_and_active_node_boundary")
print("lost=walk_identity_path_simplicity_ratio_endpoints_and_physical_chronology")
print("scope=no_tournament_no_response_composition_no_GMC_FC_or_LRC_consequence")
print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")
print("ALL INDEPENDENT EXACT CHECKS PASSED")
