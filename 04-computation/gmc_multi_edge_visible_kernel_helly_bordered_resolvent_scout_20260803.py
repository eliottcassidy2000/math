#!/usr/bin/env python3
"""Finite-exact multi-edge visibility and bordered-resolvent scout.

Close THM-3305's eleven one-edge extensions under every pair and under the
full Boolean lattice of edge choices.  The scout finds the inclusion-minimal
sets whose adjacency kernel is visible to the all-ones observer, distinguishes
one persistent circuit from five fragile arity-three circuits, and checks a
bordered Woodbury coefficient compiler on every two-edge update.

Everything here concerns symmetric static relation-walk graphs.  It is not a
tournament, chronology, response-composition, FC/GMC, or LRC result.
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


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
THM3305 = ROOT / (
    "01-canon/theorems/"
    "THM-3305-single-edge-visible-zero-mode-and-rank-two-resolvent-update.md"
)
SINGLE_SOURCE = ROOT / (
    "04-computation/"
    "gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.py"
)
SINGLE_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.out"
)
DEPENDENCIES = {
    THM3305:
        "fd9c4feaa9a39db14c754a3ff7267ad5f71e2cd71953119c07df53fafb96bab6",
    SINGLE_SOURCE:
        "5704ee8d6641a12e50d6a6179ccf58f911ba782cfe369ac70a32db51887565b2",
    SINGLE_OUTPUT:
        "b52db8c18d3182fb79ed5a230d9f75775155f438c7f32e974ab1ff3fcb15e556",
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
    single = runpy.run_path(str(SINGLE_SOURCE))
require(
    captured.getvalue().encode("utf-8") == lf_bytes(SINGLE_OUTPUT),
    "THM-3305 scout replay differs from stored output",
)

half = single["half"]
source = half["relation_source"]
selected_tree = tuple(single["selected_tree"])
core_edges = tuple(single["core_edges"])
extra_edges = tuple(single["extra_edges"])
observer_data = single["observer_data"]
zero_projection = single["zero_projection"]
require(
    len(selected_tree) == len(extra_edges) == 11
    and len(core_edges) == 22,
    "tree/core edge census drift",
)


def signature(data):
    return (
        len(data["nodes"]), len(data["arrows"]), data["order"],
        data["even_order"], data["odd_order"],
        data["kernel_dimension"], data["zero_mass"],
    )


# Pairwise closure.  No pair creates a new visible circuit: the ten visible
# pairs are exactly those containing THM-3305's special edge (7,18).
pair_records = []
base_node_set = frozenset(observer_data(selected_tree)["nodes"])
for first, second in combinations(extra_edges, 2):
    data = observer_data(selected_tree + (first, second))
    pair_records.append((first, second, signature(data), data))

pair_signature_counts = Counter(record[2] for record in pair_records)
EXPECTED_PAIR_SIGNATURE_COUNTS = Counter({
    (23, 44, 12, 6, 6, 11, sp.Integer(0)): 1,
    (23, 44, 12, 6, 6, 9, sp.Integer(0)): 2,
    (23, 44, 14, 7, 7, 9, sp.Integer(0)): 12,
    (23, 44, 15, 8, 7, 9, sp.Rational(3, 5)): 6,
    (23, 50, 12, 6, 6, 11, sp.Integer(0)): 1,
    (23, 50, 14, 7, 7, 9, sp.Integer(0)): 5,
    (23, 50, 15, 8, 7, 9, sp.Rational(3, 5)): 1,
    (24, 44, 12, 6, 6, 8, sp.Integer(0)): 1,
    (24, 44, 14, 7, 7, 10, sp.Integer(0)): 12,
    (24, 44, 14, 7, 7, 8, sp.Integer(0)): 4,
    (24, 44, 15, 8, 7, 10, sp.Rational(3, 5)): 2,
    (24, 44, 15, 8, 7, 8, sp.Rational(3, 5)): 1,
    (24, 44, 16, 8, 8, 8, sp.Integer(0)): 1,
    (24, 50, 14, 7, 7, 10, sp.Integer(0)): 2,
    (24, 50, 14, 7, 7, 8, sp.Integer(0)): 1,
    (25, 44, 14, 7, 7, 11, sp.Integer(0)): 1,
    (25, 44, 14, 7, 7, 9, sp.Integer(0)): 2,
})
require(
    pair_signature_counts == EXPECTED_PAIR_SIGNATURE_COUNTS,
    "two-edge observer-signature census drift",
)
visible_pairs = tuple(
    (first, second)
    for first, second, pair_signature, _data in pair_records
    if pair_signature[-1] != 0
)
require(
    len(visible_pairs) == 10
    and all((7, 18) in pair for pair in visible_pairs),
    "pairwise visibility is not exactly inherited from edge (7,18)",
)
pair_serialized = "".join(
    f"{first};{second};" + ";".join(str(value) for value in pair_signature) + "\n"
    for first, second, pair_signature, _data in pair_records
).encode("utf-8")
require(
    len(pair_serialized) == 1934
    and hashlib.sha256(pair_serialized).hexdigest()
        == "86badaa9186021f92c026a591b29436790335d76468db2c18fa260c7feef36f0",
    "two-edge serialized census digest drift",
)


# The full 2^11 edge lattice is small enough for an exact census.  Variable
# active-node boundaries are retained rather than padding before projection.
lattice_records = []
mass_by_mask = {}
for mask in range(1 << len(extra_edges)):
    chosen = tuple(
        extra_edges[index]
        for index in range(len(extra_edges))
        if mask & (1 << index)
    )
    nodes, arrows, adjacency = half["rebuilt_graph"](selected_tree + chosen)
    _projection, mass, kernel_dimension = zero_projection(adjacency)
    record = (
        mask, mask.bit_count(), len(nodes), len(arrows),
        kernel_dimension, mass,
    )
    lattice_records.append(record)
    mass_by_mask[mask] = mass

lattice_serialized = "".join(
    ";".join(str(value) for value in record) + "\n"
    for record in lattice_records
).encode("utf-8")
require(
    len(lattice_serialized) == 37894
    and hashlib.sha256(lattice_serialized).hexdigest()
        == "53f740b8fd3e8b1baa88a460e210662d51b6aa69ff8a157fcdfc9fd5d8510ad0",
    "Boolean-lattice serialized census digest drift",
)

size_summary = {}
records_by_size = defaultdict(list)
for record in lattice_records:
    records_by_size[record[1]].append(record)
for size, records in records_by_size.items():
    nonzero_masses = tuple(record[-1] for record in records if record[-1] != 0)
    size_summary[size] = (
        len(records), len(records) - len(nonzero_masses), len(nonzero_masses),
        len(frozenset(nonzero_masses)),
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
require(size_summary == EXPECTED_SIZE_SUMMARY, "edge-lattice size summary drift")


def edge_set(mask):
    return tuple(
        extra_edges[index]
        for index in range(len(extra_edges))
        if mask & (1 << index)
    )


minimal_visible_masks = []
for mask in range(1, 1 << len(extra_edges)):
    if mass_by_mask[mask] == 0:
        continue
    submask = (mask - 1) & mask
    has_visible_proper_subset = False
    while submask:
        if mass_by_mask[submask] != 0:
            has_visible_proper_subset = True
            break
        submask = (submask - 1) & mask
    if not has_visible_proper_subset:
        minimal_visible_masks.append(mask)

minimal_visible_sets = tuple(edge_set(mask) for mask in minimal_visible_masks)
EXPECTED_MINIMAL_VISIBLE_SETS = frozenset({
    ((7, 18),),
    ((2, 10), (13, 18), (17, 22)),
    ((2, 10), (13, 22), (17, 22)),
    ((10, 22), (13, 18), (17, 22)),
    ((10, 22), (13, 22), (17, 22)),
    ((13, 22), (17, 22), (19, 22)),
})
require(
    len(minimal_visible_sets) == 6
    and frozenset(minimal_visible_sets) == EXPECTED_MINIMAL_VISIBLE_SETS,
    "inclusion-minimal visible edge sets drift",
)


def primitive_projection(data):
    denominators = tuple(sp.denom(value) for value in data["projection"])
    scale = sp.ilcm(*denominators)
    integer_values = tuple(int(scale * value) for value in data["projection"])
    divisor = sp.igcd(*tuple(abs(value) for value in integer_values if value))
    primitive = sp.Matrix(tuple(value // divisor for value in integer_values))
    require(
        data["adjacency"] * primitive == sp.zeros(len(primitive), 1)
        and sum(primitive) > 0,
        "primitive projection direction failed",
    )
    return primitive


z, t = sp.symbols("z t")
EXPECTED_CIRCUITS = {
    ((7, 18),): {
        "mass": sp.Rational(3, 5),
        "signature": (23, 42, 15, 8, 7, 9),
        "support": (
            ((18, "Q+1"), -3), ((18, "Q+2"), 1),
            ((18, "Q+4"), 1), ((18, "Q+5"), 1),
            ((22, "Q+1"), 1), ((22, "Q+2"), 1),
            ((22, "Q+5"), 1),
        ),
        "compatible": tuple(edge for edge in extra_edges if edge != (7, 18)),
        "characteristic": z * (z**2 - 2) * (z**4 - 5*z**2 + 3)
            * (z**8 - 14*z**6 + 62*z**4 - 94*z**2 + 30),
        "odd_polynomial": (t - 2) * (t**2 - 5*t + 3)
            * (t**4 - 14*t**3 + 62*t**2 - 94*t + 30),
        "head": (-108, 19440, -180),
    },
    ((2, 10), (13, 18), (17, 22)): {
        "mass": sp.Rational(12, 35),
        "signature": (23, 52, 13, 7, 6, 11),
        "support": (
            ((2, "Q+4"), -12),
            ((11, "Q+1"), 4), ((11, "Q+2"), 4),
            ((11, "Q+5"), 4),
            ((18, "Q+1"), 3), ((18, "Q+2"), 3),
            ((18, "Q+4"), 3), ((18, "Q+5"), 3),
            ((22, "Q+1"), -4), ((22, "Q+2"), -4),
            ((22, "Q+4"), 12), ((22, "Q+5"), -4),
        ),
        "compatible": ((10, 16), (16, 21)),
        "characteristic": z * (z**2 - 2)
            * (z**10 - 24*z**8 + 198*z**6 - 696*z**4 + 992*z**2 - 420),
        "odd_polynomial": (t - 2)
            * (t**5 - 24*t**4 + 198*t**3 - 696*t**2 + 992*t - 420),
        "head": (288, 241920, 840),
    },
    ((2, 10), (13, 22), (17, 22)): {
        "mass": sp.Rational(12, 35),
        "signature": (23, 46, 13, 7, 6, 11),
        "support": (
            ((2, "Q+4"), -12),
            ((11, "Q+1"), 4), ((11, "Q+2"), 4),
            ((11, "Q+5"), 4),
            ((18, "Q+1"), 3), ((18, "Q+2"), 3),
            ((18, "Q+4"), 3), ((18, "Q+5"), 3),
            ((22, "Q+1"), -4), ((22, "Q+2"), -4),
            ((22, "Q+4"), 12), ((22, "Q+5"), -4),
        ),
        "compatible": ((10, 16), (16, 21)),
        "characteristic": z * (z**2 - 2)
            * (z**10 - 21*z**8 + 158*z**6 - 526*z**4 + 742*z**2 - 315),
        "odd_polynomial": (t - 2)
            * (t**5 - 21*t**4 + 158*t**3 - 526*t**2 + 742*t - 315),
        "head": (216, 136080, 630),
    },
    ((10, 22), (13, 18), (17, 22)): {
        "mass": sp.Rational(12, 35),
        "signature": (23, 52, 13, 7, 6, 11),
        "support": (
            ((2, "Q+4"), 12),
            ((11, "Q+1"), 4), ((11, "Q+2"), 4),
            ((11, "Q+5"), 4),
            ((18, "Q+1"), -3), ((18, "Q+2"), -3),
            ((18, "Q+4"), -3), ((18, "Q+5"), -3),
            ((22, "Q+1"), 4), ((22, "Q+2"), 4),
            ((22, "Q+4"), -12), ((22, "Q+5"), 4),
        ),
        "compatible": ((10, 16), (16, 21)),
        "characteristic": z * (z**2 - 2)
            * (z**10 - 24*z**8 + 199*z**6 - 695*z**4 + 990*z**2 - 420),
        "odd_polynomial": (t - 2)
            * (t**5 - 24*t**4 + 199*t**3 - 695*t**2 + 990*t - 420),
        "head": (288, 241920, 840),
    },
    ((10, 22), (13, 22), (17, 22)): {
        "mass": sp.Rational(12, 35),
        "signature": (23, 46, 13, 7, 6, 11),
        "support": (
            ((2, "Q+4"), 12),
            ((11, "Q+1"), 4), ((11, "Q+2"), 4),
            ((11, "Q+5"), 4),
            ((18, "Q+1"), -3), ((18, "Q+2"), -3),
            ((18, "Q+4"), -3), ((18, "Q+5"), -3),
            ((22, "Q+1"), 4), ((22, "Q+2"), 4),
            ((22, "Q+4"), -12), ((22, "Q+5"), 4),
        ),
        "compatible": ((10, 16), (16, 21)),
        "characteristic": z * (z**2 - 2)
            * (z**10 - 21*z**8 + 158*z**6 - 525*z**4 + 741*z**2 - 315),
        "odd_polynomial": (t - 2)
            * (t**5 - 21*t**4 + 158*t**3 - 525*t**2 + 741*t - 315),
        "head": (216, 136080, 630),
    },
    ((13, 22), (17, 22), (19, 22)): {
        "mass": sp.Rational(3, 7),
        "signature": (23, 46, 13, 7, 6, 11),
        "support": (
            ((2, "Q+4"), 3),
            ((22, "Q+1"), 1), ((22, "Q+2"), 1),
            ((22, "Q+4"), -3), ((22, "Q+5"), 1),
        ),
        "compatible": ((7, 18), (10, 16), (11, 13), (13, 18), (16, 21)),
        "characteristic": z * (z**2 - 2) * (z**4 - 5*z**2 + 3)
            * (z**6 - 16*z**4 + 68*z**2 - 84),
        "odd_polynomial": (t - 2) * (t**2 - 5*t + 3)
            * (t**3 - 16*t**2 + 68*t - 84),
        "head": (216, 108864, 504),
    },
}
require(
    frozenset(EXPECTED_CIRCUITS) == EXPECTED_MINIMAL_VISIBLE_SETS,
    "circuit expectation keys drift",
)


full_nodes, _full_arrows, _full_adjacency = half["rebuilt_graph"](core_edges)
full_index = {node: index for index, node in enumerate(full_nodes)}


def embedded_graph(edges):
    nodes, arrows, adjacency = half["rebuilt_graph"](edges)
    embedded = sp.zeros(len(full_nodes), len(full_nodes))
    for row, node in enumerate(nodes):
        for column, other in enumerate(nodes):
            if adjacency[row, column] != 0:
                embedded[full_index[node], full_index[other]] = adjacency[row, column]
    return frozenset(nodes), arrows, embedded


base_active, base_arrows, base_adjacency = embedded_graph(selected_tree)
circuit_records = []
for chosen, expectation in EXPECTED_CIRCUITS.items():
    data = observer_data(selected_tree + chosen)
    observed_signature = signature(data)
    require(
        observed_signature[:-1] == expectation["signature"]
        and observed_signature[-1] == expectation["mass"],
        ("minimal circuit signature drift", chosen),
    )
    primitive = primitive_projection(data)
    support = tuple(
        (half["node_label"](node), int(value))
        for node, value in zip(data["nodes"], primitive)
        if value != 0
    )
    require(support == expectation["support"], ("circuit support drift", chosen))

    characteristic = half["relation_polynomial"](
        data["order"], data["coefficients"], z
    )
    require(
        sp.factor(characteristic.as_expr() - expectation["characteristic"]) == 0,
        ("circuit characteristic drift", chosen),
    )
    half_transfer = data["adjacency"] * data["adjacency"]
    odd_start = data["adjacency"] * data["all_ones"]
    odd_order, odd_coefficients, _odd_basis = half["krylov_relation"](
        half_transfer, odd_start
    )
    odd_polynomial = half["relation_polynomial"](
        odd_order, odd_coefficients, t
    )
    require(
        sp.factor(odd_polynomial.as_expr() - expectation["odd_polynomial"]) == 0,
        ("circuit odd-half polynomial drift", chosen),
    )
    head = half["matrix_polynomial_on_vector"](
        odd_polynomial, half_transfer, data["all_ones"]
    )
    head_stats = (
        int((data["all_ones"].T * head)[0]),
        int((head.T * head)[0]),
        int(odd_polynomial.nth(0)),
    )
    require(
        head_stats == expectation["head"]
        and data["adjacency"] * head == sp.zeros(len(data["nodes"]), 1)
        and head / odd_polynomial.nth(0) == data["projection"],
        ("circuit recurrence head drift", chosen),
    )

    active, _arrows, adjacency = embedded_graph(selected_tree + chosen)
    primitive_global = sp.zeros(len(full_nodes), 1)
    for node, value in zip(data["nodes"], primitive):
        primitive_global[full_index[node]] = value
    require(
        adjacency * primitive_global == sp.zeros(len(full_nodes), 1),
        ("embedded circuit direction drift", chosen),
    )
    remaining = tuple(edge for edge in extra_edges if edge not in chosen)
    compatible = []
    for edge in remaining:
        _next_active, _next_arrows, next_adjacency = embedded_graph(
            selected_tree + chosen + (edge,)
        )
        if (next_adjacency - adjacency) * primitive_global == sp.zeros(
            len(full_nodes), 1
        ):
            compatible.append(edge)
    require(
        tuple(compatible) == expectation["compatible"],
        ("circuit one-edge compatibility drift", chosen),
    )
    breaking = tuple(edge for edge in remaining if edge not in compatible)
    circuit_records.append((
        chosen, expectation["mass"], observed_signature[:-1], support,
        int(sum(primitive)), int((primitive.T * primitive)[0]),
        tuple(compatible), breaking, sp.factor(characteristic.as_expr()),
        sp.factor(odd_polynomial.as_expr()), head_stats,
    ))


# The special direction is annihilated by every remaining edge delta, hence it
# stays in every later kernel.  Its Rayleigh contribution gives a uniform
# all-ones projection lower bound 3/5 throughout that upper interval.
special_index = extra_edges.index((7, 18))
special_bit = 1 << special_index
require(
    all(
        mass_by_mask[mask] >= sp.Rational(3, 5)
        for mask in range(1 << len(extra_edges))
        if mask & special_bit
    ),
    "persistent-circuit upper-interval mass lower bound failed",
)


# A hostile fourth-edge probe for the first fragile triple.  Three additions
# annihilate its visible kernel projection entirely, while other additions
# preserve or change its mass.  Visibility is therefore not monotone.
fragile = ((2, 10), (13, 18), (17, 22))
EXPECTED_FRAGILE_FOURTH_EDGE_MASSES = {
    (3, 22): sp.Rational(12, 59),
    (7, 18): sp.Rational(14, 15),
    (10, 16): sp.Rational(12, 35),
    (10, 22): sp.Integer(0),
    (11, 13): sp.Rational(12, 37),
    (13, 22): sp.Integer(0),
    (16, 21): sp.Rational(12, 35),
    (19, 22): sp.Integer(0),
}
fragile_fourth_edge_masses = {
    edge: observer_data(selected_tree + fragile + (edge,))["zero_mass"]
    for edge in extra_edges if edge not in fragile
}
require(
    fragile_fourth_edge_masses == EXPECTED_FRAGILE_FOURTH_EDGE_MASSES,
    "fragile fourth-edge mass probe drift",
)


# Universal bordered Woodbury compiler.  Work in the 26-node full-core
# universe and pad the selected tree by isolated future vertices.  If b new
# observer vertices are born and Delta=U C U^T is the added adjacency, then
#
#   G_new(x)-G_tree(x) = b + x s(x)^T (C-x K(x))^{-1} s(x),
#   s=U^T(I-xB)^(-1)w,  K=U^T(I-xB)^(-1)U.
#
# Coefficients require only powers of the fixed base B and a local formal
# inverse of dimension twice the number of new undirected decorated edges.
base_indicator = sp.Matrix([
    1 if node in base_active else 0 for node in full_nodes
])
COMPILER_TERMS = 20
base_sequence = half["sequence"](
    base_adjacency, base_indicator, base_indicator.T, COMPILER_TERMS
)
update_shape_counts = Counter()
max_local_dimension = 0
for first, second, _pair_signature, pair_data in pair_records:
    active, arrows, adjacency = embedded_graph(
        selected_tree + (first, second)
    )
    delta = adjacency - base_adjacency
    undirected_updates = tuple(
        (row, column)
        for row in range(len(full_nodes))
        for column in range(row + 1, len(full_nodes))
        if delta[row, column] != 0
    )
    require(
        all(delta[row, column] == 1 for row, column in undirected_updates)
        and delta == delta.T,
        ("pair adjacency is not a pure edge addition", first, second),
    )
    local_dimension = 2 * len(undirected_updates)
    births = len(active - base_active)
    update_shape_counts[(len(undirected_updates), local_dimension, births)] += 1
    max_local_dimension = max(max_local_dimension, local_dimension)

    u_matrix = sp.zeros(len(full_nodes), local_dimension)
    c_matrix = sp.zeros(local_dimension, local_dimension)
    for update_index, (left, right) in enumerate(undirected_updates):
        first_column = 2 * update_index
        second_column = first_column + 1
        u_matrix[left, first_column] = 1
        u_matrix[right, second_column] = 1
        c_matrix[first_column, second_column] = 1
        c_matrix[second_column, first_column] = 1
    require(
        u_matrix * c_matrix * u_matrix.T == delta
        and c_matrix * c_matrix == sp.eye(local_dimension),
        ("pair local update factorization failed", first, second),
    )

    observer = sp.Matrix([1 if node in active else 0 for node in full_nodes])
    s_coefficients = []
    k_coefficients = []
    observer_power = observer
    u_power = u_matrix
    for _degree in range(COMPILER_TERMS):
        s_coefficients.append(u_matrix.T * observer_power)
        k_coefficients.append(u_matrix.T * u_power)
        observer_power = base_adjacency * observer_power
        u_power = base_adjacency * u_power

    q_coefficients = [c_matrix]
    for degree in range(1, COMPILER_TERMS):
        convolution = sp.zeros(local_dimension, local_dimension)
        for index in range(degree):
            convolution += k_coefficients[index] * q_coefficients[degree - 1 - index]
        q_coefficients.append(c_matrix * convolution)

    corrections = [sp.Integer(0)] * COMPILER_TERMS
    for degree in range(1, COMPILER_TERMS):
        value = sp.Integer(0)
        target = degree - 1
        for left_degree in range(target + 1):
            for middle_degree in range(target - left_degree + 1):
                right_degree = target - left_degree - middle_degree
                value += (
                    s_coefficients[left_degree].T
                    * q_coefficients[middle_degree]
                    * s_coefficients[right_degree]
                )[0]
        corrections[degree] = sp.expand(value)

    compiled = tuple(
        int(base_sequence[degree]
            + (births if degree == 0 else 0)
            + corrections[degree])
        for degree in range(COMPILER_TERMS)
    )
    direct = half["sequence"](
        adjacency, observer, observer.T, COMPILER_TERMS
    )
    require(
        compiled == direct
        and direct[:15] == half["sequence"](
            pair_data["adjacency"], pair_data["all_ones"],
            pair_data["all_ones"].T, 15
        ),
        ("bordered pair compiler mismatch", first, second),
    )

EXPECTED_UPDATE_SHAPE_COUNTS = Counter({
    (2, 4, 0): 21,
    (2, 4, 1): 21,
    (2, 4, 2): 3,
    (5, 10, 0): 7,
    (5, 10, 1): 3,
})
require(
    update_shape_counts == EXPECTED_UPDATE_SHAPE_COUNTS
    and max_local_dimension == 10,
    "bordered-compiler local shape census drift",
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

print("MULTI-EDGE VISIBLE-KERNEL HELLY AND BORDERED-RESOLVENT FINITE-EXACT SCOUT")
print("status=FINITE-EXACT;no_theorem_id;static_symmetric_relation_walks_only")
print("dependency_hashes=BEGIN")
for dependency, expected_hash in DEPENDENCIES.items():
    print(f"{expected_hash}  {dependency.relative_to(ROOT)}")
print("dependency_hashes=END")
print(
    "pair_census=55;visible_pairs=10;new_pair_minimal_circuits=0;"
    "signature_types=17;serialized_bytes=1934;"
    "sha256=86badaa9186021f92c026a591b29436790335d76468db2c18fa260c7feef36f0"
)
print(
    "edge_lattice=2048;serialized_bytes=37894;"
    "sha256=53f740b8fd3e8b1baa88a460e210662d51b6aa69ff8a157fcdfc9fd5d8510ad0"
)
print(f"size_summary={tuple(sorted(size_summary.items()))}")
print("minimal_visible_circuits=6;arity_profile=(one:1,two:0,three:5)")
for record in circuit_records:
    (
        chosen, mass, observed_signature, support, primitive_sum,
        primitive_norm, compatible, breaking, characteristic,
        odd_polynomial, head_stats,
    ) = record
    print(
        f"circuit={chosen};mass={mass};signature={observed_signature};"
        f"primitive=(sum:{primitive_sum},norm:{primitive_norm},"
        f"support_size:{len(support)});compatible_next_edges={compatible};"
        f"breaking_next_edges={breaking};head_stats={head_stats}"
    )
print(
    "persistent_circuit=((7,18),);compatible_with_all_other_edges=true;"
    "upper_interval_zero_mass_lower_bound=3/5"
)
print(f"fragile_triple={fragile};fourth_edge_masses={fragile_fourth_edge_masses}")
print(
    f"bordered_compiler=pairs:55;terms:{COMPILER_TERMS};"
    f"local_dimension_max:{max_local_dimension};"
    f"update_shape_counts={tuple(sorted(update_shape_counts.items()))};exact_match=true"
)
print(
    "mechanism=visible_kernel_Helly_arity_and_local_formal_Woodbury_inverse;"
    "persistent_and_fragile_zero_modes_separate"
)
print(
    "preserved=exact_all_ones_static_walk_counts_active_node_births_and_kernel_mass;"
    "lost=walk_identity_chronology_tournament_orientation_and_bit_complexity"
)
print("scope=no_tournament_no_response_composition_no_GMC_FC_or_LRC_consequence")
print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")
print("ALL EXACT CHECKS PASSED")
