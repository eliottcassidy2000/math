#!/usr/bin/env python3
"""Exact aggregate-circuit reduction for static decorated relation walks.

The earlier multi-edge scouts enumerated the full 2^11 lattice of missing
row-edge choices and found six inclusion-minimal observer-visible zero modes.
This script gives a structural derivation.  The selector bipartition writes
the adjacency as A=(0,N;N^T,0).  Twelve column equations N^T xi=0 collapse to
nine aggregate equations; a two-case Boolean analysis then produces the six
minimal edge supports without enumerating decorated graphs.

The previous independent 2^11 audit is pinned but not executed.  It remains
an independent verification of the structural derivation here.  Everything
in this file concerns finite static symmetric relation-walk graphs.  There is
no tournament, chronology, response-composition, FC/GMC, LRC, asymptotic, or
bit-complexity conclusion.
"""

from __future__ import annotations

import ast
import hashlib
import io
import runpy
from contextlib import redirect_stdout
from itertools import product
from pathlib import Path

import sympy as sp


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
HALF_SOURCE = ROOT / (
    "04-computation/gmc_selector_paired_half_transfer_partial_scout_20260803.py"
)
HALF_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "gmc_selector_paired_half_transfer_partial_scout_20260803.out"
)
THM3305 = ROOT / (
    "01-canon/theorems/"
    "THM-3305-single-edge-visible-zero-mode-and-rank-two-resolvent-update.md"
)
AUDIT_SOURCE = ROOT / (
    "04-computation/"
    "gmc_multi_edge_visible_kernel_bordered_resolvent_independent_audit_20260803.py"
)
AUDIT_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "gmc_multi_edge_visible_kernel_bordered_resolvent_independent_audit_20260803.out"
)
DEPENDENCIES = {
    HALF_SOURCE:
        "da47f9db80558343bd6442196c65150417241acbed65385b876ca09a66f0490f",
    HALF_OUTPUT:
        "ef276b6a99a3ec1ac59d5cc352a3aa9a50fa64440a6615695d684a68a6683208",
    THM3305:
        "fd9c4feaa9a39db14c754a3ff7267ad5f71e2cd71953119c07df53fafb96bab6",
    AUDIT_SOURCE:
        "75a7afacc5056976efff5686b3afe44534e8c64f1e2a9fee18d26b1c70f03e52",
    AUDIT_OUTPUT:
        "ce7f9239b9ba68813ceb919662458b1c2ec32dea98d41507edab341a5417f651",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected_hash in DEPENDENCIES.items():
    require(
        hashlib.sha256(lf_bytes(dependency)).hexdigest() == expected_hash,
        ("dependency or independent-audit hash drift", dependency.name),
    )

# Replay only the older half-transfer source.  The 2^11 audit is a hash-pinned
# comparison target, not an executable dependency of this structural proof.
captured = io.StringIO()
with redirect_stdout(captured):
    half = runpy.run_path(str(HALF_SOURCE))
require(
    captured.getvalue().encode("utf-8") == lf_bytes(HALF_OUTPUT),
    "selector-paired half-transfer replay differs from stored output",
)

relation = half["relation_source"]
selected_tree = tuple(relation["selected_tree"])
core_edges = tuple(relation["core_edges"])
selected_set = frozenset(selected_tree)
extra_edges = tuple(edge for edge in core_edges if edge not in selected_set)
require(
    len(selected_tree) == len(extra_edges) == 11 and len(core_edges) == 22,
    "selected-tree/full-core edge census drift",
)


# Fix the full-core selector universe.  Every small-side node is already born
# in the selected tree.  Three full-side Q-8 nodes are born by optional edges.
full_core_nodes, _full_core_arrows, full_core_adjacency = half["rebuilt_graph"](
    core_edges
)
small_nodes = tuple(
    node for node in full_core_nodes if node[0] in half["selector_small"]
)
full_nodes = tuple(
    node for node in full_core_nodes if node[0] in half["selector_full"]
)
small_labels = tuple(half["node_label"](node) for node in small_nodes)
full_labels = tuple(half["node_label"](node) for node in full_nodes)
EXPECTED_SMALL_LABELS = (
    (2, "Q+4"),
    (11, "Q+1"), (11, "Q+2"), (11, "Q+4"), (11, "Q+5"),
    (16, "Q+1"),
    (18, "Q+1"), (18, "Q+2"), (18, "Q+4"), (18, "Q+5"),
    (22, "Q+1"), (22, "Q+2"), (22, "Q+4"), (22, "Q+5"),
)
EXPECTED_FULL_LABELS = (
    (3, "Q-8"), (3, "Q-1"), (7, "Q-7"),
    (10, "Q-8"), (10, "Q-1"),
    (13, "Q-8"), (13, "Q-1"),
    (17, "Q-8"), (17, "Q-1"), (19, "Q-1"),
    (21, "Q-8"), (21, "Q-1"),
)
require(
    small_labels == EXPECTED_SMALL_LABELS
    and full_labels == EXPECTED_FULL_LABELS,
    "full-core selector-node ordering drift",
)
small_index = {node: index for index, node in enumerate(small_nodes)}
full_index = {node: index for index, node in enumerate(full_nodes)}


def selector_block(edges):
    """Return the fixed-universe small-to-full block N for these row edges."""
    nodes, _arrows, adjacency = half["rebuilt_graph"](edges)
    local_index = {node: index for index, node in enumerate(nodes)}
    block = sp.zeros(len(small_nodes), len(full_nodes))
    for small_node, row in small_index.items():
        if small_node not in local_index:
            continue
        for full_node, column in full_index.items():
            if full_node in local_index:
                block[row, column] = adjacency[
                    local_index[small_node], local_index[full_node]
                ]
    return block


base_block = selector_block(selected_tree)
full_core_index = {
    node: index for index, node in enumerate(full_core_nodes)
}
core_block = selector_block(core_edges)
require(
    all(
        full_core_adjacency[
            full_core_index[left], full_core_index[right]
        ] == 0
        for side in (small_nodes, full_nodes)
        for left in side for right in side
    )
    and all(
        full_core_adjacency[
            full_core_index[left], full_core_index[right]
        ] == core_block[small_index[left], full_index[right]]
        for left in small_nodes for right in full_nodes
    ),
    "full-core adjacency does not have the selector block form",
)


# Each indicator records one missing row edge.  The letters are chosen to
# keep the aggregate equations legible; e,h,j,k will disappear from the
# visibility criterion but are retained in the exact block identity.
p, e, w, h, q, j, r, s, k, u, v = sp.symbols(
    "p e w h q j r s k u v"
)
indicator_by_edge = {
    (2, 10): p,
    (3, 22): e,
    (7, 18): w,
    (10, 16): h,
    (10, 22): q,
    (11, 13): j,
    (13, 18): r,
    (13, 22): s,
    (16, 21): k,
    (17, 22): u,
    (19, 22): v,
}
require(
    tuple(indicator_by_edge) == extra_edges,
    "optional-edge indicator order drift",
)
symbolic_block = base_block.copy()
for edge in extra_edges:
    delta = selector_block(selected_tree + (edge,)) - base_block
    require(
        all(value in (0, 1) for value in delta),
        ("optional relation is not a pure decorated-edge addition", edge),
    )
    symbolic_block += indicator_by_edge[edge] * delta


# Name the fourteen small-side coefficients and their row aggregates.
x, a1, a2, a4, a5, b, c1, c2, c4, c5, d1, d2, d4, d5 = sp.symbols(
    "x a1 a2 a4 a5 b c1 c2 c4 c5 d1 d2 d4 d5"
)
xi = sp.Matrix((x, a1, a2, a4, a5, b, c1, c2, c4, c5, d1, d2, d4, d5))
A_sum = a1 + a2 + a4 + a5
C_sum = c1 + c2 + c4 + c5
D_sum = d1 + d2 + d4 + d5
column_equations = symbolic_block.T * xi
EXPECTED_COLUMN_EQUATIONS = sp.Matrix((
    b,
    a4 + e*d4,
    w*c1 + D_sum,
    h*b,
    A_sum + p*x + q*d4,
    j*a1,
    x + r*C_sum + s*d4,
    b,
    x + u*d4,
    x + C_sum + v*d4,
    k*b,
    x + d4,
))
require(
    sp.simplify(column_equations - EXPECTED_COLUMN_EQUATIONS)
        == sp.zeros(len(full_nodes), 1),
    "twelve exact selector-column equations drift",
)


# The full-side kernel is universally invisible.  If E=N z denotes the small
# row equations, the active Q-8 observer and the active non-Q-8 observer are
# explicit row combinations.  Optional h,j,k create the three future Q-8
# nodes, so their observer coefficients are the same edge indicators.
z_variables = sp.symbols("z0:12")
z_vector = sp.Matrix(z_variables)
row_equations = symbolic_block * z_vector
row_by_label = {
    label: row_equations[index] for index, label in enumerate(small_labels)
}
active_full_observer = sp.Matrix((1, 1, 1, h, 1, j, 1, 1, 1, 1, k, 1))
q8_observer = (
    row_by_label[(16, "Q+1")]
    + row_by_label[(11, "Q+1")]
    - row_by_label[(11, "Q+2")]
)
non_q8_observer = (
    row_by_label[(11, "Q+4")]
    + row_by_label[(22, "Q+1")]
    + row_by_label[(2, "Q+4")]
    - p*row_by_label[(11, "Q+2")]
)
require(
    sp.expand(
        (active_full_observer.T * z_vector)[0]
        - q8_observer - non_q8_observer
    ) == 0,
    "full-side active observer is not in the row span of N",
)


# After the permanent columns force b=0 and d4=-x, the visibility-relevant
# system is exactly the following nine scalar constraints:
#
#   b=0, j a1=0, a4=e x, w c1+D=0,
#   A+(p-q)x=0, (1-s)x+rC=0,
#   (1-u)x=0, C+(1-v)x=0.
#
# The ninth datum is the observer sum S=x+A+C+D.  The lift formulas below
# show that every aggregate solution is realized by actual node coefficients.
X, AA, CC, DD = sp.symbols("X AA CC DD")
lift_w1 = sp.Matrix((
    X, 0, AA-e*X, e*X, 0, 0,
    -DD, CC+DD, 0, 0,
    DD+X, 0, -X, 0,
))
lift_w0 = sp.Matrix((
    X, 0, AA-e*X, e*X, 0, 0,
    0, CC, 0, 0,
    X, 0, -X, 0,
))
residual_w1 = sp.simplify((symbolic_block.T * lift_w1).subs(w, 1))
residual_w0 = sp.simplify((symbolic_block.T * lift_w0).subs(w, 0))
EXPECTED_RESIDUAL_W1 = sp.Matrix((
    0, 0, 0, 0,
    AA + (p-q)*X,
    0,
    (1-s)*X + r*CC,
    0,
    (1-u)*X,
    CC + (1-v)*X,
    0, 0,
))
EXPECTED_RESIDUAL_W0 = EXPECTED_RESIDUAL_W1
require(
    sp.simplify(residual_w1 - EXPECTED_RESIDUAL_W1)
        == sp.zeros(len(full_nodes), 1)
    and sp.simplify(residual_w0 - EXPECTED_RESIDUAL_W0)
        == sp.zeros(len(full_nodes), 1),
    "aggregate-to-node lifting formulas drift",
)
require(
    sp.expand(sum(lift_w1) - (X + AA + CC + DD)) == 0
    and sp.expand(sum(lift_w0) - (X + AA + CC)) == 0,
    "lifted observer sum differs from aggregate observer sum",
)


# Structural case split.
#
# If w=1, take X=AA=CC=0 and DD != 0: the singleton w is visible and every
# superset remains visible.  If w=0, then DD=0.  A visible solution must have
# X != 0, and division by X leaves only
#
#   u=1,  1-s=r(1-v),  q-p+v != 0.
#
# We use the 2^6 reduced truth table only as a hostile check of the written
# r=0/r=1 case split; no decorated graph or 2^11 lattice is enumerated.
relevant_symbols = (p, q, r, s, u, v)
symbol_to_edge = {
    p: (2, 10), q: (10, 22), r: (13, 18),
    s: (13, 22), u: (17, 22), v: (19, 22),
}
visible_reduced_supports = []
for bits in product((0, 1), repeat=len(relevant_symbols)):
    values = dict(zip(relevant_symbols, bits))
    visible = (
        values[u] == 1
        and 1-values[s] == values[r]*(1-values[v])
        and values[q]-values[p]+values[v] != 0
    )
    if visible:
        visible_reduced_supports.append(frozenset(
            symbol for symbol, bit in values.items() if bit
        ))

minimal_reduced_supports = tuple(sorted(
    (
        support for support in visible_reduced_supports
        if not any(other < support for other in visible_reduced_supports)
    ),
    key=lambda support: tuple(str(symbol) for symbol in sorted(
        support, key=str
    )),
))
minimal_reduced_serialized = tuple(
    tuple(sorted(str(symbol) for symbol in support))
    for support in minimal_reduced_supports
)
derived_minimal_sets = frozenset({((7, 18),)} | {
    tuple(sorted(symbol_to_edge[symbol] for symbol in support))
    for support in minimal_reduced_supports
})
EXPECTED_MINIMAL_SETS = frozenset({
    ((7, 18),),
    ((2, 10), (13, 18), (17, 22)),
    ((2, 10), (13, 22), (17, 22)),
    ((10, 22), (13, 18), (17, 22)),
    ((10, 22), (13, 22), (17, 22)),
    ((13, 22), (17, 22), (19, 22)),
})
require(
    derived_minimal_sets == EXPECTED_MINIMAL_SETS,
    "aggregate case split does not yield the six audited minimal circuits",
)


# For just the six structurally derived circuits, project the small-side
# all-ones vector onto ker(N^T).  The full-side projection is zero by the row
# identity above.  This recovers the primitive directions and exact masses
# found by the independent exhaustive audit.
EXPECTED_CIRCUITS = {
    ((7, 18),): (
        sp.Rational(3, 5),
        (
            ((18, "Q+1"), -3), ((18, "Q+2"), 1),
            ((18, "Q+4"), 1), ((18, "Q+5"), 1),
            ((22, "Q+1"), 1), ((22, "Q+2"), 1),
            ((22, "Q+5"), 1),
        ),
    ),
    ((2, 10), (13, 18), (17, 22)): (
        sp.Rational(12, 35),
        (
            ((2, "Q+4"), -12),
            ((11, "Q+1"), 4), ((11, "Q+2"), 4),
            ((11, "Q+5"), 4),
            ((18, "Q+1"), 3), ((18, "Q+2"), 3),
            ((18, "Q+4"), 3), ((18, "Q+5"), 3),
            ((22, "Q+1"), -4), ((22, "Q+2"), -4),
            ((22, "Q+4"), 12), ((22, "Q+5"), -4),
        ),
    ),
    ((2, 10), (13, 22), (17, 22)): (
        sp.Rational(12, 35),
        (
            ((2, "Q+4"), -12),
            ((11, "Q+1"), 4), ((11, "Q+2"), 4),
            ((11, "Q+5"), 4),
            ((18, "Q+1"), 3), ((18, "Q+2"), 3),
            ((18, "Q+4"), 3), ((18, "Q+5"), 3),
            ((22, "Q+1"), -4), ((22, "Q+2"), -4),
            ((22, "Q+4"), 12), ((22, "Q+5"), -4),
        ),
    ),
    ((10, 22), (13, 18), (17, 22)): (
        sp.Rational(12, 35),
        (
            ((2, "Q+4"), 12),
            ((11, "Q+1"), 4), ((11, "Q+2"), 4),
            ((11, "Q+5"), 4),
            ((18, "Q+1"), -3), ((18, "Q+2"), -3),
            ((18, "Q+4"), -3), ((18, "Q+5"), -3),
            ((22, "Q+1"), 4), ((22, "Q+2"), 4),
            ((22, "Q+4"), -12), ((22, "Q+5"), 4),
        ),
    ),
    ((10, 22), (13, 22), (17, 22)): (
        sp.Rational(12, 35),
        (
            ((2, "Q+4"), 12),
            ((11, "Q+1"), 4), ((11, "Q+2"), 4),
            ((11, "Q+5"), 4),
            ((18, "Q+1"), -3), ((18, "Q+2"), -3),
            ((18, "Q+4"), -3), ((18, "Q+5"), -3),
            ((22, "Q+1"), 4), ((22, "Q+2"), 4),
            ((22, "Q+4"), -12), ((22, "Q+5"), 4),
        ),
    ),
    ((13, 22), (17, 22), (19, 22)): (
        sp.Rational(3, 7),
        (
            ((2, "Q+4"), 3),
            ((22, "Q+1"), 1), ((22, "Q+2"), 1),
            ((22, "Q+4"), -3), ((22, "Q+5"), 1),
        ),
    ),
}


def primitive_direction(vector):
    denominators = tuple(sp.denom(value) for value in vector)
    scale = sp.ilcm(*denominators)
    integers = tuple(int(scale*value) for value in vector)
    divisor = sp.igcd(*tuple(abs(value) for value in integers if value))
    primitive = sp.Matrix(tuple(value // divisor for value in integers))
    if sum(primitive) < 0:
        primitive = -primitive
    return primitive


circuit_records = []
primitive_by_circuit = {}
for chosen in sorted(EXPECTED_CIRCUITS):
    block = selector_block(selected_tree + chosen)
    kernel_basis = block.T.nullspace()
    require(kernel_basis, ("small-side kernel vanished", chosen))
    basis = sp.Matrix.hstack(*kernel_basis)
    observer = sp.ones(len(small_nodes), 1)
    projection = basis * (basis.T*basis).inv() * basis.T * observer
    require(
        block.T*projection == sp.zeros(len(full_nodes), 1)
        and projection != sp.zeros(len(small_nodes), 1),
        ("derived minimal circuit is not observer-visible", chosen),
    )
    mass = sp.factor(sum(projection))
    require(
        sp.factor((projection.T*projection)[0]) == mass,
        ("orthogonal projection mass identity failed", chosen),
    )
    primitive = primitive_direction(projection)
    support = tuple(
        (label, int(value))
        for label, value in zip(small_labels, primitive) if value != 0
    )
    expected_mass, expected_support = EXPECTED_CIRCUITS[chosen]
    require(
        mass == expected_mass and support == expected_support,
        ("circuit mass or primitive support differs from prior audit", chosen),
    )
    primitive_by_circuit[chosen] = primitive
    circuit_records.append((chosen, mass, support, int(sum(primitive)),
                            int((primitive.T*primitive)[0])))

# The singleton's audited primitive direction is stronger than mere
# visibility: every other optional-edge block delta annihilates it.  It
# therefore remains a fixed kernel direction on all 2^10 supersets, and its
# Rayleigh contribution gives the uniform mass lower bound 3/5 without a
# superset census.
singleton = ((7, 18),)
singleton_primitive = primitive_by_circuit[singleton]
require(
    all(
        (
            selector_block(selected_tree + (edge,)) - base_block
        ).T * singleton_primitive == sp.zeros(len(full_nodes), 1)
        for edge in extra_edges if edge != singleton[0]
    )
    and sp.Rational(
        int(sum(singleton_primitive))**2,
        int((singleton_primitive.T*singleton_primitive)[0]),
    ) == sp.Rational(3, 5),
    "singleton circuit is not a persistent mass-3/5 direction",
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


print("VISIBLE-KERNEL AGGREGATE CIRCUIT REDUCTION FINITE-EXACT AUDIT")
print("status=FINITE-EXACT;structural_derivation_not_full_lattice_enumeration")
print("dependency_hashes=BEGIN")
for dependency, expected_hash in DEPENDENCIES.items():
    print(f"{expected_hash}  {dependency.relative_to(ROOT)}")
print("dependency_hashes=END")
print(
    f"selector_block_shape={symbolic_block.shape};"
    f"small_labels={small_labels};full_labels={full_labels}"
)
print(
    "column_reduction="
    "b=0;j*a1=0;a4=e*x;w*c1+D=0;A+(p-q)*x=0;"
    "(1-s)*x+r*C=0;(1-u)*x=0;C+(1-v)*x=0;"
    "observer_sum=x+A+C+D"
)
print(
    "full_side_invisibility="
    "Q8=E16+E11Q1-E11Q2;"
    "nonQ8=E11Q4+E22Q1+E2-p*E11Q2;PASS"
)
print(
    "case_w=1=singleton_persistent_via_x=A=C=0,D_nonzero;"
    "case_w=0=requires_x_nonzero,u=1,1-s=r*(1-v),q-p+v_nonzero"
)
print(
    "singleton_fixed_direction=annihilated_by_all_10_remaining_edge_deltas;"
    "all_superset_zero_masses_at_least_3/5"
)
print(
    f"reduced_boolean_hostile_cases={2**len(relevant_symbols)};"
    f"visible_cases={len(visible_reduced_supports)};"
    f"minimal_reduced_supports={minimal_reduced_serialized}"
)
for chosen, mass, support, primitive_sum, primitive_norm in circuit_records:
    print(
        f"circuit={chosen};mass={mass};sum={primitive_sum};"
        f"norm={primitive_norm};support={support}"
    )
print(
    "mechanism=origin_in_affine_hull_of_active_small_row_neighborhood_vectors;"
    "minimal_optional_edge_support_is_an_edge_parameter_origin_affine_circuit"
)
print(
    "preserved=exact_static_zero_visibility_and_primitive_projection_direction;"
    "lost=walk_chronology_response_dynamics_and_general_graph_asymptotics"
)
print(
    "scope=static_symmetric_relation_walks_only;"
    "no_tournament_no_FC_GMC_or_LRC_consequence_no_bit_complexity_claim"
)
print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")
print("ALL EXACT CHECKS PASSED")
