#!/usr/bin/env python3
"""Exact rational series for the THM-3287 maximizing-witness automata.

The decorated vertices are active pairs (response row, reset-link witness).
An arrow (u,s)->(v,t) exists precisely when u,v are adjacent in the selected
row graph and (s,t) belongs to the oriented static dominance relation.  Thus
the coefficient a_n counts witness-decorated walks of length n.  It does not
count chronological response trajectories.
"""

from __future__ import annotations

import ast
import io
import runpy
from contextlib import redirect_stdout
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / (
    "04-computation/gmc_backbone_maximizing_witness_section_scout_20260803.py"
)
SOURCE_OUTPUT = ROOT / (
    "05-knowledge/results/gmc_backbone_maximizing_witness_section_scout_20260803.out"
)
DEPENDENCIES = {
    SOURCE: "aed89d67ec7acabfe5b4feae4a83f7c57b78053928be44a6cc319d81fa4a9cc6",
    SOURCE_OUTPUT: "89200bc6cff7284dd33f636352f4f7f56294d90bcd2902fa96092d9f967f5fe3",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected in DEPENDENCIES.items():
    require(sha256(lf_bytes(dependency)).hexdigest() == expected,
            ("dependency drift", dependency.name))


captured = io.StringIO()
with redirect_stdout(captured):
    source = runpy.run_path(str(SOURCE))
require(captured.getvalue().encode("utf-8") == lf_bytes(SOURCE_OUTPUT),
        "THM-3287 candidate replay differs from stored output")


def lifted_automaton(edges):
    nodes = set()
    arrows = set()
    for left_row, right_row in edges:
        relation = source["oriented_relation"](left_row, right_row)
        for left_state, right_state in relation:
            left = (left_row, left_state)
            right = (right_row, right_state)
            nodes.update((left, right))
            arrows.update(((left, right), (right, left)))
    nodes = tuple(sorted(nodes))
    index = {node: position for position, node in enumerate(nodes)}
    adjacency = [set() for _ in nodes]
    for left, right in arrows:
        adjacency[index[left]].add(index[right])
    adjacency = tuple(tuple(sorted(row)) for row in adjacency)
    require(sum(map(len, adjacency)) == len(arrows), "duplicate automaton arrow")
    require(all(
        left in adjacency[right]
        for right, row in enumerate(adjacency) for left in row
    ), "lifted automaton lost orientation reversal")
    return nodes, tuple(sorted(arrows)), adjacency


def walk_sequence(adjacency, terms):
    vector = [1] * len(adjacency)
    answer = []
    for _ in range(terms):
        answer.append(sum(vector))
        following = [0] * len(adjacency)
        for left, count in enumerate(vector):
            for right in adjacency[left]:
                following[right] += count
        vector = following
    return tuple(answer)


def berlekamp_massey_rational(sequence):
    connection = [Fraction(1)]
    previous = [Fraction(1)]
    length = 0
    shift = 1
    previous_discrepancy = Fraction(1)
    for position in range(len(sequence)):
        discrepancy = Fraction(sequence[position])
        for index in range(1, length + 1):
            discrepancy += connection[index] * sequence[position - index]
        if discrepancy == 0:
            shift += 1
            continue
        saved = connection[:]
        multiplier = -discrepancy / previous_discrepancy
        needed = len(previous) + shift
        if len(connection) < needed:
            connection += [Fraction(0)] * (needed - len(connection))
        for index, coefficient in enumerate(previous):
            connection[index + shift] += multiplier * coefficient
        if 2 * length <= position:
            length = position + 1 - length
            previous = saved
            previous_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1
    coefficients = tuple(-connection[index] for index in range(1, length + 1))
    return length, coefficients


EXPECTED = {
    "backbone": {
        "vertices": 18,
        "arrows": 32,
        "order": 10,
        "coefficients": (0, 16, 0, -94, 0, 244, 0, -263, 0, 93),
        "initial": (18, 32, 84, 160, 412, 798, 2060, 4012, 10426,
                    20362, 53286, 104272, 274628, 538258, 1425644),
    },
    "selected_tree": {
        "vertices": 23,
        "arrows": 40,
        "order": 14,
        "coefficients": (
            0, 20, 0, -158, 0, 630, 0, -1343, 0, 1493, 0, -774, 0, 144,
        ),
        "initial": (23, 40, 104, 212, 536, 1174, 2950, 6700, 16850,
                    39058, 98502, 231232, 584846, 1384502, 3510318),
    },
    "full_core": {
        "vertices": 26,
        "arrows": 68,
        "order": 15,
        "coefficients": (
            0, 34, 0, -404, 0, 2285, 0, -6701, 0, 10096, 0, -7156,
            0, 1808, 0,
        ),
        "initial": (26, 68, 274, 1018, 4134, 16400, 67270, 270718,
                    1114054, 4498562, 18529428, 74885666, 308530344,
                    1247186562, 5138801640),
    },
}


def audit_family(name, edges):
    nodes, arrows, adjacency = lifted_automaton(edges)
    sequence = walk_sequence(adjacency, 4 * len(nodes) + 40)
    order, coefficients = berlekamp_massey_rational(sequence)
    expected = EXPECTED[name]
    require(len(nodes) == expected["vertices"], (name, "vertex census"))
    require(len(arrows) == expected["arrows"], (name, "arrow census"))
    require(order == expected["order"], (name, "recurrence order"))
    require(coefficients == expected["coefficients"],
            (name, "recurrence coefficients"))
    require(sequence[:15] == expected["initial"], (name, "initial terms"))
    for position in range(order, len(sequence)):
        predicted = sum(
            coefficients[lag - 1] * sequence[position - lag]
            for lag in range(1, order + 1)
        )
        require(sequence[position] == predicted,
                (name, "withheld recurrence failure", position))

    hankel = sp.Matrix(
        order, order,
        lambda row, column: sequence[row + column],
    )
    extended_hankel = sp.Matrix(
        order + 1, order + 1,
        lambda row, column: sequence[row + column],
    )
    require(hankel.rank() == order, (name, "minimal Hankel rank"))
    require(extended_hankel.rank() == order, (name, "stable Hankel rank"))

    x, z = sp.symbols("x z")
    denominator = sp.expand(
        1 - sum(
            sp.Rational(coefficient.numerator, coefficient.denominator) * x**lag
            for lag, coefficient in enumerate(coefficients, 1)
        )
    )
    numerator_coefficients = []
    for position in range(order):
        coefficient = Fraction(sequence[position])
        for lag in range(1, min(position, order) + 1):
            coefficient -= coefficients[lag - 1] * sequence[position - lag]
        numerator_coefficients.append(coefficient)
    numerator = sp.expand(sum(
        sp.Rational(coefficient.numerator, coefficient.denominator) * x**position
        for position, coefficient in enumerate(numerator_coefficients)
    ))
    characteristic = sp.expand(
        z**order - sum(
            sp.Rational(coefficient.numerator, coefficient.denominator)
            * z ** (order - lag)
            for lag, coefficient in enumerate(coefficients, 1)
        )
    )
    require(sp.Poly(denominator, x).LC() != 0, (name, "zero denominator"))

    # Coefficientwise multiplication by the denominator recovers the printed
    # numerator and vanishes thereafter, including the full-core n=14 atom.
    product_coefficients = []
    denominator_coefficients = sp.Poly(denominator, x).all_coeffs()[::-1]
    denominator_degree = sp.degree(denominator, x)
    for position in range(len(sequence)):
        value = sum(
            denominator_coefficients[lag] * sequence[position - lag]
            for lag in range(min(position, denominator_degree) + 1)
        )
        product_coefficients.append(value)
    numerator_poly = sp.Poly(numerator, x)
    require(all(
        product_coefficients[position] == numerator_poly.nth(position)
        for position in range(sp.degree(numerator, x) + 1)
    ), (name, "generating numerator mismatch"))
    require(all(
        value == 0
        for value in product_coefficients[sp.degree(numerator, x) + 1:]
    ), (name, "generating tail mismatch"))

    print(
        f"family={name};active_vertices={len(nodes)};directed_arrows={len(arrows)};"
        f"minimal_scalar_order={order};hankel_rank={order}"
    )
    print(f"family={name};a_0_to_a_14={sequence[:15]}")
    print(f"family={name};recurrence_coefficients_lag1_to_{order}={coefficients}")
    print(
        f"family={name};characteristic_factor="
        f"{sp.factor(characteristic)}"
    )
    print(
        f"family={name};generating_denominator_factor="
        f"{sp.factor(denominator)}"
    )
    print(f"family={name};generating_numerator={sp.factor(numerator)}")


audit_family("backbone", source["backbone"])
audit_family("selected_tree", source["selected_tree"])
audit_family("full_core", source["core_edges"])


syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assert_nodes == 0 and float_literals == 0, "source AST gate")

print("MAXIMIZING-WITNESS LIFTED-WALK RATIONAL SERIES")
print(f"dependency_hash_checks={len(DEPENDENCIES)};source_replay=PASS")
print("definition=a_n_counts_length_n_walks_in_active_(row,static_witness)_automaton")
print("mechanism=finite_symmetric_bipartite_transfer_matrix;parity_recurrences_are_exact")
print("scope=static_dominance_witness_sequences_only;not_simple_row_paths_not_chronological_response_dynamics")
print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")
print("all_exact_checks=PASS")
