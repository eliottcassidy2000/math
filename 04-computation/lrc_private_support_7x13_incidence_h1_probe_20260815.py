#!/usr/bin/env python3
"""Exact 7-by-13 incidence/H1 sidecar for the THM-3473 support graph.

This is a finite structural probe, not an LRC(14) verifier.  It builds the
bidirected two-section underlying the proved private-sheet atlas, checks
that the reduced incidence matrix has shape 7x13 and full row rank, proves
the unweighted Gram determinant 256 by Cauchy--Binet, and constructs the
six-dimensional cycle space.  The decisive hostile is also explicit: every
owner-potential current lies in the coboundary image and pairs trivially with
all six cycles, so incidence rank alone supplies no nonzero H^1 flux.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from json import dumps
from pathlib import Path
import sys


PACKETS = ((0, 3, 4, 5), (1, 2, 4, 7), (4, 6))
HUB = 4
EXPECTED_EDGES = (
    (0, 3), (0, 4), (0, 5), (1, 2), (1, 4), (1, 7), (2, 4),
    (2, 7), (3, 4), (3, 5), (4, 5), (4, 6), (4, 7),
)
EXPECTED_SEMANTIC_SHA256 = "648523d2d701d19ba35567fd2c96aea65cc4a7ceb65b5e200abfcac8bc8d0a82"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def determinant(matrix: list[list[int | Fraction]]) -> Fraction:
    n = len(matrix)
    require(all(len(row) == n for row in matrix), "square determinant")
    work = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for col in range(n):
        pivot = next((row for row in range(col, n) if work[row][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            result = -result
        value = work[col][col]
        result *= value
        for j in range(col, n):
            work[col][j] /= value
        for row in range(col + 1, n):
            scale = work[row][col]
            if scale:
                for j in range(col, n):
                    work[row][j] -= scale * work[col][j]
    return result


def transpose(matrix: list[list[int | Fraction]]) -> list[list[int | Fraction]]:
    return [list(column) for column in zip(*matrix)]


def multiply(left: list[list[int | Fraction]], right: list[list[int | Fraction]]) -> list[list[Fraction]]:
    require(len(left[0]) == len(right), "matrix product")
    right_t = transpose(right)
    return [
        [sum(Fraction(a) * Fraction(b) for a, b in zip(row, column)) for column in right_t]
        for row in left
    ]


def rref_nullspace(matrix: list[list[int]]) -> tuple[tuple[Fraction, ...], ...]:
    rows = len(matrix)
    cols = len(matrix[0])
    work = [[Fraction(value) for value in row] for row in matrix]
    pivots = []
    pivot_row = 0
    for col in range(cols):
        pivot = next((row for row in range(pivot_row, rows) if work[row][col]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][col]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(rows):
            if row != pivot_row and work[row][col]:
                scale = work[row][col]
                work[row] = [a - scale * b for a, b in zip(work[row], work[pivot_row])]
        pivots.append(col)
        pivot_row += 1
        if pivot_row == rows:
            break
    free = [col for col in range(cols) if col not in pivots]
    basis = []
    for free_col in free:
        vector = [Fraction(0) for _ in range(cols)]
        vector[free_col] = 1
        for row, pivot_col in enumerate(pivots):
            vector[pivot_col] = -work[row][free_col]
        basis.append(tuple(vector))
    return tuple(basis)


def security_report(source: Path) -> tuple[str, ...]:
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert present")
    forbidden = {
        "eval", "exec", "compile", "open", "system", "popen", "run", "Popen",
        "write_text", "write_bytes", "unlink", "remove", "rename",
    }
    called = {
        node.func.id for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    called.update(
        node.func.attr for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    )
    require(not called & forbidden, ("forbidden calls", sorted(called & forbidden)))
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".")[0])
    allowed = {"__future__", "ast", "fractions", "hashlib", "itertools", "json", "pathlib", "sys"}
    require(imports <= allowed, ("unexpected imports", sorted(imports - allowed)))
    return tuple(sorted(imports))


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    edges = set()
    for packet in PACKETS:
        edges.update(combinations(packet, 2))
    edges = tuple(sorted(edges))
    require(edges == EXPECTED_EDGES, edges)
    require(len(edges) == 13, edges)

    full_incidence = [[0 for _edge in edges] for _vertex in range(8)]
    for column, (left, right) in enumerate(edges):
        full_incidence[left][column] = -1
        full_incidence[right][column] = 1
    require(all(sum(full_incidence[row][col] for row in range(8)) == 0 for col in range(13)),
            "column boundary")

    row_vertices = tuple(vertex for vertex in range(8) if vertex != HUB)
    reduced = [full_incidence[vertex] for vertex in row_vertices]
    require((len(reduced), len(reduced[0])) == (7, 13), "7x13 shape")
    gram = multiply(reduced, transpose(reduced))
    gram_det = determinant(gram)
    require(gram_det == 256, gram_det)

    nonzero_minors = []
    sign_hist = {-1: 0, 1: 0}
    cauchy_binet = Fraction(0)
    for columns in combinations(range(13), 7):
        minor = [[row[column] for column in columns] for row in reduced]
        value = determinant(minor)
        cauchy_binet += value * value
        if value:
            require(value in (-1, 1), (columns, value))
            nonzero_minors.append(columns)
            sign_hist[int(value)] += 1
    require(len(nonzero_minors) == 256, len(nonzero_minors))
    require(cauchy_binet == gram_det, (cauchy_binet, gram_det))
    require(256 == 16 * 16, "two K4 factorization")

    cycles = rref_nullspace(reduced)
    require(len(cycles) == 6, len(cycles))
    for cycle in cycles:
        require(
            all(sum(Fraction(reduced[row][col]) * cycle[col] for col in range(13)) == 0
                for row in range(7)),
            ("cycle kernel", cycle),
        )

    # Every owner-potential edge current is B^T f, hence has zero pairing
    # with every cycle.  This is the exact H1 hostile.
    coboundary = transpose(reduced)
    cycle_pairing = multiply([list(cycle) for cycle in cycles], coboundary)
    require(all(value == 0 for row in cycle_pairing for value in row), cycle_pairing)

    # Deleting the hub exposes two rooted K4 blocks and one bridge.  Their
    # unweighted reduced Laplacians have determinants 16,16,1.
    packet_determinants = []
    for packet in PACKETS:
        packet_edges = tuple(combinations(packet, 2))
        root = HUB
        vertices = tuple(vertex for vertex in packet if vertex != root)
        incidence = [[0 for _edge in packet_edges] for _vertex in vertices]
        for column, (left, right) in enumerate(packet_edges):
            for row, vertex in enumerate(vertices):
                incidence[row][column] = int(vertex == right) - int(vertex == left)
        packet_gram = multiply(incidence, transpose(incidence))
        packet_determinants.append(determinant(packet_gram))
    require(tuple(packet_determinants) == (16, 16, 1), packet_determinants)

    semantic_payload = {
        "packets_one_based": tuple(tuple(vertex + 1 for vertex in packet) for packet in PACKETS),
        "edges_one_based": tuple((left + 1, right + 1) for left, right in edges),
        "hub_one_based": HUB + 1,
        "shape": (7, 13),
        "gram": tuple(tuple(int(value) for value in row) for row in gram),
        "gram_det": int(gram_det),
        "nonzero_maximal_minors": len(nonzero_minors),
        "minor_sign_hist": sign_hist,
        "packet_tree_factors": tuple(int(value) for value in packet_determinants),
        "cycle_dimension": len(cycles),
        "cycle_basis": cycles,
        "coboundary_cycle_pairing_zero": True,
        "security": security_report(Path(__file__)),
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":"), default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, semantic_hash)

    print("PRIVATE-SUPPORT 7x13 INCIDENCE / H1 EXACT PROBE")
    print("STATUS: FINITE-EXACT STRUCTURAL SIDECAR TO PROVED THM-3473")
    print(f"PACKETS_ONE_BASED: {semantic_payload['packets_one_based']}")
    print(f"EDGE_COUNT_AND_REDUCED_SHAPE: {(len(edges), (7, 13))}")
    print(f"GRAM_MATRIX: {semantic_payload['gram']}")
    print(f"GRAM_DETERMINANT: {gram_det}")
    print(f"NONZERO_7x7_MINORS_AND_SIGN_HIST: {(len(nonzero_minors), sign_hist)}")
    print(f"TREE_FACTORIZATION_TWO_K4_AND_BRIDGE: {tuple(packet_determinants)}")
    print(f"CYCLE_SPACE_DIMENSION: {len(cycles)}")
    print("H1_HOSTILE: every owner-potential current B^T f pairs zero with all six cycles")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print("VERDICT: canonical full-rank 7x13 incidence carrier exists, but static Boolean currents are coboundaries; a phase/holonomy sidecar is required for nonzero H1 flux")


if __name__ == "__main__":
    main()
