#!/usr/bin/env python3
"""Exact companion for independently audited THM-3482.

The owner order is the proved THM-3473 order.  Every coactivity edge is
oriented from the smaller to the larger owner index.  The private-sheet count
vector is used as an owner potential, and its coboundary supplies thirteen
signed edge weights.  The program derives the two symbolic K4 tree sums,
checks the weighted 7x7 reduced-Laplacian determinant directly, and verifies
that the same spectrally nonzero edge vector has zero pairing with every graph
cycle.  It is not an LRC current realization or an LRC(14) verifier.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from json import dumps
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = (
    (
        "THM-3473",
        "01-canon/theorems/THM-3473-three-times-p-eight-owner-private-sheet-partition-and-irredundancy.md",
        "05236a7338a4c92b443b2798508d407e2e4f82d942111d022064f6ec9fe86ca2",
    ),
    (
        "7x13-incidence-sidecar",
        "04-computation/lrc_private_support_7x13_incidence_h1_probe_20260815.py",
        "6efa87aa9f9b50d57d7a2db3c282ad216057b41c10af08f648e3f3398e457b91",
    ),
)
PACKET_A = (0, 3, 4, 5)
PACKET_B = (1, 2, 4, 7)
BRIDGE = (4, 6)
PACKETS = (PACKET_A, PACKET_B, BRIDGE)
HUB = 4
EXPECTED_EDGES = (
    (0, 3), (0, 4), (0, 5), (1, 2), (1, 4), (1, 7), (2, 4),
    (2, 7), (3, 4), (3, 5), (4, 5), (4, 6), (4, 7),
)
EXPECTED_TAU_A = {
    0: (0, 0, 0, -64),
    1: (0, 0, 0, -64),
    2: (8, -48, 96, -64),
}
EXPECTED_TAU_B = {
    0: (-8, 16, 32, -64),
    1: (8, 16, -96, -64),
    2: (8, -48, 96, -64),
}
EXPECTED_BRIDGE = {0: (0, -4), 1: (0, -4), 2: (0, -4)}
EXPECTED_ZERO_EDGES = {0: 4, 1: 4, 2: 6}
EXPECTED_SEMANTIC_SHA256 = "98cba04620048b6c6f8fab03518ab39e0125623f033423eda29245c9b53a0162"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def trim(poly: tuple[int, ...]) -> tuple[int, ...]:
    values = list(poly)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def poly_add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    size = max(len(left), len(right))
    return trim(tuple(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ))


def poly_neg(poly: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-value for value in poly)


def poly_sub(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return poly_add(left, poly_neg(right))


def poly_mul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    values = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            values[i + j] += a * b
    return trim(tuple(values))


def poly_eval(poly: tuple[int, ...], value: int) -> int:
    total = 0
    for coefficient in reversed(poly):
        total = total * value + coefficient
    return total


def lf_digest(relative: str) -> str:
    data = (ROOT / relative).read_bytes().replace(b"\r\n", b"\n")
    return sha256(data).hexdigest()


def dependency_hashes() -> tuple[tuple[str, str, str], ...]:
    rows = []
    for label, relative, expected in DEPENDENCIES:
        observed = lf_digest(relative)
        require(observed == expected, (label, observed, expected))
        rows.append((label, relative, observed))
    return tuple(rows)


def private_count_polynomials(state: int) -> tuple[tuple[int, ...], ...]:
    e0 = int(state == 0)
    e1 = int(state == 1)
    e2 = int(state == 2)
    return (
        (0, 4),
        (-2 * e1, 4),
        (-2 * e0, 4),
        (0, 4),
        (-2 * e2, 8),
        (0, 4),
        (-2 * e2, 4),
        (0, 4),
    )


def private_counts(k: int) -> tuple[int, ...]:
    return tuple(poly_eval(poly, k) for poly in private_count_polynomials(k % 3))


def packet_edges(packet: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(combinations(packet, 2))


def is_tree(packet: tuple[int, ...], edges: tuple[tuple[int, int], ...]) -> bool:
    if len(edges) != len(packet) - 1:
        return False
    adjacency = {vertex: set() for vertex in packet}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    reached = {packet[0]}
    stack = [packet[0]]
    while stack:
        vertex = stack.pop()
        for neighbour in adjacency[vertex]:
            if neighbour not in reached:
                reached.add(neighbour)
                stack.append(neighbour)
    return len(reached) == len(packet)


def tree_polynomial(
    packet: tuple[int, ...], potentials: tuple[tuple[int, ...], ...]
) -> tuple[int, ...]:
    total = (0,)
    edges = packet_edges(packet)
    for chosen in combinations(edges, len(packet) - 1):
        if not is_tree(packet, chosen):
            continue
        product = (1,)
        for left, right in chosen:
            product = poly_mul(product, poly_sub(potentials[right], potentials[left]))
        total = poly_add(total, product)
    return total


def determinant(matrix: list[list[int]]) -> int:
    size = len(matrix)
    require(size and all(len(row) == size for row in matrix), "square matrix")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        pivot_row = next(
            (row for row in range(pivot_index, size) if work[row][pivot_index]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = work[row][column] * pivot - work[row][pivot_index] * work[pivot_index][column]
                require(numerator % previous == 0, (pivot_index, numerator, previous))
                work[row][column] = numerator // previous
        previous = pivot
        for row in range(pivot_index + 1, size):
            work[row][pivot_index] = 0
    return sign * work[-1][-1]


def full_incidence(edges: tuple[tuple[int, int], ...]) -> list[list[int]]:
    matrix = [[0 for _edge in edges] for _vertex in range(8)]
    for column, (left, right) in enumerate(edges):
        matrix[left][column] = -1
        matrix[right][column] = 1
    return matrix


def reduced_incidence(edges: tuple[tuple[int, int], ...]) -> list[list[int]]:
    full = full_incidence(edges)
    return [row for vertex, row in enumerate(full) if vertex != HUB]


def edge_weights(potentials: tuple[int, ...], edges: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    return tuple(potentials[right] - potentials[left] for left, right in edges)


def weighted_laplacian_determinant(potentials: tuple[int, ...]) -> int:
    incidence = reduced_incidence(EXPECTED_EDGES)
    weights = edge_weights(potentials, EXPECTED_EDGES)
    matrix = [
        [
            sum(incidence[row][edge] * weights[edge] * incidence[column][edge]
                for edge in range(len(EXPECTED_EDGES)))
            for column in range(7)
        ]
        for row in range(7)
    ]
    return determinant(matrix)


def rref_nullspace(matrix: list[list[int]]) -> tuple[tuple[Fraction, ...], ...]:
    rows = len(matrix)
    columns = len(matrix[0])
    work = [[Fraction(value) for value in row] for row in matrix]
    pivots = []
    pivot_row = 0
    for column in range(columns):
        selected = next((row for row in range(pivot_row, rows) if work[row][column]), None)
        if selected is None:
            continue
        work[pivot_row], work[selected] = work[selected], work[pivot_row]
        value = work[pivot_row][column]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(rows):
            if row != pivot_row and work[row][column]:
                scale = work[row][column]
                work[row] = [a - scale * b for a, b in zip(work[row], work[pivot_row])]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == rows:
            break
    free = [column for column in range(columns) if column not in pivots]
    basis = []
    for free_column in free:
        vector = [Fraction(0) for _ in range(columns)]
        vector[free_column] = 1
        for row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -work[row][free_column]
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

    edges = tuple(sorted({edge for packet in PACKETS for edge in packet_edges(packet)}))
    require(edges == EXPECTED_EDGES, edges)
    dependencies = dependency_hashes()
    incidence = reduced_incidence(edges)
    cycles = rref_nullspace(incidence)
    require(len(cycles) == 6, len(cycles))

    symbolic_rows = []
    for state in range(3):
        potentials = private_count_polynomials(state)
        tau_a = tree_polynomial(PACKET_A, potentials)
        tau_b = tree_polynomial(PACKET_B, potentials)
        bridge = poly_sub(potentials[BRIDGE[1]], potentials[BRIDGE[0]])
        require(tau_a == EXPECTED_TAU_A[state], (state, "A", tau_a))
        require(tau_b == EXPECTED_TAU_B[state], (state, "B", tau_b))
        require(bridge == EXPECTED_BRIDGE[state], (state, "bridge", bridge))
        full = poly_mul(poly_mul(tau_a, tau_b), bridge)

        symbolic_weights = tuple(
            poly_sub(potentials[right], potentials[left]) for left, right in edges
        )
        for cycle in cycles:
            pairing = (0,)
            for coefficient, weight in zip(cycle, symbolic_weights):
                require(coefficient.denominator == 1, (state, coefficient))
                term = tuple(coefficient.numerator * value for value in weight)
                pairing = poly_add(pairing, term)
            require(pairing == (0,), (state, cycle, pairing))

        zero_edges = sum(weight == (0,) for weight in symbolic_weights)
        require(zero_edges == EXPECTED_ZERO_EDGES[state], (state, zero_edges))
        symbolic_rows.append((state, tau_a, tau_b, bridge, full, zero_edges))

    direct_rows = []
    for k in range(1, 3001):
        state = k % 3
        potentials = private_counts(k)
        weights = edge_weights(potentials, edges)
        tau_a = poly_eval(EXPECTED_TAU_A[state], k)
        tau_b = poly_eval(EXPECTED_TAU_B[state], k)
        bridge = poly_eval(EXPECTED_BRIDGE[state], k)
        factored = tau_a * tau_b * bridge
        direct = weighted_laplacian_determinant(potentials)
        require(direct == factored, (k, direct, factored))
        require(tau_a < 0 and tau_b < 0 and bridge < 0 and direct < 0,
                (k, tau_a, tau_b, bridge, direct))
        require(sum(weight == 0 for weight in weights) == EXPECTED_ZERO_EDGES[state],
                (k, weights))
        direct_rows.append((k, state, tau_a, tau_b, bridge, direct))

    # Sharp hostiles: a constant potential kills every edge, and merely
    # equalizing owners five and seven kills the forced bridge.
    constant_hostile = weighted_laplacian_determinant((1,) * 8)
    require(constant_hostile == 0, constant_hostile)
    bridge_hostile = list(private_counts(11))
    bridge_hostile[6] = bridge_hostile[4]
    require(weighted_laplacian_determinant(tuple(bridge_hostile)) == 0, bridge_hostile)

    direct_digest = sha256(
        dumps(direct_rows, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    semantic_payload = {
        "orientation": "owner index increasing",
        "packets_one_based": tuple(tuple(vertex + 1 for vertex in packet) for packet in PACKETS),
        "edges_one_based": tuple((left + 1, right + 1) for left, right in edges),
        "symbolic_rows": tuple(symbolic_rows),
        "cycle_dimension": len(cycles),
        "all_symbolic_cycle_pairings_zero": True,
        "direct_window": (1, 3000),
        "direct_digest": direct_digest,
        "first_three": tuple(direct_rows[:3]),
        "constant_hostile": constant_hostile,
        "bridge_hostile": 0,
        "dependencies": dependencies,
        "security": security_report(Path(__file__)),
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":"), default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, semantic_hash)

    print("THM-3482 PRIVATE-COUNT GRADIENT WEIGHTED SPECTRAL EXACT COMPANION")
    print("STATUS: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    for label, relative, digest in dependencies:
        print(f"DEPENDENCY: {label} {digest} {relative}")
    print("ORIENTATION_GAUGE: each edge is oriented from smaller to larger THM-3473 owner index")
    print(f"PACKETS_ONE_BASED: {semantic_payload['packets_one_based']}")
    print(f"SYMBOLIC_ROWS_STATE_TA_TB_BRIDGE_FULL_ZEROEDGES: {tuple(symbolic_rows)}")
    print("FACTORED_FULL_DETERMINANTS:")
    print("  k=0 mod3: -2048*k^4*(2k-1)^2*(2k+1)")
    print("  k=1 mod3: -2048*k^4*(8k^3+12k^2-2k-1)")
    print("  k=2 mod3: -256*k*(2k-1)^6")
    print(f"DIRECT_7x7_AUDIT_K1_TO3000_SHA256: {direct_digest}")
    print(f"FIRST_THREE_K_TA_TB_BRIDGE_DET: {tuple(direct_rows[:3])}")
    print("H1_PAIRING: all six symbolic cycle pairings vanish identically in every ternary state")
    print("HOSTILES: constant owner potential det=0; equalized (5,7) bridge det=0")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print("VERDICT: the canonical private-count coboundary has zero absolute graph H1 class but nonzero owner-order weighted 7x13 determinant for every k>=1; no LRC current or LRC(14) conclusion follows")


if __name__ == "__main__":
    main()
