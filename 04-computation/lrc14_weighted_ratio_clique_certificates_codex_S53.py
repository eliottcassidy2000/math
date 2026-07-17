#!/usr/bin/env python3
"""Generate and replay exact weighted-ratio clique certificates for LRC(14)."""

from __future__ import annotations

import argparse
import json
import time
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path


TAU = {
    2: Fraction(6, 637),
    3: Fraction(2, 441),
    4: Fraction(5, 2646),
    5: Fraction(1, 1764),
    6: Fraction(1, 3136),
    7: Fraction(5, 37632),
    8: Fraction(1, 9996),
    9: Fraction(1, 14112),
}

TWO_PATH_CUTOFF = Fraction(4, 539)
ONE_PATH_CUTOFF = Fraction(5, 588)

CERTIFICATE_SPECS = tuple((level, TAU[level], level - 1) for level in range(3, 10))
SCHEMA = "lrc14-weighted-ratio-clique-v1"


def bernoulli2_residue(residue: int) -> Fraction:
    value = Fraction(residue % 14, 14)
    return value * value - value + Fraction(1, 6)


def pair_covariance(numerator: int, denominator: int) -> Fraction:
    return (
        bernoulli2_residue(numerator - denominator)
        - bernoulli2_residue(numerator + denominator)
    ) / (numerator * denominator)


def check_covariance_envelope() -> None:
    residue_differences = (
        abs(
            bernoulli2_residue(first - second)
            - bernoulli2_residue(first + second)
        )
        for first in range(14)
        for second in range(14)
    )
    if max(residue_differences) != Fraction(12, 49):
        raise AssertionError("pair-covariance envelope mismatch")


def strict_ratio_graph(threshold: Fraction) -> tuple[list[tuple[int, int]], list[int]]:
    envelope = Fraction(12, 49) / threshold
    product_cap = (envelope.numerator - 1) // envelope.denominator
    oriented_ratios: set[tuple[int, int]] = set()
    for numerator in range(1, isqrt(product_cap) + 1):
        for denominator in range(numerator + 1, product_cap // numerator + 1):
            if gcd(numerator, denominator) != 1:
                continue
            covariance = pair_covariance(numerator, denominator)
            if covariance < 0 and -covariance > threshold:
                oriented_ratios.add((numerator, denominator))
                oriented_ratios.add((denominator, numerator))
    vertices = sorted(oriented_ratios)
    vertex_set = set(vertices)
    adjacency = [0] * len(vertices)
    for right, (a, b) in enumerate(vertices):
        for left in range(right):
            c, d = vertices[left]
            quotient_numerator = a * d
            quotient_denominator = b * c
            common = gcd(quotient_numerator, quotient_denominator)
            quotient = (
                quotient_numerator // common,
                quotient_denominator // common,
            )
            if quotient in vertex_set:
                adjacency[right] |= 1 << left
                adjacency[left] |= 1 << right
    return vertices, adjacency


def greedy_color_order(candidate_mask: int, adjacency: list[int]) -> tuple[list[int], list[int]]:
    vertices: list[int] = []
    color_bounds: list[int] = []
    uncolored = candidate_mask
    color = 0
    while uncolored:
        color += 1
        available = uncolored
        while available:
            vertex_bit = available & -available
            vertex = vertex_bit.bit_length() - 1
            vertices.append(vertex)
            color_bounds.append(color)
            uncolored &= ~vertex_bit
            available &= ~vertex_bit
            available &= ~adjacency[vertex]
    return vertices, color_bounds


def generate_proof(adjacency: list[int], target_size: int) -> tuple[int, list[list[int]]]:
    memo: dict[tuple[int, int], int] = {}
    nodes: list[list[int] | None] = []

    def prove(candidate_mask: int, target: int) -> int:
        key = (candidate_mask, target)
        if key in memo:
            return memo[key]
        node_id = len(nodes)
        memo[key] = node_id
        nodes.append(None)
        if target <= 1:
            if candidate_mask:
                raise RuntimeError("false no-clique claim reached")
            nodes[node_id] = []
            return node_id
        order, bounds = greedy_color_order(candidate_mask, adjacency)
        working = candidate_mask
        branches: list[int] = []
        for position in range(len(order) - 1, -1, -1):
            if bounds[position] < target:
                break
            vertex = order[position]
            vertex_bit = 1 << vertex
            if working & vertex_bit:
                child = prove(working & adjacency[vertex], target - 1)
                branches.extend((vertex, child))
                working &= ~vertex_bit
        nodes[node_id] = branches
        return node_id

    root_mask = (1 << len(adjacency)) - 1
    root = prove(root_mask, target_size)
    return root, [node for node in nodes if node is not None]


def replay_proof(adjacency: list[int], target_size: int, root: int, nodes: list[list[int]]) -> None:
    contexts: dict[int, tuple[int, int]] = {}
    checked: set[int] = set()
    active: set[int] = set()

    def replay(node_id: int, candidate_mask: int, target: int) -> None:
        if not 0 <= node_id < len(nodes):
            raise AssertionError("node id out of range")
        context = (candidate_mask, target)
        previous = contexts.setdefault(node_id, context)
        if previous != context:
            raise AssertionError("proof node reused in a different state")
        if node_id in checked:
            return
        if node_id in active:
            raise AssertionError("cyclic proof")
        active.add(node_id)
        branches = nodes[node_id]
        if len(branches) % 2:
            raise AssertionError("malformed branch list")
        if target <= 1:
            if candidate_mask or branches:
                raise AssertionError("invalid size-one leaf")
        else:
            order, bounds = greedy_color_order(candidate_mask, adjacency)
            working = candidate_mask
            branch_position = 0
            for position in range(len(order) - 1, -1, -1):
                if bounds[position] < target:
                    break
                vertex = order[position]
                vertex_bit = 1 << vertex
                if working & vertex_bit:
                    if branch_position + 1 >= len(branches):
                        raise AssertionError("missing branch")
                    stored_vertex = branches[branch_position]
                    child = branches[branch_position + 1]
                    if stored_vertex != vertex:
                        raise AssertionError("branch vertex mismatch")
                    replay(child, working & adjacency[vertex], target - 1)
                    branch_position += 2
                    working &= ~vertex_bit
            if branch_position != len(branches):
                raise AssertionError("extraneous branch")
        active.remove(node_id)
        checked.add(node_id)

    replay(root, (1 << len(adjacency)) - 1, target_size)
    if len(checked) != len(nodes):
        raise AssertionError("unreachable proof nodes")


def certificate_path(directory: Path, level: int) -> Path:
    return directory / f"strict_above_tau{level}.json"


def generate_certificates(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for level, threshold, target in CERTIFICATE_SPECS:
        vertices, adjacency = strict_ratio_graph(threshold)
        root, nodes = generate_proof(adjacency, target)
        artifact = {
            "schema": SCHEMA,
            "level": level,
            "threshold": [threshold.numerator, threshold.denominator],
            "target_ratio_clique": target,
            "vertex_count": len(vertices),
            "edge_count": sum(mask.bit_count() for mask in adjacency) // 2,
            "root": root,
            "nodes": nodes,
        }
        certificate_path(directory, level).write_text(
            json.dumps(artifact, separators=(",", ":"), sort_keys=True) + "\n",
            encoding="utf-8",
        )


def check_certificate(path: Path) -> dict[str, int | str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("schema") != SCHEMA:
        raise AssertionError("certificate schema mismatch")
    level = artifact["level"]
    threshold = Fraction(*artifact["threshold"])
    target = artifact["target_ratio_clique"]
    expected = (level, TAU[level], level - 1)
    if (level, threshold, target) != expected:
        raise AssertionError("certificate specification mismatch")
    vertices, adjacency = strict_ratio_graph(threshold)
    edge_count = sum(mask.bit_count() for mask in adjacency) // 2
    if len(vertices) != artifact["vertex_count"] or edge_count != artifact["edge_count"]:
        raise AssertionError("graph fingerprint mismatch")
    replay_proof(adjacency, target, artifact["root"], artifact["nodes"])
    return {
        "file": path.name,
        "threshold": str(threshold),
        "vertices": len(vertices),
        "edges": edge_count,
        "nodes": len(artifact["nodes"]),
        "branches": sum(len(node) // 2 for node in artifact["nodes"]),
        "bytes": path.stat().st_size,
    }


def turan_edges(vertex_count: int, clique_cap: int) -> int:
    quotient, remainder = divmod(vertex_count, clique_cap)
    parts = [quotient + 1] * remainder + [quotient] * (clique_cap - remainder)
    return (vertex_count * vertex_count - sum(part * part for part in parts)) // 2


def pair_bound() -> tuple[Fraction, Fraction]:
    bound = 78 * TAU[9]
    for clique_cap in range(8, 2, -1):
        bound += turan_edges(13, clique_cap) * (TAU[clique_cap] - TAU[clique_cap + 1])
    bound += turan_edges(13, 2) * (TWO_PATH_CUTOFF - TAU[3])
    bound += 24 * (ONE_PATH_CUTOFF - TWO_PATH_CUTOFF)
    bound += 12 * (TAU[2] - ONE_PATH_CUTOFF)
    return bound, Fraction(13, 30) - bound


def check_directory(directory: Path) -> list[dict[str, int | str]]:
    check_covariance_envelope()
    summaries = [check_certificate(certificate_path(directory, level)) for level in range(3, 10)]
    above_maximum_vertices, _ = strict_ratio_graph(TAU[2])
    one_path_vertices, _ = strict_ratio_graph(ONE_PATH_CUTOFF)
    two_path_vertices, _ = strict_ratio_graph(TWO_PATH_CUTOFF)
    if above_maximum_vertices:
        raise AssertionError("negative pair above tau_2")
    if one_path_vertices != [(1, 13), (13, 1)]:
        raise AssertionError("one-path ratio-color classification mismatch")
    if two_path_vertices != [(1, 12), (1, 13), (12, 1), (13, 1)]:
        raise AssertionError("two-path ratio-color classification mismatch")
    bound, margin = pair_bound()
    if bound != Fraction(176738453, 411675264):
        raise AssertionError("pair-bound arithmetic mismatch")
    if margin != Fraction(8270807, 2058376320) or margin <= 0:
        raise AssertionError("pair-bound margin mismatch")
    return summaries


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("generate", "check"))
    parser.add_argument("certificate_directory", type=Path)
    arguments = parser.parse_args()
    started = time.perf_counter()
    if arguments.mode == "generate":
        generate_certificates(arguments.certificate_directory)
    summaries = check_directory(arguments.certificate_directory)
    for summary in summaries:
        print(json.dumps(summary, sort_keys=True))
    bound, margin = pair_bound()
    print(f"negative_pair_mass_bound={bound} ({float(bound):.15f})")
    print(f"margin_to_13_over_30={margin} ({float(margin):.15f})")
    print(f"runtime_seconds={time.perf_counter() - started:.6f}")


if __name__ == "__main__":
    main()
