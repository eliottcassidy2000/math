#!/usr/bin/env python3
"""Independent switching-class and Held--Karp audit for THM-4084.

Unlike the primary audit, this path never materializes a list of cycles.  It
counts negative Hamilton cycles on every vertex subset by a parity-refined
Held--Karp recurrence, exhausts every labelled switching class for n=6,7,
computes frustration by all cuts, and classifies the complete index-at-most-two
universe.  It separately checks every matching layer through n=10.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from math import comb, factorial


EXPECTED_SEMANTIC_SHA256 = "9c6bf209688101bcb18b79f1f23986939a8dc807c84145eb318b696d520673ff"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def edge(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def edge_index(n: int) -> dict[tuple[int, int], int]:
    return {item: index for index, item in enumerate(itertools.combinations(range(n), 2))}


def mask_from_edges(items: tuple[tuple[int, int], ...], indices: dict[tuple[int, int], int]) -> int:
    mask = 0
    for a, b in items:
        mask |= 1 << indices[edge(a, b)]
    return mask


def cut_masks(n: int, indices: dict[tuple[int, int], int]) -> list[int]:
    masks: list[int] = []
    for subset_bits in range(1 << (n - 1)):
        chosen = {vertex for vertex in range(1, n) if subset_bits & (1 << (vertex - 1))}
        mask = 0
        for a, b in indices:
            if (a in chosen) != (b in chosen):
                mask |= 1 << indices[(a, b)]
        masks.append(mask)
    require(len(set(masks)) == 1 << (n - 1), f"cut masks collided n={n}")
    return masks


def edge_parity(signing: int, a: int, b: int, indices: dict[tuple[int, int], int]) -> int:
    return (signing >> indices[edge(a, b)]) & 1


def negative_cycles_dp(signing: int, vertices: tuple[int, ...], indices: dict[tuple[int, int], int]) -> int:
    root = vertices[0]
    tail = vertices[1:]
    width = len(tail)
    states: dict[tuple[int, int, int], int] = {}
    for index, vertex in enumerate(tail):
        parity = edge_parity(signing, root, vertex, indices)
        states[(1 << index, index, parity)] = 1

    for used_size in range(1, width):
        next_states = dict(states)
        for (used, last, parity), count in states.items():
            if used.bit_count() != used_size:
                continue
            for nxt in range(width):
                if used & (1 << nxt):
                    continue
                new_used = used | (1 << nxt)
                new_parity = parity ^ edge_parity(signing, tail[last], tail[nxt], indices)
                key = (new_used, nxt, new_parity)
                next_states[key] = next_states.get(key, 0) + count
        states = next_states

    full = (1 << width) - 1
    negative_oriented = 0
    for last in range(width):
        for parity in (0, 1):
            if parity ^ edge_parity(signing, tail[last], root, indices):
                negative_oriented += states.get((full, last, parity), 0)
    require(negative_oriented % 2 == 0, f"cycle reversal parity failed vertices={vertices}")
    return negative_oriented // 2


def negative_profile_dp(signing: int, n: int, indices: dict[tuple[int, int], int], lengths: tuple[int, ...]) -> dict[int, int]:
    profile: dict[int, int] = {}
    for length in lengths:
        profile[length] = sum(
            negative_cycles_dp(signing, vertices, indices)
            for vertices in itertools.combinations(range(n), length)
        )
    return profile


def profile_formula(n: int, length: int, size: int) -> int:
    return sum(
        (-4) ** (chosen - 1)
        * comb(size, chosen)
        * factorial(length - chosen - 1)
        * comb(n - 2 * chosen, length - 2 * chosen)
        for chosen in range(1, min(size, length // 2) + 1)
    )


def edge_profile(n: int, length: int) -> int:
    return factorial(n - 2) // factorial(n - length)


def canonical_matching(size: int) -> tuple[tuple[int, int], ...]:
    return tuple((2 * index, 2 * index + 1) for index in range(size))


def edge_shape(items: tuple[tuple[int, int], ...]) -> str:
    degrees: dict[int, int] = {}
    for a, b in items:
        degrees[a] = degrees.get(a, 0) + 1
        degrees[b] = degrees.get(b, 0) + 1
    signature = tuple(sorted(degrees.values()))
    names = {
        (1, 1): "single",
        (1, 1, 1, 1): "disjoint",
        (1, 1, 2): "adjacent",
        (1, 1, 1, 1, 1, 1): "matching3",
        (1, 1, 1, 1, 2): "path2_plus_edge",
        (1, 1, 2, 2): "path3",
        (1, 1, 1, 3): "star3",
        (2, 2, 2): "triangle",
    }
    require(signature in names, f"unknown edge shape {signature}")
    return names[signature]


def cycles_containing_edges(n: int, length: int, items: tuple[tuple[int, int], ...]) -> int:
    vertices = {vertex for item in items for vertex in item}
    degrees = {vertex: 0 for vertex in vertices}
    neighbors = {vertex: set() for vertex in vertices}
    for a, b in items:
        degrees[a] += 1
        degrees[b] += 1
        neighbors[a].add(b)
        neighbors[b].add(a)
    if max(degrees.values()) > 2:
        return 0

    components = 0
    unseen = set(vertices)
    while unseen:
        components += 1
        stack = [unseen.pop()]
        while stack:
            vertex = stack.pop()
            for nxt in neighbors[vertex]:
                if nxt in unseen:
                    unseen.remove(nxt)
                    stack.append(nxt)

    edge_count = len(items)
    if all(degree == 2 for degree in degrees.values()):
        return 1 if components == 1 and length == edge_count else 0
    vertex_count = len(vertices)
    if length < vertex_count:
        return 0
    return (
        2 ** (components - 1)
        * factorial(length - edge_count - 1)
        * comb(n - vertex_count, length - vertex_count)
    )


def inclusion_exclusion_profile(n: int, length: int, items: tuple[tuple[int, int], ...]) -> int:
    total = 0
    for size in range(1, len(items) + 1):
        for chosen in itertools.combinations(items, size):
            total += (-2) ** (size - 1) * cycles_containing_edges(n, length, chosen)
    return total


def exhaustive_low_frustration(n: int) -> tuple[dict[int, int], int, int, int]:
    indices = edge_index(n)
    gauge_edges = tuple(item for item in indices if 0 not in item)
    cuts = cut_masks(n, indices)
    class_counts: dict[int, int] = {}
    low_classes = 0
    dp_cycle_gates = 0
    d5_equality_classes = 0
    edge_value = sum(edge_profile(n, length) for length in range(3, 7))

    for gauge_bits in range(1 << len(gauge_edges)):
        signing = mask_from_edges(
            tuple(gauge_edges[index] for index in range(len(gauge_edges)) if gauge_bits & (1 << index)),
            indices,
        )
        switched = [(signing ^ cut).bit_count() for cut in cuts]
        frustration = min(switched)
        class_counts[frustration] = class_counts.get(frustration, 0) + 1
        if frustration > 3 or frustration == 0:
            continue

        low_classes += 1
        minimum_representatives = {signing ^ cut for cut, weight in zip(cuts, switched) if weight == frustration}
        shapes: set[str] = set()
        for representative in minimum_representatives:
            negative_edges = [item for item, position in indices.items() if representative & (1 << position)]
            shapes.add(edge_shape(tuple(negative_edges)))
        require(len(shapes) == 1, f"minimum shape was not intrinsic n={n},mask={signing},shapes={shapes}")
        shape = next(iter(shapes))

        profile = negative_profile_dp(signing, n, indices, (3, 4, 5, 6))
        dp_cycle_gates += sum(comb(n, length) for length in range(3, 7))
        if shape == "single":
            expected = {length: edge_profile(n, length) for length in range(3, 7)}
        elif shape == "disjoint":
            expected = {length: profile_formula(n, length, 2) for length in range(3, 7)}
        elif shape == "adjacent":
            expected = {
                length: 2 * edge_profile(n, length) - 2 * factorial(n - 3) // factorial(n - length)
                for length in range(3, 7)
            }
        else:
            representative_edges = tuple(
                item for item, position in indices.items() if next(iter(minimum_representatives)) & (1 << position)
            )
            expected = {
                length: inclusion_exclusion_profile(n, length, representative_edges)
                for length in range(3, 7)
            }
        require(profile == expected, f"Held--Karp low-frustration profile failed n={n},mask={signing},shape={shape}")
        value = sum(profile.values())
        require(value >= edge_value, f"low-frustration D5 value failed n={n},mask={signing}")
        if value == edge_value:
            require(shape == "single", f"extra low-frustration equality n={n},mask={signing}")
            d5_equality_classes += 1

    require(d5_equality_classes == comb(n, 2), f"single-edge equality multiplicity changed n={n}")
    return class_counts, low_classes, dp_cycle_gates, d5_equality_classes


def all_layer_matching_dp() -> tuple[int, list[dict[str, int]]]:
    dp_gates = 0
    rows: list[dict[str, int]] = []
    for n in range(6, 11):
        indices = edge_index(n)
        cuts = cut_masks(n, indices)
        row_profiles = 0
        for size in range(1, n // 2 + 1):
            signing = mask_from_edges(canonical_matching(size), indices)
            require(min((signing ^ cut).bit_count() for cut in cuts) == size, f"DP matching frustration changed n={n},r={size}")
            profile = negative_profile_dp(signing, n, indices, tuple(range(3, n + 1)))
            for length, value in profile.items():
                require(value == profile_formula(n, length, size), f"DP all-layer formula failed n={n},r={size},k={length}")
                dp_gates += comb(n, length)
            row_profiles += 1
        rows.append({"n": n, "matching_sizes": row_profiles})
    return dp_gates, rows


def boundary_control() -> tuple[int, int]:
    n = 4
    indices = edge_index(n)
    signing = mask_from_edges(((0, 1), (2, 3)), indices)
    profile = negative_profile_dp(signing, n, indices, (3, 4))
    value = sum(profile.values())
    edge_value = sum(edge_profile(n, length) for length in (3, 4))
    require(profile == {3: 4, 4: 0}, f"K4 matching boundary profile changed: {profile}")
    require(value == edge_value == 4, "K4 D3 boundary tie changed")
    return value, edge_value


def main() -> None:
    exhaustive_rows: list[dict[str, object]] = []
    total_low_classes = 0
    total_exhaustive_dp_gates = 0
    for n in (6, 7):
        class_counts, low_classes, dp_gates, equalities = exhaustive_low_frustration(n)
        exhaustive_rows.append(
            {
                "n": n,
                "class_counts": {str(key): value for key, value in sorted(class_counts.items())},
                "low_nontrivial_classes": low_classes,
                "d5_equalities": equalities,
            }
        )
        total_low_classes += low_classes
        total_exhaustive_dp_gates += dp_gates

    all_layer_dp_gates, matching_rows = all_layer_matching_dp()
    boundary_value, boundary_edge = boundary_control()
    payload = {
        "all_layer_dp_gates": all_layer_dp_gates,
        "boundary": [boundary_value, boundary_edge],
        "exhaustive_dp_gates": total_exhaustive_dp_gates,
        "exhaustive_rows": exhaustive_rows,
        "matching_rows": matching_rows,
        "total_low_classes": total_low_classes,
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()
    require(digest == EXPECTED_SEMANTIC_SHA256, f"semantic digest changed: {digest}")
    print("THM-4084 independent switching-class Held-Karp audit")
    print("exhaustive_class_range=6,7 all_layer_matching_range=6..10")
    print(f"low_frustration_classes={total_low_classes} exhaustive_subset_cycle_gates={total_exhaustive_dp_gates}")
    print(f"all_layer_subset_cycle_gates={all_layer_dp_gates}")
    print(f"exhaustive_rows={json.dumps(exhaustive_rows, sort_keys=True, separators=(',', ':'))}")
    print(f"matching_rows={json.dumps(matching_rows, sort_keys=True, separators=(',', ':'))}")
    print(f"n4_D3_boundary=matching:{boundary_value},edge:{boundary_edge}")
    print(f"semantic_sha256={digest}")
    print("PASS: exhaustive frustration<=3 firewall, independent all-layer matching profile, and sharp boundary")


if __name__ == "__main__":
    main()
