#!/usr/bin/env python3
"""Exact audit for THM-4113.

The first path transcribes Appendix B.1.2 of arXiv:2608.22870v1.  The
second path generates every noncrossing partial matching and tests maximality
directly.  The third path uses the proved 2--3-block coefficient formula.
All arithmetic is integral and the checks remain active under ``python -O``.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from functools import lru_cache
from itertools import combinations
from math import factorial


Pair = tuple[int, int]
Matching = tuple[Pair, ...]
Block = tuple[int, ...]
Partition = tuple[Block, ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def canon_matching(pairs: tuple[Pair, ...] | list[Pair]) -> Matching:
    return tuple(sorted((min(a, b), max(a, b)) for a, b in pairs))


def canon_partition(blocks: tuple[Block, ...] | list[Block]) -> Partition:
    return tuple(sorted(tuple(sorted(block)) for block in blocks))


def crossing(first: Pair, second: Pair) -> bool:
    a, b = first
    c, d = second
    return (a < c < b < d) or (c < a < d < b)


def noncrossing(matching: Matching) -> bool:
    return all(not crossing(p, q) for p, q in combinations(matching, 2))


@lru_cache(maxsize=None)
def catalan_matchings(q: int) -> tuple[Matching, ...]:
    """The GetValidParens(q) family on [0,2q)."""
    if q == 0:
        return ((),)
    result: set[Matching] = set()
    for partner in range(1, 2 * q, 2):
        left_q = (partner - 1) // 2
        right_q = q - 1 - left_q
        for left in catalan_matchings(left_q):
            for right in catalan_matchings(right_q):
                pairs: list[Pair] = [(0, partner)]
                pairs.extend((a + 1, b + 1) for a, b in left)
                pairs.extend((a + partner + 1, b + partner + 1) for a, b in right)
                result.add(canon_matching(pairs))
    return tuple(sorted(result))


@lru_cache(maxsize=None)
def paper_half_kempe_atlas(n: int) -> tuple[Matching, ...]:
    """Literal set-valued transcription of GetPlanarHalfKempes(n)."""
    require(n >= 0, "negative boundary size")
    if n <= 1:
        return ((),)
    all_generated: set[Matching] = set()
    result: set[Matching] = set()
    universe = tuple(range(n))
    for q in range(n // 2, 0, -1):
        for matched_vertices in combinations(universe, 2 * q):
            matched_set = set(matched_vertices)
            unmatched = tuple(v for v in universe if v not in matched_set)
            unmatched_set = set(unmatched)
            if any((v + 1) % n in unmatched_set for v in unmatched):
                continue
            for base in catalan_matchings(q):
                matching = canon_matching(
                    [(matched_vertices[a], matched_vertices[b]) for a, b in base]
                )
                all_generated.add(matching)
                redundant = any(
                    canon_matching(list(matching) + [(x, y)]) in all_generated
                    for x, y in combinations(unmatched, 2)
                )
                if not redundant:
                    result.add(matching)
    return tuple(sorted(result))


@lru_cache(maxsize=None)
def all_noncrossing_partial_matchings(n: int) -> tuple[Matching, ...]:
    """Independent Motzkin recursion, with no Appendix-B filtering."""
    require(n >= 0, "negative matching size")
    if n == 0:
        return ((),)
    result: set[Matching] = set()
    for rest in all_noncrossing_partial_matchings(n - 1):
        result.add(canon_matching([(a + 1, b + 1) for a, b in rest]))
    for partner in range(1, n):
        for inside in all_noncrossing_partial_matchings(partner - 1):
            for outside in all_noncrossing_partial_matchings(n - partner - 1):
                pairs: list[Pair] = [(0, partner)]
                pairs.extend((a + 1, b + 1) for a, b in inside)
                pairs.extend((a + partner + 1, b + partner + 1) for a, b in outside)
                result.add(canon_matching(pairs))
    return tuple(sorted(result))


def unmatched_vertices(n: int, matching: Matching) -> tuple[int, ...]:
    used = {v for pair in matching for v in pair}
    return tuple(v for v in range(n) if v not in used)


def is_maximal_noncrossing(n: int, matching: Matching) -> bool:
    require(noncrossing(matching), "maximality test received a crossing matching")
    unmatched = unmatched_vertices(n, matching)
    return not any(
        all(not crossing((x, y), pair) for pair in matching)
        for x, y in combinations(unmatched, 2)
    )


def deepest_parent(vertex: int, matching: Matching) -> Pair | None:
    containers = [(b - a, a, b) for a, b in matching if a < vertex < b]
    if not containers:
        return None
    _, a, b = min(containers)
    return (a, b)


def encode_matching(n: int, matching: Matching) -> tuple[str, object]:
    """Map a maximal matching to C or x*C^2 2--3-block data."""
    require(is_maximal_noncrossing(n, matching), "encoding a nonmaximal matching")
    by_parent: dict[Pair | None, list[int]] = defaultdict(list)
    for vertex in unmatched_vertices(n, matching):
        by_parent[deepest_parent(vertex, matching)].append(vertex)
    require(all(len(vertices) <= 1 for vertices in by_parent.values()), "two holes in one region")

    blocks: list[Block] = []
    for pair in matching:
        holes = by_parent.get(pair, [])
        blocks.append(tuple(sorted(pair + tuple(holes))))
    root_holes = by_parent.get(None, [])
    if not root_holes:
        return ("C", canon_partition(blocks))

    require(len(root_holes) == 1, "more than one root-region singleton")
    root = root_holes[0]
    require(all(not (a < root < b) for a, b in matching), "root singleton is enclosed")
    left = canon_partition([block for block in blocks if max(block) < root])
    right = canon_partition([block for block in blocks if min(block) > root])
    require(sum(map(len, left)) + sum(map(len, right)) + 1 == n, "root split lost vertices")
    return ("X", (root, left, right))


def matching_from_blocks(blocks: Partition) -> Matching:
    pairs: list[Pair] = []
    for block in blocks:
        require(len(block) in (2, 3), "invalid block size")
        pairs.append((block[0], block[-1]))
    return canon_matching(pairs)


def decode_encoding(encoded: tuple[str, object]) -> Matching:
    kind, payload = encoded
    if kind == "C":
        return matching_from_blocks(payload)  # type: ignore[arg-type]
    require(kind == "X", "unknown encoding kind")
    _, left, right = payload  # type: ignore[misc]
    return canon_matching(list(matching_from_blocks(left)) + list(matching_from_blocks(right)))


def shift_partition(partition: Partition, offset: int) -> Partition:
    return canon_partition([tuple(v + offset for v in block) for block in partition])


@lru_cache(maxsize=None)
def partitions_23(n: int) -> tuple[Partition, ...]:
    """All noncrossing partitions of [0,n) into blocks of size 2 or 3."""
    require(n >= 0, "negative partition size")
    if n == 0:
        return ((),)
    result: set[Partition] = set()
    for second in range(1, n):
        for inside in partitions_23(second - 1):
            for outside in partitions_23(n - second - 1):
                blocks = [(0, second)]
                blocks.extend(shift_partition(inside, 1))
                blocks.extend(shift_partition(outside, second + 1))
                result.add(canon_partition(blocks))
    for second in range(1, n):
        for third in range(second + 1, n):
            for first_gap in partitions_23(second - 1):
                for second_gap in partitions_23(third - second - 1):
                    for outside in partitions_23(n - third - 1):
                        blocks = [(0, second, third)]
                        blocks.extend(shift_partition(first_gap, 1))
                        blocks.extend(shift_partition(second_gap, second + 1))
                        blocks.extend(shift_partition(outside, third + 1))
                        result.add(canon_partition(blocks))
    return tuple(sorted(result))


def formula_count(n: int, singletons: int) -> int:
    """Coefficient of u^singletons x^n in C+u*x*C^2."""
    require(n >= 0 and singletons >= 0, "negative coefficient index")
    total = 0
    first_numerator = n - 3 * singletons
    if first_numerator >= 0 and first_numerator % 2 == 0:
        pairs = first_numerator // 2
        rest = n + 1 - pairs - singletons
        total += factorial(n) // (
            factorial(pairs) * factorial(singletons) * factorial(rest)
        )
    second_numerator = n + 2 - 3 * singletons
    if singletons >= 1 and second_numerator >= 0 and second_numerator % 2 == 0:
        pairs = second_numerator // 2
        rest = n + 2 - pairs - singletons
        total += 2 * factorial(n) // (
            factorial(pairs) * factorial(singletons - 1) * factorial(rest)
        )
    return total


def connected_component_sizes(vertices: tuple[Pair, ...], edges: set[tuple[Pair, Pair]]) -> list[int]:
    adjacency = {vertex: set() for vertex in vertices}
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    unseen = set(vertices)
    sizes: list[int] = []
    while unseen:
        root = min(unseen)
        queue = deque([root])
        unseen.remove(root)
        size = 0
        while queue:
            vertex = queue.popleft()
            size += 1
            for neighbour in adjacency[vertex]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    queue.append(neighbour)
        sizes.append(size)
    return sorted(sizes)


def main() -> None:
    max_n = 14
    direct_max_n = 12
    paper_formula_gates = 0
    direct_gates = 0
    encode_gates = 0
    type_a_gates = 0
    type_x_gates = 0

    print("THM-4113 maximal noncrossing half-Kempe atlas exact audit")
    print("n  total  singleton_profile")
    totals: list[int] = []
    for n in range(max_n + 1):
        paper = set(paper_half_kempe_atlas(n))
        profile = Counter(len(unmatched_vertices(n, matching)) for matching in paper)
        predicted = {
            singletons: formula_count(n, singletons)
            for singletons in range(n + 1)
            if formula_count(n, singletons)
        }
        require(dict(sorted(profile.items())) == predicted, f"formula mismatch at n={n}")
        paper_formula_gates += len(profile) + 1
        totals.append(len(paper))

        for matching in paper:
            require(is_maximal_noncrossing(n, matching), f"paper output nonmaximal at n={n}")
            encoded = encode_matching(n, matching)
            require(decode_encoding(encoded) == matching, f"encoding failed at n={n}")
            encode_gates += 2

        if n <= direct_max_n:
            direct = {
                matching
                for matching in all_noncrossing_partial_matchings(n)
                if is_maximal_noncrossing(n, matching)
            }
            require(paper == direct, f"direct maximal atlas mismatch at n={n}")
            direct_gates += len(all_noncrossing_partial_matchings(n)) + 1

        type_a = set(partitions_23(n))
        require(len(type_a) == sum(1 for m in paper if encode_matching(n, m)[0] == "C"),
                f"type-C count mismatch at n={n}")
        for partition in type_a:
            matching = matching_from_blocks(partition)
            require(is_maximal_noncrossing(n, matching), f"type-C inverse nonmaximal at n={n}")
            require(encode_matching(n, matching) == ("C", partition), f"type-C inverse failed at n={n}")
            type_a_gates += 2

        type_x_count = 0
        for root in range(n):
            for left_relative in partitions_23(root):
                left = shift_partition(left_relative, 0)
                for right_relative in partitions_23(n - root - 1):
                    right = shift_partition(right_relative, root + 1)
                    encoded = ("X", (root, left, right))
                    matching = decode_encoding(encoded)
                    require(is_maximal_noncrossing(n, matching), f"type-X inverse nonmaximal at n={n}")
                    require(encode_matching(n, matching) == encoded, f"type-X inverse failed at n={n}")
                    type_x_count += 1
                    type_x_gates += 2
        require(type_x_count == sum(1 for m in paper if encode_matching(n, m)[0] == "X"),
                f"type-X count mismatch at n={n}")

        rendered = ",".join(f"{s}:{count}" for s, count in sorted(profile.items()))
        print(f"{n:2d} {len(paper):6d}  {rendered}")

    require(totals[:11] == [1, 1, 1, 3, 4, 10, 20, 42, 98, 210, 492],
            "initial sequence changed")

    pair_vertices = tuple(combinations(range(5), 2))
    petersen_edges: set[tuple[Pair, Pair]] = set()
    crossing_edges: set[tuple[Pair, Pair]] = set()
    planar_edges: set[tuple[Pair, Pair]] = set()
    for first, second in combinations(pair_vertices, 2):
        if set(first).isdisjoint(second):
            edge = tuple(sorted((first, second)))
            petersen_edges.add(edge)
            if crossing(first, second):
                crossing_edges.add(edge)
            else:
                planar_edges.add(edge)
    require(len(petersen_edges) == 15, "K(5,2) edge count")
    require(len(crossing_edges) == 5 and len(planar_edges) == 10, "pentagon split")
    atlas5_edges = {
        tuple(sorted(matching))
        for matching in paper_half_kempe_atlas(5)
    }
    require(atlas5_edges == planar_edges, "five-boundary atlas is not the planar Petersen edge set")
    crossing_degrees = Counter(vertex for edge in crossing_edges for vertex in edge)
    require(sorted(crossing_degrees.values()) == [2] * 5, "crossing subgraph is not C5")
    require(connected_component_sizes(tuple(crossing_degrees), crossing_edges) == [5],
            "crossing active subgraph is disconnected")
    planar_degrees = Counter(vertex for edge in planar_edges for vertex in edge)
    require(sorted(planar_degrees.values()) == [1] * 5 + [3] * 5,
            "planar Petersen degree profile")
    require(connected_component_sizes(pair_vertices, planar_edges) == [10],
            "planar Petersen survivor is disconnected")

    print(f"paper_formula_gates={paper_formula_gates}")
    print(f"direct_motzkin_gates={direct_gates} through_n={direct_max_n}")
    print(f"matching_encode_decode_gates={encode_gates}")
    print(f"type_C_bijection_gates={type_a_gates}")
    print(f"type_X_bijection_gates={type_x_gates}")
    print("petersen_vertices=10 disjoint_edges=15 crossing_C5_edges=5 planar_atlas_edges=10")
    print("petersen_planar_degree_profile=1^5,3^5 components=10")
    print("result=PASS")


if __name__ == "__main__":
    main()
