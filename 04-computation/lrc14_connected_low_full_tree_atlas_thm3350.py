#!/usr/bin/env python3
"""Exact THM-3350 referee for the connected intrinsic-low six-level atlas.

This file is intentionally self-contained and writes nothing.  It performs
two complete enumerations with different duplicate control:

  A. ordinary layer BFS followed by set deduplication;
  B. canonical-parent reverse search (no cross-parent duplicates).

It then checks the graph statistics and the upper-median low-only zero
witness on the hostile reflected body H=(1,2,3,4,6,12).
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm

LOW_UP = frozenset((F(4, 3), F(3, 2), F(2), F(5, 2), F(3), F(4), F(5), F(6)))
LOW_ORIENTED = tuple(sorted(LOW_UP | {1 / x for x in LOW_UP}))
EXPECTED_LAYER_COUNTS = (1, 8, 94, 1295, 19389, 305909)
EXPECTED_EDGE_HIST = (
    (5, 171119), (6, 88196), (7, 32244), (8, 10046), (9, 3131),
    (10, 894), (11, 210), (12, 58), (13, 9), (14, 2),
)
EXPECTED_HIGH_COMPONENT_HIST = (
    ((2, 1, 1, 1, 1), 2), ((2, 2, 1, 1), 1),
    ((3, 1, 1, 1), 8), ((3, 3), 4), ((4, 1, 1), 210),
    ((4, 2), 8), ((5, 1), 3914), ((6,), 301762),
)
EXPECTED_ORBIT_DIGEST = "6d0b14b4608ba696bb107689afecf985a8ad70f5e5e3ada8698b70e5bede8681"
EXPECTED_LABELLED_LEDGER_DIGEST = "357926e866c6fd1b2705175a2caf57b30f44e04729e157f16ed66c54238b5a43"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def low_adjacent(x: F | int, y: F | int) -> bool:
    a, b = sorted((F(x), F(y)))
    require(a < b, (x, y))
    return b / a in LOW_UP


def connected(row: tuple[F, ...]) -> bool:
    unseen = set(row[1:])
    frontier = {row[0]}
    while frontier:
        x = frontier.pop()
        hit = {y for y in unseen if low_adjacent(x, y)}
        unseen -= hit
        frontier |= hit
    return not unseen


def primitive(row: tuple[F, ...]) -> tuple[int, ...]:
    denominator = lcm(*(x.denominator for x in row))
    values = tuple(int(x * denominator) for x in row)
    divisor = gcd(*values)
    return tuple(x // divisor for x in values)


def enumerate_set_bfs():
    """Completeness by arbitrary connected growth; global set deduplication."""
    layer = {(F(1),)}
    counts = [1]
    for _size in range(2, 7):
        children = set()
        for row in layer:
            have = set(row)
            for x in row:
                for ratio in LOW_ORIENTED:
                    y = x * ratio
                    if y >= 1 and y not in have:
                        children.add(tuple(sorted((*row, y))))
        layer = children
        counts.append(len(layer))
    require(tuple(counts) == EXPECTED_LAYER_COUNTS, ("BFS counts", counts))
    require(all(row[0] == 1 and connected(row) for row in layer), "BFS invariant")
    return layer, tuple(counts)


def canonical_parent(row: tuple[F, ...]) -> tuple[F, ...]:
    """Delete the greatest removable non-root vertex.

    A connected graph has at least two non-cut vertices, so one is not the
    distinguished root 1.  This makes the parent total and unique.
    """
    for value in reversed(row[1:]):
        candidate = tuple(x for x in row if x != value)
        if connected(candidate):
            return candidate
    raise RuntimeError(("no canonical parent", row))


def enumerate_reverse_search():
    """Duplicate-free search: accept a child only from its canonical parent."""
    layer = [(F(1),)]
    counts = [1]
    for _size in range(2, 7):
        children = []
        for row in layer:
            have = set(row)
            candidates = sorted(
                {x * r for x in row for r in LOW_ORIENTED if x * r >= 1 and x * r not in have}
            )
            for y in candidates:
                child = tuple(sorted((*row, y)))
                if canonical_parent(child) == row:
                    children.append(child)
        require(len(children) == len(set(children)), ("reverse duplicate", _size))
        layer = children
        counts.append(len(layer))
    require(tuple(counts) == EXPECTED_LAYER_COUNTS, ("reverse counts", counts))
    return layer, tuple(counts)


def orbit_digests(rows: list[tuple[int, ...]]):
    orbit = sha256()
    labelled = sha256()
    for row in rows:
        line = ",".join(map(str, row))
        orbit.update((line + "\n").encode())
        labelled.update((line + ";orbit=720\n").encode())
    return orbit.hexdigest(), labelled.hexdigest()


def edge_mask(row: tuple[int, ...]) -> int:
    mask = 0
    for bit, (i, j) in enumerate(combinations(range(6), 2)):
        if low_adjacent(row[i], row[j]):
            mask |= 1 << bit
    return mask


def component_sizes(mask: int, complement: bool = False) -> tuple[int, ...]:
    adjacency = [0] * 6
    for bit, (i, j) in enumerate(combinations(range(6), 2)):
        on = bool(mask & (1 << bit))
        if complement:
            on = not on
        if on:
            adjacency[i] |= 1 << j
            adjacency[j] |= 1 << i
    unseen = (1 << 6) - 1
    sizes = []
    while unseen:
        frontier = unseen & -unseen
        seen = 0
        while frontier:
            seen |= frontier
            unseen &= ~frontier
            neighbours = 0
            work = frontier
            while work:
                bit = work & -work
                work -= bit
                neighbours |= adjacency[bit.bit_length() - 1]
            frontier = neighbours & ~seen
        sizes.append(seen.bit_count())
    return tuple(sorted(sizes, reverse=True))


def clip(left: F, right: F):
    left, right = max(F(0), left), min(F(1), right)
    return (left, right) if left < right else None


def limiting_arcs(multiplier: int, label: int, cell: int, ruler: int):
    """The exact homogeneous clause ||P*x-(e*j mod L)/L||<1/14."""
    phase = (label * cell) % ruler
    radius = F(1, 14 * multiplier)
    arcs = []
    for k in range(-1, multiplier + 2):
        centre = F(phase, ruler * multiplier) + F(k, multiplier)
        piece = clip(centre - radius, centre + radius)
        if piece is not None:
            arcs.append(piece)
    return tuple(arcs)


def fracpart(value: F) -> F:
    return value - value.numerator // value.denominator


def bernoulli3(value: F) -> F:
    value = fracpart(value)
    return value**3 - F(3, 2) * value**2 + F(1, 2) * value


def homogeneous_limit(ruler: int, cell: int, e: int, p: int, f: int, q: int) -> F:
    """Exact common-dilation limit, including transverse body drift."""
    divisor = gcd(p, q)
    P, Q = p // divisor, q // divisor
    if P > Q:
        P, Q, e, f = Q, P, f, e
    R, S = (e * cell) % ruler, (f * cell) % ruler
    cross = Q * e - P * f
    if cross == 0:
        return intersection_mass(
            limiting_arcs(P, e, cell, ruler),
            limiting_arcs(Q, f, cell, ruler),
        )
    determinant = Q * R - P * S
    a, b = F(P, 14), F(Q, 14)
    u, v = F(determinant + cross, ruler), F(-determinant, ruler)
    psi = (
        bernoulli3(u + a - b) + bernoulli3(u - a + b)
        + bernoulli3(v + a - b) + bernoulli3(v - a + b)
        - bernoulli3(u + a + b) - bernoulli3(u - a - b)
        - bernoulli3(v + a + b) - bernoulli3(v - a - b)
    )
    return F(1, 49) + F(ruler, 6 * P * Q * cross) * psi


def reflected_level_arcs(ruler: int, label: int, level: int, cell: int):
    """Exact physical scale-one reflected clause, independent of repo code."""
    z = level * ruler - label
    phase = (label * cell) % ruler
    lo, hi = (-phase) // ruler - 1, (z - phase) // ruler + 1
    radius = F(ruler, 14 * z)
    arcs = []
    for k in range(lo, hi + 1):
        piece = clip(F(phase + k * ruler, z) - radius, F(phase + k * ruler, z) + radius)
        if piece is not None:
            arcs.append(piece)
    return tuple(sorted(set(arcs)))


def intersection_mass(first, second) -> F:
    i = j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        total += max(F(0), min(first[i][1], second[j][1]) - max(first[i][0], second[j][0]))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def safe_cell_ranges(body: tuple[int, ...]):
    ruler = 14 * lcm(*body)
    danger = []
    for label in body:
        radius, period = ruler // (14 * label), ruler // label
        for k in range(label + 1):
            centre = k * period
            danger.append((max(0, centre - radius), min(ruler, centre + radius)))
    merged = []
    for left, right in sorted(danger):
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    safe, cursor = [], 0
    for left, right in merged:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < ruler:
        safe.append((cursor, ruler))
    return ruler, tuple(safe)


def reduced_channel(x: int, y: int):
    divisor = gcd(x, y)
    return x // divisor, y // divisor


def kruskal(weights):
    parent = list(range(6))

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    total, tree = F(0), []
    for weight, edge, channel in sorted(weights, key=lambda x: (x[0], x[1]), reverse=True):
        a, b = map(root, edge)
        if a == b:
            continue
        parent[a] = b
        total += weight
        tree.append((weight, edge, channel))
        if len(tree) == 5:
            return total, tuple(tree)
    return None


def witness_audit():
    body = (1, 2, 3, 4, 6, 12)
    word = (1, 3, 9, 27, 108, 432)
    ruler, ranges = safe_cell_ranges(body)
    cells = tuple(j for left, right in ranges for j in range(left, right))
    require((ruler, len(cells), cells[len(cells) // 2]) == (168, 88, 90), (ruler, len(cells)))
    cell = 90
    low_edges = tuple(
        (i, j) for i, j in combinations(range(6), 2)
        if sum(reduced_channel(word[i], word[j])) <= 7
    )
    path = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5))
    require(low_edges == path, low_edges)

    limiting = []
    physical = []
    for edge in path:
        i, j = edge
        P, Q = reduced_channel(word[i], word[j])
        limiting.append(homogeneous_limit(ruler, cell, body[i], word[i], body[j], word[j]))
        physical.append(
            intersection_mass(
                reflected_level_arcs(ruler, body[i], word[i], cell),
                reflected_level_arcs(ruler, body[j], word[j], cell),
            )
        )
    require(tuple(limiting) == (F(0),) * 5, limiting)
    require(tuple(physical) == (F(0),) * 5, physical)
    debt = sum((F(e, 7 * (q * ruler - e)) for e, q in zip(body, word)), F(0))
    require(debt == F(5824887360027324, 3056780423523449161), debt)

    # A positive control: the same unique low tree at the first safe cell.
    alternate = []
    for i, j in path:
        P, Q = reduced_channel(word[i], word[j])
        alternate.append(homogeneous_limit(ruler, 12, body[i], word[i], body[j], word[j]))
    require(sum(alternate, F(0)) == F(29, 224), alternate)

    # The high complement is connected and repairs the low-only zero at j=90.
    high_limit, high_physical = [], []
    for edge in combinations(range(6), 2):
        if edge in path:
            continue
        i, j = edge
        P, Q = reduced_channel(word[i], word[j])
        channel = (P, Q)
        high_limit.append((
            homogeneous_limit(ruler, cell, body[i], word[i], body[j], word[j]),
            edge, channel,
        ))
        high_physical.append((
            intersection_mass(
                reflected_level_arcs(ruler, body[i], word[i], cell),
                reflected_level_arcs(ruler, body[j], word[j], cell),
            ), edge, channel,
        ))
    limit_tree = kruskal(high_limit)
    physical_tree = kruskal(high_physical)
    require(
        limit_tree is not None and limit_tree[0] == F(1453183, 12108096),
        limit_tree,
    )
    require(
        physical_tree is not None
        and physical_tree[0] == F(315022568681158, 2614867770336569),
        physical_tree,
    )
    margin = physical_tree[0] - debt
    require(margin == F(362436495428246378, 3056780423523449161), margin)
    return {
        "body": body, "ruler": ruler, "cell": cell, "word": word,
        "path": path, "limits": tuple(limiting), "physical": tuple(physical),
        "debt": debt, "alternate_cell": 12,
        "alternate_low_tree_limit": sum(alternate, F(0)),
        "high_limit_tree": limit_tree,
        "high_physical_tree": physical_tree,
        "high_physical_margin": margin,
    }


def main() -> None:
    bfs, bfs_counts = enumerate_set_bfs()
    bfs_rows = sorted(primitive(row) for row in bfs)
    require(len(bfs_rows) == len(set(bfs_rows)), "primitive normalization collision")

    reverse, reverse_counts = enumerate_reverse_search()
    reverse_rows = sorted(primitive(row) for row in reverse)
    require(bfs_rows == reverse_rows, "independent enumerators disagree")
    rows = bfs_rows

    orbit_digest, labelled_digest = orbit_digests(rows)
    require(orbit_digest == EXPECTED_ORBIT_DIGEST, orbit_digest)
    require(labelled_digest == EXPECTED_LABELLED_LEDGER_DIGEST, labelled_digest)
    require(len(rows) * 720 == 220254480, len(rows) * 720)
    max_height = max(row[-1] for row in rows)
    max_rows = tuple(row for row in rows if row[-1] == max_height)
    require((max_height, max_rows) == (7776, ((1, 6, 36, 216, 1296, 7776),)), max_rows)

    masks = tuple(edge_mask(row) for row in rows)
    edge_hist = tuple(sorted(Counter(mask.bit_count() for mask in masks).items()))
    require(edge_hist == EXPECTED_EDGE_HIST, edge_hist)
    fourteen = tuple(row for row, mask in zip(rows, masks) if mask.bit_count() == 14)
    require(
        fourteen == ((1, 2, 3, 4, 6, 12), (2, 3, 4, 6, 8, 12)),
        fourteen,
    )
    high_hist = tuple(sorted(Counter(component_sizes(mask, True) for mask in masks).items()))
    require(high_hist == EXPECTED_HIGH_COMPONENT_HIST, high_hist)
    require(sum(count for shape, count in high_hist if shape == (6,)) == 301762, high_hist)
    require(sum(count for shape, count in high_hist if shape != (6,)) == 4147, high_hist)
    coordinate_alphabet = len({x for row in rows for x in row})
    graph_masks = len(set(masks))
    require((coordinate_alphabet, graph_masks) == (126, 9233), (coordinate_alphabet, graph_masks))

    witness = witness_audit()
    print("LRC14 CONNECTED INTRINSIC-LOW SIX-LEVEL ATLAS")
    print("status=FINITE-EXACT;two independent complete enumerators;read-only exact referee")
    print(f"layer_counts={bfs_counts};reverse_layer_counts={reverse_counts}")
    print(f"unlabelled_orbits={len(rows)};labelled_vectors={len(rows)*720}")
    print(f"max_height={max_height};unique_max_orbit={max_rows[0]}")
    print(f"orbit_digest={orbit_digest}")
    print(f"labelled_orbit_ledger_digest={labelled_digest}")
    print(f"edge_count_hist={edge_hist}")
    print(f"fourteen_low_edge_orbits={fourteen}")
    print(f"coordinate_alphabet={coordinate_alphabet};sorted_position_graph_masks={graph_masks}")
    print(f"high_complement_component_hist={high_hist}")
    print("high_connected_orbits=301762;high_disconnected_orbits=4147")
    print(f"zero_witness={witness}")


if __name__ == "__main__":
    main()
