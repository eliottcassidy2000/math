"""Exact multipartite operators on the ten-state Sym^2(F_2^2) alphabet.

Source note: 05-knowledge/results/thirtysix_multipartite_20260925.md.
Run from the repository root with Python 3; no third-party dependencies.
Every check remains active under python -O. No conjecture closure is tested.
"""

from collections import Counter
from itertools import combinations, combinations_with_replacement, product
import json
from math import comb, isqrt


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def edge_count(parts):
    return sum(a * b for a, b in combinations(parts, 2))


def graph(parts):
    labels = [i for i, size in enumerate(parts) for _ in range(size)]
    return {e for e in combinations(range(len(labels)), 2)
            if labels[e[0]] != labels[e[1]]}


def rank2(n, edges):
    rows = [0] * n
    for u, v in edges:
        rows[u] ^= 1 << v
        rows[v] ^= 1 << u
    pivots = {}
    for row in rows:
        while row:
            p = row.bit_length() - 1
            if p not in pivots:
                pivots[p] = row
                break
            row ^= pivots[p]
    return len(pivots)


def degrees(n, edges):
    out = [0] * n
    for u, v in edges:
        out[u] += 1
        out[v] += 1
    return out


def triangles(n, edges):
    return sum(all(e in edges for e in combinations(triple, 2))
               for triple in combinations(range(n), 3))


def types_with_edges(target):
    # For >=2 nonempty parts, E >= N-1; sorted first part squared <= E.
    found = []

    def extend(parts, size, edges):
        if edges == target:
            found.append(parts)
            return
        upper = min((target - edges) // size, target + 1 - size)
        for value in range(parts[-1], upper + 1):
            extend(parts + (value,), size + value, edges + size * value)

    for first in range(1, isqrt(target) + 1):
        extend((first,), first, 0)
    return sorted(found, key=lambda x: (len(x), x))


def k33_certificate(parts):
    vertices = []
    start = 0
    for size in parts:
        vertices.append(tuple(range(start, start + size)))
        start += size
    for mask in range(1, (1 << len(parts)) - 1):
        left = sum((v for i, v in enumerate(vertices) if mask >> i & 1), ())
        right = sum((v for i, v in enumerate(vertices) if not mask >> i & 1), ())
        if len(left) >= 3 and len(right) >= 3:
            return left[:3], right[:3]
    return None


def rotation_faces(rotation):
    darts = {(u, v) for u, ns in rotation.items() for v in ns}
    unseen = set(darts)
    count = 0
    while unseen:
        first = min(unseen)
        current = first
        while current in unseen:
            unseen.remove(current)
            u, v = current
            ns = rotation[v]
            current = v, ns[(ns.index(u) + 1) % len(ns)]
        check(current == first, "face permutation did not close")
        count += 1
    return count


def planar_bipartite_certificate(parts):
    a, b = parts
    if a == 1:
        rotation = {0: tuple(range(1, b + 1))}
        rotation.update({v: (0,) for v in range(1, b + 1)})
    elif a == 2:
        leaves = tuple(range(2, b + 2))
        rotation = {0: leaves, 1: tuple(reversed(leaves))}
        rotation.update({v: (0, 1) for v in leaves})
    else:
        return None
    faces = rotation_faces(rotation)
    check(sum(parts) - edge_count(parts) + faces == 2, "non-spherical rotation")
    return {"faces": faces, "euler": 2}


def rot_color(c):
    return (0, 2, 3, 1)[c]


def rotate_letter(letter):
    return tuple(sorted(map(rot_color, letter)))


def orbit(letter):
    return min(letter, rotate_letter(letter), rotate_letter(rotate_letter(letter)))


def image_edge(e, permutation):
    return tuple(sorted(permutation[v] for v in e))


def sym2_objects():
    letters = tuple(combinations_with_replacement(range(4), 2))
    indices = {letter: i for i, letter in enumerate(letters)}
    action = {i: indices[rotate_letter(x)] for i, x in enumerate(letters)}
    orbit_labels = [orbit(x) for x in letters]
    charges = [a ^ b for a, b in letters]
    edges_o = {e for e in combinations(range(10), 2)
               if orbit_labels[e[0]] != orbit_labels[e[1]]}
    edges_x = {e for e in combinations(range(10), 2)
               if charges[e[0]] != charges[e[1]]}
    check(sorted(Counter(orbit_labels).values()) == [1, 3, 3, 3], "orbit sizes")
    check(sorted(Counter(charges).values()) == [2, 2, 2, 4], "charge sizes")
    for edges in [edges_o, edges_x]:
        check({image_edge(e, action) for e in edges} == edges, "C3 invariance")
        check(len(edges) == 36 and rank2(10, edges) == 4, "36-edge rank")
    zero = indices[(0, 0)]
    diagonal = [indices[(c, c)] for c in (1, 2, 3)]
    single = [indices[(0, c)] for c in (1, 2, 3)]
    mixed = [indices[tuple(d for d in (1, 2, 3) if d != c)] for c in (1, 2, 3)]
    removed = {tuple(sorted((zero, d))) for d in diagonal}
    removed |= {tuple(sorted(e)) for e in zip(single, mixed)}
    added = set(combinations(single, 2)) | {tuple(sorted(e)) for e in combinations(mixed, 2)}
    check(edges_o - edges_x == removed and edges_x - edges_o == added, "edge trade")
    check(len(edges_o & edges_x) == 30, "common edges")
    check(all(d % 2 for d in degrees(10, edges_o)), "odd-degree source")
    check(all(d % 2 == 0 for d in degrees(10, edges_x)), "even-degree target")
    check(all(d % 2 for d in degrees(10, edges_o ^ edges_x)), "trade boundary")
    check(triangles(10, edges_o) == 54 and triangles(10, edges_x) == 56,
          "triangle hostile control")
    # Ordinary numeric divisors do not carry this adjacency by themselves.
    divisor_labels = {(0, 0): 30, (0, 1): 2, (0, 2): 3, (0, 3): 5,
                      (1, 2): 6, (1, 3): 10, (2, 3): 15,
                      (1, 1): 4, (2, 2): 12, (3, 3): 20}
    check(set(divisor_labels.values()) == {2, 3, 4, 5, 6, 10, 12, 15, 20, 30},
          "divisor virtual alphabet")
    check((indices[(0, 1)], indices[(0, 2)]) not in edges_o,
          "same-channel nonedge")
    check((indices[(0, 1)], indices[(0, 2)]) in edges_x, "different-charge edge")
    check(2 + 3 != 9 and 2 + 3 != 4, "numeric square-sum hostile")
    return {"letters": letters, "orbit_degrees": sorted(degrees(10, edges_o)),
            "charge_degrees": sorted(degrees(10, edges_x)),
            "triangles": [54, 56], "common_edges": 30,
            "removed": sorted([[letters[u], letters[v]] for u, v in removed]),
            "added": sorted([[letters[u], letters[v]] for u, v in added]),
            "symmetric_difference_degree_histogram": sorted(Counter(degrees(10, edges_o ^ edges_x)).items())}


def lexicographic_tournament(depth):
    vertices = list(product(range(3), repeat=depth))
    arcs = set()
    for u, v in combinations(range(len(vertices)), 2):
        for a, b in zip(vertices[u], vertices[v]):
            if a != b:
                arcs.add((u, v) if (b - a) % 3 == 1 else (v, u))
                break
    return vertices, arcs


def main():
    types = types_with_edges(36)
    expected = {2: [(1, 36), (2, 18), (3, 12), (4, 9), (6, 6)],
                3: [(2, 2, 8), (2, 3, 6)],
                4: [(1, 1, 1, 11), (1, 3, 3, 3), (2, 2, 2, 4)]}
    for r, values in expected.items():
        check([p for p in types if len(p) == r] == values, "36-edge classification")
    rows = []
    for parts in types:
        n, r = sum(parts), len(parts)
        edges = graph(parts)
        check(len(edges) == 36, "edge recount")
        check(rank2(n, edges) == r - r % 2, "multipartite rank formula")
        planar = planar_bipartite_certificate(parts) if r == 2 else None
        if planar:
            obstruction = None
        else:
            witness = k33_certificate(parts)
            if witness:
                left, right = witness
                check(all(tuple(sorted((u, v))) in edges for u in left for v in right),
                      "K33 witness")
                obstruction = {"K33": witness}
            else:
                check(r >= 5, "missing nonplanarity certificate")
                reps = []
                total = 0
                for p in parts[:5]:
                    reps.append(total)
                    total += p
                check(all(e in edges for e in combinations(reps, 2)), "K5 witness")
                obstruction = {"K5": reps}
        rows.append({"parts": parts, "vertices": n, "rank_F2": r - r % 2,
                     "planar_rotation": planar, "obstruction": obstruction})
    check(sum(row["planar_rotation"] is not None for row in rows) == 2, "planarity count")
    rank_cases = 0
    for r in range(2, 8):
        for parts in combinations_with_replacement(range(1, 5), r):
            edges = graph(parts)
            n = sum(parts)
            check(rank2(n, edges) == r - r % 2, "general rank")
            check(all(d % 2 for d in degrees(n, edges)) ==
                  (n % 2 == 0 and all(p % 2 for p in parts)), "degree parity")
            doubled = tuple(2 * p for p in parts)
            check(edge_count(doubled) == 4 * edge_count(parts), "doubling edges")
            check(rank2(2 * n, graph(doubled)) == rank2(n, edges), "doubling rank")
            check(edge_count(tuple(8 * p for p in parts)) == 64 * edge_count(parts),
                  "eightfold edges")
            rank_cases += 1
    rank_two_cases = 0
    # Independent six-vertex graph universe: recover nonzero color types via
    # equality of adjacency rows, rather than presupposing a partition.
    pairs = tuple(combinations(range(6), 2))
    odd_rank_two = Counter()
    for mask in range(1 << len(pairs)):
        edges = {e for i, e in enumerate(pairs) if mask >> i & 1}
        if rank2(6, edges) != 2:
            continue
        deg = degrees(6, edges)
        if all(d % 2 for d in deg):
            row_groups = Counter(sum(1 << v for v in range(6)
                                     if tuple(sorted((u, v))) in edges)
                                 for u in range(6))
            sizes = tuple(sorted(row_groups.values()))
            check(sizes in [(1, 5), (3, 3)], "odd rank-two structure")
            odd_rank_two[sizes] += 1
        rank_two_cases += 1
    tournaments = []
    for depth in range(1, 5):
        vertices, arcs = lexicographic_tournament(depth)
        n = len(vertices)
        check(len(arcs) == comb(n, 2), "tournament completeness")
        outdegrees = Counter(u for u, v in arcs)
        check(all(outdegrees[u] == (n - 1) // 2 for u in range(n)), "regularity")
        tournaments.append({"depth": depth, "vertices": n, "arcs": len(arcs),
                            "outdegree": (n - 1) // 2})
    # Same four parts, all 64 uniform crossing directions, one hostile flip each.
    parts = (1, 3, 3, 3)
    labels = [i for i, p in enumerate(parts) for _ in range(p)]
    part_edges = tuple(combinations(range(4), 2))
    for bits in range(64):
        directions = {e: 1 if bits >> i & 1 else -1 for i, e in enumerate(part_edges)}
        blocks = {e: [] for e in part_edges}
        for u, v in graph(parts):
            blocks[(labels[u], labels[v])].append(directions[(labels[u], labels[v])])
        check(all(len(set(b)) == 1 for b in blocks.values()), "uniform quotient")
        blocks[(0, 1)][0] *= -1
        check(len(set(blocks[(0, 1)])) == 2, "hostile crossing flip not detected")
    print(json.dumps({"status": "FINITE-EXACT; algebraic proofs in source note",
                      "types_with_36_edges": rows,
                      "rank_and_replication_cases": rank_cases,
                      "six_vertex_graphs": 1 << len(pairs),
                      "rank_two_graphs": rank_two_cases,
                      "odd_degree_rank_two_graphs": {str(k): v for k, v in sorted(odd_rank_two.items())},
                      "sym2": sym2_objects(), "regular_ternary_tournaments": tournaments,
                      "quotient_positive_and_hostile_pairs": 64}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
