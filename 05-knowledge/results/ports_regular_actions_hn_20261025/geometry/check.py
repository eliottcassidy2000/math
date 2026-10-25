"""Exact finite controls for the plane cube theorem. Standard library only."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def binary_rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def row_basis(rows):
    """Literal rational elimination; independent of the Walsh diagonalization."""
    pivots = {}
    for source in rows:
        row = list(map(Fraction, source))
        for p, old in sorted(pivots.items()):
            if row[p]:
                scale = row[p]
                row = [x - scale * y for x, y in zip(row, old)]
        nonzero = next((i for i, x in enumerate(row) if x), None)
        if nonzero is not None:
            scale = row[nonzero]
            pivots[nonzero] = [x / scale for x in row]
    return pivots


def span_contains(basis, source):
    row = list(map(Fraction, source))
    for p, old in sorted(basis.items()):
        scale = row[p]
        row = [x - scale * y for x, y in zip(row, old)]
    return not any(row)


def colour_number(directions):
    adjacency = [{x ^ s for s in directions} for x in range(8)]
    labels = [-1] * 8

    def extend(k, assigned):
        if assigned == 8:
            return True
        todo = [v for v in range(8) if labels[v] < 0]
        v = max(todo, key=lambda v: (
            len({labels[w] for w in adjacency[v] if labels[w] >= 0}),
            len(adjacency[v]), -v))
        forbidden = {labels[w] for w in adjacency[v]}
        for c in range(k):
            if c not in forbidden:
                labels[v] = c
                if extend(k, assigned + 1):
                    return True
        labels[v] = -1
        return False

    for k in range(1, 9):
        labels[:] = [-1] * 8
        if extend(k, 0):
            return k
    raise RuntimeError("Eight vertices require more than eight colours")


def scalar_product(a, b):
    return (a & b).bit_count() & 1


def direction_sets():
    colour_histogram = Counter()
    dimension_histogram = Counter()
    independent_count = 0
    collapsed_checks = 0
    for mask in range(128):
        directions = [i + 1 for i in range(7) if mask >> i & 1]
        rank = binary_rank(directions)
        coloops = [a for a in directions
                   if binary_rank([b for b in directions if b != a]) < rank]
        rows = []
        for a, b in combinations(directions, 2):
            for x in range(8):
                row = [0] * 8
                for y, sign in ((x, 1), (x ^ a ^ b, 1),
                                (x ^ a, -1), (x ^ b, -1)):
                    row[y] += sign
                rows.append(row)
        basis = row_basis(rows)
        expected = 2 ** (3 - rank) * (1 + len(coloops))
        require(8 - len(basis) == expected, (mask, "square nullity"))
        allowed = [a for a in range(8)
                   if sum(scalar_product(a, s) for s in directions) <= 1]
        require(len(allowed) == expected, (mask, "Walsh nullity"))
        for a in range(8):
            char = [(-1) ** scalar_product(a, x) for x in range(8)]
            actual = all(sum(x * y for x, y in zip(row, char)) == 0
                         for row in rows)
            require(actual == (a in allowed), (mask, a, "Walsh membership"))
        for a in set(directions) - set(coloops):
            for x in range(8):
                row = [0] * 8
                row[x ^ a], row[x] = 1, -1
                require(span_contains(basis, row), (mask, a, "collapse"))
                collapsed_checks += 1
        chromatic = colour_number(directions)
        line = any((a ^ b) in directions for a, b in combinations(directions, 2))
        predicted = 1 if not directions else 8 if len(directions) == 7 else 4 if line else 2
        require(chromatic == predicted, (mask, "chromatic", chromatic))
        colour_histogram[chromatic] += 1
        dimension_histogram[expected] += 1
        independent_count += rank == len(directions)
    require(colour_histogram == {1: 1, 2: 63, 4: 63, 8: 1}, "colour census")
    require(independent_count == 57, "independent direction census")
    return {"sets": 128, "chromatic_histogram": dict(colour_histogram),
            "square_flat_dimension_histogram": dict(dimension_histogram),
            "embedding_eligible_sets": independent_count,
            "noncoloop_collapse_rows": collapsed_checks}


def rational_cube():
    units = [(Fraction(1), Fraction(0)), (Fraction(3, 5), Fraction(4, 5)),
             (Fraction(5, 13), Fraction(12, 13))]
    points = [tuple(sum(((-1) ** ((x >> j) & 1)) * units[j][axis] / 2
                        for j in range(3)) for axis in range(2)) for x in range(8)]
    require(len(set(points)) == 8, "cube injectivity")
    edges = 0
    for a, b in combinations(range(8), 2):
        dist2 = sum((points[a][i] - points[b][i]) ** 2 for i in range(2))
        require((dist2 == 1) == ((a ^ b).bit_count() == 1), "induced cube")
        edges += dist2 == 1
    matrix = [[-2 if x == y else 2 if (x ^ y).bit_count() == 1 else 0
               for y in range(8)] for x in range(8)]
    for axis in range(2):
        require(all(sum(matrix[x][y] * points[y][axis] for y in range(8)) == 0
                    for x in range(8)), "coordinate kernel")
    spectrum = Counter()
    for a in range(8):
        char = [(-1) ** scalar_product(a, x) for x in range(8)]
        eigenvalue = 4 - 4 * a.bit_count()
        require(all(sum(matrix[x][y] * char[y] for y in range(8))
                    == eigenvalue * char[x] for x in range(8)), "K3 spectrum")
        spectrum[eigenvalue] += 1
    # Dropping injectivity defeats the square lemma, even though every edge is unit.
    collapse = [x.bit_count() % 2 for x in range(8)]
    require(all(abs(collapse[x] - collapse[x ^ s]) == 1
                for x in range(8) for s in (1, 2, 4, 7)), "noninjective unit control")
    require(collapse[0] + collapse[3] != collapse[1] + collapse[2],
            "noninjective square hostile")
    tetra = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    require(all(sum((a[j] - b[j]) ** 2 for j in range(3)) == 8
                for a, b in combinations(tetra, 2)), "unit tetrahedron after rescaling")
    require(any(tetra[0][j] + tetra[2][j] != tetra[1][j] + tetra[3][j]
                for j in range(3)), "dimension hostile")
    return {"vertices": 8, "unit_edges": edges, "distance_checks": 28,
            "intersection_spectrum": dict(spectrum), "hostile_controls": 2}


def fano_actions():
    histogram = Counter()
    for labels in permutations(range(1, 8)):
        c = labels[0]
        differences = [labels[1] ^ labels[2], labels[3] ^ labels[4],
                       labels[5] ^ labels[6]]
        require(differences[0] ^ differences[1] ^ differences[2] == c,
                "central sum")
        bad = differences.count(c)
        rank = binary_rank(differences)
        require((bad, rank) in ((0, 3), (1, 2), (3, 1)), "rank-crossing relation")
        translations = {sum(scalar_product(a, d) << j for j, d in enumerate(differences))
                        for a in range(8)}
        orbits = {tuple(sorted(x ^ t for t in translations)) for x in range(8)}
        require(len(orbits) == 2 ** (3 - rank), "orbit count")
        require(all(len(o) == 2 ** rank for o in orbits), "orbit size")
        require((len(translations) == 8) == (bad == 0), "free/regular criterion")
        histogram[(bad, rank)] += 1
    require(histogram == {(0, 3): 16 * 168, (1, 2): 12 * 168, (3, 1): 2 * 168},
            "Fano census")
    marked_directions = []
    for quotient in (2, 4):
        pairs = ((1, 6), (quotient, (-quotient) % 7),
                 ((1 + quotient) % 7, (-1 - quotient) % 7))
        for marking in ((1, 2, 4, 3, 6, 7, 5), (1, 2, 4, 5, 7, 3, 6)):
            differences = [marking[x] ^ marking[y] for x, y in pairs]
            require(binary_rank(differences) == 3, "torus-selected regular action")
            marked_directions.append(differences)
    require(marked_directions == [[7, 3, 5], [4, 7, 2], [7, 5, 3], [4, 2, 7]],
            "marked torus directions")
    return {"labellings": 5040, "torus_selected_directions": marked_directions,
            "bad_count_rank_labellings": {f"{b},{r}": n for (b, r), n in histogram.items()},
            "quotients_by_GL3": {f"{b},{r}": n // 168 for (b, r), n in histogram.items()}}


def main():
    print(json.dumps({"status": "PASS", "directions": direction_sets(),
                      "cube": rational_cube(), "fano": fano_actions()}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
