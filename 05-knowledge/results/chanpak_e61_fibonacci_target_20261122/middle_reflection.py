"""Contract the single flipped arc and enumerate its constrained run coefficients."""
from collections import Counter
from pathlib import Path
import sys

RESULTS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(RESULTS / "chanpak_support_birth_grammar_20261122"))
from check_perfect_port import H, LOWER, LPORTS, inflate


def contracted_skeletons(first, last):
    atom = 8
    # Vertex 7 belongs to the atom; the other vertices of port 6 still occur.
    labels = tuple(v for v in range(9) if v != 7)
    caps = {v: 3 if v in LPORTS else 1 for v in labels}
    counts = [0] * 9
    result = Counter()

    def contact(i, j):
        if i == atom:
            return "above" if j == last else bool(LOWER[last] >> j & 1)
        if j == atom:
            return "below" if i == first else bool(LOWER[i] >> first & 1)
        return bool(LOWER[i] >> j & 1)

    def extend(i, condition):
        if all(counts[v] for v in labels):
            result[tuple(counts[v] for v in LPORTS), condition] += 1
        for j in labels:
            if i == j or counts[j] == caps[j]:
                continue
            edge = contact(i, j)
            if edge:
                counts[j] += 1
                extend(j, edge if isinstance(edge, str) else condition)
                counts[j] -= 1

    for i in labels:
        counts[i] = 1
        extend(i, None)
        counts[i] = 0
    return result


def midpoint_derivative(port_width, left_width, right_width, rank):
    y, z = (1 << left_width) - 2, (1 << right_width) - 2
    b = 15 + 11*z + 5*y + 4*y*z
    d = 10 + 8*z + 5*y + 4*y*z
    remaining = port_width - 1
    if rank == 0:
        return -(1 << remaining) * b
    if rank == remaining:
        return ((1 << remaining) - 2) * b + d
    return b * ((1 << rank) - (1 << (remaining - rank)))


def main():
    delta = contracted_skeletons(6, 7)
    delta.subtract(contracted_skeletons(7, 6))
    nonzero = {key: value for key, value in delta.items() if value}
    expected = {}
    for j, k, b, d in ((1, 1, 15, 10), (1, 2, 11, 8), (2, 1, 5, 5), (2, 2, 4, 4)):
        expected[(1, j, k), "above"] = -2*b
        expected[(1, j, k), "below"] = d
        expected[(2, j, k), "above"] = -b
        expected[(2, j, k), "below"] = b
    assert nonzero == expected
    for widths in ((3, 2, 2), (4, 1, 2), (5, 1, 1), (7, 2, 2)):
        rows = inflate(LOWER, LPORTS, widths)
        base = H(rows)
        blocks = dict(zip(LPORTS, widths))
        port_start = sum(blocks.get(i, 1) for i in range(6))
        vertex7 = len(rows) - 1
        for rank in range(widths[0]):
            v = port_start + rank
            flipped = list(rows)
            flipped[v] ^= 1 << vertex7
            flipped[vertex7] ^= 1 << v
            assert H(flipped) - base == midpoint_derivative(*widths, rank)
    print("16 signed coefficients, four inflation profiles, all tested ranks exact")


if __name__ == "__main__":
    main()
