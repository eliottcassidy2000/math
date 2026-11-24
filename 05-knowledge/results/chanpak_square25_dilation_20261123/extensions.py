"""Replay the twelve fixed marked tournaments for all orders, using integer DP."""

from math import comb

# The same four labelled ten-vertex tournaments supply all the contact boards.
BOARDS = {
    "A": (766, 1020, 1016, 496, 992, 960, 896, 768, 513, 8),
    "B": (542, 1020, 1016, 1008, 864, 705, 385, 785, 33, 320),
    "C": (750, 924, 1016, 1008, 737, 322, 898, 800, 17, 288),
    "D": (1022, 1004, 216, 880, 930, 964, 912, 776, 516, 4),
}

# Board, outgoing mask of vertex 10, marked old vertex, nonzero run counts.
RAILS = (
    ("A", 512, 1, (1489, 712)),
    ("D", 1, 8, (15339, 15976, 5382)),
    ("A", 1008, 1, (2261, 1098)),
    ("D", 773, 8, (7783, 4210)),
    ("A", 962, 2, (5625, 4466, 1072)),
    ("D", 1021, 8, (1859, 1026)),
    ("A", 128, 1, (1525, 724)),
    ("B", 8, 1, (9327, 4982)),
    ("A", 1009, 1, (2345, 1230)),
    ("C", 32, 7, (7147, 6390, 1872)),
    ("A", 993, 1, (2541, 1298)),
    ("B", 1014, 1, (6191, 5156, 1218)),
)


def members(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def extend(rows, contact):
    new = 1 << len(rows)
    return tuple(row if contact & (1 << i) else row | new
                 for i, row in enumerate(rows)) + (contact,)


def reaches_all(rows):
    reached = 1
    todo = 1
    while todo:
        new = 0
        for i in members(todo):
            new |= rows[i]
        todo = new & ~reached
        reached |= todo
    return reached == (1 << len(rows)) - 1


def precedence_order(rows):
    """Parity of exact prefix/suffix path decompositions at each vertex."""
    full = (1 << len(rows)) - 1
    ends = [0] * (full + 1)
    starts = [0] * (full + 1)
    for mask in range(1, full + 1):
        for i in members(mask):
            bit = 1 << i
            rest = mask ^ bit
            if not rest or (ends[rest] & ~rows[i]).bit_count() % 2:
                ends[mask] |= bit
            if not rest or (starts[rest] & rows[i]).bit_count() % 2:
                starts[mask] |= bit
    predecessors = [0] * len(rows)
    for prefix in range(1, full + 1):
        suffix = full ^ prefix
        for j in members(ends[prefix]):
            if not suffix or (starts[suffix] & rows[j]).bit_count() % 2:
                predecessors[j] ^= prefix ^ (1 << j)
    order = sorted(range(len(rows)), key=lambda i: predecessors[i].bit_count())
    seen = 0
    for i in order:
        assert predecessors[i] == seen
        seen |= 1 << i
    return order


def triangle_matching(rows, mark):
    """Largest number of directed triangles disjoint outside the mark."""
    ins = ((1 << len(rows)) - 1) ^ rows[mark] ^ (1 << mark)
    used = {0}
    for i in members(rows[mark]):
        used |= {mask | (1 << j) for mask in used
                 for j in members(rows[i] & ins & ~mask)}
    return max(mask.bit_count() for mask in used)


def run_counts(rows, mark):
    """c[k] by last letter, exterior subset, and number of visits to mark."""
    n = len(rows)
    exterior = ((1 << n) - 1) ^ (1 << mark)
    masks = [s for s in range(1 << n) if not s & (1 << mark)]
    incoming = [exterior & ~(r | (1 << i)) for i, r in enumerate(rows)]
    previous = None
    coeff = []
    for k in range(n + 1):
        paths = [[0] * n for _ in range(1 << n)]
        for mask in masks:
            if k:
                paths[mask][mark] = int(k == 1 and not mask) + sum(
                    previous[mask][j] for j in members(mask & incoming[mark]))
            for i in members(mask):
                rest = mask ^ (1 << i)
                value = int(not k and not rest)
                if rows[mark] & (1 << i):
                    value += paths[rest][mark]
                value += sum(paths[rest][j] for j in members(rest & incoming[i]))
                paths[mask][i] = value
        coeff.append(sum(paths[exterior]))
        previous = paths
    while coeff and coeff[-1] == 0:
        coeff.pop()
    return tuple(coeff[1:])


def height(coeff, width):
    return sum(c * sum((-1) ** (k-j) * comb(k, j) * j ** width
                       for j in range(k+1))
               for k, c in enumerate(coeff, 1))


def choose(order, target):
    """Return the child construction and exact height for a residue and order."""
    if order < 11 or target % 2 == 0:
        raise ValueError("Choose an order at least 11 and an odd target.")
    width = order - 10
    phase = width if width <= 2 else 3 + (width % 2 == 0)
    for board, contact, mark, coeff in RAILS:
        if height(coeff, phase) % 24 == target % 24:
            return board, contact, mark, width, height(coeff, width)
    raise AssertionError("incomplete selector")


def main():
    full = (1 << 11) - 1
    columns = [set() for _ in range(4)]
    for board, contact, mark, coeff in RAILS:
        rows = extend(BOARDS[board], contact)
        for i, row in enumerate(rows):
            assert not row & (1 << i)
            for j in range(i):
                assert (row >> j ^ rows[j] >> i) & 1
        assert reaches_all(rows)
        assert reaches_all([full ^ r ^ (1 << i) for i, r in enumerate(rows)])
        order = precedence_order(rows)
        sides = [(rows[mark] >> j) & 1 for j in order if j != mark]
        assert sides == sorted(sides)
        assert run_counts(rows, mark) == coeff
        assert triangle_matching(rows, mark) == len(coeff) - 1
        residues = [height(coeff, w) % 24 for w in range(1, 5)]
        for col, r in zip(columns, residues):
            col.add(r)
        print(board, contact, mark, coeff, residues)
    assert all(c == set(range(1, 24, 2)) for c in columns)
    for n in range(11, 25):
        for r in range(1, 24, 2):
            assert choose(n, r)[-1] % 24 == r
    even = {c[0] % 24: height(c, 2) % 24 for *_, c in RAILS}
    odd = {c[0] % 24: height(c, 3) % 24 for *_, c in RAILS}
    assert all(odd[odd[r]] == r and odd[even[odd[even[r]]]] == r for r in even)
    for orbit in ((1, 9, 13, 21), (3, 11, 23, 15, 19, 7), (5, 17)):
        assert all(even[r] == orbit[(i+1) % len(orbit)] for i,r in enumerate(orbit))
    print("Twelve fixed marked tournaments cover every odd residue at each order ≥ 11.")


if __name__ == "__main__":
    main()
