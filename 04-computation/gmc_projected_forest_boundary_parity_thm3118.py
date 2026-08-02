#!/usr/bin/env python3
"""Exact controls for THM-3118's projected forest-boundary parity law."""

from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rgs_partitions(n, blocks):
    """Restricted-growth strings for set partitions with exactly blocks blocks."""
    if n == 0:
        if blocks == 0:
            yield ()
        return
    word = [0] * n

    def rec(i, top):
        if i == n:
            if top + 1 == blocks:
                yield tuple(word)
            return
        # Even using a new label at every remaining position must reach blocks.
        if top + 1 > blocks or top + 1 + (n - i) < blocks:
            return
        for value in range(min(top + 1, blocks - 1) + 1):
            word[i] = value
            yield from rec(i + 1, max(top, value))

    yield from rec(1, 0)


def component_word(n, edges, omitted=None):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for j, (a, b) in enumerate(edges):
        if j == omitted:
            continue
        ra, rb = find(a), find(b)
        if ra == rb:
            return None
        parent[rb] = ra

    labels = {}
    answer = []
    for x in range(n):
        root = find(x)
        if root not in labels:
            labels[root] = len(labels)
        answer.append(labels[root])
    return tuple(answer)


def projected_rank(n, r):
    rows = list(rgs_partitions(n, n - r))
    row_index = {row: i for i, row in enumerate(rows)}
    target = len(rows)
    expected = target if r % 2 == 0 else target - 1
    edges = list(combinations(range(n), 2))
    pivots = {}
    read = 0

    for edge_ids in combinations(range(len(edges)), r + 1):
        forest = tuple(edges[j] for j in edge_ids)
        if component_word(n, forest) is None:
            continue
        read += 1
        column = 0
        for omitted in range(r + 1):
            row = component_word(n, forest, omitted)
            column ^= 1 << row_index[row]
        require(column.bit_count() == r + 1,
                "forest deletion flats were not distinct")
        if r % 2:
            require(column.bit_count() % 2 == 0,
                    "odd-rank column lost its augmentation parity")

        reduced = column
        while reduced:
            pivot = reduced.bit_length() - 1
            old = pivots.get(pivot)
            if old is None:
                pivots[pivot] = reduced
                break
            reduced ^= old
        if len(pivots) == expected:
            break

    require(len(pivots) == expected,
            f"rank mismatch at {(n, r, target, len(pivots))}")
    return target, expected, read


def cut_key(mask, m):
    full = (1 << m) - 1
    other = full ^ mask
    return min(mask, other)


def tree_splits(m, edges):
    answer = []
    for omitted in range(len(edges)):
        parent = list(range(m))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for j, (a, b) in enumerate(edges):
            if j == omitted:
                continue
            ra, rb = find(a), find(b)
            require(ra != rb, "double-star control contained a cycle")
            parent[rb] = ra
        root0 = find(0)
        mask = sum(1 << x for x in range(m) if find(x) == root0)
        answer.append(cut_key(mask, m))
    return sorted(answer)


def double_star_controls(m):
    full = (1 << m) - 1
    checked = 0
    for mask in range(1, full):
        complement = full ^ mask
        if mask > complement:
            continue
        a = next(x for x in range(m) if mask & (1 << x))
        b = next(x for x in range(m) if complement & (1 << x))
        edges = [(a, b)]
        edges += [(a, x) for x in range(m) if x != a and mask & (1 << x)]
        edges += [(b, x) for x in range(m) if x != b and complement & (1 << x)]
        require(len(edges) == m - 1,
                "double-star control has the wrong edge count")
        expected = [cut_key(mask, m)]
        expected += [cut_key(1 << x, m) for x in range(m) if x not in (a, b)]
        require(tree_splits(m, edges) == sorted(expected),
                "double-star fundamental cuts do not match the lemma")
        checked += 1
    return checked


def main():
    cases = [
        (3, 0), (4, 0),
        (3, 1), (4, 1), (5, 1),
        (4, 2), (5, 2), (6, 2),
        (5, 3), (6, 3), (7, 3),
        (6, 4), (7, 4), (8, 4),
        (7, 5), (8, 5),
    ]
    print("projected_boundary_rank_over_F2")
    print("N r rows predicted_rank exact_rank columns_read")
    for n, r in cases:
        rows, rank, read = projected_rank(n, r)
        print(n, r, rows, rank, rank, read)

    print("double_star_controls")
    total = 0
    for m in range(3, 9):
        count = double_star_controls(m)
        total += count
        print(m, count)
    print("total", total)
    print("all_exact_checks_pass")


if __name__ == "__main__":
    main()
