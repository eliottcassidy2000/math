"""Replay the mixed difference and positive unit-square obstruction."""
from functools import lru_cache
from itertools import combinations, permutations, product


def orient(rows, i, j):
    rows[i] |= 1 << j
    rows[j] &= ~(1 << i)


def completions(board, holes):
    values = []
    for mask in range(1 << len(holes)):
        rows = list(board)
        for k, (i, j) in enumerate(holes):
            source, target = (i, j) if mask >> k & 1 else (j, i)
            orient(rows, source, target)
        values.append(height(tuple(rows)))
    return values


@lru_cache(maxsize=30000)
def height(rows):
    @lru_cache(None)
    def follow(i, todo):
        if not todo:
            return 1
        return sum(follow(j, todo ^ (1 << j)) for j in range(len(rows))
                   if (rows[i] & todo) >> j & 1)
    all_bits = (1 << len(rows)) - 1
    return sum(follow(i, all_bits ^ (1 << i)) for i in range(len(rows)))


def mixed(heights):
    return heights[3] - heights[1] - heights[2] + heights[0]


def disjoint():
    for r in range(3):
        cases = units = 0
        outside = list(combinations(range(r), 2))
        for bits in product(range(2), repeat=2*r + len(outside)):
            board = [0] * (4 + r)
            host = [0] * (2 + r)
            for k in range(r):
                for group in range(2):
                    if bits[2*k + group]:
                        orient(host, group, k+2)
                        for i in range(2*group, 2*group+2):
                            orient(board, i, k+4)
                    else:
                        orient(host, k+2, group)
                        for i in range(2*group, 2*group+2):
                            orient(board, k+4, i)
            for edge, flag in zip(outside, bits[2*r:]):
                i, j = edge if flag else edge[::-1]
                orient(host, i+2, j+2)
                orient(board, i+4, j+4)
            minus, plus = completions(host, [(0, 1)])
            for cross in product(range(2), repeat=4):
                ac, ad, bc, bd = cross
                for (i, j), bit in zip(product((0, 1), (2, 3)), cross):
                    source, target = (i, j) if bit else (j, i)
                    orient(board, source, target)
                h = completions(board, [(0, 1), (2, 3)])
                assert mixed(h) == (bc-bd-ac+ad)*(plus-minus)
                unit = mixed(h) == 0 and abs(h[1]-h[0]) == abs(h[2]-h[0]) == 2
                padding = plus == minus == 1 and sum(cross) in (1, 3)
                assert unit == padding
                if unit:
                    assert sorted(h) == [1, 3, 3, 5]
                cases += 1
                units += unit
        print('disjoint', 'order', r+4, 'marked exteriors and cross matrices',
              cases, 'balanced unit squares', units)


def incident():
    for n in range(3, 6):
        # Relabeling fixes the apex at 0 and its two marked neighbors at 1 and 2.
        other = [e for e in combinations(range(n), 2) if e not in ((0, 1), (0, 2))]
        checked = 0
        for flags in product(range(2), repeat=len(other)):
            board = [0] * n
            for (i, j), bit in zip(other, flags):
                source, target = (i, j) if bit else (j, i)
                orient(board, source, target)
            h = completions(board, [(0, 1), (0, 2)])
            mass = sum(all(a == 1 and b == 2 or a == 2 and b == 1 or board[a] >> b & 1
                           for a, b in zip(w, w[1:]))
                       and abs(w.index(1)-w.index(2)) == 1
                       for w in permutations(range(1, n)))
            assert mixed(h) == -mass <= -2 and mass % 2 == 0
            checked += 1
        print('incident', 'order', n, 'completions of fixed marked holes', checked)


if __name__ == '__main__':
    disjoint()
    incident()
