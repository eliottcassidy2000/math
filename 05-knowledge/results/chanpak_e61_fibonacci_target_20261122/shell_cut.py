"""Choose the next Fibonacci receiver at a specified lower height and right address by minimum cuts."""

import argparse
from itertools import combinations
from pathlib import Path
import sys

from check import extend_reader, graft_chart
from replay import p

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "chanpak_minimal_core_grammar_20261122"))
from grammar import minimum_completion


def height(child, tail, y, mask):
    c, h, k0, k1 = graft_chart(extend_reader(child, tail, mask))
    return (k1 if y % 2 else k0) + c*y + h*p(y)


def search(child, tail, y, target):
    """All choices of the free tail contacts hitting target in this shell."""
    r = max(0, tail - 1)
    base = height(child, tail, y, 0)
    linear = [height(child, tail, y, 1 << i) - base for i in range(r)]
    weights = [[0] * (r+2) for _ in range(r+2)]
    source, sink = r, r+1
    for i, j in combinations(range(r), 2):
        capacity = base + linear[i] + linear[j] - height(child, tail, y, (1 << i) | (1 << j))
        assert capacity > 0
        weights[i][j] = capacity
    offset = base
    for i, value in enumerate(linear):
        unary = value - sum(weights[i])
        if unary >= 0:
            weights[i][sink] = unary
        else:
            weights[source][i] = -unary
            offset += unary

    def descend(i, yes, no):
        cheapest = offset + minimum_completion(weights, yes | (1 << source), no | (1 << sink))
        if cheapest > target:
            return
        if i == r:
            if cheapest == target:
                yield yes, extend_reader(child, tail, yes)
            return
        yield from descend(i+1, yes, no | (1 << i))
        yield from descend(i+1, yes | (1 << i), no)

    yield from descend(0, 0, 0)


def checks():
    for child, tail, y in [((0,), 4, 3), ((2, 4, 1), 4, 5), ((2, 4, 9, 3), 3, 11)]:
        outputs = {}
        for mask in range(1 << (tail-1)):
            outputs.setdefault(height(child, tail, y, mask), set()).add(mask)
        for h, expected in outputs.items():
            assert {mask for mask, _ in search(child, tail, y, h)} == expected
    print('Minimum-cut receiver selection agrees with every contact word in three complete shells.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--child', default='0', help='comma-separated masks of a normalized rooted receiver')
    parser.add_argument('--tail', type=int, default=4)
    parser.add_argument('--right', type=int, default=3)
    parser.add_argument('--lower', type=int)
    args = parser.parse_args()
    if args.lower is None:
        checks()
    else:
        child = tuple(int(s) for s in args.child.split(','))
        if args.tail < 0 or args.right < 0 or (len(child) == 1 and not args.tail):
            parser.error('choose a valid extension and a nonnegative address')
        for mask, receiver in search(child, args.tail, args.right, args.lower):
            print('tail contacts', mask, 'receiver', receiver)
