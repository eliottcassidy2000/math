"""Zero-based ranking of tournament isotypes with a prescribed Hamiltonian path count.

The alphabet consists of SCCs in increasing (order, canonical arc code) order,
including the singleton. Ranks are lexicographic in the ordered component word.
With no --order, counts and unranking use the smallest order; H=1 includes the
order-zero unit. Ranking supplied rows or words defaults to their actual order.
"""
import argparse
import json
from collections import Counter

from fibers import bounded_classes, canon, catalog, paths, strong_cap
from solve import validate


def key(rows):
    validate(rows)
    n = len(rows)
    return n, canon.canonical(n, canon.encode(rows)).code


def components(rows):
    """Sorted score prefixes with no arc leaving them are exactly the SCC suffixes."""
    if not rows:
        return []
    validate(rows)
    vertices = sorted(range(len(rows)), key=lambda v: rows[v].bit_count())
    end = total = 0
    blocks = []
    for k, v in enumerate(vertices, 1):
        total += rows[v].bit_count()
        if total == k * (k - 1) // 2:
            block = vertices[end:k]
            blocks.append(canon.rows(len(block), canon.encode(rows, block)))
            end = k
    return blocks[::-1]


def representative(word):
    n = sum(map(len, word))
    result = []
    offset = 0
    for block in word:
        later = (1 << n) - (1 << (offset + len(block)))
        result.extend((r << offset) | later for r in block)
        offset += len(block)
    return tuple(result)


class Fiber:
    def __init__(self, height, order=None, max_core_order=8):
        if height < 1 or height % 2 == 0:
            raise ValueError('height must be a positive odd integer')
        if (order is not None and order < 0) or max_core_order < 0:
            raise ValueError('order and the strong enumeration cap must be nonnegative')
        cap = strong_cap(height)
        needed = cap if order is None else min(cap, order)
        if order is not None and needed > max_core_order:
            raise ValueError(f'order {order} requires --max-core-order {needed} for a complete fiber')
        stop = min(needed, max_core_order)
        self.height = height
        self.alphabet = {(1, 0): 1}
        for n, rows, d in bounded_classes(stop, height):
            if height % d == 0:
                self.alphabet[key(rows)] = d
        self.alphabet = dict(sorted(self.alphabet.items()))
        kinds = Counter((n, d) for (n, _), d in self.alphabet.items())
        # Only quotients reachable using available nontrivial component heights matter.
        divisors = {1, height}
        frontier = [height]
        heights = {d for _, d in kinds if d > 1}
        while frontier:
            t = frontier.pop()
            for d in heights:
                if t % d == 0 and t // d not in divisors:
                    divisors.add(t // d)
                    frontier.append(t // d)
        choices = {t: [(n, t // d, mult) for (n, d), mult in kinds.items() if t % d == 0]
                   for t in divisors}
        if order is None:
            minima = {1: 0}
            for t in sorted(divisors - {1}):
                costs = [n + minima[q] for n, q, _ in choices[t] if q < t and q in minima]
                if costs:
                    minima[t] = min(costs)
            order = minima.get(height)
            if stop < cap and (order is None or order > stop):
                raise ValueError(f'birth unresolved through order {stop}; raise --max-core-order (full cap {cap})')
        self.order = order
        self.dp = [{t: int(t == 1) for t in divisors}]
        for size in range(1, (order or 0) + 1):
            self.dp.append({t: sum(mult * self.dp[size - n][q] for n, q, mult in choices[t] if n <= size)
                            for t in divisors})
        self.count = 0 if order is None else self.dp[order][height]

    def branches(self, remaining, target):
        for component, d in self.alphabet.items():
            n, _ = component
            if n > remaining:
                break
            if target % d == 0:
                q = target // d
                yield component, remaining - n, q, self.dp[remaining - n][q]

    def unrank(self, index):
        if not 0 <= index < self.count:
            raise ValueError(f'rank must be in [0, {self.count})')
        remaining, target = self.order, self.height
        word = []
        while remaining:
            for component, n, q, mass in self.branches(remaining, target):
                if index < mass:
                    word.append(canon.rows(*component))
                    remaining, target = n, q
                    break
                index -= mass
        return word

    def rank(self, word):
        if sum(map(len, word)) != self.order or not self.count:
            raise ValueError('component word does not belong to this fiber')
        remaining, target = self.order, self.height
        index = 0
        for rows in word:
            requested = key(rows)
            for component, n, q, mass in self.branches(remaining, target):
                if component == requested and mass:
                    remaining, target = n, q
                    break
                index += mass
            else:
                raise ValueError('component word does not belong to this fiber')
        return index


def check():
    for height, order, expected in ((1, None, 1), (1, 10, 1), (7, None, 0),
                                    (17, None, 2), (3, 6, 4), (15, 8, 20)):
        fiber = Fiber(height, order)
        assert fiber.count == expected
        for i in range(fiber.count):
            word = fiber.unrank(i)
            assert fiber.rank(word) == i
            rows = representative(word)
            permuted = canon.rows(len(rows), canon.encode(rows, tuple(reversed(range(len(rows))))))
            assert fiber.rank(components(permuted)) == i
    assert Fiber(3, 100).count == 98
    for n in range(3, 7):
        census = {}
        for rows in catalog(n):
            census.setdefault(paths(rows), []).append(rows)
        for height, classes in census.items():
            fiber = Fiber(height, n)
            assert sorted(fiber.rank(components(rows)) for rows in classes) == list(range(fiber.count))
    print('Round trips passed; every catalog class at orders 3 through 6 has a distinct exhaustive rank.')


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--height', type=int)
    cli.add_argument('--order', type=int, help='exact order, including singleton padding')
    cli.add_argument('--max-core-order', type=int, default=8, help='strong enumeration cap (default 8)')
    cli.add_argument('--rows', action='store_true', help='also emit the joined representative rows')
    action = cli.add_mutually_exclusive_group()
    action.add_argument('--unrank', type=int, metavar='K')
    action.add_argument('--rank-word', type=json.loads, metavar='JSON', help='ordered list of strong row arrays')
    action.add_argument('--rank-rows', type=json.loads, metavar='JSON', help='one tournament row array in any labelling')
    action.add_argument('--check', action='store_true')
    args = cli.parse_args()
    if args.check:
        check()
        return
    if args.height is None:
        cli.error('provide --height')
    try:
        word = args.rank_word if args.rank_rows is None else components(args.rank_rows)
        order = args.order if args.order is not None or word is None else sum(map(len, word))
        fiber = Fiber(args.height, order, args.max_core_order)
        output = {'height': args.height, 'order': fiber.order, 'classes': fiber.count}
        index = fiber.rank(word) if word is not None else args.unrank
        if index is not None:
            word = fiber.unrank(index)
            output.update(rank=index, word=word)
            if args.rows:
                output['rows'] = representative(word)
        print(json.dumps(output))
    except (ValueError, TypeError) as err:
        cli.error(str(err))


if __name__ == '__main__':
    main()
