"""Invert both the normalized Fibonacci receiver and its right address from an odd lower target."""

import argparse
from check import CORE, flipped, graft_chart, ordinal, profile, reader_extensions, strong_readers
from replay import inverse_one_sided


def address_window(target, chart, width):
    """Necessary prefix, trailing ones, and upper endpoint at fixed positive bit width."""
    c, h, k0, k1 = chart
    first = 1 << (width - 1)
    room = target - min(k0, k1) - c * first
    if room < 0:
        return None
    upper = min(2*first - 1, first + room // c)
    leading_zeros = 0
    while leading_zeros + 1 < width and upper - first < first >> (leading_zeros + 1):
        leading_zeros += 1
    required_ones = 0
    while required_ones + 1 < width and first >> (required_ones + 1) > room // h:
        required_ones += 1
    return first, upper, leading_zeros, required_ones


def fixed_size_addresses(target, chart, width):
    if width < 0:
        return
    c, h, k0, k1 = chart
    if width == 0:
        if target == k0:
            yield 0
        return
    window = address_window(target, chart, width)
    if window is None:
        return
    first, last, _, t = window
    scale = 1 << t
    if t:
        intercept = k1 + c * (scale - 1)
        if target == intercept and first <= scale - 1 <= last:
            yield scale - 1
        for parity in (0, 1):
            for z in inverse_one_sided(target, intercept, h, scale*c, parity):
                y = scale*(z+1) - 1
                if first <= y <= last:
                    yield y
    else:
        for parity, k in enumerate((k0, k1)):
            for y in inverse_one_sided(target, k, h, c, parity):
                if first <= y <= last:
                    yield y


def select(target, max_order=None, vertices=None):
    """Yield (rooted receiver, right address) in the six-vertex unit-flip construction."""
    if target < 1 or target % 2 == 0:
        return
    budget = 2 * (target + 10) // 28
    if vertices is not None:
        if vertices < 8 or vertices - 7 > budget.bit_length() - 1:
            return
        root_cap = max(1, 4 * (target + 10) // (28 << (vertices - 7)))
        max_order = root_cap if max_order is None else min(root_cap, max_order)

    def descend(reader):
        data = profile(reader, 1)
        slots = data[1]
        if slots > budget:
            return
        chart = graft_chart(reader, data)
        if vertices is not None:
            for y in fixed_size_addresses(target, chart, vertices - 7 - len(reader)):
                yield reader, y
        else:
            c, h, even, odd = chart
            if target == even:
                yield reader, 0
            for parity, intercept in enumerate((even, odd)):
                for right in inverse_one_sided(target, intercept, h, c, parity):
                    yield reader, right
        extra = 1
        while slots << extra <= budget and (max_order is None or len(reader) + extra <= max_order):
            for parent in reader_extensions(reader, extra - 1):
                yield from descend(parent)
            extra += 1

    yield from descend((0,))


def checks():
    from itertools import product
    from replay import p, star_rows, strong, subpaths
    for target in (189, 407, 873, 945, 1367, 2415, 3759):
        brute = set()
        for reader in strong_readers(5):
            c, h, even, odd = graft_chart(reader)
            for y in range(1 + target // c):
                if (odd if y % 2 else even) + c*y + h*p(y) == target:
                    brute.add((reader, y))
        assert set(select(target, max_order=5)) == brute
    for width, r in product(range(1, 5), (0, 1)):
        for reader in strong_readers(width):
            c, h, even, odd = graft_chart(reader)
            assert (reader, r) in select((odd if r else even) + c*r, max_order=width)
    limit = 3000
    width = ((limit + 10) // 14).bit_length() - 1
    presentations = {}
    for reader in strong_readers(width):
        c, h, even, odd = graft_chart(reader)
        shallow = sorted((row >> 1).bit_count() for row in reader[1:]) == list(range(len(reader) - 1))
        for y in range(1 + limit // c):
            value = (odd if y % 2 else even) + c*y + h*p(y)
            if value <= limit:
                presentations.setdefault(value, []).append((shallow, reader, y))
    shallow_values = {n for n, v in presentations.items() if any(shallow for shallow, _, _ in v)}
    assert (len(presentations), len(shallow_values), min(presentations.keys() - shallow_values)) == (50, 45, 859)
    first_positive = min(n for n, v in presentations.items()
                         if all(not shallow and y > 0 for shallow, _, y in v))
    assert first_positive == 1955
    for target, y in ((859, 0), (1955, 1)):
        reader = (2, 4, 9, 3)
        assert list(select(target)) == [(reader, y)]
        assert list(select(target, vertices=7 + len(reader) + y.bit_length())) == [(reader, y)]
        deletion = ordinal(CORE, reader)
        graphs = [star_rows(g, 1 | (1 << len(CORE)), 0, y)
                  for g in (deletion, flipped(deletion, 1, 3))]
        assert all(map(strong, graphs)) and flipped(graphs[0], 1, 3) == graphs[1]
        assert [subpaths(g)[-1] for g in graphs] == [target, target+2]
    n61 = 2**60 * (2**61 - 1) - 1
    top = []
    for reader in strong_readers(4):
        window = address_window(n61, graft_chart(reader), 116-len(reader))
        if window:
            top.append((reader, *window[2:]))
    assert top == [((0,), 0, 0), ((2, 4, 1), 2, 0), ((2, 12, 9, 1), 8, 5)]
    print('All lower outputs through 3000: 50 for all receivers, 45 with transitive root deletion.')
    print('Receiver recursion agrees with direct chart inversion on the bounded controls.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lower', type=int, help='prescribed odd lower count')
    parser.add_argument('--vertices', type=int, help='optional exact order of the normalized construction')
    parser.add_argument('--max-order', type=int, help='optional maximum normalized receiver order')
    args = parser.parse_args()
    if args.lower is None:
        checks()
    elif args.lower < 1 or args.lower % 2 == 0 or (args.max_order is not None and args.max_order < 1):
        parser.error('choose a positive odd lower count and an optional positive order cap')
    else:
        count = 0
        for reader, y in select(args.lower, args.max_order, args.vertices):
            count += 1
            print('receiver rows', reader, 'tag 0, right address', y, flush=True)
        scope = 'all receiver orders' if args.max_order is None else f'receiver orders through {args.max_order}'
        if args.vertices is not None:
            scope += f', total order {args.vertices}'
        print(count, 'solutions for lower count', args.lower, 'in', scope)
