"""Check the exponent-31 edges and the two independent binary-address formula."""
import argparse
import sys
from itertools import combinations
from math import factorial
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'chanpak_gap_star_ports_20261122'))
from replay import count_from_profile, flipped, inverse_one_sided, p, profile, star_rows, strong, subpaths
from check_perfect_band import remove_common_extremes

from e31_certificates import CORE, NEARBY, WITNESSES, ordinal


def word_reader(v):
    """The tagged apex starts its unique full path on the initial spine vertex."""
    word = bin(v)[2:]
    transitive = tuple((1 << len(word)) - (2 << i) for i in range(len(word)))
    tag = int(word[::-1], 2)
    return star_rows(transitive, tag, 0, 0), 1 << len(word)


def two_slots(v):
    """Assign spine vertices to before-apex, after-apex and the other slot."""
    states = {(True, False): 1}
    for bit in bin(v)[2:]:
        next_states = {}
        for (valid_before, after_used), count in states.items():
            options = [(valid_before, after_used), (bit == '0', after_used)]
            if after_used or bit == '1':
                options.append((valid_before, True))
            for option in options:
                next_states[option] = next_states.get(option, 0) + count
        states = next_states
    return 2 * sum(count for (valid_before, _), count in states.items() if valid_before)


def word_count(v, y):
    width = v.bit_length()
    untagged = (1 << width) - int(bin(v)[2:][::-1], 2)
    return (25 * ((1 + p(v)) * (1 + p(y)) + (v + 1) * (y + 1)
                  - int(v % 2 == 0 and y % 2 == 0))
            + 28 * (two_slots(v) * (y + 1) - untagged) - 10)


def right_addresses(center, v):
    """Exact inverse on the entire right flank for a prescribed center and reader word."""
    target = center - 1
    slope = 28 * two_slots(v) + 25 * (v + 1)
    deletion_count = 25 * (1 + p(v))
    found = [0] if target == word_count(v, 0) else []
    for parity in (0, 1):
        intercept = word_count(v, parity) - slope * parity
        found.extend(inverse_one_sided(target, intercept, deletion_count, slope, parity))
    return sorted(found)


def check_pair(d, tags, edge, vrows, vtags, y, expected, order):
    deletion = ordinal(d, vrows)
    tag = tags | (vtags << len(d))
    twin = flipped(deletion, *edge)
    reads = [profile(g, tag) for g in (deletion, twin)]
    counts = [count_from_profile(data, 0, y) for data in reads]
    assert counts == [expected - 1, expected + 1]
    lower, upper = [star_rows(g, tag, 0, y) for g in (deletion, twin)]
    assert len(lower) == len(upper) == order and strong(lower) and strong(upper)
    assert flipped(lower, *edge) == upper
    print(order, 'vertices:', *counts, 'at right address', y)
    return reads


def perfect_shape(r):
    return (1 << (r - 1)) * ((1 << r) - 1)


def receiver_census():
    for n, expected in enumerate((1, 2, 5, 13, 34), start=1):
        edges = tuple(combinations(range(n), 2))
        total = 0
        for code in range(1 << len(edges)):
            r = [0] * n
            for bit, (i, j) in enumerate(edges):
                if code >> bit & 1:
                    i, j = j, i
                r[i] |= 1 << j
            starts = [subpaths(r, 1 << i)[-1] for i in range(n)]
            # Exactly one tagged start of multiplicity one; all zero-mass tags are free.
            total += starts.count(1) * (1 << starts.count(0))
        assert total == factorial(n) * expected
    print('All labelled receivers through order 5:', [1, 4, 30, 312, 4080])


def main():
    receiver_census()
    for d, tags, edge, vrows, vtags, y, order in WITNESSES:
        check_pair(d, tags, edge, vrows, vtags, y, perfect_shape(31), order)
    for offset, d, tags, edge, vrows, vtags, y, order in NEARBY:
        deletion = ordinal(d, vrows)
        tags |= vtags << len(d)
        graphs = (deletion, flipped(deletion, *edge))
        center = perfect_shape(31) + offset
        assert [count_from_profile(profile(g, tags), 0, y) for g in graphs] == [center - 1, center + 1]
        pair = [star_rows(g, tags, 0, y) for g in graphs]
        assert flipped(pair[0], *edge) == pair[1]
        a, b = remove_common_extremes(*pair)
        assert len(a) == len(b) == order and strong(a) and strong(b)
        assert sum((u ^ v).bit_count() for u, v in zip(a, b)) == 2
        print('Strong unit edge at E31 offset', offset, 'on', len(a), 'vertices')
    offsets = {0, *(entry[0] for entry in NEARBY)}
    endpoints = {offset + delta for offset in offsets for delta in (-1, 1)}
    assert set(range(-5, 22, 2)) <= endpoints
    for exponent, v, y in ((23, 250, 184567259), (31, 6, 1973896132951836)):
        target = perfect_shape(exponent)
        assert word_count(v, y) == target - 1
        reader, tag = word_reader(v)
        check_pair(CORE, 1, (1, 3), reader, tag, y, target,
                   len(CORE) + len(reader) + y.bit_length() + 1)
    for v, y in ((1, 2), (2, 3), (3, 1), (4, 1)):
        reader, tag = word_reader(v)
        data = profile(reader, tag)
        assert data[0] == 1 + p(v) and data[1] == two_slots(v)
        assert sum(data[2]) == v + 1 and data[2][0] - data[2][1] == -(v % 2 == 0)
        deletion = ordinal(CORE, reader)
        tags = 1 | (tag << len(CORE))
        for increment, g in enumerate((deletion, flipped(deletion, 1, 3))):
            literal = star_rows(g, tags, 0, y)
            assert subpaths(literal)[-1] == word_count(v, y) + 2 * increment


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exponent', type=int, help='choose the center 2^(r-1)(2^r-1)')
    parser.add_argument('--reader', type=int, help='positive binary reader address v')
    args = parser.parse_args()
    if args.exponent is None and args.reader is None:
        main()
    elif args.exponent is None or args.exponent < 2 or args.reader is None or args.reader < 1:
        parser.error('choose both a positive reader and an exponent at least two')
    else:
        center = perfect_shape(args.exponent)
        print('right addresses for center', center, 'and reader', args.reader,
              right_addresses(center, args.reader))
