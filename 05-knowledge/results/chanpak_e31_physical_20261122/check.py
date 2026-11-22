"""Check the exponent-31 edges and the two independent binary-address formula."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'chanpak_gap_star_ports_20261122'))
from replay import count_from_profile, flipped, p, profile, star_rows, strong, subpaths

CORE = (4, 5, 24, 19, 35, 15)
# D, tags, flip, V, tags, right address, order.
WITNESSES = (
    (CORE, 1, (1, 3), (2, 8, 3, 5), 8, 1973896132951836, 62),
    (CORE, 1, (1, 3), (0, 5, 17, 7, 11), 17, 986948066475918, 62),
    ((96, 13, 25, 49, 67, 22, 46), 4, (0, 4), (0, 17, 11, 19, 37, 15), 16,
     42773862359093, 60),
)


def ordinal(left, right):
    later = ((1 << len(right)) - 1) << len(left)
    return tuple(r | later for r in left) + tuple(r << len(left) for r in right)


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


def main():
    for d, tags, edge, vrows, vtags, y, order in WITNESSES:
        check_pair(d, tags, edge, vrows, vtags, y, perfect_shape(31), order)
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
    main()
