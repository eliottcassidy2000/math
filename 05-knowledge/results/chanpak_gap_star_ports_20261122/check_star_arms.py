"""Check the local determinant data and literal ordinal/star compositions."""

from replay import flipped, p, profile, star_rows, subpaths


def ordinal(first, second):
    offset = len(first)
    suffix = ((1 << len(second)) - 1) << offset
    return tuple(row | suffix for row in first) + tuple(row << offset for row in second)


def right_reads(rows, tags):
    h, _, (even, odd), _, k = profile(rows, tags)
    assert (h - k) % 2 == 0
    return h, even + (h - k) // 2, even + odd


def star(left, left_tags, right, right_tags, width):
    transitive = tuple((1 << width) - (1 << (i + 1)) for i in range(width))
    mask = (1 << len(left)) - 1
    converse = tuple(mask ^ row ^ (1 << i) for i, row in enumerate(left))
    deletion = ordinal(ordinal(right, transitive), converse)
    tags = right_tags | ((mask ^ left_tags) << (len(right) + width))
    return star_rows(deletion, tags, 0, 0)


def assert_edge(first, second, counts):
    assert len(first) == len(second)
    changed = [i for i, (a, b) in enumerate(zip(first, second)) if a != b]
    i, j = changed
    assert flipped(first, i, j) == second
    assert tuple(subpaths(row)[-1] for row in (first, second)) == counts


def main():
    rays = [
        ((34, 136, 19, 133, 43, 78, 31, 117), 164, (5, 7), 2230,
         (431, 215), (435, 218)),
        ((8, 5, 9, 18, 7), 3, (0, 1), 32, (11, 4), (13, 4)),
    ]
    for rows, tags, endpoints, slope, first, second in rays:
        twin = flipped(rows, *endpoints)
        assert right_reads(rows, tags) == (*first, slope)
        assert right_reads(twin, tags) == (*second, slope)
        for m in (1, 2):
            expected = tuple(g - 2 * j + (slope << m) for g, j in (first, second))
            assert_edge(star((), 0, rows, tags, m), star((), 0, twin, tags, m), expected)

    left = (20, 29, 24, 17, 0)
    left_twin = flipped(left, 0, 1)
    for rows, reads in ((left, (3, 29, 30)), (left_twin, (5, 31, 28))):
        h, _, masses, _, k = profile(rows, 10)
        assert (h, k, sum(masses)) == reads
    right = (24, 1, 3, 6, 14)
    for tags, mass in ((6, 29), (12, 37), (24, 47)):
        assert right_reads(right, tags) == (11, 5, mass)
        m = 1
        center = 290 + (33 * mass << m)
        assert_edge(star(left, 10, right, tags, m), star(left_twin, 10, right, tags, m),
                    (center - 1, center + 1))

    block = (6, 68, 136, 131, 15, 31, 61, 115)
    twin = flipped(block, 0, 3)
    for source, h in ((block, 313), (twin, 307)):
        assert profile(source, 0)[:2] == (h, 5008)
    # Start with an intercept-one star with a nonempty left arm.
    left, right = (0,), (0,)
    all_suffix_tags = ((1 << len(block)) - 1) << len(left)
    m = 2
    assert_edge(star(ordinal(left, block), all_suffix_tags, right, 1, m),
                star(ordinal(left, twin), all_suffix_tags, right, 1, m),
                (313 + (5008 << m), 307 + (5008 << m)))
    def joined(b, e):
        return b * e + (1 + p(b - 1)) * (1 + p(e - 1)) - (b & e & 1)

    left_factor = joined(2, 159198469)
    right_factor = joined(2, 1996165791)
    perfect = (1 << 30) * ((1 << 31) - 1)
    assert (left_factor, right_factor) == (457256543, 5042777503)
    assert left_factor * right_factor == perfect + 1
    assert joined(17, 113331563) == 7 * left_factor
    print('Checked two prescribed rays, three marked unit-gap transfers, and the sharp six-count edge.')


if __name__ == '__main__':
    main()
