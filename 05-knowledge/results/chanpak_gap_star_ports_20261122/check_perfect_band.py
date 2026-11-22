"""Replay each perfect-center certificate, including its strong reduced endpoints."""


# (perfect exponent, center offset, lower outgoing masks, flip, apex tags, left, right)
CERTIFICATES = (
    (17, -8, (4, 5, 40, 67, 15, 27, 55, 127), (5, 6), 5, 0, 60455853),
    (17, 0, (4, 5, 48, 7, 11, 27, 63, 127), (2, 3), 1, 0, 233563022),
    (17, 2, (64, 129, 11, 35, 15, 23, 62, 125), (5, 7), 65, 0, 8485613),
    (17, 4, (64, 37, 177, 23, 99, 137, 46, 91), (2, 5), 194, 17, 79709),
    (17, 6, (2, 4, 17, 7, 11, 31, 63, 127), (2, 3), 2, 0, 542468540),
    (17, 8, (6, 20, 56, 19, 33, 75, 31, 127), (2, 4), 42, 0, 18697743),
    (17, 10, (2, 4, 33, 7, 15, 27, 63, 127), (2, 3), 22, 0, 87183271),
    (17, 12, (2, 8, 27, 49, 35, 135, 63, 95), (2, 4), 9, 0, 38191461),
    (17, 14, (2, 40, 67, 69, 15, 29, 179, 63), (3, 4), 3, 0, 20695423),
    (17, 16, (0, 1, 19, 71, 171, 15, 55, 111), (4, 5), 90, 1316, 5097),
    (19, -14, (8, 69, 33, 22, 39, 75, 29, 127), (1, 2), 8, 0, 697898218),
    (19, -12, (2, 4, 17, 7, 11, 31, 63, 127), (2, 3), 2, 0, 9333359802),
    (19, -6, (2, 64, 3, 7, 15, 31, 61, 127), (1, 3), 1, 0, 2152318488),
    (19, -4, (32, 5, 9, 19, 39, 78, 31, 127), (1, 3), 1, 0, 966051136),
    (19, 0, (34, 80, 75, 19, 69, 30, 41, 127), (2, 6), 1, 0, 553528478),
    (19, 2, (34, 80, 75, 19, 69, 30, 41, 127), (2, 6), 1, 0, 474437876),
    (19, 6, (2, 4, 17, 7, 11, 31, 63, 127), (0, 2), 19, 0, 2902022354),
    (19, 8, (4, 5, 24, 19, 35, 15, 63, 127), (1, 3), 1, 0, 2773805764),
    (19, 12, (6, 32, 10, 35, 15, 85, 31, 127), (0, 1), 11, 0, 548890443),
    (19, 16, (4, 17, 10, 19, 37, 207, 31, 95), (5, 6), 4, 0, 749218164),
)


def remove_common_extremes(left, right):
    """A common source or sink is forced onto the end of every spanning path."""
    active = (1 << len(left)) - 1
    while True:
        for v in range(len(left)):
            if not active >> v & 1:
                continue
            before = left[v] & active
            if before == right[v] & active and before in (0, active ^ (1 << v)):
                active ^= 1 << v
                break
        else:
            break
    keep = [v for v in range(len(left)) if active >> v & 1]

    def restrict(rows):
        return tuple(sum(1 << i for i, w in enumerate(keep) if rows[v] >> w & 1)
                     for v in keep)

    return restrict(left), restrict(right)


def main():
    from replay import count_from_profile, flipped, profile, star_rows, strong

    offsets = {17: [], 19: []}
    for exponent, offset, lower, edge, tags, x, y in CERTIFICATES:
        center = 2**(exponent-1)*(2**exponent-1) + offset
        cores = (lower, flipped(lower, *edge))
        counts = [count_from_profile(profile(core, tags), x, y) for core in cores]
        assert counts == [center-1, center+1]
        graphs = [star_rows(core, tags, x, y) for core in cores]
        a, b = remove_common_extremes(*graphs)
        assert sum((u^v).bit_count() for u,v in zip(a,b)) == 2
        assert strong(a) and strong(b)
        offsets[exponent].append(offset)
        print(exponent, offset, len(a), *counts)
    assert all(x in offsets[17] for x in range(0, 17, 2))


if __name__ == "__main__":
    main()
