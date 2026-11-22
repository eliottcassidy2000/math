"""Small deletion cores, tags and addresses for the exponent-31 unit edges."""

CORE = (4, 5, 24, 19, 35, 15)
# D, tags, flip, V, tags, right address, order.
WITNESSES = (
    (CORE, 1, (1, 3), (2, 8, 3, 5), 8, 1973896132951836, 62),
    (CORE, 1, (1, 3), (0, 5, 17, 7, 11), 17, 986948066475918, 62),
    ((96, 13, 25, 49, 67, 22, 46), 4, (0, 4), (0, 17, 11, 19, 37, 15), 16,
     42773862359093, 60),
)

# Offset from E31, D, tags, flip, V, tags, y; forced outer endpoints may be removed.
NEARBY = (
    (-4, (24, 37, 137, 50, 70, 213, 15, 91), 227, (5, 6), (), 0, 524479531053167, 58),
    (12, (8, 33, 35, 22, 135, 89, 31, 111), 3, (0, 1), (), 0, 3403860470244866, 61),
    (20, (4, 9, 10, 97, 15, 87, 151, 63), 4, (3, 4), (), 0, 11698292815714334, 63),
    (-20, (0, 17, 139, 67, 45, 79, 151, 59), 17, (2, 3), (), 0, 6236156999269463, 62),
    (-2, (0, 5, 49, 39, 75, 147, 47, 95), 12, (4, 5), (), 0, 4848041412776110, 62),
    (4, (0, 5, 129, 7, 15, 31, 63, 123), 3, (2, 4), (), 0, 20706137851495925, 64),
    (16, (0, 17, 139, 67, 45, 79, 151, 59), 17, (2, 3), (), 0, 5546712681070622, 62),
    (-16, (2, 64, 3, 7, 15, 31, 61), 1, (1, 3), (8, 5, 9, 18, 7), 1, 566408789592188, 63),
    (-12, (2, 4, 17, 7, 11), 2, (2, 3), (2, 4, 17, 7, 11), 1, 2544946116269434, 63),
    (8, CORE, 1, (1, 3), (2, 4, 1), 1, 5242259031785913, 63),
    (10, (2, 16, 11, 35, 13, 215, 31, 95), 2, (5, 6), (2, 8, 3, 5), 1, 791860274676473, 63),
)


def ordinal(left, right):
    later = ((1 << len(right)) - 1) << len(left)
    return tuple(r | later for r in left) + tuple(r << len(left) for r in right)


def stored_pairs():
    """Yield the lower endpoint and the input for constructing both tournaments."""
    midpoint = (1 << 30) * ((1 << 31) - 1)
    centers = ((0, *record) for record in sorted(WITNESSES, key=lambda rec: rec[-1]))
    for offset, core, tags, edge, receiver, receiver_tags, y, _ in (*centers, *NEARBY):
        yield (midpoint + offset - 1, ordinal(core, receiver), edge,
               tags | (receiver_tags << len(core)), 0, y)
