"""Recover the four-core tripling equations, the two exact exits and all core edits."""
from itertools import combinations

from direct_fork import address_height, direct, flipped, read, star_rows, tagged_data


CENTER = (6, 4, 8, 3)
TAGS = 5
SOURCE = flipped(CENTER, 2, 3)
UPPER = flipped(CENTER, 0, 3)
READINGS = (
    (1, 16, 5, 5, 1, 1, 3),
    (5, 32, 14, 14, 2, 2, 15),
    (3, 24, 15, 12, 1, 0, 13),
)


def run():
    assert tuple(tagged_data(t, TAGS) for t in (SOURCE, CENTER, UPPER)) == READINGS
    # Twice the change for each edit in the four-core and its four apex tags
    # equals slope*y+offset on either ray. The desired lower edit gives -4.
    for i, j in combinations(range(5), 2):
        chart = (tagged_data(flipped(CENTER, i, j), TAGS) if j < 4
                 else tagged_data(CENTER, TAGS ^ (1 << i)))
        h, _, c, _, e, _, k = chart
        slope, offset = h + 2*c - 33, 2*e + 2*k - 5*h - 9
        assert ((slope == 0 and offset != -4) or
                (slope > 0 and slope*13 + offset > -4) or
                (slope < 0 and slope*13 + offset < -4))
    for y in range(2048):
        source, center, upper = [address_height(d, 0, y) for d in READINGS]
        assert center - 3*source == 2*read(y) - y - y%2 + 6
        assert upper - center == -2*read(y) + y - y%2 - 2
        if center == 3*source and upper == center + 2:
            w, r = y - 1, 0
            while w % 2 == 0:
                w //= 2
                r += 1
            assert r >= 2 and w in (3, 5)
    for r in range(2, 20):
        for w in (3, 5):
            y = (w << r) + 1
            source, center, upper = [address_height(d, 0, y) for d in READINGS]
            assert source == 11*w*2**(r-1) + 7
            assert (center, upper) == (3*source, 3*source+2)
            lower = [k for k in range(y.bit_length())
                     if address_height(READINGS[1], 0, y ^ (1 << k)) == center-2]
            assert lower == ([1] if y in (21, 25) else [])
    for y in (21, 25):
        c = star_rows(CENTER, TAGS, 0, y)
        assert flipped(c, 2, 3) == star_rows(SOURCE, TAGS, 0, y)
        assert flipped(c, 0, 3) == star_rows(UPPER, TAGS, 0, y)
        assert flipped(c, 5, len(c)-1) == star_rows(CENTER, TAGS, 0, y ^ 2)
        h = address_height(READINGS[0], 0, y)
        heights = [direct(t) for t in (flipped(c, 2, 3), c,
                                        flipped(c, 5, len(c)-1), flipped(c, 0, 3))]
        assert heights == [h, 3*h, 3*h-2, 3*h+2]
        print("address", y, "heights", *heights)
    print("two four-core rays; core and flank-tag tests passed")


if __name__ == "__main__":
    run()
