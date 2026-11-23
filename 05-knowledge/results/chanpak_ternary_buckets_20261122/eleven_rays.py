"""Recover the four-core rays and classify all individual lower exits."""
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



def arc_changes(rows):
    """Find every one-arc change by cutting paths into two ordered pieces."""
    n = len(rows)
    full = (1 << n) - 1
    def first_paths(t):
        dp = [[0]*n for _ in range(full+1)]
        for i in range(n):
            dp[1 << i][i] = 1
        for mask in range(1, full+1):
            bits = mask
            while bits:
                bit = bits & -bits
                i = bit.bit_length()-1
                prev = dp[mask ^ bit]
                choices = t[i] & mask
                while choices:
                    first = choices & -choices
                    dp[mask][i] += prev[first.bit_length()-1]
                    choices ^= first
                bits ^= bit
        return [[(j, c) for j, c in enumerate(row) if c] for row in dp]
    start = first_paths(rows)
    end = first_paths([full ^ r ^ (1 << i) for i, r in enumerate(rows)])
    cuts = [[0]*n for _ in rows]
    for left in range(1, full):
        for i, count in end[left]:
            row = cuts[i]
            for j, other in start[full ^ left]:
                row[j] += count*other
    return {(i,j): (cuts[j][i]-cuts[i][j]) * (1 if (rows[i] >> j) & 1 else -1)
            for i,j in combinations(range(n), 2)}


def twice_coefficients(values):
    """Twice the coefficients in 1, 2^m-1, 3^m-2*2^m+1 from m=0,1,2."""
    a, b, c = values
    return [2*a, 2*(b-a), c-3*b+2*a]


def strict_sign(grid):
    coeff = [v for row in grid for v in row]
    return coeff[0] != 0 and all(v*coeff[0] >= 0 for v in coeff)


def all_ordinary_edits():
    for w, expected in ((3, 49), (5, 60)):
        changes = {}
        for s in range(1, 8):
            changes[s] = arc_changes(star_rows(CENTER, TAGS, 0, (w << (s+1))+1))
        outside = range(5 + w.bit_length())  # four core, low digit, high digits
        def index(i, s): return i if i < 5 else i+s
        def at_mark(i, a, b):
            s = a+1+b
            edge = tuple(sorted((index(i,s), 5+a)))
            return changes[s][edge]+2
        total = 0
        for i, j in combinations(outside, 2):
            values = [changes[s][index(i,s),index(j,s)]+2 for s in (1,2,3)]
            assert strict_sign([twice_coefficients(values)])
            total += 1
        for i in outside:
            for a_range in ((0,), (1,2,3)):
                for b_range in ((0,), (1,2,3)):
                    grid = [[at_mark(i,a,b) for b in b_range] for a in a_range]
                    if len(b_range) == 3:
                        grid = [twice_coefficients(row) for row in grid]
                    if len(a_range) == 3:
                        grid = list(zip(*[twice_coefficients(col) for col in zip(*grid)]))
                    assert strict_sign(grid)
                    total += 1
        assert total == expected
        print("odd prefix", w, "strict-sign charts", total)


if __name__ == "__main__":
    run()
    all_ordinary_edits()
