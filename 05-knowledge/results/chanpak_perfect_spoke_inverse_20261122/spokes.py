"""Invert a prescribed marked-spoke flip over two independently addressed flanks."""

import argparse
from math import gcd, isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'chanpak_gap_star_ports_20261122'))
from replay import count_from_profile, profile, star_rows, subpaths


def ceil_div(num, den):
    return -((-num) // den)


def exact(rows, tags, vertex, target, change):
    """All positive left/right addresses for this ordered pair of counts."""
    if tags >> vertex & 1:
        raise ValueError('the contact to be added must initially be untagged')
    low = profile(rows, tags)
    high = profile(rows, tags | (1 << vertex))
    h, q = low[:2]
    left = sum(low[2])
    right = sum(low[3])
    A = sum(high[2]) - left
    B = right - sum(high[3])
    g = gcd(A, B)
    a, b = A // g, B // g
    answers = []
    for i in (0, 1):
        for j in (0, 1):
            constants = [d[4] - (d[3][0] - d[3][1]) * i - (d[2][0] - d[2][1]) * j
                         for d in (low, high)]
            rhs = change - (constants[1] - constants[0]) - A * j + B * i
            if rhs % (2 * g):
                continue
            c = rhs // (2 * g)
            first = (-c * pow(b, -1, a)) % a if a > 1 else 0
            second = (c + b * first) // a
            shift = max(0, ceil_div(1 - (2 * first + i), 2 * a),
                        ceil_div(1 - (2 * second + j), 2 * b))
            x0, y0 = 2 * (first + a * shift) + i, 2 * (second + b * shift) + j
            dx, dy = 2 * a, 2 * b
            C = constants[0] + q * x0 * y0 + right * x0 + left * y0
            if target < C:
                continue
            linear = q * (x0 * dy + y0 * dx) + right * dx + left * dy
            quadratic = q * dx * dy
            t = (isqrt(linear**2 + 4 * quadratic * (target - C)) - linear) // (2 * quadratic)
            x, y = x0 + dx * t, y0 + dy * t
            if count_from_profile(low, x, y) == target:
                assert count_from_profile(high, x, y) == target + change
                answers.append((x, y))
    return answers


def ordinal(left, right):
    later = ((1 << len(right)) - 1) << len(left)
    return tuple(row | later for row in left) + tuple(row << len(left) for row in right)


def dyadic_markers(rows):
    n = len(rows)
    full = (1 << n) - 1
    paths = subpaths(rows)
    paths[0] = 1
    q = sum(paths[s] * paths[full ^ s] for s in range(full + 1))
    if q & (q - 1):
        return [], paths[-1]
    reverse = tuple(row ^ full ^ (1 << i) for i, row in enumerate(rows))
    output = []
    for u in range(n):
        masses = []
        for graph in (rows, reverse):
            rooted = subpaths(graph, 1 << u)
            masses.append(sum(rooted[s] * paths[full ^ s] for s in range(1, full + 1)))
        A, B = masses
        common = gcd(A, B)
        a, b = A // common, B // common
        if a & (a - 1) or b & (b - 1):
            continue
        if isqrt(2 * q * a * b) ** 2 == 2 * q * a * b:
            output.append((u, A, B, q))
    return output, paths[-1]


def catalogue():
    catalogue_path = HERE.parents[2] / 'session212' / 'cycle-levels'
    for n in range(1, 9):
        classes = marked = nontransitive = 0
        for line in (catalogue_path / f'classes{n}.rows').read_text().splitlines():
            rows = tuple(map(int, line.split()))
            markers, h = dyadic_markers(rows)
            classes += 1
            marked += len(markers)
            if h > 1:
                nontransitive += len(markers)
        print(f'n={n} classes={classes} dyadically aligned markers={marked} nontransitive={nontransitive}')
        assert marked == n and nontransitive == 0


def checks():
    for rows, tags, u in [((0,), 0, 0), ((2, 4, 1), 1, 1), ((2, 8, 3, 5), 9, 1)]:
        for x in range(1, 9):
            for y in range(1, 9):
                ss = [tags, tags | 1 << u]
                values = [count_from_profile(profile(rows, mask), x, y) for mask in ss]
                assert (x, y) in exact(rows, tags, u, values[0], values[1] - values[0])
                if max(x, y) <= 3:
                    assert [subpaths(star_rows(rows, mask, x, y))[-1] for mask in ss] == values

    core = (4, 65, 50, 39, 75, 147, 45, 95)
    E31 = (1 << 30) * ((1 << 31) - 1)
    assert (25441035, 21176353) in exact(core, 244, 3, E31 - 78364 + 1, -2)
    sandwich = ordinal(ordinal((2, 8, 3, 5), (0,)), (2, 8, 3, 5))
    markers, h = dyadic_markers(sandwich)
    assert h == 25 and (4, 160, 160, 2048) in markers
    signed_terms = []
    for side, scalar in (((2, 4, 1), 5), ((0, 1, 11, 35, 15, 23), 3)):
        weights = [scalar * profile(side, 1 << i)[2][1] for i in range(len(side))]
        sums = {0}
        for weight in weights:
            sums |= {old + weight for old in sums}
        signed_terms.append(sums)
    differences = {a-b for a in signed_terms[0] for b in signed_terms[1]}
    assert gcd(*differences) == 1 and min(abs(n) for n in differences if n) == 5
    print('Checked 192 inverse queries, literal small hosts, and the nine-vertex aligned separator.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', action='store_true')
    args = parser.parse_args()
    checks()
    if args.catalogue:
        catalogue()
