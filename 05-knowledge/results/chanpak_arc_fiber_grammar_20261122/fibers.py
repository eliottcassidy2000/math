"""Ordinary and single-arc size series for a prescribed positive odd height."""

import argparse
import gzip
from collections import Counter
from functools import cache
from itertools import combinations
from math import comb, factorial, gcd
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "chanpak_birth_atlas_metagraph_20261117" / "factor"))
from solve import paths, strong, strong_cap
sys.path.insert(0, str(HERE.parent / "chanpak_minimal_core_grammar_20261122"))
from grammar import canon, strong_classes


def automorphisms(rows):
    """Count degree-preserving vertex permutations preserving every oriented pair."""
    degrees = [r.bit_count() for r in rows]
    candidates = [[j for j, d in enumerate(degrees) if d == degrees[i]] for i in range(len(rows))]
    order = sorted(range(len(rows)), key=lambda i: len(candidates[i]))

    def extend(images, used):
        if len(images) == len(rows):
            return 1
        i = order[len(images)]
        total = 0
        for j in candidates[i]:
            if used >> j & 1:
                continue
            if all((rows[i] >> k & 1) == (rows[j] >> l & 1) for k, l in images):
                total += extend(images + [(i, j)], used | 1 << j)
        return total

    return extend([], 0)


def catalog(n):
    source = ROOT / "session212" / "cycle-levels" / f"classes{n}.rows"
    for line in source.read_text().splitlines():
        yield tuple(map(int, line.split()))


def bounded_classes(stop, cap):
    if stop <= 8:
        for n in range(3, stop + 1):
            for rows in catalog(n):
                if strong(rows):
                    h = paths(rows)
                    if h <= cap:
                        yield n, rows, h
    else:
        for n, code, h in strong_classes(cap, max_order=stop):
            yield n, tuple(canon.rows(n, code)), h


def full_catalog(stop, cap, target=None):
    """Count strong scalar classes and central marked upward unit pairs."""
    atoms, physical = Counter(), Counter()
    for n, rows, h in bounded_classes(stop, cap):
        if target is not None and target % h:
            continue
        atoms[h, n] += 1
        if h <= 3 or (target is not None and h != target):
            continue
        good = 0
        for i, j in combinations(range(n), 2):
            other = list(rows)
            other[i] ^= 1 << j
            other[j] ^= 1 << i
            if strong(other) and paths(other) == h + 2:
                good += 1
        if good:
            group = automorphisms(rows)
            assert good % group == 0
            physical[h, n] += good // group
    return atoms, physical


def scalar_terms(h, atoms):
    """Map (total strong order, number of components) to ordered word counts."""
    @cache
    def words(t):
        if t == 1:
            return Counter({(0, 0): 1})
        output = Counter()
        for (d, n), multiplicity in atoms.items():
            if t % d:
                continue
            for (s, k), rest in words(t // d).items():
                output[s + n, k + 1] += multiplicity * rest
        return output

    return words(h)


def edge_terms(h, physical):
    if h in (1, 3):
        return {3: 1} if h == 1 else {4: 2}
    return {n: v for (d, n), v in physical.items() if d == h and v}


def coefficients(h, n, words, edges):
    scalar = sum(w * comb(n-s+k, k) for (s, k), w in words.items() if s <= n)
    marked = sum(w * (n-s+1) for s, w in edges.items() if s <= n)
    return scalar, marked


def polynomial(items):
    return ' + '.join(f'{value} z^{degree}' for degree, value in sorted(items.items()) if value) or '0'


def census_rows(n):
    bank = HERE.parent / "path_spectrum_two_angles_metagraph_20261106" / "data"
    path = bank / f"census_{n}.txt"
    with (path.open() if n < 11 else gzip.open(str(path) + ".gz", "rt")) as stream:
        for line in stream:
            if line.startswith("T "):
                yield tuple(map(int, line.split()[1:]))


def ordinary_query(h, order):
    if order is not None and 3 <= order <= 11:
        return sum(count for height, _, _, count in census_rows(order) if height == h)
    stop = strong_cap(h) if order is None else min(strong_cap(h), order)
    atoms = Counter()
    if stop <= 11:
        for n in range(3, stop + 1):
            for height, _, is_strong, count in census_rows(n):
                if is_strong and h % height == 0:
                    atoms[height, n] += count
    else:
        for n, _, height in bounded_classes(stop, h):
            if h % height == 0:
                atoms[height, n] += 1
    words = scalar_terms(h, atoms)
    if order is None:
        return ' + '.join(f'{w} z^{s}/(1-z)^{k+1}' for (s, k), w in sorted(words.items())) or '0'
    return coefficients(h, order, words, {})[0]


def fermat_births():
    for n, h in ((3, 3), (4, 5), (6, 17), (8, 257), (11, 65537)):
        levels = {t: count for height, t, strong_flag, count in census_rows(n)
                  if height == h and strong_flag}
        print(h, 'first order', n, 'classes', sum(levels.values()), 'by cyclic triangles', levels)


def print_query(h, order):
    stop = max(0, strong_cap(h))
    if order is not None:
        stop = min(stop, order)
    atoms, physical = full_catalog(stop, h, target=h)
    words, edges = scalar_terms(h, atoms), edge_terms(h, physical)
    if order is not None:
        scalar, marked = coefficients(h, order, words, edges)
        print(f'order {order}, H={h}: {scalar} classes; H={h} to {h+2}: {marked} marked flips')
        return
    print('F:', ' + '.join(f'{w} z^{s}/(1-z)^{k+1}' for (s, k), w in sorted(words.items())) or '0')
    print('U:', f'({polynomial(edges)})/(1-z)^2')
    if words:
        first = min(s for s, k in words)
        print('first order and number of ordinary classes:', first,
              sum(w for (s, k), w in words.items() if s == first))


def check_padding():
    atoms, physical = full_catalog(6, 45)
    for h in (1, 3, 5, 7, 9, 11, 13, 15):
        words, edges = scalar_terms(h, atoms), edge_terms(h, physical)
        for n in (3, 4, 5, 6):
            ordinary = marked = 0
            for rows in catalog(n):
                if paths(rows) != h:
                    continue
                ordinary += 1
                matches = 0
                for i, j in combinations(range(n), 2):
                    other = list(rows)
                    other[i] ^= 1 << j
                    other[j] ^= 1 << i
                    matches += paths(other) == h + 2
                group = automorphisms(rows)
                assert matches % group == 0
                marked += matches // group
            assert coefficients(h, n, words, edges) == (ordinary, marked)
    print('Ordinary and marked counts agree with direct class/orbit counts through order six.')


def e61_bounds():
    e = 2**60 * (2**61 - 1)
    cap = strong_cap(e - 1)
    # Split a permutation on 38 vertices into ordered blocks of sizes 8,10,10,10.
    bound = factorial(38) * 661 * 15745**3 // (factorial(8) * factorial(10)**3)
    assert bound < e - 1 and cap == strong_cap(e + 1) == 157
    assert gcd(e + 1, factorial(cap)) == 1
    print('E61 unit-pair core order: 39..157. Every upper class is rigid.')
    print('Neighbors:', e - 1, e + 1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fermat-births', action='store_true', help='read the five exact minimum fibers')
    parser.add_argument('--ordinary', action='store_true', help='only ordinary isomorphism classes, using the census through order 11')
    parser.add_argument('--height', type=int, help='positive odd lower count')
    parser.add_argument('--order', type=int, help='compute only the specified coefficient')
    parser.add_argument('--small-table', action='store_true', help='all odd lower heights through 73')
    parser.add_argument('--e61', action='store_true', help='evaluate the rigorous order and rigidity tests')
    args = parser.parse_args()
    if args.fermat_births:
        fermat_births()
    elif args.small_table:
        atoms, marked = full_catalog(8, 73)
        for h in range(1, 74, 2):
            words = scalar_terms(h, atoms)
            first = min((s for s, k in words), default=None)
            number = sum(w for (s, k), w in words.items() if s == first)
            print(h, 'minimum', first, 'classes', number, 'edge numerator', polynomial(edge_terms(h, marked)))
    elif args.e61:
        e61_bounds()
    elif args.height is None and args.order is None:
        check_padding()
    elif args.height is None or args.height < 1 or args.height % 2 == 0 or (args.order is not None and args.order < 0):
        parser.error('choose a positive odd --height and optionally a nonnegative --order')
    elif args.ordinary:
        print(ordinary_query(args.height, args.order))
    else:
        print_query(args.height, args.order)
