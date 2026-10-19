#!/usr/bin/env python3
"""Exact 8-vertex endpoint obstruction; standalone standard-library proof control."""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from fractions import Fraction
import json
NAMES = ('b', 'x', 'z', 'y', 'a', 't1', 't2', 't3')
EDGES = (('b', 'z'), ('z', 'a'), ('x', 'a'), ('b', 'y'), ('t1', 't2'), ('t2', 't3')) + tuple(((t, q) for t in ('t1', 't2', 't3') for q in ('x', 'y', 'a')))
CORE_ORDERS = ((0, 1, 2, 3, 4), (0, 1, 2, 4, 3), (0, 1, 3, 2, 4), (0, 2, 1, 3, 4), (0, 2, 1, 4, 3), (0, 2, 3, 1, 4), (0, 3, 1, 2, 4), (0, 3, 2, 1, 4), (1, 0, 2, 3, 4), (1, 0, 2, 4, 3), (1, 0, 3, 2, 4))
CORE_PAIRS = ((0, 1), (1, 2), (1, 3), (2, 3), (3, 4))

def main():
    n = len(NAMES)
    p = [0] * n
    for x, y in EDGES:
        p[NAMES.index(y)] |= 1 << NAMES.index(x)
    changed = True
    while changed:
        changed = False
        for j in range(n):
            z = p[j]
            for i in range(n):
                if p[j] >> i & 1:
                    z |= p[i]
            if z != p[j]:
                p[j] = z
                changed = True
    words = []

    def rec(word, mask):
        if len(word) == n:
            words.append(word)
            return
        for i in range(n):
            if not mask >> i & 1 and (not p[i] & ~mask):
                rec(word + (i,), mask | 1 << i)
    rec((), 0)
    E = len(words)
    hist = Counter((tuple((i for i in w if i < 5)) for w in words))
    counts = [hist[o] for o in CORE_ORDERS]
    if not (E == 53 and counts == [4, 4, 4, 10, 10, 10, 4, 4, 1, 1, 1]):
        raise AssertionError('exact certificate check failed')
    pos = [{v: j for j, v in enumerate(w)} for w in words]

    def flag(w, pair):
        return int(w[pair[0]] < w[pair[1]])

    def quadrant(first, second):
        a = [[0, 0], [0, 0]]
        for w in pos:
            a[flag(w, first)][flag(w, second)] += 1
        return a
    rows = []
    for first in CORE_PAIRS:
        true = sum((flag(w, first) for w in pos))
        heavy = int(2 * true >= E)
        mass = true if heavy else E - true
        maxcells = [max(quadrant(first, second)[heavy]) for second in CORE_PAIRS]
        rows.append(dict(first=[NAMES[i] for i in first], heavy_true=bool(heavy), heavy_mass=mass, best_second_largest_cell=min(maxcells), all_second_largest_cells=maxcells, best_heavy_branch_deficit=E - mass - min(maxcells)))
    if not [r['best_heavy_branch_deficit'] for r in rows] == [-29, -1, -2, -12, -5]:
        raise AssertionError('exact certificate check failed')
    inc = [(i, j) for i, j in combinations(range(n), 2) if not (p[j] >> i & 1 or p[i] >> j & 1)]
    degrees = [sum((i in pair for pair in inc)) for i in range(n)]
    if not max(degrees) == 5:
        raise AssertionError('exact certificate check failed')
    best = (-100, None, None, None)
    for first in inc:
        for second in inc:
            tab = quadrant(first, second)
            margin = min((E - sum(tab[i]) - tab[i][j] for i in (0, 1) for j in (0, 1)))
            if margin > best[0]:
                best = (margin, first, second, tab)
    margin, first, second, tab = best
    if not margin >= 0:
        raise AssertionError('exact certificate check failed')
    return dict(status='PASS', vertices=n, maximum_incomparability_degree=max(degrees), extensions=E, core_order_counts=counts, core_first_comparison_certificates=rows, best_core_adaptive_worst_path_deficit=-1, unrestricted_fixed_strategy=dict(first=[NAMES[i] for i in first], second=[NAMES[i] for i in second], quadrants=tab, minimum_deficit=margin))
if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
