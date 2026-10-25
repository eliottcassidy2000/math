"""Independent literal control for the signed bivariate inventory.

Only the answer function is imported. Canonicalization uses simultaneous
permutations of literal arcs and nonnegative integer vertex colors.
"""
import importlib.util
from itertools import combinations, permutations
import json
from pathlib import Path


spec = importlib.util.spec_from_file_location('candidate_excess', Path(__file__).with_name('check_excess.py'))
candidate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(candidate)


def colors(total, n):
    if n == 0:
        if total == 0:
            yield ()
        return
    for first in range(total + 1):
        for tail in colors(total - first, n - 1):
            yield (first,) + tail


def literal_colored_classes(n, excess):
    pairs = list(combinations(range(n), 2))
    indices = {pair: i for i, pair in enumerate(pairs)}
    perms = list(permutations(range(n)))
    representatives = set()
    for mask in range(1 << len(pairs)):
        def arc(i, j):
            if i < j:
                return (mask >> indices[(i, j)]) & 1
            return 1 - arc(j, i)
        renamed = [
            (perm, sum(arc(perm[a], perm[b]) << i
                       for i, (a, b) in enumerate(pairs)))
            for perm in perms]
        for coloring in colors(excess, n):
            key = min((tuple(coloring[v] for v in perm), arcs)
                      for perm, arcs in renamed)
            representatives.add(key)
    return len(representatives)


def main():
    checked = 0
    for actual in range(5):
        for excess in range(6):
            literal = literal_colored_classes(actual, excess)
            proposed = candidate.value(actual, excess)
            if proposed != (-1) ** excess * literal:
                raise RuntimeError(('colored inventory mismatch', actual, excess, proposed, literal))
            checked += 1
    sign_checks = 0
    for actual in range(21):
        for excess in range(21 - actual):
            proposed = candidate.value(actual, excess)
            if proposed.denominator != 1 or (-1) ** excess * proposed < 0:
                raise RuntimeError(('integrality/sign', actual, excess, proposed))
            sign_checks += 1
    print(json.dumps({
        'status': 'PASS',
        'literal_colored_inventories': checked,
        'actual_degree_range': [0, 4],
        'excess_degree_range': [0, 5],
        'sign_integrality_coefficients': sign_checks,
        'sign_integrality_total_degree_range': [0, 20],
        'scope': 'signed colored source tournaments; not vertex counts of resulting cores'
    }, indent=2))


if __name__ == '__main__':
    main()
