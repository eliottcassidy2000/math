"""Small exact controls for the actual-degree/excess-degree gauge repair."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial, gcd, prod
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def two(e):
    return Fraction(1 << e) if e >= 0 else Fraction(1, 1 << -e)


def beta(k, excess):
    if k == 0:
        return int(excess == 0)
    return (-1) ** excess * comb(k + excess - 1, excess)


def partitions(n, low=1):
    if n == 0:
        yield ()
    for a in range(low, n + 1, 2):
        for tail in partitions(n - a, a):
            yield (a,) + tail


def allocation_coefficient(p, excess):
    """Direct product of the finite needed pieces of (1+z^d)^(-k)."""
    coeff = [1] + [0] * excess
    for d, k in Counter(p).items():
        nxt = [0] * (excess + 1)
        for old, c in enumerate(coeff):
            for e in range((excess - old) // d + 1):
                nxt[old + d * e] += c * beta(k, e)
        coeff = nxt
    return coeff[excess]


def direct(actual, excess, ell, defect):
    total = Fraction(0)
    for p in partitions(actual):
        if p and max(p) > ell:
            continue
        coefficient = allocation_coefficient(p, excess)
        exponent = sum((d - 1) // 2 for d in p)
        exponent += sum(gcd(a, b) for a, b in combinations(p, 2))
        exponent -= sum(defect[(d - 1) // 2] for d in p)
        denominator = prod(d ** k * factorial(k) for d, k in Counter(p).items())
        total += coefficient * two(exponent) / denominator
    return total


@lru_cache(None)
def gauge(ell, actual, excess, defect):
    if actual == 0:
        return Fraction(excess == 0)
    if ell < 1:
        return Fraction(0)
    if ell > actual:
        return gauge(ell - 2, actual, excess,
                     defect[:(min(ell - 2, actual) + 1) // 2])
    if ell == 1:
        return beta(actual, excess) * two(actual * (actual - 1) // 2) / factorial(actual)
    total = Fraction(0)
    for k in range(actual // ell + 1):
        rest = actual - ell * k
        selected = defect[(ell - 1) // 2]
        exponent = k * (actual - selected) - (ell * k * k + k) // 2
        nxt = tuple(defect[(d - 1) // 2] + k * (d - gcd(ell, d))
                    for d in range(1, min(ell - 2, rest) + 1, 2))
        for e in (range(excess // ell + 1) if k else (0,)):
            total += (beta(k, e) * two(exponent) / (ell ** k * factorial(k))
                      * gauge(ell - 2, rest, excess - ell * e, nxt))
    return total


def value(actual, excess):
    if any(type(x) is not int or x < 0 for x in (actual, excess)):
        raise ValueError('actual and excess degrees must be nonnegative integers')
    ell = actual if actual % 2 else actual - 1
    return gauge(ell, actual, excess, (0,) * ((ell + 1) // 2))


def literal_pairfree_classes(n):
    pairs = list(combinations(range(n), 2))
    indices = {edge: i for i, edge in enumerate(pairs)}
    representatives = set()
    for mask in range(1 << len(pairs)):
        def edge(a, b):
            return ((mask >> indices[(a, b)]) & 1) if a < b else 1 - edge(b, a)
        if any(all(edge(a, c) == edge(b, c) for c in range(n) if c not in (a, b))
               for a, b in pairs):
            continue
        representatives.add(min(sum(edge(p[a], p[b]) << i for i, (a, b) in enumerate(pairs))
                                for p in permutations(range(n))))
    return len(representatives)


def main():
    bivariate_checks = 0
    core_counts = []
    for n in range(15):
        total = Fraction(0)
        for actual in range(n + 1):
            excess = n - actual
            ell = actual if actual % 2 else actual - 1
            defect = (0,) * ((ell + 1) // 2)
            result = value(actual, excess)
            need(result == direct(actual, excess, ell, defect), ('direct substitution', actual, excess))
            need(result.denominator == 1, ('colored-inventory integrality', actual, excess))
            need((-1) ** excess * result >= 0, ('colored-inventory sign', actual, excess))
            total += result
            bivariate_checks += 1
        need(total.denominator == 1, ('core integrality', n))
        core_counts.append(total.numerator)
        if n <= 5:
            need(total == literal_pairfree_classes(n), ('literal pair-free classes', n))

    minimality_checks = 0
    for actual in range(1, 9):
        for excess in range(9):
            for ell in (1, 3, 5):
                ds = list(range(1, min(ell, actual) + 1, 2))
                live = [i for i, d in enumerate(ds) if d < actual or excess % actual == 0]
                supported = [p for p in partitions(actual)
                             if max(p) <= ell and allocation_coefficient(p, excess)]
                profiles = [(0,) + q for q in product((0, 1), repeat=len(ds) - 1)]
                for a, b in combinations(profiles, 2):
                    observed_a = tuple(sum(a[(d - 1) // 2] for d in p) for p in supported)
                    observed_b = tuple(sum(b[(d - 1) // 2] for d in p) for p in supported)
                    need((observed_a == observed_b) == all(a[i] == b[i] for i in live),
                         ('exact live-coordinate support', actual, excess, ell, a, b))
                    minimality_checks += 1

    # Same D, different cycle counts: no common scalar merges both core futures.
    need(Fraction(2 ** 1, 2 ** 3) != Fraction(2 ** 2, 2 ** 6), 'ungraded scalar-gauge obstruction')
    need(direct(3, 2, 3, (0, 0)) == direct(3, 2, 3, (0, 8)), 'invisible top coordinate')
    need(direct(3, 3, 3, (0, 0)) != direct(3, 3, 3, (0, 8)), 'visible top coordinate restored')
    for bad in (-1, True, 1.5, '3'):
        for args in ((bad, 0), (0, bad)):
            try:
                value(*args)
            except ValueError:
                pass
            else:
                raise RuntimeError('bad degree accepted')
    print(json.dumps({'status': 'PASS', 'bivariate_coefficients_checked': bivariate_checks,
                      'total_degree_range': [0, 14], 'literal_pairfree_class_range': [0, 5],
                      'individual_continuation_minimality_checks': minimality_checks,
                      'pairfree_counts': core_counts,
                      'scope': 'exact graded repair; no timing or graph-size improvement claimed'}, indent=2))


if __name__ == '__main__':
    main()
