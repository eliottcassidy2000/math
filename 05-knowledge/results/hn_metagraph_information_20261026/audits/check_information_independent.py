"""Independent path-subcube census, exact Fourier audit, and historical probes."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations
import json


def need(value, message):
    if not value:
        raise RuntimeError(message)


def census(n):
    edges = list(combinations(range(n), 2))
    idx = {e: i for i, e in enumerate(edges)}
    size = 1 << len(edges)
    h = [0] * size
    # Each directed Hamiltonian path contributes 1 on a subcube. This is
    # independent of the lane's per-tournament dynamic program.
    for p in permutations(range(n)):
        one = zero = 0
        for a, b in zip(p, p[1:]):
            bit = 1 << idx[tuple(sorted((a, b)))]
            if a < b:
                one |= bit
            else:
                zero |= bit
        free = (size - 1) ^ (one | zero)
        sub = free
        while True:
            h[one | sub] += 1
            if not sub:
                break
            sub = (sub - 1) & free
    return edges, h


def walsh(values):
    a = values[:]
    stride = 1
    while stride < len(a):
        for start in range(0, len(a), 2 * stride):
            for j in range(start, start + stride):
                x, y = a[j], a[j + stride]
                a[j], a[j + stride] = x + y, x - y
        stride *= 2
    return a


def outdegrees(n, edges, mask):
    d = [0] * n
    for k, (a, b) in enumerate(edges):
        d[a if (mask >> k) & 1 else b] += 1
    return tuple(d)


def main():
    rows = []
    historical = {}
    for n in range(2, 7):
        edges, h = census(n)
        m, size = len(edges), len(h)
        ft = walsh(h)
        need(all(h[t] == h[t ^ (size - 1)] for t in range(size)), 'converse symmetry')
        need(all(v == 0 for i, v in enumerate(ft) if i.bit_count() % 2), 'odd Fourier vanish')
        adjacent = [(i, j) for i, j in combinations(range(m), 2) if set(edges[i]) & set(edges[j])]
        from math import factorial
        expected = Q(factorial(n - 2), 2**(n - 2))
        need(all(abs(Q(ft[(1 << i) | (1 << j)], size)) == expected for i, j in adjacent),
             'adjacent degree-two coefficient magnitude')
        # Build an actual spanning tree of the arc line graph, not all pairs.
        seen, tree = {0}, []
        while len(seen) < m:
            for i, j in adjacent:
                if (i in seen) != (j in seen):
                    tree.append((i, j))
                    seen.update((i, j))
                    break
            else:
                raise RuntimeError('disconnected arc line graph')
        signatures = Counter(tuple(((t >> i) ^ (t >> j)) & 1 for i, j in tree) for t in range(size))
        need(len(signatures) == size // 2 and set(signatures.values()) == {2}, 'sharp binary fibre')
        direct_classes = None
        if n <= 4:
            direct = defaultdict(list)
            for t in range(size):
                direct[tuple(h[t ^ f] for f in range(size))].append(t)
            need(len(direct) == size // 2 and all(v == [v[0], v[0] ^ (size - 1)] for v in direct.values()),
                 'literal full-response fibres are converse pairs')
            direct_classes = len(direct)
        rows.append({'n': n, 'labelled_tournaments': size, 'adjacent_arc_pairs': len(adjacent),
                     'tree_correlations': len(tree), 'response_classes': len(signatures),
                     'literal_response_classes': direct_classes,
                     'quadratic_coefficient_magnitude': str(expected)})
        if n == 5:
            mean = Q(sum(h), size)
            variance = sum((Q(v) - mean)**2 for v in h) / size
            third = sum((Q(v) - mean)**3 for v in h) / size
            need(third == Q(165, 32), 'historical Fisher hostile')
            score_fibres = defaultdict(list)
            for t, v in enumerate(h):
                score_fibres[outdegrees(n, edges, t)].append(v)
            conditional_variance = sum(len(vs) * (Q(sum(vs), len(vs)) - mean)**2
                                       for vs in score_fibres.values()) / size
            degree2 = sum(Q(v, size)**2 for i, v in enumerate(ft) if i.bit_count() == 2)
            need(conditional_variance != degree2, 'conditional score projection differs from Fourier truncation')
            historical = {'n': 5, 'mean_H': str(mean), 'variance_H': str(variance),
                          'Fisher_derivative_at_beta_zero': str(third),
                          'score_conditional_variance': str(conditional_variance),
                          'degree_two_Walsh_energy': str(degree2),
                          'OCR': str(conditional_variance / variance),
                          'degree_two_fraction': str(degree2 / variance)}
    print(json.dumps({'status': 'PASS', 'method': 'Hamiltonian path subcube accumulation + integer FWHT',
                      'response_controls': rows, 'historical_counterexamples': historical}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
