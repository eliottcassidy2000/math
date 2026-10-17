#!/usr/bin/env python3
"""Are the best Paley rankings found by local search close to affine (arithmetic) orders? (S40, HYP-9287 test 3)

For each witness order stored by paley_ranking_local_search_20261017.py, compute the normalised Kendall-tau distance
to the nearest affine order i -> a i + b (a in (Z/p)^*, b in Z/p; both the order and its reverse), as a fraction of the
p(p-1)/2 pairs (a uniformly random order has expected distance 1/2 to any fixed order, and its minimum over the
p(p-1) affine orders is somewhat smaller).
Usage: python paley_best_orders_vs_affine_20261017.py
Session opus-2026-10-17-S40.
"""
import os, json
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', '..', '05-knowledge', 'results')


@njit(cache=True)
def kendall(pos_a, order_b):
    """Number of discordant pairs between the order with position map pos_a and the order listed by order_b
    (inversions of the sequence pos_a[order_b[i]], counted with a Fenwick tree in O(n log n))."""
    n = order_b.shape[0]
    tree = np.zeros(n + 1, dtype=np.int64)
    d = 0
    for i in range(n - 1, -1, -1):
        v = pos_a[order_b[i]] + 1
        # count already inserted values smaller than v (they come later in order_b but earlier in order a)
        j = v - 1
        s = 0
        while j > 0:
            s += tree[j]
            j -= j & (-j)
        d += s
        j = v
        while j <= n:
            tree[j] += 1
            j += j & (-j)
    return d


@njit(cache=True)
def nearest_affine(order, p):
    best = p * p
    pos = np.empty(p, dtype=np.int64)
    for a in range(1, p):
        for b in range(p):
            for i in range(p):
                pos[(a * i + b) % p] = i          # position of vertex in the affine order
            d = kendall(pos, order)
            m = p * (p - 1) // 2
            d = min(d, m - d)                     # the reversed affine order
            if d < best:
                best = d
    return best


def main():
    rows = []
    for name in ('paley_ranking_local_search_20261017.json', 'paley_ranking_local_search_p132_500_20261017.json',
                 'paley_ranking_local_search_p501_1000_20261017.json'):
        path = os.path.join(RES, name)
        if not os.path.exists(path):
            continue
        for r in json.load(open(path)):
            p = r['p']
            order = np.array(r['order'], dtype=np.int64)
            d = nearest_affine(order, p)
            rows.append(dict(p=p, L=r['L'], beaten=r['beaten'], nearest_affine_fraction=d / (p * (p - 1) / 2)))
    rng = np.random.default_rng(1)
    for p in (43, 163, 499, 991):
        ctrl = [nearest_affine(rng.permutation(p).astype(np.int64), p) / (p * (p - 1) / 2) for _ in range(3)]
        print('control p = %d: random orders, nearest-affine fraction %s' % (p, ['%.3f' % c for c in ctrl]))
    print('    p     L    natural beaten   nearest-affine distance of the best order (fraction of pairs)')
    for r in sorted(rows, key=lambda r: r['p']):
        print('%5d %6.3f   %-6s   %.3f' % (r['p'], r['L'], 'YES' if r['beaten'] else 'no', r['nearest_affine_fraction']))
    json.dump(rows, open(os.path.join(RES, 'paley_best_orders_vs_affine_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
