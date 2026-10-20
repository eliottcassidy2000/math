"""Exact unit-gap boundary with a weight-preserving canonical quotient."""
from collections import Counter
from itertools import permutations, product
from math import comb
import json
import zero_heavy_mass as Z


def canonical(t, w):
    # Weight >1 is intrinsic. Its induced tournament is transitive; fix it
    # by internal scores, then allow every permutation inside each unit cell.
    big = [i for i, x in enumerate(w) if x > 1]
    fixed_big = tuple(big)
    big.sort(key=lambda i: -sum(t[i] >> j & 1 for j in fixed_big))
    Z.need(all(t[a] >> b & 1 for i, a in enumerate(big) for b in big[i+1:]),
           'intrinsic big-vertex order')
    cells = {}
    for i, x in enumerate(w):
        if x == 1:
            signature = tuple(t[i] >> b & 1 for b in big)
            cells.setdefault(signature, []).append(i)
    best = None
    for blocks in product(*(permutations(cells[s]) for s in sorted(cells))):
        order = big + [i for block in blocks for i in block]
        candidate = (tuple(w[i] for i in order),
                     tuple(sum((t[i] >> j & 1) << p for p, j in enumerate(order))
                           for i in order))
        if best is None or candidate < best:
            best = candidate
    return best, max(map(len, cells.values()), default=0)


def matching_gap_bijection(k):
    images = set()
    for m in Z.matchings(tuple(range(2*k))):
        for gap in range(2*k+1):
            shift = lambda i: i + (i >= gap)
            pairs = tuple(sorted(tuple(sorted((shift(a), shift(b)))) for a, b in m))
            image = tuple(sorted(pairs + ((gap, 2*k+1),)))
            images.add(image)
            p = next(pair for pair in image if 2*k+1 in pair)
            Z.need(p[0] == gap, 'new last endpoint retains marked gap')
            inv = lambda i: i - (i > gap)
            restored = tuple(sorted((inv(a), inv(b)) for a, b in image if (a, b) != p))
            Z.need(restored == tuple(sorted(m)), 'inverse of matching-gap bijection')
    expected = set(tuple(sorted(m)) for m in Z.matchings(tuple(range(2*k+2))))
    Z.need(images == expected, 'bijection onto all next matchings')
    return len(images)


def main():
    layers = []; relabel_controls = 0
    for k in (2, 3, 4):
        delta = k-1
        ms = list(Z.matchings(tuple(range(2*k))))
        for isstrong in (False, True):
            g = 2*k-1 if isstrong else 2*k+1
            seen = set(); cells = Counter(); count = 0
            for m in ms:
                if isstrong and not Z.indecomp(m, k):
                    continue
                for deficits in Z.compositions(delta, g):
                    t, w = Z.model(k, m, tuple(k-d for d in deficits), isstrong)
                    Z.need(len(Z.old.G.clone_chains(t)) == len(t), 'pair-free unit-gap boundary')
                    if isstrong:
                        Z.need(Z.strong(t), 'retained strong predicate')
                    key, cell = canonical(t, w)
                    if count < 25:
                        orders = [tuple(reversed(range(len(t)))),
                                  tuple(range(1, len(t)))+(0,)]
                        for order in orders:
                            tt = tuple(sum((t[i] >> j & 1) << p
                                           for p, j in enumerate(order)) for i in order)
                            ww = tuple(w[i] for i in order)
                            Z.need(canonical(tt, ww)[0] == key, 'literal relabelling invariant')
                            relabel_controls += 1
                    Z.need(key not in seen, 'no quotient collision between records')
                    seen.add(key); cells[cell] += 1; count += 1
            base = sum(Z.indecomp(m, k) for m in ms) if isstrong else len(ms)
            Z.need(count == base*comb(g+delta-1, delta), 'boundary count')
            layers.append(dict(k=k, delta=delta, strong=isstrong,
                               classes=count, max_unit_cell_histogram=dict(cells)))
    print(json.dumps(dict(status='FINITE-EXACT unit-gap quotient controls',
                          layers=layers,
                          relabel_controls=relabel_controls,
                          marked_gap_bijection={k: matching_gap_bijection(k) for k in range(2,5)}),
                     indent=2))


if __name__ == '__main__':
    main()
