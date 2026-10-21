#!/usr/bin/env python3
"""Exact retained-vertex records at two boundaries; no removable assertions."""
from collections import Counter
from itertools import combinations, permutations, product
from math import comb
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def matchings(xs):
    if not xs:
        yield ()
        return
    a = xs[0]
    for b in xs[1:]:
        for rest in matchings(tuple(x for x in xs[1:] if x != b)):
            yield ((a, b),) + rest


def compositions(total, parts, cap):
    if parts == 0:
        if total == 0:
            yield ()
        return
    if not 0 <= total <= parts * cap:
        return
    for a in range(max(0, total-(parts-1)*cap), min(cap, total)+1):
        for rest in compositions(total-a, parts-1, cap):
            yield (a,) + rest


def p3_graphs(k):
    xs = tuple(range(2*k-1))
    for triple in combinations(xs, 3):
        for centre in triple:
            star = tuple(tuple(sorted((centre, v))) for v in triple if v != centre)
            for rest in matchings(tuple(v for v in xs if v not in triple)):
                yield tuple(sorted(star+rest))


def indecomp(edges, n):
    return all(any(a < c <= b for a, b in edges) for c in range(1, n))


def model(n, edges, gaps, strong):
    if strong:
        gaps = (0,) + gaps + (0,)
    need(len(gaps) == n+1, 'one retained slot per actual gap')
    slots = []; w = []
    for i, g in enumerate(gaps):
        if g:
            w.append(g)
        if i < n:
            slots.append(len(w)); w.append(1)
    q = len(w)
    t = [((1 << q)-1) ^ ((1 << (i+1))-1) for i in range(q)]
    for a, b in edges:
        i, j = slots[a], slots[b]
        t[i] ^= 1 << j; t[j] ^= 1 << i
    return tuple(t), tuple(w)


def clones(t):
    full = (1 << len(t))-1
    return [(i, j) for i, j in combinations(range(len(t)), 2)
            if not ((t[i] ^ t[j]) & (full ^ ((1 << i) | (1 << j))))]


def strongly_connected(t):
    def closure(rows):
        seen = 1
        while True:
            new = seen
            for i, row in enumerate(rows):
                if seen >> i & 1:
                    new |= row
            if new == seen:
                return seen == (1 << len(rows))-1
            seen = new
    rev = tuple(sum((t[j] >> i & 1) << j for j in range(len(t))) for i in range(len(t)))
    return closure(t) and closure(rev)


def canonical(t, w):
    big = tuple(i for i, a in enumerate(w) if a > 1)
    ordered = sorted(big, key=lambda i: -sum(t[i] >> j & 1 for j in big))
    need(all(t[i] >> j & 1 for p, i in enumerate(ordered) for j in ordered[p+1:]),
         'large-vertex order is truly transitive')
    cells = {}
    for i, a in enumerate(w):
        if a == 1:
            sig = tuple(t[i] >> j & 1 for j in ordered)
            cells.setdefault(sig, []).append(i)
    best = None
    for blocks in product(*(permutations(cells[s]) for s in sorted(cells))):
        order = ordered + [i for block in blocks for i in block]
        key = (tuple(w[i] for i in order), tuple(sum((t[i] >> j & 1) << p
               for p, j in enumerate(order)) for i in order))
        if best is None or key < best:
            best = key
    return best, max(map(len, cells.values()), default=0)


def fas(t, w):
    q = len(t); dp = [0] + [10**9]*((1 << q)-1)
    for mask in range(1, 1 << q):
        for i in range(q):
            if mask >> i & 1:
                rem = mask ^ (1 << i)
                cost = sum(w[i]*w[j] for j in range(q) if rem >> j & 1 and t[i] >> j & 1)
                dp[mask] = min(dp[mask], dp[rem]+cost)
    return dp[-1]


def records(k, delta, strong):
    n = 2*k; m = n-1 if strong else n+1
    for edges in matchings(tuple(range(n))):
        if strong and not indecomp(edges, n):
            continue
        for ds in compositions(delta, m, k-1):
            yield 'A', edges, tuple(k-d for d in ds), n
        for hole in range(m):
            for ds in compositions(delta-k, m-1, k-1):
                gaps = list(k-d for d in ds); gaps.insert(hole, 0)
                yield 'B', edges, tuple(gaps), n
    if delta == k+1:
        n = 2*k-1; m = n-1 if strong else n+1
        for edges in p3_graphs(k):
            if not strong or indecomp(edges, n):
                yield 'C', edges, (k,)*m, n


def main():
    F = [1]; I = [0]; P = [0, 0]
    for k in range(1, 7):
        F.append(F[-1]*(2*k-1))
        I.append(F[k]-sum(I[j]*F[k-j] for j in range(1, k)))
        if k >= 2:
            P.append((k-1)*F[k]-sum(P[j]*sum(F[a]*F[k-j-a] for a in range(k-j+1))
                                      for j in range(2, k)))
    independent_graph_checks = []
    for k in range(2, 6):
        ms = list(matchings(tuple(range(2*k))))
        ps = list(p3_graphs(k))
        need(len(ms) == F[k], 'matching universe')
        need(sum(indecomp(e, 2*k) for e in ms) == I[k], 'irreducible recurrence')
        need(len(ps) == (k-1)*F[k], 'P3 universe')
        need(sum(indecomp(e, 2*k-1) for e in ps) == P[k], 'P3 irreducible recurrence')
        j = sum(sum(b == a+1 for a, b in e) for e in ms if indecomp(e, 2*k))
        need(j == (2*k-3)*I[k-1], 'adjacent-pair incidence')
        independent_graph_checks.append(dict(k=k, matchings=F[k], irreducible=I[k],
                                            p3_graphs=len(ps), irreducible_p3=P[k], adjacent=j))
    rows = []; relabels = 0; dp_checks = 0
    for k in (2, 3, 4):
        for delta in (k, k+1):
            for strong in (False, True):
                seen = set(); accepted = Counter(); rejected = Counter(); cell_sizes = Counter()
                for kind, edges, gaps, n in records(k, delta, strong):
                    t, w = model(n, edges, gaps, strong)
                    pairs = clones(t)
                    if kind == 'B':
                        h = gaps.index(0)+(1 if strong else 0)
                        forbidden = h > 0 and h < n and (h-1, h) in edges
                        need(bool(pairs) == forbidden, 'omitted-gap clone iff')
                        if forbidden:
                            rejected[kind] += 1
                            continue
                    need(not pairs, 'all retained classes pair free')
                    expected_mass = (2*k*k+k if strong else 2*k*k+3*k)-delta
                    need(sum(w) == expected_mass, 'actual total mass')
                    if strong:
                        need(strongly_connected(t), 'literal strong predicate')
                    key, size = canonical(t, w)
                    need(key not in seen, 'different records do not collide')
                    seen.add(key); accepted[kind] += 1; cell_sizes[size] += 1
                    if sum(accepted.values()) <= 3:
                        order = tuple(reversed(range(len(t))))
                        tt = tuple(sum((t[i] >> j & 1) << p for p, j in enumerate(order)) for i in order)
                        ww = tuple(w[i] for i in order)
                        need(canonical(tt, ww)[0] == key, 'relabel covariance'); relabels += 1
                    if k == 2 and accepted[kind] <= 2:
                        need(fas(t, w) == k, 'independent weighted feedback optimum'); dp_checks += 1
                if delta == k:
                    expected = I[k]*comb(3*k-2, k)-(2*k-3)*I[k-1] if strong else F[k]*(comb(3*k, k)-1)
                else:
                    expected = (I[k]*(comb(3*k-1, k+1)-(2*k-1))-(2*k-2)*(2*k-3)*I[k-1]+P[k]
                                if strong else F[k]*(comb(3*k+1, k+1)-3*k-2))
                need(len(seen) == expected, 'claimed consequence count')
                rows.append(dict(k=k, delta=delta, strong=strong, classes=len(seen),
                                 types=dict(accepted), clone_rejections=dict(rejected),
                                 max_unit_cell_histogram=dict(cell_sizes)))
    print(json.dumps(dict(status='FINITE-EXACT quotient and hostile controls; proofs separate',
                         graph_checks=independent_graph_checks, boundary_layers=rows,
                         relabel_controls=relabels, weighted_subset_dp_checks=dp_checks,
                         small_k2_full_compiler_reference={'delta2':[42,11], 'delta3':[81,15]}),
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
