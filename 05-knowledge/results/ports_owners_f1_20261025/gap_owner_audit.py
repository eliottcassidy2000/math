"""Independent degree/weight quotient of generated next-band records.

Imports the disclosed record generators, never their cell-owner canonicalizer
or strongness checker. This audit partitions by full weighted isomorphism
using weight and full outdegree only, then exhausts every cell permutation.
It also evaluates homogeneous pairs by literal third-vertex comparisons.
"""
from collections import defaultdict, Counter
from itertools import combinations, permutations, product
import json
import gap_check as G


CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def clone_pairs(rows):
    n = len(rows)
    return tuple((a, b) for a, b in combinations(range(n), 2)
                 if all(((rows[a] >> v) & 1) == ((rows[b] >> v) & 1)
                        for v in range(n) if v not in (a, b)))


def canonical(rows, weights):
    groups = defaultdict(list)
    n = len(rows)
    size = (n+7)//8
    for i, w in enumerate(weights):
        groups[(w, sum((rows[i] >> j) & 1 for j in range(n)))].append(i)
    best = None
    for blocks in product(*(permutations(groups[x]) for x in sorted(groups))):
        word = tuple(v for block in blocks for v in block)
        out = bytes(weights[v] for v in word)
        for v in word:
            r = 0
            for j, u in enumerate(word):
                if rows[v] >> u & 1:
                    r |= 1 << j
            out += r.to_bytes(size, 'little')
        if best is None or out < best:
            best = out
    return best


def strong(rows):
    n = len(rows)
    for converse in (False, True):
        reached = {0}
        todo = [0]
        while todo:
            a = todo.pop()
            for b in range(n):
                arc = (rows[b] >> a) & 1 if converse else (rows[a] >> b) & 1
                if arc and b not in reached:
                    reached.add(b)
                    todo.append(b)
        if len(reached) != n:
            return False
    return True


def records(k, delta, isstrong, boundary_only=False):
    for overlap in (False, True):
        if boundary_only and overlap:
            continue
        residual = delta-(k+1 if overlap else 0)
        if residual < 0:
            continue
        t = 2*k-1 if overlap else 2*k
        edgesets = G.overlap_graphs(k) if overlap else G.matchings(tuple(range(t)))
        g = t-1 if isstrong else t+1
        for edges in edgesets:
            if isstrong and not G.indecomp(edges, t):
                continue
            if boundary_only:
                ds = []
                for zero in range(g):
                    for unit in range(g):
                        if zero == unit:
                            continue
                        x = [0]*g
                        x[zero], x[unit] = k, k-1
                        ds.append(tuple(x))
            else:
                ds = G.compositions(residual, g, k-1 if overlap else k)
            for deficits in ds:
                rows, weights, gaps = G.make_record(k, edges, deficits, isstrong, overlap)
                yield rows, weights, edges, gaps, overlap


def run(k, delta, isstrong, boundary_only=False):
    seen = set()
    counts = Counter()
    for rows, weights, edges, gaps, overlap in records(k, delta, isstrong, boundary_only):
        counts['raw'] += 1
        cp = clone_pairs(rows)
        t = 2*k-1 if overlap else 2*k
        bad = not overlap and any(gaps[i] == 0 and (i-1, i) in edges for i in range(1, t))
        need(bool(cp) == bad, 'literal clone iff, independent third-vertex code')
        if cp:
            counts['rejected'] += 1
            continue
        if isstrong:
            need(strong(rows), 'independent directed reachability')
        key = canonical(rows, weights)
        need(key not in seen, 'no full weighted-isomorphism collision')
        seen.add(key)
        if len(seen) <= 8:
            order = tuple(reversed(range(len(rows))))
            rr, ww = G.relabel(rows, weights, order)
            need(canonical(rr, ww) == key, 'independent quotient label covariance')
    if boundary_only:
        g = 2*k-1 if isstrong else 2*k+1
        I, _ = G.sequences(k)
        expected = ((I[k]*g*(g-1)-(2*k-3)*I[k-1]*(g-1)) if isstrong
                    else G.df(k)*(g-1)*(g-1))
    else:
        expected = G.formula(k, delta, isstrong)
    need(len(seen) == expected, 'independent exact class count')
    return {'k': k, 'delta': delta, 'strong': isstrong,
            'scope': 'one missing gap and one unit gap' if boundary_only else 'complete structural layer',
            'classes': len(seen), **dict(counts)}


def main():
    layers = [run(3, d, s) for d in (4, 5) for s in (False, True)]
    boundaries = [run(4, 7, s, boundary_only=True) for s in (False, True)]
    print(json.dumps({'status': 'FINITE-EXACT independent weighted-isomorphism and clone audit',
                      'layers': layers, 'boundary_cells': boundaries, 'checks': CHECKS}, indent=2))


if __name__ == '__main__':
    main()
