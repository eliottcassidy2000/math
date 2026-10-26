"""Independent literal-census audit using relabelling orbits, not canonical labels."""
from collections import deque
from itertools import combinations
from math import factorial
from pathlib import Path
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def orbit_partition(n):
    edges = list(combinations(range(n), 2))
    index = {e: k for k, e in enumerate(edges)}
    generators = []
    for r in range(n - 1):
        p = list(range(n))
        p[r], p[r + 1] = p[r + 1], p[r]
        generators.append(tuple((index[tuple(sorted((p[i], p[j])))], int(p[i] > p[j]))
                                for i, j in edges))
    owner = [-1] * (1 << len(edges))
    sizes = []
    for seed in range(len(owner)):
        if owner[seed] >= 0:
            continue
        ident = len(sizes)
        owner[seed] = ident
        queue = deque([seed])
        count = 0
        while queue:
            t = queue.popleft()
            count += 1
            for g in generators:
                u = sum((((t >> k) & 1) ^ flip) << j for j, (k, flip) in enumerate(g))
                if owner[u] < 0:
                    owner[u] = ident
                    queue.append(u)
        sizes.append(count)
    return edges, owner, sizes


def rows_for(n, edges, t):
    rows = [0] * n
    for k, (i, j) in enumerate(edges):
        a, b = (i, j) if (t >> k) & 1 else (j, i)
        rows[a] |= 1 << b
    return rows


def maximal_transitive_modules(rows):
    n = len(rows)
    full = (1 << n) - 1
    modules = []
    for mask in range(1, full + 1):
        inside = [v for v in range(n) if mask >> v & 1]
        outside = [v for v in range(n) if not mask >> v & 1]
        if not all(rows[v] & mask in (0, mask) for v in outside):
            continue
        if sorted((rows[v] & mask).bit_count() for v in inside) != list(range(len(inside))):
            continue
        modules.append(mask)
    return [s for s in modules if not any(s != t and s & t == s for t in modules)]


def main():
    here = Path(__file__).resolve().parent
    bank = json.loads((here.parent / 'enumeration/representatives.json').read_text())['members']
    cached = {n: orbit_partition(n) for n in range(7)}
    result = []
    for n in range(7):
        edges, owner, sizes = cached[n]
        records = bank[str(n)]
        observed = [owner[r['code']] for r in records]
        need(len(observed) == len(set(observed)) == len(sizes), 'one representative per full relabelling orbit')
        for r in records:
            need(r['aut'] * sizes[owner[r['code']]] == factorial(n), 'independent automorphism mass')
            rows = rows_for(n, edges, r['code'])
            modules = maximal_transitive_modules(rows)
            need(sorted(s.bit_count() for s in modules) == sorted(r['weights']), 'independent chain lengths')
            q = len(modules)
            need(q == r['core_order'], 'independent core order')
            representatives = [(s & -s).bit_length() - 1 for s in modules]
            core = sum(((rows[representatives[i]] >> representatives[j]) & 1) << k
                       for k, (i, j) in enumerate(cached[q][0]))
            need(cached[q][1][core] == cached[q][1][r['core']], 'independent quotient isomorphism')
        result.append({'n': n, 'labelled_universe': len(owner), 'classes': len(sizes),
                       'module_subsets_checked_per_representative': (1 << n) - 1})
    print(json.dumps({'status': 'PASS', 'method': 'adjacent-vertex-transposition orbit BFS and all-subset transitive modules',
                      'scope': 'complete independent literal class and core audit through order six',
                      'rows': result}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
