"""Independent literal word/reversal audit; imports no producer checker."""
from itertools import combinations
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def rev(q):
    return dict(zip(q, reversed(q)))


def edge(a, b):
    return tuple(sorted((a, b)))


def edges(q):
    return {edge(a, b) for a, b in zip(q, q[1:])}


def main():
    antipodal = endcaps = 0
    for n in range(12, 513):
        p = tuple(range(n))
        rp = rev(p)
        d2 = (2, 0) + tuple(range(4, n)) + (1, 3)
        r = rev(d2)
        g = {v: r[rp[v]] for v in range(n)}
        expected = {v: v + 4 for v in range(n - 4)}
        expected.update({n - 4: 2, n - 3: 3, n - 2: 0, n - 1: 1})
        need(g == expected, 'explicit reversal product')
        if n % 2 == 0:
            unseen = set(range(n))
            cycles = []
            while unseen:
                cycle = []
                v = min(unseen)
                while v in unseen:
                    unseen.remove(v)
                    cycle.append(v)
                    v = g[v]
                cycles.append(cycle)
            defects = {edge(0, 1), edge(3, 4), edge(0, 4), edge(1, n - 1)}
            need(all(edge(a, b) not in defects or
                     (cycle.index(a) - cycle.index(b)) % len(cycle) in (1, len(cycle) - 1)
                     for cycle in cycles for a, b in combinations(cycle, 2)),
                 'internal defects adjacent')
            cycle = next(c for c in cycles if len(c) % 2 == 0 and len(c) >= 4)
            a, b = cycle[0], cycle[len(cycle) // 2]
            start = (a, b)
            for _ in range(len(cycle) // 2):
                need(edge(rp[a], rp[b]) not in defects, 'reflected antipodes avoid defect')
                a, b = g[a], g[b]
                need(edge(a, b) not in defects, 'antipodes avoid defect')
            need((a, b) == tuple(reversed(start)), 'half orbit reverses endpoints')
            antipodal += 1
        w = tuple(range(4, n - 3))
        right = (n - 3, n - 2, n - 1)
        words = ((0, 1) + w + right + (2, 3),
                 (1,) + w + (0,) + right + (2, 3),
                 (1, 2, 3) + w + (0,) + right)
        need(edges(words[0]) ^ edges(words[1]) ^ edges(words[2]) == edges(p),
             'endcap exact target')
        witnesses = ((0, 2, 1, 3), (1, 2, 4, 3), (0, n - 1, 1, 4))
        for q, (a, b, c, d) in zip(words, witnesses):
            r = rev(q)
            need((r[b], r[a]) == (c, d), 'literal endcap witness positions')
        endcaps += 1
    print(json.dumps({'status': 'PASS', 'independent_even_antipodal_orders': antipodal,
                      'endcap_word_orders': endcaps, 'order_range': [12, 512]}, indent=2))


if __name__ == '__main__':
    main()
