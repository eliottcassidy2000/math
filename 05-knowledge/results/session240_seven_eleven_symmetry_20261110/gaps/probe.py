#!/usr/bin/env python3
"""Exact cyclic universe, two-cover identity, independent non-DP path DFS.
All checks are explicit and survive python -O. No external dependencies.
"""
from collections import Counter
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def path_table(out):
    n = len(out)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    h = [0] * (full + 1)
    h[0] = 1
    for v in range(n):
        dp[1 << v][v] = 1
    for s in range(1, full + 1):
        for v in range(n):
            c = dp[s][v]
            if not c:
                continue
            h[s] += c
            rem = out[v] & (full ^ s)
            while rem:
                b = rem & -rem
                rem -= b
                dp[s | b][b.bit_length() - 1] += c
    return h


def dfs_count(out):
    """Enumerate directed simple paths directly: no subset-count DP or cache."""
    n = len(out)
    neighbors = [[w for w in range(n) if out[v] >> w & 1] for v in range(n)]
    used = [False] * n
    def visit(v, depth):
        if depth == n:
            return 1
        count = 0
        for w in neighbors[v]:
            if not used[w]:
                used[w] = True
                count += visit(w, depth + 1)
                used[w] = False
        return count
    total = 0
    for v in range(n):
        used[v] = True
        total += visit(v, 1)
        used[v] = False
    return total


def validate(out):
    n = len(out)
    for x in range(n):
        check(not (out[x] >> x & 1), 'loop')
        for y in range(x):
            check(((out[x] >> y) & 1) + ((out[y] >> x) & 1) == 1, 'not tournament')
    for root in range(n):
        seen, pending = {root}, [root]
        while pending:
            v = pending.pop()
            for w in range(n):
                if out[v] >> w & 1 and w not in seen:
                    seen.add(w)
                    pending.append(w)
        check(len(seen) == n, 'not strong')


def circulant(p, mask):
    D = [i if not mask >> (i - 1) & 1 else p - i for i in range(1, (p + 1) // 2)]
    return [sum(1 << ((x + d) % p) for d in D) for x in range(p)], D


def inflate(out):
    p = len(out)
    # R -> x -> y -> R
    return [a | (1 << p) for a in out] + [1 << (p + 1), (1 << p) - 1]


def run():
    result = {}
    expected = {
        3: {(3, 6): 2}, 5: {(15, 90): 4},
        7: {(175, 1638): 6, (189, 1638): 2},
        11: {(92411, 1571526): 10, (93027, 1574958): 10,
             (93467, 1577598): 10, (95095, 1586310): 2},
    }
    for p in [3, 5, 7, 11]:
        rows = []
        for mask in range(1 << ((p - 1) // 2)):
            out, D = circulant(p, mask)
            validate(out)
            for x in range(p):
                for y in range(p):
                    check((out[x] >> y & 1) == (out[(x + 1) % p] >> ((y + 1) % p) & 1), 'rotation')
            hs = path_table(out)
            full = (1 << p) - 1
            h = hs[-1]
            f2 = sum(hs[s] * hs[full ^ s] for s in range(1, full))
            inf = inflate(out)
            validate(inf)
            g = path_table(inf)[-1]
            check(g == 3 * h + f2, f'two-cover identity p={p} mask={mask}')
            check(h % (2 * p) == p and g % (2 * p) == p, 'odd and p-divisible')
            rows.append({'mask': mask, 'D': D, 'H': h, 'F2': f2, 'inflation_H': g})
        check(Counter((r['H'], r['F2']) for r in rows) == expected[p], 'exact census changed')
        # Multiplication x -> a*x gives explicit isomorphisms. No assumption
        # that multipliers characterize isomorphism is needed or made.
        by_D = {tuple(sorted(r['D'])): r['mask'] for r in rows}
        remaining = set(range(len(rows)))
        multiplier_orbits = []
        while remaining:
            rep = min(remaining)
            orbit = sorted({by_D[tuple(sorted(a * d % p for d in rows[rep]['D']))]
                            for a in range(1, p)})
            check(set(orbit) <= remaining, 'multiplier orbit overlap')
            remaining.difference_update(orbit)
            r = rows[rep]
            check(all((rows[m]['H'], rows[m]['F2']) == (r['H'], r['F2']) for m in orbit), 'relabel invariant')
            out, _ = circulant(p, rep)
            check(dfs_count(out) == r['H'], 'independent DFS core')
            multiplier_orbits.append({'representative': rep, 'masks': orbit, 'core_DFS': r['H']})
        mu = min(r['H'] for r in rows)
        gamma = min(r['inflation_H'] for r in rows)
        minimizers = [r['mask'] for r in rows if r['inflation_H'] == gamma]
        check(any(o['masks'] == minimizers for o in multiplier_orbits), 'minimum not a single multiplier orbit')
        result[str(p)] = {'orientations': len(rows), 'mu': mu, 'gamma': gamma,
            'multi_orbit_floor': mu * mu + 2 * p,
            'histogram': dict(sorted(Counter(r['H'] for r in rows).items())),
            'minimizers': minimizers, 'multiplier_orbits': multiplier_orbits, 'rows': rows}
    out, _ = circulant(11, 1)
    check(dfs_count(inflate(out)) == 1848759, 'independent sharp inflation DFS')
    # Hostile: non-strong padding retains H, including at arbitrarily large
    # order by ordinal-sum multiplicativity. This checks one padding vertex.
    for p, mask in [(7, 0), (11, 1)]:
        out, _ = circulant(p, mask)
        check(path_table([a | (1 << p) for a in out] + [0])[-1] == result[str(p)]['mu'], 'padding hostile')
    check(path_table([6, 4, 0])[-1] == 1, 'transitive control')
    check(path_table([2, 4, 1])[-1] == 3, 'C3 control')
    print(json.dumps({'status': 'PASS', 'checks': {'independent_minimum_13_vertex_DFS': 1848759,
        'core_DFS': 'one representative per multiplier orbit at p=3,5,7,11',
        'padding_hostiles': [7, 11], 'two_cover_identity': 'all 46 orientations'}, 'primes': result},
        indent=2, sort_keys=True))


if __name__ == '__main__':
    run()
