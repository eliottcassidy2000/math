"""Exact fixed-skeleton boundary controls; no tournament census.

Certifies complete run signatures, literal transfer/recharge, and the
small radius-two obstruction to a uniform 18h boundary repair.
"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def need(test, message):
    if not test:
        raise RuntimeError(message)


def tt(n):
    return tuple(sum(1 << j for j in range(i + 1, n)) for i in range(n))


def ordinal(A, B):
    a, b = len(A), len(B)
    return tuple(row | (((1 << b) - 1) << a) for row in A) + tuple(row << a for row in B)


@lru_cache(None)
def h_subsets(R):
    n = len(R)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(dp[mask]):
            todo = R[v] & ~mask
            while count and todo:
                bit = todo & -todo
                todo -= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return tuple([1] + [sum(dp[mask]) for mask in range(1, 1 << n)])


def H(R):
    return h_subsets(tuple(R))[-1]


def slots(R, k):
    h = h_subsets(tuple(R))
    p = [0] * len(h)
    p[0] = 1
    for _ in range(k):
        q = []
        for mask in range(len(h)):
            value, sub = 0, mask
            while True:
                value += h[sub] * p[mask ^ sub]
                if sub == 0:
                    break
                sub = (sub - 1) & mask
            q.append(value)
        p = q
    return p[-1]


def skeleton():
    R = list(tt(9))
    for u, v in ((4, 0), (5, 0), (7, 2), (8, 2)):
        R[u] ^= 1 << v
        R[v] ^= 1 << u
    return tuple(R)


def inflate(R, modules):
    owners = [(u, j) for u, M in enumerate(modules) for j in range(len(M))]
    out = [0] * len(owners)
    for i, (u, a) in enumerate(owners):
        for j, (v, b) in enumerate(owners):
            if (modules[u][a] >> b) & 1 if u == v else (R[u] >> v) & 1:
                out[i] |= 1 << j
    return tuple(out)


def pi(M, k):
    modules = [tt(1) for _ in range(9)]
    modules[1], modules[3] = tt(k), M
    return inflate(skeleton(), modules)


def gamma(M):
    return pi(M, 0)


def relabel(R, order):
    return tuple(sum(((R[u] >> v) & 1) << j for j, v in enumerate(order)) for u in order)


def signature(R, root):
    n = len(R)
    other = [v for v in range(n) if v != root]
    full = (1 << len(other)) - 1
    cap = 1 + sum((R[v] >> root) & 1 for v in other)
    result = []
    for visits in range(1, cap + 1):
        @lru_cache(None)
        def count(mask, k, last):
            if mask == full and k == visits:
                return 1
            value = 0
            if k < visits and last != root and (last < 0 or (R[last] >> root) & 1):
                value += count(mask, k + 1, root)
            for bit, v in enumerate(other):
                if not (mask >> bit) & 1 and (last < 0 or (R[last] >> v) & 1):
                    value += count(mask | (1 << bit), k, v)
            return value
        result.append(count(0, 0, -1))
    return result


def main():
    G = gamma(tt(1))
    sig = signature(G, 2)
    need(G == (230, 60, 248, 241, 225, 192, 130, 2), 'literal Gamma rows')
    need(sig == [163, 180, 60], 'complete Gamma run signature')
    modules = [('TT1', tt(1)), ('TT2', tt(2)), ('TT3', tt(3)),
               ('C3', (2, 4, 1)), ('strong4', (2, 4, 9, 3))]
    module_records = []
    for name, M in modules:
        h, p = H(M), slots(M, 3)
        before = pi(M, 1)
        changed = list(before)
        changed[1] ^= 1 << 2
        changed[2] ^= 1 << 1
        after = relabel(changed, [0, 2, 1] + list(range(3, len(changed))))
        enlarged = ordinal(tt(1), M)
        need(after == gamma(enlarged), 'literal last-clone boundary graph')
        need(H(G if M == tt(1) else gamma(M)) == 60 * p - 17 * h,
             'arbitrary-module Gamma identity')
        need(H(before) == 120 * p + h, 'inherited Pi identity')
        need(H(after) == 180 * p - 17 * h, 'exact boundary count')
        desired = (3 * H(before) - h) // 2
        need(desired - H(after) == 18 * h, 'boundary debt exactly 18h')
        # Fresh v is the ordinary root-child birth at a (vertex0).
        n = len(after)
        recharge = tuple(row | ((1 << n) if u == 0 else 0)
                         for u, row in enumerate(after)) + (((1 << n) - 1) ^ 1,)
        canonical = relabel(recharge, [0, n] + list(range(1, n)))
        need(canonical == pi(enlarged, 1), 'literal native recharge at a')
        need(H(canonical) == 2 * H(after) + 35 * h == 3 * H(before) - 2 * h,
             'recharge and complete circuit laws')
        module_records.append(dict(module=name, h=h, P3=p, before=H(before),
                                   after=H(after), desired=desired, debt=18 * h,
                                   recharged=H(canonical)))

    fits = []
    counts = {}
    pairs = list(combinations(range(8), 2))
    for width in (1, 2):
        examined, scalar_fits, full_fits = 0, 0, 0
        for chosen in combinations(pairs, width):
            examined += 1
            R = list(G)
            for u, v in chosen:
                R[u] ^= 1 << v
                R[v] ^= 1 << u
            R = tuple(R)
            if H(R) != 181:
                continue
            scalar_fits += 1
            s = signature(R, 2)
            target = [181, 180, 60] + [0] * max(0, len(s) - 3)
            matches = len(s) >= 3 and s == target
            full_fits += matches
            fits.append(dict(flips=chosen, rows=R, complete_signature=s,
                             universal_repair=matches))
        counts[str(width)] = dict(examined=examined, scalar_fits=scalar_fits,
                                  universal_fits=full_fits)
    need(counts['1']['scalar_fits'] == 0 and counts['2']['scalar_fits'] == 2,
         'named fixed-skeleton scalar fits')
    need(all(row['universal_fits'] == 0 for row in counts.values()),
         'no uniform repair in the stated radius-two universe')
    result = dict(status='PASS', Gamma_rows=G, complete_run_signature=sig,
                  module_controls=module_records, repair_universe=counts,
                  singleton_only_fits=fits,
                  scope='Fixed 8-vertex Gamma and named modules only; no census of all clocks or all tournaments.')
    (HERE / 'clock_boundary_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
