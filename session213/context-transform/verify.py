#!/usr/bin/env python3
"""Self-contained exact controls for the cycle-code path-density proof."""
import itertools
import json


def require(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def comb(m, q, tau, plus):
    n = m + 4
    rows = [0] * n
    a, b = n - 2, n - 1
    def edge(i, j):
        rows[i] |= 1 << j
    for i in range(m + 2):
        for j in range(i + 1, m + 2):
            edge(i, j)
        edge(a, i) if i == 1 else edge(i, a)
        inbound = 0 if i == 0 else (tau if i == 1 else q >> (i - 2) & 1)
        edge(i, b) if inbound else edge(b, i)
    edge(a, b) if plus else edge(b, a)
    return rows


def induced(rows, vertices):
    return [sum(((rows[i] >> j) & 1) << k for k, j in enumerate(vertices)) for i in vertices]


def all_subset_paths(rows, root=None):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        if root is None or v == root:
            dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            ways = dp[mask][v]
            if not ways:
                continue
            available = rows[v] & ~mask
            while available:
                bit = available & -available
                available -= bit
                dp[mask | bit][bit.bit_length() - 1] += ways
    counts = [sum(row) for row in dp]
    counts[0] = 1
    return counts, dp[-1]


def hc(rows):
    _, ends = all_subset_paths(rows, 0)
    return sum(z for v, z in enumerate(ends) if rows[v] & 1)


def g(v):
    if v == 0:
        return 1
    u = v // 2
    return g(u) if v % 2 else g(u) + u + u % 2


def is_path(rows, word):
    return len(word) == len(rows) and set(word) == set(range(len(rows))) and all(rows[i] >> j & 1 for i, j in zip(word, word[1:]))


def path_classes(m, q, tau):
    a, b = m + 2, m + 3
    first, second = set(), set()
    for labels in itertools.product(range(3), repeat=m):
        I, J, K = ([i + 2 for i, lab in enumerate(labels) if lab == v] for v in range(3))
        if not I or q >> (I[-1] - 2) & 1:
            first.add(tuple(I + [b, 0] + J + [a, 1] + K))
        if (tau if not J else q >> (J[-1] - 2) & 1):
            second.add(tuple(I + [a, 1] + J + [b, 0] + K))
    return first, second


def budget_count(X):
    """Four presentations per (m,q), using only the theorem's necessary budget."""
    n, r = 0, 0
    while 2 ** (r + 1) <= X:
        j = 0
        while 2 ** (r + 1) * 3 ** j <= X:
            n += 4 * 2 ** j
            j += 1
        r += 1
    return n


def run():
    records = []
    totals = dict(dp_objects=0, exact_minus_readouts=0, literal_path_objects=0,
                  literal_paths=0, permutation_objects=0, endpoint_deletions=0)
    observed = []
    for m in range(1, 8):
        vals = []
        for q in range(1 << m):
            V = sum((q >> i & 1) * 3 ** i * 2 ** (m - 1 - i) for i in range(m))
            R = 1 + sum((1 - (q >> i & 1)) * 2 ** (m - 1 - i) for i in range(m))
            e = q >> (m - 1) & 1
            for tau in (0, 1):
                for plus in (False, True):
                    rows = comb(m, q, tau, plus)
                    require(all(not (rows[i] >> i & 1) for i in range(len(rows))), 'loop')
                    require(all(((rows[i] >> j) ^ (rows[j] >> i)) & 1 for i in range(len(rows)) for j in range(i)), 'tournament')
                    paths, ends = all_subset_paths(rows)
                    h = paths[-1]
                    require(h >= (1 + tau) * 2 ** m + 2 * V, ('lower bound', m, q, tau, plus, h))
                    require(h % 2 == 1, 'Redei control')
                    require(hc(rows) == q + tau + plus, ('cycle', m, q, tau, plus))
                    if q:
                        l = q.bit_length()
                        require(h >= 2 * 3 ** (l - 1) * 2 ** (m - l), 'highest bit')
                    if m <= 4:
                        one, two = path_classes(m, q, tau)
                        require(not one & two, 'class intersection')
                        require(len(one) == 2 ** m + V and len(two) == tau * 2 ** m + V, 'class cardinality')
                        require(all(is_path(rows, p) for p in one | two), 'invalid literal path')
                        totals['literal_path_objects'] += 1
                        totals['literal_paths'] += len(one) + len(two)
                    if m <= 2:
                        brute = sum(is_path(rows, p) for p in itertools.permutations(range(len(rows))))
                        require(brute == h, ('permutations', m, q, tau, plus, h, brute))
                        totals['permutation_objects'] += 1
                    if not plus:
                        U = induced(rows, list(range(m + 2)) + [m + 3])
                        W = induced(rows, list(range(2, m + 2)) + [m + 3])
                        ups, uends = all_subset_paths(U)
                        wps, _ = all_subset_paths(W)
                        full = (1 << len(W)) - 1
                        p2 = sum(wps[s] * wps[full ^ s] for s in range(full + 1))
                        require(p2 % 2 == 0, 'cover parity')
                        require(uends[1] == (2 - tau) * e, 'ending at x1')
                        predicted = ups[-1] - (2 - tau) * e + (1 + tau) * 2 ** m + 2 * V + p2 - (2 - tau) * R
                        require(h == predicted, ('exact minus', m, q, tau, h, predicted))
                        totals['exact_minus_readouts'] += 1
                        if tau:
                            require(ups[-1] == 2 * (q + 1) + g(q), 'deletion law')
                            require(2 * q + 3 <= ups[-1] <= 3 * (q + 1), 'deletion bounds')
                            totals['endpoint_deletions'] += 1
                        if tau and q == (1 << m) - 1:
                            require(h == 2 * 3 ** m + 2 ** (m + 2) - 1, 'inherited all-one ray')
                        if tau and not q:
                            require(h == 3 * 2 ** m + 3, 'zero ray')
                    vals.append(h)
                    observed.append((m, q, h))
                    totals['dp_objects'] += 1
        records.append(dict(m=m, order=m+4, parameter_cases=len(vals), distinct_H=len(set(vals)), min_H=min(vals), max_H=max(vals)))
    grids = []
    for X in (2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152):
        direct = 0
        # Exact independent enumeration by m then q bit-length, including padding.
        for m in range(1, X.bit_length() + 1):
            for l in range(1, m + 1):
                if 2 * 3 ** (l - 1) * 2 ** (m - l) <= X:
                    direct += 4 * 2 ** (l - 1)
        require(direct == budget_count(X), 'budget summation mismatch')
        # Avoid floating alpha in any pass/fail test: the exact integer region is authoritative.
        require(sum(h <= X and q > 0 for m, q, h in observed) <= direct, 'finite cases exceed infinite budget')
        grids.append(dict(X=X, exact_nonzero_parameter_budget=direct, zero_parameter_budget=4*(X.bit_length()-1)))
    return dict(status='FINITE-EXACT controls; all-order density follows from PROOF.md',
                totals=totals, levels=records, budget_grid=grids,
                exact_deletion_hostile=dict(target=11, candidates=[2*k+g(k-1) for k in range(1,6)]))


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
