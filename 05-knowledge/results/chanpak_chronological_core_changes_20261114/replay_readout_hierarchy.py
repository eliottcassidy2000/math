"""Literal controls for the chronological budget / 7-adic strip."""
from functools import lru_cache
from math import isqrt


def insist(ok, *ctx):
    if not ok:
        raise AssertionError(ctx)


def covers(rows, r):
    n = len(rows)
    full = (1 << n) - 1

    @lru_cache(None)
    def end(s, i):
        t = s ^ (1 << i)
        return 1 if not t else sum(end(t, j) for j in range(n)
                                   if (t >> j) & (rows[j] >> i) & 1)

    @lru_cache(None)
    def visit(s, k):
        if k == 1:
            return sum(end(s, i) for i in range(n) if s >> i & 1) if s else 1
        t = s
        ans = 0
        while True:
            ans += visit(t, 1) * visit(s ^ t, k - 1)
            if not t:
                return ans
            t = (t - 1) & s

    return visit(full, r)


def strip(R, m):
    n = 4 + len(R) + m
    rows = [0] * n
    bb = n - 2
    for u in range(n):
        for v in range(u + 1, n):
            if u < 2 and v >= bb:
                rows[v] |= 1 << u
            elif 2 <= u < v < 2 + len(R):
                if R[u - 2] >> (v - 2) & 1:
                    rows[u] |= 1 << v
                else:
                    rows[v] |= 1 << u
            else:
                rows[u] |= 1 << v
    return tuple(rows)


def arcs(rows):
    return {(i, j) for i, row in enumerate(rows) for j in range(len(rows)) if row >> j & 1}


def certify(rows, rank, paths):
    E = arcs(rows)
    insist(sorted(rank) == list(range(len(rows))), rank)
    loc = {x: i for i, x in enumerate(rank)}
    backward = {(x, y) for x, y in E if loc[x] > loc[y]}
    cut = set()
    for p in paths:
        insist(sorted(p) == list(range(len(rows))), p)
        e = set(zip(p, p[1:]))
        insist(e <= E, e - E)
        cut |= e
    insist(backward <= cut, backward - cut)


def control_strip(R, rp, rank_r, m):
    rows = strip(R, m)
    n = len(rows)
    r = len(R)
    rank = (0, 1) + tuple(2 + x for x in rank_r) + tuple(range(2 + r, n))
    paths = []
    for s in (0, 1):
        paths.append((n - 2, s) + tuple(2 + x for x in rp) +
                     (n - 1, 1 - s) + tuple(range(2 + r, n - 2)))
    certify(rows, rank, paths)
    insist(min(sum(row >> i & 1 for row in rows) for i in range(n)) >= 2)
    h, p2, p3 = [covers(R, k) for k in (1, 2, 3)]
    answer = covers(rows, 1)
    insist(answer == 4 * p3 * 3**m + 5 * p2 * 2**m - 11 * h, R, m, answer)
    return answer


def main():
    R = (18, 8, 3, 5, 14)
    p = (0, 4, 1, 3, 2)
    ranking = (4, 2, 3, 0, 1)
    certify(R, ranking, [p])
    insist([covers(R, k) for k in (1, 2, 3)] == [13, 104, 471])
    insist([control_strip(R, p, ranking, m) for m in (1, 2, 3)][1] == 7 * 2699)
    prime = 2699
    insist(all(prime % j for j in range(2, isqrt(prime) + 1)))
    for m in range(1, 7):
        insist(control_strip((0,), (0,), (0,), m) == 4 * 3**(m + 1) + 5 * 2**(m + 1) - 11)
    for s in range(1, 6):
        q = 7**s
        mod = 7 * q
        def M(m):
            return (1884 * pow(3, m, mod) + 520 * pow(2, m, mod) - 143) % mod
        insist({M(6*u) // 7 for u in range(q)} == set(range(q)), s)
        insist(all(M(6*u + 2) % 49 == 28 for u in range(q)), s)
    print('Direct Hamilton/cover DP, two-path certificates, and 7-adic controls: PASS')


if __name__ == '__main__':
    main()
