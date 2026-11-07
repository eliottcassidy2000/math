#!/usr/bin/env python3
import functools
import itertools
import json


def require(b, s):
    if not b:
        raise RuntimeError(s)


def graph(m, q, t, plus):
    n = m + 4
    a, b = n - 2, n - 1
    r = [0] * n
    for i, j in itertools.combinations(range(n), 2):
        if j < a:
            u, v = i, j
        elif j == a:
            u, v = (a, i) if i == 1 else (i, a)
        elif i == a:
            u, v = (a, b) if plus else (b, a)
        else:
            into = 0 if i == 0 else (t if i == 1 else q >> (i - 2) & 1)
            u, v = (i, b) if into else (b, i)
        r[u] |= 1 << v
    return r


def path(r, p):
    return all(r[u] >> v & 1 for u, v in zip(p, p[1:]))


def hp(r):
    full = (1 << len(r)) - 1
    @functools.lru_cache(None)
    def completions(last, remaining):
        if not remaining:
            return 1
        z, available = 0, r[last] & remaining
        while available:
            b = available & -available
            available -= b
            z += completions(b.bit_length()-1, remaining ^ b)
        return z
    return sum(completions(v, full ^ (1 << v)) for v in range(len(r)))


def run():
    counts = dict(adjacent_run_cases=0, dp_pairs=0, permutation_pairs=0)
    levels = []
    for m in range(1, 9):
        deltas = []
        for q in range(1 << m):
            rev = sum((q >> i & 1) << (m - 1 - i) for i in range(m))
            e0, et = q & 1, q >> (m - 1) & 1
            for t in (0, 1):
                plus, minus = graph(m,q,t,True), graph(m,q,t,False)
                a, b = m + 2, m + 3
                ab = ba = 0
                for mask in range(1 << (m + 2)):
                    left = [i for i in range(m+2) if mask >> i & 1]
                    right = [i for i in range(m+2) if not mask >> i & 1]
                    ab += path(plus, left+[a,b]+right)
                    ba += path(minus, left+[b,a]+right)
                require(ab == (4-t)*2**m-2-rev+e0, ('ab',m,q,t,ab))
                require(ba == q+et, ('ba',m,q,t,ba))
                delta = (4-t)*2**m-2-q-rev+e0-et
                require(ab-ba == delta, 'difference')
                require((2-t)*2**m <= delta <= (4-t)*2**m-2, 'bounds')
                require((delta == (2-t)*2**m) == (q == 2**m-1), 'lower equality')
                require((delta == (4-t)*2**m-2) == (q == 0), 'upper equality')
                if m <= 6:
                    require(hp(plus)-hp(minus) == delta, ('DP',m,q,t))
                    counts['dp_pairs'] += 1
                if m <= 2:
                    vals = [sum(path(r,p) for p in itertools.permutations(range(m+4))) for r in (plus,minus)]
                    require(vals[0]-vals[1] == delta, 'brute')
                    counts['permutation_pairs'] += 1
                counts['adjacent_run_cases'] += 1
                deltas.append(delta)
        levels.append(dict(m=m, min_delta=min(deltas), max_delta=max(deltas), parameter_pairs=len(deltas)))
    return dict(status='PASS: exact family arc-flip tax',counts=counts,levels=levels)


if __name__ == '__main__':
    print(json.dumps(run(),indent=2,sort_keys=True))
