from itertools import permutations
from fractions import Fraction as F
from collections import Counter

def wrap(n):
    edges = {(i, j) for i in range(n) for j in range(i + 1, n)}
    edges |= {(i, n) for i in range(n)} | {(n, n + 1)} | {(n + 1, i) for i in range(n)}
    return [p for p in permutations(range(n + 2)) if all(((a, b) in edges for a, b in zip(p, p[1:])))]
for n in range(2, 7):
    paths = wrap(n)
    v = [F((-1) ** i * (i + 1), i + 3) for i in range(n)]
    incidence = [sum((v[i] for i in range(n) if p.index(i) < p.index(n)), F(0)) for p in paths]
    c = sum(v)
    sq = sum((x * x for x in v))
    h = 2 ** n + 1
    e = sum(incidence) / h
    e2 = sum((x * x for x in incidence)) / h
    if not e == F(2 ** (n - 1) + 1, h) * c:
        raise AssertionError()
    if not e2 == F(2 ** (n - 2), h) * (c * c + sq) + c * c / h:
        raise AssertionError()
    counter = Counter((tuple((i for i in range(n) if p.index(i) < p.index(n))) for p in paths))
    tv = sum((abs(F(counter.get(s, 0), h) - F(1, 2 ** n)) for s in (tuple((i for i in range(n) if m >> i & 1)) for m in range(2 ** n)))) / 2
    if not tv == (1 - F(1, 2 ** n)) / h:
        raise AssertionError()
    for i in range(n):
        deleted = [tuple((x - (x > i) for x in p if x != i)) for p in paths if p.index(i) < p.index(n)]
        if not (set(deleted) == set(wrap(n - 1)) and len(deleted) == len(set(deleted))):
            raise AssertionError()
    print('n', n, 'E', e, 'E2', e2, 'TV subset', tv)
Q = 10 ** 9
N = 200
ps = []
for i in range(2, N + 3):
    T = i * Q ** i
    lo, hi = (Q, 2 * Q)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** i <= T:
            lo = mid
        else:
            hi = mid
    if not lo ** i <= T < hi ** i:
        raise AssertionError()
    ps.append(lo - Q)
D = sum(((-1) ** i * ps[i - 2] for i in range(2, N + 1)))
lo = 2 * D - 198 - 2 * ps[199] - 2 + ps[200]
hi = 2 * D + 200 - ps[199]
if not (D, ps[199], ps[200], lo, hi) == (201254860, 26735756, 26626878, 375664886, 375774164):
    raise AssertionError()
print('root integer audit', D, ps[199], ps[200], F(lo, 2 * Q), F(hi, 2 * Q))
