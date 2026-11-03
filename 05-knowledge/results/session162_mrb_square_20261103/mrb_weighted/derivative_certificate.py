from itertools import permutations
from fractions import Fraction
from math import prod

def paths(out):
    n = len(out)
    return [p for p in permutations(range(n)) if all((out[x] >> y & 1 for x, y in zip(p, p[1:])))]

def wrap(n):
    B, C = (n, n + 1)
    out = [0] * (n + 2)
    for i in range(n):
        for j in range(i + 1, n):
            out[i] |= 1 << j
        out[i] |= 1 << B
        out[C] |= 1 << i
    out[B] |= 1 << C
    return out
for n in range(1, 7):
    ps = paths(wrap(n))
    v = [Fraction((-1) ** i * (i + 1), 2 * i + 7) for i in range(n)]
    counts = [sum((p.index(i) < p.index(n) for p in ps)) for i in range(n)]
    if not counts == [2 ** (n - 1) + 1] * n:
        raise AssertionError()
    incidence = sum((sum((v[i] for i in range(n) if p.index(i) < p.index(n)), Fraction()) for p in ps), Fraction())
    if not incidence == (2 ** (n - 1) + 1) * sum(v):
        raise AssertionError()
    if n >= 2:
        subpaths = paths(wrap(n - 1))
        for i in range(n):
            projected = [tuple((x - (x > i) for x in p if x != i)) for p in ps if p.index(i) < p.index(n)]
            if not (len(projected) == len(set(projected)) and set(projected) == set(subpaths)):
                raise AssertionError()
    print('n', n, 'P', len(ps), 'mark', counts[0], 'directional derivative', incidence)
