"""Exact checks for the wrapped-transitive partition gadget and the 145 seed."""
from itertools import permutations
from fractions import Fraction
from math import prod

def paths(out):
    n = len(out)
    for p in permutations(range(n)):
        if all((out[x] >> y & 1 for x, y in zip(p, p[1:]))):
            yield p

def build_wrap(n):
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
    out = build_wrap(n)
    w = [Fraction((-1) ** i * (i + 1), 2 * i + 7) for i in range(n)]
    ps = list(paths(out))
    Z = sum((prod((w[i] for i in p[:p.index(n)] if i < n)) for p in ps), Fraction(0))
    rhs = prod((1 + x for x in w)) + prod(w)
    if not (len(ps) == 2 ** n + 1 and Z == rhs):
        raise AssertionError()
    print('wrap', n, 'count', len(ps), 'weighted=', Z)
O = [82, 40, 43, 65, 14, 89, 22]
if not all(((O[i] >> j & 1) + (O[j] >> i & 1) == 1 for i in range(7) for j in range(i + 1, 7))):
    raise AssertionError()
ps = list(paths(O))
if not len(ps) == 145:
    raise AssertionError()
cycle = (0, 1, 3, 6, 4, 2, 5, 0)
if not all((O[i] >> j & 1 for i, j in zip(cycle, cycle[1:]))):
    raise AssertionError()
print('145 seed endpoint counts:', [sum((p[-1] == i for p in ps)) for i in range(7)])
print('Hamilton cycle:', cycle)
