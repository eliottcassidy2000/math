"""Exact integer replay for the simultaneous residue/height selector."""
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt, prod
from pathlib import Path
from runpy import run_path

bank = run_path(str(Path(__file__).with_name("third-slot-bank.py")))["data"]
selected = [bank[i] for i in (0, 2, 6, 11, 13)]
lo, hi = zip(*[(u[0], v[0]) for _, _, u, v in selected])
ratios = tuple(F(v, u) for u, v in zip(lo, hi))
moves = ((0, 0, 2, -1, 0), (1, 0, -1, 1, -1), (3, 0, -1, -10, 9),
         (5, 0, -5, -1, 1), (8, 1, -14, -6, -8))
budget = (9, 7, 11, 5, 11)
t = 240
base = (0, 0, 211, 196, 103)
lower, upper = F(10065, 10000), F(10070, 10000)
eps = F(164, 10**9)

def P(a):
    return prod(r**n for r, n in zip(ratios, a))

def mod(r, M):
    return r.numerator * pow(r.denominator, -1, M) % M

def reading(x, M=None):
    if M:
        return prod(pow(u, t - n, M) * pow(v, n, M)
                    for u, v, n in zip(lo, hi, x)) % M
    return prod(u**(t - n) * v**n for u, v, n in zip(lo, hi, x))

def corners(D, Q, R):
    return (11*D+42*Q+20*R, 71*D+258*Q+116*R,
            23*D+90*Q+44*R, 77*D+294*Q+140*R)

G = tuple(map(P, moves))
assert tuple(map(lambda r: mod(r, 77), ratios)) == (41, 37, 47, 53, 2)
assert all(g > 1 and mod(g, 77) == 1 for g in G)
assert all(g**(b + 1) > prev for prev, g, b in zip((upper,) + G, G, budget))
assert G[-1] < 1 + eps
budget_both = budget[:-1] + (budget[-1] + 1,)
for a, col in zip(base, zip(*moves)):
    assert a + sum(min(0, v)*b for v, b in zip(col, budget_both)) >= 0
    assert a + 4 + sum(max(0, v)*b for v, b in zip(col, budget_both)) <= t

units = {s for s in range(77) if gcd(s, 77) == 1}
anchors = {}
Rmod = reading(base, 77)
for a in product(range(5), repeat=5):
    p = P(a)
    if 1 <= p <= lower:
        anchors.setdefault(Rmod * mod(p, 77) % 77, (a, p))
assert set(anchors) == units and len(units) == 60
Dmod = prod(pow(d, t, 2310) for d, _, _, _ in selected) % 2310
Qmod = prod(pow(q, t, 2310) for _, q, _, _ in selected) % 2310
assert (Dmod, Qmod) == (385, 1156)
assert all(u % 2 == v % 2 == 1 and u % 3 == v % 3 == 0
           and u % 5 == v % 5 != 0 for u, v in zip(lo, hi))
assert reading(base, 30) == 21
good = {s for s in units if all(gcd(c, 2310) == 1 for c in
        corners(Dmod, Qmod, next(r for r in range(s, 2310, 77) if r % 30 == 21))[:3])}
assert good == {s for s in units if s % 7 not in (2, 4) and s % 11 not in (1, 10)}
assert len(good) == 32

for y in (lower, F(10066, 10000), F(10068, 10000), upper):
    for s, (a, p) in anchors.items():
        digits = []
        for g, b in zip(G, budget):
            j = 0
            while p*g <= y:
                p *= g
                j += 1
            assert j <= b
            digits.append(j)
        x = tuple(u + v + sum(j*w for j, w in zip(digits, col))
                  for u, v, col in zip(base, a, zip(*moves)))
        xx = tuple(u + v for u, v in zip(x, moves[-1]))
        assert all(0 <= u <= t for u in (*x, *xx))
        assert reading(x, 77) == reading(xx, 77) == s
        assert P(tuple(u-v for u, v in zip(x, base))) == p
        assert P(tuple(u-v for u, v in zip(xx, base))) == p*G[-1]
        assert p <= y < p*G[-1] and p*(G[-1]-1) < eps*y

readings = [m[0] for _, _, u, v in bank for m in (u, v)]
markers = (79, 7499, 3469, 757, 2957, 587, 8117, 29789, 85597, 28559, 1223,
           20807, 99623, 557, 523, 307, 92269, 331, 479, None, 6113, 3671, 457,
           7283, 67, 43, 35267, 937)
for i, p in enumerate(markers):
    if p is None:
        assert readings[i] == 66825 and readings[i] > 1
    else:
        assert all(p % d != 0 for d in range(2, isqrt(p)+1))
        assert [j for j, r in enumerate(readings) if r % p == 0] == [i]
assert all(d <= q and q**7 < r**6
           for d, q, u, v in bank for r in (u[0], v[0]))
print("Verified both-sided 77-selector (60 phases, 32 guarded), and rank-28 capacity certificate.")
