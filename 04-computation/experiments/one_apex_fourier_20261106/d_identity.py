# S83 (D): exact checks for the two-block decomposition and the X^(2/3) counting family.
# g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2).
import sys, random, math
from functools import lru_cache
sys.setrecursionlimit(10000)
def g(v):
    r = 1; extra = 0
    # iterative: process bits from LSB
    acc = 0
    while v > 0:
        if v & 1:
            v >>= 1
        else:
            w = v >> 1
            acc += w + (w & 1)
            v = w
    return 1 + acc
# sanity: recursion
def g_rec(v):
    if v == 0: return 1
    if v & 1: return g_rec(v >> 1)
    w = v >> 1
    return g_rec(w) + w + (w & 1)
for v in range(1, 5000):
    assert g(v) == g_rec(v), v
def Cm(l, m):
    return sum(((1 - ((l >> i) & 1)) << (m - 1 - i)) for i in range(m))
# (*) g(h 2^m + l) = g(h) + g(l) - 1 + h*C_m(l) + (h&1)*[l < 2^(m-1)]
random.seed(1)
bad = 0
for _ in range(200000):
    m = random.randint(1, 20); h = random.randint(1, 1 << random.randint(1, 30)); l = random.randrange(1 << m)
    lhs = g((h << m) + l)
    rhs = g(h) + g(l) - 1 + h * Cm(l, m) + (h & 1) * (1 if l < (1 << (m - 1)) else 0)
    if lhs != rhs: bad += 1
print("(*) violations in 200000 random cases:", bad)
# C_m bijection, g(l) <= l+1
for m in range(1, 13):
    assert sorted(Cm(l, m) for l in range(1 << m)) == list(range(1 << m))
assert all(g(v) <= v + 1 for v in range(1 << 18))
print("C_m bijective for m<=12; g(v)<=v+1 for v<2^18")
