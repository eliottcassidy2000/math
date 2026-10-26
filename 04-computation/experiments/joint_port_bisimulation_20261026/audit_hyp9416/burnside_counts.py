# Independent audit of HYP-9416, check A: numbers of rooted two-port states [C; s, t] (s != t) and of strong ones,
# WITHOUT enumerating any tournament.
#   (1) Burnside over S_n: a permutation of cycle type lam fixes a tournament iff every cycle is odd, and then fixes
#       2^(sum (l_i - 1)/2 + sum_{i<j} gcd(l_i, l_j)) tournaments.  Orbits of (T), (T, v), (T, s, t) are weighted by
#       1, m_1, m_1 (m_1 - 1) (m_1 = number of fixed points).
#   (2) Strong-component calculus: every tournament is uniquely an ordered sequence of strong components, and its
#       automorphism group is the product of the components' groups.  With A, B, Cc the generating functions of strong
#       classes, vertex-rooted strong classes and ordered-pair-rooted strong classes, and T, V, R those of all classes:
#           T = 1/(1-A),   V = T^2 B,   R = T^2 Cc + 2 T^3 B^2.
#       Hence A, B, Cc (in particular the strong two-port counts) follow from the Burnside numbers alone.
# Audit of mac-mini-2026-10-25-pyth's HYP-9416; written independently (imports nothing from the lane).
from fractions import Fraction
from math import gcd, factorial
from collections import Counter

NMAX = 9

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield (); return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest

def z(lam):
    c = Counter(lam); r = 1
    for part, mult in c.items(): r *= part ** mult * factorial(mult)
    return r

def fixT(lam):
    if any(l % 2 == 0 for l in lam): return 0
    e = sum((l - 1) // 2 for l in lam) + sum(gcd(lam[i], lam[j]) for i in range(len(lam)) for j in range(i + 1, len(lam)))
    return 2 ** e

t = [Fraction(1)] + [Fraction(0)] * NMAX; v = [Fraction(0)] * (NMAX + 1); r = [Fraction(0)] * (NMAX + 1)
for n in range(1, NMAX + 1):
    for lam in partitions(n):
        f = fixT(lam)
        if not f: continue
        m1 = lam.count(1)
        t[n] += Fraction(f, z(lam)); v[n] += Fraction(f * m1, z(lam)); r[n] += Fraction(f * m1 * (m1 - 1), z(lam))
assert all(x.denominator == 1 for x in t + v + r)
t = [int(x) for x in t]; v = [int(x) for x in v]; r = [int(x) for x in r]

def mul(a, b):
    c = [0] * (NMAX + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(NMAX + 1 - i): c[i + j] += x * b[j]
    return c
def inv(a):   # a[0] = 1
    b = [0] * (NMAX + 1); b[0] = 1
    for n in range(1, NMAX + 1): b[n] = -sum(a[k] * b[n - k] for k in range(1, n + 1))
    return b

Tinv = inv(t)
A = [0] + [-x for x in Tinv[1:]]              # A = 1 - 1/T
T2inv = mul(Tinv, Tinv)
B = mul(v, T2inv)                             # B = V / T^2
T3 = mul(mul(t, t), t)
Cc = mul([ri - 2 * x for ri, x in zip(r, mul(T3, mul(B, B)))], T2inv)

print("n | tournament classes | vertex-rooted | two-port states (all) | strong classes | strong vertex-rooted | strong two-port states")
for n in range(1, NMAX + 1):
    print(f"{n} | {t[n]} | {v[n]} | {r[n]} | {A[n]} | {B[n]} | {Cc[n]}")
claimed_all = {3: 8, 4: 32, 5: 176, 6: 1408, 7: 17824}
claimed_strong = {3: 2, 4: 12, 5: 92, 6: 946, 7: 14162}
ok = all(r[n] == claimed_all[n] and Cc[n] == claimed_strong[n] for n in claimed_all)
print("HYP-9416 counts 8/2, 32/12, 176/92, 1408/946, 17824/14162 reproduced:", ok)
assert [t[n] for n in range(1, 9)] == [1, 1, 2, 4, 12, 56, 456, 6880]       # OEIS A000568
assert [A[n] for n in range(1, 9)] == [1, 0, 1, 1, 6, 35, 353, 6008]        # OEIS A051337
