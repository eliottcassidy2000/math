#!/usr/bin/env python3
"""audit G1 (independent): two-step reciprocals of THM-4700(2), the set S, and the Pythagorean-tree remark.
Pure integer arithmetic; no session code is imported."""
from math import isqrt
from sympy import factorint, primerange

def sqfree(n):
    out = 1
    for p, e in factorint(n).items():
        if e % 2: out *= p
    return out

# --- S and Gamma cap Q
S = [p for p in primerange(2, 1100) if p % 15 == 14 and pow(p, 5, 11) == 10]
print('S up to 1100:', S)
print('primes < 29 in S:', [p for p in primerange(2, 29) if p in S], '; 29 in S:', 29 in S)

# --- two-step reciprocals: L(1/n) = 2 iff sqrt(-(4n^2 - 1)) lies in F iff sqfree(4n^2 - 1) in {3, 11, 15, 55}
Dset = {3, 11, 15, 55}
hits = [(n, sqfree(4 * n * n - 1)) for n in range(2, 12001) if sqfree(4 * n * n - 1) in Dset]
print('n <= 12000 with L(1/n) = 2 (n, D):', hits)
# Pell generation independently
def pell_even_x(D, xmax):
    # fundamental solution by brute force, then powers; keep even x
    y = 1
    while True:
        x2 = 1 + D * y * y; x = isqrt(x2)
        if x * x == x2: break
        y += 1
    x1, y1 = x, y; out = []; x, y = x1, y1
    while x <= xmax:
        if x % 2 == 0: out.append((x // 2, y))
        x, y = x1 * x + D * y1 * y, x1 * y + y1 * x
    return (x1, y1), out
for D in sorted(Dset):
    fund, sols = pell_even_x(D, 24000)
    print('D = %2d: fundamental solution %s ; n = x/2 <= 12000 with x even: %s' % (D, fund, [n for n, y in sols if n >= 2]))
# parity of x for D = 55 along the first 30 powers
x1, y1 = 89, 12; x, y = x1, y1; par = []
for _ in range(30):
    par.append(x % 2); x, y = x1 * x + 55 * y1 * y, x1 * y + y1 * x
print('D = 55: x mod 2 along the first 30 solutions:', set(par))
# the walk for each hit: v = (1 + y sqrt-D)/(2n) has norm (1 + D y^2)/(4n^2) = 1
for n, D in hits:
    y2 = (4 * n * n - 1) // D; y = isqrt(y2); assert y * y == y2 and 1 + D * y * y == 4 * n * n
print('every hit has an integral y with 1 + D y^2 = 4 n^2: True')

# --- Pythagorean tree (Berggren/Barning/Hall matrices A, B, C)
A = ((1, -2, 2), (2, -1, 2), (2, -2, 3)); B = ((1, 2, 2), (2, 1, 2), (2, 2, 3)); Cm = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
def ap(M, v): return tuple(sum(M[i][k] * v[k] for k in range(3)) for i in range(3))
tree = {(3, 4, 5): ''}
frontier = [((3, 4, 5), '')]
for depth in range(1, 10):
    nxt = []
    for v, word in frontier:
        for name, M in (('A', A), ('B', B), ('C', Cm)):
            w = ap(M, v); tree[w] = word + name; nxt.append((w, word + name))
    frontier = nxt
print('tree nodes to depth 9:', len(tree))
def near60(t):
    a, b, c = t
    return a == (c + 1) // 2 and (c + 1) % 2 == 0 or a == (c - 1) // 2 and (c - 1) % 2 == 0 or b == (c + 1) // 2 and (c + 1) % 2 == 0 or b == (c - 1) // 2 and (c - 1) % 2 == 0
n60 = sorted(((t, tree[t]) for t in tree if near60(t)), key=lambda z: z[0][2])
print('near-60-degree triples (a leg (c +- 1)/2) to depth 9:', n60)
print('their words are prefixes of CACACA...:', all('CACACACACA'.startswith(wd) for _, wd in n60))
# Pell check X = (3c -+ 1)/2, Y = other leg
for (a, b, c), wd in n60:
    leg = a if a in ((c + 1) // 2, (c - 1) // 2) and abs(2 * a - c) == 1 else b
    other = b if leg == a else a
    X = (3 * c - 1) // 2 if 2 * leg == c + 1 else (3 * c + 1) // 2
    assert X * X - 3 * other * other == 1, (a, b, c)
print('X^2 - 3 Y^2 = 1 with X = (3c -+ 1)/2 (sign as leg = (c +- 1)/2) and Y = other leg: True for all')
# the converse: Pell solutions (2 + sqrt3)^k, k >= 2, give primitive near-60 triples on the CA branch
x, y = 2, 1; got = []
for k in range(2, 12):
    x, y = 2 * x + 3 * y, x + 2 * y
    for s in (1, -1):
        if (2 * x + s) % 3 == 0:
            c = (2 * x + s) // 3
            leg = (c + 1) // 2 if s == 1 else (c - 1) // 2
            if c % 2 == 1 and leg * leg + y * y == c * c: got.append((k, c, leg, y))
print('Pell solutions (2 + sqrt3)^k -> triples (k, c, leg (c +- 1)/2, other leg):', got)
# 'symmetric three-step walk' 1/c = delta + conj(delta) -+ 1 with delta = (a + b i)/c
for (a, b, c), wd in n60:
    for leg in (a, b):
        if 2 * leg == c + 1: assert 2 * leg - c == 1           # 1/c = delta + conj delta - 1
        if 2 * leg == c - 1: assert -2 * leg + c == 1          # 1/c = -(delta + conj delta) + 1 with delta = (leg + i other)/c
print('three-step walks verified: 1/c = (2 leg)/c - 1 or 1 - (2 leg)/c')
# Eisenstein triangles behind the delta walks of THM-4700(2): the special angle is opposite the side named "opp"
for (p, q, opp, ang) in ((3, 5, 7, 120), (5, 16, 19, 120), (7, 15, 13, 60), (8, 15, 13, 60)):
    twocos = {120: -1, 60: 1}[ang]
    print('triangle with sides %d, %d and %d opposite a %d-degree angle: %s' % (p, q, opp, ang, opp * opp == p * p + q * q - twocos * p * q))
# the Eisenstein integers behind them: 13 v = 8 + 15 zeta3 for v = (1 + 15 sqrt-3)/26; g^2 for g = 2 + 3 zeta3, 3 + 5 zeta3
N = lambda a, b: a * a - a * b + b * b          # norm of a + b zeta3
print('N(8 + 15 zeta3) = %d = 13^2 ; N(2 + 3 zeta3) = %d ; N(3 + 5 zeta3) = %d' % (N(8, 15), N(2, 3), N(3, 5)))
