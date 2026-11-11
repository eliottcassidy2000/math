#!/usr/bin/env python3
"""Input for halfpaths.c: the Paley tournament on GF(q), q = 3 (mod 4) prime or q = 27, with psi(x) = -x, v0 = 0.
Vertex 2i is a representative a_i of the i-th pair, vertex 2i+1 is -a_i.  Usage: paley_input.py q"""
import sys
q = int(sys.argv[1])
if q == 27:
    # GF(27) = F_3[x]/(x^3 - x - 1) (x^3 = x + 1); elements as integers c0 + 3 c1 + 9 c2
    def dig(a): return [a % 3, (a // 3) % 3, (a // 9) % 3]
    def und(d): return d[0] + 3 * d[1] + 9 * d[2]
    def add(a, b): return und([(x + y) % 3 for x, y in zip(dig(a), dig(b))])
    def neg(a): return und([(-x) % 3 for x in dig(a)])
    def mul(a, b):
        A, B = dig(a), dig(b); c = [0] * 5
        for i in range(3):
            for j in range(3): c[i + j] = (c[i + j] + A[i] * B[j]) % 3
        # x^4 = x^2 + x, x^3 = x + 1
        c[2] = (c[2] + c[4]) % 3; c[1] = (c[1] + c[4]) % 3
        c[1] = (c[1] + c[3]) % 3; c[0] = (c[0] + c[3]) % 3
        return und(c[:3])
    elems = list(range(27))
else:
    add = lambda a, b: (a + b) % q
    neg = lambda a: (-a) % q
    mul = lambda a, b: (a * b) % q
    elems = list(range(q))
sq = {mul(a, a) for a in elems if a != 0}
assert len(sq) == (q - 1) // 2 and all(neg(a) not in sq for a in sq)
arc = lambda a, b: add(b, neg(a)) in sq          # a -> b iff b - a is a non-zero square
reps = sorted(sq)                                 # one representative per pair: the squares
verts = []
for a in reps: verts += [a, neg(a)]
idx = {v: i for i, v in enumerate(verts)}
print(len(reps))
for v in verts:
    mask = 0
    for w in verts:
        if w != v and arc(v, w): mask |= 1 << idx[w]
    print(mask, 1 if arc(0, v) else 0)
