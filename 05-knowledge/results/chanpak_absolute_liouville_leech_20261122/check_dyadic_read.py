"""Replay the new exact polynomial binary-read algorithm."""
from fractions import Fraction as Q
from functools import cache
from random import Random


@cache
def read(n):
    if not n:
        return 0
    q, b = divmod(n, 2)
    return read(q) + (0 if b else q + (q & 1))


def reverse(n, k):
    return sum(((n >> i) & 1) << (k - 1 - i) for i in range(k))


def rad(n):
    k = n.bit_length()
    return Q(reverse(n, k), 2**k)


def add(*polys):
    out = [sum(a[j] for a in polys if j < len(a))
           for j in range(max(map(len, polys)))]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def multiple(poly, c):
    return [c * a for a in poly]


def binary_read_poly(Y):
    if len(Y) == 1:
        return [read(Y[0])]
    w, *U = Y
    if w >= 0:
        return add(binary_read_poly(U), [0] + multiple(U, 1 - rad(w)),
                   multiple(U, -1), [read(w) + (U[0] & 1)])
    B = 2**(-w).bit_length()
    v = B + w
    R = reverse(-w-1, B.bit_length() - 1)
    return add(binary_read_poly(add(U, [-1])),
               [0] + multiple(U, rad(-w-1)),
               [read(v) - R + (v < B//2)])


def eval_poly(co, x):
    return sum(a*x**j for j, a in enumerate(co))


def check():
    rng = Random(31523)
    cases = [[w, u] for w in range(-24, 25) for u in (1, 2, 9)]
    cases += [[rng.randrange(-200, 201) for _ in range(d)] + [rng.randrange(1, 40)]
              for d in range(2, 7) for _ in range(160)]
    cases += [[-1, -1, 1], [0, -1, 1], [-1, 0, 1], [0, 0, 1]]
    for Y in cases:
        P = binary_read_poly(Y)
        k = 1-rad(Y[0]) if Y[0] >= 0 else rad(-Y[0]-1)
        degree = len(Y)-1
        assert (P[degree] if degree < len(P) else 0) == k*Y[-1]
        for m in (12, 19, 31, 73):
            q = 1 << m
            val = eval_poly(Y, q)
            assert val > 0 and eval_poly(P, q) == read(val), (Y, P, m)
    print(f"Exact new binary-read polynomial: {len(cases)} templates, four widths each.")


if __name__ == "__main__":
    check()
