#!/usr/bin/env python3
"""Subset-DP helpers excerpted from Session243 constructions/verify.py at c6b359b9c0ca. Stdlib only."""
from itertools import permutations
from math import gcd, isqrt, lcm
import json
import sys


def need(test, message):
    if not test:
        raise RuntimeError(message)


def validate(a):
    n = len(a)
    need(all(0 <= x < 1 << n and not (x >> i & 1) for i, x in enumerate(a)), 'invalid row')
    need(all(((a[i] >> j) ^ (a[j] >> i)) & 1 for i in range(n) for j in range(i + 1, n)), 'not tournament')


def tables(a):
    validate(a)
    n = len(a)
    incoming = [sum(1 << u for u in range(n) if a[u] >> v & 1) for v in range(n)]
    end = [[0] * n for _ in range(1 << n)]
    h = [0] * (1 << n)
    h[0] = 1
    for s in range(1, 1 << n):
        vs = s
        while vs:
            bit = vs & -vs
            v = bit.bit_length() - 1
            t = s ^ bit
            if t == 0:
                end[s][v] = 1
            else:
                choices = t & incoming[v]
                total = 0
                while choices:
                    bitu = choices & -choices
                    total += end[t][bitu.bit_length() - 1]
                    choices ^= bitu
                end[s][v] = total
            vs ^= bit
        h[s] = sum(end[s])
    return h, end


def reverse(a):
    full = (1 << len(a)) - 1
    return [full ^ x ^ (1 << i) for i, x in enumerate(a)]


def pair(a):
    h, _ = tables(a)
    full = (1 << len(a)) - 1
    return h[full], sum(h[s] * h[full ^ s] for s in range(full + 1))


def state(a, root):
    h, end = tables(a)
    _, start = tables(reverse(a))
    full = (1 << len(a)) - 1
    return (h[full], end[full][root], sum(start[s][root] * h[full ^ s] for s in range(1, full + 1)), h[full ^ (1 << root)])


def ab(z):
    D, E, U, d = z
    return D + 2 * d - E - U, U - d


def tt(n):
    return [sum(1 << j for j in range(i + 1, n)) for i in range(n)]


def endpoint(k):
    a = tt(k + 2)
    a[0] ^= 1 << (k + 1)
    a[k + 1] |= 1
    return a


def join(a, b):
    n, m = len(a), len(b)
    return [x | (((1 << m) - 1) << n) for x in a] + [x << n for x in b]


def graft(a, root, b):
    n, m = len(a), len(b)
    new = [x if i == root else x | (((1 << m) - 1) << n) for i, x in enumerate(a)]
    return new + [(x << n) | (1 << root) for x in b]


def strong(a):
    full = (1 << len(a)) - 1
    for r in range(len(a)):
        seen = 1 << r
        while True:
            new = seen
            for u in range(len(a)):
                if seen >> u & 1:
                    new |= a[u]
            if new == seen:
                break
            seen = new
        if seen != full:
            return False
    return True


def labelled(n, code):
    a=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if code&1:a[i]|=1<<j
            else:a[j]|=1<<i
            code >>= 1
    return a
