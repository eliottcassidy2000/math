#!/usr/bin/env python3
"""Auditor bridges-note, script A5 (independent code).

The note (section 10, and next obligation (5) of section 12) says: "the Cayley-Dickson sign tournaments are probably
the June doubling tower" and calls it "an observation for kind-pasteur's audit, not tested here".  The audited
THM-5133 (G) in the same worktree states: CD_3 = Paley P_7, CD_4 = T_15 of THM-448, CD_5 is NOT T_31 of THM-448.
This script tests that statement directly with cheap invariants:

  CD_n : sign tournament of the Cayley-Dickson algebra of dimension 2^n, doubling
         (a, b)(c, d) = (a c - conj(d) b, d a + b conj(c)), i -> j iff e_i e_j = + e_(i xor j)
  T_(2^k - 1) : core of the skew Hadamard tower of THM-448,
         S[i, j] = (-1)^(popcount((i & j) >> (v + 1)) + bit_v(i)), v = lsb(i xor j)

Invariants: doubly regular, H (orders 7, 15), multisets of H over the out-links and the in-links (order 31),
triple-intersection spectrum.
"""
import itertools
from collections import Counter
import numpy as np
from numba import njit


def conj(v):
    if len(v) == 1:
        return list(v)
    h = len(v) // 2
    return conj(v[:h]) + [-t for t in v[h:]]


def add(x, y):
    return [a + b for a, b in zip(x, y)]


def neg(x):
    return [-a for a in x]


def mul(x, y):
    n = len(x)
    if n == 1:
        return [x[0] * y[0]]
    h = n // 2
    a, b = x[:h], x[h:]
    c, d = y[:h], y[h:]
    first = add(mul(a, c), neg(mul(conj(d), b)))
    second = add(mul(d, a), mul(b, conj(c)))
    return first + second


def cd_tournament(n):
    N = 1 << n
    arcs = set()
    for i in range(1, N):
        ei = [0] * N
        ei[i] = 1
        for j in range(1, N):
            if i == j:
                continue
            ej = [0] * N
            ej[j] = 1
            p = mul(ei, ej)
            k = i ^ j
            assert sum(abs(t) for t in p) == 1 and abs(p[k]) == 1
            if p[k] == 1:
                arcs.add((i - 1, j - 1))
    # antisymmetry
    m = N - 1
    assert all(((i, j) in arcs) != ((j, i) in arcs) for i in range(m) for j in range(i + 1, m))
    return m, arcs


def tower_tournament(k):
    N = 1 << k
    arcs = set()
    for i in range(1, N):
        for j in range(1, N):
            if i == j:
                continue
            x = i ^ j
            v = (x & -x).bit_length() - 1
            e = bin((i & j) >> (v + 1)).count("1") + ((i >> v) & 1)
            if e % 2 == 0:
                arcs.add((i - 1, j - 1))
    m = N - 1
    assert all(((i, j) in arcs) != ((j, i) in arcs) for i in range(m) for j in range(i + 1, m))
    return m, arcs


@njit(cache=True)
def hcount(adj, n):
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, full + 1):
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            out = adj[v] & ~mask
            w = 0
            while out:
                if out & 1:
                    dp[mask | (1 << w), w] += c
                out >>= 1
                w += 1
    s = 0
    for v in range(n):
        s += dp[full, v]
    return s


def H_of(vs, arcs):
    idx = {v: t for t, v in enumerate(vs)}
    n = len(vs)
    adj = np.zeros(n, np.int64)
    for u in vs:
        for w in vs:
            if u != w and (u, w) in arcs:
                adj[idx[u]] |= (1 << idx[w])
    return int(hcount(adj, n))


def describe(name, m, arcs):
    out = {v: [w for w in range(m) if (v, w) in arcs] for v in range(m)}
    inn = {v: [w for w in range(m) if (w, v) in arcs] for v in range(m)}
    reg = all(len(out[v]) == (m - 1) // 2 for v in range(m))
    dr = all(len(set(out[u]) & set(out[v])) == (m - 3) // 4 for u in range(m) for v in range(u + 1, m))
    info = {"regular": reg, "doubly regular": dr}
    if m <= 15:
        info["H"] = H_of(list(range(m)), arcs)
    if m >= 15:
        info["out-link H"] = sorted(Counter(H_of(out[v], arcs) for v in range(m)).items())
        info["in-link H"] = sorted(Counter(H_of(inn[v], arcs) for v in range(m)).items())
    trip = Counter()
    for u, v, w in itertools.combinations(range(m), 3):
        trip[len(set(out[u]) & set(out[v]) & set(out[w]))] += 1
    info["triple out-intersections"] = sorted(trip.items())
    print(name)
    for k, v in info.items():
        print("    %s: %s" % (k, v))
    return info


def main():
    res = {}
    for n in (3, 4, 5):
        m, arcs = cd_tournament(n)
        res[("CD", n)] = describe("CD_%d (%d vertices)" % (n, m), m, arcs)
        m, arcs = tower_tournament(n)
        res[("T", n)] = describe("T_%d of the THM-448 tower" % m, m, arcs)

    def same(a, b):
        # isomorphism invariants up to converse
        direct = all(a[k] == b[k] for k in a if k in b)
        if direct:
            return True
        sw = dict(a)
        if "out-link H" in sw:
            sw["out-link H"], sw["in-link H"] = a["in-link H"], a["out-link H"]
        return all(sw[k] == b[k] for k in sw if k in b and k != "triple out-intersections")
    for n in (3, 4, 5):
        print("n = %d: invariants of CD_%d and T_%d agree (up to converse): %s" % (n, n, (1 << n) - 1, same(res[("CD", n)], res[("T", n)])))


if __name__ == "__main__":
    main()
