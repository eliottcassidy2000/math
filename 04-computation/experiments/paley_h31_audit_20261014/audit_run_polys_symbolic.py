#!/usr/bin/env python3
"""AUDIT: exact (symbolic in q) injective pattern sums for unions of paths, valid for every q = 3 mod 4.

Sigma*(G) = sum over set partitions pi of V(G) of mu(0,pi) Sigma(G/pi),   mu = prod_B (-1)^(|B|-1) (|B|-1)!
Sigma(H) (unrestricted) is evaluated symbolically by the class-F calculus of THM-4652 part 5 extended by reduction
rules: loops (T(0)), isolated vertices (q), pendant vertices (sum_d T(d) = gamma q + beta), series reduction of degree-2
vertices (convolution), parallel merging (pointwise product), and, at an irreducible core, expansion of one label
alpha chi + gamma + beta delta into its three parts (delete / contract).  A core whose labels are all pure chi is
evaluated as 0 if it is a complete graph (THM-4652 part 4, unrestricted = injective there); any other pure-chi core is
reported (it would carry non-polynomial arithmetic, e.g. W_4).
Outputs exact polynomials for A_2, A_4, A_6, A_8 (single runs) and the joint A_(2,2), A_(2,4); validated against brute
force.  Usage: python audit_run_polys_symbolic.py
"""
import sys, itertools, time
from fractions import Fraction
from math import factorial
import numpy as np
import numba as nb

sys.setrecursionlimit(100000)


# ---------------------------------------------------------------- polynomials in q: dict deg -> int
def padd(a, b):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + v
    return {k: v for k, v in r.items() if v}


def pmul(a, b):
    r = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            r[k1 + k2] = r.get(k1 + k2, 0) + v1 * v2
    return {k: v for k, v in r.items() if v}


def pscale(a, c):
    return {k: c * v for k, v in a.items() if c * v}


ONE = {0: 1}
Q = {1: 1}
ZERO = {}


def pstr(a):
    if not a:
        return '0'
    return ' + '.join(f'{a[k]}*q^{k}' for k in sorted(a, reverse=True))


def peval(a, q):
    return sum(v * q ** k for k, v in a.items())


# label = (alpha, gamma, beta): T(d) = alpha chi(d) + gamma + beta [d = 0]
CHI = (ONE, ZERO, ZERO)


def rev(T):
    a, g, b = T
    return (pscale(a, -1), g, b)


def conv(F, G):
    a1, g1, b1 = F
    a2, g2, b2 = G
    alpha = padd(pmul(a1, b2), pmul(b1, a2))
    gamma = padd(padd(pmul(a1, a2), pmul(pmul(g1, g2), Q)), padd(pmul(g1, b2), pmul(b1, g2)))
    beta = padd(pscale(pmul(pmul(a1, a2), Q), -1), pmul(b1, b2))
    return (alpha, gamma, beta)


def pprod(F, G):
    a1, g1, b1 = F
    a2, g2, b2 = G
    alpha = padd(pmul(a1, g2), pmul(g1, a2))
    gamma = padd(pmul(a1, a2), pmul(g1, g2))
    beta = padd(padd(pscale(pmul(a1, a2), -1), pmul(g1, b2)), padd(pmul(b1, g2), pmul(b1, b2)))
    return (alpha, gamma, beta)


def is_zero_label(T):
    return not T[0] and not T[1] and not T[2]


CORES = {}


def unrestricted(verts, edges):
    """verts: set; edges: list of (u, v, T) meaning T(x_v - x_u). Returns polynomial Sigma."""
    verts = set(verts)
    edges = [e for e in edges]
    factor = ONE
    while True:
        changed = False
        # loops
        ne = []
        for (u, v, T) in edges:
            if u == v:
                factor = pmul(factor, padd(T[1], T[2]))       # T(0) = gamma + beta
                changed = True
            else:
                ne.append((u, v, T))
        edges = ne
        if not factor:
            return ZERO
        # parallel merge
        groups = {}
        for (u, v, T) in edges:
            key = (min(u, v), max(u, v))
            TT = T if u == key[0] else rev(T)
            if key in groups:
                groups[key] = pprod(groups[key], TT)
                changed = True
            else:
                groups[key] = TT
        edges = [(k[0], k[1], T) for k, T in groups.items() if not is_zero_label(T)]
        if any(is_zero_label(T) for T in groups.values()):
            # a zero label kills everything
            return ZERO
        deg = {v: 0 for v in verts}
        for (u, v, T) in edges:
            deg[u] += 1
            deg[v] += 1
        # isolated
        for v in list(verts):
            if deg[v] == 0:
                factor = pmul(factor, Q)
                verts.discard(v)
                changed = True
        if changed:
            continue
        # pendant
        done = False
        for v in list(verts):
            if deg[v] == 1:
                for i, (a, b, T) in enumerate(edges):
                    if a == v or b == v:
                        factor = pmul(factor, padd(pmul(T[1], Q), T[2]))     # sum_d T(d) = gamma q + beta
                        edges.pop(i)
                        break
                verts.discard(v)
                done = True
                break
        if done:
            if not factor:
                return ZERO
            continue
        # series
        for w in list(verts):
            if deg[w] == 2:
                inc = [(i, e) for i, e in enumerate(edges) if e[0] == w or e[1] == w]
                (i1, e1), (i2, e2) = inc
                # orient e1 as (u -> w), e2 as (w -> v)
                u = e1[0] if e1[1] == w else e1[1]
                T1 = e1[2] if e1[1] == w else rev(e1[2])
                v = e2[1] if e2[0] == w else e2[0]
                T2 = e2[2] if e2[0] == w else rev(e2[2])
                newT = conv(T1, T2)
                edges = [e for j, e in enumerate(edges) if j not in (i1, i2)]
                edges.append((u, v, newT))
                verts.discard(w)
                done = True
                break
        if done:
            continue
        break
    if not edges:
        return pmul(factor, {0: 1} if not verts else {len(verts): 1})
    # irreducible core: expand a non-pure label if any
    for i, (u, v, T) in enumerate(edges):
        a, g, b = T
        if g or b:
            rest = edges[:i] + edges[i + 1:]
            tot = ZERO
            if a:
                tot = padd(tot, unrestricted(verts, rest + [(u, v, (a, ZERO, ZERO))]))
            if g:
                tot = padd(tot, pmul(g, unrestricted(verts, rest)))
            if b:
                # contract v into u
                rest2 = []
                for (x, y, S) in rest:
                    x2 = u if x == v else x
                    y2 = u if y == v else y
                    rest2.append((x2, y2, S))
                tot = padd(tot, pmul(b, unrestricted(verts - {v}, rest2)))
            return pmul(factor, tot)
    # pure chi core
    nv = len(verts)
    pairs = {(min(u, v), max(u, v)) for (u, v, T) in edges}
    key = (nv, len(edges))
    CORES[key] = CORES.get(key, 0) + 1
    if len(pairs) == nv * (nv - 1) // 2 and len(edges) == len(pairs):
        return ZERO                     # complete pattern: chi(Vandermonde) symmetry (THM-4652 part 4)
    raise RuntimeError(f'non-complete pure-chi core: {nv} vertices, edges {[(u, v) for (u, v, T) in edges]}')


def partitions_no_adjacent(n, forbid):
    """set partitions of range(n) (restricted growth), with forbid[i] = list of j < i that may not share i's block"""
    blocks = []

    def rec(i):
        if i == n:
            yield [list(b) for b in blocks]
            return
        for b in blocks:
            if not any(j in b for j in forbid[i]):
                b.append(i)
                yield from rec(i + 1)
                b.pop()
        blocks.append([i])
        yield from rec(i + 1)
        blocks.pop()
    yield from rec(0)


def injective_sum_paths(lengths):
    """Sigma* of disjoint directed paths with the given numbers of edges, as a polynomial in q (q = 3 mod 4)"""
    verts = []
    edges = []
    forbid = {}
    off = 0
    for L in lengths:
        for i in range(L + 1):
            forbid[off + i] = [off + i - 1] if i > 0 else []
        for i in range(L):
            edges.append((off + i, off + i + 1))
        off += L + 1
    n = off
    total = ZERO
    cnt = 0
    for part in partitions_no_adjacent(n, forbid):
        cnt += 1
        blk = {}
        mu = 1
        for bi, b in enumerate(part):
            for x in b:
                blk[x] = bi
            mu *= (-1) ** (len(b) - 1) * factorial(len(b) - 1)
        qverts = set(range(len(part)))
        qedges = [(blk[a], blk[b], CHI) for a, b in edges]
        s = unrestricted(qverts, qedges)
        if s:
            total = padd(total, pscale(s, mu))
    return total, cnt


# ---------------------------------------------------------------- brute force (validation)
@nb.njit(cache=True)
def brute_paths(q, chi, lens):
    # injective sum over disjoint paths with lengths lens (<= 2 paths, total points <= 7)
    npts = 0
    for L in lens:
        npts += L + 1
    v = np.zeros(npts, dtype=np.int64)
    total = 1
    for _ in range(npts):
        total *= q
    tot = 0
    for code in range(total):
        c = code
        for i in range(npts):
            v[i] = c % q
            c //= q
        ok = True
        for i in range(npts):
            for j in range(i + 1, npts):
                if v[i] == v[j]:
                    ok = False
                    break
            if not ok:
                break
        if not ok:
            continue
        p = 1
        off = 0
        for L in lens:
            for i in range(L):
                p *= chi[(v[off + i + 1] - v[off + i]) % q]
            off += L + 1
        tot += p
    return tot


def chi_table(q):
    chi = np.full(q, -1, dtype=np.int64)
    chi[0] = 0
    for x in range(1, q):
        chi[(x * x) % q] = 1
    return chi


def factor_str(poly):
    import sympy
    qq = sympy.symbols('q')
    expr = sum(sympy.Integer(v) * qq ** k for k, v in poly.items())
    return str(sympy.factor(expr))


def main():
    results = {}
    for lens in ([2], [3], [4], [5], [6], [2, 2], [2, 4], [8]):
        t0 = time.time()
        poly, cnt = injective_sum_paths(lens)
        results[tuple(lens)] = poly
        print(f'Sigma*(paths {lens}) = {factor_str(poly)}    [{cnt} partitions, {time.time() - t0:.1f} s]', flush=True)
    print('cores met (vertices, edges): count ->', CORES)
    # validation against brute force at small q
    for lens in ([2], [4], [6], [2, 2]):
        for q in (7, 11):
            if sum(L + 1 for L in lens) > q:
                continue
            b = int(brute_paths(q, chi_table(q), np.array(lens, dtype=np.int64)))
            s = peval(results[tuple(lens)], q)
            print(f'   brute force {lens} at q = {q}: {b}; symbolic {s}; agree {b == s}')
    # comparisons with the claims
    import sympy
    qq = sympy.symbols('q')
    claims = {(2,): qq * (qq - 1), (4,): 2 * qq * (qq - 1) * (qq - 3), (6,): qq * (qq - 1) * (qq - 3) * (5 * qq - 29),
              (2, 2): qq * (qq - 1) * (qq - 3) * (qq - 9)}
    for k, c in claims.items():
        expr = sum(sympy.Integer(v) * qq ** d for d, v in results[k].items())
        print(f'claim for {list(k)}: {sympy.factor(c)} -> identical: {sympy.expand(expr - c) == 0}')
    # A_8 against the June ADDENDUM-8 values
    a8 = results[(8,)]
    print('A_8 at q = 19, 23, 31:', [peval(a8, q) for q in (19, 23, 31)], '(THM-438 ADDENDUM-8: 10243584, 34327040, 202487040)')


if __name__ == '__main__':
    main()
