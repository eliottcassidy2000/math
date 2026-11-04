#!/usr/bin/env python3
"""auditb_common.py -- shared helpers for audit B of THM-4871 statement 7 (2026-11-04).

Independent of cpq_defect_census.py: own circulant classes, own brute-force canonical forms for the 4-coloured
objects on the k <= 4 fixed points (min over S_k, no nauty colour partitions), own object builder.
Sign convention: sgn_X(g) = (-1)^#{edges {u<v} of X with g(u) > g(v)}.
"""
import itertools, math
import pynauty


def units(m):
    return [a for a in range(1, m) if math.gcd(a, m) == 1]


def symmetric_sets(m):
    """all C with C = -C, 0 not in C (2^((m-1)/2) of them)."""
    half = list(range(1, (m + 1) // 2))
    out = []
    for bits in itertools.product([0, 1], repeat=len(half)):
        out.append(frozenset(x for d, b in zip(half, bits) if b for x in (d, m - d)))
    return out


def skew_sets(m):
    half = list(range(1, (m + 1) // 2))
    return [frozenset((d if b == 0 else m - d) for d, b in zip(half, bits))
            for bits in itertools.product([0, 1], repeat=len(half))]


def dilation_classes(sets, m):
    seen, reps = set(), []
    for C in sets:
        if C in seen:
            continue
        orb = {frozenset((a * c) % m for c in C) for a in units(m)}
        seen |= orb
        reps.append(C)
    return reps


def reversal_parity(g, edges):
    """parity of #{(u<v) in edges : g[u] > g[v]}."""
    return sum(1 for (u, v) in edges if g[u] > g[v]) & 1


def circulant_edges(C, m):
    return sorted({(min(i, (i + c) % m), max(i, (i + c) % m)) for i in range(m) for c in C})


def multiplier_stabiliser(C, m):
    return [a for a in units(m) if frozenset((a * c) % m for c in C) == C]


def multiplier_even(C, m):
    edges = circulant_edges(C, m)
    for a in multiplier_stabiliser(C, m):
        g = [(a * i) % m for i in range(m)]
        if reversal_parity(g, edges):
            return False
    return True


def perm_sign(g):
    """sign of a permutation given as a list."""
    n = len(g); seen = [False] * n; s = 1
    for i in range(n):
        if not seen[i]:
            j = i; L = 0
            while not seen[j]:
                seen[j] = True; j = g[j]; L += 1
            if L % 2 == 0:
                s = -s
    return s


# ---------------------------------------------------------------- 4-coloured objects on k points, brute force
def coloured_objects(k, directed, ncol=4):
    """4-coloured graphs (directed=False) or tournaments (directed=True) on k points up to colour-preserving
    isomorphism, by brute force over S_k.  Returns list of (edge_set, colour_tuple, even) where for graphs
    even = every colour-preserving automorphism has even reversal parity (colour-sign-even), and for
    tournaments even = True.  edge_set: frozenset of (u, v) with u < v for graphs, arcs (u, v) for tournaments."""
    if k == 0:
        return [(frozenset(), (), True)]
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    perms = list(itertools.permutations(range(k)))
    out = {}
    for bits in itertools.product([0, 1], repeat=len(pairs)):
        if directed:
            E = frozenset((u, v) if b else (v, u) for (u, v), b in zip(pairs, bits))
        else:
            E = frozenset(pr for pr, b in zip(pairs, bits) if b)
        for col in itertools.product(range(ncol), repeat=k):
            best = None; auts = []
            for s in perms:
                col2 = [None] * k
                for u in range(k):
                    col2[s[u]] = col[u]
                col2 = tuple(col2)
                if directed:
                    E2 = frozenset((s[u], s[v]) for (u, v) in E)
                else:
                    E2 = frozenset((min(s[u], s[v]), max(s[u], s[v])) for (u, v) in E)
                key = (col2, tuple(sorted(E2)))
                if best is None or key < best:
                    best = key
                if col2 == col and E2 == E:
                    auts.append(s)
            if best in out:
                continue
            if directed:
                even = True
            else:
                even = all(reversal_parity(s, sorted(E)) == 0 for s in auts)
            out[best] = (E, col, even)
    return list(out.values())


def burnside_b_k4(k):
    """number of isomorphism types of 4-coloured tournaments on k points = sum over odd cycle types."""
    from fractions import Fraction
    total = Fraction(0)
    for parts in partitions_odd(k):
        # z_lambda
        from collections import Counter
        cnt = Counter(parts); z = 1
        for part, mult in cnt.items():
            z *= part ** mult * math.factorial(mult)
        c1 = len(parts)
        c2 = sum((part - 1) // 2 for part in parts)  # pair orbits inside cycles
        c2 += sum(math.gcd(a, b) for i, a in enumerate(parts) for b in parts[i + 1:])
        total += Fraction(4 ** c1 * 2 ** c2, z)
    return total


def partitions_odd(k, maxpart=None):
    if maxpart is None:
        maxpart = k if k % 2 == 1 else k - 1
    if k == 0:
        yield []
        return
    for part in range(min(k, maxpart), 0, -1):
        if part % 2 == 1:
            for rest in partitions_odd(k - part, part):
                yield [part] + rest


# ---------------------------------------------------------------- the H-fixed objects
def build_graph(p, q, k, CP, CQ, cross, EF, col):
    """H-fixed graph on n = p+q+k vertices: P = 0..p-1, Q = p..p+q-1, F = p+q..n-1; returns sorted edge list."""
    n = p + q + k
    E = set()
    for i in range(p):
        for c in CP:
            j = (i + c) % p
            E.add((min(i, j), max(i, j)))
    for i in range(q):
        for c in CQ:
            j = (i + c) % q
            E.add((min(p + i, p + j), max(p + i, p + j)))
    if cross:
        for i in range(p):
            for j in range(q):
                E.add((i, p + j))
    for t in range(k):
        f = p + q + t
        x, y = col[t] & 1, col[t] >> 1
        if x:
            for i in range(p):
                E.add((i, f))
        if y:
            for j in range(q):
                E.add((p + j, f))
    for (u, v) in EF:
        E.add((p + q + u, p + q + v))
    return sorted(E)


def build_tournament(p, q, k, SP, SQ, cross, AF, col):
    """H-fixed tournament; returns list of arcs (u, v)."""
    n = p + q + k
    A = []
    for i in range(p):
        for s in SP:
            A.append((i, (i + s) % p))
    for i in range(q):
        for s in SQ:
            A.append((p + i, p + (i + s) % q))
    for i in range(p):
        for j in range(q):
            A.append((i, p + j) if cross == 0 else (p + j, i))
    for t in range(k):
        f = p + q + t
        x, y = col[t] & 1, col[t] >> 1
        for i in range(p):
            A.append((i, f) if x else (f, i))
        for j in range(q):
            A.append((p + j, f) if y else (f, p + j))
    for (u, v) in AF:
        A.append((p + q + u, p + q + v))
    return A


def graph_cert_and_gens(n, edges):
    adj = {v: [] for v in range(n)}
    for (u, v) in edges:
        adj[u].append(v); adj[v].append(u)
    g = pynauty.Graph(n, directed=False, adjacency_dict=adj)
    gens, grpsize1, grpsize2, _, _ = pynauty.autgrp(g)
    order = grpsize1 * 10 ** grpsize2        # nauty encodes |Aut| as mantissa * 10^exponent (float)
    return pynauty.certificate(g), gens, order


def tournament_cert(n, arcs):
    adj = {v: [] for v in range(n)}
    for (u, v) in arcs:
        adj[u].append(v)
    g = pynauty.Graph(n, directed=True, adjacency_dict=adj)
    return pynauty.certificate(g)


def sign_even(gens, edges):
    """sgn_X is a homomorphism, so X is sign-even iff every generator of Aut X is even."""
    return all(reversal_parity(g, edges) == 0 for g in gens)
