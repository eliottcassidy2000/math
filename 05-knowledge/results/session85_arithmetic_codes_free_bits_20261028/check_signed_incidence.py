#!/usr/bin/env python3
"""Exact signed-incidence, hypergraph and auxiliary-graph controls; stdlib only."""
from itertools import combinations, product
from math import gcd
from fractions import Fraction
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def determinant(a):
    a = [list(r) for r in a]
    n = len(a)
    if not n:
        return 1
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        p = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * p - a[i][k] * a[k][j]) // previous
        for i in range(k + 1, n):
            a[i][k] = 0
        previous = p
    return sign * a[-1][-1]


def smith_by_minors(a, n):
    """Nonzero Smith factors via gcds of every minor; no SNF oracle."""
    m, previous, factors = len(a), 1, []
    for k in range(1, min(m, n) + 1):
        delta = 0
        for rows in combinations(range(m), k):
            for cols in combinations(range(n), k):
                d = determinant([[a[i][j] for j in cols] for i in rows])
                delta = gcd(delta, abs(d))
                if delta == 1:
                    break
            if delta == 1:
                break
        if not delta:
            break
        require(delta % previous == 0, 'determinantal divisibility')
        factors.append(delta // previous)
        previous = delta
    require(all(y % x == 0 for x, y in zip(factors, factors[1:])), 'Smith divisibility')
    return factors


def rank_mod(a, n, p):
    a = [[x % p for x in r] for r in a]
    r = 0
    for j in range(n):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][j], -1, p)
        a[r] = [x * inv % p for x in a[r]]
        for i in range(len(a)):
            if i != r:
                c = a[i][j]
                a[i] = [(x - c * y) % p for x, y in zip(a[i], a[r])]
        r += 1
    return r


def signed_matrix(n, edges):
    # (u,v,sigma) means x_v - sigma*x_u = q_e.
    a = []
    for u, v, s in edges:
        row = [0] * n
        row[u] -= s
        row[v] += 1
        a.append(row)
    return a


def components_balance(n, edges):
    adj = [[] for _ in range(n)]
    for u, v, s in edges:
        adj[u].append((v, s))
        adj[v].append((u, s))
    signs, c, b = {}, 0, 0
    for root in range(n):
        if root in signs:
            continue
        c += 1
        balanced = True
        signs[root] = 1
        todo = [root]
        while todo:
            u = todo.pop()
            for v, s in adj[u]:
                if v not in signs:
                    signs[v] = s * signs[u]
                    todo.append(v)
                elif signs[v] != s * signs[u]:
                    balanced = False
        b += balanced
    return c, b


def graph_controls():
    count, profiles, switch_classes = 0, {}, {}
    for n in range(1, 5):
        pairs = list(combinations(range(n), 2))
        for state in product((0, -1, 1), repeat=len(pairs)):
            edges = [(u, v, s) for (u, v), s in zip(pairs, state) if s]
            a = signed_matrix(n, edges)
            c, b = components_balance(n, edges)
            expected = [1] * (n - c) + [2] * (c - b)
            require(smith_by_minors(a, n) == expected, ('signed Smith', n, edges))
            require(rank_mod(a, n, 2) == n - c, 'mod2 rank')
            require(rank_mod(a, n, 3) == n - b, 'mod3 rank')
            require(rank_mod(a, n, 5) == n - b, 'mod5 rank')
            # Canonical switching: make every spanning-forest edge positive.
            epsilon = {}
            for root in range(n):
                if root in epsilon:
                    continue
                epsilon[root] = 1
                todo = [root]
                while todo:
                    u = todo.pop()
                    for a0, v0, sign in edges:
                        v = v0 if a0 == u else a0 if v0 == u else None
                        if v is not None and v not in epsilon:
                            epsilon[v] = sign * epsilon[u]
                            todo.append(v)
            canonical = tuple(epsilon[u] * sign * epsilon[v] for u, v, sign in edges)
            underlying = (n, tuple((u, v) for u, v, sign in edges))
            switch_classes.setdefault(underlying, set()).add(canonical)
            require((all(z == 1 for z in canonical)) == (b == c), 'balance iff trivial holonomy')
            count += 1
            profiles[str((n, c, b))] = profiles.get(str((n, c, b)), 0) + 1
    explicit_orbit_count = 0
    for (n, pairs), classes in switch_classes.items():
        c, _ = components_balance(n, [(u, v, 1) for u, v in pairs])
        expected_classes = 2 ** (len(pairs) - n + c)
        require(len(classes) == expected_classes, 'switching quotient size')
        # Independent orbit construction: act by every vertex-sign assignment,
        # without using the spanning-forest normalization above.
        orbits = set()
        for signs in product((-1, 1), repeat=len(pairs)):
            orbit = frozenset(
                tuple(eta[u] * sign * eta[v] for (u, v), sign in zip(pairs, signs))
                for eta in product((-1, 1), repeat=n)
            )
            require(len(orbit) == 2 ** (n - c), 'explicit switching orbit size')
            orbits.add(orbit)
        require(len(orbits) == expected_classes, 'explicit switching orbit count')
        require(len(set().union(*orbits)) == 2 ** len(pairs), 'orbit coverage')
        require(sum(map(len, orbits)) == 2 ** len(pairs), 'orbit disjointness')
        require(all(len(orbit.intersection(classes)) == 1 for orbit in orbits),
                'one canonical forest representative per orbit')
        explicit_orbit_count += len(orbits)
    return {'all_signed_simple_graphs_through_n4': count, 'profiles': profiles,
            'underlying_graph_switching_quotients': len(switch_classes),
            'explicit_vertex_switching_orbits': explicit_orbit_count}


def matvec(a, x):
    return [sum(u * v for u, v in zip(row, x)) for row in a]


def petersen_controls():
    outer = [(i, (i + 1) % 5) for i in range(5)]
    inner = [(i + 5, (i + 2) % 5 + 5) for i in range(5)]
    spokes = [(i, i + 5) for i in range(5)]
    edges = outer + inner + spokes
    b = signed_matrix(10, [(u, v, -1) for u, v in edges])
    require(smith_by_minors(b, 10) == [1] * 9 + [2], 'Petersen Smith')
    # Exhaustive minors of the 15x10 matrix are unnecessary after first unit
    # factors: the above gcd procedure nevertheless checks them exactly.
    k5 = outer + [(u - 5, v - 5) for u, v in inner]
    bk = signed_matrix(5, [(u, v, -1) for u, v in k5])
    require(smith_by_minors(bk, 5) == [1] * 4 + [2], 'K5 Smith')
    for j in range(10):
        x = [int(i == j) for i in range(10)]
        q = matvec(b, x)
        adjusted = q[:5] + [q[10 + u] + q[10 + v] - q[5 + i]
                             for i, (u, v) in enumerate(k5[5:])]
        require(matvec(bk, x[:5]) == adjusted, 'spoke elimination matrix identity')
        require(x[5:] == [q[10 + i] - x[i] for i in range(5)], 'spoke label reconstruction')
    keep = [v for v in range(10) if v != 0]
    old = [(keep.index(u), keep.index(v), -1) for u, v in edges if u and v]
    require(smith_by_minors(signed_matrix(9, old), 9) == [1] * 8 + [2], 'Petersen minus0 Smith')
    remaining = [2, 3, 6, 7, 8, 9]
    ears = [(2, 1, 6), (3, 4, 9), (7, 5, 8)]
    quotient = [(u, v, -1) for u, v in edges if u in remaining and v in remaining]
    quotient += [(u, v, 1) for u, w, v in ears]
    marked = [(remaining.index(u), remaining.index(v), s) for u, v, s in quotient]
    bare = [(u, v, -1) for u, v, s in marked]
    require(len(quotient) == 9, 'K33 size')
    require(smith_by_minors(signed_matrix(6, marked), 6) == [1] * 5 + [2], 'marked K33 Smith')
    require(smith_by_minors(signed_matrix(6, bare), 6) == [1] * 5, 'bare K33 Smith')
    for j in range(10):
        x = [int(i == j) for i in range(10)]
        for u, w, v in ears:
            require(x[v] - x[u] == (x[w] + x[v]) - (x[u] + x[w]), 'ear signed matrix identity')
    return {'Petersen': [1] * 9 + [2], 'K5': [1] * 4 + [2],
            'Petersen_minus0': [1] * 8 + [2], 'marked_K33': [1] * 5 + [2],
            'bare_K33': [1] * 5, 'marked_edges': [list(t) for t in ears]}


def hypergraph_controls():
    simplex = {}
    for k in range(1, 7):
        a = [[int(i != j) for j in range(k + 1)] for i in range(k + 1)]
        expected = [1] * k + ([k] if k != 1 else [1])
        require(smith_by_minors(a, k + 1) == expected, 'J-I Smith')
        simplex[str(k)] = expected
    # Large k: exact determinant + explicit image/solution tests, not all minors.
    for k in (15,):
        a = [[int(i != j) for j in range(k + 1)] for i in range(k + 1)]
        require(abs(determinant(a)) == k, 'large simplex determinant')
        q = list(range(k + 1))
        x = [Fraction(sum(q), k) - z for z in q]
        require(matvec(a, x) == q, 'large simplex inverse')
    squares = [64, 81, 100, 121]
    x = [sum(squares) // 3 - q for q in squares]
    require(x == [58, 41, 22, 1], 'positive distinct square triple witness')
    require(matvec([[int(i != j) for j in range(4)] for i in range(4)], x) == squares, 'triple square witness')
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a in range(1, 8) for b in range(a + 1, 8)})
    n = [[int(p in line) for p in range(1, 8)] for line in lines]
    snf = smith_by_minors(n, 7)
    require(snf == [1, 1, 1, 1, 2, 2, 6], 'Fano point-line Smith')
    require(abs(determinant(n)) == 24, 'Fano determinant')
    require([rank_mod(n, 7, p) for p in (2, 3, 5)] == [4, 6, 7], 'Fano finite-field ranks')
    for i in range(7):
        for j in range(7):
            require(sum(n[i][k] * n[j][k] for k in range(7)) == (3 if i == j else 1), 'Fano Gram')
    tested = 0
    for q in product(range(3), repeat=7):
        total = sum(q)
        sums = [sum(n[i][p] * q[i] for i in range(7)) for p in range(7)]
        x = [Fraction(3 * s - total, 6) for s in sums]
        require(matvec(n, x) == list(q), 'Fano exact inverse')
        integral = all(t.denominator == 1 for t in x)
        criterion = total % 3 == 0 and all((s - total // 3) % 2 == 0 for s in sums)
        require(integral == criterion, 'Fano integral image')
        tested += 1
    return {'J_minus_I_Smith': simplex, 'k15_determinant': 15,
            'positive_distinct_square_triple_labels': [58, 41, 22, 1],
            'Fano_lines': lines, 'Fano_Smith': snf, 'Fano_inverse_rhs_tests': tested}


def fano_auxiliary_controls():
    h = [[(j >> i) & 1 for j in range(1, 8)] for i in range(3)]
    code = [x for x in product((0, 1), repeat=7) if all(v % 2 == 0 for v in matvec(h, x))]
    require(len(code) == 16, 'Hamming dimension')
    require(smith_by_minors(h, 7) == [1, 1, 1], 'H integral Smith')
    incidences = [(v, 7 + i) for i in range(3) for v in range(7) if h[i][v]]
    cycle_count, factor_count = 0, 0
    for y in product((0, 1), repeat=12):
        parity = [0] * 10
        for value, (u, v) in zip(y, incidences):
            parity[u] ^= value
            parity[v] ^= value
        cycle_count += not any(parity)
        equality = all(len({y[j] for j, (u, v) in enumerate(incidences) if u == i}) == 1 for i in range(7))
        if equality and not any(parity[7:]):
            factor_count += 1
            decoded = tuple(next(y[j] for j, (u, v) in enumerate(incidences) if u == i) for i in range(7))
            require(decoded in code, 'factor graph decode')
    require((cycle_count, factor_count) == (8, 16), 'Tanner graph law distinction')
    generated = set()
    for a, b, c, d in product((0, 1), repeat=4):
        generated.add((a ^ b ^ d, a ^ c ^ d, a, b ^ c ^ d, b, c, d))
    require(generated == set(code), 'four cycle bits plus decoder')
    return {'H_Smith': [1, 1, 1], 'Hamming_words': 16, 'Tanner_vertices': 10,
            'Tanner_edges': 12, 'ordinary_even_subgraphs': cycle_count,
            'equality_parity_assignments': factor_count, 'decoder_words': len(generated)}


def projection_controls():
    # Every simple graph on four labelled vertices, every retained/hidden split.
    pairs, tested = list(combinations(range(4), 2)), 0
    for roles in product((0, 1, 2), repeat=6):
        s = [e for e, r in zip(pairs, roles) if r == 1]
        t = [e for e, r in zip(pairs, roles) if r == 2]
        parent = list(range(4))
        def find(v):
            while parent[v] != v:
                v = parent[v]
            return v
        for u, v in t:
            parent[find(u)] = find(v)
        projected = set()
        for y in product((0, 1), repeat=len(s) + len(t)):
            parity = [0] * 4
            for bit, (u, v) in zip(y, s + t):
                parity[u] ^= bit
                parity[v] ^= bit
            if not any(parity):
                projected.add(y[:len(s)])
        contracted = set()
        for y in product((0, 1), repeat=len(s)):
            parity = [0] * 4
            for bit, (u, v) in zip(y, s):
                parity[find(u)] ^= bit
                parity[find(v)] ^= bit
            if not any(parity):
                contracted.add(y)
        require(projected == contracted, 'auxiliary cycle projection theorem')
        tested += 1
    return {'all_graph_retained_hidden_assignments_on_four_vertices': tested}


def guard_controls():
    triangle = signed_matrix(3, [(0, 1, -1), (1, 2, -1), (0, 2, -1)])
    rational_labels = [Fraction(1, 2)] * 3
    require(matvec(triangle, rational_labels) == [1, 1, 1], 'rational odd square hostile')
    require(all(v.denominator != 1 for v in rational_labels), 'integral odd square rejection')
    require(determinant(triangle) != 0, 'unique rational triangle solution')
    good = [6, 19, 30]
    require(sorted(good[i] + good[j] for i, j in combinations(range(3), 2)) == [25, 36, 49], 'triangle square sums')
    path_labels = [0, 1, 2]
    oriented = [path_labels[0] - path_labels[1], path_labels[2] - path_labels[1]]
    require(oriented == [-1, 1] and len(set(map(abs, oriented))) == 1, 'mirror factor hostile')
    return {'square_triangle_labels': good, 'graceful_mirror_hostile': oriented}


def main():
    result = {'signed_graphs': graph_controls(), 'petersen': petersen_controls(),
              'hypergraphs': hypergraph_controls(), 'auxiliaries': fano_auxiliary_controls(),
              'projection': projection_controls(), 'guards': guard_controls()}
    print(json.dumps(result, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
