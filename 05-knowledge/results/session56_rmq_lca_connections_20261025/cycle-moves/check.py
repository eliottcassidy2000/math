#!/usr/bin/env python3
"""Exact portable checks; Python standard library only. Bits 1 orient low->high."""
from itertools import combinations, permutations
from collections import defaultdict, Counter
from pathlib import Path
import json

def edges(n):
    return [(i, j) for i, j in combinations(range(n), 2) if j != i + 1]

def score(n, E, m):
    d = [0] * n
    for k, (i, j) in enumerate(E):
        d[i if m >> k & 1 else j] += 1
    return tuple(d)

def arc(E, idx, m, u, v):
    return bool(m >> idx[tuple(sorted((u, v)))] & 1) == (u < v)

def cycles(n, E, induced=False):
    S = set(E)
    out = []
    for k in range(3, n + 1):
        for vs in combinations(range(n), k):
            for p in permutations(vs[1:]):
                c = (vs[0],) + p
                if c[1] > c[-1]:
                    continue
                pairs = {tuple(sorted((c[i], c[(i + 1) % k]))) for i in range(k)}
                if not pairs <= S:
                    continue
                if induced and any((e in S - pairs for e in combinations(vs, 2))):
                    continue
                out.append(c)
    return out

def reverse_cycle(E, idx, m, c, moves):
    """Recursively implement reversal using only induced directed cycles."""
    if not all((arc(E, idx, m, c[i], c[(i + 1) % len(c)]) for i in range(len(c)))):
        raise AssertionError('Check failed at original line 29')
    for i in range(len(c)):
        for j in range(i + 2, len(c)):
            if i == 0 and j == len(c) - 1:
                continue
            if tuple(sorted((c[i], c[j]))) not in idx:
                continue
            a = c[i:j + 1]
            b = c[j:] + c[:i + 1]
            first, second = (b, a) if arc(E, idx, m, c[i], c[j]) else (a, b)
            m = reverse_cycle(E, idx, m, first, moves)
            return reverse_cycle(E, idx, m, second, moves)
    mask = sum((1 << idx[tuple(sorted((c[i], c[(i + 1) % len(c)])))] for i in range(len(c))))
    moves.append(list(c))
    return m ^ mask

def route(n, E, start, target):
    idx = {e: k for k, e in enumerate(E)}
    m = start
    moves = []
    q = 0
    if not score(n, E, m) == score(n, E, target):
        raise AssertionError('Check failed at original line 44')
    while m != target:
        out = defaultdict(list)
        for k, (i, j) in enumerate(E):
            if (m ^ target) >> k & 1:
                out[i if m >> k & 1 else j].append(j if m >> k & 1 else i)
        v = min(out)
        walk = []
        seen = {}
        while v not in seen:
            seen[v] = len(walk)
            walk.append(v)
            v = min(out[v])
        c = tuple(walk[seen[v]:])
        before = m
        m = reverse_cycle(E, idx, m, c, moves)
        q += 1
        if not m ^ before == sum((1 << idx[tuple(sorted((c[i], c[(i + 1) % len(c)])))] for i in range(len(c)))):
            raise AssertionError('Check failed at original line 54')
        if not (m ^ target).bit_count() < (before ^ target).bit_count():
            raise AssertionError('Check failed at original line 55')
    if not m == target:
        raise AssertionError('Check failed at original line 56')
    if not len(moves) <= (start ^ target).bit_count() - 2 * q:
        raise AssertionError('Check failed at original line 57')
    m = start
    s = score(n, E, m)
    for c in moves:
        if not len(c) <= 4:
            raise AssertionError('Check failed at original line 61')
        if not all((arc(E, idx, m, c[i], c[(i + 1) % len(c)]) for i in range(len(c)))):
            raise AssertionError('Check failed at original line 62')
        cm = 0
        for i in range(len(c)):
            cm ^= 1 << idx[tuple(sorted((c[i], c[(i + 1) % len(c)])))]
        if not sum((tuple(sorted(e)) in idx for e in combinations(c, 2))) == len(c):
            raise AssertionError('Check failed at original line 65')
        m ^= cm
        if not score(n, E, m) == s:
            raise AssertionError('Check failed at original line 66')
    if not m == target:
        raise AssertionError('Check failed at original line 67')
    return moves

def census(n):
    E = edges(n)
    idx = {e: k for k, e in enumerate(E)}
    C = cycles(n, E, True)
    if not all((len(c) <= 4 for c in C)):
        raise AssertionError('Check failed at original line 72')
    templates = []
    for c in C:
        mask = pat = 0
        for i, u in enumerate(c):
            v = c[(i + 1) % len(c)]
            b = 1 << idx[tuple(sorted((u, v)))]
            mask |= b
            if u < v:
                pat |= b
        templates.append((mask, pat))
    N = 1 << len(E)
    par = list(range(N))
    fibres = defaultdict(list)
    move_edges = 0

    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a
    for m in range(N):
        fibres[score(n, E, m)].append(m)
        for mask, pat in templates:
            if m & mask in (pat, mask ^ pat):
                t = m ^ mask
                if m < t:
                    par[find(m)] = find(t)
                    move_edges += 1
    for f in fibres.values():
        if not len({find(m) for m in f}) == 1:
            raise AssertionError('Check failed at original line 90')
    route_count = max_moves = 0
    longest = None
    for f in fibres.values():
        for m in f[1:]:
            moves = route(n, E, m, f[0])
            route_count += 1
            if len(moves) > max_moves:
                max_moves = len(moves)
                longest = {'start': m, 'target': f[0], 'moves': moves}
    return {'n': n, 'residual_edges': len(E), 'orientations': N, 'fibres': len(fibres), 'fibre_size_histogram': dict(sorted(Counter(map(len, fibres.values())).items())), 'induced_cycle_counts': dict(Counter(map(len, C))), 'undirected_move_edges': move_edges, 'constructive_routes_checked': route_count, 'maximum_constructed_route_length': max_moves, 'longest_example': longest}

def hp(n, E, m):
    idx = {e: k for k, e in enumerate(E)}

    def a(u, v):
        return u + 1 == v if abs(u - v) == 1 else arc(E, idx, m, u, v)
    return sum((all((a(p[i], p[i + 1]) for i in range(n - 1))) for p in permutations(range(n))))
report = {'status': 'FINITE-EXACT', 'universe': 'All orientations of K_n minus the path 0-1-...-(n-1), for every n=1..7; labelled residual outdegree fibres. Fixed tournament path arcs all forward.', 'census': []}
for n in range(1, 8):
    r = census(n)
    report['census'].append(r)
    print(n, r['orientations'], r['fibres'], r['maximum_constructed_route_length'], flush=True)
n = 5
E = edges(n)
allones = (1 << len(E)) - 1
m = allones
for e in [(0, 3), (1, 4), (2, 4)]:
    m ^= 1 << E.index(e)
f = [x for x in range(1 << len(E)) if score(n, E, x) == score(n, E, m)]
if not len(f) == 2:
    raise AssertionError('Check failed at original line 110')
w = []
for x in f:
    directed_triangles = [c for c in cycles(n, E) if len(c) == 3 and (all((arc(E, {e: k for k, e in enumerate(E)}, x, c[i], c[(i + 1) % 3]) for i in range(3))) or all((not arc(E, {e: k for k, e in enumerate(E)}, x, c[i], c[(i + 1) % 3]) for i in range(3))))]
    if not not directed_triangles:
        raise AssertionError('Check failed at original line 114')
    w.append({'mask': x, 'reversed_forward_edges': [list(e) for k, e in enumerate(E) if not x >> k & 1], 'Hamilton_paths': hp(n, E, x)})
if not sorted((v['Hamilton_paths'] for v in w)) == [11, 13]:
    raise AssertionError('Check failed at original line 116')
report['sharpness'] = {'n': 5, 'residual_edges': E, 'residual_score': score(n, E, m), 'full_score': [d + int(i < 4) for i, d in enumerate(score(n, E, m))], 'entire_fibre': w, 'route': route(n, E, f[0], f[1]), 'directed_triangle_moves': 0, 'minimal_n': True}
generic = []
for n in range(1, 6):
    allE = list(combinations(range(n), 2))
    tested = states = witnesses = 0
    for support in range(1 << len(allE)):
        E = [e for k, e in enumerate(allE) if support >> k & 1]
        idx = {e: k for k, e in enumerate(E)}
        C = cycles(n, E, True)
        N = 1 << len(E)
        par = list(range(N))
        fibres = defaultdict(list)

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a
        templates = []
        for c in C:
            mask = pat = 0
            for i, u in enumerate(c):
                v = c[(i + 1) % len(c)]
                bit = 1 << idx[tuple(sorted((u, v)))]
                mask |= bit
                if u < v:
                    pat |= bit
            templates.append((mask, pat))
        for m in range(N):
            fibres[score(n, E, m)].append(m)
            for mask, pat in templates:
                if m & mask in (pat, mask ^ pat):
                    par[find(m)] = find(m ^ mask)
        if not all((len({find(m) for m in f}) == 1 for f in fibres.values())):
            raise AssertionError('Check failed at original line 139')
        if C:
            c = max(C, key=len)
            inside = set(c)
            m = 0
            for k, (u, v) in enumerate(E):
                if u in inside and v in inside:
                    forward = any((c[i] == u and c[(i + 1) % len(c)] == v for i in range(len(c))))
                elif (u in inside) != (v in inside):
                    forward = u not in inside
                else:
                    forward = True
                if forward:
                    m |= 1 << k
            if not len(fibres[score(n, E, m)]) == 2:
                raise AssertionError('Check failed at original line 148')
            if not (fibres[score(n, E, m)][0] ^ fibres[score(n, E, m)][1]).bit_count() == len(c):
                raise AssertionError('Check failed at original line 149')
            witnesses += 1
        tested += 1
        states += N
    generic.append({'n': n, 'graphs': tested, 'orientation_states': states, 'sharp_induced_cycle_fibres': witnesses})
report['all_graph_checks'] = generic
Path(__file__).with_name('certificate.json').write_text(json.dumps(report, indent=2) + '\n')
print('PASS')
