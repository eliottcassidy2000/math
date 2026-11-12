#!/usr/bin/env python3
"""Exact finite controls for the all-order theorem in PROOF.md. Standard library."""
import itertools as it
import json
from fractions import Fraction
from math import factorial
from pathlib import Path

def words(k):
    return [w for m in range(k + 1) for w in it.permutations(range(k), m)]

def swap(w, i):
    return w[:i] + (w[i + 1], w[i]) + w[i + 2:]

def neighbors(w, k):
    out = [swap(w, i) for i in range(len(w) - 1)]
    if w:
        out.append(w[:-1])
    out.extend((w + (a,) for a in range(k) if a not in w))
    return out

def tubing(w, k):
    unused = set(range(k)) - set(w)
    tubes = [frozenset([a]) for a in unused]
    for j in range(1, len(w) + 1):
        tubes.append(frozenset({k} | set(range(k)) - set(w[:j])))
    return frozenset(tubes)

def compatible(a, b, k):
    if a <= b or b <= a:
        return True
    return not a & b and k not in a | b

def physical(core, signatures, w):
    r = len(core)
    names = tuple(range(r)) + tuple((r + a for a in w))
    arcs = {(i, j) for i in range(r) for j in range(r) if core[i] >> j & 1}
    for a in w:
        for i in range(r):
            arcs.add((r + a, i) if signatures[a] >> i & 1 else (i, r + a))
    arcs.update(((r + a, r + b) for i, a in enumerate(w) for b in w[i + 1:]))
    return (names, frozenset(arcs))

def hamiltonian_paths(graph):
    names, arcs = graph
    return sum((all(((p[i], p[i + 1]) in arcs for i in range(len(p) - 1))) for p in it.permutations(names)))

def main():
    reports = []
    for k in range(1, 6):
        ws = words(k)
        if not len(ws) == sum((factorial(k) // factorial(k - m) for m in range(k + 1))):
            raise AssertionError('exact control failed')
        if not len(ws) == k * len(words(k - 1)) + 1:
            raise AssertionError('exact control failed')
        edges = {tuple(sorted((w, v))) for w in ws for v in neighbors(w, k)}
        if not all((len(neighbors(w, k)) == k for w in ws)):
            raise AssertionError('exact control failed')
        if not len(edges) * 2 == k * len(ws):
            raise AssertionError('exact control failed')
        ts = {w: tubing(w, k) for w in ws}
        if not len(set(ts.values())) == len(ws):
            raise AssertionError('exact control failed')
        if not all((len(t) == k for t in ts.values())):
            raise AssertionError('exact control failed')
        if not all((all((compatible(a, b, k) for a, b in it.combinations(t, 2))) for t in ts.values())):
            raise AssertionError('exact control failed')
        if k <= 4:
            if not edges == {tuple(sorted((w, v))) for w, v in it.combinations(ws, 2) if len(ts[w] & ts[v]) == k - 1}:
                raise AssertionError('exact control failed')
        ordered = sorted(edges)
        weights = {e: Fraction((j * j + 7 * j + 3) % 31 - 15, j % 5 + 1) for j, e in enumerate(ordered)}

        def d(u, v):
            e = tuple(sorted((u, v)))
            return weights[e] if (u, v) == e else -weights[e]

        def circ(cycle):
            return sum((d(cycle[j], cycle[(j + 1) % len(cycle)]) for j in range(len(cycle))), Fraction())
        F = {w: sum((d(w[:j], w[:j + 1]) for j in range(len(w))), Fraction()) for w in ws}
        checks = 0
        for w in ws:
            for i in range(len(w) - 1):
                v = swap(w, i)
                p = w[:i]
                a, b = w[i:i + 2]
                base = p + (a, b)
                other = p + (b, a)
                total = circ([p, p + (a,), base, other, p + (b,)])
                for c in w[i + 2:]:
                    total -= circ([base, other, other + (c,), base + (c,)])
                    base, other = (base + (c,), other + (c,))
                if not total == d(w, v) - (F[v] - F[w]):
                    raise AssertionError('exact control failed')
                checks += 1
        reports.append({'k': k, 'vertices': len(ws), 'edges': len(edges), 'oriented_swap_filling_checks': checks})
    k = 3
    special = ((0, 1, 2), (1, 0, 2))

    def defect(u, v):
        return 1 if (u, v) == special else -1 if (v, u) == special else 0

    def sumdef(c):
        return sum((defect(c[j], c[(j + 1) % len(c)]) for j in range(len(c))))
    pentagons = []
    for p in words(k):
        for a, b in it.combinations(set(range(k)) - set(p), 2):
            cyc = [p, p + (a,), p + (a, b), p + (b, a), p + (b,)]
            if not sumdef(cyc) == 0:
                raise AssertionError('exact control failed')
            pentagons.append(cyc)
    square = [(0, 1), (1, 0), (1, 0, 2), (0, 1, 2)]
    if not (len(pentagons) == 6 and sumdef(square) == -1):
        raise AssertionError('exact control failed')
    hexagon = [(0, 1, 2), (1, 0, 2), (1, 2, 0), (2, 1, 0), (2, 0, 1), (0, 2, 1)]
    if not sumdef(hexagon) == 1:
        raise AssertionError('exact control failed')
    core = [0, 1]
    sigs = (1, 2, 2)
    ws = words(3)
    graphs = {w: physical(core, sigs, w) for w in ws}
    for w in ws:
        for v in neighbors(w, 3):
            nw, aw = graphs[w]
            nv, av = graphs[v]
            if len(w) == len(v):
                if not (set(nw) == set(nv) and len(aw ^ av) == 2):
                    raise AssertionError('exact control failed')
                diff = list(aw ^ av)
                if not diff[0] == diff[1][::-1]:
                    raise AssertionError('exact control failed')
            else:
                common = set(nw) & set(nv)
                if not {(a, b) for a, b in aw if a in common and b in common} == {(a, b) for a, b in av if a in common and b in common}:
                    raise AssertionError('exact control failed')
    hs = {w: hamiltonian_paths(graphs[w]) for w in ws}
    if not hs[0, 1] == hs[1, 0] == 5:
        raise AssertionError('exact control failed')
    if not (hs[0, 1, 2] == 11 and hs[1, 0, 2] == 13):
        raise AssertionError('exact control failed')
    report = {'grammar': reports, 'pentagon_only_hostile': {'all_six_pentagon_defects': 0, 'square_defect': sumdef(square), 'hexagon_defect': sumdef(hexagon)}, 'scalar_quotient_hostile': {'core_rows': core, 'tail_beaten_core_masks': sigs, 'AB': 5, 'BA': 5, 'ABC': 11, 'BAC': 13}, 'all_k3_path_counts': {' '.join(map(str, w)) or 'empty': hs[w] for w in ws}}
    print(json.dumps(report, sort_keys=True, indent=2))
if __name__ == '__main__':
    main()
