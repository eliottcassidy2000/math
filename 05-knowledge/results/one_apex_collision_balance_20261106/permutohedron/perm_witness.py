"""S84 (THM-4891): explicit witnesses for partially coherent permutohedron orientations, with an independent check.

A witness is a table kappa: class -> which element of the pair comes first at the sink end.
  L1 (square-coherent): class = (a, b, S), S = set of elements before the pair.
  L2 (hexagon-coherent): class = (a, b, sigma), sigma = relative order of the other elements.
The independent check rebuilds the oriented Cayley graph of S_n (generators = adjacent transpositions)
directly from kappa, counts vertices of out-degree 0, and runs a topological sort.
"""
import itertools, json, sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from perm_relax import edges_and_classes
from perm_relax3 import solve


def table_from_bits(n, level, bits):
    perms, cls, lits = edges_and_classes(n, level)
    tab = {}
    for key, c in cls.items():
        pair = key if level == 'L0' else key[0]
        a, b = pair
        first = a if bits[c] == 1 else b
        if level == 'L1':
            ctx = sorted(key[1])
        elif level == 'L2':
            ctx = list(key[1])
        else:
            ctx = []
        tab[(a, b, tuple(ctx))] = first
    return tab


def check(n, level, tab):
    perms = list(itertools.permutations(range(n)))
    index = {p: k for k, p in enumerate(perms)}
    out = {p: [] for p in perms}
    for p in perms:
        for i in range(n - 1):
            a, b = p[i], p[i + 1]
            lo, hi = min(a, b), max(a, b)
            if level == 'L1':
                ctx = tuple(sorted(p[:i]))
            elif level == 'L2':
                ctx = tuple(x for x in p if x not in (a, b))
            first = tab[(lo, hi, ctx)]
            q = list(p); q[i], q[i + 1] = q[i + 1], q[i]
            q = tuple(q)
            # edge p--q points into the arrangement where 'first' precedes the other element of the pair
            if p[i] != first:
                out[p].append(q)        # p -> q
    sinks = [p for p in perms if not out[p]]
    # Kahn
    indeg = {p: 0 for p in perms}
    for p in perms:
        for q in out[p]:
            indeg[q] += 1
    stack = [p for p in perms if indeg[p] == 0]
    seen = 0
    while stack:
        p = stack.pop(); seen += 1
        for q in out[p]:
            indeg[q] -= 1
            if indeg[q] == 0:
                stack.append(q)
    # coherence checks: squares and hexagons
    def oriented(p, i):
        q = list(p); q[i], q[i + 1] = q[i + 1], q[i]
        return tuple(q) in out[p]
    sq_ok = hex_ok = True
    for p in perms:
        for i in range(n - 1):
            for j in range(n - 1):
                if abs(i - j) >= 2:
                    pj = list(p); pj[j], pj[j + 1] = pj[j + 1], pj[j]; pj = tuple(pj)
                    if oriented(p, i) != oriented(pj, i):
                        sq_ok = False
            if i + 2 <= n - 1:
                # hexagon on positions i, i+1, i+2: the pair at (i,i+1) with c at i+2 versus c at i, pair at (i+1,i+2)
                a, b, c = p[i], p[i + 1], p[i + 2]
                r = list(p); r[i], r[i + 1], r[i + 2] = c, a, b; r = tuple(r)
                if oriented(p, i) != oriented(r, i + 1):
                    hex_ok = False
    return len(sinks), seen == len(perms), sq_ok, hex_ok, sinks


def main():
    res = {}
    # n = 4, L2, 7 sinks, acyclic
    ok, it, bits = solve(4, 'L2', 7, acyclic=True)
    tab = table_from_bits(4, 'L2', bits)
    ns, ac, sq, hx, sinks = check(4, 'L2', tab)
    res['n4_L2_7'] = {'sinks': ns, 'acyclic': ac, 'square_coherent': sq, 'hexagon_coherent': hx,
                      'sink_list': [''.join(map(str, s)) for s in sinks],
                      'table': {'%d%d|%s' % (a, b, ''.join(map(str, ctx))): f for (a, b, ctx), f in sorted(tab.items())}}
    for k in (7, 21):
        ok, it, bits = solve(5, 'L1', k, acyclic=True)
        tab = table_from_bits(5, 'L1', bits)
        ns, ac, sq, hx, sinks = check(5, 'L1', tab)
        res['n5_L1_%d' % k] = {'sinks': ns, 'acyclic': ac, 'square_coherent': sq, 'hexagon_coherent': hx,
                               'sink_list': [''.join(map(str, s)) for s in sinks],
                               'table': {'%d%d|{%s}' % (a, b, ''.join(map(str, ctx))): f for (a, b, ctx), f in sorted(tab.items())}}
        ok, it, bits = solve(5, 'L2', k, acyclic=True)
        tab = table_from_bits(5, 'L2', bits)
        ns, ac, sq, hx, sinks = check(5, 'L2', tab)
        res['n5_L2_%d' % k] = {'sinks': ns, 'acyclic': ac, 'square_coherent': sq, 'hexagon_coherent': hx,
                               'sink_list': [''.join(map(str, s)) for s in sinks]}
    for name, r in res.items():
        print(name, {k: v for k, v in r.items() if k not in ('table', 'sink_list')}, 'sinks:', ' '.join(r['sink_list']))
    print('table n4_L2_7 (pair|order of others -> first element):', json.dumps(res['n4_L2_7']['table']))
    print('table n5_L1_7 (pair|{set before} -> first element):', json.dumps(res['n5_L1_7']['table']))


if __name__ == '__main__':
    main()
