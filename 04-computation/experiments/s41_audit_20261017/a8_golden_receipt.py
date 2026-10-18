#!/usr/bin/env python3
"""S41 audit A8 (independent referee code): golden_additive_boundary_20261017.md section 2, from clock_duality's receipt
04-computation/experiments/root_positive_transport_20261017g/golden_shared.json only.

1. Edge classification of the 30,856 receipt edges: Cartesian axis families (A x C x C, A x C x C) and the 56 mixed edges.
2. The 570 patterns are exactly the canonical proper 6-colourings of the induced port graph.
3. For every pattern, an additive colouring c_A(a) + c_C(b) + c_C(c) mod 6 with ONE coefficient colouring c_C (c_1 = c_2)
   is found with an encoding different from the session's (forbidden-sum clauses enumerated explicitly), and the full
   5400-vertex colouring is checked against every receipt edge and the literal port colours.
4. The two-distance graph on A (unit edges + mixed first-coordinate steps): the five pentagon ports span K_5, and the
   number of patterns extending to a proper 6-colouring of it (claimed 105 = 570 - 465).
Usage: python a8_golden_receipt.py
"""
import os, sys, json, itertools, time
sys.dont_write_bytecode = True
from pysat.solvers import Solver

HERE = os.path.dirname(os.path.abspath(__file__))
REC = os.path.join(os.path.dirname(HERE), 'root_positive_transport_20261017g', 'golden_shared.json')


def canonical_colourings(n, edges, k=6):
    out = []
    adj = {i: set() for i in range(n)}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)

    def rec(col):
        i = len(col)
        if i == n:
            out.append(tuple(col))
            return
        top = max(col) + 1 if col else 0
        for c in range(min(top, k - 1) + 1):
            if all(col[j] != c for j in adj[i] if j < i):
                rec(col + [c])
    rec([])
    return out


def main():
    t0 = time.time()
    R = json.load(open(REC))
    addr = [tuple(x) for x in R['source_addresses']]
    edges = [tuple(e) for e in R['edges']]
    mixed = set(tuple(e) for e in R['mixed_edges'])
    nA, nC = R['first_axis_vertices'], R['coefficient_vertices']
    pats = [tuple(p) for p in R['patterns']]
    res = dict(vertices=len(addr), edges=len(edges), mixed=len(mixed))
    assert len(set(addr)) == len(addr) == nA * nC * nC
    assert all(addr[p] == (p, 0, 0) for p in range(8)), 'ports are (p, 0, 0)'
    assert all(x == '0' for x in R['spindle_coordinates'][0]), 'coefficient index 0 is the zero element'
    # 1. classification
    EA, EU, EV, MX = set(), set(), set(), []
    famA, famU, famV = {}, {}, {}
    for e in edges:
        (a, b, c), (a2, b2, c2) = addr[e[0]], addr[e[1]]
        ch = (a != a2, b != b2, c != c2)
        if e in mixed or tuple(reversed(e)) in mixed:
            MX.append(((a, b, c), (a2, b2, c2)))
            assert ch == (True, True, True), 'mixed edges change all three coordinates'
        elif ch == (True, False, False):
            EA.add(frozenset((a, a2)))
            famA.setdefault((b, c), set()).add(frozenset((a, a2)))
        elif ch == (False, True, False):
            EU.add(frozenset((b, b2)))
            famU.setdefault((a, c), set()).add(frozenset((b, b2)))
        elif ch == (False, False, True):
            EV.add(frozenset((c, c2)))
            famV.setdefault((a, b), set()).add(frozenset((c, c2)))
        else:
            raise RuntimeError('unclassified edge %s' % (e,))
    cart = (len(famA) == nC * nC and all(v == EA for v in famA.values()) and len(famU) == nA * nC
            and all(v == EU for v in famU.values()) and len(famV) == nA * nC and all(v == EV for v in famV.values()))
    res.update(EA=len(EA), EU=len(EU), EV=len(EV), same_coefficient_graph=(EU == EV), cartesian_axis_families=cart,
               mixed_classified=len(MX))
    print('1. A edges %d, u-axis %d, v-axis %d (same graph %s), every axis family Cartesian: %s; mixed %d (all change '
          'three coordinates)' % (len(EA), len(EU), len(EV), EU == EV, cart, len(MX)))
    assert len(EA) * nC * nC + 2 * len(EU) * nA * nC + len(MX) == len(edges)
    # 2. patterns
    port_edges = sorted(tuple(sorted(e)) for e in edges if max(e) < 8)
    can = canonical_colourings(8, port_edges)
    res['patterns_equal_canonical'] = (sorted(can) == sorted(pats))
    print('2. induced port graph %s; canonical proper 6-colourings: %d; receipt patterns identical: %s'
          % (port_edges, len(can), sorted(can) == sorted(pats)))
    # 3. additive colourings with c_1 = c_2
    def var(kind, i, c):
        return (0 if kind == 'A' else nA * 6) + 6 * i + c + 1
    base = []
    for kind, n in (('A', nA), ('C', nC)):
        for i in range(n):
            base.append([var(kind, i, c) for c in range(6)])
            for c1, c2 in itertools.combinations(range(6), 2):
                base.append([-var(kind, i, c1), -var(kind, i, c2)])
    for e in EA:
        a, a2 = tuple(e)
        for c in range(6):
            base.append([-var('A', a, c), -var('A', a2, c)])
    for e in EU:
        b, b2 = tuple(e)
        for c in range(6):
            base.append([-var('C', b, c), -var('C', b2, c)])
    base.append([var('C', 0, 0)])
    # mixed: forbid every assignment of the involved variables whose signed sum is 0 mod 6
    for (a, b, c), (a2, b2, c2) in MX:
        terms = [('A', a2, 1), ('A', a, -1), ('C', b2, 1), ('C', b, -1), ('C', c2, 1), ('C', c, -1)]
        coef = {}
        for kind, i, s in terms:
            coef[(kind, i)] = coef.get((kind, i), 0) + s
        keys = [k for k, s in coef.items() if s % 6 != 0]
        for cols in itertools.product(range(6), repeat=len(keys)):
            if sum(coef[k] * cc for k, cc in zip(keys, cols)) % 6 == 0:
                base.append([-var(k[0], k[1], cc) for k, cc in zip(keys, cols)])
    found = 0
    allok = True
    with Solver(name='cd15', bootstrap_with=base) as S:
        for pat in pats:
            ok = S.solve(assumptions=[var('A', p, pat[p]) for p in range(8)])
            if not ok:
                allok = False
                continue
            m = set(x for x in S.get_model() if x > 0)
            cA = [next(c for c in range(6) if var('A', i, c) in m) for i in range(nA)]
            cC = [next(c for c in range(6) if var('C', i, c) in m) for i in range(nC)]
            col = [(cA[a] + cC[b] + cC[c]) % 6 for (a, b, c) in addr]
            good = all(col[x] != col[y] for x, y in edges) and tuple(col[:8]) == pat
            allok &= good
            found += good
    res['additive_shared_found_and_verified'] = found
    print('3. patterns with a verified additive colouring with c_1 = c_2 = c_C: %d of %d (every receipt edge checked) [%.0f s]'
          % (found, len(pats), time.time() - t0))
    # 4. two-distance graph on A
    TD = set(EA) | set(frozenset((x[0], y[0])) for x, y in MX)
    k5 = all(frozenset((i, j)) in TD for i, j in itertools.combinations(range(5), 2))
    base2 = []
    for i in range(nA):
        base2.append([6 * i + c + 1 for c in range(6)])
    for e in TD:
        a, a2 = tuple(e)
        for c in range(6):
            base2.append([-(6 * a + c + 1), -(6 * a2 + c + 1)])
    ext = 0
    with Solver(name='cd15', bootstrap_with=base2) as S:
        for pat in pats:
            ext += S.solve(assumptions=[6 * p + pat[p] + 1 for p in range(8)])
    res.update(two_distance_edges=len(TD), pentagon_K5=k5, two_distance_extending=ext, two_distance_blocked=len(pats) - ext)
    print('4. two-distance graph on A: %d edges; ports 0..4 span K_5: %s; patterns extending %d, blocked %d'
          % (len(TD), k5, ext, len(pats) - ext))
    json.dump(res, open(os.path.join(HERE, 'a8_golden_receipt.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
