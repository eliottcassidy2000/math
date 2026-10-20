#!/usr/bin/env python3
"""Literal labelled packing/stabilizer audit, independent of quotient verifier.

Exhaust every graph type through n=6 by full S_n orbits. Enumerate actual
disjoint twin-triple subsets, actual automorphism stabilizers, and the actual
contracted graph signs. Compare with a separately implemented quotient formula.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
import json

CHECKS = 0

def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)

def moved(mask, edge_images):
    result = 0
    while mask:
        bit = mask & -mask
        result |= edge_images[bit.bit_length()-1]
        mask -= bit
    return result

def original_sign(edges, p):
    return (-1) ** sum(p[u] > p[v] for u, v in edges)

def quotient_D(n, edges):
    neighbors = [frozenset(v if u == i else u for u, v in edges if i in (u, v))
                 for i in range(n)]
    classes = []
    for neighborhood in dict.fromkeys(neighbors):
        classes.append(tuple(i for i in range(n) if neighbors[i] == neighborhood))
    c = list(map(len, classes))
    E = {(i, j) for i, j in combinations(range(len(c)), 2)
         if classes[j][0] in neighbors[classes[i][0]]}
    group = [p for p in permutations(range(len(c)))
             if all(c[i] == c[p[i]] for i in range(len(c)))
             and {tuple(sorted((p[i], p[j]))) for i, j in E} == E]
    total = 0
    for p in group:
        value = (-1) ** sum(c[i]*c[j] for i, j in E if p[i] > p[j])
        used = set()
        for i in range(len(c)):
            if i in used:
                continue
            j, cycle = i, []
            while j not in used:
                used.add(j)
                cycle.append(j)
                j = p[j]
            k = c[i] // 3
            value *= sum((-1) ** (len(cycle)*a) for a in range(k+1))
        total += value
    return Fraction(total, len(group))

def literal_packing(n, edges, auts):
    nb = [frozenset(v if u == i else u for u, v in edges if i in (u, v))
          for i in range(n)]
    triples = [t for t in combinations(range(n), 3)
               if nb[t[0]] == nb[t[1]] == nb[t[2]]]
    all_packings = []
    for j in range(1, n//3 + 1):
        all_packings.extend(tuple(P) for P in combinations(triples, j)
                            if len(set(v for t in P for v in t)) == 3*j)
    def action(P, p):
        return tuple(sorted(tuple(sorted(p[v] for v in t)) for t in P))
    seen, b, checked_stabilizers = set(), 0, 0
    for P in all_packings:
        if P in seen:
            continue
        seen.update(action(P, p) for p in auts)
        covered = {v for t in P for v in t}
        blocks = list(P) + [(v,) for v in range(n) if v not in covered]
        owner = {v: i for i, B in enumerate(blocks) for v in B}
        contracted = {tuple(sorted((owner[u], owner[v]))) for u, v in edges}
        positive = True
        for p in auts:
            if action(P, p) != P:
                continue
            alpha = tuple(owner[p[B[0]]] for B in blocks)
            sign = original_sign(contracted, alpha)
            need(sign == original_sign(edges, p), 'literal odd contraction identity')
            positive &= sign == 1
            checked_stabilizers += 1
        if positive:
            b += (-1) ** (len(P)+1)
    even = all(original_sign(edges, p) == 1 for p in auts)
    e = int(even and bool(triples))
    D = quotient_D(n, edges)
    need(D.denominator == 1, 'quotient integrality')
    need(b == int(even) - D, 'literal packing differs from quotient formula')
    return b, e, checked_stabilizers

def graph_audit(n):
    pairs = list(combinations(range(n), 2))
    positions = {e: i for i, e in enumerate(pairs)}
    perms = list(permutations(range(n)))
    images = [[1 << positions[tuple(sorted((p[u], p[v])))] for u, v in pairs]
              for p in perms]
    seen = set()
    count = btotal = etotal = stabilizers = 0
    defects = Counter()
    for mask in range(1 << len(pairs)):
        if mask in seen:
            continue
        orbit = [moved(mask, image) for image in images]
        seen.update(orbit)
        auts = [p for p, transformed in zip(perms, orbit) if transformed == mask]
        edges = {e for i, e in enumerate(pairs) if mask >> i & 1}
        b, e, s = literal_packing(n, edges, auts)
        count += 1
        btotal += b
        etotal += e
        stabilizers += s
        defects[b-e] += 1
    need(count == [0, 1, 2, 4, 11, 34, 156][n], 'unlabelled graph count control')
    need(btotal == etotal, 'global defect sum')
    return dict(n=n, all_graph_types=count, sum_b=btotal, sum_e=etotal,
                defect_histogram=dict(sorted(defects.items())),
                explicit_packing_stabilizers_checked=stabilizers)

def four_quotient_audit():
    pairs = list(combinations(range(4), 2))
    rows = []
    for mask in range(64):
        E = {e for i, e in enumerate(pairs) if mask >> i & 1}
        nb = [frozenset(v if u == i else u for u, v in E if i in (u, v))
              for i in range(4)]
        if len(set(nb)) != 4:
            continue
        expanded = {(3*i+a, 3*j+b) for i, j in E for a in range(3) for b in range(3)}
        D = quotient_D(12, expanded)
        deg = sorted(map(len, nb))
        expected = -2 if deg == [1, 1, 2, 2] else 1 if deg == [1, 1, 1, 1] else 0
        need(D == expected, 'four-class boundary classification')
        rows.append(dict(mask=mask, degrees=deg, D=int(D)))
    return rows

def main():
    rows = [graph_audit(n) for n in range(1, 7)]
    boundary = four_quotient_audit()
    print(json.dumps(dict(status='PASS', checks=CHECKS, graph_rows=rows,
                          four_class_weight_three_boundary=boundary),
                     indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
