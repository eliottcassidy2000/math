#!/usr/bin/env python3
"""Exact sharpness control for the full-row-rank four-color ceiling.

Standard library only. Constructs the Moser spindle over F_17(theta), theta^2=3,
checks the entire induced unit graph, computes its rigidity rank independently,
and solves the first two carry systems to exhibit lifts modulo 17^2 and 17^3.
"""
from itertools import combinations, product
import json
from pathlib import Path

P, D = 17, 3
VERTICES = [(0, 0), (1, 0), (9, 2), (10, 2), (15, 1), (5, 5), (3, 6)]
EXPECTED_EDGES = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3),
                  (2, 3), (3, 6), (4, 5), (4, 6), (5, 6)]


def norm(v):
    return v[0]**2 - D*v[1]**2


def diff(a, b):
    return (a[0]-b[0], a[1]-b[1])


def solve(rows, rhs):
    """RREF modulo P, free variables zero; return rank and one exact solution."""
    a = [[x % P for x in row] + [b % P] for row, b in zip(rows, rhs)]
    pivots = []
    n = len(rows[0])
    for j in range(n):
        h = len(pivots)
        pivot = next((i for i in range(h, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[h], a[pivot] = a[pivot], a[h]
        inv = pow(a[h][j], -1, P)
        a[h] = [(x*inv) % P for x in a[h]]
        for i in range(len(a)):
            if i != h:
                c = a[i][j]
                a[i] = [(x-c*y) % P for x, y in zip(a[i], a[h])]
        pivots.append(j)
        if len(pivots) == len(a):
            break
    assert all(any(row[:n]) or row[n] == 0 for row in a)
    answer = [0]*n
    for i, j in enumerate(pivots):
        answer[j] = a[i][n]
    assert all((sum(x*y for x, y in zip(row, answer))-b) % P == 0
               for row, b in zip(rows, rhs))
    return len(pivots), answer


def main():
    assert pow(D, (P-1)//2, P) == P-1
    assert len(set(VERTICES)) == 7
    edges = [(i, j) for i, j in combinations(range(7), 2)
             if norm(diff(VERTICES[j], VERTICES[i])) % P == 1]
    assert edges == EXPECTED_EDGES
    rows = []
    for i, j in edges:
        dx, dy = diff(VERTICES[j], VERTICES[i])
        row = [0]*14
        row[2*i], row[2*j] = -2*dx, 2*dx
        row[2*i+1], row[2*j+1] = 2*D*dy, -2*D*dy
        rows.append(row)
    rank, _ = solve(rows, [0]*len(edges))
    assert rank == len(edges) == 11
    proper = lambda colors: all(colors[i] != colors[j] for i, j in edges)
    # This universe has exactly 3^7=2187 labeled color assignments.
    assert not any(proper(c) for c in product(range(3), repeat=7))
    colors = [0, 1, 2, 0, 1, 2, 3]
    assert proper(colors)
    # Every induced subgraph obeys the hereditary sparsity inequality.
    checked = 0
    for mask in range(1 << 7):
        selected = [i for i in range(7) if mask >> i & 1]
        if len(selected) < 2:
            continue
        m = sum(i in selected and j in selected for i, j in edges)
        assert m <= 2*len(selected)-3
        checked += 1
    lift = list(VERTICES)
    lifts = []
    for r in [1, 2]:
        modulus = P**r
        errors = [norm(diff(lift[j], lift[i]))-1 for i, j in edges]
        assert all(e % modulus == 0 for e in errors)
        rhs = [-(e//modulus) % P for e in errors]
        _, h = solve(rows, rhs)
        lift = [(x+modulus*h[2*i], y+modulus*h[2*i+1])
                for i, (x, y) in enumerate(lift)]
        assert all((norm(diff(lift[j], lift[i]))-1) % (P*modulus) == 0
                   for i, j in edges)
        assert [(x % P, y % P) for x, y in lift] == VERTICES
        lifts.append({'depth': r+1, 'modulus': P*modulus, 'vertices': lift})
    out = {'p': P, 'quadratic_nonsquare': D, 'vertices': VERTICES, 'edges': edges,
           'rigidity_rank': rank, 'stress_dimension': len(edges)-rank,
           'chromatic_number': 4, 'three_color_assignments_checked': 3**7,
           'four_coloring': colors, 'induced_sparsity_checks': checked,
           'exact_lifts': lifts,
           'scope': 'This is a sharp four-chromatic positive control, not a five- or six-chromatic witness.'}
    Path(__file__).with_name('rigidity-barrier-verified.json').write_text(json.dumps(out, indent=2)+'\n')
    print('EXACT RIGIDITY BARRIER CONTROL PASSED')
    print('Moser spindle over F_17^2: 7 distinct vertices, exactly 11 unit edges, rank 11.')
    print('All 2187 labeled three-color assignments fail; explicit four-coloring passes.')
    print('All 120 induced subgraphs with at least two vertices satisfy m <= 2v-3.')
    print('Exact carry solutions give distinct-vertex lifts modulo 17^2 and 17^3.')


if __name__ == '__main__':
    main()
