#!/usr/bin/env python3
"""Independent finite audit; stdlib only, no imports/execution of producer code.

Usage: python3 independent_checks.py [path/to/exact_checks.json]
Checks quotient returns by multinomial counts, voltage returns by group-ring
dynamic programming, characteristic polynomials by Leibniz determinant, and
Hamiltonian paths by depth-first extension. Explicit failures survive -O.
"""
import json
import sys
from collections import Counter
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT = HERE.parent / 'exact_checks.json'
if not DEFAULT.exists():
    DEFAULT = HERE.parent / 'session16-geometric-forcing' / 'exact_checks.json'
INPUT = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
claimed = json.loads(INPUT.read_text())
checks = 0

def require(condition, name):
    global checks
    checks += 1
    if not condition:
        raise ValueError(name)

def norm(p):
    x, y = p
    return x*x + x*y + y*y

def residue(p):
    return (p[0] + 3*p[1]) % 7

S = tuple((x, y) for x, y in product(range(-2, 3), repeat=2) if norm((x, y)) == 1)
require(len(S) == 6, 'derive all six unit vectors from norm')
require(set(S) == {tuple(s) for s in claimed['unit_steps']}, 'producer unit vectors')
require({residue(s) for s in S} == set(range(1, 7)), 'quotient locally bijective')
require(all((x-y) % 3 for x, y in S), 'three-color differences')
require(all(norm((x, y)) % 7 == 0 for x, y in product(range(-7, 8), repeat=2)
            if residue((x, y)) == 0), 'kernel norm divisibility finite control')
require(min(norm(p) for p in product(range(-3, 4), repeat=2)
            if p != (0, 0) and residue(p) == 0) == 7, 'shortest period control')
require(claimed.get('triangular_lattice_chromatic_number') == 3,
        'lattice chromatic value is explicitly typed')
require('plane_chromatic_number' not in claimed, 'no unqualified plane=3 claim')

def kernel_coords(x, y):
    require((y-2*x) % 7 == 0 and (x+3*y) % 7 == 0, 'voltage belongs to kernel')
    return (y-2*x)//7, (x+3*y)//7

step_for = {residue(s): s for s in S}
voltage = {}
for r in range(7):
    for s in range(7):
        if r != s:
            x, y = step_for[(s-r) % 7]
            voltage[r, s] = kernel_coords(r+x-s, y)

triangles = []
gauge = {r: (r*r-3*r, 2*r+1) for r in range(7)}
for i, j, k in combinations(range(7), 3):
    edges = ((i, j), (j, k), (k, i))
    vx = sum(voltage[e][0] for e in edges)
    vy = sum(voltage[e][1] for e in edges)
    direct = tuple(sum(step_for[(s-r) % 7][axis] for r, s in edges) for axis in (0, 1))
    require(direct == (-3*vx+vy, vx+2*vy), 'triangle voltage/displacement equality')
    changed = tuple(sum(voltage[r, s][axis]+gauge[r][axis]-gauge[s][axis]
                        for r, s in edges) for axis in (0, 1))
    require(changed == (vx, vy), 'closed-walk gauge invariance')
    triangles.append({'vertices': [i, j, k], 'winding': list(direct)})
require(triangles == claimed['triangles'], 'complete triangle records')
zero_triangles = sum(t['winding'] == [0, 0] for t in triangles)
require(zero_triangles == 14, '14 true triangles, 21 winding artifacts')

def compositions(total, parts):
    if parts == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for rest in compositions(total-first, parts-1):
                yield (first,) + rest

def multinomial_returns(steps, n):
    counts = Counter()
    for multiplicities in compositions(n, len(steps)):
        displacement = tuple(sum(m*s[axis] for m, s in zip(multiplicities, steps))
                             for axis in (0, 1))
        coefficient = factorial(n)
        divisor = 1
        for m in multiplicities:
            divisor *= factorial(m)
        counts[displacement] += coefficient // divisor
    require(sum(counts.values()) == len(steps)**n, 'multinomial universe complete')
    return counts[(0, 0)], sum(c for p, c in counts.items() if residue(p) == 0)

# Group-ring adjacency at a fixed quotient root. Its states explicitly retain
# coordinates in the period basis, independently of direct step-count formulas.
gr = Counter({(0, 0, 0): 1})
undirected = []
for n in range(9):
    lift, quotient = multinomial_returns(S, n)
    gr_zero = gr[(0, 0, 0)]
    gr_total = sum(c for (r, _, _), c in gr.items() if r == 0)
    require((gr_zero, gr_total) == (lift, quotient), 'group ring vs multinomial')
    require(quotient == (6**n + 6*(-1)**n)//7, 'K7 spectral formula')
    record = claimed['return_moments'][n]
    require((lift, quotient) == (record['lattice_returns_per_vertex'], record['quotient_returns_per_vertex']),
            'producer undirected moments')
    undirected.append({'length': n, 'lattice': lift, 'quotient': quotient})
    nxt = Counter()
    for (r, a, b), count in gr.items():
        for s in range(7):
            if r != s:
                u, v = voltage[r, s]
                nxt[s, a+u, b+v] += count
    gr = nxt

def determinant_charpoly(A):
    n = len(A)
    ascending = [0]*(n+1)
    for p in permutations(range(n)):
        off = [i for i in range(n) if p[i] != i]
        if any(A[i][p[i]] == 0 for i in off):
            continue
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        ascending[n-len(off)] += (-1)**(inversions+len(off))
    return list(reversed(ascending))

def hamiltonian_from(A, start):
    paths = cycles = 0
    endpoints = Counter()
    cycle_winding = Counter()
    def extend(vertices):
        nonlocal paths, cycles
        if len(vertices) == len(A):
            paths += 1
            endpoints[vertices[-1]] += 1
            cycles += A[vertices[-1]][start]
            if A[vertices[-1]][start]:
                closed = vertices + (start,)
                displacement = tuple(sum(step_for[(s-r) % 7][axis]
                                         for r, s in zip(closed, closed[1:])) for axis in (0, 1))
                cycle_winding[displacement] += 1
            return
        for w in range(len(A)):
            if w not in vertices and A[vertices[-1]][w]:
                extend(vertices+(w,))
    extend((start,))
    return paths, cycles, dict(sorted(endpoints.items())), cycle_winding

orientation_results = {}
for name, steps in {'carousel': ((1, 0), (0, 1), (-1, 1)),
                    'paley': ((1, 0), (0, -1), (-1, 1))}.items():
    residues = {residue(s) for s in steps}
    A = [[int((s-r) % 7 in residues) for s in range(7)] for r in range(7)]
    require(all(A[r][s]+A[s][r] == int(r != s) for r in range(7) for s in range(7)),
            'orientation is tournament')
    counts = [hamiltonian_from(A, r) for r in range(7)]
    total_paths = sum(t[0] for t in counts)
    cycles = counts[0][1]
    require(len({(t[0], t[1]) for t in counts}) == 1, 'cyclic rooted-count symmetry')
    require(total_paths == claimed['orientations'][name]['paths'], 'Hamiltonian paths')
    require(cycles == claimed['orientations'][name]['cycles'], 'cycles modulo rotation, not reversal')
    require(counts[0][3][(0, 0)] == 0, 'neither quotient Hamiltonian cycle lifts closed')
    charpoly = determinant_charpoly(A)
    require(charpoly == claimed['orientations'][name]['characteristic_polynomial_descending'],
            'Leibniz determinant independent of Newton identities')
    moments = []
    for n in range(13):
        lift, quotient = multinomial_returns(steps, n)
        expected = int(n == 0) if name == 'carousel' else (
            factorial(n)//factorial(n//3)**3 if n % 3 == 0 else 0)
        require(lift == expected, 'orientation closed-walk formula')
        row = claimed['orientations'][name]['moments'][n]
        require((lift, quotient) == (row['zero_winding_per_vertex'], row['quotient_per_vertex']),
                'producer directed moments')
        moments.append({'length': n, 'lattice': lift, 'quotient': quotient})
    orientation_results[name] = {'paths': total_paths, 'cycles_mod_rotation': cycles,
        'paths_from_zero': counts[0][0], 'endpoint_counts_from_zero': counts[0][2],
        'cycle_winding_census': [{'displacement': list(p), 'count': c}
                                for p, c in sorted(counts[0][3].items())],
        'characteristic_polynomial_descending': charpoly, 'moments': moments}

# Equality at eight distinct integer inputs certifies the degree-seven factor
# identity independently of the producer's coefficient multiplication.
paley_polynomial = orientation_results['paley']['characteristic_polynomial_descending']
for t in range(-3, 5):
    val = 0
    for c in paley_polynomial:
        val = val*t+c
    require(val == (t-3)*(t*t+t+2)**3, 'Paley characteristic factorization')
u = [2, -1]
for n in range(2, 13):
    u.append(-u[-1]-2*u[-2])
for n, row in enumerate(orientation_results['paley']['moments']):
    require(3**n+3*u[n] == 7*row['quotient'], 'Paley recurrence vs multinomial count')
require(u == claimed['paley_recurrence']['u_through_12'], 'producer recurrence values')

require(all(x+2*y > 0 for x, y in ((1, 0), (0, 1), (-1, 1))), 'carousel strict height')
require(tuple(sum(s[a] for s in ((1, 0), (0, -1), (-1, 1))) for a in (0, 1)) == (0, 0),
        'Paley local zero relation')
result = {'status': 'PASS', 'checks': checks,
          'scope': 'Independent finite controls only; no new plane bound or Lean build.',
          'triangle_count': 35, 'zero_winding_triangles': zero_triangles,
          'zero_winding_trace_length_3': 7*undirected[3]['lattice'],
          'quotient_trace_length_3': 7*undirected[3]['quotient'],
          'undirected_moments': undirected, 'orientations': orientation_results}
out = json.dumps(result, indent=2)+'\n'
HERE.joinpath('independent_checks.json').write_text(out)
print(out, end='')
