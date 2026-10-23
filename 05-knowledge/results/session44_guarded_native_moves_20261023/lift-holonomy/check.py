#!/usr/bin/env python3
"""Exact supplied-cover compiler and independently authenticated native control.

No numeric order20 electrical potential is invented. Native phase preservation
is proved by explicit source/partner isomorphisms, for any supplied class potential.
"""
from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path
import argparse
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = next((p for p in HERE.parents if (p/'.git').exists()), HERE)
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def points(r):
    return list(product(range(3), repeat=r))


def eye(r):
    return tuple(tuple(int(i == j) for j in range(r)) for i in range(r))


def mv(A, x):
    return tuple(sum(a*b for a, b in zip(row, x)) % 3 for row in A)


def add(x, y):
    return tuple((a+b) % 3 for a, b in zip(x, y))


def mm(A, B):
    n = len(A)
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(n)) % 3
                       for j in range(n)) for i in range(n))


def matrix_inverse(A):
    n = len(A)
    rows = [list(row)+list(eye(n)[i]) for i, row in enumerate(A)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if rows[i][k]), None)
        if pivot is None:
            raise ValueError('singular affine transport')
        rows[k], rows[pivot] = rows[pivot], rows[k]
        inverse = pow(rows[k][k], -1, 3)
        rows[k] = [(x*inverse) % 3 for x in rows[k]]
        for i in range(n):
            if i != k:
                c = rows[i][k]
                rows[i] = [(x-c*y) % 3 for x, y in zip(rows[i], rows[k])]
    return tuple(tuple(row[n:]) for row in rows)


def identity(r):
    return eye(r), (0,)*r


def translation(t):
    return eye(len(t)), tuple(t)


def apply(f, x):
    return add(mv(f[0], x), f[1])


def compose(f, g):
    """f after g."""
    return mm(f[0], g[0]), add(mv(f[0], g[1]), f[1])


def inverse(f):
    A = matrix_inverse(f[0])
    return A, tuple(-x % 3 for x in mv(A, f[1]))


def partition(vertices, adjacency):
    remaining = set(vertices)
    parts = []
    while remaining:
        root = min(remaining)
        seen, queue = {root}, [root]
        while queue:
            u = queue.pop()
            for v in adjacency[u]:
                if v not in seen:
                    seen.add(v)
                    queue.append(v)
        remaining.difference_update(seen)
        parts.append(frozenset(seen))
    return frozenset(parts)


def validate_map(f, r):
    if len(f[0]) != r or any(len(row) != r for row in f[0]) or len(f[1]) != r:
        raise ValueError('wrong fibre dimensions')
    if any(type(x) is not int or x not in (0, 1, 2) for row in f[0] for x in row) or any(type(x) is not int or x not in (0, 1, 2) for x in f[1]):
        raise ValueError('use canonical integer F3 entries')
    inverse(f)


def compile_cover(vertices, r, edges, internal=()):
    """Topology compiler ONLY. Native legality must be authenticated separately.

    Every listed undirected edge means a complete matching on the full fibre;
    every internal generator acts on the full fibre. Reverse moves are included.
    """
    vertices = tuple(vertices)
    if not vertices or len(set(vertices)) != len(vertices):
        raise ValueError('need nonempty distinct base vertices')
    V = points(r)
    adjacency = {u: [] for u in vertices}
    for u, v, f in edges:
        if u not in adjacency or v not in adjacency:
            raise ValueError('edge outside base')
        validate_map(f, r)
        rev = inverse(f)
        adjacency[u].append((v, f))
        adjacency[v].append((u, rev))
    root = vertices[0]
    chart = {root: identity(r)}
    queue = deque([root])
    while queue:
        u = queue.popleft()
        for v, f in adjacency[u]:
            if v not in chart:
                chart[v] = compose(f, chart[u])
                queue.append(v)
    if len(chart) != len(vertices):
        raise ValueError('disconnected base: compile each component separately')
    holonomy = [compose(inverse(chart[v]), compose(f, chart[u])) for u, v, f in edges]
    for v, f in internal:
        if v not in chart:
            raise ValueError('invalid internal generator vertex')
        validate_map(f, r)
        holonomy.append(compose(inverse(chart[v]), compose(f, chart[v])))
    orbit_adjacency = {a: set() for a in V}
    for a in V:
        for h in holonomy:
            orbit_adjacency[a].add(apply(h, a))
            orbit_adjacency[a].add(apply(inverse(h), a))
    orbits = partition(V, orbit_adjacency)
    predicted = frozenset(frozenset((v, apply(chart[v], a)) for v in vertices for a in orbit)
                          for orbit in orbits)
    return chart, holonomy, orbits, predicted


def literal_components(vertices, r, edges, internal=()):
    states = [(v, a) for v in vertices for a in points(r)]
    adjacency = {s: set() for s in states}
    for u, v, f in list(edges)+[(v, v, f) for v, f in internal]:
        for a in points(r):
            s, t = (u, a), (v, apply(f, a))
            adjacency[s].add(t)
            adjacency[t].add(s)
    return partition(states, adjacency)


def check_cover(vertices, r, edges, internal=()):
    result = compile_cover(vertices, r, edges, internal)
    actual = literal_components(vertices, r, edges, internal)
    need(result[3] == actual, 'root holonomy orbit partition versus literal graph')
    return result


def rank(rows, r):
    A = [list(row) for row in rows]
    p = 0
    for j in range(r):
        pivot = next((i for i in range(p, len(A)) if A[i][j]), None)
        if pivot is None:
            continue
        A[p], A[pivot] = A[pivot], A[p]
        z = pow(A[p][j], -1, 3)
        A[p] = [z*x % 3 for x in A[p]]
        for i in range(len(A)):
            if i != p:
                z = A[i][j]
                A[i] = [(x-z*y) % 3 for x, y in zip(A[i], A[p])]
        p += 1
    return p


def affine_bank(r):
    result = []
    for values in product(range(3), repeat=r*r):
        A = tuple(tuple(values[i*r:(i+1)*r]) for i in range(r))
        try:
            matrix_inverse(A)
        except ValueError:
            continue
        result.extend((A, c) for c in points(r))
    return result


def algebra_controls():
    A1, A2 = affine_bank(1), affine_bank(2)
    need(len(A1) == 6 and len(A2) == 432, 'complete affine banks')
    count = 0
    # Every assignment of three AGL(1,3) edge maps on a triangle.
    for f, g, h in product(A1, repeat=3):
        check_cover((0, 1, 2), 1, [(0, 1, f), (1, 2, g), (2, 0, h)])
        count += 1
    # All AGL(2,3) holonomies, plus nontrivial vertex gauges and transported internals.
    for i, h in enumerate(A2):
        edges = [(0, 1, identity(2)), (1, 2, identity(2)), (2, 0, h)]
        internal = [(1, A2[(37*i+11) % len(A2)])]
        base = check_cover((0, 1, 2), 2, edges, internal)
        gauges = [A2[(13*i+j*113) % len(A2)] for j in range(3)]
        changed = [(u, v, compose(gauges[v], compose(f, inverse(gauges[u])))) for u, v, f in edges]
        new_internal = [(v, compose(gauges[v], compose(f, inverse(gauges[v])))) for v, f in internal]
        moved = check_cover((0, 1, 2), 2, changed, new_internal)
        expected = frozenset(frozenset((v, apply(gauges[v], a)) for v, a in part) for part in base[3])
        need(moved[3] == expected, 'affine gauge/converse conjugacy')
        count += 2
    translation_cases = 0
    for r in (1, 2, 3):
        for s, t in product(points(r), repeat=2):
            _, generators, orbits, _ = check_cover((0, 1), r,
                [(0, 1, identity(r)), (0, 1, translation(s))], [(1, translation(t))])
            d = rank([h[1] for h in generators], r)
            need(len(orbits) == 3**(r-d), 'translation rank component count')
            need({len(orbit) for orbit in orbits} == {3**d}, 'translation equal orbit sizes')
            annihilator = [k for k in points(r) if all(sum(x*y for x, y in zip(k, h[1])) % 3 == 0 for h in generators)]
            need(len(annihilator) == len(orbits), 'dual zero-frequency-sector obstruction count')
            translation_cases += 1
    _, _, orbits, _ = check_cover((0, 1), 1,
        [(0, 1, identity(1)), (0, 1, (((2,),), (0,)))])
    need(sorted(map(len, orbits)) == [1, 2], 'affine orbit sizes need not be ternary powers')
    flat = check_cover((0, 1), 1, [(0, 1, identity(1))])
    need(len(flat[2]) == 3, 'connected base with three disconnected lifts')
    for args in [((0, 1), 1, [], []), ((0, 1), 1, [(0, 1, (((0,),), (0,)))], [])]:
        try:
            compile_cover(*args)
        except ValueError:
            pass
        else:
            raise RuntimeError('invalid cover accepted')
    need(len(check_cover((0, 1), 0, [(0, 1, identity(0))])[2]) == 1, 'zero-dimensional fibre boundary')
    need(len(check_cover((0,), 2, [])[2]) == 9, 'single base vertex without invented internal moves')
    translated = check_cover((0, 1, 2), 1, [(0, 1, identity(1)), (1, 2, identity(1)), (2, 0, translation((1,)))])
    need(len(translated[2]) == 1 and not any(apply(translation((1,)), a) == a for a in points(1)), 'connected lift need not admit parallel section')
    # A guard may leave only one statewise edge, not a full matching.
    states = [(v, a) for v in (0, 1) for a in points(1)]
    adjacency = {s: set() for s in states}
    adjacency[(0, (0,))].add((1, (0,)))
    adjacency[(1, (0,))].add((0, (0,)))
    need(len(partition(states, adjacency)) == 5, 'partial-edge hostile has five components, not three')
    return dict(affine_cover_cases=count, translation_cover_cases=translation_cases,
                affine_orbit_size_hostile=[1, 2], partial_edge_actual_components=5,
                incorrectly_completed_matching_components=3)


def matmul(A, B):
    return [[sum(x*y for x, y in zip(row, column)) for column in zip(*B)] for row in A]


def charpoly(A):
    """Exact Newton trace recurrence, coefficients highest degree first."""
    n = len(A)
    current = [[int(i == j) for j in range(n)] for i in range(n)]
    traces = [0]
    for _ in range(n):
        current = matmul(current, A)
        traces.append(sum(current[i][i] for i in range(n)))
    coefficients = [Fraction(1)]
    for k in range(1, n+1):
        coefficients.append(-sum(coefficients[k-i]*traces[i] for i in range(1, k+1))/k)
    return coefficients


def polymul(f, g):
    out = [0]*(len(f)+len(g)-1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            out[i+j] += a*b
    return out


def root_of_unity(power):
    return ((1, 0), (0, 1), (-1, -1))[power % 3]


def cyclotomic_add(x, y):
    return x[0]+y[0], x[1]+y[1]


def cyclotomic_multiply(x, y):
    a, b = x
    c, d = y
    return a*c-b*d, a*d+b*c-b*d


def character(k, a):
    return root_of_unity(sum(x*y for x, y in zip(k, a)))


def character_controls():
    # Independent pair arithmetic in Q[w]/(w^2+w+1), never floating roots.
    affine_checks = 0
    for f in affine_bank(2):
        for k, a in product(points(2), repeat=2):
            transpose_k = tuple(sum(f[0][i][j]*k[i] for i in range(2)) % 3 for j in range(2))
            need(character(k, apply(f, a)) == cyclotomic_multiply(character(k, f[1]), character(transpose_k, a)),
                 'exact affine character permutation and phase factor')
            affine_checks += 1
    outgoing = {v: [(v, (0, 0), Fraction(1, 2))] for v in range(3)}
    for u, v, t in [(0, 1, (1, 0)), (1, 2, (0, 1)), (2, 0, (1, 1))]:
        outgoing[u].append((v, t, Fraction(1, 4)))
        outgoing[v].append((u, tuple(-x % 3 for x in t), Fraction(1, 4)))
    block_checks = 0
    for k in points(2):
        for target_vertex in range(3):
            for u, a in product(range(3), points(2)):
                direct, coefficient = (0, 0), (0, 0)
                for v, t, weight in outgoing[u]:
                    if v == target_vertex:
                        z = character(k, add(a, t))
                        direct = cyclotomic_add(direct, (weight*z[0], weight*z[1]))
                        z = character(k, t)
                        coefficient = cyclotomic_add(coefficient, (weight*z[0], weight*z[1]))
                need(direct == cyclotomic_multiply(character(k, a), coefficient), 'complete translation block action on basis')
                block_checks += 1
    return dict(affine_character_identities=affine_checks, translation_block_basis_entries=block_checks,
                coefficient_ring='Q[w]/(w^2+w+1)', floating_arithmetic=False)


def spectral_controls():
    B = [[Fraction(5, 9), Fraction(4, 9)], [Fraction(4, 9), Fraction(5, 9)]]
    need(charpoly(B) == [1, Fraction(-10, 9), Fraction(1, 9)], 'base gap exactly eight ninths')
    outputs = []
    for epsilon in (Fraction(0), Fraction(1, 27), Fraction(1, 81), Fraction(1, 243)):
        states = [(v, a) for v in range(2) for a in range(3)]
        P = [[Fraction(0) for _ in states] for _ in states]
        index = {s: i for i, s in enumerate(states)}
        for (v, a), row in zip(states, P):
            for w in range(2):
                row[index[(w, a)]] += B[v][w]
            row[index[(v, a)]] -= epsilon
            for delta in (1, 2):
                row[index[(v, (a+delta) % 3)]] += epsilon/2
        need(all(sum(row) == 1 and min(row) >= 0 for row in P), 'exact reversible stochastic lift')
        eigenvalues = [Fraction(1), Fraction(1, 9)]+[1-3*epsilon/2]*2+[Fraction(1, 9)-3*epsilon/2]*2
        expected = [1]
        for eigenvalue in eigenvalues:
            expected = polymul(expected, [1, -eigenvalue])
        need(charpoly(P) == expected, 'full lift exact characteristic polynomial')
        adjacency = {i: {j for j, value in enumerate(row) if value} for i, row in enumerate(P)}
        components = len(partition(range(6), adjacency))
        need(components == (3 if not epsilon else 1), 'spectral hostile topology')
        gap = 1-sorted(eigenvalues, reverse=True)[1]
        need(gap == 3*epsilon/2, 'global gap independent of fixed base gap')
        outputs.append(dict(epsilon=str(epsilon), components=components, gap=str(gap),
                            eigenvalues=[str(x) for x in eigenvalues]))
    return dict(base_eigenvalues=['1', '1/9'], fixed_base_gap='8/9', lifts=outputs)


def native_controls():
    # Source singleton 0; C9 cells at 1..9 and 10..18; sink singleton19.
    cells = [list(range(1, 10)), list(range(10, 19))]
    T = set()
    for u in range(20):
        for v in range(u+1, 20):
            if u == 0 or v == 19:
                forward = True
            elif (u-1)//9 == (v-1)//9:
                forward = 1 <= (v-u) % 9 <= 4
            else:
                forward = (u-1) % 3 == (v-10) % 3
            T.add((u, v) if forward else (v, u))
    def path(a, b):
        return (0,)+tuple(1+(i+a) % 9 for i in range(9))+tuple(10+(2+i+b) % 9 for i in range(9))+(19,)
    def valid(P, arcs=T):
        return len(P) == 20 and len(set(P)) == 20 and all((u, v) in arcs for u, v in zip(P, P[1:]))
    def partner(arcs, P):
        keep = set(zip(P, P[1:]))
        return {(u, v) if (u, v) in keep else (v, u) for u, v in arcs}
    def permutation(a, b):
        return (0,)+tuple(1+(i+a) % 9 for i in range(9))+tuple(10+(i+b) % 9 for i in range(9))+(19,)
    def transport(arcs, perm):
        return {(perm[u], perm[v]) for u, v in arcs}
    Topp = {(v, u) for u, v in T}
    P0 = path(0, 0)
    Q0 = partner(T, P0)
    records = {}
    for a, b in product(range(9), repeat=2):
        P, perm = path(a, b), permutation(a, b)
        need(valid(P) == ((a-b) % 3 == 0), 'literal native boundary iff residue compatibility')
        need((transport(T, perm) == T) == ((a-b) % 3 == 0), 'exhaustive supplied C9 chart language')
        if (a-b) % 3:
            continue
        Q = partner(T, P)
        need(tuple(perm[x] for x in P0) == P, 'path transporter owner by owner')
        need(transport(Q0, perm) == Q, 'actual partner isomorphism, hence every class potential preserved')
        need(partner(Q, P) == T, 'native partner involution')
        need(P[0] == 0 and P[-1] == 19 and perm[0] == 0 and perm[19] == 19, 'literal indexed endpoints pinned')
        reversed_path = tuple(reversed(P))
        need(valid(reversed_path, Topp), 'converse actual path')
        need(partner(Topp, reversed_path) == {(v, u) for u, v in Q}, 'actual converse partner square')
        q = a % 3
        records[(q, ((a-q)//3, (b-q)//3))] = (a, b, P)
    need(len(records) == 27 and len({P for a, b, P in records.values()}) == 27, 'all twenty-seven actual path lifts retained')
    edges = [(0, 1, identity(2)), (1, 2, identity(2)), (2, 0, translation((1, 1)))]
    compiled = check_cover((0, 1, 2), 2, edges)
    repaired = check_cover((0, 1, 2), 2, edges, [(q, translation((1, 0))) for q in range(3)])
    for q, coordinate in records:
        a, b, P = records[(q, coordinate)]
        for da, db in ((1, 1), (-1, -1), (3, 0), (-3, 0)):
            aa, bb = (a+da) % 9, (b+db) % 9
            qq = aa % 3
            target = (qq, ((aa-qq)//3, (bb-qq)//3))
            PP = records[target][2]
            perm = permutation(da % 9, db % 9)
            need(tuple(perm[x] for x in P) == PP and valid(PP), 'every proposed native move authenticated')
            need(transport(T, perm) == T and perm[0] == 0 and perm[19] == 19, 'move fixes source and actual ports')
            need(transport(partner(T, P), perm) == partner(T, PP), 'move preserves native target class')
        qq = (q+1) % 3
        expected = coordinate if q != 2 else add(coordinate, (1, 1))
        need(records[(qq, expected)][:2] == ((a+1) % 9, (b+1) % 9), 'actual wrap carry exactly affine holonomy')
    need(sorted(map(len, compiled[3])) == [9, 9, 9], 'native shift graph three nine-cycles')
    need(len(repaired[3]) == 1, 'authenticated kernel move joins all native lifts')
    labels = {(q, uv): (uv[1]-uv[0]) % 3 for q, uv in records}
    need(all(len({labels[state] for state in component}) == 1 for component in compiled[3]), 'missing kernel-difference component label')
    A = [[11*int(i == j)+8*int((j-i) % 9 in (1, 8)) for j in range(9)] for i in range(9)]
    expected = [1]
    for factor in ([1, -27], [1, -3], [1, -3], [1, -33, 171, 1293], [1, -33, 171, 1293]):
        expected = polymul(expected, factor)
    need(charpoly(A) == expected, 'exact native C9 component characteristic polynomial scaled by27')
    base = [[11 if i == j else 8 for j in range(3)] for i in range(3)]
    need(charpoly(base) == polymul([1, -27], polymul([1, -3], [1, -3])), 'native quotient gap eight ninths')
    digest = sha256(json.dumps(sorted(T), separators=(',', ':')).encode()).hexdigest()
    return dict(order=20, arcs=len(T), source_arc_sha256=digest,
                supplied_chart_candidates=81, compatible_native_lifts=27,
                base_vertices=3, fibre_dimension=2, fibre_size=9,
                fixed_actual_endpoint_ports=[0, 19], holonomy_translation=[1, 1],
                holonomy_rank=1, native_component_sizes=[9, 9, 9],
                component_label='(b-a)/3 mod3', additional_kernel_translation=[1, 0],
                repaired_holonomy_rank=2, repaired_components=1,
                base_gap='8/9', component_gap='16/27 * (1-cos(2*pi/9))',
                scaled_component_charpoly=[int(x) for x in charpoly(A)],
                scaled_component_factorization='(z-27)(z-3)^2(z^3-33z^2+171z+1293)^2',
                electrical_scope='All class potentials preserved by explicit source/partner isomorphisms; no numeric order20 potential asserted',
                native_edge_scope='Specified actual automorphism moves, not all native metagraph edges')


def source_pins(repo):
    data = json.loads((HERE/'source_pins.json').read_text())
    for record in data['sources']:
        raw = subprocess.check_output(['git', 'show', data['commit']+':'+record['path']], cwd=repo)
        need(sha256(raw).hexdigest() == record['sha256'], 'immutable inherited source pin '+record['path'])
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=ROOT, help='Repository holding pinned source commits')
    parser.add_argument('--write', action='store_true', help='Explicitly regenerate the frozen receipt')
    args = parser.parse_args()
    pins = source_pins(args.repo)
    out = dict(status='PROVED scoped statements; FINITE-EXACT controls',
               source_commit=pins['commit'], algebra=algebra_controls(),
               characters=character_controls(), spectral=spectral_controls(), native=native_controls())
    out['exact_check_count'] = CHECKS
    serialized = json.dumps(out, indent=2, sort_keys=True)+'\n'
    if args.write:
        (HERE/'receipt.json').write_text(serialized)
    else:
        need((HERE/'receipt.json').read_text() == serialized, 'frozen receipt replay')
    print(serialized, end='')


if __name__ == '__main__':
    main()
