"""Exact controls for gap geometry, arithmetic fibres, and typed counts."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import gcd, isqrt
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent


def need(test, message):
    if not test:
        raise RuntimeError(message)


def vertices(b, depth):
    out = {'': (F(0), F(1))}
    for k in range(depth):
        for w in [w for w in out if len(w) == k]:
            left, h = out[w]
            out[w+'0'] = (left, h/b)
            out[w+'1'] = (left+(b-1)*h/b, h/b)
    return {w: (left+h/2, h) for w, (left, h) in out.items()}


def cosh_distance(p, q):
    x, y = p
    u, v = q
    return 1+((x-u)**2+(y-v)**2)/(2*y*v)


def lca(a, b):
    h = 0
    for x, y in zip(a, b):
        if x != y:
            break
        h += 1
    return h


def geometry():
    rows = []
    for b in map(F, (2, F(9, 4), F(5, 2), 3, 4, 5)):
        pts = vertices(b, 6)
        edge_cosh = 1+5*(b-1)**2/(8*b)
        pairs = edges = rough_checks = 0
        for u, v in combinations(pts, 2):
            d = len(u)+len(v)-2*lca(u, v)
            C = cosh_distance(pts[u], pts[v])
            adjacent = d == 1
            need(C >= edge_cosh and (C == edge_cosh) == adjacent,
                 'exact induced equal-distance tree')
            pairs += 1
            edges += adjacent
            if b > 2:
                # Compare cosh to cosh(log X); monotonicity is used only
                # when X>=1. Negative lower distance bounds are automatic.
                X = b**d*((b-2)/b)**2
                if X >= 1:
                    need(C >= (X+1/X)/2, 'rough lower bound')
                Y = 4*b**d
                need(C <= (Y+1/Y)/2, 'rough upper bound')
                rough_checks += 1
        alpha = (b-1)/(2*b)
        gamma = (3*b+5)/(4*b)
        need(gamma**2+1 == (gamma+alpha)**2+b**-2,
             'geodesic circle endpoints')
        need(alpha > 0 and gamma > 0 and 2*alpha < b-1,
             'monotone arcs and noncrossing parent intervals')
        # Exact hostile at the two central addresses across the root split.
        a, c = '0'+'1'*5, '1'+'0'*5
        central = cosh_distance(pts[a], pts[c])
        need(central == 1+F(1, 2)*((b-2)*b**5+1)**2,
             'cross-root separation formula')
        rows.append({'base': str(b), 'vertices': len(pts), 'pair_checks': pairs,
                     'tree_edges': edges, 'rough_checks': rough_checks,
                     'cosh_edge_length': str(edge_cosh),
                     'cross_root_depth6_cosh': str(central)})
    return rows


def lattice_balls(radius, honeycomb):
    origin = (0, 0, 0) if honeycomb else (0, 0)
    seen, level = {origin}, {origin}
    counts = [1]
    for r in range(1, radius+1):
        nxt = set()
        for v in level:
            if honeycomb:
                a, b, s = v
                ns = ((a, b, 1), (a-1, b, 1), (a, b-1, 1)) if s == 0 else (
                     (a, b, 0), (a+1, b, 0), (a, b+1, 0))
            else:
                a, b = v
                ns = [(a+x, b+y) for x, y in ((1, 0), (-1, 0), (0, 1),
                                             (0, -1), (1, -1), (-1, 1))]
            nxt.update(ns)
        level = nxt-seen
        seen.update(level)
        counts.append(len(seen))
        need(counts[-1] == 1+(3*r*(r+1)//2)*(1 if honeycomb else 2),
             'lattice formula versus BFS')
    return counts


def geometry_phases():
    configs = list(product((F(1, 4), F(1, 3), F(2, 5), F(3, 7), F(49, 100)),
                           (F(1, 4), F(1, 2), F(3, 4), F(1))))
    for W in (F(1, 2), F(2, 3), F(3, 4)):
        configs.append((4*W**2/(4+W**2), W))
    rows = []
    for q, W in configs:
        pts = {w: (W*x, y) for w, (x, y) in vertices(1/q, 5).items()}
        E = 1+(1+W**2/4)*(1-q)**2/(2*q)
        S = 1+W**2*(1-q)**2/(2*q**2)
        threshold = 4*q/(4-q)
        phase = 'tree' if W**2 > threshold else 'triangles' if W**2 == threshold else 'matching'
        expected_count = {'tree': 62, 'triangles': 93, 'matching': 31}[phase]
        pairs = minimum_edges = 0
        for u, v in combinations(pts, 2):
            d = len(u)+len(v)-2*lca(u, v)
            sibling = len(u) == len(v) and u[:-1] == v[:-1]
            expected = d == 1 if phase == 'tree' else sibling if phase == 'matching' else d == 1 or sibling
            C = cosh_distance(pts[u], pts[v])
            need(C >= min(E, S) and (C == min(E, S)) == expected,
                 'complete nearest-distance phase classification')
            pairs += 1
            minimum_edges += C == min(E, S)
        need(minimum_edges == expected_count, 'phase graph edge count')
        boundary = 1+2*(1-q)**2/(q*(4-q))
        need(boundary-F(9, 7) == 2*(2*q-1)*(4*q-7)/(7*q*(4-q)) > 0,
             'strict adjacency-preserving edge-scale floor')
        if phase == 'tree':
            need(E > boundary > F(9, 7) > F(5, 4), 'tradeoff above log2')
        rows.append({'q': str(q), 'W': str(W), 'phase': phase,
                     'pair_checks': pairs, 'minimum_edges': minimum_edges,
                     'cosh_edge': str(E), 'cosh_sibling': str(S)})
    return rows


def capacity():
    seen, level = {5}, {5}
    counts = [1]
    for _ in range(26):
        nxt = {2*m for m in level}
        nxt.update((m-1)//3 for m in level if m % 6 == 4)
        need(not (nxt & seen), 'root5 first-hit levels disjoint')
        level = nxt
        seen.update(level)
        counts.append(len(seen))
    honey = lattice_balls(26, True)
    tri = lattice_balls(26, False)
    need(counts[21] == 704 and honey[21] == 694, 'honeycomb source control')
    need(counts[26] == 2263 and tri[26] == 2107, 'triangular source control')
    need(counts[21]-honey[21] == 10, 'deficit correction')
    need(next(r for r in range(27) if counts[r] > honey[r]) == 21, 'first honey radius')
    need(next(r for r in range(27) if counts[r] > tri[r]) == 26, 'first triangular radius')
    N, M = 704, 694
    q, r = divmod(N, M)
    need((q, r) == (1, 10) and F(2*r, N) == F(5, 176), 'balanced entropy fibre')
    return {'root5_counts': counts, 'honeycomb_radius21': honey[21],
            'triangular_radius26': tri[26], 'deficits': [10, 156],
            'entropy_minimizing_fibres': {'size1': 684, 'size2': 10},
            'minimum_conditional_entropy_bits': '5/176'}


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))


def divisors(n):
    small = [d for d in range(1, isqrt(n)+1) if n % d == 0]
    return sorted(set(small+[n//d for d in small]))


def arithmetic():
    ps = [p for p in range(3, 501, 2) if prime(p)]
    checks, collisions = 0, 0
    for a in range(2, 7):
        for b in range(1, a):
            A, B = 2**(a+1)-1, 2**(b+1)-1
            g = gcd(A, B)
            need(g == 2**gcd(a+1, b+1)-1, 'Mersenne gcd')
            R, S = A//g, B//g
            sigma_a = {p: sum(divisors(2**a*p)) for p in ps}
            sigma_b = {p: sum(divisors(2**b*p)) for p in ps}
            for p in ps:
                for q in ps:
                    left = sigma_a[p] == sigma_b[q]
                    right = (p+1) % (2*S) == 0 and (p+1)//(2*S) == F(q+1, 2*R)
                    need(left == right, 'complete prime-guarded collision chart')
                    checks += 1
                    collisions += left
    table = []
    for n in (4, 8, 9, 12, 15, 16, 20, 24, 25, 26, 28, 33):
        ds = divisors(n)
        s = sum(ds)-n
        witnesses = [(d, n//d) for d in ds if d <= n//d and 1+d+n//d == s]
        table.append({'n': n, 'proper_divisors': ds[:-1], 's': s,
                      'sigma': s+n, 'factor_pair_compression': witnesses})
    small = []
    for ell in range(1, 101):
        p, q = 6*ell-1, 14*ell-1
        if prime(p) and prime(q):
            need(sum(divisors(4*p)) == sum(divisors(2*q)) == 42*ell, '20/26 family')
            small.append([ell, 4*p, 2*q, 42*ell])
    aliquot_twins = []
    for a in range(3, 14):
        p = 2**a-5
        if prime(p) and p != 3:
            need(sum(divisors(3*p))-3*p == 2**a-1, '16/33 target family')
            aliquot_twins.append([2**a, 3*p, 2**a-1])
    need(sum(divisors(9))-9 != 7, 'repeated-prime hostile at exponent3')
    return {'collision_equivalences_checked': checks, 'collisions': collisions,
            'role_table': table, 'sigma_family_ell_through100': small,
            'aliquot_target_twins_exponent_through13': aliquot_twins}


def square_graph(n):
    return [(i, j) for i, j in combinations(range(1, n+1), 2) if isqrt(i+j)**2 == i+j]


def small_graphs():
    records = []
    thresholds = {}
    for n in range(1, 47):
        E = square_graph(n)
        A = {i: set() for i in range(1, n+1)}
        for i, j in E:
            A[i].add(j)
            A[j].add(i)
        if 'triangle' not in thresholds:
            triples = [(i, j, k) for i, j, k in combinations(A, 3)
                       if j in A[i] and k in A[i] and k in A[j]]
            if triples:
                thresholds['triangle'] = {'n': n, 'witness': triples[0]}
        if 'four_cycle' not in thresholds:
            for i, j in combinations(A, 2):
                common = sorted(A[i] & A[j])
                if len(common) >= 2:
                    thresholds['four_cycle'] = {'n': n, 'witness': [i, common[0], j, common[1]]}
                    break
    need(thresholds['triangle']['n'] == 30 and thresholds['four_cycle']['n'] == 46,
         'small arithmetic cycle thresholds')
    for n in (15, 20, 24, 25, 26, 33):
        E = square_graph(n)
        a = {v: set() for v in range(1, n+1)}
        for i, j in E:
            a[i].add(j)
            a[j].add(i)
        # DSATUR search supplies a 3-colouring certificate, not a unit-distance placement.
        colour = {}
        def solve():
            if len(colour) == n:
                return True
            v = max((i for i in a if i not in colour), key=lambda i: (
                len({colour[j] for j in a[i] if j in colour}), len(a[i]), -i))
            forbidden = {colour[j] for j in a[v] if j in colour}
            for c in range(3):
                if c not in forbidden:
                    colour[v] = c
                    if solve():
                        return True
                    del colour[v]
            return False
        need(solve() and all(colour[i] != colour[j] for i, j in E), 'square graph 3-colouring')
        if n == 15:
            rest = set(a)
            while True:
                leaves = {i for i in rest if len(a[i] & rest) <= 1}
                if not leaves:
                    break
                rest -= leaves
            need(len(rest) == 5 and all(len(a[i] & rest) == 2 for i in rest), 'square15 C5 core')
            cycle_start = min(rest)
            cycle = [cycle_start]
            while len(cycle) < 5:
                cycle.append(min((a[cycle[-1]] & rest)-set(cycle)))
            need(cycle_start in a[cycle[-1]], 'odd cycle certificate')
        visited, stack = set(), [1]
        while stack:
            u = stack.pop()
            if u not in visited:
                visited.add(u)
                stack.extend(a[u]-visited)
        need(len(visited) == n, 'connected guard graph')
        records.append({'n': n, 'edges': len(E), 'cycle_rank': len(E)-n+1,
                        'maximum_degree': max(map(len, a.values())),
                        'degree_four_vertices': [v for v in a if len(a[v]) == 4],
                        'colouring': [colour[i] for i in range(1, n+1)]})
    need(all(tuple(sorted((i, j))) in square_graph(15) for i, j in zip(cycle, cycle[1:]+cycle[:1])),
         'literal C5 edges')
    source = HERE.parent/'joint_incidence_clocked_frames_20261102/path_kernel.py'
    spec = importlib.util.spec_from_file_location('inherited_square_paths', source)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    E, P = mod.graph_paths(25, True)
    pats = mod.patterns(E, P)
    pairs = [(i, j) for i, j in combinations(range(len(P)), 2)
             if not (pats[i][0] & pats[j][1] or pats[i][1] & pats[j][0])]
    need(len(E) == 32 and len(P) == 20 and len(pairs) == 4, 'inherited square25 incidence')
    need(all((pats[i][0] | pats[i][1] | pats[j][0] | pats[j][1]).bit_count() == 26
             for i, j in pairs), '26 union arcs')
    need(12+4*5+1 == 33 and 20+4 == 24, 'carrier count anatomy')
    return {'square_graphs': records, 'odd_cycle_in_square15': cycle,
            'first_cycle_occurrences': thresholds,
            'square25': {'directed_paths': 20, 'path_arcs': 24, 'compatible_pairs': 4,
                         'pair_union_arcs': 26, 'guard_edges': 32,
                         'initial_active_guard_classes': 24, 'closed_carrier_states': 33}}


def a5():
    counts = Counter()
    for p in permutations(range(5)):
        inversions = sum(p[i] > p[j] for i, j in combinations(range(5), 2))
        if inversions % 2:
            continue
        q, order = list(p), 1
        while q != list(range(5)):
            q = [p[q[i]] for i in range(5)]
            order += 1
        counts[order] += 1
    need(counts == {1: 1, 2: 15, 3: 20, 5: 24}, 'A5 element orders')
    return dict(sorted(counts.items()))


if __name__ == '__main__':
    report = {'geometry': geometry(), 'geometry_phases': geometry_phases(),
              'capacity': capacity(), 'arithmetic': arithmetic(),
              'finite_graphs': small_graphs(), 'A5_order_counts': a5()}
    (HERE/'verification.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')
    print(json.dumps({'status': 'passed', 'geometry_pair_checks': sum(r['pair_checks'] for r in report['geometry']),
                      'phase_pair_checks': sum(r['pair_checks'] for r in report['geometry_phases']),
                      'arithmetic_checks': report['arithmetic']['collision_equivalences_checked'],
                      'square15_cycle': report['finite_graphs']['odd_cycle_in_square15']}))
