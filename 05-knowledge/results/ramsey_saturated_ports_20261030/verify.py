"""Exact native checks. No solver or graph package is needed to check witnesses.

Graph-bank completeness is CITED, not established by checking its members.
No 56/90-vertex Ramsey coloring is constructed or excluded here.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
GATES = 0


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def graph(n, predicate):
    return tuple(sum(1 << j for j in range(n) if i != j and predicate(i, j))
                 for i in range(n))


def g6(line):
    data = [v-63 for v in line.strip()]
    n = data.pop(0)
    need(n < 63, 'Small graph6 input')
    bits = [(v >> k) & 1 for v in data for k in range(5, -1, -1)]
    rows = [0]*n
    for bit, (i, j) in zip(bits, ((i, j) for j in range(n) for i in range(j))):
        if bit:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    return tuple(rows)


def sets_of_size(A, k, edge):
    return [S for S in combinations(range(len(A)), k)
            if all(((A[u] >> v) & 1) == edge for u, v in combinations(S, 2))]


def iso_maps(A, B, first=False):
    """Exact backtracking with dynamic minimum-candidate choice."""
    n = len(A)
    if len(B) != n:
        return []
    da, db = [r.bit_count() for r in A], [r.bit_count() for r in B]
    if sorted(da) != sorted(db):
        return []
    chosen, used, out = {}, set(), []
    def visit():
        if len(chosen) == n:
            out.append(tuple(chosen[u] for u in range(n)))
            return first
        choices = []
        for u in range(n):
            if u in chosen:
                continue
            candidates = [v for v in range(n) if v not in used and da[u] == db[v]
                          and all(((A[u] >> a) & 1) == ((B[v] >> b) & 1)
                                  for a, b in chosen.items())]
            if not candidates:
                return False
            choices.append((len(candidates), u, candidates))
        _, u, candidates = min(choices)
        for v in candidates:
            chosen[u] = v
            used.add(v)
            if visit():
                return True
            used.remove(v)
            del chosen[u]
        return False
    visit()
    return out


def orbit_sizes(group, n):
    todo, out = set(range(n)), []
    while todo:
        u = min(todo)
        orbit = {p[u] for p in group}
        out.append(len(orbit))
        todo -= orbit
    return sorted(out)


def banks():
    specs = [('r35_13.g6', 1, 3, 5), ('r36_17.g6', 7, 3, 6), ('r44_17.g6', 1, 4, 4)]
    out, all_graphs = {}, {}
    for name, count, s, t in specs:
        raw = (HERE/'inputs'/name).read_bytes()
        gs = [g6(line) for line in raw.splitlines() if line]
        need(len(gs) == count, 'Bank size')
        entries = []
        for i, A in enumerate(gs):
            need(not sets_of_size(A, s, 1), 'No forbidden clique')
            need(not sets_of_size(A, t, 0), 'No forbidden independent set')
            need(all(not iso_maps(A, other, True) for other in gs[:i]), 'Pairwise nonisomorphic')
            au = iso_maps(A, A)
            entries.append({'index': i, 'n': len(A), 'edges': sum(r.bit_count() for r in A)//2,
                            'degrees': dict(sorted(Counter(r.bit_count() for r in A).items())),
                            'aut_order': len(au), 'root_orbits': orbit_sizes(au, len(A)),
                            'ordered_pair_orbits': sum(sum(p[u] == u for u in range(len(A)))**2 for p in au)//len(au)})
        out[name] = {'sha256': hashlib.sha256(raw).hexdigest(), 'entries': entries}
        all_graphs[name] = gs
    need(sum(len(e['root_orbits']) for e in out['r36_17.g6']['entries']) == 67,
         'Seven unrooted 17-vertex graphs require 67 rooted types')
    need(sum(e['ordered_pair_orbits'] for e in out['r36_17.g6']['entries']) == 1065,
         'Exact joint two-root inventory')
    H = graph(13, lambda i, j: (i-j) % 13 in {1, 5, 8, 12})
    need(bool(iso_maps(H, all_graphs['r35_13.g6'][0], True)), 'Native cyclic H is bank member')
    au = iso_maps(H, H)
    affine = {tuple((a*x+b) % 13 for x in range(13)) for a in [1, 5, 8, 12] for b in range(13)}
    need(set(au) == affine, 'Aut(H)=C13:C4')
    need(all(sum(p[i] == i for i in range(13)) <= 1 for p in au if p != tuple(range(13))),
         'Two fixed points force identity on H')
    N = [1, 5, 8, 12]
    private = {a: next(x for x in range(1, 13) if x not in N and (H[x] & H[0]) == 1 << a) for a in N}
    double = {x: tuple(a for a in N if (H[x] >> a) & 1) for x in range(1, 13)
              if x not in N and (H[x] & H[0]).bit_count() == 2}
    arrows = []
    for a in N:
        opposite = next(b for b in N if b != a and (H[a] & H[b]).bit_count() == 1)
        x = next(x for x in double if (H[private[a]] >> x) & 1)
        target = next(b for b in double[x] if b != opposite)
        arrows.append((a, target))
    need(set(arrows) == {(1, 8), (8, 12), (12, 5), (5, 1)}, 'Intrinsic directed four-cycle')
    deletion_counts, addition_counts = Counter(), Counter()
    for u, v in combinations(range(13), 2):
        if (H[u] >> v) & 1:
            edited = list(H)
            edited[u] ^= 1 << v
            edited[v] ^= 1 << u
            deletion_counts[len(sets_of_size(edited, 5, 0))] += 1
        else:
            addition_counts[(H[u] & H[v]).bit_count()] += 1
    need(deletion_counts == {3: 26} and addition_counts == {1: 26, 2: 26}, 'Exact double-saturation witnesses')
    squares = {i*i % 17 for i in range(1, 17)}
    P = graph(17, lambda i, j: (i-j) % 17 in squares)
    need(bool(iso_maps(P, all_graphs['r44_17.g6'][0], True)), 'Paley17 is the separate (4,4) bank member')
    out['H13_structure'] = {'private': private, 'double': double, 'arrows': arrows,
                           'independent_sets': {k: len(sets_of_size(H, k, 0)) for k in range(1, 6)},
                           'edge_deletion_witness_histogram': dict(deletion_counts),
                           'nonedge_addition_witness_histogram': dict(addition_counts)}
    return out


def multiply(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 0b10011
    return out


def gf16_control():
    powers, x = [], 1
    for _ in range(15):
        powers.append(x)
        x = multiply(x, 2)
    need(x == 1 and len(set(powers)) == 15, 'Primitive GF16 element')
    color = {x: k % 3 for k, x in enumerate(powers)}
    frame = powers[::3]
    need(frame[0] ^ frame[1] ^ frame[2] ^ frame[3] ^ frame[4] == 0, 'Five-frame sum')
    for subset in combinations(frame, 4):
        span = {0}
        for a in subset:
            span |= {x ^ a for x in tuple(span)}
        need(len(span) == 16, 'Every four frame elements are independent')
    C = [[-1 if i == j else color[i ^ j] for j in range(16)] for i in range(16)]
    graphs = [graph(16, lambda i, j: C[i][j] == c) for c in range(3)]
    for A in graphs:
        need({r.bit_count() for r in A} == {5}, 'Each color is 5-regular')
        need(not sets_of_size(A, 3, 1), 'Each cube-coset color is triangle-free')
        need({(A[u] & A[v]).bit_count() for u, v in combinations(range(16), 2) if not ((A[u] >> v) & 1)} == {2},
             'Clebsch nonadjacent codegree two')
    aut0 = iso_maps(graphs[0], graphs[0])
    autC = [p for p in aut0 if all(C[i][j] == C[p[i]][p[j]] for i, j in combinations(range(16), 2))]
    need(len(aut0) == 1920 and len(autC) == 160, 'Forgetting two colors increases symmetry')
    for v in range(16):
        B = [u for u in range(16) if C[v][u] == 2]
        counts = [sum(C[u][w] == c for u, w in combinations(B, 2)) for c in range(3)]
        need(counts == [5, 5, 0], 'Edge-budget transfer positive control')
    return {'n': 16, 'color_classes_of_differences': [[x for x in range(1, 16) if color[x] == c] for c in range(3)],
            'single_color_aut_order': len(aut0), 'all_colors_aut_order': len(autC),
            'third_port_edge_counts': [5, 5, 0]}


def witness(name):
    data = json.loads((HERE/(name+'.json')).read_text())
    raw = (HERE/(name+'.cnf')).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == data['dimacs_sha256'], 'CNF hash')
    need(data['status'] == 'SAT', 'Stored existence witness')
    X = data['cross_colors']
    need(len(X) == 13 and all(len(row) == 13 and set(row) <= {0, 1, 2} for row in X), 'Matrix type')
    C = [[-1]*27 for _ in range(27)]
    def edge(i, j, c):
        C[i][j] = C[j][i] = c
    for i in range(13):
        edge(0, 1+i, 0)
        edge(0, 14+i, 1)
        for j in range(13):
            edge(1+i, 14+j, X[i][j])
    for i, j in combinations(range(13), 2):
        h = (i-j) % 13 in {1, 5, 8, 12}
        edge(1+i, 1+j, 1 if h else 2)
        edge(14+i, 14+j, 0 if h else 2)
    for c, k in [(0, 3), (1, 3), (2, 5)]:
        need(not any(all(C[i][j] == c for i, j in combinations(S, 2))
                     for S in combinations(range(27), k)), 'Native forbidden-clique check')
    red_rows = [row.count(0) for row in X]
    green_rows = [row.count(1) for row in X]
    red_cols = [sum(X[i][j] == 0 for i in range(13)) for j in range(13)]
    green_cols = [sum(X[i][j] == 1 for i in range(13)) for j in range(13)]
    need(set(red_rows) == set(green_cols) == {4}, 'Necessary directed port margins')
    need(sum(green_rows) == sum(red_cols) == 52, 'Mean four in remaining margins')
    if data['equitable_extra_assumption']:
        need(set(green_rows) == set(red_cols) == {4}, 'Extra equitability assumption')
    if data['forced_green_row_defect']:
        need(green_rows[0] <= 3 and len(set(green_rows)) > 1, 'Local axioms do not force uniform other margins')
    return {'n': 27, 'red_row_counts': red_rows, 'green_column_counts': green_cols,
            'green_row_counts': green_rows, 'red_column_counts': red_cols,
            'forbidden_counts': [0, 0, 0], 'scope': 'local witness only; not extendibility to 56'}


def identities():
    n, a, b, d = 56, 13, 29, 4
    e = Fraction(a*(b-2*a+1+3*d), 2)
    q = b*(b-1)//2-2*e
    need(e == 104 and q == 198 and n*q/3 == 3696, 'Third-port budget')
    need(Fraction(2*q, b) == Fraction(396, 29), 'Nonintegral blue edge-regular parameter')
    for t in range(10):
        table = [[0, 4, 8], [4, t, 9-t], [8, 9-t, 12+t]]
        need([sum(row) for row in table] == [12, 13, 29], 'Variable red-edge table')
    profiles = [(16, 17, 56), (17, 16, 56), (17, 17, 55)]
    enumerated = [(r, g, bl) for r in range(18) for g in range(18) for bl in range(57) if r+g+bl == 89]
    need(set(enumerated) == set(profiles), '90-vertex boundary has three degree profiles')
    parity_solutions = []
    for x in range(91):
        for y in range(91-x):
            z = 90-x-y
            if (1530-x) % 2 == (1530-y) % 2 == (5040-z) % 2 == 0:
                need(x % 2 == y % 2 == z % 2 == 0, 'All three profile populations even')
                parity_solutions.append([x, y, z])
    need(91*812 % 3 != 0, 'Boza cited upper-bound arithmetic')
    return {'56_degree_profile': [a, a, b], '29_port_edge_counts': [104, 104, 198],
            '29_port_deficit_totals': [24, 24, 68], 'blue_triangle_count': 3696,
            'blue_edge_codegree_mean': '396/29', 'red_edge_variable_range': [0, 9],
            '90_degree_profiles': profiles, 'even_profile_population_triples': len(parity_solutions),
            'symmetry_proof_numbers': {'local_aut': 52, 'faithful_neighborhoods_required': 1,
                                      'point_stabilizer_divides': 4,
                                      'full_color_preserving_aut_divides': 224}}


def hostile():
    path = HERE.parent/'regular_parent_symmetry_20261029'/'verify.py'
    spec = importlib.util.spec_from_file_location('inherited_hostile', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.hostile_checks()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    result = {'banks': banks(), 'GF16': gf16_control(),
              'balanced_local_witness': witness('balanced_ports'),
              'nonuniform_local_witness': witness('unequal_ports'),
              'identities': identities(), 'hostile_81_27': hostile()}
    result['native_gates'] = GATES
    data = json.dumps(result, indent=2, sort_keys=True)+'\n'
    target = HERE/'verification.json'
    if args.write:
        target.write_text(data, encoding='utf-8', newline='\n')
    else:
        need(target.read_text(encoding='utf-8') == data, 'Reproduction agrees with saved output')
    print(data, end='')


if __name__ == '__main__':
    main()
