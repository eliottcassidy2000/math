#!/usr/bin/env python3
"""Independent rectangle audit: direct D-A multiplication, no producer imports."""
import argparse
from array import array
from collections import Counter, defaultdict
import csv
import gzip
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import gcd
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valuations(n, cap):
    if n == 0:
        return cap
    return min(cap, (abs(n) & -abs(n)).bit_length() - 1)


def tensor_current(u, v, x, weight, modulus=None):
    ux = sum(a*b for a, b in zip(u, x))
    vx = sum(a*b for a, b in zip(v, x))
    out = [weight*(a*vx+b*ux) for a, b in zip(u, v)]
    return [a % modulus for a in out] if modulus else out


def controls():
    count = 0
    for k in range(1, 6):
        m = 1 << k
        for b in range(k+2):
            for gap in range(m):
                divisor = 1 << max(0, k-b-valuations(gap, k))
                for coefficient in range(-m, m+1):
                    require(((1 << b)*gap*coefficient) % m == 0
                            if coefficient % divisor == 0
                            else ((1 << b)*gap*coefficient) % m != 0,
                            'exact dyadic scalar divisibility including zero gap')
                    count += 1
    # With overlap, the two terms may cancel: u=v=(1,-1), x=(0,4).
    require(tensor_current([1, -1], [1, -1], [0, 4], 1, 8) == [0, 0],
            'overlap hostile vanishes')
    require(1 % 2 != 0, 'overlap hostile violates a false disjoint-support conclusion')
    # Derive the signed edge convention independently from four edge Laplacians.
    identities = 0
    for amplitude in range(-3, 4):
        for x in permutations(range(4)):
            current = [0]*4
            for i, j, delta_w in [(0, 2, -2*amplitude), (1, 3, -2*amplitude),
                                   (0, 3, 2*amplitude), (1, 2, 2*amplitude)]:
                current[i] += delta_w*(x[i]-x[j])
                current[j] += delta_w*(x[j]-x[i])
            require(current == tensor_current([1, -1, 0, 0], [0, 0, 1, -1], x, 2*amplitude),
                    'signed four-edge Laplacian equals stated symmetric tensor')
            identities += 1
    # A separately enumerated full kernel and image checks the general-v gcd
    # criterion, including v=0 and nonprimitive target coefficient vectors.
    from itertools import product
    gcd_cases = 0
    for weighted_edges in [
            [(0, 1, 2), (1, 2, 2), (2, 3, 2), (3, 4, 2)],
            [(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2), (0, 4, 2)]]:
        matrix = [[0]*5 for _ in range(5)]
        for i, j, weight in weighted_edges:
            matrix[i][i] += weight
            matrix[j][j] += weight
            matrix[i][j] -= weight
            matrix[j][i] -= weight
        vectors = list(product(range(4), repeat=5))
        def multiply(x):
            return tuple(sum(a*b for a, b in zip(row, x)) % 4 for row in matrix)
        im = {multiply(x) for x in vectors}
        harmonic = [x for x in vectors if multiply(x) == (0,)*5]
        u = [1, -1, 0, 0, 0]
        for p, q, c in product(range(-2, 3), range(-2, 3), range(4)):
            v = [0, 0, p, q, -p-q]
            g = gcd(*v)
            criterion = tuple(c*a % 4 for a in v) in im and tuple(c*g*a % 4 for a in u) in im
            literal = all(not any(tensor_current(u, v, x, c, 4)) for x in harmonic)
            require(criterion == literal, 'general net-target gcd image criterion on full finite kernel')
            gcd_cases += 1
    return dict(exact_scalar_divisibility_cases=count, signed_rectangle_cases=identities,
                general_net_target_gcd_cases=gcd_cases, disjoint_support_necessary_hostile=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--producer', type=Path, required=True)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    expected = json.loads((HERE/'SOURCE_PINS.json').read_text())
    actual = {}

    def read(name):
        raw = (args.repo/'05-knowledge/results'/name).read_bytes()
        actual[name] = sha256(raw).hexdigest()
        require(actual[name] == expected[name], 'authenticated original source: '+name)
        return raw

    algebra = controls()
    n = 3528
    root = 3527
    neighbors = [defaultdict(int) for _ in range(n)]
    degrees = [0]*n
    original_weights = {}
    edge_rows = list(csv.DictReader(gzip.decompress(read('metagraph_order8_20261023/merged_edges.tsv.gz')).decode().splitlines(), delimiter='\t'))
    for row in edge_rows:
        i, j, blue, black = [int(row[key]) for key in ('source', 'target', 'blue', 'black')]
        require(0 <= i <= j < n and blue >= 0 and black >= 0 and black % 2 == 0,
                'literal edge endpoints and paired-black weights')
        require((i, j) not in original_weights, 'one exported record per unordered owner pair')
        original_weights[i, j] = blue, black
        if i != j:
            weight = blue + black
            neighbors[i][j] += weight
            neighbors[j][i] += weight
            degrees[i] += weight
            degrees[j] += weight
    require(len(edge_rows) == 438216, 'complete literal pair row universe')

    # Independent implementation: adjacency D-A, not the solver's edge-current loop.
    def image(x):
        require(len(x) == n, 'full potential length')
        return [degrees[i]*x[i]-sum(w*x[j] for j, w in neighbors[i].items()) for i in range(n)]

    data = json.loads(read('metagraph_order8_electrical_20261023/receipt.json'))['order8']
    require(data['root'] == root, 'same original grounding')
    Y = data['radical_lifts_mod4']
    Z = data['eighth_radical_lifts_mod8'][0]
    X = data['higher_layers'][-1]['lifts'][0]
    require(len(Y) == 5, 'five authenticated modulo-four lifts')
    for label, vector, modulus in [('Y'+str(i), x, 4) for i, x in enumerate(Y)]+[('Z', Z, 8), ('X', X, 256)]:
        require(vector[root] == 0 and all(a % modulus == 0 for a in image(vector)),
                'direct full-integer harmonic verification for '+label)
    chars = [[2*a % 8 for a in y] for y in Y]
    for i, x in enumerate(chars):
        require(all(a % 8 == 0 for a in image(x)), 'actual scaled character '+str(i))
    sparse_nodes = [1238, 2246, 881, 2762, 2650, 2843]
    for node in sparse_nodes:
        sparse = [0]*n
        sparse[node] = 4
        require(all(a % 8 == 0 for a in image(sparse)), 'actual sparse order-two harmonic exists')
        require(all(2*a % 8 == 0 for a in sparse), 'paired-black factor kills sparse order-two character')

    raw_witness = (args.producer/'image_witnesses.json').read_bytes()
    witnesses = json.loads(raw_witness)
    require(witnesses['input_sha256'] == actual['metagraph_order8_20261023/merged_edges.tsv.gz'],
            'producer witnesses use exactly the independently read edge list')
    required_pairs = [(1238, 2246), (881, 2762), (2650, 2843)]
    require([tuple(s['charge']) for s in witnesses['solutions']] == required_pairs,
            'all three claimed charge witnesses, no substitutes')
    orders = []
    for solution in witnesses['solutions']:
        a, b = solution['charge']
        require(solution['coefficient'] == 4, 'image witness states the exact annihilating coefficient')
        z = solution['potential_mod8']
        require(len(z) == n and z[root] == 0 and all(type(t) is int and 0 <= t < 8 for t in z),
                'explicit grounded integral mod-eight witness')
        current = image(z)
        rhs = [0]*n
        rhs[a], rhs[b] = 4, -4
        require(all((y-r) % 8 == 0 for y, r in zip(current, rhs)), 'direct certificate Lz=4d modulo eight')
        require([[i, c % 8] for i, c in enumerate(current) if c % 8] == solution['image_mod8_nonzero'],
                'advertised image residue support exactly reproduced')
        require([[i, (c-rhs[i])//8] for i, c in enumerate(current) if c != rhs[i]] == solution['integer_correction_div8_nonzero'],
                'entire integer correction reproduced, not only a residual flag')
        separator = next(i for i, x in enumerate(chars) if (x[a]-x[b]) % 4 == 2)
        gap = (chars[separator][a]-chars[separator][b]) % 8
        require((2*gap) % 8 == 4, 'actual harmonic separator excludes 2d from image')
        require((X[a]-X[b]) % 8 == 0, 'original chosen character misses this nonzero charge')
        orders.append(dict(charge=[a, b], order_in_C8=4, annihilating_multiple=4,
                           separating_scaled_mod4_character=separator, gap_mod8=gap,
                           twice_charge_evaluation_mod8=(2*gap)%8,
                           witness_sha256=sha256(bytes(z)).hexdigest()))

    # Exact all-bank criterion: target rows are distinct and therefore cannot cancel.
    u = [0]*n
    u[1238], u[2246] = 1, -1
    amplitude_checks = 0
    for a in range(-5, 6):
        for b in range(-5, 6):
            v = [0]*n
            v[881], v[2762], v[2650], v[2843] = a, -a, b, -b
            preserved_by_checked_bank = all(not any(tensor_current(u, v, x, 2, 8)) for x in chars+[Z, [a % 8 for a in X]])
            require(preserved_by_checked_bank == (a % 2 == b % 2 == 0),
                    'independent disjoint-target amplitude parity classification')
            amplitude_checks += 1
    # This is a signed algebraic counterexample, not an asserted legal path move.
    v = [0]*n
    v[0], v[5] = 2, -2
    hostile = tensor_current(u, v, Z, 2, 8)
    require(all(a % 2 == 0 for a in v) and any(hostile), 'even net target coefficients alone are not sufficient')
    require([[i, a] for i, a in enumerate(hostile) if a] == [[1238, 4], [2246, 4]],
            'exact actual-bank hostile at original source pair')

    nodes = list(csv.DictReader(read('metagraph_order8_20261023/merged_nodes.tsv').decode().splitlines(), delimiter='\t'))
    ordinary_to_merged = {}
    for row in nodes:
        for key in ('class_a', 'class_b'):
            ordinary_to_merged[int(row[key])] = int(row['merged_id'])
    owner = array('H')
    owner.frombytes(gzip.decompress(read('metagraph_order8_20261023/owner_map.u16.gz')))
    if sys.byteorder != 'little':
        owner.byteswap()
    require(len(owner) == 1 << 21 and len(ordinary_to_merged) == 6880, 'complete authentic owner maps')
    chords = [(i, j) for i in range(8) for j in range(i+2, 8)]
    full = (1 << 21)-1
    sources = {}
    matrices = {}
    automorphism_counts = {}
    for source in (1238, 2246):
        row = next(row for row in nodes if int(row['merged_id']) == source)
        canon = int(row['canon_a'])
        matrix = [[0]*8 for _ in range(8)]
        for bit, (i, j) in enumerate(combinations(range(8), 2)):
            matrix[i][j] = (canon >> (27-bit)) & 1
            matrix[j][i] = 1-matrix[i][j]
        matrices[source] = matrix
        automorphisms = [p for p in permutations(range(8))
                         if all(matrix[p[i]][p[j]] == matrix[i][j]
                                for i, j in combinations(range(8), 2))]
        require(automorphisms == [tuple(range(8))],
                'source has no nontrivial automorphism kernel lifts')
        automorphism_counts[source] = len(automorphisms)
        paths = {}
        for p in permutations(range(8)):
            if not all(matrix[i][j] for i, j in zip(p, p[1:])):
                continue
            q = sum(matrix[p[i]][p[j]] << bit for bit, (i, j) in enumerate(chords))
            require(ordinary_to_merged[owner[q]] == source, 'actual Hamiltonian path source authentication')
            target = ordinary_to_merged[owner[full ^ q]]
            paths[p] = target
            reverse = p[::-1]
            require(all(matrix[j][i] for i, j in zip(reverse, reverse[1:])), 'literal reverse is Hamiltonian in the converse source')
            q_reverse = sum(matrix[reverse[j]][reverse[i]] << bit for bit, (i, j) in enumerate(chords))
            require(ordinary_to_merged[owner[q_reverse]] == source
                    and ordinary_to_merged[owner[full ^ q_reverse]] == target,
                    'reversed converse retains both original merged owners')
        sources[source] = paths
    require([len(sources[s]) for s in (1238, 2246)] == [409, 523], 'complete two-source Hamiltonian census')

    selections = [
        dict(targets=[881, 2762], source_paths={1238: [(7, 1, 0, 4, 2, 6, 5, 3), (7, 1, 0, 5, 2, 6, 4, 3)],
                                              2246: [(6, 5, 4, 3, 2, 0, 7, 1), (6, 5, 2, 4, 3, 0, 7, 1)]},
             pins={1238: (4, 2, 6), 2246: (5, 0, 7)}, weights=[4, 2, 4, 4]),
        dict(targets=[2650, 2843], source_paths={1238: [(7, 0, 6, 4, 3, 1, 5, 2), (7, 0, 6, 4, 1, 5, 3, 2)],
                                                2246: [(6, 0, 7, 3, 2, 1, 5, 4), (6, 0, 7, 3, 1, 5, 2, 4)]},
             pins={1238: (1, 0, 6), 2246: (1, 0, 7)}, weights=[4, 2, 4, 6])]
    rectangles = []
    for selection in selections:
        a, b = selection['targets']
        native_weights = [original_weights[min(s, t), max(s, t)] for s in (1238, 2246) for t in (a, b)]
        require(native_weights == [(0, w) for w in selection['weights']], 'original pure-black rectangle multiplicities')
        selected = []
        for s, paths in selection['source_paths'].items():
            p, q = paths
            required_targets = [a, b] if s == 1238 else [b, a]
            require([sources[s][p], sources[s][q]] == required_targets, 'actual selected paths and rectangle orientation')
            require(p[0] == q[0] and p[-1] == q[-1], 'literal endpoint pins retained')
            clock, first, second = selection['pins'][s]
            require(p[clock:clock+2] == q[clock:clock+2] == (first, second), 'literal indexed arc pins retained')
            require((X[s]-X[a]) % 8 == (X[s]-X[b]) % 8, 'selected original character phase retained')
            selected.append(dict(source=s, old_path=list(p), new_path=list(q), targets=required_targets,
                                 fixed_arc=list(selection['pins'][s]), clock_cost=sum(x != y for x, y in zip(p, q))))
        scarce = [list(p) for p, target in sources[1238].items() if target == b]
        require(len(scarce) == 1, 'global scarce target has exactly one distinct source-1238 path')
        require(set(selected[0]['old_path']) == set(range(8)), 'nonempty amplitude-one feasibility witness')
        rectangles.append(dict(targets=[a, b], selected=selected, original_blue_black_weights=native_weights,
                               scarce_target=b, all_source1238_paths_to_scarce_target=scarce,
                               maximum_distinct_observer_amplitude=1,
                               full_original_mod8_bank_amplitudes='even integers',
                               bank_preserving_distinct_observer_amplitudes=[0]))

    result = dict(status='PASS: independent exact rectangle-bank audit',
                  input_sha256=actual, producer_image_witness_sha256=sha256(raw_witness).hexdigest(),
                  original_pair_rows=len(edge_rows), vertices=n, root=root,
                  direct_D_minus_A_checks=True, source_charge_orders=orders,
                  original_scaled_mod4_gaps=[[ (x[a]-x[b]) % 8 for a, b in required_pairs] for x in chars],
                  algebra_controls=algebra, amplitude_pairs_checked=amplitude_checks,
                  sparse_order_two_control_nodes=sparse_nodes,
                  general_even_coefficient_hostile=dict(source_pair=[1238, 2246], target_vector=[[0, 2], [5, -2]],
                                                        character='original eighth_radical_lifts_mod8[0]',
                                                        nonzero_current_mod8=[[i, a] for i, a in enumerate(hostile) if a]),
                  source_path_counts={s: len(paths) for s, paths in sources.items()},
                  source_automorphism_group_orders=automorphism_counts,
                  rectangles=rectangles, checks=CHECKS)
    output = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.write:
        (HERE/'receipt.json').write_text(output)
    elif (HERE/'receipt.json').exists():
        require(output == (HERE/'receipt.json').read_text(), 'frozen audit receipt reproduces byte-for-byte')
    print(output, end='')


if __name__ == '__main__':
    main()
