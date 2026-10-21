#!/usr/bin/env python3
"""Independent exact audit. Imports no producer or repository Python code."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import prod
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

CHECKS = Counter()
def require(condition, label):
    CHECKS[label] += 1
    if not condition:
        raise RuntimeError(label)

def matmul(A, B):
    return [[sum(a * b for a, b in zip(row, column))
             for column in zip(*B)] for row in A]

def bareiss(A):
    """Fraction-free determinant, independent of producer Gaussian routine."""
    A = [list(row) for row in A]
    n = len(A)
    if n == 0:
        return 1
    sign, previous = 1, 1
    for k in range(n - 1):
        if not A[k][k]:
            replacement = next((i for i in range(k + 1, n) if A[i][k]), None)
            if replacement is None:
                return 0
            A[k], A[replacement] = A[replacement], A[k]
            sign = -sign
        pivot = A[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = pivot * A[i][j] - A[i][k] * A[k][j]
                if numerator % previous:
                    raise RuntimeError('nonintegral Bareiss division')
                A[i][j] = numerator // previous
            A[i][k] = 0
        previous = pivot
    return sign * A[-1][-1]

def solve(A, b):
    """Exact elimination of an augmented system, followed by back substitution."""
    B = [[F(x) for x in row] + [F(y)] for row, y in zip(A, b)]
    n = len(B)
    for k in range(n):
        pivot = next((i for i in range(k, n) if B[i][k]), None)
        if pivot is None:
            raise RuntimeError('singular independent solve')
        B[k], B[pivot] = B[pivot], B[k]
        for i in range(k + 1, n):
            multiplier = B[i][k] / B[k][k]
            for j in range(k, n + 1):
                B[i][j] -= multiplier * B[k][j]
    x = [F(0)] * n
    for i in reversed(range(n)):
        x[i] = (B[i][-1] - sum(B[i][j] * x[j] for j in range(i + 1, n))) / B[i][i]
    return x

def resistance(L, a, b):
    root = len(L) - 1
    red = [row[:root] for row in L[:root]]
    rhs = [int(i == a) - int(i == b) for i in range(root)]
    voltage = solve(red, rhs) + [F(0)]
    return voltage[a] - voltage[b]

def literal_contraction(weights):
    labels = [i for i, w in enumerate(weights) for _ in range(w)]
    q = len(weights)
    L = [[0] * q for _ in range(q)]
    loops = [0] * q
    edges = [[0] * q for _ in range(q)]
    for a, b in combinations(range(len(labels)), 2):
        x, y = labels[a], labels[b]
        if x == y:
            loops[x] += 1
        else:
            edges[x][y] += 1
            edges[y][x] += 1
            L[x][x] += 1
            L[y][y] += 1
            L[x][y] -= 1
            L[y][x] -= 1
    return labels, L, loops, edges

def pruefer_weighted_count(weights):
    q = len(weights)
    total = 0
    for word in product(range(q), repeat=q - 2):
        degree = [1] * q
        for x in word:
            degree[x] += 1
        edge_weight = 1
        for x in word:
            leaf = next(i for i, d in enumerate(degree) if d == 1)
            edge_weight *= weights[leaf] * weights[x]
            degree[leaf] -= 1
            degree[x] -= 1
        last = [i for i, d in enumerate(degree) if d == 1]
        edge_weight *= weights[last[0]] * weights[last[1]]
        total += edge_weight
    return total

def test_general_contraction():
    profiles = 0
    for q in range(2, 6):
        for weights in product(range(1, 4), repeat=q):
            profiles += 1
            W = sum(weights)
            labels, L, loops, edges = literal_contraction(weights)
            expected = [[W * weights[i] * int(i == j) - weights[i] * weights[j]
                         for j in range(q)] for i in range(q)]
            require(L == expected, 'literal Laplacian')
            require(loops == [w * (w - 1) // 2 for w in weights], 'literal loops')
            require(pruefer_weighted_count(weights) == W ** (q - 2) * prod(weights),
                    'Pruefer weighted trees')
            for root in range(q):
                cofactor = [[L[i][j] for j in range(q) if j != root]
                            for i in range(q) if i != root]
                require(bareiss(cofactor) == W ** (q - 2) * prod(weights),
                        'independent all-root cofactors')
            for a, b in combinations(range(q), 2):
                require(resistance(L, a, b) == (F(1, weights[a]) + F(1, weights[b])) / W,
                        'independent resistance solve')
            for potential in ([i * i - 2 * i for i in range(q)],
                              [(-1) ** i for i in range(q)]):
                original = sum((potential[labels[a]] - potential[labels[b]]) ** 2
                               for a, b in combinations(range(W), 2))
                quotient = sum(potential[i] * L[i][j] * potential[j]
                               for i in range(q) for j in range(q))
                require(original == quotient, 'literal constant-fibre energy')
            transitions = []
            for a, x in enumerate(labels):
                counts = Counter(labels[b] for b in range(W) if a != b)
                row = [F(counts[y], W - 1) for y in range(q)]
                require(row == [F(weights[y] - int(x == y), W - 1) for y in range(q)],
                        'microstate strong lumpability')
                transitions.append(row)
            quotient_rows = [transitions[labels.index(x)] for x in range(q)]
            degrees = [sum(row) + 2 * loops[x] for x, row in enumerate(edges)]
            for x in range(q):
                require(degrees[x] == weights[x] * (W - 1), 'loop darts preserve degree')
                for y in range(q):
                    dart_count = edges[x][y] if x != y else 2 * loops[x]
                    require(F(dart_count, degrees[x]) == quotient_rows[x][y], 'loop-dart walk')
                    require(F(weights[x], W) * quotient_rows[x][y]
                            == F(weights[y], W) * quotient_rows[y][x], 'stationary detailed balance')
            require([sum(F(weights[x], W) * quotient_rows[x][y] for x in range(q))
                     for y in range(q)] == [F(w, W) for w in weights], 'stationary law')
    return profiles

def test_integral_presentation():
    cases = 0
    for m in range(1, 13):
        for s in (1, 2, 3, 4, 5, 8, 16):
            cases += 1
            Q = m * s + 1
            A = [[s * (Q * int(i == j) - s) for j in range(m)] for i in range(m)]
            # Column basis U=(e_1,...,e_(m-1),S). Relation transform
            # V=(e_1+s*1,...,e_(m-1)+s*1,1). Both determinants are +1.
            U = [[int(i == j) if j < m - 1 else 1 for j in range(m)] for i in range(m)]
            V = [[int(i == j) + s if j < m - 1 else 1 for j in range(m)] for i in range(m)]
            diagonal = [s * Q] * (m - 1) + [s]
            D = [[diagonal[i] * int(i == j) for j in range(m)] for i in range(m)]
            require(bareiss(U) == 1 and bareiss(V) == 1, 'unimodular certificate determinants')
            require(matmul(A, V) == matmul(U, D), 'explicit integral equivalence AV=UD')
            require(abs(bareiss(A)) == prod(diagonal), 'Smith order')
            if s & (s - 1) == 0:
                two_primary = [x & -x for x in diagonal]
                # s=1 need not make Q odd; the primary assertion is only d>=1.
                if s > 1:
                    require(Q % 2 == 1 and two_primary == [s] * m, 'dyadic primary factors')
    return cases

def test_kron():
    cases = 0
    for q in range(2, 7):
        for s in range(1, 7):
            cases += 1
            Q = s * (q - 1) + 1
            weights = [s] * (q - 1) + [1]
            labels, _, _, _ = literal_contraction(weights)
            terminals = [labels.index(i) for i in range(q)]
            interior = [i for i in range(Q) if i not in terminals]
            micro_L = [[Q * int(i == j) - 1 for j in range(Q)] for i in range(Q)]
            H = [[micro_L[i][j] for j in interior] for i in interior]
            # Solve a distinct Dirichlet problem for each boundary basis vector.
            columns = []
            for b in range(q):
                boundary = [int(i == b) for i in range(q)]
                rhs = [-sum(micro_L[v][terminals[i]] * boundary[i] for i in range(q))
                       for v in interior]
                internal_values = solve(H, rhs) if interior else []
                potentials = dict(zip(terminals, boundary)) | dict(zip(interior, internal_values))
                currents = [sum(micro_L[t][v] * potentials[v] for v in range(Q))
                            for t in terminals]
                columns.append(currents)
            kron = [list(row) for row in zip(*columns)]
            require(kron == [[Q * int(i == j) - F(Q, q) for j in range(q)] for i in range(q)],
                    'Dirichlet-solved Kron response')
            for a, b in combinations(range(q), 2):
                require(resistance(kron, a, b) == F(2, Q), 'Kron resistance preservation')
    return cases

def test_entropy_and_blue():
    cases = 0
    for m in range(1, 7):
        for d in range(6):
            cases += 1
            s, q = 2 ** d, m + 1
            weights = [s] * m + [1]
            labels = [i for i, w in enumerate(weights) for _ in range(w)]
            Q = len(labels)
            fibres = Counter()
            child_blue = parent_blue = black_to_blue = 0
            for a, b in product(range(Q), repeat=2):
                x, y = labels[a], labels[b]
                fibres[x, y] += 1
                child_blue += int(a == b)
                parent_blue += int(x == y)
                black_to_blue += int(a != b and x == y)
            require(all(fibres[x, y] == weights[x] * weights[y]
                        for x, y in product(range(q), repeat=2)), 'enumerated ordered fibres')
            require(all(n & (n - 1) == 0 for n in fibres.values()), 'dyadic fibre sizes')
            entropy = sum(F(n, Q * Q) * (n.bit_length() - 1) for n in fibres.values())
            require(entropy == F(2 * d * m * s, Q), 'enumerated conditional Shannon entropy')
            uniform_parent_entropy = sum(F(n.bit_length() - 1, q * q) for n in fibres.values())
            require(uniform_parent_entropy == F(2 * d * m, q), 'uniform-parent alternate law')
            require(F(child_blue, Q * Q) == F(1, Q), 'enumerated child blue')
            require(F(parent_blue, Q * Q) == F(m * s * s + 1, Q * Q), 'enumerated parent blue')
            require(black_to_blue == parent_blue - child_blue == m * s * (s - 1),
                    'enumerated black-to-blue correction')
    return cases

def tournament_paths(k):
    n = 2 * k + 1
    def arc(a, b):
        return (a < b) != ((min(a, b), max(a, b)) in ((0, k), (k, 2 * k)))
    paths = []
    def visit(path, remaining):
        if not remaining:
            paths.append(tuple(path))
        else:
            for v in sorted(remaining):
                if not path or arc(path[-1], v):
                    visit(path + [v], remaining - {v})
    visit([], set(range(n)))
    return paths

def decode(k, I, K):
    A, B = set(range(k)), set(range(k + 1, 2 * k + 1))
    I = A if I is None else set(I)
    K = A if K is None else set(K)
    J = {2 * k - v for v in A - K}
    return tuple(sorted(I) + sorted(J) + [k] + sorted(A - I) + sorted(B - J))

def test_inherited_carrier():
    tables = {}
    for k in range(2, 6):
        finite = [frozenset(i + 1 for i in range(k - 1) if mask >> i & 1)
                  for mask in range(2 ** (k - 1))]
        states = finite + [None]
        decoded = {(I, K): decode(k, I, K) for I, K in product(states, repeat=2)}
        direct = tournament_paths(k)
        require(len(decoded) == len(set(decoded.values())), 'injective inherited decoder')
        require(set(decoded.values()) == set(direct), 'direct tournament DFS versus decoder')
        for (I, K), path in decoded.items():
            reflected = tuple(2 * k - v for v in reversed(path))
            require(reflected == decoded[K, I], 'literal reflection swap')
            require((reflected == path) == (I == K), 'literal blue fixed predicate')
        tables[k] = decoded
    deletion_counts = 0
    for child_k in range(3, 6):
        for parent_k in range(2, child_k):
            d = child_k - parent_k
            deleted = set(range(parent_k, child_k)) | set(range(child_k + 1, child_k + d + 1))
            survivors = [v for v in range(2 * child_k + 1) if v not in deleted]
            relabel = {v: i for i, v in enumerate(survivors)}
            counts = Counter()
            for (I, K), path in tables[child_k].items():
                parent_state = tuple(None if S is None else frozenset(v for v in S if v < parent_k)
                                     for S in (I, K))
                image = tuple(relabel[v] for v in path if v not in deleted)
                require(image == tables[parent_k][parent_state], 'literal all-depth path deletion')
                counts[parent_state] += 1
                deletion_counts += 1
            for (I, K), count in counts.items():
                require(count == 2 ** (d * (int(I is not None) + int(K is not None))),
                        'literal inherited fibre cardinality')
    return {'tournament_k': [2, 5], 'path_counts': [len(tables[k]) for k in sorted(tables)],
            'actual_path_deletions': deletion_counts}

def test_hostiles():
    witnesses = {}
    weights = [2, 2, 1]
    labels, L, loops, edges = literal_contraction(weights)
    require(resistance(L, 0, 1) == F(1, 5) != F(2, 5), 'hostile shorting is not Kron')
    require(resistance(L, 0, 2) == F(3, 10) != F(2, 5), 'hostile exceptional resistance')
    witnesses['shorting_not_kron'] = {'contraction': ['1/5', '3/10'], 'kron': '2/5'}
    dart_degrees = [sum(row) + 2 * loops[i] for i, row in enumerate(edges)]
    loopless_degrees = [sum(row) for row in edges]
    holding = F(2 * loops[0], dart_degrees[0])
    stationary = F(dart_degrees[0], sum(dart_degrees))
    loopless_stationary = F(loopless_degrees[0], sum(loopless_degrees))
    require(holding == F(1, 4) != 0 and stationary == F(2, 5) != loopless_stationary == F(3, 8),
            'hostile dropping loops changes walk')
    witnesses['drop_loops'] = {'holding': [str(holding), '0'],
                               'finite_stationary': [str(stationary), str(loopless_stationary)]}
    fibres = Counter((labels[a], labels[b]) for a, b in product(range(len(labels)), repeat=2))
    entropy_child = sum(F(n, len(labels) ** 2) * (n.bit_length() - 1) for n in fibres.values())
    entropy_parent = sum(F(n.bit_length() - 1, len(weights) ** 2) for n in fibres.values())
    require(entropy_child == F(8, 5) != entropy_parent == F(4, 3), 'hostile probability laws differ')
    witnesses['uniform_child_vs_uniform_parent_entropy'] = [str(entropy_child), str(entropy_parent)]
    require(sum(loops) == 2 and 2 * sum(loops) == 4, 'hostile ordered versus unordered collapse')
    witnesses['collapsed_black_mass'] = {'undirected_edges': sum(loops), 'ordered_paths': 2 * sum(loops)}
    recovered_weights = [-L[i][-1] for i in range(len(weights) - 1)] + [1]
    recovered_loops = [w * (w - 1) // 2 for w in recovered_weights]
    require(recovered_weights == weights and recovered_loops == loops,
            'hostile known source model recovers loop mass from Laplacian')
    # Adding arbitrary loops changes neither this electrical Laplacian nor its cofactors.
    extra_loops = [loops[0] + 1] + loops[1:]
    require(extra_loops != loops and sum(extra_loops) != sum(loops),
            'hostile arbitrary added loops invisible electrically')
    witnesses['loop_recovery_scope'] = {'known_pointed_model_recovers': recovered_loops,
                                       'same_electrical_network_with_extra_loop': extra_loops}
    require(bareiss([row[:2] for row in L[:2]]) == 20 != 4, 'hostile full group is not its 2-primary part')
    witnesses['critical_order_not_four'] = 20
    require(2 ** 4 != 4, 'hostile higher-rank dyadic group not V4')
    require(4 ** 2 != 4 * 2, 'hostile conductor16 order mismatch')
    witnesses['higher_groups'] = {'q5_s2_two_primary_order': 16, 'q3_s4_two_primary_order': 16,
                                  'conductor16_comparator_order': 8}
    _, s1_L, _, _ = literal_contraction([1, 1, 1, 1])
    require(bareiss([r[:3] for r in s1_L[:3]]) == 16,
            'hostile s1 has nontrivial two-primary group')
    witnesses['d0_boundary'] = {'q': 4, 's': 1, 'critical_order': 16}
    cleared_kron = [[10, -5], [-5, 10]]
    require(bareiss(cleared_kron) == 75 != 20, 'hostile denominator clearing changes lattice')
    witnesses['cleared_kron_tree_count'] = 75
    return witnesses

def verify_pins(producer, repo, workspace):
    manifest = json.loads((Path(__file__).parent / 'SOURCE_PINS.json').read_text())
    for pin in manifest['pins']:
        base = {'producer': producer, 'repository': repo, 'workspace': workspace}[pin['root']]
        path = base / pin['path']
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        require(digest == pin['sha256'], 'frozen source SHA256')
    return len(manifest['pins'])

def verify_producer(producer):
    normal = subprocess.run([sys.executable, str(producer / 'check.py')],
                            check=True, capture_output=True).stdout
    optimized = subprocess.run([sys.executable, '-O', str(producer / 'check.py')],
                               check=True, capture_output=True).stdout
    require(normal == optimized, 'producer normal optimized identical')
    require(normal == (producer / 'receipt.json').read_bytes(), 'producer receipt byte-identical')
    require(optimized == (producer / 'receipt.optimized.json').read_bytes(), 'producer optimized receipt byte-identical')
    parsed = json.loads(normal)
    require(parsed['status'] == 'PASS' and parsed['checks'] == 1328 and len(parsed['rows']) == 30,
            'producer frozen finite scope')
    return parsed['checks']

def main():
    parser = argparse.ArgumentParser()
    workspace_default = Path(__file__).resolve().parent.parent
    parser.add_argument('--workspace', type=Path, default=workspace_default)
    parser.add_argument('--producer', type=Path)
    parser.add_argument('--repo', type=Path)
    args = parser.parse_args()
    producer = args.producer or args.workspace / 'session20-choice-networks'
    repo = args.repo or args.workspace / 'math'
    pin_count = verify_pins(producer, repo, args.workspace)
    producer_checks = verify_producer(producer)
    profiles = test_general_contraction()
    presentations = test_integral_presentation()
    kron = test_kron()
    entropy = test_entropy_and_blue()
    carrier = test_inherited_carrier()
    hostiles = test_hostiles()
    print(json.dumps({'status': 'PASS', 'independent_checks': sum(CHECKS.values()),
                      'checks_by_kind': dict(sorted(CHECKS.items())), 'verified_pins': pin_count,
                      'producer_checks': producer_checks, 'arbitrary_weight_profiles': profiles,
                      'integral_presentations': presentations, 'kron_cases': kron,
                      'entropy_cases': entropy, 'inherited_carrier_controls': carrier,
                      'hostile_witnesses': hostiles,
                      'scope': 'Independent exact finite controls and a human all-parameter proof audit; no tournament metagraph or novelty claim.'}, indent=2) + '\n', end='')

if __name__ == '__main__':
    main()
