#!/usr/bin/env python3
"""Independent permutation-determinant and literal cycle audit.

Does not import either producer checker or the other audit checker.
"""
import argparse
import hashlib
import itertools as it
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(condition, reason):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(reason)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source(n, bits):
    free = [(i, j) for i in range(n) for j in range(i + 2, n)]
    arcs = {(i, i + 1) for i in range(n - 1)}
    arcs.update((i, j) if bits & (1 << k) else (j, i) for k, (i, j) in enumerate(free))
    return frozenset(arcs)


def complement(arcs, path):
    fixed = set(zip(path, path[1:]))
    need(fixed <= arcs, 'path is Hamiltonian')
    return frozenset((i, j) if (i, j) in fixed else (j, i) for i, j in arcs)


def targets(n, bits):
    T = source(n, bits)
    S = T | {(n, i) for i in range(n)} | {(i, n + 1) for i in range(n)} | {(n, n + 1)}
    K = complement(T, tuple(range(n)))
    H = complement(S, (n,) + tuple(range(n)) + (n + 1,))
    return K, H


def cycles(arcs, n, length):
    result = []
    for chosen in it.combinations(range(n), length):
        # One distinguished minimum vertex removes exactly the rotation factor.
        for rest in it.permutations(chosen[1:]):
            order = (chosen[0],) + rest
            if all((order[k], order[(k + 1) % length]) in arcs for k in range(length)):
                result.append(order)
    return result


def determinant(arcs, n):
    # Leibniz expansion of det(I-zA). Diagonal entries are 1. Every
    # nonzero off-diagonal factor is -z; each permutation is one monomial.
    coefficients = [0] * (n + 1)
    for perm in it.permutations(range(n)):
        moved = [(i, perm[i]) for i in range(n) if i != perm[i]]
        if not all(edge in arcs for edge in moved):
            continue
        inversions = sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        coefficients[len(moved)] += -1 if (inversions + len(moved)) % 2 else 1
    return coefficients


def walk_data(arcs, n, limit):
    adjacency = [[j for j in range(n) if (i, j) in arcs] for i in range(n)]
    distributions = [[int(i == j) for j in range(n)] for i in range(n)]
    result = []
    for k in range(limit + 1):
        result.append([row[:] for row in distributions])
        following = [[0] * n for _ in range(n)]
        for start in range(n):
            for current, count in enumerate(distributions[start]):
                for endpoint in adjacency[current]:
                    following[start][endpoint] += count
        distributions = following
    return result


def product(left, right, limit):
    return [sum(left[j] * right[k - j] for j in range(k + 1)
                if j < len(left) and k - j < len(right)) for k in range(limit + 1)]


def direct_cycle_control(n, bits):
    K, H = targets(n, bits)
    K2 = [[sum((i, v) in K and (v, j) in K for v in range(n)) for j in range(n)] for i in range(n)]
    a_term = sum(K2[0][j] for j in range(1, n))
    b_term = sum(K2[i][n - 1] for i in range(n - 1))
    both_term = int((0, n - 1) in K)
    cycle4 = cycles(H, n + 2, 4)
    boundary = [order for order in cycle4 if n in order or n + 1 in order]
    only_a = sum(n in order and n + 1 not in order for order in boundary)
    only_b = sum(n + 1 in order and n not in order for order in boundary)
    both = sum(n in order and n + 1 in order for order in boundary)
    need((only_a, only_b, both) == (a_term, b_term, both_term), 'literal four-cycle boundary split')
    need(len(boundary) == a_term + b_term + both_term, 'literal four-cycle formula')
    need(len(cycle4) - len(cycles(K, n, 4)) == len(boundary), 'literal induced four-cycle readout')
    return K, H


def formal_control(n, bits, K, H):
    limit = n + 8
    kw, hw = walk_data(K, n, limit), walk_data(H, n + 2, limit)
    aa = [sum(M[0][j] for j in range(1, n)) for M in kw]
    bb = [M[0][n - 1] for M in kw]
    cc = [sum(M[i][j] for i in range(n - 1) for j in range(1, n)) for M in kw]
    dd = [sum(M[i][n - 1] for i in range(n - 1)) for M in kw]
    need(aa[0] == bb[0] == dd[0] == 0, 'zero constant response terms')
    q = [int(k == 0) for k in range(limit + 1)]
    ad, bc = product(aa, dd, limit), product(bb, cc, limit)
    for k in range(2, limit + 1):
        q[k] -= aa[k - 2] + dd[k - 2]
    for k in range(3, limit + 1):
        q[k] -= bb[k - 3]
    for k in range(4, limit + 1):
        q[k] += ad[k - 4] - bc[k - 4]
    dk, dh = determinant(K, n), determinant(H, n + 2)
    need(product(dk, q, limit) == dh + [0] * (limit + 1 - len(dh)),
         'Schur response versus independent Leibniz determinant polynomial')
    # Formal division of -z Q' by Q, without using a trace recurrence.
    logarithmic_derivative = [0] * (limit + 1)
    for k in range(1, limit + 1):
        logarithmic_derivative[k] = -k * q[k] - sum(q[j] * logarithmic_derivative[k - j] for j in range(1, k + 1))
        literal = sum(hw[k][i][i] for i in range(n + 2)) - sum(kw[k][i][i] for i in range(n))
        need(logarithmic_derivative[k] == literal, 'formal logarithmic derivative versus literal rooted-walk count')
    boundary3 = len(cycles(H, n + 2, 3)) - len(cycles(K, n, 3))
    boundary4 = len(cycles(H, n + 2, 4)) - len(cycles(K, n, 4))
    need(logarithmic_derivative[3] == 3 * boundary3 and logarithmic_derivative[4] == 4 * boundary4,
         'rooted versus unrooted cycle factors')
    return {'n': n, 'bits': bits, 'detK': dk, 'detH': dh,
            'series_through_degree': limit, 'boundary_triangle_count': boundary3,
            'boundary_four_cycle_count': boundary4}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    pins = json.loads((HERE / 'GENERATING_INPUT_PINS.json').read_text())
    for row in pins['inputs']:
        p = HERE / row['snapshot']
        need(p.is_file() and digest(p) == row['sha256'] and p.stat().st_size == row['bytes'], 'source authentication')
    exhaustive_cycles = exhaustive_formal = 0
    witnesses = []
    for n in range(2, 7):
        for bits in range(1 << ((n - 1) * (n - 2) // 2)):
            K, H = direct_cycle_control(n, bits)
            exhaustive_cycles += 1
            if n <= 5:
                row = formal_control(n, bits, K, H)
                exhaustive_formal += 1
                if n == 5 and bits in (2, 3):
                    witnesses.append(row)
    extra_formal = []
    for n, bits in ((6, 0), (6, 277), (6, 1023)):
        K, H = targets(n, bits)
        extra_formal.append(formal_control(n, bits, K, H))
    need(exhaustive_cycles == 1099 and exhaustive_formal == 75, 'specified source universes')
    result = {'status': 'PASS', 'explicit_checks': CHECKS, 'python_optimization': sys.flags.optimize,
              'checker_sha256': digest(Path(__file__)), 'input_pins_sha256': digest(HERE / 'GENERATING_INPUT_PINS.json'),
              'exhaustive_direct_four_cycle_cases': exhaustive_cycles,
              'exhaustive_leibniz_series_cases': exhaustive_formal,
              'series_depth': 'n+8, two coefficients beyond producer controls',
              'witnesses': witnesses, 'extra_order6_formal_cases': extra_formal,
              'elapsed_seconds': round(time.monotonic() - started, 3),
              'scope': 'Directed simple four-cycles with chords allowed; rooted walks at larger lengths are not simple cycles.'}
    output = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
