"""Regular-seven cells: exact kernel, marking and refinement controls.

Independent standard-library checker. All mathematical checks survive -O.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import hashlib
import json
import random

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise ValueError(message)


def code(A):
    return sum(A[i][j] << k for k, (i, j) in enumerate(combinations(range(len(A)), 2)))


def regular_sevens():
    A = [[0] * 7 for _ in range(7)]
    degrees = [0] * 7

    def visit(v):
        if v == 7:
            yield tuple(tuple(row) for row in A)
            return
        future = range(v + 1, 7)
        wanted = 3 - degrees[v]
        if not 0 <= wanted <= 6 - v:
            return
        for outgoing in combinations(future, wanted):
            S = set(outgoing)
            for w in future:
                A[v][w] = int(w in S)
                A[w][v] = 1 - A[v][w]
                degrees[w] += A[w][v]
            if all(degrees[w] <= 3 <= degrees[w] + 5 - v for w in future):
                yield from visit(v + 1)
            for w in future:
                degrees[w] -= A[w][v]

    return list(visit(0))


@lru_cache(None)
def autos(A):
    return tuple(p for p in permutations(range(len(A)))
                 if all(A[i][j] == A[p[i]][p[j]]
                        for i, j in combinations(range(len(A)), 2)))


def cycles(p):
    seen = set()
    result = []
    for v in range(len(p)):
        if v not in seen:
            cycle = []
            while v not in seen:
                seen.add(v)
                cycle.append(v)
                v = p[v]
            result.append(tuple(cycle))
    return result


def refine(A, cells):
    colour = [0] * len(A)
    for i, D in enumerate(cells):
        for v in D:
            colour[v] = i
    rounds = 0
    while True:
        names = sorted(set(colour))
        rows = [(colour[v], tuple(sum(A[v][u] for u in range(len(A))
                                    if colour[u] == c) for c in names))
                for v in range(len(A))]
        ids = {x: i for i, x in enumerate(sorted(set(rows)))}
        new = [ids[x] for x in rows]
        if len(set(new)) == len(set(colour)):
            return tuple(tuple(v for v in range(len(A)) if new[v] == c)
                         for c in sorted(set(new))), rounds
        rounds += 1
        colour = new


def marked_cells(n, marks):
    return [(v,) for v in marks] + [tuple(v for v in range(n) if v not in marks)]


def join(A, B, M):
    n = len(A)
    T = [[0] * (2 * n) for _ in range(2 * n)]
    for i, j in product(range(n), repeat=2):
        T[i][j] = A[i][j]
        T[n + i][n + j] = B[i][j]
        T[i][n + j] = M[i][j]
        T[n + j][i] = 1 - M[i][j]
    return T


def kernel(M, G, column=True):
    if column:
        return [p for p in G if all(row[j] == row[p[j]] for row in M for j in range(7))]
    return [p for p in G if all(M[i][j] == M[p[i]][j] for i, j in product(range(7), repeat=2))]


def block_bank():
    rng = random.Random(700310)
    bank = {}
    for k in range(1, 7):
        M = [[int((j - i) % 7 < k) for j in range(7)] for i in range(7)]
        for step in range(1200):
            a, b = rng.sample(range(7), 2)
            c, d = rng.sample(range(7), 2)
            if M[a][c] == M[b][d] and M[a][d] == M[b][c] and M[a][c] != M[a][d]:
                for i, j in ((a, c), (a, d), (b, c), (b, d)):
                    M[i][j] ^= 1
            if step % 12 == 0:
                U = tuple(tuple(row) for row in M)
                bank[U] = k
    return sorted(bank.items())


def all_partitions(items):
    if not items:
        yield ()
        return
    v, *tail = items
    for partition in all_partitions(tail):
        yield ((v,),) + partition
        for k in range(len(partition)):
            yield partition[:k] + ((v,) + partition[k],) + partition[k + 1:]


def equitable(A, P):
    return all(len({sum(A[v][w] for w in E) for v in D}) == 1 for D in P for E in P)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    regular = regular_sevens()
    need(len(regular) == 2640, 'regular-seven labelled count')
    need(len({code(A) for A in regular}) == 2640, 'duplicate regular tournament')
    by_code = {code(A): A for A in regular}
    remaining = set(by_code)
    reps, orbit_sizes = [], []
    while remaining:
        A = by_code[min(remaining)]
        orbit = {code(tuple(tuple(A[p[i]][p[j]] for j in range(7)) for i in range(7)))
                 for p in permutations(range(7))}
        need(orbit <= set(by_code), 'relabelled regular tournament absent')
        need(orbit <= remaining, 'orbits overlap')
        remaining -= orbit
        reps.append(A)
        orbit_sizes.append(len(orbit))
    need(len(reps) == 3, 'regular-seven class count')
    census = []
    invariant_rows = 0
    for A, orbit_size in zip(reps, orbit_sizes):
        G = autos(A)
        need(len(G) * orbit_size == 5040, 'orbit-stabilizer mismatch')
        shapes = Counter(tuple(sorted(map(len, cycles(p)))) for p in G)
        need(set(shapes) <= {(1,) * 7, (1, 3, 3), (7,)}, 'forbidden automorphism shape')
        one = []
        for v in range(7):
            P, _ = refine(A, marked_cells(7, (v,)))
            stabilizer = [p for p in G if p[v] == v]
            need((len(P) == 7) == (len(stabilizer) == 1), 'one-mark refinement/base discrepancy')
            one.append(len(P) == 7)
        for u, v in combinations(range(7), 2):
            P, _ = refine(A, marked_cells(7, (u, v)))
            need(len(P) == 7, 'two owner refinement failure')
            need(sum(p[u] == u and p[v] == v for p in G) == 1, 'two owner base failure')
        partition_shapes = Counter()
        for P in all_partitions(list(range(7))):
            if equitable(A, P):
                shape = tuple(sorted(map(len, P)))
                partition_shapes[shape] += 1
                need(shape in {(7,), (1, 3, 3), (1,) * 7}, 'unexpected equitable partition')
                if shape == (1, 3, 3):
                    cell_groups = [p for p in G if all(set(p[x] for x in D) == set(D) for D in P)]
                    need(len(cell_groups) == 3, 'proper nondiscrete partition lost C3')
        # Exhaust ALL matrices invariant under every nonidentity representative
        # automorphism, row-regular with any degree 0..7. Orbit-union words
        # give at most two choices in each fixed degree for these shapes.
        for p in G:
            if p == tuple(range(7)):
                continue
            orbits = cycles(p)
            words = []
            for chosen in product((0, 1), repeat=len(orbits)):
                row = [0] * 7
                for flag, orbit in zip(chosen, orbits):
                    for v in orbit:
                        row[v] = flag
                words.append(tuple(row))
            for k in range(8):
                choices = [row for row in words if sum(row) == k]
                for M in product(choices, repeat=7):
                    invariant_rows += 1
                    if all(sum(M[i][j] for i in range(7)) == k for j in range(7)):
                        need(k in (0, 7), 'nonuniform regular block has one-sided kernel')
        cyclic_out = sum(all(sum(A[x][y] for y in range(7) if A[v][y]) == 1
                             for x in range(7) if A[v][x]) for v in range(7))
        actual_paths = [P for P in permutations(range(7)) if all(A[u][v] for u, v in zip(P, P[1:]))]
        for P in actual_paths:
            survivors = list(G)
            chosen = []
            for v in P:
                if any(p[v] != v for p in survivors):
                    chosen.append(v)
                    survivors = [p for p in survivors if p[v] == v]
            need(len(chosen) == (2 if len(G) == 21 else 1), 'actual path greedy base budget')
            Q, _ = refine(A, marked_cells(7, chosen))
            need(len(Q) == 7, 'actual path owners do not refine regular-seven')
        census.append({'mask': code(A), 'labelled_orbit': orbit_size, 'automorphisms': len(G),
                       'automorphism_cycle_shapes': {str(k): v for k, v in sorted(shapes.items())},
                       'one_owner_discrete_choices': sum(one), 'two_owner_discrete_choices': 21,
                       'actual_hamiltonian_path_mark_controls': len(actual_paths),
                       'cyclic_outneighbourhoods': cyclic_out,
                       'equitable_partition_shapes': {str(k): v for k, v in sorted(partition_shapes.items())}})
    bank = block_bank()
    repeated_examples = []
    max_rounds = 0
    cases = 0
    for M, degree in bank:
        need(all(sum(row) == degree for row in M), 'row regularity lost')
        need(all(sum(M[i][j] for i in range(7)) == degree for j in range(7)), 'column regularity lost')
        repeated = len(set(M)) < 7 or len(set(zip(*M))) < 7
        if repeated and len(repeated_examples) < 3:
            repeated_examples.append({'degree': degree, 'rows': [''.join(map(str, row)) for row in M],
                                      'distinct_rows': len(set(M)), 'distinct_columns': len(set(zip(*M)))})
        for A, B in product(reps, repeat=2):
            need(len(kernel(M, autos(A), False)) == 1, 'row kernel')
            need(len(kernel(M, autos(B))) == 1, 'column kernel')
            T = join(A, B, M)
            P, rounds = refine(T, [(v,) for v in range(7)] + [tuple(range(7, 14))])
            need(len(P) == 14, 'discreteness failed to propagate')
            max_rounds = max(max_rounds, rounds)
            # Retain cell names and two actual owners in one cell.
            P, rounds = refine(T, [(0,), (1,), tuple(range(2, 7)), tuple(range(7, 14))])
            need(len(P) == 14, 'two root marks fail across attachment')
            # The converse uses the same owners and complement cross block.
            C = [[0 if i == j else 1 - T[i][j] for j in range(14)] for i in range(14)]
            P, _ = refine(C, [(0,), (1,), tuple(range(2, 7)), tuple(range(7, 14))])
            need(len(P) == 14, 'converse propagation failure')
            cases += 1
    need(bool(repeated_examples), 'distinct-row hypothesis not actually weakened')
    rectangle_cases = 0
    for b in range(20):
        M, degree = bank[b]
        # Repeat a regular block to make an equitable 21-by-7 attachment.
        # The source is the regular cyclic tournament of order21.
        rectangular = M + M + M
        for B in reps:
            T = [[int(i != j and (j - i) % 21 in range(1, 11)) for j in range(21)]
                 + [0] * 7 for i in range(21)]
            T += [[0] * 21 + list(row) for row in B]
            for i, j in product(range(21), range(7)):
                T[i][21 + j] = rectangular[i][j]
                T[21 + j][i] = 1 - rectangular[i][j]
            Q, _ = refine(T, [(v,) for v in range(21)] + [tuple(range(21, 28))])
            need(len(Q) == 28, 'rectangular source naming failed to propagate')
            need(len(kernel(rectangular, autos(B))) == 1, 'rectangular one-sided kernel')
            rectangle_cases += 1
    # Uniform attachments are genuinely nonfunctional and cannot propagate.
    U = tuple((1,) * 7 for _ in range(7))
    for A, B in product(reps, repeat=2):
        T = join(A, B, U)
        P, _ = refine(T, [(v,) for v in range(7)] + [tuple(range(7, 14))])
        need(len(P) == 8, 'uniform hostile unexpectedly resolved')
        need(len(kernel(U, autos(B))) == len(autos(B)), 'uniform kernel lost')
    out = {'status': 'FINITE-EXACT checks of proved scoped rules', 'labelled_regular_sevens': 2640,
           'census': census, 'all_invariant_row_regular_matrices_checked': invariant_rows,
           'sampled_regular_blocks': len(bank), 'distinct_internal_type_pair_controls': cases,
           'rectangular_twenty_one_to_seven_controls': rectangle_cases,
           'repeated_row_or_column_examples': repeated_examples,
           'maximum_neighbour_refinement_rounds_in_controls': max_rounds,
           'checks': CHECKS,
           'source_sha256_lf': hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest()}
    encoded = json.dumps(out, indent=2, sort_keys=True) + '\n'
    if args.write:
        (HERE / 'receipt.json').write_text(encoded, encoding='utf8')
    else:
        need(json.loads((HERE / 'receipt.json').read_text(encoding='utf8')) == out, 'receipt differs')
    print(json.dumps({'regular_sevens': 2640, 'classes': census, 'blocks': len(bank),
                      'attachment_cases': cases, 'checks': CHECKS}, sort_keys=True))


if __name__ == '__main__':
    main()
