#!/usr/bin/env python3
"""Exact order12 repair: whole affine families, no tournament/HP census.

Standalone Python standard library. Affine matrices use bit0 as the
constant and bit1,... as independent forward-arc variables.
"""
from itertools import combinations, product
import json

N = 12
P = tuple(range(N))
PORTS = ((2, 7, 9), (6, 7, 9))
OLD = (
    (0, 1, 10, 11, 2, 3, 4, 5, 6, 7, 8, 9),
    (1, 10, 0, 11, 2, 3, 4, 5, 6, 7, 8, 9),
    (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 11),
)
NEW = (
    (0, 1, 2, 3, 4, 9, 10, 5, 6, 7, 8, 11),
    (0, 1, 2, 3, 4, 9, 10, 5, 7, 8, 11, 6),
    (0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 6),
)
GUARDS = ((4, 9), (10, 5), (7, 5))
HOSTILE_ROWS = (
    (1402, 3748, 1945, 2674, 1698, 2372, 3478, 2857, 2586, 1121, 2472, 533),
    (1342, 3748, 2008, 2674, 1698, 836, 3475, 2857, 2586, 3137, 2472, 53),
)
FALLBACK = {
    'open': (
        (0, 1, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7),
        (1, 8, 0, 9, 10, 11, 2, 3, 4, 5, 6, 7),
        (1, 10, 11, 2, 3, 4, 5, 6, 7, 8, 0, 9),
    ),
    'closed_zero': (
        (8, 10, 11, 2, 3, 4, 5, 6, 7, 1, 9, 0),
        (9, 8, 10, 11, 2, 3, 4, 5, 6, 7, 1, 0),
        (0, 9, 1, 10, 11, 2, 3, 4, 5, 6, 7, 8),
    ),
    'closed_one': (
        (8, 1, 9, 0, 10, 11, 2, 3, 4, 5, 6, 7),
        (9, 8, 1, 0, 10, 11, 2, 3, 4, 5, 6, 7),
        (0, 9, 1, 10, 11, 2, 3, 4, 5, 6, 7, 8),
    ),
}


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def edges(q):
    return {tuple(sorted((u, v))) for u, v in zip(q, q[1:])}


def reversal(q):
    return dict(zip(q, q[::-1]))


def eliminate(equations):
    basis = {}
    for e in equations:
        while e > 1:
            pivot = e.bit_length()-1
            if pivot in basis:
                e ^= basis[pivot]
            else:
                basis[pivot] = e
                break
        if e == 1:
            return None
    return basis


def reduce_form(form, basis):
    answer = form & 1
    form &= ~1
    while form:
        pivot = form.bit_length()-1
        if pivot in basis:
            form ^= basis[pivot]
            answer ^= form & 1
            form &= ~1
        else:
            bit = 1 << pivot
            answer ^= bit
            form ^= bit
    return answer


def raw_matrix():
    A = [[0]*N for _ in P]
    for i, (u, v) in enumerate(combinations(P, 2), 1):
        A[u][v] = 1 << i
        A[v][u] = (1 << i) ^ 1
    return A


def blue_equations(A, q):
    r = reversal(q)
    return [A[u][v] ^ A[r[v]][r[u]] for u, v in combinations(P, 2)]


def k_move(A, q):
    E = edges(q)
    return [[0 if u == v else A[u][v] ^ int(tuple(sorted((u, v))) not in E)
             for v in P] for u in P]


def specialize(A, equations):
    basis = eliminate(equations)
    need(basis is not None, 'source family must be nonempty')
    return [[0 if u == v else reduce_form(A[u][v], basis)
             for v in P] for u in P]


def source(ports):
    A = raw_matrix()
    equations = [A[i][i+1] ^ 1 for i in range(N-1)]
    equations += [A[11][0] ^ 1] + blue_equations(A, P)
    equations += [A[0][v] ^ int(v not in ports) for v in range(1, 11)]
    equations += [A[1][10] ^ 1, A[2][10] ^ 1]
    return specialize(A, equations)


def variables(A):
    mask = 0
    for row in A:
        for form in row:
            mask |= form & ~1
    return [i for i in range(1, mask.bit_length()) if mask >> i & 1]


def symbolic_replay(A, words, target_word=P):
    initial = A
    xor = set()
    black = []
    states = []
    for q in words:
        need(sorted(q) == list(P), 'not Hamiltonian word')
        need(all(A[u][v] == 1 for u, v in zip(q, q[1:])),
             ('nonuniform join', q))
        states.append(A)
        black.append(eliminate(blue_equations(A, q)) is None)
        xor ^= edges(q)
        A = k_move(A, q)
    need(xor == edges(target_word), 'path-edge XOR')
    need(A == k_move(initial, target_word), 'literal affine target')
    return black, states


def literal_model(A, assignment=0):
    return tuple(sum(((A[u][v] & 1) ^
                      ((A[u][v] & assignment).bit_count() & 1)) << v
                     for v in P if u != v) for u in P)


def literal_source(rows, ports):
    bit = lambda u, v: (rows[u] >> v) & 1
    need(all(not bit(u, u) for u in P), 'source loops')
    need(all(bit(u, v) ^ bit(v, u) == 1
             for u, v in combinations(P, 2)), 'source tournament')
    need(all(bit(i, i+1) for i in range(N-1)) and bit(11, 0),
         'source closed marked path')
    need(all(bit(u, v) == bit(11-v, 11-u)
             for u, v in combinations(P, 2)), 'source P anti')
    need({v for v in range(1, 11) if bit(v, 0)} == set(ports),
         'source exact port set')
    need(rows[0].bit_count() == 7 and rows[11].bit_count() == 4,
         'source endpoint scores')
    need(bit(1, 10) and bit(2, 10), 'source two branch guards')


def literal_replay(rows, words, target_word=P):
    def bit(R, u, v):
        return (R[u] >> v) & 1
    def blue(R, q):
        r = reversal(q)
        return all(bit(R, u, v) == bit(R, r[v], r[u])
                   for u, v in combinations(P, 2))
    def move(R, q):
        E = edges(q)
        return tuple(sum((bit(R, u, v) ^
                          int(tuple(sorted((u, v))) not in E)) << v
                         for v in P if u != v) for u in P)
    target = move(rows, target_word)
    flags = []
    for q in words:
        need(sorted(q) == list(P) and
             all(bit(rows, u, v) for u, v in zip(q, q[1:])),
             'literal word legality')
        flags.append(blue(rows, q))
        rows = move(rows, q)
    need(rows == target, 'literal target')
    return flags


def main():
    X = {(4, 5), (8, 9), (10, 11), (4, 9), (5, 10), (8, 11)}
    Y = {(5, 6), (6, 7), (5, 7), (6, 11)}
    need(edges(NEW[0]) == edges(P) ^ X, 'first edge difference')
    need(edges(NEW[1]) == edges(P) ^ X ^ Y, 'second edge difference')
    need(edges(NEW[2]) == edges(P) ^ Y, 'third edge difference')
    rH = reversal(OLD[0])
    g = {v: rH[11-v] for v in P}
    expected_cycles = ((0, 6, 2), (1, 7, 3), (4, 10, 8), (5, 11, 9))
    expected_g = {cycle[i]: cycle[(i+1) % 3]
                  for cycle in expected_cycles for i in range(3)}
    need(g == expected_g, 'literal reflection-product cycles')

    reports = []
    total_models = 0
    for ports, stored in zip(PORTS, HOSTILE_ROWS):
        A = source(ports)
        need(len(variables(A)) == 18, 'full source dimension')
        old_flags, _ = symbolic_replay(A, OLD)
        need(old_flags == [False, True, True], 'old exact blackness scope')

        G = specialize(A, [A[u][v] ^ 1 for u, v in GUARDS])
        need(len(variables(G)) == 15, 'guard family dimension')
        new_flags, states = symbolic_replay(G, NEW)
        need(new_flags == [True]*3, 'guarded full-affine blackness')
        need(reversal(NEW[0])[9] == 10 and reversal(NEW[0])[0] == 11,
             'first reflected pair')
        need(states[0][0][9] == 0 and states[0][10][11] == 1,
             'first black witness')
        need(reversal(NEW[1])[4] == 5 and reversal(NEW[1])[0] == 6,
             'second reflected pair')
        need(states[1][0][4] == 0 and states[1][5][6] == 1,
             'second black witness')
        B = k_move(G, P)
        need(reversal(NEW[2])[10] == 2 and reversal(NEW[2])[1] == 11,
             'third reflected pair')
        need(B[1][10] == 0 and B[2][11] == 1, 'third black witness')

        H = specialize(A, blue_equations(A, OLD[0]))
        need(len(variables(H)) == 3, 'blue first-word dimension')
        need(all(H[u][v] == H[g[u]][g[v]]
                 for u, v in combinations(P, 2)), 'source automorphism')
        need(all(H[u][v] == 1 for u, v in GUARDS),
             'failed-word symmetry forces all new joins')
        need(H[8][11] == H[3][11] == H[6][11] == 1,
             'fixed outgoing-port exclusions')
        need(symbolic_replay(H, NEW)[0] == [True]*3, 'whole blue family')

        literal_count = 0
        for values in product((0, 1), repeat=3):
            assignment = sum(value << var for value, var
                             in zip(values, variables(H)))
            rows = literal_model(H, assignment)
            literal_source(rows, ports)
            need(literal_replay(rows, OLD) == [True, False, False],
                 'complete old-blue family')
            need(literal_replay(rows, NEW) == [False]*3, 'literal new repair')
            literal_count += 1
        literal_source(stored, ports)
        need(literal_replay(stored, OLD) == [True, False, False],
             'inherited hostile retained')
        need(literal_replay(stored, NEW) == [False]*3, 'inherited hostile repaired')
        total_models += literal_count

        fallback = []
        for kind, guards in (
            ('open', [(1, 8, 1)]),
            ('closed_zero', [(1, 8, 0), (10, 8, 0)]),
            ('closed_one', [(1, 8, 0), (10, 8, 1)]),
        ):
            F = specialize(H, [H[u][v] ^ value for u, v, value in guards])
            replacement = FALLBACK[kind]
            need(symbolic_replay(F, replacement, OLD[0])[0] == [True]*3,
                 'five-step replacement prefix')
            five = replacement + OLD[1:]
            need(symbolic_replay(F, five)[0] == [True]*5,
                 'independent five-black fallback')
            fallback.append({'branch': kind, 'dimension': len(variables(F))})

        # Each guard is necessary for this particular new ledger.
        omitted = []
        for index, missing in enumerate(GUARDS):
            hostile = specialize(A, [A[u][v] ^ int((u, v) != missing)
                                     for u, v in GUARDS])
            stage = hostile if index < 2 else k_move(hostile, NEW[0])
            pair = missing if index < 2 else (5, 7)
            need(stage[pair[0]][pair[1]] == 0, 'actual illegal guarded join')
            rows = literal_model(hostile)
            literal_source(rows, ports)
            need(literal_replay(rows, OLD) == [False]*3,
                 'old repair handles this complementary source')
            omitted.append({'omitted_guard': missing, 'source_rows': rows,
                            'illegal_word': 1 if index < 2 else 2,
                            'illegal_join': pair})
        reports.append({
            'ports': ports, 'full_dimension': 18, 'new_guard_dimension': 15,
            'blue_old_first_dimension': 3, 'literal_blue_sources': literal_count,
            'old_all_completion_black': old_flags,
            'new_guard_all_completion_black': new_flags,
            'inherited_hostile_rows': stored,
            'omitted_guard_hostiles': omitted,
            'five_black_fallback_cells': fallback,
        })
    print(json.dumps({
        'status': 'PASS', 'source_cells': reports,
        'old_first_reflection_product_cycles': expected_cycles,
        'new_words': NEW, 'new_guards': GUARDS,
        'complete_literal_old_blue_sources': total_models,
        'method': 'standalone full-affine Gaussian proof gates plus all16 literal symmetric sources',
        'scope': 'both complete order12 failed base cells have three-black repairs; no quotient or clone contraction',
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
