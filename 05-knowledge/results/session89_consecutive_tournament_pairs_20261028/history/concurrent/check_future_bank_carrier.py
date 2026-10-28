#!/usr/bin/env python3
"""Exact standalone future-decision carrier checks; no producer imports."""
from collections import Counter
from itertools import permutations
import json

checks = 0


def require(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)


def rank_basis(rows):
    basis = {}
    for row in rows:
        insert(row, basis)
    return basis


def insert(row, basis):
    while row:
        pivot = row.bit_length() - 1
        if pivot not in basis:
            basis[pivot] = row
            return
        row ^= basis[pivot]


def rectangle_rows(n):
    return [1 ^ (1 << (i*n)) ^ (1 << j) ^ (1 << (i*n+j))
            for i in range(1, n) for j in range(1, n)]


def syndrome(word, rows):
    return sum(((word & row).bit_count() & 1) << k
               for k, row in enumerate(rows))


def coordinate_action(p):
    """Output-to-input coordinates of P_p and inversion offset of alpha_p."""
    n = len(p)
    inv = [p.index(i) for i in range(n)]
    sources = [None] * (n*n)
    offset = 0
    for i in range(n):
        sources[i*n+i] = inv[i]*n + inv[i]
        for j in range(i+1, n):
            lo, hi = sorted((inv[i], inv[j]))
            sources[i*n+j] = lo*n+hi
            sources[j*n+i] = hi*n+lo
            if inv[i] > inv[j]:
                offset ^= (1 << (i*n+j)) ^ (1 << (j*n+i))
    return sources, offset


def move(word, action):
    sources, offset = action
    return offset ^ sum(((word >> source) & 1) << target
                        for target, source in enumerate(sources))


def pullback(row, sources):
    answer = 0
    while row:
        bit = row & -row
        answer ^= 1 << sources[bit.bit_length()-1]
        row ^= bit
    return answer


def span(basis):
    values = {0}
    for vector in basis:
        values |= {word ^ vector for word in values}
    return values


def symmetric_cut_basis(n):
    return [(1 << (n*n))-1] + [
        sum(int((i == k) != (j == k)) << (i*n+j)
            for i in range(n) for j in range(n))
        for k in range(1, n)]


def normalize(word, n):
    """Unique first-row-zero representative modulo symmetric cuts + J."""
    bit = lambda i, j: (word >> (i*n+j)) & 1
    return sum((bit(i,j) ^ bit(0,i) ^ bit(0,j) ^ bit(0,0)) << (i*n+j)
               for i in range(n) for j in range(n))


def canonicalize(word, n):
    """Compact n(n-1)-bit carrier; the theorem's n>=3 regime only."""
    if n < 3:
        raise ValueError('Use the separately stated n=1,2 exceptional carriers.')
    return normalize(word, n) >> n


def carrier_toggle(state, cell, n):
    return canonicalize((state << n) ^ (1 << cell), n)


def carrier_relabel(state, permutation):
    n = len(permutation)
    return canonicalize(move(state << n, coordinate_action(permutation)), n)


def carrier_predicate(state, n):
    """The decision kernel consists exactly of constant remaining rows."""
    rowmask = (1 << n)-1
    return all(((state >> (r*n)) & rowmask) in (0, rowmask)
               for r in range(n-1))


def carrier_syndrome(state, n):
    answer = 0
    for r in range(n-1):
        row = (state >> (r*n)) & ((1 << n)-1)
        for j in range(1, n):
            answer |= (((row >> j) ^ row) & 1) << (r*(n-1)+j-1)
    return answer


def main():
    global checks
    checks = 0
    rank_results = []
    for n in range(1, 9):
        rows = rectangle_rows(n)
        require(len(rank_basis(rows)) == (n-1)**2, ('current syndrome rank', n))
        all_basis = {}
        permutation_count = 0
        for p in permutations(range(n)):
            sources, _ = coordinate_action(p)
            require(sorted(sources) == list(range(n*n)), ('P coordinate permutation', n, p))
            for row in rows:
                insert(pullback(row, sources), all_basis)
            permutation_count += 1
        expected_rank = 0 if n == 1 else 1 if n == 2 else n*n-n
        require(len(all_basis) == expected_rank, ('all-permutation constraint rank', n))
        if n >= 3:
            dbasis = symmetric_cut_basis(n)
            require(len(rank_basis(dbasis)) == n, ('symmetric cut dimension', n))
            require(all((row & vector).bit_count() % 2 == 0
                        for row in all_basis.values() for vector in dbasis),
                    ('all constraints annihilate claimed D', n))
        rank_results.append({'n': n, 'permutations': permutation_count,
                             'transported_rectangle_rows': permutation_count * len(rows),
                             'current_syndrome_rank': len(rows),
                             'future_constraint_rank': len(all_basis),
                             'common_kernel_dimension': n*n-len(all_basis)})

    # Full n=3 universe: independent truth contexts, linear syndrome profiles,
    # claimed kernel, explicit carrier coordinates, and deterministic transitions.
    n = 3
    rows = rectangle_rows(n)
    perms = list(permutations(range(n)))
    actions = [coordinate_action(p) for p in perms]
    cuts = {word for word in range(512) if syndrome(word, rows) == 0}
    require(len(cuts) == 32, 'n3 current code size')
    D = span(symmetric_cut_basis(n))
    require(len(D) == 8, 'n3 D size')
    future_basis = rank_basis(pullback(row, sources)
                              for sources, _ in actions for row in rows)
    require(len(future_basis) == 6, 'n3 independent future checks')
    basis_profiles = [syndrome(word, list(future_basis.values())) for word in range(512)]
    all_syndrome_profiles = [sum(syndrome(move(word, action), rows) << (4*g)
                                for g, action in enumerate(actions)) for word in range(512)]
    truth_profiles = []
    for word in range(512):
        profile = 0
        for g, action in enumerate(actions):
            moved = move(word, action)
            for toggle in range(512):
                if moved ^ toggle in cuts:
                    profile |= 1 << (512*g+toggle)
        truth_profiles.append(profile)
    require(Counter(Counter(basis_profiles).values()) == {8: 64}, '64 basis profiles, fibres8')
    require(Counter(Counter(all_syndrome_profiles).values()) == {8: 64}, '64 full syndrome profiles, fibres8')
    require(Counter(Counter(truth_profiles).values()) == {8: 64}, '64 full truth profiles, fibres8')
    require(Counter(Counter(syndrome(word, rows) for word in range(512)).values()) == {32: 16},
            '16 current syndromes, fibres32')
    for x in range(512):
        require(normalize(x, n) & 7 == 0, ('normal first row zero', x))
        require(x ^ normalize(x, n) in D, ('normalization difference in D', x))
        for y in range(512):
            expected = x ^ y in D
            require((truth_profiles[x] == truth_profiles[y]) == expected,
                    ('truth classes are D cosets', x, y))
            require((basis_profiles[x] == basis_profiles[y]) == expected,
                    ('basis classes are D cosets', x, y))
            require((all_syndrome_profiles[x] == all_syndrome_profiles[y]) == expected,
                    ('affine syndrome classes are D cosets', x, y))
            require((normalize(x, n) == normalize(y, n)) == expected,
                    ('explicit carrier classes are D cosets', x, y))
        representative = normalize(x, n)
        state = canonicalize(x, n)
        require(state == representative >> n, ('compact carrier', x))
        require(carrier_predicate(state, n) == (x in cuts), ('compact predicate', x))
        require(carrier_syndrome(state, n) == syndrome(x, rows), ('compact syndrome', x))
        for cell in range(9):
            require(normalize(x ^ (1 << cell), n) == normalize(representative ^ (1 << cell), n),
                    ('toggle descends', x, cell))
            require(carrier_toggle(state, cell, n) == canonicalize(x ^ (1 << cell), n),
                    ('compact toggle API', x, cell))
        for p, action in zip(perms, actions):
            require(normalize(move(x, action), n) == normalize(move(representative, action), n),
                    ('affine action descends', x, action[1]))
            require(carrier_relabel(state, p) == canonicalize(move(x, action), n),
                    ('compact relabel API', x, p))
        require((x in cuts) == (representative in cuts), ('decision descends', x))

    # In the compact 2x3 rectangle, each column is a chosen vector in F2^2.
    # The current-decision kernel is the diagonal three-column subspace.
    diagonal_columns = {sum(((u >> r) & 1) << (3*r+c)
                            for r in range(2) for c in range(3)) for u in range(4)}
    require({state for state in range(64) if carrier_predicate(state, 3)} == diagonal_columns,
            'current decision kernel is diagonal F2^2 in three columns')
    one_column_states = {sum(((u >> r) & 1) << (3*r+c) for r in range(2))
                         for c in range(3) for u in (1, 2, 3)}
    require(len(one_column_states) == 9, 'nine nonzero one-column states')
    require(len({carrier_syndrome(state, 3) for state in one_column_states}) == 9,
            'four-bit syndrome separates the nine selected states')
    require(all(carrier_syndrome(state, 3) != 0 for state in one_column_states),
            'selected states are outside current decision kernel')

    # Compare the exact nine-port map, in the prior bank's named cell coordinates.
    # The listed six-bit masks are followed there by an invertible GF64 basis map.
    port_columns = (1, 4, 5, 34, 2, 32, 16, 24, 8)

    def port_map(word):
        answer = 0
        for k, column in enumerate(port_columns):
            if word >> k & 1:
                answer ^= column
        return answer

    field_kernel = {word for word in range(512) if port_map(word) == 0}
    row_triples = {0b111 << (3*i) for i in range(3)}
    require(field_kernel == span(row_triples), 'actual F64 port kernel is row triples')
    require(len({port_map(word) for word in range(512)}) == 64, 'port map rank6')
    D_weights = dict(sorted(Counter(word.bit_count() for word in D).items()))
    field_weights = dict(sorted(Counter(word.bit_count() for word in field_kernel).items()))
    require(D_weights == {0: 1, 4: 3, 5: 3, 9: 1}, 'future kernel enumerator')
    require(field_weights == {0: 1, 3: 3, 6: 3, 9: 1}, 'field kernel enumerator')
    require(D_weights != field_weights, 'no cell permutation identifies the two kernels')

    print(json.dumps({'status': 'PASS', 'always_active_checks': checks,
                      'all_permutation_ranks': rank_results,
                      'n3': {'universe': 512, 'full_decision_contexts_per_word': 3072,
                             'distinct_truth_profiles': len(set(truth_profiles)),
                             'distinct_basis_profiles': len(set(basis_profiles)),
                             'profile_fibre_size': 8, 'compact_carrier_bits': 6,
                             'one_column_states_separated_by_current_syndrome': 9,
                             'future_kernel_weight_enumerator': D_weights,
                             'prior_field_port_kernel_weight_enumerator': field_weights},
                      'scope': 'Declared cell-toggle and natural-relabel decision contexts only; no pair-isoclass or canonical field identification.'}, indent=2))


if __name__ == "__main__":
    main()
