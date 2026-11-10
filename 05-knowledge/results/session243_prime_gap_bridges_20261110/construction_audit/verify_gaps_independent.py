#!/usr/bin/env python3
"""Independent physical replay of the compact gap and endpoint-suffix claims."""
import json
import sys
from hashlib import sha256
from verify_independent import require, profile, graft, ordinal, strong, decode, transitive


def move_state(s, letter):
    D, E, U, d = s
    if letter == 'O': return D, 0, 2 * U, d
    if letter == 'Z': return D + U - E, d, 2 * U - d, d
    if letter == 'J': return D + U - E, D - E, D + U, D
    raise ValueError(letter)


def word_state(s, word):
    for letter in word: s = move_state(s, letter)
    return s


def word_graph(rows, root, word):
    for letter in word:
        old_n = len(rows)
        if letter == 'O':
            rows = ordinal(rows, (0,))
        else:
            rows = graft(rows, root, (0,))
            if letter == 'J': root = old_n
    return rows, root


def Uword(target):
    require(target >= 2, 'U must be at least two')
    reverse = []
    while target != 2:
        if target % 2:
            reverse.append('Z'); target = (target + 1) // 2
        else:
            reverse.append('O'); target //= 2
    return ''.join(reversed(reverse))


def compact_words(N):
    common = '' if N == 1 else Uword(N) + 'O'
    return common + 'Z', common + 'OZ'


def endpoint_source(N, extra=0):
    # Interior endpoint weights are 1,2,4,...,2^(n-3).
    n = max(4, N.bit_length() + 2) + extra
    rows = list(transitive(n))
    rows[0] ^= 1 << (n - 1)
    rows[n - 1] |= 1
    rows.append(0)
    for i in range(n):
        if 1 <= i <= n - 2 and N >> (i - 1) & 1:
            rows[i] |= 1 << n
        else:
            rows[n] |= 1 << i
    return tuple(rows), n


def main():
    if '--fail' in sys.argv: require(False, 'intentional failure')
    seed = (1, 0, 2, 1)
    report = {'status': 'PASS', 'compact_samples': [], 'endpoint_samples': []}
    receipt = sha256()
    for N in range(1, 10001):
        lower, upper = compact_words(N)
        ls, us = word_state(seed, lower), word_state(seed, upper)
        require(us[0] - ls[0] == 2 * N, 'compact exact gap')
        require(2 * N < ls[0] <= 3 * N and us[0] <= 5 * N, 'compact magnitude')
        require(max(len(lower), len(upper)) + 2 <= (N - 1).bit_length() + 4, 'compact order')
        receipt.update(repr((N, lower, upper, ls, us)).encode())
        if N <= 128:
            small, _ = word_graph((2, 0), 0, lower)
            large, _ = word_graph((2, 0), 0, upper)
            deleted = len(large) - 2
            keep = [v for v in range(len(large)) if v != deleted]
            induced = tuple(sum(1 << j for j, v in enumerate(keep) if large[u] >> v & 1) for u in keep)
            require(induced == small, 'literal induced one-vertex extension')
            for word, state in [(lower, ls), (upper, us)]:
                graph, root = word_graph((2, 0), 0, word)
                require(profile(graph, root) == state, 'literal compact state')
                require(strong(graph), 'literal compact strongness')
        if N in [1, 2, 3, 4, 5, 6, 7, 8, 9, 16, 17, 31, 32, 33, 127, 128, 10000]:
            report['compact_samples'].append({'N': N, 'words': [lower, upper],
                'counts': [ls[0], us[0]], 'orders': [len(lower) + 2, len(upper) + 2]})
    report['compact_exact_arithmetic_targets'] = 10000
    report['compact_literal_targets'] = 128
    report['compact_receipt_sha256'] = receipt.hexdigest()

    w1, w2 = 'JOJOZ', 'JOOOZ'
    basis = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    r1 = [word_state(s, w1)[0] for s in basis]
    r2 = [word_state(s, w2)[0] for s in basis]
    require(r1 == [9, -3, 9, 0] and r2 == [9, -1, 9, 0], 'suffix exact rows')
    report['suffix_rows'] = [r1, r2]
    tests = 0
    for n in range(1, 5):
        for k in range(1 << (n * (n - 1) // 2)):
            rows = decode(n, k)
            for root in range(n):
                before = profile(rows, root)
                values = []
                for word in [w1, w2]:
                    graph, newroot = word_graph(rows, root, word)
                    after = profile(graph, newroot)
                    require(after == word_state(before, word), 'literal suffix state')
                    values.append(after[0])
                    if strong(rows) and n >= 3:
                        require(strong(graph), 'literal suffix strongness')
                require(values[1] - values[0] == 2 * before[1], 'literal endpoint gap')
                tests += 1
    report['suffix_all_roots_sources_through_order4'] = tests
    for N in range(1, 13):
        previous = 0
        for extra in range(3):
            rows, root = endpoint_source(N, extra)
            state = profile(rows, root)
            require(state[1] == N and strong(rows), 'endpoint selector')
            values = []
            for word in [w1, w2]:
                graph, newroot = word_graph(rows, root, word)
                # Literal full-state enumeration for these controls (order <=14).
                after = profile(graph, newroot)
                require(after == word_state(state, word) and strong(graph), 'endpoint-pair physical')
                values.append(after[0])
            require(values[1] - values[0] == 2 * N, 'endpoint repeated gap')
            require(values[0] > previous, 'finite growing witness')
            previous = values[0]
            report['endpoint_samples'].append({'N': N, 'extra_base_vertices': extra,
                'source_order': len(rows), 'E': state[1], 'counts': values,
                'common_output_order': len(rows) + 5})
    print(json.dumps(report, sort_keys=True, indent=2))

if __name__ == '__main__': main()
