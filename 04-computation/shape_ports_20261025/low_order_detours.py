#!/usr/bin/env python3
"""Independent exhaustive endpoint-three audit; no repository imports.

Universe: every labelled tournament on [n], 1 <= n <= 6. Retain exactly
standard directed Hamilton words P with z->a, blue reversal and endpoint
scores (n-4,3). Search black first/last paths, then replay with directed sets.
"""
from itertools import combinations, permutations
import argparse
import hashlib
import json
from pathlib import Path


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def audit_order(n):
    pairs = list(combinations(range(n), 2))
    index = {pair: i for i, pair in enumerate(pairs)}
    full = (1 << len(pairs)) - 1

    def arc(t, a, b):
        return bool((t >> index[tuple(sorted((a, b)))]) & 1) == (a < b)

    def edges(word):
        return sum(1 << index[tuple(sorted(e))] for e in zip(word, word[1:]))

    def directed(t, word):
        return all(arc(t, a, b) for a, b in zip(word, word[1:]))

    def blue(t, word):
        reverse = {word[i]: word[n - 1 - i] for i in range(n)}
        return all(arc(t, a, b) != arc(t, reverse[a], reverse[b]) for a, b in pairs)

    def rows(t):
        return [sum(1 << b for b in range(n) if a != b and arc(t, a, b)) for a in range(n)]

    def path_list(t, only_black):
        return [(p, edges(p)) for p in permutations(range(n))
                if directed(t, p) and (not only_black or not blue(t, p))]

    def words_of_edge_set(mask):
        adjacency = [[] for _ in range(n)]
        if mask.bit_count() != n - 1:
            return []
        for i, (a, b) in enumerate(pairs):
            if mask >> i & 1:
                adjacency[a].append(b)
                adjacency[b].append(a)
        ends = [i for i in range(n) if len(adjacency[i]) == 1]
        if len(ends) != 2 or any(not a or len(a) > 2 for a in adjacency):
            return []
        word = [ends[0]]
        while len(word) < n:
            options = [x for x in adjacency[word[-1]] if len(word) == 1 or x != word[-2]]
            if len(options) != 1 or options[0] in word:
                return []
            word.append(options[0])
        return [tuple(word), tuple(reversed(word))]

    def independent_replay(source_rows, words, expected_rows):
        arcs = {(a, b) for a in range(n) for b in range(n) if source_rows[a] >> b & 1}
        witnesses = []
        for word in words:
            check(sorted(word) == list(range(n)), 'not a permutation')
            check(all((a, b) in arcs for a, b in zip(word, word[1:])), 'illegal path')
            reverse = dict(zip(word, reversed(word)))
            hostile = next(((a, b) for a, b in sorted(arcs)
                            if (reverse[a], reverse[b]) in arcs), None)
            check(hostile is not None, 'blue step in alleged black certificate')
            witnesses.append(hostile)
            path_pairs = {frozenset((a, b)) for a, b in zip(word, word[1:])}
            arcs = {(a, b) if frozenset((a, b)) in path_pairs else (b, a) for a, b in arcs}
        actual = [sum(1 << b for a, b in arcs if a == v) for v in range(n)]
        check(actual == expected_rows, 'nonliteral target')
        return witnesses

    base = tuple(range(n))
    base_edges = edges(base)
    retained = []
    for source in range(1 << len(pairs)):
        if n < 2 or not directed(source, base) or not arc(source, n - 1, 0) or not blue(source, base):
            continue
        if sum(arc(source, 0, b) for b in range(1, n)) != n - 4:
            continue
        if sum(arc(source, n - 1, b) for b in range(n - 1)) != 3:
            continue
        target = source ^ full ^ base_edges
        first, last = path_list(source, True), path_list(target, True)
        certificate = None
        for q1, e1 in first:
            middle_source = source ^ full ^ e1
            for q3, e3 in last:
                for q2 in words_of_edge_set(base_edges ^ e1 ^ e3):
                    if directed(middle_source, q2) and not blue(middle_source, q2):
                        certificate = (q1, q2, q3)
                        break
                if certificate:
                    break
            if certificate:
                break
        source_rows, target_rows = rows(source), rows(target)
        witnesses = independent_replay(source_rows, certificate, target_rows) if certificate else None
        check(certificate is not None or not last,
              'unresolved search: not certified impossible by absence of target black paths')
        retained.append({'source_rows': source_rows, 'target_rows': target_rows,
                         'black_source_paths': len(first), 'black_target_paths': len(last),
                         'target_all_paths': len(path_list(target, False)),
                         'certificate': certificate, 'anti_failure_witnesses': witnesses})
    return {'n': n, 'labelled_universe': 1 << len(pairs), 'retained': len(retained),
            'three_black_detours': sum(r['certificate'] is not None for r in retained),
            'rows': retained}


def odd_overlap_witness():
    n = 7
    incoming = {2, 4}
    reflected = {n - 1 - v for v in incoming}
    arcs = {(a, b) for a in range(1, n - 1) for b in range(a + 1, n - 1)}
    arcs.add((n - 1, 0))
    for v in range(1, n - 1):
        arcs.add((v, 0) if v in incoming else (0, v))
        arcs.add((n - 1, v) if v in reflected else (v, n - 1))
    check(incoming == reflected, 'not overlap two')
    check(all((i, i + 1) in arcs for i in range(n - 1)), 'not a path')
    check(all((n - 1 - b, n - 1 - a) in arcs for a, b in arcs), 'not blue')
    check(sum(a == 0 for a, _ in arcs) == 3 and sum(a == n - 1 for a, _ in arcs) == 3,
          'wrong endpoint layer')
    return {'n': n, 'I': sorted(incoming), 'J': sorted(reflected),
            'rows': [sum(1 << b for a, b in arcs if a == v) for v in range(n)]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    orders = [audit_order(n) for n in range(1, 7)]
    check([r['retained'] for r in orders] == [0, 0, 0, 0, 2, 12], 'universe mismatch')
    check([r['three_black_detours'] for r in orders] == [0, 0, 0, 0, 0, 11], 'detour count mismatch')
    output = {'status': 'FINITE-EXACT complete labelled universes n<=6', 'orders': orders,
              'odd_overlap_two': odd_overlap_witness()}
    payload = json.dumps(output, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(payload)
    print(json.dumps({'retained': [r['retained'] for r in orders],
                      'detours': [r['three_black_detours'] for r in orders],
                      'sha256': hashlib.sha256(payload.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
