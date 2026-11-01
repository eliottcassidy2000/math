#!/usr/bin/env python3
"""All-labelled n<=6 audit of fixed Hamiltonian-count insertion increments.

No inherited atlas or third-party package is imported. Checks use explicit
exceptions and remain live with Python -O. Counts refer to directed vertex
words, not automorphism orbits.
"""
from collections import Counter
from hashlib import sha256
from itertools import permutations
import json
from math import factorial
from pathlib import Path


def require(ok, message, data=None):
    if not ok:
        raise RuntimeError((message, data))


def from_code(n, code):
    rows = [0] * n
    bit = 0
    for u in range(n):
        for v in range(u + 1, n):
            if (code >> bit) & 1:
                rows[u] |= 1 << v
            else:
                rows[v] |= 1 << u
            bit += 1
    return rows


def path_counts(rows):
    n = len(rows)
    total = 1 << n
    incoming = [sum(1 << u for u in range(n) if rows[u] >> v & 1)
                for v in range(n)]
    table = [[0] * n for _ in range(total)]
    H = [1] + [0] * (total - 1)
    for mask in range(1, total):
        bits = mask
        while bits:
            b = bits & -bits
            bits -= b
            v = b.bit_length() - 1
            old = mask ^ b
            value = int(old == 0)
            pred = old & incoming[v]
            while pred:
                a = pred & -pred
                pred -= a
                value += table[old][a.bit_length() - 1]
            table[mask][v] = value
            H[mask] += value
    return H


def components(rows):
    n = len(rows)
    reach = [(1 << u) | rows[u] for u in range(n)]
    for mid in range(n):
        for u in range(n):
            if reach[u] >> mid & 1:
                reach[u] |= reach[mid]
    unseen = (1 << n) - 1
    answer = []
    while unseen:
        b = unseen & -unseen
        u = b.bit_length() - 1
        part = sum(1 << v for v in range(n)
                   if reach[u] >> v & 1 and reach[v] >> u & 1)
        answer.append(part)
        unseen ^= part
    return answer


def literal_count(rows, mask):
    vertices = [u for u in range(len(rows)) if mask >> u & 1]
    return sum(all(rows[u] >> v & 1 for u, v in zip(word, word[1:]))
               for word in permutations(vertices))


def odd_cycles(rows):
    """Literal directed cycles modulo rotation only, retaining multiplicity."""
    n = len(rows)
    answer = []
    for length in range(3, n + 1, 2):
        for root in range(n):
            for tail in permutations(range(root + 1, n), length - 1):
                word = (root,) + tail
                if all(rows[word[i]] >> word[(i + 1) % length] & 1
                       for i in range(length)):
                    answer.append(sum(1 << u for u in word))
    return answer


def main():
    checks = Counter()
    pairs = {1: set(), 2: set()}
    witnesses = {}
    max_core = {2: 0, 4: 0}
    digest = sha256()
    inventory = []
    for n in range(1, 7):
        strong_count = 0
        for code in range(1 << (n * (n - 1) // 2)):
            rows = from_code(n, code)
            H = path_counts(rows)
            full = (1 << n) - 1
            comps = components(rows)
            strong_count += len(comps) == 1 and n >= 3
            F = H[full]
            core = sum(c for c in comps if c.bit_count() > 1)
            require(H[core] == F, 'canonical core H preservation')
            require(F % 2 == 1, 'Redei parity')
            product = 1
            for part in comps:
                product *= H[part]
            require(product == F, 'SCC product')
            checks['tournaments'] += 1
            cycles = odd_cycles(rows) if n <= 5 else None
            if n <= 5:
                require(F == literal_count(rows, full), 'literal path check')
                checks['literal_full_path_checks'] += 1
            for v in range(n):
                bit = 1 << v
                h = H[full ^ bit]
                delta = F - h
                require(delta >= 0 and delta % 2 == 0, 'even monotone deletion')
                part = next(c for c in comps if c & bit)
                require((delta == 0) == (part == bit), 'zero iff singleton SCC')
                checks['one_vertex_deletions'] += 1
                if cycles is not None:
                    weighted = sum(H[full ^ c] for c in cycles if c & bit)
                    require(delta == 2 * weighted, 'literal OCF deletion identity')
                    checks['literal_cycle_deletion_checks'] += 1
                if not delta:
                    continue
                K = delta // 2
                a = F // H[part]
                localK = (H[part] - H[part ^ bit]) // 2
                require(a * localK == K and K % a == 0, 'exact local scaling')
                require(1 <= a <= K, 'outside budget')
                require(part.bit_count() <= 2 * localK + 2, 'Moon affected-core bound')
                require(core.bit_count() <= 2 * K + 2, 'whole canonical core gain bound')
                require(2 * h <= factorial(2 * K + 1), 'sharpened fixed-increment old endpoint bound')
                require(F <= factorial(2 * K + 1) // 2 + 2 * K, 'sharpened final endpoint bound')
                checks['positive_increment_bounds'] += 1
                if K in pairs:
                    pair = (h, F)
                    pairs[K].add(pair)
                    key = f'{h}->{F}'
                    witnesses.setdefault(key, {'rows': rows, 'removed_vertex': v,
                                                'old_H': h, 'new_H': F})
                digest.update(json.dumps([n, code, v, h, F, a, localK],
                                         separators=(',', ':')).encode() + b'\n')
            # Checking all induced subsets is cheap once every subset count exists.
            for mask in range(1 << n):
                delta = F - H[mask]
                require(delta >= 0 and delta % 2 == 0, 'induced monotonicity')
                require(H[mask & core] == H[mask], 'simultaneous actual-pair kernel preserves old H')
                if delta:
                    require(core.bit_count() <= delta + 2, 'arbitrary induced gap core bound')
                    require(2 * H[mask] <= factorial(delta + 1), 'arbitrary induced gap endpoint bound')
                    checks['actual_induced_pair_kernel'] += 1
                if delta in (2, 4):
                    max_core[delta] = max(max_core[delta], core.bit_count())
                    expected = ({(1, 3), (3, 5)} if delta == 2 else
                                {(1, 5), (5, 9)})
                    require((H[mask], F) in expected, 'exact induced small-gap classification')
                    checks[f'induced_gap_{delta}'] += 1
        inventory.append({'n': n, 'labelled': 1 << (n * (n - 1) // 2),
                          'strong_labelled': strong_count})
    require(max_core == {2: 4, 4: 5}, 'exact small-gap maximum canonical core')
    require(pairs[1] == {(1, 3), (3, 5)}, 'complete +2 table')
    require(pairs[2] == {(1, 5), (5, 9)}, 'complete +4 table')
    # A same-count nonisomorphic representation cannot be freely inserted:
    # C3 ordinal C3 has H9, but no H11 tournament contains any H9 tournament.
    # This is a theorem consequence, not a finite search claim.
    require(path_counts([2 | 56, 4 | 56, 1 | 56, 16, 32, 8])[-1] == 9,
            'ordinal product positive control')
    for key, record in witnesses.items():
        rows = record['rows']
        F = literal_count(rows, (1 << len(rows)) - 1)
        v = record['removed_vertex']
        h = literal_count(rows, ((1 << len(rows)) - 1) ^ (1 << v))
        require((h, F) == (record['old_H'], record['new_H']), 'retained actual witnesses')
    # Explicit hostile mutation: every purported +2 pair above5 is rejected.
    def accepted_unit_pair(old, new):
        require((old, new) in {(1, 3), (3, 5)}, 'invalid unit-increment claim')
    try:
        accepted_unit_pair(9, 11)
    except RuntimeError:
        checks['hostile_mutation_rejected'] += 1
    else:
        raise RuntimeError('hostile mutation accepted')
    print(json.dumps({'status': 'PASS', 'universe': 'all labelled tournaments of orders1..6; every deletion and induced subset',
                      'checks': dict(sorted(checks.items())), 'inventory': inventory,
                      'one_vertex_pairs_by_half_gain': {str(k): sorted(v) for k, v in pairs.items()},
                      'exact_maximum_core_by_gain': max_core,
                      'witnesses': witnesses, 'semantic_sha256': digest.hexdigest(),
                      'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
                      'scope': 'Complete all-order +2/+4 pair classification uses the proved affected-SCC cutoff; higher bounded increments are not enumerated.'},
                     indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
