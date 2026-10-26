#!/usr/bin/env python3
"""Finite controls for REPORT.md; standard library only; deterministic JSON.

This is not a halting decider or a proof of any famous conjecture.
All guards remain active under python -O.
"""
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def codec(n, edges):
    """edges are ordered (source, target, complete-label) records."""
    m = len(edges)
    N = n + 2 * m
    rows = [0] * N
    types = ['V'] * n + ['P'] * m + ['Q'] * m
    labels = [None] * n + [e[2] for e in edges] + [None] * m
    for u in range(N):
        for v in range(u + 1, N):
            if types[u] == types[v]:
                forward = True
            elif types[u] == 'V' and types[v] == 'P':
                forward = u == edges[v - n][0]
            elif types[u] == 'V' and types[v] == 'Q':
                forward = u != edges[v - n - m][1]
            else:
                require(types[u] == 'P' and types[v] == 'Q', 'bad types')
                forward = u - n == v - n - m
            a, b = (u, v) if forward else (v, u)
            rows[a] |= 1 << b
    return {'n': n, 'm': m, 'types': types, 'labels': labels, 'rows': rows}


def decode(c):
    """Validate and decode a canonical ordered presentation, independently."""
    n, m = c['n'], c['m']
    rows, types, labels = c['rows'], c['types'], c['labels']
    N = n + 2 * m
    require(n >= 1 and m >= 0, 'sizes')
    require(len(rows) == len(types) == len(labels) == N, 'lengths')
    require(types == ['V'] * n + ['P'] * m + ['Q'] * m, 'type partition')
    V, P, Q = range(n), range(n, n + m), range(n + m, N)
    def arc(a, b):
        return (rows[a] >> b) & 1
    for a in range(N):
        require(0 <= rows[a] < (1 << N) and not arc(a, a), 'loop/row')
        for b in range(a + 1, N):
            require(arc(a, b) + arc(b, a) == 1, 'not tournament')
            if types[a] == types[b]:
                require(arc(a, b), 'noncanonical filler')
    paired = set()
    out = []
    for p in P:
        qmatches = [q for q in Q if arc(p, q)]
        require(len(qmatches) == 1, 'P/Q matching')
        q = qmatches[0]
        require(q not in paired, 'noninjective matching')
        paired.add(q)
        # Q order is part of this canonical codec, not a graph invariant.
        require(q - n - m == p - n, 'noncanonical pairing')
        sources = [v for v in V if arc(v, p)]
        targets = [v for v in V if arc(q, v)]
        require(len(sources) == len(targets) == 1, 'source/target incidence')
        out.append((sources[0], targets[0], labels[p]))
    require(len(paired) == m, 'incomplete matching')
    return out


def all_codec_controls():
    counts = []
    total_vertices = 0
    max_vertices = 0
    digest = hashlib.sha256()
    for n in range(1, 4):
        count = 0
        for targets in itertools.product(range(-1, n), repeat=2 * n):
            edges = []
            for port, t in enumerate(targets):
                if t >= 0:
                    source, read = divmod(port, 2)
                    # Full labels deliberately differ even with same endpoints.
                    label = f'read={read};write={(source + read) % 2};move={(-1, 1)[read]}'
                    edges.append((source, t, label))
            c = codec(n, edges)
            require(decode(c) == edges, 'roundtrip')
            N = len(c['rows'])
            require(N == n + 2 * len(edges) and N <= 5 * n, 'overhead')
            digest.update(json.dumps(c, sort_keys=True).encode())
            total_vertices += N
            max_vertices = max(max_vertices, N)
            count += 1
        require(count == (n + 1) ** (2 * n), 'enumeration count')
        counts.append({'states': n, 'tables': count})
    require(sum(x['tables'] for x in counts) == 4181, 'total')
    sample_edges = [(0, 0, 'loop:a'), (0, 0, 'loop:b'), (0, 1, 'a'), (1, 0, 'b')]
    c = codec(2, sample_edges)
    require(decode(c) == sample_edges, 'loops/parallel')
    rejected = 0
    for a in range(len(c['rows'])):
        for b in range(a + 1, len(c['rows'])):
            bad = {**c, 'rows': c['rows'][:]}
            bad['rows'][a] ^= 1 << b
            bad['rows'][b] ^= 1 << a
            try:
                decode(bad)
            except RuntimeError:
                rejected += 1
            else:
                raise RuntimeError(f'single flipped arc survived: {a},{b}')
    require(rejected == 45, 'mutation count')
    return {'by_state_count': counts, 'total': 4181, 'total_encoded_vertices': total_vertices,
            'maximum_encoded_vertices': max_vertices, 'canonical_single_arc_mutations_rejected': rejected,
            'corpus_sha256': digest.hexdigest()}


def tournament_controls():
    all_rows = []
    for n in range(1, 6):
        pairs = list(itertools.combinations(range(n), 2))
        permutations = list(itertools.permutations(range(n)))
        count = 0
        aut_sum = 0
        aut_sizes = set()
        for mask in range(1 << len(pairs)):
            edges = {(a, b) if mask & (1 << i) else (b, a) for i, (a, b) in enumerate(pairs)}
            autos = [p for p in permutations if all((p[a], p[b]) in edges for a, b in edges)]
            aut_sizes.add(len(autos))
            aut_sum += len(autos)
            require(len(autos) % 2 == 1, 'even automorphism group')
            for p in autos:
                involution = all(p[p[a]] == a for a in range(n))
                identity = all(p[a] == a for a in range(n))
                require(not involution or identity, 'nontrivial tournament involution')
            count += 1
        all_rows.append({'n': n, 'labelled_tournaments': count, 'automorphism_sizes': sorted(aut_sizes),
                         'automorphism_order_sum': aut_sum,
                         'unlabelled_orbits_by_Burnside': aut_sum // len(permutations)})
    require(sum(x['labelled_tournaments'] for x in all_rows) == 1099, 'tournament total')
    require([x['unlabelled_orbits_by_Burnside'] for x in all_rows] == [1, 1, 2, 4, 12], 'orbit control')
    return all_rows


def simulate(table, max_steps=100):
    q, h, tape = 'A', 0, {}
    seen = {}
    trace = []
    for step in range(max_steps):
        configuration = (q, h, tuple(sorted(tape.items())))
        if configuration in seen:
            return {'outcome': 'repeated_configuration', 'step': step, 'first_seen': seen[configuration], 'trace': trace}
        seen[configuration] = step
        read = tape.get(h, 0)
        inst = table[(q, read)]
        trace.append({'step': step, 'q': q, 'head': h, 'ones': sorted(tape), 'read': read, 'instruction': inst})
        if inst is None:
            return {'outcome': 'halt', 'executed_instructions': step + 1, 'trace': trace}
        write, move, q = inst
        if write:
            tape[h] = write
        else:
            tape.pop(h, None)
        h += move
    return {'outcome': 'bounded_prefix_only', 'steps': max_steps, 'trace': trace}


def machine_hostile():
    a = {('A', 0): (0, 1, 'B'), ('A', 1): None,
         ('B', 0): (0, -1, 'A'), ('B', 1): (1, -1, 'A')}
    b = {**a, ('A', 0): (1, 1, 'B')}
    def skeleton(t):
        return sorted((q, read, 'HALT' if inst is None else inst[2]) for (q, read), inst in t.items())
    require(skeleton(a) == skeleton(b), 'skeleton not equal')
    r0, r1 = simulate(a), simulate(b)
    require(r0['outcome'] == 'repeated_configuration' and r0['step'] == 2, 'zero writer')
    require(r1['outcome'] == 'halt' and r1['executed_instructions'] == 3, 'one writer')
    return {'same_source_read_target_skeleton': skeleton(a), 'write_zero_machine': r0, 'write_one_machine': r1}


def lrc_exact(v, eps=None):
    """Finite closed rational intervals for one tuple, no conjecture used."""
    eps = Fraction(1, len(v) + 1) if eps is None else eps
    surviving = [(Fraction(0), Fraction(1))]
    for x in v:
        require(x > 0, 'positive speed')
        good = [(Fraction(j, x) + eps / x, Fraction(j + 1, x) - eps / x) for j in range(x)]
        new = []
        for a, b in surviving:
            for c, d in good:
                lo, hi = max(a, c), min(b, d)
                if lo <= hi:
                    new.append((lo, hi))
        surviving = new
    return surviving


def arithmetic_controls():
    root_pairs = []
    for k in range(1, 5):
        c = Fraction(1, 2 ** k - 3)
        h = Fraction(3 ** (k - 1), 3 ** k - 2 ** k)
        require((3 * c + 1) / 2 ** k == c, 'Collatz fixed equation')
        require(3 ** (k - 1) * (3 * h - 1) / 2 ** k == h, 'twisted fixed equation')
        root_pairs.append({'k': k, 'collatz': str(c), 'twisted': str(h)})
    require(root_pairs[0]['collatz'] == '-1' and root_pairs[1]['twisted'] == '3/5', 'domain hostile')
    tuples = []
    for k in range(1, 6):
        speeds = tuple(range(1, k + 1))
        intervals = lrc_exact(speeds)
        require(intervals, 'AP positive control')
        t = intervals[0][0]
        min_distance = min(min((x * t) % 1, 1 - (x * t) % 1) for x in speeds)
        require(min_distance >= Fraction(1, k + 1), 'exact witness')
        tuples.append({'speeds': list(speeds), 'first_witness': str(t), 'min_distance': str(min_distance),
                       'closed_components': [[str(a), str(b)] for a, b in intervals]})
    # The ordinary 1/2 threshold for speeds 1 and 2 is impossible; this does
    # not use the conjectural 1/3 threshold of lrc_exact.
    require(not lrc_exact((1, 2), Fraction(1, 2)), 'strict threshold hostile')
    return {'constant_valuation_fixed_points': root_pairs, 'finite_LRC_positive_controls': tuples,
            'impossible_half_threshold_for_speeds_1_2': lrc_exact((1, 2), Fraction(1, 2))}


def source_integrity():
    root = Path(__file__).resolve().parent
    data = json.loads((root / 'SOURCE_PINS.json').read_text())
    for row in data['files']:
        got = hashlib.sha256((root / row['snapshot']).read_bytes()).hexdigest()
        require(got == row['sha256'], 'source hash mismatch: ' + row['snapshot'])
    return {'commit': data['commit'], 'verified_source_snapshots': len(data['files'])}


def main():
    print(json.dumps({'status': 'PASS_FINITE_CONTROLS_ONLY', 'codec': all_codec_controls(),
                      'tournaments': tournament_controls(), 'write_label_hostile': machine_hostile(),
                      'arithmetic': arithmetic_controls(), 'source_integrity': source_integrity()},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
