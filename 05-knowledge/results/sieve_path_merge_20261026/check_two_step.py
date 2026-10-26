#!/usr/bin/env python3
"""Independent replay of incoming two-step decoder on every labelled T, n<=6.

Only subset-DP used_arcs and tournament construction are imported from the
previous independent audit; its main is not run. No incoming producer imports.
"""
import argparse
from collections import Counter
from itertools import combinations
from pathlib import Path
import hashlib
import json
import runpy

GATES = 0


def need(value, message):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(message)


def reachable(rows, start):
    seen, pending = {start}, [start]
    while pending:
        u = pending.pop()
        for v in range(len(rows)):
            if rows[u] >> v & 1 and v not in seen:
                seen.add(v)
                pending.append(v)
    return seen


def restore(used):
    n = len(used)
    reach = [reachable(used, u) for u in range(n)]
    remaining, blocks = set(range(n)), []
    while remaining:
        u = min(remaining)
        block = {v for v in remaining if v in reach[u] and u in reach[v]}
        blocks.append(block)
        remaining -= block
    owner = {v: i for i, block in enumerate(blocks) for v in block}
    answer, records = list(used), []
    for x, y in combinations(range(n), 2):
        if used[x] >> y & 1 or used[y] >> x & 1:
            continue
        if owner[x] != owner[y]:
            xy, yx = y in reach[x], x in reach[y]
            need(xy != yx, 'promised support components must be linearly ordered')
            u, v = (x, y) if xy else (y, x)
            records.append({'arc': [u, v], 'kind': 'outer'})
        else:
            block = blocks[owner[x]]
            xy = [w for w in block if used[x] >> w & 1 and used[w] >> y & 1]
            yx = [w for w in block if used[y] >> w & 1 and used[w] >> x & 1]
            need((len(xy) >= 2) != (len(yx) >= 2), 'exactly one missing-pair two-step count must reach two')
            u, v, good, bad = (x, y, xy, yx) if len(xy) >= 2 else (y, x, yx, xy)
            need(len(good) >= 2 and len(bad) <= 1, 'direct decoder bounds')
            records.append({'arc': [u, v], 'kind': 'internal', 'forward_middle_vertices': good, 'reverse_middle_vertices': bad})
        answer[u] |= 1 << v
    return tuple(answer), records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dp-source', type=Path,
        default=Path(__file__).resolve().parents[1] / 'growing_sieve_path_20261026/audit_paths/check.py')
    args = parser.parse_args()
    source = args.dp_source.read_bytes()
    expected = hashlib.sha256(source).hexdigest()
    module = runpy.run_path(str(args.dp_source), run_name='independent_subset_dp_support')
    used_arcs, tournament = module['used_arcs'], module['tournament']
    need(hashlib.sha256(args.dp_source.read_bytes()).hexdigest() == expected, 'DP source unchanged during load')
    census, first_internal = [], None
    histogram = Counter()
    total_tournaments = 0
    for n in range(1, 7):
        total = 1 << (n*(n-1)//2)
        inner = outer = nonfull = 0
        for code in range(total):
            actual = tournament(n, code)
            support = used_arcs(actual)
            restored, certificates = restore(support)
            need(restored == actual, ('two-step inverse disagrees with literal T', n, code))
            total_tournaments += 1
            nonfull += support != actual
            for cert in certificates:
                x, y = cert['arc']
                need(actual[x] >> y & 1, ('certificate direction oracle', n, code, x, y))
                if cert['kind'] == 'internal':
                    inner += 1
                    histogram[len(cert['forward_middle_vertices']), len(cert['reverse_middle_vertices'])] += 1
                    if first_internal is None:
                        first_internal = {'n': n, 'code': code, 'rows': actual, 'support': support, 'certificate': cert}
                else:
                    outer += 1
        census.append({'n': n, 'labelled_tournaments': total, 'nonfull_supports': nonfull,
                       'missing_internal_arcs': inner, 'missing_outer_arcs': outer})
    need(total_tournaments == 33867, 'entire labelled universe size')
    need(first_internal['n'] == 5 and first_internal['code'] == 10, 'first strong internal boundary')
    need(first_internal['certificate'] == {'arc': [4, 1], 'kind': 'internal', 'forward_middle_vertices': [2, 3], 'reverse_middle_vertices': [0]}, 'sharp threshold-two witness')

    transitive = (6, 4, 0)
    support = used_arcs(transitive)
    forward = sum(bool(support[0] >> w & 1 and support[w] >> 2 & 1) for w in range(3))
    need(forward == 1 and restore(support)[0] == transitive, 'outer missing arc requires component-order rule')
    paley = tuple(sum(1 << v for v in range(7) if (v-u) % 7 in (1, 2, 4)) for u in range(7))
    p_support = used_arcs(paley)
    c01 = sum(bool(p_support[0] >> w & 1 and p_support[w] >> 1 & 1) for w in range(7))
    c10 = sum(bool(p_support[1] >> w & 1 and p_support[w] >> 0 & 1) for w in range(7))
    need(p_support == paley and paley[0] >> 1 & 1, 'Paley observed arc control')
    need((c01, c10) == (1, 2), 'two-step rule must not be applied to already observed arcs')
    need(restore(p_support)[0] == paley, 'existing arcs retained')
    print(json.dumps({'status': 'FINITE-EXACT PASS', 'gates': GATES,
        'universe': 'Every labelled tournament at orders1..6, without quotient/filter; one explicit Paley7 scope hostile',
        'labelled_tournaments': total_tournaments, 'dp_source_sha256': expected,
        'census': census, 'missing_internal_count_histogram': [dict(forward=a, reverse=b, arcs=v) for (a,b),v in sorted(histogram.items())],
        'first_internal': first_internal,
        'scope_hostiles': {'transitive3_missing_outer_forward_count': forward,
                          'Paley7_existing_arc_0_to_1_counts': [c01, c10]},
        'method': 'Runpy imports only prior independent subset-endpoint DP; fresh traversal SCCs and direct two-step decoder; literal tournament rows are oracle.'}, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
