#!/usr/bin/env python3
"""Independent strong-tournament presentation cover and literal path counting.

Uses no producer imports, no subset DP, and no inherited census. Every strong
n-vertex tournament has a Hamilton cycle; relabel it 0->1->...->n-1->0. All
orientations of remaining pairs therefore cover every strong isomorphism class.
This is a presentation cover, not a count of isomorphism classes.
"""
from itertools import permutations
from hashlib import sha256
from pathlib import Path
import json


def check(ok, label):
    if not ok:
        raise RuntimeError(label)


def count(rows, vertices):
    result = 0
    for word in permutations(vertices):
        valid = True
        for u, v in zip(word, word[1:]):
            if not rows[u] & (1 << v):
                valid = False
                break
        result += valid
    return result


def main():
    pairs = {1: set(), 2: set()}
    counts = []
    total_deletions = 0
    digest = sha256()
    for n in range(3, 7):
        cycle_rows = [1 << ((u + 1) % n) for u in range(n)]
        free = [(u, v) for u in range(n) for v in range(u + 1, n)
                if not(cycle_rows[u] & (1 << v) or cycle_rows[v] & (1 << u))]
        per_n = 1 << len(free)
        for code in range(per_n):
            rows = list(cycle_rows)
            for j, (u, v) in enumerate(free):
                if code & (1 << j):
                    rows[u] |= 1 << v
                else:
                    rows[v] |= 1 << u
            F = count(rows, tuple(range(n)))
            values = []
            for v in range(n):
                h = count(rows, tuple(u for u in range(n) if u != v))
                check(F > h and (F - h) % 2 == 0, 'strict/even strong deletion')
                K = (F - h) // 2
                check(n <= 2 * K + 2, 'Moon affected-SCC cutoff')
                if K in pairs:
                    pairs[K].add((h, F))
                values.append(h)
                total_deletions += 1
            digest.update(json.dumps([n, code, F, values], separators=(',', ':')).encode() + b'\n')
        counts.append({'n': n, 'Hamilton_cycle_presentations': per_n})
    check(pairs[1] == {(1, 3), (3, 5)}, 'complete unit pairset')
    check(pairs[2] == {(1, 5), (5, 9)}, 'complete double pairset')
    hostile = [14, 4, 8, 18, 39, 79, 31]
    hostile_F = count(hostile, tuple(range(7)))
    hostile_h = count(hostile, tuple(range(1, 7)))
    check((hostile_h, hostile_F) == (17, 27), 'linear refinement hostile literal count')
    print(json.dumps({'linear_refinement_hostile': {'rows': hostile, 'deleted_vertex': 0, 'old_H': hostile_h, 'new_H': hostile_F, 'half_gain': 5},
                      'status': 'PASS', 'universe': 'all strong tournament isomorphism classes of orders3..6 via fixed Hamilton-cycle presentations',
                      'literal_deletion_checks': total_deletions, 'presentation_counts': counts,
                      'one_vertex_pairs_by_half_gain': {str(k): sorted(v) for k, v in pairs.items()},
                      'semantic_sha256': digest.hexdigest(),
                      'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
                      'scope': 'All-order +2/+4 spectra require the proof of the SCC cutoff and odd external factor, not finite extrapolation.'},
                     indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
