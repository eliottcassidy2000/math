#!/usr/bin/env python3
"""
CORRECTION CHECK: Is the hereditary property true for ALL maximizers
at odd n, or only REGULAR (vertex-transitive) maximizers?

Previous claim: "At odd n, every vertex deletion from maximizer gives
the (n-1)-maximizer."

Current data shows n=5 maximizer with del_hs=[3,5,5,5,3], where
H=3 ≠ max H(4)=5. This contradicts the claim!

Let's exhaustively check ALL maximizers at each n.

kind-pasteur-2026-03-06-S18g
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count

MAX_H = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189}

def score_seq(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

def delete_vertex(T, v):
    n = len(T)
    verts = [i for i in range(n) if i != v]
    return [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]

for n in range(3, 8):
    m = n * (n - 1) // 2
    max_h = MAX_H[n]
    max_h_prev = MAX_H[n - 1]

    print(f"\n{'='*60}")
    print(f"n={n}: H_max={max_h}, H_max({n-1})={max_h_prev}")
    print(f"{'='*60}")

    hereditary_count = 0
    non_hereditary_count = 0
    total_maximizers = 0
    by_score = {}

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h != max_h:
            continue

        total_maximizers += 1
        s = score_seq(T)
        del_hs = [hamiltonian_path_count(delete_vertex(T, v)) for v in range(n)]
        all_max = all(dh == max_h_prev for dh in del_hs)

        if all_max:
            hereditary_count += 1
        else:
            non_hereditary_count += 1

        if s not in by_score:
            by_score[s] = {'total': 0, 'hereditary': 0, 'non_hereditary': 0,
                           'example_del_hs': None}
        by_score[s]['total'] += 1
        if all_max:
            by_score[s]['hereditary'] += 1
        else:
            by_score[s]['non_hereditary'] += 1
            if by_score[s]['example_del_hs'] is None:
                by_score[s]['example_del_hs'] = del_hs

    print(f"Total maximizers: {total_maximizers}")
    print(f"  Hereditary: {hereditary_count}")
    print(f"  Non-hereditary: {non_hereditary_count}")

    for s in sorted(by_score.keys()):
        info = by_score[s]
        is_regular = len(set(s)) == 1
        print(f"\n  Score {s}{' (regular)' if is_regular else ''}:")
        print(f"    Total: {info['total']}, hereditary: {info['hereditary']}, "
              f"non-hereditary: {info['non_hereditary']}")
        if info['example_del_hs']:
            print(f"    Example non-hereditary del_hs: {sorted(info['example_del_hs'], reverse=True)}")

print("\nDone.")
