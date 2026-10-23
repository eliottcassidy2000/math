"""Independent static rooted-motif controls; no imports from the frozen producer."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
from math import comb
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def build(g, phase):
    n = 3 * g
    A = [[0] * n for _ in range(n)]
    for i in range(g):
        for x in range(3):
            A[3 * i + x][3 * i + (x + 1) % 3] = 1
    for i, j in combinations(range(g), 2):
        for x in range(3):
            for y in range(3):
                if j == i + 1:
                    bit = int(y == x)
                elif i == 0 and j == g - 1:
                    bit = int(y == (x + phase) % 3)
                else:
                    bit = 1
                A[3 * i + x][3 * j + y] = bit
                A[3 * j + y][3 * i + x] = 1 - bit
    return A


def encode(A, order):
    word = 0
    for i, j in combinations(range(len(order)), 2):
        word = 2 * word + A[order[i]][order[j]]
    return word


@lru_cache(None)
def rooted_canon(k, word):
    pairs = list(combinations(range(k), 2))
    A = [[0] * k for _ in range(k)]
    for bit_index, (i, j) in enumerate(reversed(pairs)):
        A[i][j] = (word >> bit_index) & 1
        A[j][i] = 1 - A[i][j]
    return min(encode(A, (0,) + tail) for tail in permutations(range(1, k)))


def profile(A, root, cutoff):
    n = len(A)
    others = [u for u in range(n) if u != root]
    output = []
    for k in range(1, cutoff + 1):
        counts = Counter(rooted_canon(k, encode(A, (root,) + S))
                         for S in combinations(others, k - 1))
        need(sum(counts.values()) == comb(n - 1, k - 1), 'rooted subset universe incomplete')
        output.append(tuple(sorted(counts.items())))
    return tuple(output)


def refine_once(A, colors):
    palette = sorted(set(colors))
    records = [(colors[u], tuple(sum(A[u][v] for v in range(len(A)) if colors[v] == color)
                                 for color in palette)) for u in range(len(A))]
    ids = {record: i for i, record in enumerate(sorted(set(records)))}
    return [ids[record] for record in records], records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    results = []
    for g in (3, 4, 5):
        graphs = [build(g, delta) for delta in (1, 2)]
        profiles = [[profile(A, root, g) for root in range(3 * g)] for A in graphs]
        for root in range(3 * g):
            need(profiles[0][root] == profiles[1][root], 'pointwise singly-rooted profile differs')
        for P in profiles:
            for i in range(g):
                need(P[3 * i] == P[3 * i + 1] == P[3 * i + 2], 'profile not uniform inside a cell')
        palette = {P: i for i, P in enumerate(sorted(set(profiles[0])))}
        colors = [[palette[P] for P in records] for records in profiles]
        need(colors[0] == colors[1], 'initial corresponding colors differ')
        need(len(set(colors[0])) == g, 'full rooted profiles do not identify intrinsic cells')
        updated = [refine_once(A, c) for A, c in zip(graphs, colors)]
        need(updated[0][1] == updated[1][1], 'ordinary neighbor records differ')
        need(updated[0][0] == colors[0] and updated[1][0] == colors[1], 'equitable profile partition is not stable')
        # Independently authenticate the point-fixing map on every full transversal.
        rooted_transversals = 0
        for root in range(3 * g):
            r, x = divmod(root, 3)
            for S in combinations([v for v in range(3 * g) if v != root], g - 1):
                vertices = (root,) + S
                if len({v // 3 for v in vertices}) != g:
                    continue
                mapped = tuple(3 * (v // 3) + (-v % 3 + 2 * x) % 3 for v in vertices)
                need(mapped[0] == root, 'transversal map moved the actual root')
                need(all(graphs[0][u][v] == graphs[1][uu][vv]
                         for u, uu in zip(vertices, mapped) for v, vv in zip(vertices, mapped)),
                     'root-fixed transversal inversion is not an induced isomorphism')
                rooted_transversals += 1
        row = {'cells': g, 'order': 3 * g, 'cutoff': g,
               'actual_roots_checked_in_each_graph': 3 * g,
               'all_rooted_induced_motif_profiles_equal_pointwise': True,
               'profile_cell_sizes': sorted(Counter(colors[0]).values()),
               'ordinary_refinement_is_already_stable': True,
               'root_fixed_transversal_checks': rooted_transversals,
               'profiles_at_cell_representatives': [profiles[0][3 * i] for i in range(g)],
               'pointwise_profile_sha256': [hashlib.sha256(json.dumps(P).encode()).hexdigest()
                                            for P in profiles[0]]}
        results.append(row)
    output = {'scope': 'Full singly-rooted uncolored induced motif counts at every actual vertex, through cutoff g, on both nonzero-phase 3g-vertex tournaments for g3,4,5; ordinary refinement is not individualization.',
              'controls': results, 'rooted_canonical_cache': list(rooted_canon.cache_info()),
              'checks_before_receipt_comparison': CHECKS}
    output = json.loads(json.dumps(output))
    receipt = HERE / 'rooted_cutoff_receipt.json'
    if args.write:
        receipt.write_text(json.dumps(output, indent=2) + '\n')
    else:
        need(json.loads(receipt.read_text()) == output, 'frozen rooted cutoff receipt differs')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
