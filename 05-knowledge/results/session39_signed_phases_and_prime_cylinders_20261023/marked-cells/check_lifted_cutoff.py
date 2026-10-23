"""Literal controls of the two-vertex lift and static rooted-profile transfer."""
from collections import Counter
from pathlib import Path
import argparse
import hashlib
import json
import check_rooted_cutoff as rooted

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def lift(A):
    n = len(A)
    B = [row + [1, 0] for row in A]
    B += [[0] * n + [0, 1], [1] * n + [0, 0]]
    need([sum(row) for row in B].count(1) == 1, 'new minimum degree not unique')
    need([sum(row) for row in B].count(n) == 1, 'new maximum degree not unique')
    need(sum(B[n]) == 1 and sum(B[n + 1]) == n, 'new extreme labels incorrect')
    recovered = [[B[u][v] for v in range(n)] for u in range(n)]
    need(recovered == A, 'deleting intrinsic extremes does not recover old tournament')
    return B


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    output = []
    for g in (3, 4, 5):
        graphs = [rooted.build(g, phase) for phase in (1, 2)]
        anti = [3 * (g - 1 - i) + (-x) % 3 for i in range(g) for x in range(3)]
        for r in range(3):
            n = len(graphs[0])
            u, v, w = 0, 1, 3 * (g - 1)
            need(not (graphs[0][u][v] == graphs[0][v][w] == graphs[0][w][u])
                 and graphs[1][u][v] == graphs[1][v][w] == graphs[1][w][u],
                 'indexed three-subset hostile failed')
            profiles = [[rooted.profile(A, u, g) for u in range(n)] for A in graphs]
            need(profiles[0] == profiles[1], 'lifted static rooted profiles differ')
            ids = {P: i for i, P in enumerate(sorted(set(profiles[0])))}
            colors = [[ids[P] for P in rec] for rec in profiles]
            expected_sizes = [1] * (2 * r) + [3] * g
            need(sorted(Counter(colors[0]).values()) == expected_sizes, 'lifted profile partition differs')
            records = [rooted.refine_once(A, cs) for A, cs in zip(graphs, colors)]
            need(records[0] == records[1], 'lifted neighbor records differ')
            need(records[0][0] == colors[0], 'lifted profile partition is not stable')
            for A in graphs:
                need(all(A[u][v] == A[anti[v]][anti[u]] for u in range(n) for v in range(n)),
                     'extended anti-isomorphism failed')
                for start in range(n):
                    reached = {start}
                    queue = [start]
                    for u in queue:
                        for v in range(n):
                            if A[u][v] and v not in reached:
                                reached.add(v)
                                queue.append(v)
                    need(len(reached) == n, 'lift not strong')
            output.append({'base_cells': g, 'lift_iterations': r, 'order': n, 'cutoff': g,
                           'all_actual_root_profiles_equal': True,
                           'indexed_three_subset_differs': [0, 1, 3 * (g - 1)],
                           'stable_profile_cell_sizes': expected_sizes,
                           'self_converse_and_strong': True,
                           'profile_sha256_by_vertex': [hashlib.sha256(json.dumps(P).encode()).hexdigest()
                                                       for P in profiles[0]]})
            if r < 2:
                graphs = [lift(A) for A in graphs]
                anti += [n + 1, n]
    arithmetic = []
    for n in range(11, 101):
        r = next(r for r in range(3) if (n - 2 * r) % 3 == 0)
        g = (n - 2 * r) // 3
        need(g >= 3 and n == 3 * g + 2 * r, 'all-order representation failed')
        need(3 * (g + 1) >= n - 1, 'uniform cutoff inequality failed')
        arithmetic.append([n, g, r, g + 1])
    result = {'scope': 'Actual two-vertex cyclic lifts preserve strongness, self-converse maps, intrinsic recovery, and every static singly-rooted profile through the base cutoff in the stated finite controls.',
              'dependency_sha256': hashlib.sha256((HERE / 'check_rooted_cutoff.py').read_bytes()).hexdigest(),
              'controls': output, 'arithmetic_controls_n_g_r_required_cutoff': arithmetic,
              'checks_before_receipt_comparison': CHECKS + rooted.CHECKS}
    result = json.loads(json.dumps(result))
    dest = HERE / 'lifted_cutoff_receipt.json'
    if args.write:
        dest.write_text(json.dumps(result, indent=2) + '\n')
    else:
        need(json.loads(dest.read_text()) == result, 'frozen lifted cutoff receipt differs')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
