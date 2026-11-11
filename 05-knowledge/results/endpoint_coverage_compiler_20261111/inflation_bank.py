"""Exact endpoint-cover coefficients for replacing an unselected vertex by TT_a.

Coefficients multiply P_j(a)=j! S(a,j), ordered nonempty path covers of TT_a.
The finite quotient DP counts visits to the contracted vertex, not large clones.
"""
from pathlib import Path
import sys
import json
import argparse
from itertools import combinations

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'modular_centers_collatz_pell_20261111'))
import content_checks as cc
import tournament_checks as tc


def visit_ends(rows, marked):
    n = len(rows)
    old = [v for v in range(n) if v != marked]
    order = old + [marked]
    adjacency = [[bool(rows[a] >> b & 1) for b in order] for a in order]
    limit = n
    dp = [[[0] * n for _ in range(limit + 1)] for _ in range(1 << (n - 1))]
    for mask in range(1 << (n - 1)):
        for visits in range(limit + 1):
            for v in range(n - 1):
                if not mask >> v & 1:
                    continue
                rest = mask ^ (1 << v)
                if not rest and visits == 0:
                    dp[mask][visits][v] = 1
                else:
                    dp[mask][visits][v] = sum(dp[rest][visits][u] for u in range(n) if adjacency[u][v])
            if visits:
                dp[mask][visits][-1] = (1 if mask == 0 and visits == 1 else 0)
                dp[mask][visits][-1] += sum(dp[mask][visits - 1][u] for u in range(n - 1) if adjacency[u][-1])
    return order, dp


def coefficients(rows, marked):
    n = len(rows)
    order, end = visit_ends(rows, marked)
    reverse = [sum(1 << b for b in range(n) if rows[b] >> a & 1) for a in range(n)]
    _, start = visit_ends(reverse, marked)
    full = (1 << (n - 1)) - 1
    h = [0] + [sum(end[full][j]) for j in range(1, n + 1)]
    u = [[0] * n for _ in range(n + 1)]
    w = [[[0] * n for _ in range(n)] for _ in range(n + 1)]
    for j in range(1, n + 1):
        directed = [[0] * n for _ in range(n)]
        for b in range(n - 1):
            u[j][order[b]] = start[full][j][b]  # empty first path
        for mask in range(full + 1):
            for first_j in range(j + 1):
                left, right = end[mask][first_j], start[full ^ mask][j - first_j]
                left_sum = sum(left)
                if not left_sum:
                    continue
                for b in range(n - 1):
                    if not right[b]:
                        continue
                    u[j][order[b]] += left_sum * right[b]
                    for a in range(n - 1):
                        directed[order[a]][order[b]] += left[a] * right[b]
        for a in range(n):
            for b in range(n):
                w[j][a][b] = directed[a][b] + directed[b][a]
    return h, u, w


def onto(a, j):
    # Independent ordinary set-partition recurrence for j! S(a,j).
    row = [1] + [0] * j
    for _ in range(a):
        row = [0] + [k * (row[k] + row[k - 1]) for k in range(1, j + 1)]
    return row[j]


def inflate(rows, marked, size):
    labels = [v for v in range(len(rows)) if v != marked] + [marked] * size
    out = [0] * len(labels)
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            forward = True if labels[i] == labels[j] else bool(rows[labels[i]] >> labels[j] & 1)
            out[i if forward else j] |= 1 << (j if forward else i)
    return labels, out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=6, choices=range(1, 8))
    args = parser.parse_args()
    bank, counts, records = [[]], [], []
    marked_checked, controls, presentations, negative = 0, 0, 0, 0
    defect_types, isolated_hits = {}, []
    for n in range(1, args.max_n + 1):
        codes = {cc.canonical(cc.child(rows, mask)) for rows in bank for mask in range(1 << (n - 1))}
        bank = [cc.decode(n, code) for code in sorted(codes)]
        counts.append(len(bank))
        for code, rows in zip(sorted(codes), bank):
            base_h, base_u, base_w, _, _, _ = tc.path_tables(rows)
            unit_sets = []
            for selection in range(1, 1 << n):
                old_selected = [v for v in range(n) if selection >> v & 1]
                if base_h - sum(base_w[a][b] for a, b in combinations(old_selected, 2)) == 1:
                    unit_sets.append((selection, old_selected))
            for marked in range(n):
                eligible = [(mask, selected) for mask, selected in unit_sets if not mask >> marked & 1]
                if not eligible:
                    continue
                h, u, w = coefficients(rows, marked)
                marked_checked += 1
                tc.require(h[1] == base_h, 'base H coefficient')
                for mask, selected in eligible:
                    delta = [h[j] - sum(w[j][a][b] for a, b in combinations(selected, 2)) for j in range(1, n + 1)]
                    slope = [sum(u[j][v] for v in selected) for j in range(1, n + 1)]
                    tc.require(delta[0] == 1 and slope[0] == sum(base_u[v] for v in selected), 'base seed coefficient')
                    stable = delta == [1] + [0] * (n - 1)
                    presentations += 1
                    negative += int(any(x < 0 for x in delta))
                    signature = tuple(delta[:max(i for i,x in enumerate(delta) if x)+1])
                    if signature not in defect_types:
                        defect_types[signature] = dict(n=n, code=code, rows=rows, selected=mask,
                                                      inflated_vertex=marked, coefficients=list(signature), count=0)
                    defect_types[signature]['count'] += 1
                    if not stable:
                        successes = [size for size in range(2,21)
                                     if sum(d*onto(size,j) for j,d in enumerate(delta,1)) == 1]
                        if successes:
                            isolated_hits.append(dict(n=n,code=code,selected=mask,inflated_vertex=marked,
                                                      defect_coefficients=delta,sizes_2_to_20=successes))
                    # Literal inflated-graph control, including hostile nonconstant defects.
                    if n <= 4 or stable:
                        for size in (2, 3):
                            labels, expanded = inflate(rows, marked, size)
                            selected_expanded = [i for i, v in enumerate(labels) if v in selected]
                            eh, eu, ew, _, _, _ = tc.path_tables(expanded)
                            actual_delta = eh - sum(ew[a][b] for a, b in combinations(selected_expanded, 2))
                            actual_slope = sum(eu[v] for v in selected_expanded)
                            weights = [onto(size, j) for j in range(1, n + 1)]
                            tc.require(actual_delta == sum(a*b for a,b in zip(delta, weights)), 'inflated defect')
                            tc.require(actual_slope == sum(a*b for a,b in zip(slope, weights)), 'inflated slope')
                            controls += 1
                    if stable:
                        last = max(i for i, x in enumerate(slope) if x)
                        coeff = slope[:last + 1]
                        record = dict(n=n, code=code, rows=rows, selected=mask, inflated_vertex=marked,
                                      defect_coefficients=delta, slope_cover_coefficients=coeff)
                        if len(coeff) <= 2:
                            a = coeff[1] if len(coeff) == 2 else 0
                            record['slope_A_2power_plus_B'] = [a, coeff[0] - 2*a]
                        records.append(record)
    useful = [r for r in records if r.get('slope_A_2power_plus_B', [0, 0])[0] and r['slope_A_2power_plus_B'][1]]
    result = {'ordinary_base_counts': counts, 'marked_bases_checked': marked_checked,
              'presentations_checked': presentations, 'negative_defect_presentations': negative,
              'literal_inflated_seed_controls': controls, 'stable_presentations': len(records),
              'nonzero_affine_constant_presentations': len(useful), 'records': records,
              'defect_types': list(defect_types.values()), 'nonstable_hits_sizes_2_to_20': isolated_hits}
    name = 'inflation_bank.json' if args.max_n == 6 else f'inflation_bank_{args.max_n}.json'
    (HERE / name).write_text(json.dumps(result, indent=2) + '\n', encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','defect_types','nonstable_hits_sizes_2_to_20')}, indent=2))
    print(json.dumps({'nonstable_hits': isolated_hits[:3]}, indent=2))


if __name__ == '__main__':
    main()
