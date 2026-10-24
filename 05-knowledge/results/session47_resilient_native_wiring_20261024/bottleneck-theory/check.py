#!/usr/bin/env python3
"""Independent finite controls for exact weighted bottlenecks and compression.

No upstream checker is imported. Every check remains active under python -O.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
INF = math.inf
CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def components(n, edges):
    adjacency = [[] for _ in range(n)]
    for i, j in edges:
        adjacency[i].append(j)
        adjacency[j].append(i)
    result = [-1] * n
    for start in range(n):
        if result[start] != -1:
            continue
        result[start] = start
        stack = [start]
        while stack:
            v = stack.pop()
            for w in adjacency[v]:
                if result[w] == -1:
                    result[w] = start
                    stack.append(w)
    return result


def minimax(matrix):
    distance = [row[:] for row in matrix]
    n = len(matrix)
    for k in range(n):
        for i in range(n):
            for j in range(n):
                distance[i][j] = min(distance[i][j],
                                     max(distance[i][k], distance[k][j]))
    return distance


def cut_value(matrix, mask):
    return min((matrix[i][j] for i in range(len(matrix))
                if mask >> i & 1 for j in range(len(matrix))
                if not mask >> j & 1), default=INF)


def kruskal(matrix):
    n = len(matrix)
    labels = list(range(n))
    selected = []
    for w, i, j in sorted((matrix[i][j], i, j)
                          for i in range(n) for j in range(i + 1, n)
                          if math.isfinite(matrix[i][j])):
        a, b = labels[i], labels[j]
        if a != b:
            selected.append((i, j))
            labels = [a if t == b else t for t in labels]
    return selected


def tree_matrix(matrix, edges):
    n = len(matrix)
    result = [[0 if i == j else INF for j in range(n)] for i in range(n)]
    for i, j in edges:
        result[i][j] = result[j][i] = matrix[i][j]
    return result


def graph_controls():
    rows = []
    for n in range(1, 5):
        pairs = list(itertools.combinations(range(n), 2))
        trees = [edges for edges in itertools.combinations(pairs, n - 1)
                 if len(set(components(n, edges))) == 1]
        count = connected = mst_count = 0
        for values in itertools.product((0, 1, 3, INF), repeat=len(pairs)):
            matrix = [[0 if i == j else INF for j in range(n)] for i in range(n)]
            for (i, j), value in zip(pairs, values):
                matrix[i][j] = matrix[j][i] = value
            beta = minimax(matrix)
            cuts = {mask: cut_value(matrix, mask)
                    for mask in range(1, (1 << n) - 1)}
            for i, j in pairs:
                expected = max(v for mask, v in cuts.items()
                               if mask >> i & 1 and not mask >> j & 1)
                need(beta[i][j] == expected, 'Pair max-cut formula failed')
            for i, j, k in itertools.product(range(n), repeat=3):
                need(beta[i][k] <= max(beta[i][j], beta[j][k]),
                     'Ultrametric inequality failed')
            for budget in (0, 1, 3):
                cc = components(n, [(i, j) for i, j in pairs
                                    if matrix[i][j] <= budget])
                need(all((cc[i] == cc[j]) == (beta[i][j] <= budget)
                         for i in range(n) for j in range(n)),
                     'Independent threshold connectivity failed')
            global_b = max((beta[i][j] for i, j in pairs), default=0)
            need(global_b == max(cuts.values(), default=0), 'Global cut maximum failed')
            finite_trees = [t for t in trees
                            if all(math.isfinite(matrix[i][j]) for i, j in t)]
            best_tree_peak = min((max((matrix[i][j] for i, j in t), default=0)
                                  for t in finite_trees), default=INF)
            need(global_b == best_tree_peak, 'Exhaustive spanning-tree peak failed')
            forest = kruskal(matrix)
            need(minimax(tree_matrix(matrix, forest)) == beta,
                 'Kruskal forest did not preserve all minimax distances')
            if finite_trees:
                connected += 1
                total = min(sum(matrix[i][j] for i, j in t) for t in finite_trees)
                for t in finite_trees:
                    if sum(matrix[i][j] for i, j in t) == total:
                        mst_count += 1
                        need(minimax(tree_matrix(matrix, t)) == beta,
                             'An exhaustive minimum-total tree failed')
            # Vertex-supply deletion cannot improve fixed-endpoint routes.
            for removed in range(n):
                keep = [i for i in range(n) if i != removed]
                restricted = minimax([[matrix[i][j] for j in keep] for i in keep])
                need(all(beta[i][j] <= restricted[a][b]
                         for a, i in enumerate(keep) for b, j in enumerate(keep)),
                     'Fixed-endpoint supply monotonicity failed')
            count += 1
        rows.append({'vertices': n, 'weighted_graphs': count,
                     'connected_graphs': connected, 'minimum_total_trees_checked': mst_count})
    return {'edge_alphabet': [0, 1, 3, 'absent'], 'rows': rows,
            'total_graphs': sum(row['weighted_graphs'] for row in rows)}


class Poset:
    def __init__(self, n, edges):
        self.n = n
        self.all = (1 << n) - 1
        relation = [[i == j for j in range(n)] for i in range(n)]
        for i, j in edges:
            relation[i][j] = True
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    relation[i][j] |= relation[i][k] and relation[k][j]
        self.relation = relation
        self.upsets = [m for m in range(1 << n) if self.up(m) == m]
        self.intervals = [(a, b) for a in self.upsets for b in self.upsets if a & ~b == 0]

    def up(self, mask):
        return sum(1 << j for j in range(self.n)
                   if any(mask >> i & 1 and self.relation[i][j] for i in range(self.n)))

    def down(self, mask):
        return sum(1 << i for i in range(self.n)
                   if any(mask >> j & 1 and self.relation[i][j] for j in range(self.n)))

    def members(self, interval):
        low, high = interval
        return [m for m in self.upsets if low & ~m == 0 and m & ~high == 0]

    def route(self, source, target):
        current = source
        route = [current]
        while current & ~target:
            remaining = [i for i in range(self.n) if (current & ~target) >> i & 1]
            x = next(i for i in remaining if not any(
                self.relation[j][i] for j in remaining if j != i))
            current ^= 1 << x
            route.append(current)
        while target & ~current:
            remaining = [i for i in range(self.n) if (target & ~current) >> i & 1]
            x = next(i for i in remaining if not any(
                self.relation[i][j] for j in remaining if j != i))
            current ^= 1 << x
            route.append(current)
        return route


def posets(n):
    pairs = list(itertools.combinations(range(n), 2))
    unique = {}
    for flags in itertools.product((False, True), repeat=len(pairs)):
        p = Poset(n, [pair for pair, flag in zip(pairs, flags) if flag])
        key = tuple(tuple(row) for row in p.relation)
        unique[key] = p
    return list(unique.values())


def weight_table(weights):
    return [sum(w for i, w in enumerate(weights) if m >> i & 1)
            for m in range(1 << len(weights))]


def boundary(a, b, wt):
    la, ra = a
    lb, rb = b
    return wt[la & ~rb] + wt[lb & ~ra], la | (lb & ra), lb | (la & rb)


def interval_controls():
    rows = []
    for n in range(1, 5):
        closure_count = pair_count = route_count = 0
        pp = posets(n)
        for p in pp:
            for forced_one, forced_zero in itertools.product(range(1 << n), repeat=2):
                low, high = p.up(forced_one), p.all ^ p.down(forced_zero)
                literal = [m for m in p.upsets if forced_one & ~m == 0 and not forced_zero & m]
                compiled = [] if low & ~high else p.members((low, high))
                need(literal == compiled, 'Forced interval compilation failed')
                closure_count += 1
            members = {interval: p.members(interval) for interval in p.intervals}
            for interval, states in members.items():
                for source, target in itertools.product(states, repeat=2):
                    route = p.route(source, target)
                    need(route[0] == source and route[-1] == target, 'Route endpoints failed')
                    need(all(m in states for m in route), 'Internal route left interval')
                    changes = [a ^ b for a, b in zip(route, route[1:])]
                    need(all(m.bit_count() == 1 for m in changes), 'Route changed multiple components')
                    need(len(changes) == (source ^ target).bit_count(), 'Route was not conformal')
                    route_count += 1
            for weights in ((1,) * n, tuple(i % 3 for i in range(n)),
                            tuple(i + 1 for i in range(n))):
                wt = weight_table(weights)
                for a, b in itertools.product(p.intervals, repeat=2):
                    gap, x, y = boundary(a, b, wt)
                    literal = min(wt[u ^ v] for u in members[a] for v in members[b])
                    need(gap == literal, 'Boundary formula not the exact minimum')
                    need(x in members[a] and y in members[b], 'Boundary witness infeasible')
                    need(wt[x ^ y] == gap, 'Boundary witness cost mismatch')
                    pair_count += 1
        rows.append({'components': n, 'naturally_labelled_posets': len(pp),
                     'forced_constraint_pairs': closure_count,
                     'weighted_interval_pairs': pair_count,
                     'internal_routes': route_count})
    return rows


def charge(mask, labels):
    return sum(q for i, q in enumerate(labels) if mask >> i & 1) % 2


def quotient_instance(p, helper_intervals, helper_charges, labels, weights):
    wt = weight_table(weights)
    states = [(h, m) for h, interval in enumerate(helper_intervals)
              for m in p.members(interval)
              if (helper_charges[h] + charge(m, labels)) % 2 == 0]
    cost = [[(h != g) + wt[u ^ v] for g, v in states] for h, u in states]
    full_beta = minimax(cost)
    previous_projection = None
    budget_checks = 0
    for budget in range(sum(weights) + 2):
        visible = sum(1 << i for i in range(p.n) if labels[i] or weights[i] > budget)
        neutral = p.all ^ visible
        quotient = []
        for h, (lh, rh) in enumerate(helper_intervals):
            for t in range(1 << p.n):
                if t & ~visible or (helper_charges[h] + charge(t, labels)) % 2:
                    continue
                low = p.up(lh | t)
                high = p.all ^ p.down((p.all ^ rh) | (visible ^ t))
                if low & ~high:
                    continue
                need(low & visible == t and high & visible == t, 'Visible closure mismatch')
                quotient.append((h, t, low & neutral, high & neutral))
        lookup = {(h, t): j for j, (h, t, lo, hi) in enumerate(quotient)}
        projection = [lookup[(h, m & visible)] for h, m in states]
        need(set(projection) == set(range(len(quotient))), 'Quotient vertices not exactly feasible')
        groups = [[i for i, q in enumerate(projection) if q == j] for j in range(len(quotient))]
        for j, (h, t, lo, hi) in enumerate(quotient):
            literal = {states[i][1] & neutral for i in groups[j]}
            compiled = {m for m in range(1 << p.n) if not m & ~neutral
                        and not lo & ~m and not m & ~hi
                        and all(not (m >> x & 1 and neutral >> y & 1)
                                or (m >> y & 1) for x in range(p.n) for y in range(p.n)
                                if p.relation[x][y])}
            need(literal == compiled, 'Neutral fibre membership mismatch')
            cc = components(len(groups[j]), [(a, b) for a in range(len(groups[j]))
                            for b in range(a + 1, len(groups[j]))
                            if cost[groups[j][a]][groups[j][b]] <= budget])
            need(len(set(cc)) == 1, 'Projected fibre not internally budget-connected')
        edges = []
        for j in range(len(quotient)):
            h, t, la, ra = quotient[j]
            for k in range(j + 1, len(quotient)):
                g, u, lb, rb = quotient[k]
                gap, wa, wb = boundary((la, ra), (lb, rb), wt)
                crossing = (h != g) + wt[t ^ u] + gap
                literal = min(cost[a][b] for a in groups[j] for b in groups[k])
                need(crossing == literal, 'Quotient crossing minimum mismatch')
                need((h, t | wa) in states and (g, u | wb) in states,
                     'Explicit full crossing witness not feasible')
                if crossing <= budget:
                    edges.append((j, k))
        cc = components(len(quotient), edges)
        need(all((cc[projection[i]] == cc[projection[j]]) == (full_beta[i][j] <= budget)
                 for i in range(len(states)) for j in range(len(states))),
             'Exact quotient lifting failed')
        if previous_projection is not None:
            old_keys, old_projection = previous_projection
            for j, (h, t) in enumerate(old_keys):
                reduced = lookup[(h, t & visible)]
                need(all(projection[i] == reduced for i, v in enumerate(old_projection) if v == j),
                     'Budget projections did not commute')
        previous_projection = ([(h, t) for h, t, lo, hi in quotient], projection)
        budget_checks += 1
    return budget_checks


def quotient_controls():
    rows = []
    for n in (1, 2):
        instances = budgets = empty = 0
        for p in posets(n):
            for intervals in itertools.product(p.intervals, repeat=2):
                for hcharges in itertools.product((0, 1), repeat=2):
                    for labels in itertools.product((0, 1), repeat=n):
                        for weights in itertools.product((0, 1, 2), repeat=n):
                            budgets += quotient_instance(p, intervals, hcharges, labels, weights)
                            instances += 1
        rows.append({'components': n, 'scope': 'exhaustive two-helper intervals, C2 charges and weights 0/1/2',
                     'instances': instances, 'budget_checks': budgets})
    # All seven naturally labelled three-element posets and every ordered pair
    # of helper intervals; two contrasting charge patterns, two helper-charge pairs.
    instances = budgets = 0
    for p in posets(3):
        for intervals in itertools.product(p.intervals, repeat=2):
            for hcharges in ((0, 0), (0, 1)):
                for labels in ((0, 0, 0), (0, 1, 0)):
                    budgets += quotient_instance(p, intervals, hcharges, labels, (1, 2, 3))
                    instances += 1
    rows.append({'components': 3, 'scope': 'all helper-interval pairs; weights 1/2/3; stated C2 charge patterns',
                 'instances': instances, 'budget_checks': budgets})
    return rows


def counterexamples():
    unit = weight_table((1,))
    off, free, on = (0, 0), (0, 1), (1, 1)
    a = boundary(off, free, unit)[0]
    b = boundary(free, on, unit)[0]
    c = boundary(off, on, unit)[0]
    need((a, b, c) == (0, 0, 1) and c > a + b, 'Nonmetric boundary control failed')
    need(set((0,)).isdisjoint((1,)) and 1 + c == 2, 'Intersection-only control failed')
    zero_gap = boundary(off, on, weight_table((0,)))[0]
    need(zero_gap == 0, 'Zero cost disjoint intervals control failed')
    heavy = [[0, 2], [2, 0]]
    need(len(set(components(2, []))) == 2 and minimax(heavy)[0][1] == 2,
         'Heavy-neutral lower bound failed')
    need(components(1, []) == [0], 'Erroneous collapse control failed')
    old = [[0, 1], [1, 0]]
    enlarged = [[0, 1, 9], [1, 0, 9], [9, 9, 0]]
    need(max(map(max, minimax(old))) == 1 and max(map(max, minimax(enlarged))) == 9,
         'Universal-supply nonmonotonicity control failed')
    need(minimax(enlarged)[0][1] == 1, 'Old endpoint was not preserved')
    uniform = sum((2,) * 11)
    mixed = sum((2,) * 10 + (3,))
    need((uniform, mixed) == (22, 23), 'Clock scaling counterexample failed')
    need(sum((2, 2, 2, 2, 3, 3, 4, 4)) == 22, 'Eight-slot cost sum failed')
    def recode(v):
        return v ^ (0b1010 if v & 1 else 0)
    need(recode(11) == 1, 'Parity 11-to-1 recoding failed')
    need(len({recode(v) for v in range(16)}) == 16, 'Parity recoding not invertible')
    need(all(recode(a ^ b) == recode(a) ^ recode(b)
             for a in range(16) for b in range(16)), 'Parity recoding not additive')
    diagonal_rows = []
    for m in range(1, 13):
        # Abstract records induced by inherited diagonal C3 description only.
        states = [tuple(tuple((phase + j) % 3 for j in range(3))
                        for _ in range(m + 1)) for phase in range(3)]
        matrix = [[sum(a != b for ra, rb in zip(s, t) for a, b in zip(ra, rb))
                   for t in states] for s in states]
        need(all(matrix[i][j] == (0 if i == j else 3 * (m + 1))
                 for i in range(3) for j in range(3)), 'Diagonal C3 cost failed')
        need(cut_value(matrix, 1) == 3 * (m + 1), 'Diagonal C3 cut failed')
        need(max(map(max, minimax(matrix))) == 3 * (m + 1), 'Diagonal C3 bottleneck failed')
        diagonal_rows.append({'old_slots': m, 'total_slots': m + 1,
                              'clock_bottleneck': 3 * (m + 1)})
    return {'interval_gap_triangle': [a, b, c], 'heavy_neutral_weight': 2,
            'heavy_neutral_budget': 1, 'disjoint_interval_valid_bridge': 2,
            'zero_weight_disjoint_interval_gap': zero_gap,
            'enlarged_universal_budget': [1, 9],
            'same_11_slot_threshold_clock_budgets': [uniform, mixed],
            'eight_slot_22_weights': [2, 2, 2, 2, 3, 3, 4, 4],
            'parity_recoding': {'input': 11, 'output': recode(11), 'linear_bijection': True},
            'diagonal_C3_record_models_not_native_replays': diagonal_rows}


def supply_controls():
    deletion_tests = helper_pair_tests = 0
    for r in range(1, 4):
        for t in range(4):
            packing = [set(range(j * r, (j + 1) * r)) for j in range(t + 1)]
            universe = range((t + 1) * r)
            for size in range(t + 1):
                for deleted in itertools.combinations(universe, size):
                    need(any(not block.intersection(deleted) for block in packing),
                         'Disjoint packing failed bounded forbidden-U coverage')
                    deletion_tests += 1
            if t:
                insufficient = packing[:t]
                obstruction = {min(block) for block in insufficient}
                need(all(block & obstruction for block in insufficient),
                     'Packing-count sharpness control failed')
    for r in range(1, 5):
        states = [set(s) for s in itertools.combinations(range(2 * r + 1), r)]
        for left in states:
            for right in states:
                support = len(left ^ right)
                need(support <= 2 * r and 2 * support <= 4 * r,
                     'Helper replacement-cost bound failed')
                helper_pair_tests += 1
        need(max(len(left ^ right) for left in states for right in states) == 2 * r,
             'Maximum helper replacement bound was not attained')
    for m in range(1, 13):
        wt = weight_table((m,))
        gap = boundary((0, 0), (1, 1), wt)[0]
        need(gap == m and 1 + gap > m, 'Missing neutral interface-gap control failed')
    return {'forbidden_U_deletion_cases': deletion_tests,
            'complete_r_selected_helper_pairs': helper_pair_tests,
            'packing_scope': 'r=1..3, t=0..3, all deletions of at most t in (t+1) disjoint r-sets',
            'helper_pair_scope': 'r=1..4, all ordered r-subset pairs of a 2r+1-slot bank',
            'missing_boundary_cases': 12}


def verify_sources():
    pins = json.loads((HERE / 'SOURCEPINS.json').read_text())
    for row in pins['local_inputs']:
        path = HERE / row['snapshot']
        need(path.is_file(), f'Missing source snapshot: {path}')
        need(digest(path) == row['sha256'], f'Source hash mismatch: {path}')
        need(path.stat().st_size == row['bytes'], f'Source size mismatch: {path}')
    return {'snapshots_verified': len(pins['local_inputs']),
            'sourcepins_sha256': digest(HERE / 'SOURCEPINS.json')}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default=None)
    parser.add_argument('--section', choices=('all', 'graphs', 'intervals', 'quotients', 'counterexamples', 'supply'), default='all')
    args = parser.parse_args()
    start = time.monotonic()
    receipt = {'status': 'PASS', 'python_optimization': sys.flags.optimize,
               'checker_sha256': digest(Path(__file__)), 'source_inputs': verify_sources()}
    jobs = [('graphs', graph_controls), ('intervals', interval_controls),
            ('quotients', quotient_controls), ('counterexamples', counterexamples),
            ('supply', supply_controls)]
    for name, job in jobs:
        if args.section in ('all', name):
            receipt[name] = job()
            print(f'{name}: PASS', file=sys.stderr, flush=True)
    receipt['explicit_checks'] = CHECKS
    receipt['elapsed_seconds'] = round(time.monotonic() - start, 3)
    serialized = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    if args.output:
        (HERE / args.output).write_text(serialized)
    print(serialized, end='')


if __name__ == '__main__':
    main()
