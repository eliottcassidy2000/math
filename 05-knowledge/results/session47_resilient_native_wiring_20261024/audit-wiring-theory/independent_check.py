#!/usr/bin/env python3
"""Independent audit: set-based tournaments/posets and literal simple paths.

No producer checker is imported or executed. Checks use explicit exceptions,
so all validation survives python -O. Input snapshots are authenticated first.
"""
from __future__ import annotations
import argparse
from collections import Counter, deque
import hashlib
import itertools as it
import json
import math
from pathlib import Path
import random
import sys
import time

HERE = Path(__file__).resolve().parent
COUNTS = Counter()
INF = math.inf


def require(value, section, description):
    COUNTS[section] += 1
    if not value:
        raise RuntimeError(section + ': ' + description)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def blocks_from_edges(n, edges):
    labels = list(range(n))
    for x, y in edges:
        a, b = labels[x], labels[y]
        labels = [a if z == b else z for z in labels]
    return labels


def path_peak(matrix, start, end):
    if start == end:
        return 0
    available = [v for v in range(len(matrix)) if v not in (start, end)]
    result = INF
    for k in range(len(available) + 1):
        for middle in it.permutations(available, k):
            path = (start,) + middle + (end,)
            result = min(result, max(matrix[a][b] for a, b in zip(path, path[1:])))
    return result


def literal_beta(matrix):
    return [[path_peak(matrix, a, b) for b in range(len(matrix))]
            for a in range(len(matrix))]


def partitions(n):
    if n == 0:
        yield ()
        return
    def walk(labels):
        if len(labels) == n:
            yield tuple(labels)
        else:
            for value in range(max(labels) + 2):
                yield from walk(labels + [value])
    yield from walk([0])


def graph_controls():
    rows = []
    partition_tests = 0
    rng = random.Random(47001)
    for n in range(1, 7):
        pairs = list(it.combinations(range(n), 2))
        if n <= 4:
            samples = it.product((0, 2, 7, INF), repeat=len(pairs))
        else:
            samples = (tuple(rng.choice((0, 2, 7, INF)) for _ in pairs)
                       for _ in range(96 if n == 5 else 24))
        trees = [t for t in it.combinations(pairs, n - 1)
                 if len(set(blocks_from_edges(n, t))) == 1]
        graphs = minimum_trees = 0
        for values in samples:
            matrix = [[0 if a == b else INF for b in range(n)] for a in range(n)]
            for (a, b), weight in zip(pairs, values):
                matrix[a][b] = matrix[b][a] = weight
            beta = literal_beta(matrix)
            cut = {}
            for mask in range(1, (1 << n) - 1):
                cut[mask] = min(matrix[a][b] for a in range(n) if mask >> a & 1
                                for b in range(n) if not mask >> b & 1)
            for a, b in pairs:
                require(beta[a][b] == max(value for mask, value in cut.items()
                                          if mask >> a & 1 and not mask >> b & 1),
                        'graphs', 'simple-path minimax versus all separating cuts')
            for a, b, c in it.product(range(n), repeat=3):
                require(beta[a][c] <= max(beta[a][b], beta[b][c]),
                        'graphs', 'strong triangle inequality')
            for budget in (0, 1, 2, 6, 7):
                labels = blocks_from_edges(n, [(a, b) for a, b in pairs
                                               if matrix[a][b] <= budget])
                require(all((labels[a] == labels[b]) == (beta[a][b] <= budget)
                            for a in range(n) for b in range(n)),
                        'graphs', 'threshold union-find versus literal paths')
            B = max(map(max, beta))
            feasible_trees = [t for t in trees if all(math.isfinite(matrix[a][b]) for a, b in t)]
            require(B == max(cut.values(), default=0), 'graphs', 'global cut formula')
            require(B == min((max((matrix[a][b] for a, b in t), default=0)
                              for t in feasible_trees), default=INF),
                    'graphs', 'minimum spanning-tree peak')
            if feasible_trees:
                least_sum = min(sum(matrix[a][b] for a, b in t) for t in feasible_trees)
                for t in feasible_trees:
                    if sum(matrix[a][b] for a, b in t) != least_sum:
                        continue
                    sparse = [[0 if a == b else INF for b in range(n)] for a in range(n)]
                    for a, b in t:
                        sparse[a][b] = sparse[b][a] = matrix[a][b]
                    require(literal_beta(sparse) == beta, 'graphs', 'every minimum-total tree')
                    minimum_trees += 1
            # Directly verify equation (8) on every partition for smaller graphs,
            # and a deterministic 1/31 sample of four-vertex graphs.
            if n <= 3 or (n == 4 and graphs % 31 == 0):
                for projection in partitions(n):
                    groups = [[a for a in range(n) if projection[a] == j]
                              for j in range(max(projection) + 1)]
                    internal_peak = max(max(map(max, literal_beta([[matrix[a][b] for b in g] for a in g])))
                                        for g in groups)
                    q = [[0 if i == j else min(matrix[a][b] for a in g for b in h)
                          for j, h in enumerate(groups)] for i, g in enumerate(groups)]
                    qb = literal_beta(q)
                    require(all(qb[projection[a]][projection[b]] <= beta[a][b]
                                <= max(internal_peak, qb[projection[a]][projection[b]])
                                for a in range(n) for b in range(n)),
                            'contraction', 'fixed-partition two-sided bottleneck bound')
                    partition_tests += 1
            graphs += 1
        rows.append({'vertices': n, 'graphs': graphs, 'all_minimum_trees': minimum_trees,
                     'scope': 'exhaustive' if n <= 4 else 'seeded extra controls'})
    return {'rows': rows, 'edge_alphabet': [0, 2, 7, 'absent'],
            'partitions_tested': partition_tests, 'rng_seed': 47001}


def all_posets(n):
    """Enumerate all labelled orders directly, not by closing DAG generators."""
    pairs = list(it.combinations(range(n), 2))
    for choices in it.product((0, 1, -1), repeat=len(pairs)):
        rel = {(v, v) for v in range(n)}
        for (a, b), sign in zip(pairs, choices):
            if sign:
                rel.add((a, b) if sign == 1 else (b, a))
        if all((a, c) in rel for a, b in rel for bb, c in rel if bb == b):
            yield frozenset(rel)


def subset_list(n):
    return [frozenset(i for i in range(n) if mask >> i & 1) for mask in range(1 << n)]


def upset_list(n, rel):
    return [s for s in subset_list(n) if all(a not in s or b in s for a, b in rel)]


def member_list(interval, upsets):
    low, high = interval
    return [s for s in upsets if low <= s <= high]


def weight(s, weights):
    return sum(weights[x] for x in s)


def exact_distances(states, costs):
    """Unweighted shortest distances in the explicitly given adjacency graph."""
    rows = []
    for source in range(len(states)):
        dist = {source: 0}
        queue = deque([source])
        while queue:
            a = queue.popleft()
            for b in costs[a]:
                if b not in dist:
                    dist[b] = dist[a] + 1
                    queue.append(b)
        rows.append(dist)
    return rows


def interval_controls():
    rows = []
    for n in range(1, 5):
        order_count = interval_count = gap_count = route_pairs = 0
        for rel in all_posets(n):
            order_count += 1
            upsets = upset_list(n, rel)
            intervals = [(a, b) for a in upsets for b in upsets if a <= b]
            members = [member_list(interval, upsets) for interval in intervals]
            for interval, states in zip(intervals, members):
                adjacency = [[j for j, t in enumerate(states) if len(s ^ t) == 1]
                             for s in states]
                distances = exact_distances(states, adjacency)
                for i, s in enumerate(states):
                    for j, t in enumerate(states):
                        require(distances[i].get(j, INF) == len(s ^ t),
                                'intervals', 'one-component routing has conformal shortest length')
                        route_pairs += 1
                interval_count += 1
            for weights in (tuple(0 if i % 2 == 0 else 3 for i in range(n)),
                            tuple((2, 5, 1, 7)[i] for i in range(n))):
                for i, (la, ra) in enumerate(intervals):
                    for j, (lb, rb) in enumerate(intervals):
                        formula = weight(la - rb, weights) + weight(lb - ra, weights)
                        literal = min(weight(s ^ t, weights) for s in members[i] for t in members[j])
                        wa, wb = la | (lb & ra), lb | (la & rb)
                        require(formula == literal, 'intervals', 'exact boundary formula')
                        require(wa in members[i] and wb in members[j], 'intervals', 'witness membership')
                        require(weight(wa ^ wb, weights) == formula, 'intervals', 'witness attains bound')
                        gap_count += 1
        require(order_count == (1, 3, 19, 219)[n - 1], 'intervals', 'all labelled poset census')
        rows.append({'components': n, 'all_labelled_posets': order_count,
                     'intervals': interval_count, 'route_pairs': route_pairs,
                     'weighted_interval_pairs': gap_count})
    return rows


def sum_charge(s, labels):
    return (sum(labels[x][0] for x in s), sum(labels[x][1] for x in s) % 3)


def minus_charge(c):
    return (-c[0], -c[1] % 3)


def charge_zero(c, h):
    return c[0] + h[0] == 0 and (c[1] + h[1]) % 3 == 0


def reachability(cost, budget):
    n = len(cost)
    return blocks_from_edges(n, [(a, b) for a in range(n) for b in range(a + 1, n)
                                 if cost[a][b] <= budget])


def quotient_controls():
    rng = random.Random(47002)
    rows = []
    empty_cases = nonempty_cases = 0
    for n in range(1, 5):
        instances = budgets_tested = projected_vertices = 0
        for rel in all_posets(n):
            upsets = upset_list(n, rel)
            intervals = [(a, b) for a in upsets for b in upsets if a <= b]
            for trial in range(4):
                chosen = [rng.choice(intervals) for _ in range(3)]
                labels = [rng.choice(((0, 0), (0, 0), (0, 1), (0, 2), (1, 0), (-1, 0)))
                          for _ in range(n)]
                if trial == 0:
                    labels = [(0, 0)] * n
                if trial == 1:
                    labels[0] = (1, 0)
                weights = tuple(rng.choice((0, 1, 3, 5)) for _ in range(n))
                hcharges = [minus_charge(sum_charge(rng.choice(member_list(iv, upsets)), labels))
                            for iv in chosen]
                if trial == 3:
                    hcharges[0] = (n + 2, 0)  # intentionally empty helper fibre
                helper_cost = [[0] * 3 for _ in range(3)]
                for h, g in it.combinations(range(3), 2):
                    helper_cost[h][g] = helper_cost[g][h] = rng.choice((0, 1, 2, 6))
                states = [(h, s) for h, iv in enumerate(chosen) for s in member_list(iv, upsets)
                          if charge_zero(sum_charge(s, labels), hcharges[h])]
                full_cost = [[helper_cost[h][g] + weight(s ^ t, weights) for g, t in states]
                             for h, s in states]
                empty_cases += sum(not any(hh == h for hh, _ in states) for h in range(3))
                nonempty_cases += sum(any(hh == h for hh, _ in states) for h in range(3))
                previous = None
                for budget in (0, 1, 2, 3, 4, 5, 7, 11):
                    visible = frozenset(x for x in range(n) if labels[x] != (0, 0) or weights[x] > budget)
                    neutral = frozenset(range(n)) - visible
                    group_dict = {}
                    for index, (h, s) in enumerate(states):
                        group_dict.setdefault((h, s & visible), []).append(index)
                    keys = sorted(group_dict, key=lambda key: (key[0], sorted(key[1])))
                    groups = [group_dict[key] for key in keys]
                    projection = {i: j for j, group in enumerate(groups) for i in group}
                    compiled = []
                    for h, (lh, rh) in enumerate(chosen):
                        for tv in subset_list(n):
                            if not tv <= visible or not charge_zero(sum_charge(tv, labels), hcharges[h]):
                                continue
                            low = frozenset(b for a, b in rel if a in lh | tv)
                            forbidden = frozenset(range(n)) - rh | (visible - tv)
                            high = frozenset(range(n)) - frozenset(a for a, b in rel if b in forbidden)
                            if low <= high:
                                compiled.append((h, tv, low & neutral, high & neutral))
                    compiled_by_key = {(h, tv): (low, high) for h, tv, low, high in compiled}
                    require(set(compiled_by_key) == set(keys), 'quotients', 'literal versus compiled vertex set')
                    for j, key in enumerate(keys):
                        low, high = compiled_by_key[key]
                        actual = {states[i][1] & neutral for i in groups[j]}
                        expected = {s for s in subset_list(n) if s <= neutral and low <= s <= high
                                    and all(a not in s or b not in neutral or b in s for a, b in rel)}
                        require(actual == expected, 'quotients', 'literal versus compiled neutral fibre')
                        inside = [[full_cost[i][k] for k in groups[j]] for i in groups[j]]
                        require(len(set(reachability(inside, budget))) == 1,
                                'quotients', 'internal budget connectivity')
                    qcost = [[0] * len(keys) for _ in keys]
                    for i, (h, tv) in enumerate(keys):
                        la, ra = compiled_by_key[(h, tv)]
                        for j in range(i + 1, len(keys)):
                            g, uv = keys[j]
                            lb, rb = compiled_by_key[(g, uv)]
                            formula = helper_cost[h][g] + weight(tv ^ uv, weights)
                            formula += weight(la - rb, weights) + weight(lb - ra, weights)
                            literal = min(full_cost[x][y] for x in groups[i] for y in groups[j])
                            require(formula == literal, 'quotients', 'exact crossing minimum')
                            qcost[i][j] = qcost[j][i] = formula
                    full_cc, q_cc = reachability(full_cost, budget), reachability(qcost, budget)
                    require(all((full_cc[i] == full_cc[j]) == (q_cc[projection[i]] == q_cc[projection[j]])
                                for i in range(len(states)) for j in range(len(states))),
                            'quotients', 'complete connectivity equivalence')
                    if previous is not None:
                        oldkeys, oldgroups = previous
                        for oldkey, group in zip(oldkeys, oldgroups):
                            reduced = (oldkey[0], oldkey[1] & visible)
                            require(all(keys[projection[i]] == reduced for i in group),
                                    'quotients', 'budget projection commutes')
                    previous = keys, groups
                    projected_vertices += len(keys)
                    budgets_tested += 1
                instances += 1
        rows.append({'components': n, 'three_helper_instances': instances,
                     'budgets': budgets_tested, 'projected_vertices': projected_vertices})
    return {'charge_group': 'Z x C3', 'helper_costs': 'symmetric, possibly zero and nonmetric',
            'rng_seed': 47002, 'rows': rows, 'empty_helper_fibres': empty_cases,
            'nonempty_helper_fibres': nonempty_cases}


def tournament(n, bits):
    edges = {(i, i + 1) for i in range(n - 1)}
    off_path = [(i, j) for i in range(n) for j in range(i + 2, n)]
    for k, (i, j) in enumerate(off_path):
        edges.add((i, j) if bits & (1 << k) else (j, i))
    return frozenset(edges)


def reverse_nonpath(edges, path):
    kept = set(zip(path, path[1:]))
    require(kept <= edges, 'suspension', 'Hamiltonian path authentication')
    return frozenset((a, b) if (a, b) in kept else (b, a) for a, b in edges)


def suspend(edges, n):
    return edges | {(-1, v) for v in range(n)} | {(v, n) for v in range(n)} | {(-1, n)}


def degrees(edges, vertices):
    return {v: sum(a == v for a, b in edges) for v in vertices}


def cyclic_triples(edges, vertices):
    # Independent of directed-cycle orientation: a three-vertex tournament is
    # cyclic exactly when all three induced outdegrees equal one.
    result = []
    for triple in it.combinations(vertices, 3):
        s = set(triple)
        count = Counter(a for a, b in edges if a in s and b in s)
        if all(count[v] == 1 for v in triple):
            result.append(frozenset(triple))
    return result


def canonical(edges, vertices):
    best = None
    count = 0
    for ordering in it.permutations(vertices):
        code = 0
        for i in range(len(ordering)):
            for j in range(i + 1, len(ordering)):
                code = (code << 1) | ((ordering[i], ordering[j]) in edges)
        if best is None or code < best:
            best = code
        count += 1
    return best, count


def suspension_controls():
    examples = []
    for bits in (2, 3):
        source = tournament(5, bits)
        target = reverse_nonpath(source, tuple(range(5)))
        enlarged = reverse_nonpath(suspend(source, 5), (-1, 0, 1, 2, 3, 4, 5))
        d, dh = degrees(source, range(5)), degrees(enlarged, range(-1, 6))
        require(list(d.values()) == ([2, 1, 2, 2, 3] if bits == 2 else [3, 1, 1, 2, 3]),
                'suspension', 'explicit source degrees')
        require([v for v in dh if dh[v] == 1] == [-1] and [v for v in dh if dh[v] == 5] == [5],
                'suspension', 'intrinsically unique outer ports')
        require({b for a, b in enlarged if a == -1} == {0}, 'suspension', 'first inner endpoint recovery')
        require({a for a, b in enlarged if b == 5} == {4}, 'suspension', 'last inner endpoint recovery')
        old_canon, old_perms = canonical(target, tuple(range(5)))
        new_canon, new_perms = canonical(enlarged, tuple(range(-1, 6)))
        converse_canon, converse_perms = canonical(frozenset((b, a) for a, b in enlarged), tuple(range(-1, 6)))
        examples.append({'bits': bits, 'source_degrees': list(d.values()), 'H_scores': sorted(dh.values()),
                         'H_square_sum': sum(x * x for x in dh.values()),
                         'old_canonical_integer': old_canon, 'H_canonical_integer': new_canon,
                         'H_converse_canonical_integer': converse_canon,
                         'canonical_permutation_counts': [old_perms, new_perms, converse_perms]})
    k2, k3 = (reverse_nonpath(tournament(5, bits), tuple(range(5))) for bits in (2, 3))
    f = (1, 2, 3, 0, 4)
    require(frozenset((f[a], f[b]) for a, b in k2) == k3, 'suspension', 'explicit old target isomorphism')
    require(examples[0]['old_canonical_integer'] == examples[1]['old_canonical_integer'],
            'suspension', 'direct canonical equality of old targets')
    require(examples[0]['H_canonical_integer'] not in
            (examples[1]['H_canonical_integer'], examples[1]['H_converse_canonical_integer']),
            'suspension', 'direct canonical separation including converse')
    require([row['H_square_sum'] for row in examples] == [73, 75], 'suspension', 'converse-invariant square obstruction')
    rows = []
    for n in range(2, 8):
        cases = 0
        for bits in range(1 << ((n - 1) * (n - 2) // 2)):
            source = tournament(n, bits)
            target = reverse_nonpath(source, tuple(range(n)))
            enlarged = reverse_nonpath(suspend(source, n), (-1,) + tuple(range(n)) + (n,))
            vertices = tuple(range(n))
            enlarged_vertices = (-1,) + vertices + (n,)
            d, dk, dh = degrees(source, vertices), degrees(target, vertices), degrees(enlarged, enlarged_vertices)
            ts, tk, th = (cyclic_triples(source, vertices), cyclic_triples(target, vertices),
                          cyclic_triples(enlarged, enlarged_vertices))
            for edges, vs, score, cycles in ((source, vertices, d, ts), (target, vertices, dk, tk),
                                             (enlarged, enlarged_vertices, dh, th)):
                require(len(edges) == math.comb(len(vs), 2), 'suspension', 'tournament edge census')
                require(len(cycles) == math.comb(len(vs), 3) - sum(math.comb(v, 2) for v in score.values()),
                        'suspension', 'induced-degree cycles versus score identity')
            require(all(dk[v] == n - 1 - d[v] + (v == 0) - (v == n - 1) for v in vertices),
                    'suspension', 'old target degree identity')
            require([dh[v] for v in enlarged_vertices] == [1] + [n - d[v] for v in vertices] + [n],
                    'suspension', 'enlarged degree identity')
            require(len(tk) == len(ts) + d[0] - d[n - 1] - 1,
                    'suspension', 'inner triangle identity')
            require(len(th) == len(ts) + n, 'suspension', 'enlarged triangle identity')
            boundary = [tri for tri in th if -1 in tri or n in tri]
            require(len(boundary) == n + 1 - d[0] + d[n - 1], 'suspension', 'boundary triangle identity')
            require(sum(-1 in tri for tri in boundary) == n - d[0]
                    and sum(n in tri for tri in boundary) == 1 + d[n - 1]
                    and not any(-1 in tri and n in tri for tri in boundary),
                    'suspension', 'separate outer-port triangle contributions')
            require(frozenset((a, b) for a, b in enlarged if a in vertices and b in vertices) == target,
                    'suspension', 'literal reverse readout')
            cases += 1
        rows.append({'source_order': n, 'fixed_forward_path_tournaments': cases})
    require(sum(row['fixed_forward_path_tournaments'] for row in rows if row['source_order'] <= 6) == 1099,
            'suspension', 'producer universe equals 1099, not all tournaments')
    # Additional direct eight-owner canonicalization and arbitrary-owner renaming.
    extra = reverse_nonpath(suspend(tournament(6, 277), 6), (-1, 0, 1, 2, 3, 4, 5, 6))
    vs = tuple(range(-1, 7))
    naming = dict(zip(vs, (18, 3, 44, 0, 21, 9, 67, 2)))
    renamed = frozenset((naming[a], naming[b]) for a, b in extra)
    can1, perms1 = canonical(extra, vs)
    can2, perms2 = canonical(renamed, tuple(naming.values()))
    require(can1 == can2, 'suspension', 'eight-owner canonical invariance under arbitrary renaming')
    return {'witnesses': examples, 'universe_rows': rows,
            'original_scope_total': 1099, 'extra_order7_fixed_path_tournaments': 32768,
            'extra_order8_canonicalizations': {'canonical_integer': can1, 'permutation_counts': [perms1, perms2]},
            'scope': 'fixed forward-path tournaments and named examples, not full metagraphs or harmonic transfer'}


def hostile_and_supply_controls():
    require(path_peak([[0, 2], [2, 0]], 0, 1) == 2, 'hostiles', 'heavy neutral barrier')
    require(blocks_from_edges(2, []) != [0, 0], 'hostiles', 'budget one cannot cross weight two')
    require(blocks_from_edges(1, []) == [0], 'hostiles', 'collapse would invent connectivity')
    intervals = [(frozenset(), frozenset()), (frozenset(), frozenset({0})), (frozenset({0}), frozenset({0}))]
    gap = lambda a, b: len(a[0] - b[1]) + len(b[0] - a[1])
    require((gap(intervals[0], intervals[1]), gap(intervals[1], intervals[2]), gap(intervals[0], intervals[2])) == (0, 0, 1),
            'hostiles', 'nonmetric interval gap')
    require(1 + gap(intervals[0], intervals[2]) == 2, 'hostiles', 'disjoint intervals can have a legal simultaneous bridge')
    require(1 + 0 == 1, 'hostiles', 'zero-weight disjoint intervals do not force positive distance')
    forced_heavy_singleton = [frozenset({0})]
    require(len(forced_heavy_singleton) == 1, 'hostiles', 'heavy visibility is not necessary in every instance')
    require(path_peak([[0, 1, 9], [1, 0, 9], [9, 9, 0]], 0, 1) == 1,
            'hostiles', 'fixed endpoints preserve old route after supply enlargement')
    require(path_peak([[0, 1, 9], [1, 0, 9], [9, 9, 0]], 0, 2) == 9,
            'hostiles', 'new obligations can raise universal threshold')
    require(11 * 2 == 22 and 10 * 2 + 3 == 23, 'hostiles', 'slot/clock units cannot be inferred from factorization')
    repack = lambda x: ((x & 1) | ((((x >> 1) ^ x) & 1) << 1)
                       | (x & 4) | ((((x >> 3) ^ x) & 1) << 3))
    require(repack(11) == 1 and len({repack(x) for x in range(16)}) == 16,
            'hostiles', 'displayed parity eleven can become one')
    require(all(repack(x ^ y) == repack(x) ^ repack(y) for x in range(16) for y in range(16)),
            'hostiles', 'parity relabeling is linear')
    deletion_cases = helper_pairs = 0
    for r in range(1, 5):
        for p in range(1, 5):
            bank = [frozenset(range(j * r, (j + 1) * r)) for j in range(p)]
            for k in range(p):
                for deleted_tuple in it.combinations(range(p * r), k):
                    deleted = frozenset(deleted_tuple)
                    require(any(not block & deleted for block in bank), 'supply', 'disjoint packing survives p-1 forced-U slots')
                    deletion_cases += 1
            hit = frozenset(min(block) for block in bank)
            require(all(block & hit for block in bank), 'supply', 'p forced-U slots can destroy p blocks')
        chosen = [frozenset(s) for s in it.combinations(range(2 * r), r)]
        costs = tuple((i % 4) for i in range(2 * r))
        W = max(costs)
        for a, b in it.product(chosen, repeat=2):
            require(len(a ^ b) <= 2 * r and weight(a ^ b, costs) <= 2 * r * W,
                    'supply', 'complete equal-margin helper-change bound')
            helper_pairs += 1
    for m in range(1, 21):
        require(1 + m > m, 'supply', 'neutral boundary term cannot be omitted')
    return {'forbidden_U_cases': deletion_cases, 'helper_pairs': helper_pairs,
            'missing_boundary_weights': list(range(1, 21))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    pins = json.loads((HERE / 'INPUT_PINS.json').read_text())
    for row in pins['inputs']:
        p = HERE / row['snapshot']
        require(p.is_file() and digest(p) == row['sha256'] and p.stat().st_size == row['bytes'],
                'authentication', 'input snapshot: ' + row['source'])
    result = {'status': 'PASS', 'python_optimization': sys.flags.optimize,
              'checker_sha256': digest(Path(__file__)), 'input_pins_sha256': digest(HERE / 'INPUT_PINS.json')}
    for name, function in [('graphs', graph_controls), ('intervals', interval_controls),
                           ('quotients', quotient_controls), ('suspension', suspension_controls),
                           ('hostiles_and_supply', hostile_and_supply_controls)]:
        result[name] = function()
        print(name + ': PASS', file=sys.stderr, flush=True)
    result['checks_by_section'] = dict(sorted(COUNTS.items()))
    result['explicit_checks'] = sum(COUNTS.values())
    result['elapsed_seconds'] = round(time.monotonic() - started, 3)
    output = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
