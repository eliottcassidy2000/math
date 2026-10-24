#!/usr/bin/env python3
"""Literal controls for charge/count and exact boundary-cost tree messages."""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import argparse
import json

from compiler import Group, TreeCarrier

HERE = Path(__file__).resolve().parent
COUNT = 0


def need(ok, why):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(why)


def relation(mask, q=2):
    values = list(product(range(q), repeat=2))
    return [t for k, t in enumerate(values) if mask >> k & 1]


def count_control(carrier, tables):
    exact = carrier.literal(tables)
    answer, _ = carrier.count(tables)
    need({q: n for q, (n, w) in answer.items()} == {q: len(v) for q, v in exact.items()},
         'count messages equal independent complete global enumeration')
    for q, (n, witness) in answer.items():
        need(witness in exact[q], 'section is a literal complete record at stated charge')
    return exact


def pair_control(carrier, first, second):
    aa, bb = carrier.literal(first), carrier.literal(second)
    answers, _ = carrier.bridge(first, second)
    need(set(answers) == set(product(aa, bb)), 'all and only nonempty paired charge fibres')
    for (q, r), (cost, a, b) in answers.items():
        expected = min(sum(carrier.costs[i][x[i], y[i]] for i in carrier.slots)
                       for x in aa[q] for y in bb[r])
        need(cost == expected, 'bridge minimum equals global pair enumeration')
        need(a in aa[q] and b in bb[r], 'minimum witnesses are actual feasible endpoints')
        need(sum(carrier.costs[i][a[i], b[i]] for i in carrier.slots) == cost,
             'actual endpoint clock bill equals reported minimum')


def binary_carrier(bags, parents, group, weights):
    slots = sorted(set().union(*map(set, bags)))
    domains = {i: (0, 1) for i in slots}
    charges = {i: {b: tuple(b*w for w in weights[i]) for b in (0, 1)} for i in slots}
    costs = {i: {(a, b): (i+1)*(a != b) for a, b in product((0, 1), repeat=2)} for i in slots}
    return TreeCarrier(bags, parents, domains, charges, costs, group)


def finite_controls():
    cases = Counter()
    # Every local binary relation on both edges, every binary charge vector.
    for weights in product((0, 1), repeat=3):
        c = binary_carrier([(0, 1), (1, 2)], [-1, 0], Group((2,)),
                           {i: (weights[i],) for i in range(3)})
        for a, b in product(range(16), repeat=2):
            tables = [relation(a), relation(b)]
            count_control(c, tables)
            cases['all_path_three_count_instances'] += 1
            if weights in ((1, 1, 1), (1, 0, 1)):
                other = [relation(15-b), relation((a*5+3) % 16)]
                pair_control(c, tables, other)
                cases['path_three_two_arrival_instances'] += 1
    # Branched separators: all choices among full/equal/unequal/implication.
    masks = (15, 9, 6, 13)
    for first in product(masks, repeat=3):
        c = binary_carrier([(0,), (0, 1), (0, 2), (2, 3)], [-1, 0, 0, 2],
                           Group((0, 3)), {i: (1, i+1) for i in range(4)})
        tables = [[(0,), (1,)]] + [relation(x) for x in first]
        count_control(c, tables)
        for second in product(masks, repeat=3):
            other = [[(0,), (1,)]] + [relation(x) for x in second]
            pair_control(c, tables, other)
            cases['branched_Z_times_C3_pair_instances'] += 1
    # A shared slot must contribute once, not once per bag occurrence.
    c = binary_carrier([(0, 1), (1, 2)], [-1, 0], Group((0,)),
                       {i: (1,) for i in range(3)})
    answer, _ = c.count([relation(15), relation(15)])
    need({q[0]: n for q, (n, w) in answer.items()} == {0: 1, 1: 3, 2: 3, 3: 1},
         'separator slot counted exactly once in integer margin')
    for bad_bags, bad_parent in [([(0,), (1,), (0, 1)], [-1, 0, 1]),
                                 ([(0, 1), (0, 1)], [-1, -1])]:
        try:
            binary_carrier(bad_bags, bad_parent, Group((2,)), {0: (1,), 1: (1,)})
        except ValueError:
            need(True, 'invalid carrier rejected')
        else:
            need(False, 'invalid carrier incorrectly accepted')
    # Independent bridge minima do not compose through a disconnected fibre.
    c = TreeCarrier([(0,)], [-1], {0: (0, 1)}, {0: {0: (), 1: ()}},
                    {0: {(x, y): 2*(x != y) for x, y in product(range(2), repeat=2)}}, Group(()))
    aa, bb, cc = [[(0,)]], [[(0,), (1,)]], [[(1,)]]
    ab, _ = c.bridge(aa, bb)
    bc, _ = c.bridge(bb, cc)
    ac, _ = c.bridge(aa, cc)
    need([t[(), ()][0] for t in (ab, bc, ac)] == [0, 0, 2],
         'minimum fibre distance is not a metric and cannot certify composed routes')
    need(ab[(), ()][2] != bc[(), ()][1], 'two independent bridge minima choose different middle records')
    return dict(cases)


def coset_controls():
    group = Group((2, 3))
    elements = list(product(range(2), range(3)))
    subgroups = []
    for mask in range(1 << 6):
        h = {x for i, x in enumerate(elements) if mask >> i & 1}
        if group.zero in h and all(group.add(a, b) in h for a, b in product(h, repeat=2)):
            subgroups.append(h)
    need(len(subgroups) == 4, 'complete finite subgroup universe C2xC3')
    for h, k, a, b in product(subgroups, subgroups, elements, elements):
        left = {group.add(group.add(a, x), group.add(b, y)) for x in h for y in k}
        right = {group.add(group.add(a, b), z) for z in {group.add(x, y) for x in h for y in k}}
        need(left == right, 'conditional coset convolution law')
    # Phase-dependent actual record IDs carry private charge directions.
    domains = {0: tuple(range(3)), 1: tuple(range(6)), 2: tuple(range(9))}
    charges = {0: {x: (0, 0) for x in domains[0]},
               1: {x: (x//3, x % 3) for x in domains[1]},
               2: {x: (0, x % 3) for x in domains[2]}}
    costs = {i: {(a, b): (i+1)*(a != b) for a, b in product(d, repeat=2)}
             for i, d in domains.items()}
    c = TreeCarrier([(0,), (0, 1), (0, 2)], [-1, 0, 0], domains, charges, costs, group)
    local_x = [(z, x) for z in domains[0] for x in domains[1] if x % 3 == z]
    local_y = [(z, y) for z in domains[0] for y in domains[2] if y//3 == z]
    for z in range(3):
        tables = [[(z,)], local_x, local_y]
        exact = count_control(c, tables)
        need(set(exact) == set(elements) and all(len(v) == 1 for v in exact.values()),
             'separate C2 and C3 directions give full coverage at EACH shared phase')
        for zz in range(3):
            pair_control(c, tables, [[(zz,)], local_x, local_y])
    # The minimal incompatible-composition hostile.
    c = binary_carrier([(0, 1)], [-1], Group((2,)), {0: (1,), 1: (1,)})
    need(set(c.literal([relation(9)])) == {(0,)},
         'individually onto contributions coupled equal have only total0')
    need({(x+y) % 2 for x, y in product(range(2), repeat=2)} == {0, 1},
         'discarding the owner relation creates a false charge1 witness')
    return {'subgroup_count': len(subgroups), 'conditional_positive_phases': 3,
            'minimal_hostile': {'group': 'C2', 'slots': 2, 'factor': 'x=y', 'only_total': 0}}


def native_controls():
    base_pairs = list(combinations(range(8), 2))
    base = [[0]*8 for _ in range(8)]
    for j, (a, b) in enumerate(base_pairs):
        base[a][b] = (29558785 >> (27-j)) & 1
        base[b][a] = 1-base[a][b]
    p0, p1 = (3, 0, 6, 7, 1, 5, 4, 2), (3, 1, 6, 7, 0, 5, 4, 2)
    rows = []
    for r in range(1, 4):
        n, start = 8*r+3, 8*r
        a = [[0]*n for _ in range(n)]
        for u, v in combinations(range(n), 2):
            if u//8 == v//8 and v < start:
                value = base[u % 8][v % 8]
            elif u >= start:
                value = int((v-u) % 3 == 1)
            else:
                value = 1
            a[u][v], a[v][u] = value, 1-value
        chords = [(u, v) for u in range(n) for v in range(u+2, n)]
        refl = [chords.index((n-1-v, n-1-u)) for u, v in chords]
        full = (1 << len(chords))-1
        all_records, resources = {}, set()
        for address in range(1 << r):
            prefix = tuple(8*j+x for j in range(r) for x in (p1 if address >> j & 1 else p0))
            paths = [prefix + tuple(start+(k+j) % 3 for j in range(3)) for k in range(3)]
            words = []
            for p in paths:
                need(all(a[u][v] for u, v in zip(p, p[1:])), 'literal native full path')
                words.append(sum(a[p[u]][p[v]] << j for j, (u, v) in enumerate(chords)))
            need(len(set(words)) == 1, 'three C3 kernel lifts use ONE normalized resource')
            word = words[0]
            reflected = sum(((word >> j) & 1) << k for j, k in enumerate(refl))
            resource = min(word, reflected, full ^ word, full ^ reflected)
            need(resource not in resources, 'prefix addresses are distinct actual normalized units')
            resources.add(resource)
            for phase in (1, 2):
                rotate = list(range(start)) + [start+(j+phase) % 3 for j in range(3)]
                need(all(a[u][v] == a[rotate[u]][rotate[v]] for u, v in combinations(range(n), 2)),
                     'literal final-C3 source automorphism fixes all prefix owners')
                need(tuple(rotate[x] for x in paths[0]) == paths[phase], 'actual phase transport of path')
            all_records[address] = paths
        for m in range(1, (1 << r)+1):
            bags = [(0,)] + [(j-1, j) for j in range(1, m)]
            parents = [-1] + list(range(m-1))
            domains = {j: (0, 1, 2) for j in range(m)}
            charges = {j: {k: () for k in range(3)} for j in range(m)}
            costs = {j: {(k, ell): sum(x != y for x, y in zip(all_records[j][k], all_records[j][ell]))
                         for k, ell in product(range(3), repeat=2)} for j in range(m)}
            c = TreeCarrier(bags, parents, domains, charges, costs, Group(()))
            equality = [[(k, k) for k in range(3)] for _ in range(m-1)]
            tables0, tables1 = [[(0,)]] + equality, [[(1,)]] + equality
            pair_control(c, tables0, tables1)
            answer, _ = c.bridge(tables0, tables1)
            need(answer[(), ()][0] == 3*m, 'native width-one exact boundary cost3m')
            rows.append({'r': r, 'source_order': n, 'slots': m, 'separator_width': 1,
                         'exact_clock_cost': 3*m, 'resource_count': m})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    out = {'status': 'FINITE-EXACT', 'finite_universes': finite_controls(),
           'coset_controls': coset_controls(), 'native_boundary_family': native_controls(),
           'checks': COUNT,
           'source_hashes': {p.name: sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                             for p in (HERE/'compiler.py', HERE/'check.py')}}
    path = HERE/'receipt.json'
    if args.write:
        path.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    else:
        need(out == json.loads(path.read_text()), 'frozen finite compiler receipt')
    print(json.dumps({'status': 'PASS', 'checks': COUNT, 'universe': out['finite_universes']}))


if __name__ == '__main__':
    main()
