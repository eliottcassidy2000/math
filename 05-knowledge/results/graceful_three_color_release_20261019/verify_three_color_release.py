#!/usr/bin/env python3
"""Exact arithmetic audits for THREE_COLOR_RELEASE_DICHOTOMY.md; stdlib only.
Checks are explicit exceptions, so python -O has identical behavior.
"""
import itertools
import json
from collections import Counter, deque

CHECKS = 0

def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)

def partition(n, edges):
    parent = list(range(n))
    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edges:
        x, y = root(a), root(b)
        parent[x] = y
    roots = [root(x) for x in range(n)]
    names = {x: i for i, x in enumerate(dict.fromkeys(roots))}
    return tuple(names[x] for x in roots)

def rank(n, edges):
    return n - len(set(partition(n, edges)))

def color_edges(m, d):
    return tuple((i, i + d) for i in range(m + 1 - d))

def subsets(items):
    return (tuple(c) for k in range(len(items) + 1)
            for c in itertools.combinations(items, k))

def graph_components(states):
    adjacency = [[] for _ in states]
    for i, x in enumerate(states):
        for j in range(i):
            if sum(a != b for a, b in zip(x, states[j])) <= 2:
                adjacency[i].append(j)
                adjacency[j].append(i)
    unseen = set(range(len(states)))
    components = []
    while unseen:
        todo = [unseen.pop()]
        component = set(todo)
        while todo:
            for w in adjacency[todo.pop()]:
                if w in unseen:
                    unseen.remove(w)
                    component.add(w)
                    todo.append(w)
        components.append(component)
    return components

def audit(m, fixed_by_color, old_colors, expected):
    n = m + 1
    check(set(fixed_by_color).isdisjoint(old_colors), 'color disjointness')
    check(set(fixed_by_color) | set(old_colors) == set(range(1, m + 1)),
          'all difference colors present')
    check(len(old_colors) == 3, 'exactly three old colors')
    for d, (a, b) in fixed_by_color.items():
        check(0 <= a < b <= m and b - a == d, 'valid fixed edge')
    check(rank(n, fixed_by_color.values()) == len(fixed_by_color), 'fixed forest')
    old_union = sum((color_edges(m, d) for d in old_colors), ())
    rows = []
    original = None
    original_components = None
    for release in [None] + sorted(fixed_by_color):
        fixed = tuple(e for d, e in fixed_by_color.items() if d != release)
        colors = tuple(sorted(old_colors + (() if release is None else (release,))))
        states = []
        for state in itertools.product(*(color_edges(m, d) for d in colors)):
            check(all(b - a == d for d, (a, b) in zip(colors, state)),
                  'enumerated exact differences')
            if rank(n, fixed + state) == m:
                states.append(state)
        check(bool(states), 'nonempty enlarged fiber')
        components = graph_components(states)
        sizes = sorted(map(len, components))
        check(sizes == expected[release], 'expected exact component census')
        ranks = {d: rank(n, fixed + color_edges(m, d)) - len(fixed)
                 for d in colors}
        staircase = all(q >= len(colors) - i
                        for i, q in enumerate(sorted(ranks.values(), reverse=True)))
        row = {'release': release, 'colors': colors, 'states': len(states),
               'component_sizes': sizes, 'color_ranks': ranks,
               'rank_staircase': staircase}
        if release is None:
            original = states
            original_components = components
        else:
            old_ranks = {}
            tight = []
            for cs in subsets(old_colors):
                q = rank(n, fixed + sum((color_edges(m, d) for d in cs), ())) - len(fixed)
                old_ranks[','.join(map(str, cs)) or 'empty'] = q
                check(q >= len(cs), 'Rado inequalities inherited from original fiber')
                if q == len(cs):
                    tight.append(cs)
            kernel = sorted(set().union(*map(set, tight)))
            check(tuple(kernel) in tight, 'maximal tight union is tight')
            for cs in subsets(tuple(d for d in old_colors if d not in kernel)):
                if cs:
                    combined = tuple(sorted(tuple(kernel) + cs))
                    check(old_ranks[','.join(map(str, combined))] - len(kernel) > len(cs),
                          'strict Rado inequality after kernel contraction')
            old_rank = rank(n, fixed + old_union) - len(fixed)
            bridge = rank(n, fixed + old_union) == m - 1
            check(old_rank in (3, 4), 'one-edge total-rank dichotomy')
            check(bridge == (old_rank == 3), 'bridge/rank equivalence')
            buffer_full = ranks[release] == 4
            row.update(old_union_rank=old_rank, bridge=bridge,
                       buffer_full_rank=buffer_full, tight_old_subsets=tight,
                       tight_kernel=kernel, old_subset_ranks=old_ranks)
            if buffer_full and not bridge:
                check(len(kernel) <= 2, 'small kernel in nonbridge case')
                check(len(components) == 1, 'positive dichotomy conclusion')
            if bridge:
                sides = partition(n, fixed + old_union)
                crossing = [e for e in color_edges(m, release) if sides[e[0]] != sides[e[1]]]
                old_positions = [colors.index(d) for d in old_colors]
                projections = Counter(tuple(state[i] for i in old_positions) for state in states)
                check(set(projections) == set(original), 'bridge product old-state projection')
                check(set(projections.values()) == {len(crossing)}, 'bridge product multiplicity')
                original_sets = {frozenset(original[i] for i in comp)
                                 for comp in original_components}
                for comp in components:
                    projected = frozenset(tuple(states[i][j] for j in old_positions) for i in comp)
                    check(projected in original_sets, 'bridge component projects exactly')
                check(len(components) == len(original_components), 'bridge component bijection')
                row['crossing_buffer_edges'] = crossing
        rows.append(row)
    return {'m': m, 'fixed': fixed_by_color, 'old_colors': old_colors,
            'old_free_union_connected': rank(n, old_union) == m,
            'release_census': rows}

def main():
    witness12 = audit(12,
        {1:(5,6),2:(1,3),3:(8,11),4:(3,7),5:(0,5),7:(4,11),
         8:(1,9),10:(0,10),12:(0,12)}, (6,9,11),
        {None:[2,3],1:[54],2:[56],3:[24],4:[20],5:[54],7:[14],8:[20],10:[16],12:[2,3]})
    witness14 = audit(14,
        {1:(6,7),2:(1,3),3:(10,13),4:(8,12),5:(0,5),6:(3,9),7:(6,13),
         9:(4,13),10:(1,11),12:(0,12),14:(0,14)}, (8,11,13),
        {None:[3,3],1:[6,6],2:[68],3:[26],4:[26],5:[24],6:[15],7:[26],
         9:[15],10:[22],12:[24],14:[3,3]})
    check(witness12['old_free_union_connected'], 'old witness global free connectivity')
    r3 = next(x for x in witness12['release_census'] if x['release'] == 3)
    check(r3['buffer_full_rank'] and not r3['bridge'] and not r3['rank_staircase'],
          'new theorem strictly extends staircase certificate')
    r1 = next(x for x in witness14['release_census'] if x['release'] == 1)
    check(r1['buffer_full_rank'] and r1['bridge'], 'nonforced full-buffer negative example')
    check(r1['crossing_buffer_edges'] == [(6,7),(7,8)], 'two legal released-buffer realizations')
    print(json.dumps({'status':'PASS', 'checks':CHECKS,
                      'audits':[witness12,witness14]}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
