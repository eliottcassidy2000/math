#!/usr/bin/env python3
"""Exact graphic replay; analytic proof covers all represented fields.
Standard library only. Checks survive python -O.
"""
import itertools
import json
from collections import Counter

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
    names = {x:i for i,x in enumerate(dict.fromkeys(roots))}
    return tuple(names[x] for x in roots)

def rank(n, edges):
    return n - len(set(partition(n, edges)))

def edges(d):
    return tuple((i, i+d) for i in range(17-d))

def subsets(items):
    return (c for k in range(len(items)+1) for c in itertools.combinations(items,k))

def components(states):
    neighbors = [[] for _ in states]
    for i, state in enumerate(states):
        for j in range(i):
            if sum(a != b for a,b in zip(state,states[j])) <= 2:
                neighbors[i].append(j)
                neighbors[j].append(i)
    unseen = set(range(len(states)))
    result = []
    while unseen:
        todo = [unseen.pop()]
        component = set(todo)
        while todo:
            for j in neighbors[todo.pop()]:
                if j in unseen:
                    unseen.remove(j)
                    component.add(j)
                    todo.append(j)
        result.append(component)
    return result

def matrix_check():
    # All zero-diagonal binary 3x3 matrices: only two satisfy both
    # invertibility and vanishing principal 2x2 minors.
    off = [(i,j) for i in range(3) for j in range(3) if i != j]
    survivors = []
    for bits in itertools.product((0,1), repeat=6):
        a = [[0]*3 for _ in range(3)]
        for (i,j), bit in zip(off,bits):
            a[i][j] = bit
        pair_zero = all(a[i][j]*a[j][i] == 0 for i in range(3) for j in range(i))
        determinant = (a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1]) % 2
        if pair_zero and determinant:
            check(all(sum(row) == 1 for row in a), 'one nonzero in every row')
            check(all(sum(a[i][j] for i in range(3)) == 1 for j in range(3)),
                  'one nonzero in every column')
            survivors.append(a)
        check(len(a) == 3 and all(a[i][i] == 0 for i in range(3)), 'matrix coverage')
    check(len(survivors) == 2, 'exactly two directed-three-cycle matrices')
    return {'matrices_checked':64, 'surviving_matrices':survivors}

def main():
    matrix_receipt = matrix_check()
    fixed = {1:(7,8),2:(7,9),3:(3,6),4:(7,11),5:(10,15),6:(1,7),
             7:(5,12),8:(2,10),9:(3,12),11:(3,14),12:(4,16),14:(0,14),16:(0,16)}
    old_colors = (10,13,15)
    old_union = sum((edges(d) for d in old_colors), ())
    check(set(fixed) | set(old_colors) == set(range(1,17)), 'full difference coverage')
    check(set(fixed).isdisjoint(old_colors), 'color disjointness')
    for d,(a,b) in fixed.items():
        check(0 <= a < b <= 16 and b-a == d, 'valid fixed edge')
    check(rank(17,fixed.values()) == 13, 'fixed forest')
    full_partition = partition(17,fixed.values())
    support = {}
    for d in old_colors:
        support[d] = Counter(tuple(sorted((full_partition[a],full_partition[b])))
                             for a,b in edges(d) if full_partition[a] != full_partition[b])
        check(len(support[d]) == 2, 'two directions per color')
    support_union = set.union(*(set(x) for x in support.values()))
    check(len(support_union) == 3, 'three common tree directions')
    check(not set.intersection(*(set(x) for x in support.values())), 'each direction omitted once')
    for I in subsets(old_colors):
        q = rank(17,tuple(fixed.values())+sum((edges(d) for d in I),()))-13
        check(q == min(3,2*len(I)), 'necessary three-color rank pattern')
    expected = {():[1,3],(1,):[2,6],(2,):[2,6],(3,):[21],(1,2):[8,24],(1,3):[86]}
    results = []
    original_states = None
    original_components = None
    for B, expected_sizes in expected.items():
        G = tuple(e for d,e in fixed.items() if d not in B)
        colors = tuple(sorted(old_colors+B))
        states = []
        for state in itertools.product(*(edges(d) for d in colors)):
            check(all(b-a == d for d,(a,b) in zip(colors,state)), 'candidate exact differences')
            if rank(17,G+state) == 16:
                states.append(state)
        comps = components(states)
        sizes = sorted(map(len,comps))
        check(sizes == expected_sizes, 'exact release component counts')
        ranks = {d:rank(17,G+edges(d))-len(G) for d in colors}
        row = {'released':B,'states':len(states),'component_sizes':sizes,'color_ranks':ranks,
               'old_union_rank':rank(17,G+old_union)-len(G)}
        if not B:
            original_states = states
            original_components = [frozenset(states[i] for i in c) for c in comps]
        if B and set(B) <= {1,2}:
            check(row['old_union_rank'] == 3, 'bridge-only old flat remains tight')
            old_positions = [colors.index(d) for d in old_colors]
            projections = Counter(tuple(s[i] for i in old_positions) for s in states)
            check(set(projections) == set(original_states), 'bridge product old choices')
            check(len(set(projections.values())) == 1, 'uniform bridge product multiplicity')
            multiplicity = next(iter(projections.values()))
            check(multiplicity == (8 if len(B) == 2 else 2), 'buffer factor count')
            for component in comps:
                projected = frozenset(tuple(states[i][j] for j in old_positions) for i in component)
                check(projected in original_components, 'old component cannot be crossed')
            row['buffer_factor_states'] = multiplicity
        if B == (1,2):
            check(ranks[1] == ranks[2] == 5, 'two full-rank buffers still fail')
        if B == (1,3):
            nonbuffer = (3,10,13,15)
            tight = []
            for I in subsets(nonbuffer):
                q = rank(17,G+sum((edges(d) for d in I),()))-len(G)
                check(q >= len(I), 'nonbuffer Rado inequality')
                if len(I) >= 3:
                    check(q > len(I), 'no tight subset of size at least three')
                if q == len(I):
                    tight.append(I)
            kernel = sorted(set().union(*map(set,tight)))
            check(len(kernel) <= 2, 'two-release small tight kernel')
            check(ranks[1] == 5, 'color-1 connected buffer')
            row['tight_nonbuffer_subsets'] = tight
            row['tight_kernel'] = kernel
        results.append(row)
    Q = tuple(fixed.values())+old_union
    for d in (1,2):
        check(rank(17,tuple(e for c,e in fixed.items() if c != d)+old_union) == 15,
              'fixed color-1 and color-2 edges are original bridges')
    check(rank(17,tuple(e for c,e in fixed.items() if c != 3)+old_union) == 16,
          'fixed color-3 edge is an original nonbridge')
    failed_cut = partition(17,tuple(e for d,e in fixed.items() if d != 1)+old_union)
    crossing3 = [e for e in edges(3) if failed_cut[e[0]] != failed_cut[e[1]]]
    check((5,8) in crossing3, 'color 3 crosses failed color-1 cut')
    print(json.dumps({'status':'PASS','checks':CHECKS,'matrix_audit':matrix_receipt,
                      'fixed_edges':fixed,'old_colors':old_colors,
                      'contact_supports':{d:{str(k):v for k,v in x.items()} for d,x in support.items()},
                      'release_audits':results,'crossing_color3_edges':crossing3},
                     indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
