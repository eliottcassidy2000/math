#!/usr/bin/env python3
"""Exact finite controls for abstract monodromy passports, not Keller realizations.
Run normally or with python3 -O. Every check uses require, never assert.
"""
import json
from collections import deque
from itertools import permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def identity(n):
    return tuple(range(n))


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    out = [0] * len(p)
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)


def group(generators, n):
    gens = list(generators)
    seen = {identity(n)}
    queue = deque(seen)
    while queue:
        p = queue.popleft()
        for g in gens:
            q = compose(g, p)
            if q not in seen:
                seen.add(q)
                queue.append(q)
    return seen


def transposition(n, i, j):
    p = list(range(n))
    p[i], p[j] = p[j], p[i]
    return tuple(p)


def fixed(p):
    return {i for i, j in enumerate(p) if i == j}


def normal_closure(elements, G, n):
    conjugates = {compose(compose(g, h), inverse(g))
                  for g in G for h in elements}
    return group(conjugates, n)


def fixed_generated(G, n):
    # Fixed-point property is conjugacy-invariant, so generated subgroup is normal.
    return group([g for g in G if fixed(g)], n)


def lift_outer(h, m, n):
    return tuple(h[i] * n + j for i in range(m) for j in range(n))


def lift_inner(g, block, m, n):
    return tuple(i * n + (g[j] if i == block else j)
                 for i in range(m) for j in range(n))


def wreath(inner_generators, outer_generators, m, n):
    gens = [lift_outer(h, m, n) for h in outer_generators]
    gens += [lift_inner(g, i, m, n)
             for i in range(m) for g in inner_generators]
    return group(gens, m * n)


def induced_outer(p, m, n):
    out = tuple(p[i * n] // n for i in range(m))
    require(all(p[i * n + j] // n == out[i]
                for i in range(m) for j in range(n)), 'block preservation')
    return out


def fixed_formula(p, m, n):
    h = induced_outer(p, m, n)
    # Convention p(i,j)=(h(i), g_i(j)).
    return sum(sum(p[i * n + j] % n == j for j in range(n))
               for i in range(m) if h[i] == i)


def run():
    result = {'scope': 'FINITE-EXACT abstract permutation controls; no polynomial realization'}
    paths = []
    # Explicit groups are small; the companion proof covers every d>=3.
    for d in range(3, 8):
        generators = [transposition(d, i, i + 1) for i in range(d - 1)]
        G = group(generators, d)
        require(len(G) == len(list(permutations(range(d)))), 'adjacent swaps generate S_d')
        require(all(len(fixed(g)) == d - 2 for g in generators), 'fixed-sheet count')
        require(normal_closure([generators[0]], G, d) == G, 'one transposition normally generates')
        I = {identity(d), generators[0]}
        require(any(compose(compose(g, generators[0]), inverse(g)) not in I for g in G),
                'cyclic inertia is not normal')
        paths.append({'degree': d, 'group_order': len(G),
                      'path_edges': len(generators), 'fixed_per_edge': d - 2})
    result['adjacent_transposition_paths'] = paths
    c2 = [transposition(2, 0, 1)]
    s3 = [transposition(3, 0, 1), transposition(3, 1, 2)]
    W1 = wreath(c2, s3, 3, 2)  # Inner C2, outer S3.
    W2 = wreath(s3, c2, 2, 3)  # Inner S3, outer C2.
    require(len(W1) == 48 and len(W2) == 72, 'wreath orders')
    rows = []
    for name, G, m, n in [('C2 wr S3', W1, 3, 2), ('S3 wr C2', W2, 2, 3)]:
        require(all(len(fixed(p)) == fixed_formula(p, m, n) for p in G), 'wreath fixed formula')
        N = fixed_generated(G, m * n)
        rows.append({'group': name, 'degree': m * n, 'order': len(G),
                     'fixed_generated_order': len(N),
                     'fixed_generated_outer_order': len({induced_outer(p, m, n) for p in N})})
    require(fixed_generated(W1, 6) == W1, 'inner cyclic blocks can be hidden')
    N2 = fixed_generated(W2, 6)
    require(len(N2) == 36, 'regular outer quotient survives')
    require({induced_outer(p, 2, 3) for p in N2} == {identity(2)}, 'outer quotient obstruction')
    result['wreath_controls'] = rows
    # Degree-four dihedral group: a nonregular nilpotent action.
    rotation = (1, 2, 3, 0)
    reflection = (0, 3, 2, 1)
    D8 = group([rotation, reflection], 4)
    ND = fixed_generated(D8, 4)
    require(len(D8) == 8 and len(ND) == 4, 'nonregular nilpotent exclusion')
    result['nilpotent_control'] = {'group': 'D8 on square vertices', 'order': 8,
                                 'fixed_generated_order': len(ND)}
    # Marked survivors need not be all inertia-fixed sheets.
    p = lift_inner(c2[0], 0, 3, 2)
    outer_survivors = {1}  # Block 2 is fixed but omitted by the affine-chart mark.
    inner_survivors = {1: {0}}  # Other unramified sheet likewise need not be included.
    marked = {i * 2 + j for i in outer_survivors for j in inner_survivors[i]}
    require(marked == {2} and marked < fixed(p), 'marked is strictly smaller than fixed')
    result['mark_control'] = {'inertia_fixed': sorted(fixed(p)), 'marked_survivors': sorted(marked),
                              'status': 'formal marks only, not claimed geometrically realizable'}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    run()
