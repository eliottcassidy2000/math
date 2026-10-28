#!/usr/bin/env python3
"""Exact GF(2) scope-loss, LCA rank, and affine join-tree controls.
No external libraries; explicit exceptions survive python -O.
"""
from itertools import product
import json
import random

checks = 0

def need(p, note):
    global checks
    checks += 1
    if not p:
        raise RuntimeError(note)

def basis(rows):
    piv = {}
    for row in rows:
        x = row
        while x:
            p = x.bit_length() - 1
            if p in piv:
                x ^= piv[p]
            else:
                piv[p] = x
                break
    return tuple(piv[p] for p in sorted(piv, reverse=True))

def span(rows):
    out = {0}
    for x in basis(rows):
        out |= {y ^ x for y in tuple(out)}
    return frozenset(out)

def rank(rows):
    return len(basis(rows))

def dot(x, y):
    return (x & y).bit_count() & 1

def sol(rows, n, offset=0):
    return frozenset(x for x in range(1 << n)
                     if all(dot(r, x ^ offset) == 0 for r in rows))

def proj(C, scope):
    return frozenset(x & scope for x in C)

def dim(C):
    size = len(C)
    need(size > 0 and size & (size - 1) == 0, 'not an affine-size power of 2')
    return size.bit_length() - 1

def affine_projection(equations, scope, n):
    """Eliminate outside coordinates; return independent (mask,rhs) checks.

    None denotes empty support. Coefficients use the original n labels.
    """
    work = [r | (b << n) for r, b in equations]
    for p in range(n):
        if (scope >> p) & 1:
            continue
        pivot_index = next((i for i, row in enumerate(work) if (row >> p) & 1), None)
        if pivot_index is None:
            continue
        pivot = work.pop(pivot_index)
        work = [row ^ pivot if (row >> p) & 1 else row for row in work]
    full = (1 << n) - 1
    if rank(work) != rank(row & full for row in work):
        return None
    return [(row & full, (row >> n) & 1) for row in basis(work)]


def affine_message_count(bags, local_rows, n):
    """Exact chain join-tree count by elimination, no assignment enumeration."""
    child_support = []
    exponent = 0
    for v in range(len(bags) - 1, -1, -1):
        B = bags[v]
        combined = local_rows[v] + child_support
        G = affine_projection(combined, B, n)
        if G is None:
            return 0
        separator = B & bags[v - 1] if v else 0
        support = affine_projection(G, separator, n)
        if support is None:
            raise RuntimeError('nonempty affine set has empty projection')
        dim_G = B.bit_count() - rank(r for r, b in G)
        dim_S = separator.bit_count() - rank(r for r, b in support)
        exponent += dim_G - dim_S
        child_support = support
    return 1 << exponent


def subspaces(n):
    spaces = {frozenset({0})}
    for x in range(1, 1 << n):
        spaces |= {span((*W, x)) for W in tuple(spaces)}
    return sorted(spaces, key=lambda W: (len(W), tuple(sorted(W))))

def scope_test(W, scopes, n, offset):
    C = sol(W, n, offset)
    local = span(r for B in scopes for r in W if r & ~B == 0)
    predicted = sol(local, n, offset)
    projections = [proj(C, B) for B in scopes]
    joined = frozenset(x for x in range(1 << n)
                       if all(x & B in P for B, P in zip(scopes, projections)))
    d = dim(W) - dim(local)
    need(joined == predicted, 'scope annihilator theorem')
    need(len(joined) == len(C) * (1 << d), 'scope defect count')
    return C, joined, d

pair_cases = 0
for n in range(1, 5):
    full = (1 << n) - 1
    for W in subspaces(n):
        for A, B in product(range(1 << n), repeat=2):
            if A | B != full:
                continue
            for offset in (0, full):
                C, joined, d = scope_test(W, (A, B), n, offset)
                need(d == dim(proj(C, A)) + dim(proj(C, B))
                     - dim(proj(C, A & B)) - dim(C), 'two-bag rank defect')
                pair_cases += 1

rng = random.Random(20261028)
random_scope_cases = 0
for _ in range(1200):
    n = rng.randrange(1, 8)
    W = span(rng.randrange(1 << n) for _ in range(rng.randrange(0, n + 3)))
    scopes = [rng.randrange(1 << n) for _ in range(rng.randrange(1, 8))]
    scope_test(W, scopes, n, rng.randrange(1 << n))
    random_scope_cases += 1

# Complete binary leaf hierarchy: row is owned at its support's LCA.
node_scopes = {i: 1 << i for i in range(8)}
children = {}
next_id = 8
level = list(range(8))
while len(level) > 1:
    new = []
    for j in range(0, len(level), 2):
        v = next_id
        next_id += 1
        children[v] = level[j:j + 2]
        node_scopes[v] = node_scopes[level[j]] | node_scopes[level[j + 1]]
        new.append(v)
    level = new
root = level[0]
lca_cases = 0
for _ in range(1200):
    rows = [rng.randrange(1, 256) for _ in range(rng.randrange(0, 24))]
    owned = {v: [] for v in node_scopes}
    for row in rows:
        v = min((v for v, B in node_scopes.items() if row & ~B == 0),
                key=lambda v: node_scopes[v].bit_count())
        owned[v].append(row)
    accumulated = {}
    increments = {}
    for v in node_scopes:
        previous = [r for c in children.get(v, []) for r in accumulated[c]]
        accumulated[v] = basis(previous + owned[v])
        increments[v] = rank(accumulated[v]) - rank(previous)
    need(sum(increments.values()) == rank(rows), 'LCA telescoping rank')
    need(8 - sum(increments.values()) == dim(sol(rows, 8)), 'LCA entropy budget')
    lca_cases += 1

# Join-tree affine support + constant fibre exponent recurrence.
# We implement supports by explicit sets as an independent oracle for elimination.
bags = [0b000111, 0b001110, 0b011100, 0b111000]
tree_cases = 0
for _ in range(1500):
    local_rows = []
    for B in bags:
        admissible = [r for r in range(64) if r & ~B == 0]
        local_rows.append([(rng.choice(admissible), rng.randrange(2))
                           for _ in range(rng.randrange(0, 5))])
    all_rows = [eq for eqs in local_rows for eq in eqs]
    global_solutions = frozenset(x for x in range(64)
                               if all(dot(r, x) == b for r, b in all_rows))
    child_support = None
    child_separator = 0
    exponent = 0
    for v in range(3, -1, -1):
        B = bags[v]
        G = frozenset(x for x in range(64) if x & ~B == 0
                      and all(dot(r, x) == b for r, b in local_rows[v])
                      and (child_support is None or x & child_separator in child_support))
        separator = B & bags[v - 1] if v else 0
        support = proj(G, separator)
        if G:
            exponent += dim(G) - dim(support)
            fibres = [sum(x & separator == s for x in G) for s in support]
            need(len(set(fibres)) == 1, 'affine fibres not constant')
        child_support, child_separator = support, separator
    predicted = (1 << exponent) if child_support else 0
    need(predicted == len(global_solutions), 'join-tree solution count')
    need(affine_message_count(bags, local_rows, 6) == predicted, 'elimination-message count')
    # An independent augmented-rank consistency test.
    R = rank(r for r, b in all_rows)
    augmented = rank(r | (b << 6) for r, b in all_rows)
    need(predicted == (1 << (6 - R) if augmented == R else 0), 'augmented rank')
    tree_cases += 1

H = [sum(1 << (j - 1) for j in range(1, 8) if (j >> k) & 1)
     for k in range(3)]
W = span(H)
C = sol(W, 7)
lines = [sum(1 << (j - 1) for j in (a, b, a ^ b))
         for a in range(1, 8) for b in range(a + 1, 8) if a ^ b > b]
fours = [127 ^ L for L in lines]
need(len(set(lines)) == 7 and set(fours) == W - {0}, 'Fano support list')
_, joined_lines, line_defect = scope_test(W, lines, 7, 0)
_, joined_fours, four_defect = scope_test(W, fours, 7, 0)
need(line_defect == 3 and len(joined_lines) == 128, 'Fano line loss')
need(four_defect == 0 and joined_fours == C, 'Fano four recovery')
need(rank(fours) == 3 and len(C) == 16, 'seven correlated guards')
circuits = [x for x in C if x and not any(y and y != x and y & x == y for y in C)]
need(set(circuits) == set(lines + fours), 'Fano matroid circuits')
need(len(C - {0}) == 15 and 127 in C and 127 not in circuits, '15 nonzero / 14 minimal')
# Any three distinct nonzero dual checks: 7 dependent triples, 28 bases.
from itertools import combinations
triple_ranks = {r: sum(rank(t) == r for t in combinations(fours, 3)) for r in (2, 3)}
need(triple_ranks == {2: 7, 3: 28}, 'correlated-check triples')
# A hidden 3-variable parity survives a genuine join-tree skeleton.
_, parity_join, parity_defect = scope_test(span([7]), [3, 6], 3, 0)
need(parity_defect == 1 and len(parity_join) == 8, 'hidden global parity')
# Pairwise rowspace intersections miss a three-way dependency.
U, V, Z = span([3]), span([6]), span([5])
need(U & V == U & Z == V & Z == frozenset({0}) and rank([3, 6, 5]) == 2,
     'three-way matroid dependency')
# Duplicate guard cannot be charged twice; contradictory duplicate is impossible.
need(rank([3, 3]) == 1 and len(sol([3, 3], 2)) == 2, 'duplicate guard')
need(not any(dot(3, x) == 0 and dot(3, x) == 1 for x in range(4)), 'inconsistent duplicate')

print(json.dumps({
    'status': 'PASS', 'explicit_checks': checks,
    'two_bag_affine_cases': pair_cases,
    'random_arbitrary_scope_cases': random_scope_cases,
    'LCA_rank_cases': lca_cases,
    'join_tree_affine_cases': tree_cases,
    'Fano': {'H_rank': rank(H), 'code_dimension': dim(C), 'circuits': len(circuits),
             'nonzero_codewords': len(C) - 1, 'line_projection_lost_bits': line_defect,
             'complementary_four_projection_lost_bits': four_defect,
             'seven_guard_joint_probability': '1/8', 'triple_guard_ranks': triple_ranks},
    'hidden_global_parity_lost_bits': parity_defect,
}, sort_keys=True, indent=2))
