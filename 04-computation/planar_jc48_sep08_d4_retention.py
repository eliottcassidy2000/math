#!/usr/bin/env python3
"""Complete partial-retention lower bounds for the literal W(D4) model.

No central-involution or full-fixed-set filter is imposed.  Geometric
access, cusp injection and Euler integration are analytic dependencies.
"""
from collections import Counter, deque
from itertools import combinations, permutations, product
import hashlib
import json

GATES = 0


def check(condition, name):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(name)


letters = (1, 2, 3, 4, -1, -2, -3, -4)
index = {x: i for i, x in enumerate(letters)}
elements = []
for permutation in permutations((1, 2, 3, 4)):
    for signs in product((-1, 1), repeat=4):
        if signs[0] * signs[1] * signs[2] * signs[3] != 1:
            continue
        positive = tuple(permutation[i] * signs[i] for i in range(4))
        elements.append(tuple(index[x] for x in positive + tuple(-x for x in positive)))
elements.sort()
lookup = {g: i for i, g in enumerate(elements)}
check(len(elements) == len(lookup) == 192, 'exact even signed permutation universe')
multiply = [[lookup[tuple(p[q[i]] for i in range(8))] for q in elements] for p in elements]
identity = lookup[tuple(range(8))]
inverse = [next(j for j in range(192) if multiply[i][j] == identity) for i in range(192)]


def reflection(i, j, negative=False):
    images = list(letters)
    for sign in (-1, 1):
        images[index[sign * i]] = sign * j * (-1 if negative else 1)
        images[index[sign * j]] = sign * i * (-1 if negative else 1)
    return lookup[tuple(index[x] for x in images)]


a, b, c, d = reflection(1, 2), reflection(2, 3), reflection(3, 4), reflection(3, 4, True)
e = multiply[multiply[b][c]][inverse[b]]
f = multiply[multiply[inverse[e]][a]][e]
pairs = ((b, c), (b, d), (a, e), (f, b), (a, d), (c, d))


def closure(generators):
    found, todo = {identity}, [identity]
    while todo:
        x = todo.pop()
        for g in generators:
            y = multiply[x][g]
            if y not in found:
                found.add(y)
                todo.append(y)
    return frozenset(found)


check(len(closure((a, b, c, d))) == 192, 'marked reflections generate the whole model')
for g in (a, b, c, d):
    check(multiply[g][g] == identity, 'marked generators are involutions')
for g in (a, c, d):
    check(multiply[multiply[g][b]][g] == multiply[multiply[b][g]][b], 'marked braid relation')
for g, h in ((a, c), (a, d), (c, d)):
    check(multiply[g][h] == multiply[h][g], 'marked commuting leaves')

# Every transitive action having a reflection-fixed point is W/H with b in H.
# Generate every such subgroup, with no central or index cutoff.
base = closure((b,))
subgroups, queue = {base: (b,)}, deque([base])
while queue:
    subgroup = queue.popleft()
    generators = subgroups[subgroup]
    covered = set(subgroup)
    for g in range(192):
        if g in covered:
            continue
        covered.update(multiply[g][h] for h in subgroup)
        extension = closure(generators + (g,))
        if extension not in subgroups:
            subgroups[extension] = generators + (g,)
            queue.append(extension)
check(len(subgroups) == 53, 'all fifty-three reflection-containing subgroups')
for subgroup, generators in subgroups.items():
    check(closure(generators) == subgroup and b in subgroup, 'subgroup witness and fixed point')
    for g in range(192):
        check(closure(generators + (g,)) in subgroups, 'every element extension checked separately')


def bounds(degree, fixed, counts, retained):
    cusp = [max(0, retained - fixed + n, (3 * retained - degree + 1) // 2)
            for n in counts[:3]]
    node = [max(0, degree - 2 * retained, degree - 2 * fixed + n)
            for n in counts[3:]]
    return cusp, node, -2 * retained + sum(cusp) + sum(node)


rows, candidates, trials = [], [], 0
for subgroup, generators in subgroups.items():
    assignment, representatives = {}, []
    for g in range(192):
        if g in assignment:
            continue
        label = len(representatives)
        representatives.append(g)
        coset = {multiply[g][h] for h in subgroup}
        check(len(coset) == len(subgroup), 'full coset cardinality')
        for x in coset:
            assignment[x] = label
    degree = len(representatives)
    check(degree * len(subgroup) == len(assignment) == 192, 'complete coset partition')
    action = {q: tuple(assignment[multiply[q][g]] for g in representatives)
              for pair in pairs for q in pair}
    fixed = sum(action[b][i] == i for i in range(degree))
    check(fixed > 0, 'positive fixed-point entry')
    for q in action:
        check(sum(action[q][i] == i for i in range(degree)) == fixed,
              'all original access meridians have the same fixed count')
    counts = [sum(action[q][i] == i and action[r][i] == i for i in range(degree))
              for q, r in pairs]
    row = {'degree': degree, 'fixed': fixed, 'counts': counts,
           'generators': generators, 'retention_bounds': []}
    for retained in range(1, min(fixed, degree - 1) + 1):
        trials += 1
        cusp, node, lower = bounds(degree, fixed, counts, retained)
        check(all(n >= 0 for n in cusp + node), 'nonnegative local count lower bounds')
        row['retention_bounds'].append([retained, cusp, node, lower])
        if lower <= 1:
            candidates.append((degree, fixed, retained, tuple(counts), lower))
        else:
            check(lower >= 2, 'strict integer Euler obstruction')
    rows.append(row)
check(len(candidates) == 6, 'only six labelled nontrivial potential actions')
for degree, fixed, retained, counts, lower in candidates:
    check((degree, fixed, retained, lower) == (4, 2, 2, 1), 'every survivor is degree four and full retained')
    check(counts[:3] == (1, 1, 1) and sorted(counts[3:]) == [0, 0, 2], 'sharp survivor local profile')

# Exhaustive retained-subset controls in the natural eight-letter model.
# At a cusp B=(sigma*tau)A; a joint fixed point lies in A iff it lies in B.
fixed8 = {q: frozenset(i for i in range(8) if elements[q][i] == i) for pair in pairs for q in pair}
proper_control = False
for q, r in pairs[:3]:
    joint = len(fixed8[q] & fixed8[r])
    for size in range(1, len(fixed8[q]) + 1):
        for chosen in combinations(sorted(fixed8[q]), size):
            A = frozenset(chosen)
            B = frozenset(elements[q][elements[r][i]] for i in A)
            check(B <= fixed8[r], 'actual braid reaccess preserves retained fixedness')
            n = len(A & B)
            check(n >= max(0, size - len(fixed8[q]) + joint, (3 * size - 8 + 1) // 2),
                  'cusp lower bound for every natural retained subset')
            if size == 1 and n == 0:
                proper_control = True
check(proper_control, 'hostile: retained sheets can be a proper fixed subset')
for q, r in pairs[3:]:
    joint = len(fixed8[q] & fixed8[r])
    for size in range(1, len(fixed8[q]) + 1):
        for chosen_a in combinations(sorted(fixed8[q]), size):
            for chosen_b in combinations(sorted(fixed8[r]), size):
                A, B = set(chosen_a), set(chosen_b)
                omega = len((set(range(8)) - A) & (set(range(8)) - B))
                check(omega >= max(0, 8 - 2 * size, 8 - 2 * len(fixed8[q]) + joint),
                      'node bound for every natural pair of retained subsets')
natural_counts = [len(fixed8[q] & fixed8[r]) for q, r in pairs]
check(natural_counts == [2, 2, 2, 0, 0, 4], 'signed-eight fixed-count hostile')
check([bounds(8, 4, natural_counts, k)[2] for k in range(1, 5)] == [16, 8, 5, 2],
      'all natural-eight retention sizes fail Euler one')

print('Literal W(D4) partial-retention Euler obstruction PASS; actual geometric suppliers required')
print('Universe:192 signed permutations;53 reflection-containing subgroups; retention trials:', trials)
print('Only potential Euler-one actions: six labelled degree4 stabilizers, fixed=retained=2')
print('Extra nodes contribute nonnegative overlap; full-fixed and central-quotient hypotheses are unnecessary')
print('Natural degree8 lower bounds for retained sizes1..4:16,8,5,2; proper-retention hostile checked')
print('Always-active gates:', GATES)
print('Semantic SHA256:', hashlib.sha256(json.dumps(rows, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
