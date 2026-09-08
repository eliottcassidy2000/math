#!/usr/bin/env python3
"""Exact all-pair certificate for H4 permutations of type (3)(3).

Status: RESERVED pending independent audit.  A proof reduces every pair to
12 labels; no full-tuple ambient-degree census is used.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json

GATES = 0


def check(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def mul(p, q):
    return tuple(p[x] for x in q)


def inv(p):
    q = list(p)
    for i, j in enumerate(p):
        q[j] = i
    return tuple(q)


def power(p, n):
    q = tuple(range(len(p)))
    for _ in range(n):
        q = mul(q, p)
    return q


def cycle(p, block, orientation=0):
    q = list(p)
    block = tuple(block)
    if orientation:
        block = (block[0], block[2], block[1])
    for x, y in zip(block, block[1:] + block[:1]):
        q[x] = y
    return tuple(q)


def moved_cycles(p):
    seen = set()
    answer = []
    for first in range(len(p)):
        if first in seen or p[first] == first:
            continue
        row = []
        x = first
        while x not in row:
            row.append(x)
            seen.add(x)
            x = p[x]
        answer.append(tuple(row))
    return answer


def word(p, q, length):
    answer = tuple(range(len(p)))
    for j in range(length):
        answer = mul(answer, p if j % 2 == 0 else q)
    return answer


def point_word(p, q, length, x):
    # Direct rightmost-first letter action, independent of mul/word.
    for j in range(length - 1, -1, -1):
        x = (p if j % 2 == 0 else q)[x]
    return x


def flags(p, q):
    return (word(p, q, 3) == word(q, p, 3),
            word(p, q, 5) == word(q, p, 5),
            mul(p, q) == mul(q, p))


def joint_orbits(p, q):
    unseen = set(range(len(p)))
    answer = []
    while unseen:
        first = min(unseen)
        found = {first}
        queue = [first]
        for x in queue:
            for g in (p, q):
                y = g[x]
                if y not in found:
                    found.add(y)
                    queue.append(y)
        unseen -= found
        if len(found) > 1:
            answer.append(tuple(sorted(found)))
    return tuple(answer)


def matrix(p, q):
    rows = [set(c) for c in moved_cycles(p)]
    cols = [set(c) for c in moved_cycles(q)]
    entries = [[len(a & b) for b in cols] for a in rows]
    return min(tuple(entries[i][j] for i in rr for j in cc)
               for rr in ((0, 1), (1, 0)) for cc in ((0, 1), (1, 0)))


def class33(degree):
    identity = tuple(range(degree))
    for aa in combinations(range(degree), 3):
        rest = tuple(x for x in range(degree) if x not in aa)
        for bb in combinations(rest, 3):
            if aa > bb:
                continue
            for ao, bo in product(range(2), repeat=2):
                yield cycle(cycle(identity, aa, ao), bb, bo)


D = 12
identity = tuple(range(D))
sigma = cycle(cycle(identity, (0, 1, 2)), (3, 4, 5))
universe = tuple(class33(D))
allset = set(universe)
check(len(universe) == len(allset) == 36960, 'complete unfiltered conjugacy class')
check(36960 == 924 * 40, 'independent support-size universe count')

expected_tables = (
    {(0, 2, 2, 0): 270, (0, 2, 3, 0): 36, (0, 3, 3, 0): 1,
     (1, 1, 1, 1): 270, (1, 2, 2, 1): 9},
    {(0, 1, 1, 0): 3240, (0, 1, 3, 0): 180,
     (0, 3, 3, 0): 1, (1, 2, 2, 1): 18},
    {(0, 0, 0, 0): 40, (0, 0, 0, 3): 160, (0, 3, 3, 0): 4},
)
expected_components = (
    {(3, 3): 1, (3, 4): 36, (4, 4): 270, (6,): 9, (8,): 270},
    {(3, 3): 1, (3, 5): 180, (5, 5): 3240, (6,): 18},
    {(3, 3): 4, (3, 3, 3): 160, (3, 3, 3, 3): 40},
)
hist = [Counter() for _ in range(3)]
component_hist = [Counter() for _ in range(3)]
data = {}
raw = []
A = set(range(6))
a_cycles = [set(c) for c in moved_cycles(sigma)]
for tau in universe:
    cyc = moved_cycles(tau)
    check(sorted(map(len, cyc)) == [3, 3] and power(tau, 3) == identity,
          'literal type and order')
    ff = flags(sigma, tau)
    for length, flag in zip((3, 5), ff[:2]):
        direct = all(point_word(sigma, tau, length, x) ==
                     point_word(tau, sigma, length, x) for x in range(D))
        check(direct == flag, 'independent literal odd word action')
    cells = matrix(sigma, tau)
    orbits = joint_orbits(sigma, tau)
    sizes = tuple(sorted(map(len, orbits)))
    data[tau] = (ff, cells, sizes)
    raw.append((tau, ff, cells, sizes))
    for k, yes in enumerate(ff):
        if yes:
            hist[k][cells] += 1
            component_hist[k][sizes] += 1
    if ff[0]:
        check(tau == sigma or sizes in ((3, 4), (4, 4), (6,), (8,)),
              'complete ordinary joint-orbit types')
        check(all(len(set(c) & A) >= 2 for c in cyc),
              'ordinary run bound in each three-cycle')
        check(len(set(x for c in cyc for x in c) - A) <= 2,
              'ordinary support enlargement at most two')
    if ff[1]:
        check(all(set(c) & A for c in cyc),
              'fifth run bound in every three-cycle')
    if ff[2]:
        check(all(not(set(c) & A) or set(c) in a_cycles for c in cyc),
              'commuting equal-cycle blocks do not mix')

for k in range(3):
    check(dict(hist[k]) == expected_tables[k], 'complete matrix histogram')
    check(dict(component_hist[k]) == expected_components[k],
          'complete component-size histogram')
check([sum(h.values()) for h in hist] == [586, 3439, 204],
      'three complete partner totals')
check(next(t for t in universe if data[t][0][0] and
           data[t][2] == (3, 3)) == sigma, 'ordinary split diagonal is unique')

# Full centralizer C3 wreath S2 x S6; the equal blocks may be exchanged.
centralizer_gens = [
    cycle(identity, (0, 1, 2)),
    cycle(identity, (3, 4, 5)),
    (3, 4, 5, 0, 1, 2, 6, 7, 8, 9, 10, 11),
]
for x in range(6, 11):
    g = list(identity)
    g[x], g[x+1] = g[x+1], g[x]
    centralizer_gens.append(tuple(g))
for g in centralizer_gens:
    check(mul(g, sigma) == mul(sigma, g), 'declared centralizer generator')
inverses = [inv(g) for g in centralizer_gens]
unseen = set(universe)
orbit_bank = []
orbit_hist = [Counter() for _ in range(3)]
while unseen:
    rep = min(unseen)
    found = {rep}
    queue = [rep]
    for p in queue:
        for g, gi in zip(centralizer_gens, inverses):
            q = mul(mul(g, p), gi)
            check(q in allset, 'centralizer image retains full raw universe')
            if q not in found:
                found.add(q)
                queue.append(q)
    unseen -= found
    ff, cells, sizes = data[rep]
    for p in found:
        check(data[p] == data[rep], 'exact orbit predicate and unordered-cell invariance')
    orbit_bank.append((rep, len(found), ff, cells, sizes))
    for k, yes in enumerate(ff):
        if yes:
            orbit_hist[k][cells] += len(found)
check(len(orbit_bank) == 61, 'complete raw centralizer orbit count')
check(sum(row[1] for row in orbit_bank) == 36960, 'all raw partners retained')
check([sum(row[2][k] for row in orbit_bank) for k in range(3)] == [5, 4, 6],
      'relation-specific orbit counts')
check(orbit_hist == hist, 'centralizer-weighted histograms equal raw histograms')

# Minimal hostile: same support size and same unordered matrix, different
# cyclic order and different relations. Neither admits two invariant triples.
id6 = tuple(range(6))
s6 = (1, 2, 0, 4, 5, 3)
t3 = (1, 3, 4, 0, 5, 2)
t5 = (1, 3, 5, 0, 2, 4)


def generated_group(gens):
    found = {tuple(range(len(gens[0])))}
    queue = list(found)
    for p in queue:
        for g in gens:
            q = mul(p, g)
            if q not in found:
                found.add(q)
                queue.append(q)
    return found


hostiles = []
for name, tau, expected_flags, order in (
        ('ordinary', t3, (True, False, False), 12),
        ('fifth', t5, (False, True, False), 60)):
    check(flags(s6, tau) == expected_flags, 'named literal hostile relations')
    group = generated_group((s6, tau))
    check(len(group) == order and len({p[0] for p in group}) == 6,
          'named group order and actual transitivity')
    block_count = 0
    trials = 0
    for aa in combinations(range(6), 3):
        aa = frozenset(aa)
        bb = frozenset(range(6)) - aa
        if min(aa) > min(bb):
            continue
        trials += 1
        invariant = all(frozenset(g[x] for x in aa) in (aa, bb) for g in (s6, tau))
        block_count += invariant
    check(trials == 10 and block_count == 0, 'complete unordered two-triple block test')
    check(matrix(s6, tau) == (1, 2, 2, 1), 'equal cells do not retain cyclic order')
    hostiles.append((name, tau, order, block_count))

# Finite controls for the analytical centralizer argument; these do not replace
# its transitive-orbit proof or impose a full-tuple degree bound.
for length in (4, 8):
    check(length % 3 != 0, 'semiregular order-three obstruction')
check(6 // 3 == 2 and 2 < 3, 'six-orbit has only two centralizer cycles')
check(1 < 2, 'split 3+4 outside overlap cannot support an ordinary three-cycle')
check(2 < 4, 'ordinary b-c overlap forces a common a-cycle')
check(3 < 6 and 6-3 == 3, 'common cycle is a proper nontrivial invariant set')
for degree in (6, 7, 12, 25):
    ii = tuple(range(degree))
    g = cycle(cycle(ii, (0, 1, 2)), (3, 4, 5))
    check(flags(g, g) == (True, True, True), 'all-equal H4 control')
    check(len(joint_orbits(g, g)) == 2, 'equal type33 action is never transitive')

raw_hash = sha256(json.dumps(raw, separators=(',', ':')).encode()).hexdigest()
orbit_hash = sha256(json.dumps(orbit_bank, separators=(',', ':')).encode()).hexdigest()
report = {
    'universe': 'all 36960 conjugate partners of (123)(456) on 12 labels',
    'partner_totals': {'braid3': 586, 'braid5': 3439, 'commuting': 204},
    'matrix_tables': {
        name: [[list(k), v] for k, v in sorted(h.items())]
        for name, h in zip(('braid3', 'braid5', 'commuting'), hist)},
    'joint_orbit_sizes': {
        name: [[list(k), v] for k, v in sorted(h.items())]
        for name, h in zip(('braid3', 'braid5', 'commuting'), component_hist)},
    'raw_centralizer_orbits': 61,
    'admitted_orbit_counts': [5, 4, 6],
    'minimal_hostiles': hostiles,
    'raw_sha256': raw_hash,
    'orbit_sha256': orbit_hash,
    'always_active_gates': GATES,
    'scope': 'uniform H4 type(3)(3) equality via pair orbit types; no full-tuple census',
}
print(json.dumps(report, sort_keys=True, indent=2))

