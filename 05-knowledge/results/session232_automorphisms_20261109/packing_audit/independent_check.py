#!/usr/bin/env python3
"""Independent audit: permutation cycle-covers, explicit words and Burnside.
Does not import or execute the producer. Guards remain live under python -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
import hashlib
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def tournament(n, bits):
    arcs = set()
    for index, (u, v) in enumerate(combinations(range(n), 2)):
        arcs.add((u, v) if bits & (1 << index) else (v, u))
    return n, frozenset(arcs)

def circulant(n, steps):
    return n, frozenset((u, (u+s) % n) for u in range(n) for s in steps)

def order_join(A, B):
    n, a = A
    m, b = B
    return n+m, a | {(n+u, n+v) for u, v in b} | {(u, n+v) for u in range(n) for v in range(m)}

def words(T):
    n, arcs = T
    # Actual full words, not a subset counting recurrence.
    def extend(prefix, remaining):
        if not remaining:
            yield prefix
        else:
            for v in remaining:
                if not prefix or (prefix[-1], v) in arcs:
                    yield from extend(prefix+(v,), remaining-{v})
    return set(extend((), frozenset(range(n))))

def automorphisms(T):
    n, arcs = T
    # Generate permutations only within invariant score classes; check every arc.
    classes = defaultdict(list)
    for u in range(n):
        classes[sum(x == u for x, y in arcs)].append(u)
    blocks = list(classes.values())
    group = []
    for targets in product(*(permutations(block) for block in blocks)):
        p = [None] * n
        for block, images in zip(blocks, targets):
            for u, v in zip(block, images):
                p[u] = v
        p = tuple(p)
        if {(p[u], p[v]) for u, v in arcs} == arcs:
            group.append(p)
    return tuple(sorted(group))

def permutation_components(p):
    unseen = set(range(len(p)))
    cycles = []
    while unseen:
        u = min(unseen)
        cycle = []
        v = u
        while v in unseen:
            cycle.append(v)
            unseen.remove(v)
            v = p[v]
        check(v == u, 'invalid permutation decomposition')
        cycles.append(tuple(cycle))
    return cycles

def colored_covers(T):
    n, arcs = T
    p = [None] * n
    # Perfect matchings of allowed (row,image) pairs: loops or directed arcs.
    # Descending rows are merely a search order, not a component assumption.
    def match(row, unused):
        if row < 0:
            components = permutation_components(p)
            if any(len(c) % 2 == 0 for c in components):
                return
            nontrivial = [c for c in components if len(c) > 1]
            for bits in product((0, 1), repeat=len(nontrivial)):
                colors = [-1] * n
                for cycle, bit in zip(nontrivial, bits):
                    for vertex in cycle:
                        colors[vertex] = bit
                yield tuple(p), tuple(colors)
            return
        for image in unused:
            if row == image or (row, image) in arcs:
                p[row] = image
                yield from match(row-1, unused-{image})
    objects = list(match(n-1, frozenset(range(n))))
    check(len(objects) == len(set(objects)), 'cycle-cover duplication')
    return set(objects)

def relabel_cover(x, g):
    p, colors = x
    q = [None] * len(p)
    c = [None] * len(p)
    for v in range(len(p)):
        q[g[v]] = g[p[v]]
        c[g[v]] = colors[v]
    return tuple(q), tuple(c)

def cover_orbits(objects, group):
    by_rep = {}
    lookup = {}
    for x in sorted(objects):
        images = {relabel_cover(x, g) for g in group}
        check(images <= objects, 'cover action does not close')
        rep = min(images)
        lookup[x] = rep
        by_rep[rep] = images
    stabilizers = {rep: sum(relabel_cover(rep, g) == rep for g in group) for rep in by_rep}
    for rep, members in by_rep.items():
        check(len(members) * stabilizers[rep] == len(group), 'orbit-stabilizer failure')
    fixed_counts = [sum(relabel_cover(x, g) == x for x in objects) for g in group]
    check(sum(fixed_counts) == len(group) * len(by_rep), 'Burnside discrepancy')
    return by_rep, lookup, stabilizers, fixed_counts

def stats(T, group=None):
    n, arcs = T
    group = automorphisms(T) if group is None else group
    paths = words(T)
    covers = colored_covers(T)
    check(len(paths) == len(covers), 'independent OCF count mismatch')
    path_orbits = {min(tuple(g[v] for v in p) for g in group) for p in paths}
    for p in paths:
        images = {tuple(g[v] for v in p) for g in group}
        check(images <= paths and len(images) == len(group), 'path action failure')
    orbits, lookup, stabilizers, fixed = cover_orbits(covers, group)
    mass = sum((Fraction(1, s) for s in stabilizers.values()), Fraction(0))
    check(mass == len(path_orbits), 'mass does not equal explicit path quotient')
    check(mass == Fraction(len(covers), len(group)), 'mass identity failure')
    defect = sum((1-Fraction(1,s) for s in stabilizers.values()), Fraction(0))
    check(defect == len(orbits)-len(path_orbits), 'defect identity failure')
    check((defect == 0) == (len(group) == 1), 'strictness boundary failure')
    summary = {'n':n, 'H':len(paths), 'group_order':len(group),
               'path_orbits':len(path_orbits), 'coarse_packing_orbits':len(orbits),
               'mass':str(mass), 'defect':str(defect),
               'stabilizer_histogram':dict(sorted(Counter(stabilizers.values()).items())),
               'burnside_fixed_counts':fixed}
    return summary, (paths, covers, group, orbits, lookup, stabilizers, path_orbits)

def display_cover(x):
    p, colors = x
    return [{'cycle':c, 'color':colors[c[0]]} for c in permutation_components(p) if len(c)>1]

def witness(T):
    summary, data = stats(T)
    paths, covers, group, orbits, lookup, stabilizers, path_reps = data
    out = {'summary':summary, 'arcs':sorted(T[1]), 'automorphisms':group,
           'all_hamilton_words':sorted(paths), 'path_orbit_representatives':sorted(path_reps),
           'all_colored_packings':[display_cover(x) for x in sorted(covers)],
           'packing_orbits':[{'representative':display_cover(rep), 'size':len(members),
                              'stabilizer_order':stabilizers[rep]} for rep,members in sorted(orbits.items())]}
    return out, data

def restriction_check(lower, upper):
    low_paths, low_covers, low_group, low_orbits, low_lookup, low_stab, low_path_reps = lower
    high_paths, high_covers, high_group, high_orbits, high_lookup, high_stab, high_path_reps = upper
    n = len(next(iter(low_paths)))
    def restrict_cover(x):
        p,c = x
        check(all(p[v] < n for v in range(n)), 'cover crosses the restriction boundary')
        return p[:n], c[:n]
    pack_fibres = {rep: [] for rep in low_orbits}
    for rep,members in sorted(high_orbits.items()):
        images = {low_lookup[restrict_cover(x)] for x in members}
        check(len(images) == 1, 'packing quotient restriction ill-defined')
        pack_fibres[next(iter(images))].append(rep)
    check({len(v) for v in pack_fibres.values()} == {3}, 'actual packing fibre not three')
    path_fibres = {rep: [] for rep in low_path_reps}
    for rep in high_path_reps:
        # Check an entire upper orbit, not just its representative.
        images = set()
        for g in high_group:
            p = tuple(g[v] for v in rep)
            r = tuple(v for v in p if v < n)
            check(r in low_paths, 'path restriction not Hamiltonian')
            images.add(min(tuple(h[v] for v in r) for h in low_group))
        check(len(images) == 1, 'path quotient restriction ill-defined')
        path_fibres[next(iter(images))].append(rep)
    check({len(v) for v in path_fibres.values()} == {1}, 'actual path fibre not singleton')
    return {'packing_fibre_histogram':dict(sorted(Counter(map(len,pack_fibres.values())).items())),
            'path_fibre_histogram':dict(sorted(Counter(map(len,path_fibres.values())).items())),
            'packing_map': [{'target':display_cover(rep), 'preimages':[display_cover(x) for x in src]}
                             for rep,src in sorted(pack_fibres.items())]}

def main():
    all_records = []
    max_group = 0
    identities = 0
    for n in range(1,6):
        for bits in range(1 << (n*(n-1)//2)):
            T = tournament(n,bits)
            summary, data = stats(T)
            all_records.append((n,bits,summary))
            max_group = max(max_group,summary['group_order'])
            # A separate trivial-action positive control at every tournament.
            if summary['group_order'] > 1:
                stats(T, (tuple(range(n)),))
                identities += 1
    C = circulant(3, {1})
    R = circulant(5, {1,2})
    examples = {}
    tower = []
    tower_data = []
    T = R
    for k in range(3):
        examples['tower_'+str(k)], data = witness(T)
        row = examples['tower_'+str(k)]['summary']
        check((row['H'],row['group_order'],row['path_orbits'],row['coarse_packing_orbits']) ==
              (15*3**k,5*3**k,3,7*3**k), 'tower numerical formula')
        check(row['stabilizer_histogram'] == {3**k:2*3**k,5*3**k:5*3**k}, 'tower stabilizer formula')
        tower.append(row)
        tower_data.append(data)
        T = order_join(T,C)
    examples['C3'], _ = witness(C)
    examples['C3_to_R5'], _ = witness(order_join(C,R))
    check(examples['C3']['summary']['stabilizer_histogram'] == {3:3}, 'C3 histogram')
    check(examples['C3_to_R5']['summary']['coarse_packing_orbits'] == 21, 'reverse join coarse count')
    restrictions = [restriction_check(tower_data[k-1],tower_data[k]) for k in (1,2)]
    (OUT/'WITNESSES.json').write_text(json.dumps(examples,sort_keys=True,indent=2)+'\n')
    (OUT/'RESTRICTIONS.json').write_text(json.dumps(restrictions,sort_keys=True,indent=2)+'\n')
    records = json.dumps(all_records,sort_keys=True,separators=(',',':'))
    print(json.dumps({'status':'PASS', 'all_labelled_tournaments_n1_through_n5':len(all_records),
                      'extra_trivial_action_controls':identities, 'max_group_order_small':max_group,
                      'all_records_sha256':hashlib.sha256(records.encode()).hexdigest(),
                      'C3':examples['C3']['summary'], 'reverse_join':examples['C3_to_R5']['summary'],
                      'tower':tower, 'restriction_fibres':[{k:v for k,v in x.items() if k!='packing_map'} for x in restrictions],
                      'all_order_statement_is_proved_not_enumerated':True},sort_keys=True,indent=2))

if __name__ == '__main__':
    main()
