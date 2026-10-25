"""Standalone exact orbit/cancellation and literal tournament checks (stdlib only)."""
from collections import Counter
from itertools import product
from math import comb, prod
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def cycles(p):
    unseen = set(range(len(p)))
    result = []
    while unseen:
        start = min(unseen)
        orbit = []
        v = start
        while v in unseen:
            unseen.remove(v)
            orbit.append(v)
            v = p[v]
        need(v == start, 'permutation cycle closes')
        result.append(tuple(orbit))
    return result


def act(mask, p):
    return sum(1 << p[i] for i in range(len(p)) if mask >> i & 1)


def multiply(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i+j] += x*y
    return out


def cycle_polynomial(p):
    out = [1]
    for c in cycles(p):
        factor = [0] * (len(c)+1)
        factor[0] = factor[-1] = 1
        out = multiply(out, factor)
    return out


def shifted(p):
    return [sum(p[k] * comb(k, j) * (-1)**(k-j)
                for k in range(j, len(p))) for j in range(len(p))]


def orbit_components(group):
    n = len(group[0])
    unseen = set(range(n))
    result = []
    while unseen:
        o = frozenset(p[min(unseen)] for p in group)
        need(o <= unseen, 'point orbits partition')
        result.append(o)
        unseen -= o
    return result


def inspect_action(group, odd=True, cube_check=True):
    group = sorted(set(group))
    n = len(group[0])
    ident = tuple(range(n))
    need(ident in group, 'identity')
    universe = set(group)
    for p in group:
        need(sorted(p) == list(ident), 'permutation')
        for q in group:
            need(tuple(p[q[i]] for i in range(n)) in universe, 'group closure')
    reps = {min(act(mask, p) for p in group) for mask in range(1 << n)}
    counts = Counter(mask.bit_count() for mask in reps)
    polynomial = [counts[k] for k in range(n+1)]
    burnside = [0]*(n+1)
    for p in group:
        for i, a in enumerate(cycle_polynomial(p)):
            burnside[i] += a
    need(all(x % len(group) == 0 for x in burnside), 'Burnside integrality')
    need(polynomial == [x//len(group) for x in burnside], 'brute rank census equals cycle-index')
    local = shifted(polynomial)
    order = next(i for i, x in enumerate(local) if x)
    minimum = min(len(cycles(p)) for p in group)
    components = orbit_components(group)
    r = len(components)
    cycle_histogram = Counter(tuple(sorted(map(len, cycles(p)))) for p in group)
    if odd:
        need(len(group) % 2 == 1, 'group odd order')
        need(all(len(c) % 2 for p in group for c in cycles(p)), 'every cycle odd')
        need(order == minimum, 'root multiplicity equals minimum cycle count')
        numerator = sum(prod(map(len, cycles(p))) for p in group if len(cycles(p)) == minimum)
        need(numerator % len(group) == 0, 'leading coefficient is integral')
        need(local[order] == numerator//len(group) > 0, 'positive leading coefficient formula')
        need(order >= r and (order-r) % 2 == 0, 'orbit bound and parity')
        sync = any(all(sum(c[0] in o for c in cycles(p)) == 1 for o in components) for p in group)
        need((order == r) == sync, 'simultaneous cyclic transitivity criterion')
        need(len(reps) % (1 << order) == 0, 'extra integral divisibility at t=1')
        need(polynomial == polynomial[::-1], 'complement reciprocal symmetry')
        signed_moments = [sum((-1)**(n-mask.bit_count())*comb(mask.bit_count(), j)
                              for mask in reps) for j in range(order+1)]
        need(signed_moments == [0]*order + [local[order]], 'first clone signed binomial moment')
        if r == 1 and n > 1:
            need(3*order <= n, 'sharp transitive bound')
        if cube_check:
            invariant_masks = [sum(1 << x for j, o in enumerate(components)
                                   if selector >> j & 1 for x in o)
                               for selector in range(1 << r)]
            need(len(set(invariant_masks)) == 1 << r, 'all invariant translations')
            for mask in reps:
                cube = {min(act(mask ^ f, p) for p in group) for f in invariant_masks}
                need(len(cube) == 1 << r, 'Boolean toggle action free on subset orbits')
                for j, o in enumerate(components):
                    f = sum(1 << x for x in o)
                    need(((mask ^ f).bit_count()-mask.bit_count()) % 2, 'generator reverses cut parity')
    return {'points': n, 'group_order': len(group), 'point_orbits': r,
            'root_order': order, 'min_cycles': minimum, 'leading_coefficient': local[order],
            'subset_orbits': len(reps), 'polynomial': polynomial,
            'cycle_types': {','.join(map(str,k)):v for k,v in sorted(cycle_histogram.items())}}


def regular_abelian(moduli):
    points = list(product(*(range(m) for m in moduli)))
    index = {x:i for i,x in enumerate(points)}
    return [tuple(index[tuple((x[j]+a[j]) % moduli[j] for j in range(len(moduli)))]
                  for x in points) for a in points]


def projection_action():
    forms = ((1,0),(0,1),(1,1),(1,2))
    return [tuple(3*i+(c+u*a+v*b)%3 for i,(u,v) in enumerate(forms) for c in range(3))
            for a,b in product(range(3),repeat=2)]


def core(k, forms=None):
    points = list(product(range(3), repeat=k))
    if forms is None:
        forms = [tuple(int(i == j) for j in range(k)) for i in range(k)]
    m = len(forms)
    labels = [('x', a) for a in points] + [('y', i, c) for i in range(m) for c in range(3)]
    arcs = set()
    def orient(u,v):
        if u[0] == v[0] == 'x':
            first = next(j for j in range(k) if u[1][j] != v[1][j])
            return (v[1][first]-u[1][first]) % 3 == 1
        if u[0] == v[0] == 'y':
            return u[1] < v[1] if u[1] != v[1] else (v[2]-u[2]) % 3 == 1
        if u[0] == 'x':
            return sum(a*b for a,b in zip(u[1],forms[v[1]])) % 3 == v[2]
        return not orient(v,u)
    for i in range(len(labels)):
        for j in range(i+1,len(labels)):
            arcs.add((i,j) if orient(labels[i],labels[j]) else (j,i))
    return labels, frozenset(arcs)


def inflate(labels, arcs, kind):
    blocks = []
    expanded = []
    for i,u in enumerate(labels):
        block = tuple(range(len(expanded), len(expanded)+(2 if u[0] == kind else 1)))
        blocks.append(block)
        expanded.extend((u,j) for j in range(len(block)))
    out = set()
    for block in blocks:
        for i,u in enumerate(block):
            for v in block[i+1:]:
                out.add((u,v))
    for i,j in arcs:
        out.update((u,v) for u in blocks[i] for v in blocks[j])
    return expanded, frozenset(out), blocks


def homogeneous_pairs(n, arcs):
    return [(u,v) for u in range(n) for v in range(u+1,n)
            if all(((u,w) in arcs) == ((v,w) in arcs) for w in range(n) if w not in (u,v))]


def automorphisms(n, arcs):
    """All arc-preserving bijections; score pruning plus exact adjacency only."""
    adj = [set(v for u,v in arcs if u == i) for i in range(n)]
    scores = list(map(len,adj))
    domains = [{v for v in range(n) if scores[v] == scores[u]} for u in range(n)]
    result = []
    def search(mapping, used):
        if len(mapping) == n:
            result.append(tuple(mapping[i] for i in range(n)))
            return
        choices = []
        for u in range(n):
            if u in mapping:
                continue
            candidates = [v for v in domains[u]-used
                          if all((s in adj[u]) == (t in adj[v]) for s,t in mapping.items())]
            if not candidates:
                return
            choices.append((len(candidates),u,candidates))
        _,u,candidates = min(choices)
        for v in candidates:
            mapping[u] = v
            search(mapping,used|{v})
            del mapping[u]
    search({},set())
    need(all(all((p[u],p[v]) in arcs for u,v in arcs) for p in result), 'literal full automorphisms')
    return result


def inspect_tournament(k, forms=None, inflate_kind='x', full=True):
    labels, arcs = core(k, forms)
    base_pairs = homogeneous_pairs(len(labels),arcs)
    need(not base_pairs, 'base tournament pair-free')
    expanded, tarcs, blocks = inflate(labels,arcs,inflate_kind)
    n = len(expanded)
    pairs = homogeneous_pairs(n,tarcs)
    expected = [tuple(b) for b in blocks if len(b)==2]
    need(pairs == expected, 'canonical modules exactly intended clone pairs and singletons')
    group = automorphisms(n,tarcs)
    index = {g:i for i,g in enumerate(pairs)}
    gap_action = [tuple(index[(p[u],p[v])] for u,v in pairs) for p in group]
    need(len(group) == 3**k, 'full tournament automorphism order')
    need(len(set(gap_action)) == len(group), 'gap action faithful in examples')
    if full:
        answer = inspect_action(gap_action)
    else:
        answer = {'points':len(pairs), 'point_orbits':len(orbit_components(gap_action)),
                  'min_cycles':min(len(cycles(p)) for p in gap_action)}
    answer.update({'vertices':n, 'base_vertices':len(labels), 'automorphisms':len(group),
                   'base_score_classes':dict(Counter(sum(u==v for u,_ in arcs) for v in range(len(labels)))),
                   'literal_homogeneous_pairs':len(pairs)})
    return answer


def divide_one_plus_t(p):
    need(len(p) >= 2, 'nonconstant factor dividend')
    q = [p[0]]
    for x in p[1:-1]:
        q.append(x-q[-1])
    need(q[-1] == p[-1], 'exact monic factor division')
    return q


def exterior_permutation_polynomial(group):
    n = len(group[0])
    total = [0]*(n+1)
    for p in group:
        term = [1]
        for c in cycles(p):
            d = len(c)
            factor = [0]*(d+1)
            factor[0] = 1
            factor[-1] = 1 if d % 2 else -1
            term = multiply(term,factor)
        for j,a in enumerate(term):
            total[j] += a
    need(all(a % len(group) == 0 for a in total), 'exterior Molien integrality')
    result = [a//len(group) for a in total]
    need(all(a >= 0 for a in result), 'exterior invariant dimensions nonnegative')
    return result


def main():
    out = {}
    out['empty'] = inspect_action([()])
    out['trivial_five'] = inspect_action([tuple(range(5))])
    out['cyclic_three'] = inspect_action(regular_abelian((3,)))
    out['cyclic_nine'] = inspect_action(regular_abelian((9,)))
    out['regular_c3_squared'] = inspect_action(regular_abelian((3,3)))
    out['four_projection_orbits'] = inspect_action(projection_action())
    out['even_hostile'] = inspect_action([(0,1),(1,0)],odd=False)
    need(out['even_hostile']['root_order'] == 0 and out['even_hostile']['min_cycles'] == 1,
         'even order refutes unqualified root formula')
    out['tournament_24'] = inspect_tournament(2)
    need(out['tournament_24']['polynomial'] == out['regular_c3_squared']['polynomial'], 'actual regular gap action consequence')
    out['tournament_63_structure'] = inspect_tournament(3, full=False)
    need(out['tournament_63_structure']['min_cycles'] == 9, 'k=3 next-scale minimum')
    forms = ((1,0),(0,1),(1,1),(1,2))
    out['tournament_33'] = inspect_tournament(2,forms,inflate_kind='y')
    need(out['tournament_33']['polynomial'] == out['four_projection_orbits']['polynomial'], 'actual synchronization hostile consequence')
    labels, arcs = core(1)
    out['k1_boundary'] = {'vertices':len(labels), 'homogeneous_pairs':homogeneous_pairs(len(labels),arcs)}
    need(len(out['k1_boundary']['homogeneous_pairs']) == 3, 'k=1 base fails pair-free hypothesis')
    need(out['regular_c3_squared']['root_order'] == 3 and out['regular_c3_squared']['leading_coefficient'] == 24, 'strict root example')
    need(out['four_projection_orbits']['root_order'] == 6 and out['four_projection_orbits']['point_orbits'] == 4, 'synchronization failure')
    p = out['regular_c3_squared']['polynomial']
    need(exterior_permutation_polynomial(regular_abelian((3,3))) == p,
         'odd subset orbit polynomial equals exterior Molien')
    quotients = []
    for i in range(3):
        p = divide_one_plus_t(p)
        quotients.append(p)
    need(all(a >= 0 for a in quotients[0]), 'rank-r exterior factor stays coefficientwise nonnegative')
    need(any(a < 0 for a in quotients[-1]), 'rank-d positive exterior factor is impossible')
    out['regular_c3_squared_factor_quotients'] = quotients
    even_exterior = exterior_permutation_polynomial([(0,1),(1,0)])
    need(even_exterior == [1,1,0], 'even action needs signed subset stabilizers')
    out['even_hostile_exterior_polynomial'] = even_exterior
    # C3 acting on C^3 by scalar primitive cube root: invariants only degrees 0,3.
    complex_h = [comb(3,j) if j % 3 == 0 else 0 for j in range(4)]
    need(complex_h == [1,0,0,1] and shifted(complex_h)[0] == 0,
         'complex representation root-order theorem needs real hypothesis')
    out['complex_hostile'] = {'exterior_polynomial':complex_h, 'min_fixed_dimension':0, 'root_order':1}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
