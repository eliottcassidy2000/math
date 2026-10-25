"""Exact action, positive-jet, and full-automorphism hostile controls."""
from itertools import combinations, permutations, product
from math import comb
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def add(a, b):
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(max(len(a), len(b)))]


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def power(a, exponent):
    out = [1]
    for _ in range(exponent):
        out = mul(out, a)
    return out


def cycles(p):
    unseen = set(range(len(p)))
    sizes = []
    while unseen:
        v = min(unseen)
        size = 0
        while v in unseen:
            unseen.remove(v)
            size += 1
            v = p[v]
        sizes.append(size)
    return sizes


def burnside_polynomial(group):
    total = [0]
    for p in group:
        term = [1]
        for length in cycles(p):
            factor = [0] * (length + 1)
            factor[0] = factor[-1] = 1
            term = mul(term, factor)
        total = add(total, term)
    need(all(c % len(group) == 0 for c in total), 'integer orbit polynomial')
    return [c // len(group) for c in total]


def shifted(poly):
    n = len(poly) - 1
    return [sum((-1) ** (n + k) * poly[k] * comb(k, j) for k in range(j, n + 1))
            for j in range(n + 1)]


def positive_cycle_formula(group):
    total = [0]
    for p in group:
        term = [1]
        for length in cycles(p):
            need(length % 2 == 1, 'odd cycle hypothesis')
            term = mul(term, [0] + [comb(length, i) for i in range(1, length + 1)])
        total = add(total, term)
    need(all(c % len(group) == 0 for c in total), 'integer positive jet')
    return [c // len(group) for c in total]


def subset_orbits(group):
    n = len(group[0])
    representatives = {}
    for mask in range(1 << n):
        key = min(sum(1 << p[i] for i in range(n) if mask >> i & 1) for p in group)
        representatives[key] = mask.bit_count()
    out = [0] * (n + 1)
    for size in representatives.values():
        out[size] += 1
    return out


def cayley25():
    vertices = list(product(range(5), repeat=2))
    connection = {(0, 3), (0, 4), (1, 0), (3, 0), (3, 1), (3, 2),
                  (3, 3), (3, 4), (4, 1), (4, 2), (4, 3), (4, 4)}
    arcs = {(i, j) for i, x in enumerate(vertices) for j, y in enumerate(vertices)
            if tuple((b - a) % 5 for a, b in zip(x, y)) in connection}
    need(all(((a, b) in arcs) != ((b, a) in arcs) for a, b in combinations(range(25), 2)), 'Cayley tournament')
    need(all(sum((a, b) in arcs for b in range(25)) == 12 for a in range(25)), 'regular base')
    colors = [int(i == 0) for i in range(25)]
    trace = [len(set(colors))]
    while len(set(colors)) < 25:
        signatures = [(colors[i], tuple(sum((i, j) in arcs for j in range(25) if colors[j] == c)
                                        for c in range(max(colors) + 1)),
                       tuple(sum((j, i) in arcs for j in range(25) if colors[j] == c)
                             for c in range(max(colors) + 1))) for i in range(25)]
        names = {s: k for k, s in enumerate(sorted(set(signatures)))}
        new = [names[s] for s in signatures]
        need(len(set(new)) > len(set(colors)), 'rooted refinement must individualize')
        colors = new
        trace.append(len(set(colors)))
    need(trace == [2, 3, 15, 21, 25], 'refinement certificate')
    index = {v: i for i, v in enumerate(vertices)}
    group = [tuple(index[tuple((a + b) % 5 for a, b in zip(v, shift))] for v in vertices)
             for shift in vertices]
    need(all(all((p[a], p[b]) in arcs for a, b in arcs) for p in group), 'all translations preserve arcs')
    # Rooted discrete refinement proves no further stabilizer; orbit-stabilizer gives full group.
    inflated = {(2*a+i, 2*b+j) for a, b in arcs for i in range(2) for j in range(2)}
    inflated |= {(2*a, 2*a+1) for a in range(25)}
    homogeneous = {(a, b) for a, b in combinations(range(50), 2)
                   if all(((a, c) in inflated) == ((b, c) in inflated)
                          for c in range(50) if c not in (a, b))}
    need(homogeneous == {(2*a, 2*a+1) for a in range(25)}, 'exact canonical TT2 blocks')
    return group, trace


def main():
    group25, trace = cayley25()
    p25 = burnside_polynomial(group25)
    q25 = shifted(p25)
    need(q25 == positive_cycle_formula(group25), 'two positive-jet routes')
    need(all(c >= 0 for c in q25), 'nonnegative jet')
    need(next(i for i, c in enumerate(q25) if c) == 5, 'depth five')
    need(sum(p[0] == 0 for p in group25) == 1, 'regular full action')
    elementary = []
    for r in (1, 2, 3):
        q = 5 ** r
        h = q // 5
        numerator = add([comb(q, j) for j in range(q+1)],
                        [(q-1)*comb(h, j//5) if j % 5 == 0 else 0 for j in range(q+1)])
        need(all(c % q == 0 for c in numerator), 'elementary orbit coefficient integrality')
        poly = [c // q for c in numerator]
        jet = shifted(poly)
        positive = add([0] * q + [1], [(q-1)*c for c in power([0,5,10,10,5,1],h)])
        need(all(c % q == 0 for c in positive), 'elementary positive integrality')
        need(jet == [c // q for c in positive], 'elementary closed formula')
        need(next(i for i,c in enumerate(jet) if c) == h, 'elementary depth')
        elementary.append({'rank': r, 'gaps': q, 'depth': h, 'leading_coefficient': jet[h]})

    arcs7 = {(a,b) for a in range(7) for b in range(7) if (b-a)%7 in (1,2,4)}
    full7 = [p for p in permutations(range(7)) if all((p[a],p[b]) in arcs7 for a,b in arcs7)]
    need(len(full7) == 21, 'full Paley7 automorphisms')
    p7 = burnside_polynomial(full7)
    need(p7 == subset_orbits(full7), 'literal small cut orbits')
    need(shifted(p7) == positive_cycle_formula(full7), 'Paley positive jet')
    doubled = [tuple(p) + (7,) + tuple(8+p[i] for i in range(7)) for p in full7]
    p15 = burnside_polynomial(doubled)
    naive = mul([1,1], mul(p7,p7))
    need(p15[2] == 7 and naive[2] == 5, 'independent-copy multiplication hostile')
    need(next(i for i,c in enumerate(shifted(p15)) if c) == 3, 'doubled depth')

    # A regular C3^2 subgroup is insufficient: C3[C3] has a nine-cycle automorphism.
    v9 = list(product(range(3),repeat=2))
    index9 = {v:i for i,v in enumerate(v9)}
    arcs9 = {(i,j) for i,(a,b) in enumerate(v9) for j,(c,d) in enumerate(v9)
             if (a!=c and (c-a)%3==1) or (a==c and (d-b)%3==1)}
    odometer = tuple(index9[((a+1)%3,(b+int(a==2))%3)] for a,b in v9)
    need(cycles(odometer) == [9] and all((odometer[a],odometer[b]) in arcs9 for a,b in arcs9),
         'full-action subgroup hostile')
    print(json.dumps({'status':'PASS','explicit_C5_squared_TT2':{'base_vertices':25,'inflated_vertices':50,
                      'rooted_refinement_class_counts':trace,'full_automorphism_order':25,
                      'invariant_gap_toggle_rank':1,'rational_cancellation_depth':5},
                      'closed_family_formula_controls':elementary,
                      'F21_doubling_t_squared_coefficients':{'correct':p15[2],'independent_copy_error':naive[2]},
                      'C3_squared_subgroup_hostile':{'subgroup_depth':3,'full_action_depth':1},
                      'scope':'exact concrete TRR certificate; infinite existence uses cited Godsil theorem'},indent=2))


if __name__ == '__main__':
    main()
