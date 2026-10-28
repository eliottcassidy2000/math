#!/usr/bin/env python3
"""EXPLORATORY (post hoc, not pre-registered).  kind-pasteur-2026-10-27-S8.

(a) Tier analysis of the non-linear carriers n <= 7: defect delta(g) = chi(G^g) - 1 for odd-order elements and for
    involutions; is non-weak explained by an odd-order defect?
(b) Count tier: for the non-weak carriers, does the marked count still equal |Tour(X)/A| for some A-set X?
(c) The theta_p family: p directed u -> w paths of length k, plus p-1 isolated vertices (n = pk + 1): the Z_p mark
    has no matching X (proof check), so never weakly trijective.
(d) K_4 (transitive) + 3K_1: linear for the subgroup Z_3, not for Z_2 or the full S_3.
(e) Functional averaging identity: sum_T sum_r prod_{v != r} s_v(T) = F_n 2^C(n-1,2), F_n = number of rooted
    loop- and 2-cycle-free functional digraphs on [n].
"""
import json, itertools
from math import comb, factorial
from collections import Counter
from carrier_core import automorphisms, Group, signed_types, coset_action, partitions_into, marks


def elem_order(g):
    n = len(g); k = 1; h = list(g)
    while h != list(range(n)):
        h = [g[x] for x in h]; k += 1
    return k


def delta(n, edges, g):
    F = {v for v in range(n) if g[v] == v}
    return len(F) - sum(1 for e in edges if e <= F) - 1


def linear_for_group(n, arcs, elems, j=1):
    G = Group(elems)
    edges = {frozenset(a) for a in arcs}
    Pi = [(a, b) for a, b in itertools.combinations(range(n), 2) if frozenset((a, b)) not in edges]
    actV = lambda i, x: G.el[i][x]
    piT = signed_types(G, actV, Pi)
    subs = list(G.subgroups()); reps = {}
    for H in subs:
        reps.setdefault(G.canon_sub(H), H)
    parts = [(cid, len(G.el) // len(H)) for cid, H in reps.items()]
    for combo in partitions_into(n - j, parts):
        tables = []; off = 0
        for cid in combo:
            m, tab = coset_action(G, reps[cid]); tables.append((off, m, tab)); off += m
        def actX(i, x, tables=tables):
            for o, m, tab in tables:
                if o <= x < o + m:
                    return o + tab[i][x - o]
        if signed_types(G, actX, list(itertools.combinations(range(off), 2))) == piT:
            return True, [len(G.el) // len(reps[c]) for c in combo]
    return False, None


def tour_quotient_counts(n, arcs):
    """marked count, and the set of |Tour(X)/A| over all A-sets X of size n-1."""
    A = automorphisms(n, arcs); G = Group(A)
    edges = {frozenset(a) for a in arcs}
    Pi = [(a, b) for a, b in itertools.combinations(range(n), 2) if frozenset((a, b)) not in edges]
    actV = lambda i, x: G.el[i][x]
    cyc = [G.closure([i]) for i in range(len(G.el))]
    def count(act, pairs):
        tot = 0
        for K in cyc:
            m = marks(G, act, pairs, [K])[K]
            tot += 0 if m is None else 2 ** m
        return tot // len(G.el)
    marked = count(actV, Pi)
    subs = list(G.subgroups()); reps = {}
    for H in subs:
        reps.setdefault(G.canon_sub(H), H)
    parts = [(cid, len(G.el) // len(H)) for cid, H in reps.items()]
    vals = set()
    for combo in partitions_into(n - 1, parts):
        tables = []; off = 0
        for cid in combo:
            m, tab = coset_action(G, reps[cid]); tables.append((off, m, tab)); off += m
        def actX(i, x, tables=tables):
            for o, m, tab in tables:
                if o <= x < o + m:
                    return o + tab[i][x - o]
        vals.add(count(actX, list(itertools.combinations(range(off), 2))))
    return marked, sorted(vals)


def theta_p(p, k):
    u, w = 0, 1; arcs = []; nxt = 2
    for _ in range(p):
        prev = u
        for _ in range(k - 1):
            arcs.append((prev, nxt)); prev = nxt; nxt += 1
        arcs.append((prev, w))
    n = nxt + (p - 1)
    return n, arcs


def theta_mark_check(p, k):
    """Z_p rotation of the p paths: the Pi-mark exponent vs every admissible X cycle type (x fixed, c p-cycles)."""
    n, arcs = theta_p(p, k)
    edges = {frozenset(a) for a in arcs}
    v = 2 + (p - 1); e = 0                       # fixed: u, w and the p-1 isolated vertices; no fixed edges
    fixed_pairs_Pi = comb(v, 2) - e
    orbits_Pi = (comb(n - 1, 2) + (p - 1) * fixed_pairs_Pi) // p
    ok = []
    for x in range(0, n):
        if (n - 1 - x) % p:
            continue
        orbits_X = (comb(n - 1, 2) + (p - 1) * comb(x, 2)) // p
        if orbits_X == orbits_Pi:
            ok.append(x)
    return n, orbits_Pi, ok


def functional_count(n):
    m = n - 1
    s = sum((-1) ** k * factorial(m) // (factorial(k) * 2 ** k * factorial(m - 2 * k)) * m ** (m - 2 * k)
            for k in range(m // 2 + 1))
    return n * s


def functional_identity(n):
    pairs = list(itertools.combinations(range(n), 2)); tot = 0
    for mask in range(1 << len(pairs)):
        s = [0] * n
        for i, (a, b) in enumerate(pairs):
            if (mask >> i) & 1:
                s[a] += 1
            else:
                s[b] += 1
        for r in range(n):
            pr = 1
            for v in range(n):
                if v != r:
                    pr *= s[v]
            tot += pr
    return tot, functional_count(n) * 2 ** comb(n - 1, 2)


def main():
    out = {}
    census = json.load(open('carrier_census_n7.json'))
    rows = []
    for n_str, rs in census.items():
        n = int(n_str)
        for r in rs:
            if r['linear'] is False:
                arcs = [tuple(a) for a in r['arcs']]
                edges = {frozenset(a) for a in arcs}
                A = automorphisms(n, arcs)
                odd = [delta(n, edges, g) for g in A if elem_order(g) % 2 == 1 and elem_order(g) > 1]
                inv = [delta(n, edges, g) for g in A if elem_order(g) == 2]
                rows.append({'n': n, 'weak': r['weak'], 'odd_defect': any(d != 0 for d in odd),
                             'inv_defect': any(d != 0 for d in inv), 'arcs': arcs})
    tab = Counter((x['weak'], x['odd_defect'], x['inv_defect']) for x in rows)
    out['a_tiers(weak, odd_defect, involution_defect)'] = {str(k): v for k, v in tab.items()}
    print('(a) non-linear carriers n<=7 by (weak, odd-order defect, involution defect):', dict(tab), flush=True)
    nonweak = [x for x in rows if not x['weak']]
    out['b_count_tier'] = []
    for x in nonweak:
        marked, vals = tour_quotient_counts(x['n'], x['arcs'])
        out['b_count_tier'].append({'arcs': x['arcs'], 'marked': marked, 'tour_quotients': vals,
                                    'count_matches_some_X': marked in vals})
        print('(b) non-weak', x['arcs'], 'marked', marked, 'Tour(X)/A values', vals, 'match:', marked in vals, flush=True)
    out['c_theta'] = []
    for p, k in ((3, 2), (3, 3), (3, 4), (5, 2), (5, 3), (7, 2)):
        n, orb, ok = theta_mark_check(p, k)
        out['c_theta'].append({'p': p, 'k': k, 'n': n, 'Z_p_orbits_on_Pi': orb, 'admissible_x': ok})
        print(f'(c) theta_{p}, k={k}, n={n}: Z_p-orbits on Pi = {orb}; X fixed-point counts matching = {ok}', flush=True)
    # (d) K4 (transitive on 0..3) + isolated 4,5,6
    arcs = [(a, b) for a, b in itertools.combinations(range(4), 2)]
    A = automorphisms(7, arcs)
    z3 = [g for g in A if g[4:] in ((4, 5, 6), (5, 6, 4), (6, 4, 5))]
    z2 = [g for g in A if g[4:] in ((4, 5, 6), (5, 4, 6))]
    out['d_K4_3K1'] = {'full_S3': linear_for_group(7, arcs, A), 'Z3': linear_for_group(7, arcs, z3),
                       'Z2': linear_for_group(7, arcs, z2)}
    print('(d) K4+3K1 linear for full S3 / Z3 / Z2:', out['d_K4_3K1'], flush=True)
    out['e_functional_identity'] = {}
    for n in range(3, 7):
        lhs, rhs = functional_identity(n)
        out['e_functional_identity'][n] = {'sum_T_sum_r_prod': lhs, 'F_n_times_2^C(n-1,2)': rhs, 'F_n': functional_count(n)}
        print(f'(e) n={n}: sum_T sum_r prod s_v = {lhs}; F_n 2^C(n-1,2) = {rhs}; F_n = {functional_count(n)}', flush=True)
    json.dump(out, open('tiers_and_sporadics.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
