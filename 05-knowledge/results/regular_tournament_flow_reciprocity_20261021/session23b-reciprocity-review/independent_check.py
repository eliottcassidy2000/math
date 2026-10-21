#!/usr/bin/env python3
"""Independent exact regular-flow reciprocity controls; no asserts or Git writes.

Uses exhaustive tournament/group checks through order 7, exact rational ranks,
small flow slices, and symbolic rational-function identities. It does not infer
an all-order theorem from finite coefficients.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
from math import factorial, lcm
import json
import sympy as S


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def rank(rows):
    if not rows:
        return 0
    a = [[Fraction(v) for v in row] for row in rows]
    r = 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        z = a[r][c]
        a[r] = [v / z for v in a[r]]
        for i in range(r + 1, len(a)):
            if a[i][c]:
                z = a[i][c]
                a[i] = [u - z * v for u, v in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def cycles(p):
    todo = set(range(len(p)))
    ans = []
    while todo:
        i = min(todo)
        cc = []
        while i in todo:
            todo.remove(i)
            cc.append(i)
            i = p[i]
        ans.append(tuple(cc))
    return ans


def order(p):
    return lcm(*(len(c) for c in cycles(p))) if p else 1


def power(p, k):
    ans = tuple(range(len(p)))
    for _ in range(k % order(p)):
        ans = tuple(p[v] for v in ans)
    return ans


def regular_masks(n):
    edge_index = {e: k for k, e in enumerate(combinations(range(n), 2))}
    r = (n - 1) // 2

    def visit(i, degree, mask):
        if i == n:
            yield mask
            return
        need = r - degree[i]
        if not 0 <= need <= n - i - 1:
            return
        later = range(i + 1, n)
        for winners in combinations(later, need):
            winners = set(winners)
            dd = degree[:]
            mm = mask
            for j in later:
                if j in winners:
                    mm |= 1 << edge_index[i, j]
                else:
                    dd[j] += 1
            if all(dd[j] <= r for j in later):
                yield from visit(i + 1, dd, mm)

    return set(visit(0, [0] * n, 0))


def arcs_for(n, mask):
    return tuple((i, j) if mask >> k & 1 else (j, i)
                 for k, (i, j) in enumerate(combinations(range(n), 2)))


def relabel_mask(arcs, p, index):
    return sum(1 << index[p[u], p[v]] for u, v in arcs if p[u] < p[v])


def circulation(arcs, x, n):
    div = [0] * n
    for (u, v), t in zip(arcs, x):
        div[u] += t
        div[v] -= t
    return not any(div)


def move(x, ap):
    y = [0] * len(x)
    for i, j in enumerate(ap):
        y[j] = x[i]
    return tuple(y)


def arc_action(arcs, p, anti):
    index = {e: k for k, e in enumerate(arcs)}
    return tuple(index[(p[v], p[u]) if anti else (p[u], p[v])]
                 for u, v in arcs)


def fixed_equations(ap):
    out = []
    for i, j in enumerate(ap):
        row = [0] * len(ap)
        row[i] += 1
        row[j] -= 1
        out.append(row)
    return out


def small_flows(n, arcs, bound=8):
    index = {e: k for k, e in enumerate(arcs)}
    directed = []
    for length in range(3, min(n, bound) + 1):
        for seq in permutations(range(n), length):
            if seq[0] != min(seq):
                continue
            es = [(seq[i], seq[(i + 1) % length]) for i in range(length)]
            if all(e in index for e in es):
                v = [0] * len(arcs)
                for e in es:
                    v[index[e]] += 1
                directed.append((length, tuple(v)))
    flow = [set() for _ in range(bound + 1)]
    flow[0].add((0,) * len(arcs))
    for w in range(1, bound + 1):
        for length, c in directed:
            if length <= w:
                for x in flow[w - length]:
                    flow[w].add(tuple(a + b for a, b in zip(x, c)))
    # An independent edge-vector enumeration validates a smaller complete slice.
    direct_bound = min(bound, 8 if n <= 5 else 4)
    direct = [set() for _ in range(direct_bound + 1)]
    examined = 0
    for w in range(direct_bound + 1):
        for support in combinations_with_replacement(range(len(arcs)), w):
            examined += 1
            x = [0] * len(arcs)
            for i in support:
                x[i] += 1
            if circulation(arcs, x, n):
                direct[w].add(tuple(x))
        require(direct[w] == flow[w], 'direct edge/cycle enumeration')
    return flow, len(directed), direct_bound, examined


all_records = []
for n, expected_labelled, expected_classes in [(1, 1, 1), (3, 2, 1), (5, 24, 1), (7, 2640, 3)]:
    unseen = regular_masks(n)
    require(len(unseen) == expected_labelled, 'regular labelled census')
    per = list(permutations(range(n)))
    ix = {e: k for k, e in enumerate(combinations(range(n), 2))}
    full = (1 << len(ix)) - 1
    nclasses = 0
    while unseen:
        mask = min(unseen)
        arcs = arcs_for(n, mask)
        m = len(arcs)
        d = m - n + 1
        r = (n - 1) // 2
        B = [[int(u == v) - int(w == v) for u, w in arcs] for v in range(n)]
        require(rank(B) == n - 1, 'incidence rank')
        require(all(sum(row) == 0 for row in B), 'all ones circulation')
        images = [relabel_mask(arcs, p, ix) for p in per]
        au = [p for p, mm in zip(per, images) if mm == mask]
        an = [p for p, mm in zip(per, images) if mm == full ^ mask]
        orbit = set(images)
        require(len(orbit) * len(au) == factorial(n), 'orbit-stabilizer')
        require(orbit <= regular_masks(n), 'relabeling preserves regularity')
        unseen -= orbit
        nclasses += 1
        require(len(au) % 2 == 1, 'odd automorphism order')
        require(len(an) in (0, len(au)), 'anti coset cardinality')
        actions = []
        histogram = Counter()
        for anti, group in [(False, au), (True, an)]:
            for p in group:
                ap = arc_action(arcs, p, anti)
                ac = cycles(ap)
                C = [[sum(row[i] for i in cc) for cc in ac] for row in B]
                da = len(ac) - rank(C)
                require(da == m - rank(B + fixed_equations(ap)), 'two fixed dimension computations')
                require((da - d) % 2 == 0, 'fixed dimension parity')
                # Check the signed incidence intertwiner on every column.
                eps = -1 if anti else 1
                for i in range(m):
                    for v in range(n):
                        require(B[p[v]][ap[i]] == eps * B[v][i], 'signed incidence action')
                require(all(sum(row) == 0 for row in C), 'fixed-lattice ones translation')
                require(sum(map(len, ac)) == m, 'translation weight')
                o = order(p)
                # Fixed dimension also equals the exact average character.
                chars = []
                for k in range(o):
                    pp, aa = power(p, k), power(ap, k)
                    trE = sum(i == v for i, v in enumerate(aa))
                    trW = (-1 if anti and k % 2 else 1) * (sum(i == v for i, v in enumerate(pp)) - 1)
                    chars.append(trE - trW)
                require(sum(chars) == o * da, 'character fixed dimension')
                if n > 1 and anti:
                    require(o % 4 == 2, 'anti order twice odd')
                    s = o // 2
                    j, h = power(p, s), power(p, 2)
                    aj, ah = arc_action(arcs, j, True), arc_action(arcs, h, False)
                    require(order(j) == 2 and order(h) == s, 'anti odd part decomposition')
                    require(tuple(j[h[v]] for v in range(n)) == tuple(h[j[v]] for v in range(n)), 'commuting decomposition')
                    dj = m - rank(B + fixed_equations(aj))
                    require(dj == r * r, 'anti involution dimension')
                    require(sum(i == v for i, v in enumerate(j)) == 1, 'one fixed vertex')
                    require(sum(i == v for i, v in enumerate(aj)) == r, 'r fixed arcs')
                    require(m - rank(B + fixed_equations(aj) + fixed_equations(ah)) == da, 'fixed intersection')
                    require(all(len(c) == 1 or len(c) % 4 == 2 for c in cycles(p)), 'anti vertex cycle restriction')
                elif not anti:
                    require(o % 2 == 1, 'individual automorphism odd order')
                    if n >= 5 and o > 1:
                        require(da <= d - 2, 'faithful ordinary action dimension drop')
                binary = None
                if len(ac) <= 12:
                    closed, opened = set(), set()
                    for z in product(range(2), repeat=len(ac)):
                        if all(sum(c * x for c, x in zip(row, z)) == 0 for row in C):
                            closed.add(z)
                    for z in product(range(1, 3), repeat=len(ac)):
                        if all(sum(c * x for c, x in zip(row, z)) == 0 for row in C):
                            opened.add(tuple(x - 1 for x in z))
                    require(closed == opened, 'independent small-box interior translation')
                    binary = len(closed)
                histogram[('anti' if anti else 'auto', o, da, len(ac))] += 1
                actions.append((anti, p, ap, binary))
        flows, num_cycles, direct_bound, examined = small_flows(n, arcs)
        autos = [ap for anti, _, ap, _ in actions if not anti]
        antis = [ap for anti, _, ap, _ in actions if anti]
        ordinary, converse = [], []
        for fs in flows:
            reps = {min(move(x, ap) for ap in autos) for x in fs}
            fixsum = sum(move(x, ap) == x for x in fs for ap in autos)
            require(fixsum == len(autos) * len(reps), 'ordinary Burnside')
            ordinary.append(len(reps))
            if antis:
                fixed_reps = sum(min(move(move(x, antis[0]), ap) for ap in autos) == x for x in reps)
                antisum = sum(move(x, ap) == x for x in fs for ap in antis)
                require(antisum == len(autos) * fixed_reps, 'anti-coset Burnside')
                converse.append(fixed_reps)
            for x in fs:
                shifted = tuple(t + 1 for t in x)
                require(circulation(arcs, shifted, n), 'shifted sample circulation')
                for _, _, ap, _ in actions:
                    require((move(x, ap) == x) == (move(shifted, ap) == shifted), 'shift respects fixedness')
        all_records.append({
            'n': n, 'representative_mask': mask, 'arcs': arcs,
            'labelled_isomorphism_orbit_size': len(orbit), 'automorphisms': len(au),
            'anti_automorphisms': len(an), 'dimension': d, 'shift': m,
            'action_histogram': [dict(kind=k[0], order=k[1], fixed_dimension=k[2], arc_orbits=k[3], multiplicity=v)
                                 for k, v in sorted(histogram.items())],
            'binary_box_action_tests': sum(b is not None for _, _, _, b in actions),
            'direct_edge_vector_bound': direct_bound, 'direct_edge_vectors_examined': examined,
            'directed_simple_cycles': num_cycles,
            'raw_flow_coefficients_q0_to_q8': list(map(len, flows)),
            'ordinary_flow_orbit_coefficients_q0_to_q8': ordinary,
            'converse_fixed_orbit_coefficients_q0_to_q8': converse,
        })
    require(nclasses == expected_classes, 'regular isomorphism census')


q = S.symbols('q')
P = 1 + 2*q**3 + 3*q**4 + q**5 + 3*q**6 + q**7 + 3*q**8 + 2*q**9 + q**12
N = 1 + q**4 + q**5 + q**6 + q**7 + q**8 + q**12
R5 = P / ((1-q**3)**3 * (1-q**4)**2 * (1-q**5))
K5 = N / ((1-q**3)*(1-q**5)*(1-q**6)*(1-q**8))
U5 = 1 / (1-q**5)**2
F5 = (R5 + 4*U5)/5
C3 = 1 / (1-q**3)
U7 = 1 / (1-q**7)**3
U3 = (1 + q**3 + q**6 + q**9 + q**12) / ((1-q**3)*(1-q**6)**2*(1-q**9)**2)
V6 = (1 + q**15) / ((1-q**6)*(1-q**12)*(1-q**18))
rational_checks = []
for name, F, sign, shift in [('C3', C3, -1, 3), ('R5_raw', R5, 1, 10),
        ('R5_rotation', U5, 1, 10), ('R5_orbits', F5, 1, 10), ('R5_anti', K5, 1, 10),
        ('R5_extended', (F5+K5)/2, 1, 10), ('Paley7_order7_auto', U7, -1, 21),
        ('Paley7_order3_auto', U3, -1, 21), ('Paley7_order6_anti', V6, -1, 21)]:
    require(S.cancel(F.subs(q, 1/q) - sign*q**shift*F) == 0, name + ' exact reciprocity')
    rational_checks.append(name)
require(S.cancel(U3 - (1-q**15)/((1-q**3)**2*(1-q**6)**2*(1-q**9)**2)) == 0,
        'order3 geometric sum simplification')
v6_sum = (1/(1-q**12) - q**3/(1-q**18))/((1-q**6)*(1-q**3))
require(S.cancel(V6-v6_sum) == 0, 'order6 geometric sum simplification')

# Independently reconstruct both Paley quotient lattices from the named maps.
paley = tuple((u, v) for u in range(7) for v in range(7) if (v-u) % 7 in {1, 2, 4})
paley_fixed_controls = []
for multiplier, anti, expected_orbit_lengths, target in [(2, False, [3]*7, U3), (3, True, [6,6,6,3], V6)]:
    ap = arc_action(paley, tuple(multiplier*i % 7 for i in range(7)), anti)
    ac = cycles(ap)
    require(list(map(len, ac)) == expected_orbit_lengths, 'Paley named arc orbits')
    C = [[sum((paley[i][0] == v) - (paley[i][1] == v) for i in cc) for cc in ac] for v in range(7)]
    if multiplier == 2:
        expected_C = [[3,0,0,0,-3,0,0], [-1,0,1,1,0,-1,0], [-1,0,1,1,0,-1,0],
                      [0,0,-1,-1,1,1,0], [-1,0,1,1,0,-1,0], [0,0,-1,-1,1,1,0], [0,0,-1,-1,1,1,0]]
        coeff = [0]*81
        for a, f, b, g in product(range(27), repeat=4):
            w = 3*(3*a+2*f+b+g)
            if w <= 80:
                coeff[w] += a+f+1
    else:
        expected_C = [[0,0,0,0],[-1,0,2,-1],[-1,0,2,-1],[1,0,-2,1],[-1,0,2,-1],[1,0,-2,1],[1,0,-2,1]]
        coeff = [0]*81
        for a, b, c, dd in product(range(15), repeat=4):
            w = 6*a+6*b+6*c+3*dd
            if w <= 80 and a+dd == 2*c:
                coeff[w] += 1
    require(C == expected_C, 'explicit Paley quotient equations')
    # Compare via denominator convolution, avoiding expensive series expansion.
    num, den = map(S.Poly, S.fraction(S.cancel(target)))
    cc = [int(den.nth(k)) for k in range(den.degree()+1)]
    for w in range(81):
        require(sum(cc[k]*coeff[w-k] for k in range(min(w, len(cc)-1)+1)) == int(num.nth(w)),
                'Paley fixed formula bounded direct coefficient check')
    paley_fixed_controls.append({'multiplier':multiplier,'anti':anti,'orbits':[[paley[i] for i in c] for c in ac],
                                 'incidence_quotient':C,'coefficient_bound':80,'coefficients':coeff})

T4 = 1/((1-q**3)**2*(1-q**4))
T4anti = 1/((1-q**6)*(1-q**4))
T4ext = (T4+T4anti)/2
require(S.cancel(T4.subs(q,1/q) + q**10*T4) == 0, 'nonregular raw true shift')
require(S.cancel(T4.subs(q,1/q) + q**6*T4) != 0, 'nonregular edge-count shift rejected')
require(S.cancel(T4anti.subs(q,1/q) - q**10*T4anti) == 0, 'nonregular anti opposite sign')
require(all(S.cancel(T4ext.subs(q,1/q)-sgn*q**10*T4ext) != 0 for sgn in [-1,1]),
        'nonregular extended neither common sign')

alpha_raw = S.limit((1-q)**6*R5, q, 1)/factorial(6)
alpha_ordinary = alpha_raw/5
alpha_converse = S.limit((1-q)**4*K5, q, 1)/factorial(4)
require((alpha_raw,alpha_ordinary,alpha_converse) == (S.Rational(17,1555200),S.Rational(17,7776000),S.Rational(7,17280)),
        'R5 cumulative leading constants')
require(alpha_converse/alpha_ordinary == S.Rational(3150,17), 'R5 cumulative ratio constant')

print(json.dumps({'status':'PASS','checks_survive_python_O':True,
                  'regular_classes':all_records,'exact_rational_function_checks':rational_checks,
                  'paley_named_fixed_cones':paley_fixed_controls,
                  'hostile_T4':{'raw_dimension':3,'anti_fixed_dimension':2,'true_shift':10,
                               'edgecount_shift_6_rejected':True,'extended_common_sign_rejected':True},
                  'R5_cumulative_asymptotics':{'raw_leading_constant':str(alpha_raw),
                      'ordinary_leading_constant':str(alpha_ordinary),'converse_leading_constant':str(alpha_converse),
                      'ratio_leading_constant':str(alpha_converse/alpha_ordinary),'ratio_exponent':-2}},
                 sort_keys=True,indent=2))
