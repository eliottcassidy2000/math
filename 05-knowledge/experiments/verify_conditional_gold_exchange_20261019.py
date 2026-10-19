#!/usr/bin/env python3
"""Selected exact controls for conditional Gold/exchange theorems.

No poset census. Two prescribed ten-point witnesses, their two-block ordinal
sums, two tiny controls, and bounded integer tables for algebraic identities.
All guards remain active under python -O. Standard library only.
"""
from functools import lru_cache
from itertools import combinations, permutations, product
from fractions import Fraction
from collections import Counter
import json

stats = Counter()


def require(ok, message):
    if not ok:
        raise RuntimeError(message)
    stats['checks'] += 1


def close(pred):
    p = list(pred)
    for k in range(len(p)):
        for v in range(len(p)):
            if p[v] >> k & 1:
                p[v] |= p[k]
    if any(p[v] >> v & 1 for v in range(len(p))):
        return None
    return tuple(p)


def add(pred, a, b):
    p = list(pred)
    p[b] |= 1 << a
    return close(p)


@lru_cache(None)
def count(pred):
    if pred is None:
        return 0
    n, full = len(pred), (1 << len(pred))-1
    @lru_cache(None)
    def suffix(mask):
        if mask == full:
            return 1
        return sum(suffix(mask | 1 << v) for v in range(n)
                   if not mask >> v & 1 and not pred[v] & ~mask)
    return suffix(0)


def words(pred):
    n = len(pred)
    out = []
    def visit(mask, word):
        if len(word) == n:
            out.append(tuple(word))
            return
        for v in range(n):
            if not mask >> v & 1 and not pred[v] & ~mask:
                visit(mask | 1 << v, word + [v])
    visit(0, [])
    return out


def before(word, pair):
    return word.index(pair[0]) < word.index(pair[1])


def upsets(pred):
    return tuple(sum(1 << j for j in range(len(pred)) if pred[j] >> i & 1)
                 for i in range(len(pred)))


def relabel(pred, sigma):
    q = [0] * len(pred)
    for b in range(len(pred)):
        for a in range(len(pred)):
            if pred[b] >> a & 1:
                q[sigma[b]] |= 1 << sigma[a]
    return tuple(q)


def involution(n, swaps):
    sigma = list(range(n))
    for a, b in swaps:
        sigma[a], sigma[b] = b, a
    return tuple(sigma)


def join_relabel(pred, sigma):
    q = relabel(pred, sigma)
    return close(tuple(a | b for a, b in zip(pred, q)))


def good_candidates(pred):
    n, E = len(pred), count(pred)
    u = upsets(pred)
    def chain(mask):
        vertices = [i for i in range(n) if mask >> i & 1]
        return all(pred[a] >> b & 1 or pred[b] >> a & 1
                   for a, b in combinations(vertices, 2))
    out = []
    for dual, (d, v) in enumerate(((pred, u), (u, pred))):
        for a, b in permutations(range(n), 2):
            if ((d[a] | v[a]) >> b & 1 or d[a] & ~d[b]
                    or not chain(v[b] & ~v[a])):
                continue
            N = count(add(pred, b, a) if dual else add(pred, a, b))
            out.append({'dual': dual, 'a': a, 'b': b, 'numerator': N,
                        'good': 2*N <= E})
    return out


def ordinal_sum(p, q):
    n = len(p)
    return tuple(p) + tuple(((1 << n)-1) | (v << n) for v in q)


def witness(pred, label, expected_E, expected_b, expected_cells, swaps):
    pred = tuple(pred)
    require(close(pred) == pred, label + ': transitivity')
    u = upsets(pred)
    degrees = [len(pred)-1-(a | b).bit_count() for a, b in zip(pred, u)]
    require(max(degrees) == 7, label + ': seven-thin, not six-thin')
    E = count(pred)
    require(E == expected_E, label + ': extension count')
    allwords = words(pred)
    require(len(allwords) == E, label + ': independent enumeration')
    first, second = (3, 5), (2, 4)
    R = add(pred, *first)
    b = count(R)
    require(b == expected_b and 2*b >= E, label + ': heavy branch')
    rw = {w for w in allwords if before(w, first)}
    cells = [sum(before(w, first) == a and before(w, second) == z
                 for w in allwords) for a, z in [(True, True), (True, False),
                                                (False, True), (False, False)]]
    require(cells == expected_cells, label + ': four cells')
    margins = [E-b-cells[0], E-b-cells[1], b-cells[2], b-cells[3]]
    require(min(margins) > 0, label + ': strict Gold')
    eligible = good_candidates(pred)
    require(not any(e['good'] for e in eligible), label + ': no good pair')
    sigma = involution(len(pred), swaps)
    require(all(sigma[sigma[v]] == v for v in range(len(pred))), label + ': involution')
    require(sigma[2] == 4 and sigma[4] == 2, label + ': toggled pair')
    J = join_relabel(R, sigma)
    M = count(J)
    transformed = {tuple(sigma[v] for v in w) for w in rw}
    common = rw & transformed
    require(len(common) == M, label + ': union/intersection count')
    require({tuple(sigma[v] for v in w) for w in common} == common,
            label + ': exchange-invariant intersection')
    positive = sum(before(w, second) for w in common)
    require(2*positive == M, label + ': exact pairing')
    imbalance, defect, tolerance = abs(cells[0]-cells[1]), b-M, 2*E-3*b
    require(imbalance <= defect <= tolerance, label + ': exchange certificate')
    # Two-block control uses DP only; the theorem proves every block count.
    two = ordinal_sum(pred, pred)
    twoR = add(two, *first)
    twosigma = sigma + tuple(range(len(pred), 2*len(pred)))
    require(count(two) == E*E, label + ': ordinal total')
    require(count(twoR) == b*E, label + ': ordinal heavy branch')
    require(count(join_relabel(twoR, twosigma)) == M*E,
            label + ': ordinal overlap')
    return {'label': label, 'pred': pred, 'degrees': degrees, 'E': E,
            'R_pred': R, 'b': b, 'first': first, 'second': second,
            'cells_11_10_01_00': cells, 'gold_margins': margins,
            'sigma': sigma, 'overlap': M, 'imbalance': imbalance,
            'defect': defect, 'tolerance': tolerance,
            'original_signatures': [(a.bit_count(), z.bit_count())
                                    for a, z in zip(pred, u)],
            'branch_signatures': [(a.bit_count(), z.bit_count())
                                  for a, z in zip(R, upsets(R))],
            'eligible_good_orientations': eligible,
            'eligible_count': len(eligible),
            'minimum_eligible_numerator': min(e['numerator'] for e in eligible)}


def main():
    Q = (0, 0, 0, 6, 0, 19, 2, 18, 86, 151)
    Qprime = add(Q, 6, 7)
    first = witness(Q, 'packet117 Q', 11971, 6562,
                    [3281, 3281, 1899, 3510], [(2, 4), (3, 7), (5, 9)])
    require(relabel(tuple(first['R_pred']), first['sigma']) == tuple(first['R_pred']),
            'Q conditioned automorphism')
    require(first['eligible_count'] == 48 and first['minimum_eligible_numerator'] == 6547,
            'Q good-pair regression')
    second = witness(Qprime, 'Q plus 6<7', 7249, 4072,
                     [2221, 1851, 1224, 1953], [(2, 4)])
    require(second['overlap'] == 3192 and second['defect'] == 880,
            'Qprime exchange regression')
    require(second['eligible_count'] == 42 and second['minimum_eligible_numerator'] == 3804,
            'Qprime no-good regression')
    require(len(set(second['original_signatures'])) == 10, 'Qprime rigidity signatures')
    sig = second['branch_signatures']
    coincident = [(a, b) for a, b in combinations(range(10), 2) if sig[a] == sig[b]]
    require(coincident == [(2, 4)], 'Rprime only duplicate signature')
    Rprime = tuple(second['R_pred'])
    require(Rprime[3] >> 2 & 1 and not Rprime[3] >> 4 & 1,
            'fixed vertex3 separates2 and4, so Rprime rigid')

    # Small sharpness control for the exchange-defect coefficient.
    toy = (0, 0, 1)  # 0<2; vertex1 is isolated.
    tau = involution(3, [(0, 1)])
    toy_E = count(toy)
    toy_M = count(join_relabel(toy, tau))
    toy_c = count(add(toy, 0, 1))
    require((toy_E, toy_M, toy_c) == (3, 2, 2), 'sharp exchange toy')
    require(abs(2*toy_c-toy_E) == toy_E-toy_M, 'defect coefficient is sharp')
    toy_counts = [count(add(toy, a, b)) for a, b in combinations(range(3), 2)]
    require(not any(2*toy_E <= 5*n <= 3*toy_E for n in toy_counts),
            'unconditional2/5 strengthening fails on tiny poset')
    require(any(e['good'] for e in good_candidates(toy)), 'tiny hostile has good pair')

    # A symmetric heavy child with p=3/4 still cannot certify Gold.
    hostile = (0, 0, 2, 6, 15, 15)
    H = add(hostile, 1, 0)
    hsigma = involution(6, [(4, 5)])
    require(count(hostile) == 8 and count(H) == 6, 'heavy-symmetry hostile totals')
    require(relabel(H, hsigma) == H and count(add(H, 4, 5)) == 3,
            'heavy-symmetry hostile automorphism')
    require(8-6-3 == -1, 'symmetry without root-mass bound fails')
    require(count((0, 0)) == 2 and 2-1-1 == 0, 'total-child E2 base')
    require(join_relabel((0, 1), involution(2, [(0, 1)])) is None,
            'cyclic union has overlap zero')

    # Arbitrary four-cell count tables, not posets: exact Gold/covariance algebra.
    for values in product(range(5), repeat=4):
        n, k, ell, z = values
        E = sum(values)
        if E == 0:
            continue
        b = n+k
        if 2*b < E:
            continue
        m = n+ell
        delta = E*n-b*m
        interval = 2*b-E <= n <= E-b
        gold = all(E >= b+c for c in (n, k) if c)
        require(gold == interval, 'joint interval equivalence')
        require(interval == (abs(2*n-b) <= 2*E-3*b), 'imbalance equivalence')
        require(abs(b*(2*m-E)+2*delta) == E*abs(2*n-b), 'covariance identity')
        require(all(E >= E-b+c for c in (ell, z)), 'small branch automatic')
        stats['four_cell_tables'] += 1

    # Six-order integer-weight controls for the classical cyclic-triple theorem.
    orders = [(0,1,2), (1,2,0), (2,0,1), (0,2,1), (1,0,2), (2,1,0)]
    arcs = [(0,1), (1,2), (2,0)]
    for weights in product(range(4), repeat=6):
        E = sum(weights)
        if not E:
            continue
        by_order = dict(zip(orders, weights))
        masses = [sum(w for order, w in zip(orders, weights) if before(order, arc))
                  for arc in arcs]
        if not all(2*m > E for m in masses):
            continue
        index = min(range(3), key=lambda i: masses[i])
        x, y, z = index, (index+1)%3, (index+2)%3
        p = masses[index]
        a, c = by_order[(x,y,z)], by_order[(z,x,y)]
        require(3*p <= 2*E, 'cyclic minimum at most2/3')
        require(a >= 2*p-E and c >= 2*p-E, 'cyclic singleton lower bounds')
        chosen = a if a <= E-p else c
        require(2*p-E <= chosen <= E-p, 'cyclic Gold singleton')
        stats['cyclic_weight_laws'] += 1

    # Exact transfer inequality; no assertion that arbitrary children are balanced.
    for E in range(2, 31):
        for b in range((E+1)//2, 3*E//5+1):
            for c in range((b+2)//3, 2*b//3+1):
                require(max(c, b-c) <= E-b, 'conditional1/3 transfer')
                stats['conditional_transfer_cases'] += 1

    # Three already-recorded examples: cheap falsification checks, not a scan.
    prior = [
        ('packet117 twelve-point witness',
         (0,0,0,0,7,13,14,12,206,5,685,535)),
        ('previous nine-point collar witness', (0,0,1,0,8,24,3,5,189)),
    ]
    names = ('b','x','z','y','a','t1','t2','t3')
    edges = [('b','z'),('z','a'),('x','a'),('b','y'),('t1','t2'),('t2','t3')]
    edges += [(v,z) for v in ('t1','t2','t3') for z in ('x','y','a')]
    pred = [0]*len(names)
    for a,b in edges:
        pred[names.index(b)] |= 1 << names.index(a)
    prior.append(('previous eight-point rooted-five obstruction', close(pred)))
    prior_checks = []
    for name, pred in prior:
        pred = tuple(pred)
        E = count(pred)
        U = upsets(pred)
        require(max(len(pred)-1-(d|u).bit_count() for d,u in zip(pred,U)) <= 7,
                name+': seven-thin scope')
        mid = [(a,b,count(add(pred,a,b))) for a,b in combinations(range(len(pred)),2)
               if 2*E <= 5*count(add(pred,a,b)) <= 3*E]
        good = sum(row['good'] for row in good_candidates(pred))
        require(bool(mid) or good > 0, name+': does not falsify proposed dichotomy')
        prior_checks.append({'label':name,'pred':pred,'E':E,'good_orientations':good,
                             'two_fifths_comparison':mid[0] if mid else None})

    print(json.dumps({'status': 'ALL_SELECTED_EXACT_CONTROLS_PASSED',
                      'scope': 'No poset census; two fixed ten-point witnesses, two-block ordinal products, tiny hostile controls, bounded integer probability tables.',
                      'statistics': dict(stats), 'witnesses': [first, second],
                      'selected_prior_dichotomy_controls': prior_checks}, indent=2))


if __name__ == '__main__':
    main()
