#!/usr/bin/env python3
"""Exact controls for the degree-uniform marked three-cusp D4 exclusion.

The analytic parity argument, not a finite ambient census, proves the theorem.
No external packages. Checks remain active under python -O.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations
import json

GATES = 0


def check(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def mul(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inv(p):
    z = [0] * len(p)
    for i, j in enumerate(p):
        z[j] = i
    return tuple(z)


def conj(p, q):
    return mul(mul(p, q), inv(p))


def fixed(p):
    return {i for i, j in enumerate(p) if i == j}


def support(p):
    return set(range(len(p))) - fixed(p)


def image(p, s):
    return {p[i] for i in s}


def braid(p, q):
    return mul(mul(p, q), p) == mul(mul(q, p), q)


def cycle_type(p):
    todo = set(range(len(p)))
    sizes = []
    while todo:
        i = min(todo)
        j = i
        n = 0
        while j in todo:
            todo.remove(j)
            n += 1
            j = p[j]
        sizes.append(n)
    return tuple(sorted(sizes))


def partitions(n, least=1):
    if n == 0:
        yield ()
    for first in range(least, n + 1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def canonical(parts):
    p = list(range(sum(parts)))
    pos = 0
    for n in parts:
        for i in range(pos, pos + n):
            p[i] = pos + (i - pos + 1) % n
        pos += n
    return tuple(p)


def subsets(s):
    ss = sorted(s)
    for k in range(len(ss) + 1):
        for z in combinations(ss, k):
            yield set(z)


def p_add(*items):
    out = Counter()
    for scale, p in items:
        for v, coefficient in p.items():
            out[v] += scale * coefficient
    return {v: c for v, c in out.items() if c}


def affine(**kwargs):
    return dict(kwargs)


def main():
    # Formal affine coefficient arithmetic in the independent symbols F,C,k,delta.
    upper = affine(F=3, one=1, k=-1)
    union = affine(F=2, C=2, k=-1, delta=-1)
    central = affine(C=2, F=-1)
    generic = affine(F=3, C=3, k=-3, delta=3)
    check(p_add((1, upper), (-1, union)) ==
          affine(F=1, C=-2, one=1, delta=1), 'union comparison identity')
    check(p_add((1, upper), (-1, central)) ==
          affine(F=4, C=-2, one=1, k=-1), 'central comparison identity')
    check(p_add((1, upper), (-1, generic)) ==
          affine(C=-3, one=1, k=2, delta=-3), 'generic comparison identity')
    # The parity contradictions are affine identities for unrestricted integer m.
    check(p_add((1, affine(m=3, one=3)),
                (-1, affine(m=2, one=-1))) == affine(m=1, one=4),
          'odd lower k exceeds upper by m+4')
    check(p_add((1, affine(m=3, one=-2)),
                (-1, affine(m=2, one=1))) == affine(m=1, one=-3),
          'even lower k exceeds upper when m>3')

    # Exhaustive integer universe AFTER the analytic all-degree reduction m<=3.
    integer_rows = []
    for m in (1, 2, 3):
        t = 2 * m
        k_max = 2 * m + 1
        delta_max = (2 * k_max - 3 * m + 1) // 3
        for delta in range(max(0, m - 1), delta_max + 1):
            for k in range(1, k_max + 1):
                if 2 * k < 3 * m + 3 * delta - 1:
                    continue
                d = t + k + delta
                upper_w = 3 * m + 1 - k
                lower_w = max(0, 3 * t - d, 3 * m - t, 3 * (d - 2 * k))
                check(lower_w <= upper_w, 'integer-row comparison')
                integer_rows.append([t, delta, k, d, lower_w, upper_w])
    expected = [[2, 0, 1, 3, 3, 3], [2, 0, 2, 4, 2, 2],
                [2, 0, 3, 5, 1, 1], [2, 1, 3, 6, 1, 1],
                [4, 1, 4, 9, 3, 3], [4, 1, 5, 10, 2, 2],
                [6, 2, 7, 15, 3, 3]]
    check(integer_rows == expected, 'complete reduced integer universe')
    check([r for r in integer_rows if r[0] == 6] == [[6, 2, 7, 15, 3, 3]],
          'unique moved-six numerical survivor')

    # Hostile: all three set bounds can be sharp before commuting is imposed.
    bset = set(range(6))
    leaves = [{0, 1, 2, 6, 7, 8}, {2, 3, 4, 9, 10, 11},
              {0, 4, 5, 12, 13, 14}]
    intersections = [len(leaves[i] & leaves[j]) for i, j in combinations(range(3), 2)]
    check([len(x) for x in leaves] == [6, 6, 6], 'formal supports have size six')
    check([len(bset & x) for x in leaves] == [3, 3, 3], 'central half-overlaps')
    check(intersections == [1, 1, 1], 'formal node intersections are one')
    check(len(set.union(*leaves)) == 15, 'formal leaf union exhausts degree15')
    check(-2 * 7 + 3 * 4 + sum(intersections) == 1, 'formal Euler survivor')

    # Independent bounded finite-set check; repetitions among leaves are retained.
    set_rows = 0
    for d in range(2, 8):
        omega = set(range(d))
        for t in range(2, d):
            bset = set(range(t))
            family = [set(x) for x in combinations(range(d), t)
                      if len(set(x) & bset) >= (t + 1) // 2]
            for ix in combinations_with_replacement(range(len(family)), 3):
                ls = [family[i] for i in ix]
                total = sum(len(ls[i] & ls[j]) for i, j in combinations(range(3), 2))
                restricted = [bset & x for x in ls]
                total_b = sum(len(restricted[i] & restricted[j])
                              for i, j in combinations(range(3), 2))
                check(total == 3 * t - len(set.union(*ls)) + len(set.intersection(*ls)),
                      'three-set exact inclusion-exclusion')
                check(total >= 3 * t - d, 'ambient support union bound')
                check(total >= total_b >= 3 * ((t + 1) // 2) - t,
                      'central support incidence bound')
                check(set.union(*ls) <= omega, 'finite-set ambient universe')
                set_rows += 1

    # Complete conjugacy representatives sigma and all same-type tau, D=2,...,6.
    # This bank checks mixed cycles as well as single cycles; it is not the proof.
    pair_rows = braid_rows = retained_rows = commuting_rows = 0
    for d in range(2, 7):
        omega = set(range(d))
        bank = list(permutations(range(d)))
        by_type = {}
        for p in bank:
            by_type.setdefault(cycle_type(p), []).append(p)
        check(set(by_type) == set(partitions(d)), 'complete cycle-type partition universe')
        for parts in partitions(d):
            sigma = canonical(parts)
            t = len(support(sigma))
            check(cycle_type(sigma) == parts, 'canonical permutation type')
            for tau in by_type[parts]:
                pair_rows += 1
                common = support(sigma) & support(tau)
                if mul(sigma, tau) == mul(tau, sigma):
                    commuting_rows += 1
                    check(image(sigma, common) == common, 'commuting invariant intersection')
                    check(len(common) != 1, 'commuting supports cannot meet in one')
                if not braid(sigma, tau):
                    continue
                braid_rows += 1
                g = mul(sigma, tau)
                check(conj(g, sigma) == tau, 'braid conjugator')
                check(2 * len(common) >= t, 'arbitrary-cycle half-support overlap')
                fs, ft = fixed(sigma), fixed(tau)
                joint = fs & ft
                check(image(g, fs) == ft, 'full-fixed re-access')
                for aa in subsets(fs):
                    bb = image(g, aa)
                    check(bb <= ft, 'retained re-access type')
                    k, n = len(aa), len(aa & bb)
                    uu = aa - bb
                    outside = omega - (aa | bb)
                    check(image(tau, uu) <= outside, 'literal cusp injection')
                    check(2 * n >= 3 * k - d, 'actual cusp numerical injection')
                    check(n >= k - len(fs) + len(joint), 'retention-deficit bound')
                    check(n >= k - t // 2, 'actual arbitrary-cycle cusp bound')
                    # Inside-pair H^- keeps n and the deleted overlap, without
                    # asserting a new re-access equation for arbitrary subsets.
                    a_new, b_new = bb, image(inv(tau), aa)
                    check(len(a_new & b_new) == n, 'inside-pair retained count')
                    check(len((omega-a_new) & (omega-b_new)) ==
                          len((omega-aa) & (omega-bb)), 'inside-pair deleted count')
                    retained_rows += 1

    # Marked transitive S4 hostile. It satisfies the Euler ledger and survives
    # the support inequalities; the geometric degree-four supplier is essential.
    a = canonical((1, 1, 2))
    d = a
    b = (0, 2, 1, 3)
    c = (1, 0, 2, 3)
    e = conj(b, c)
    f = conj(inv(e), a)
    check(all(braid(b, leaf) for leaf in (a, c, d)), 'S4 D4 braid edges')
    check(all(mul(x, y) == mul(y, x) for x, y in ((a,c), (a,d), (c,d))),
          'S4 D4 commuting leaves')
    check((conj(e, f), conj(e, b)) == (a, c), 'marked node0 to leaf pair')
    check((conj(c, a), conj(c, e)) == (a, b), 'marked third cusp to central edge')
    node = [len(support(x) & support(y)) for x, y in ((f,b), (a,d), (c,d))]
    cusp = [len(fixed(x) & fixed(y)) for x, y in ((b,c), (b,d), (a,e))]
    check((cusp, node) == ([1,1,1], [0,2,0]), 'marked S4 count control')
    check(-2 * 2 + sum(cusp) + sum(node) == 1, 'S4 exact Euler value one')
    orbit = {0}
    while True:
        larger = orbit | set.union(*(image(p, orbit) for p in (a,b,c,d)))
        if larger == orbit:
            break
        orbit = larger
    check(orbit == set(range(4)), 'S4 transitivity')

    report = {'integer_rows': integer_rows, 'set_rows': set_rows,
              'pair_rows': pair_rows, 'braid_rows': braid_rows,
              'retained_rows': retained_rows, 'commuting_rows': commuting_rows,
              'S4_cusps': cusp, 'S4_nodes': node, 'formal_D15_intersections': intersections}
    raw = json.dumps(report, sort_keys=True, separators=(',', ':')).encode()
    print('DEGREE-UNIFORM PROOF CONTROLS: PASS')
    print('No ambient finite census is used to infer the all-degree theorem.')
    print('Formal comparison identities and parity reduction: PASS')
    print('Complete post-reduction integer rows [t,delta,k,D,Wmin,Wmax]:')
    for row in integer_rows:
        print(' ', row)
    print('Formal D15 set/Euler survivor: node intersections [1,1,1]; commuting fails.')
    print('Bounded set configurations:', set_rows)
    print('Canonical same-type permutation pairs:', pair_rows)
    print('Braided pairs / retained subsets:', braid_rows, '/', retained_rows)
    print('Commuting pairs:', commuting_rows)
    print('Marked transitive S4 hostile: cusps', cusp, 'nodes', node, 'Euler=1')
    print('Semantic SHA256:', sha256(raw).hexdigest())
    print('Always-active gates:', GATES)


if __name__ == '__main__':
    main()
