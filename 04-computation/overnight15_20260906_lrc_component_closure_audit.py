"""Independent positive-box, exact-grid, and actual component controls."""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd
import sys

Q = 91**6
gates = 0


def check(ok, message):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(message)


def distance(x):
    x %= 1
    return min(x, 1-x)


def eligible_sum(s):
    if s > 356:
        return False
    d = 2
    while d*d <= s:
        power = 0
        while s % d == 0:
            power += 1
            s //= d
        if power and (d % 3 != 2 or power > 2):
            return False
        d += 1
    return s == 1 or s % 3 == 2


def component_sizes(row):
    parts = [{i} for i in range(len(row))]
    for i, j in combinations(range(len(row)), 2):
        d = gcd(row[i], row[j])
        if eligible_sum((row[i]+row[j])//d):
            first = next(s for s in parts if i in s)
            second = next(s for s in parts if j in s)
            if first is not second:
                first.update(second)
                parts.remove(second)
    return sorted(map(len, parts))


def positive_box_controls():
    count = 0
    for q in range(2, 15):
        for k in range(2, q+1):
            for x in range(1, q*(k+1)+1):
                s = min(x//k, q)
                r = x-s*k
                check(0 <= r <= q and 0 <= s <= q and r+s*k == x, 'positive bounded division')
                count += 1
    physical = 0
    for q in range(2, 9):
        for k in range(2, q+1):
            for g in range(1, q+1):
                for t in range(1, 13):
                    if gcd(g, t) != 1:
                        continue
                    for v in range(1, 9):
                        delta = gcd(g, v)
                        check(delta == gcd(g, t*v), 'coprime scale clearing')
                        c, x = g//delta, t*v//delta
                        if x <= q*(k+1):
                            s = min(x//k, q)
                            r = x-s*k
                            check(c*t*v-r*g-s*g*k == 0 and max(c, r, s) <= q,
                                  'literal physical crossing row')
                            physical += 1
    print('positive unit-pair boxes', count, 'literal physical crossings', physical)


def grid_controls():
    count = 0
    for t in range(1, 30):
        for g in range(1, t+1):
            if gcd(t, g) != 1:
                continue
            for eta in (F(1, 4), F(3, 7)):
                points = [(g*(eta+j)/t) % 1 for j in range(t)]
                check(len(set(points)) == t, 'coprime physical phases give complete grid')
                for center in (F(1, 11), F(1, 2), F(13, 17)):
                    nearest = min(distance(y-center) for y in points)
                    check(nearest <= F(1, 2*t), 'closed spacing endpoint')
                    check(nearest < F(1, 2*t)+F(1, 2*t*(t+1)), 'strict interior with slack')
                    count += 1
    # Without coprimality the claimed t-grid can collapse.
    collapsed = {(2*(F(1, 4)+j)/4) % 1 for j in range(4)}
    check(len(collapsed) == 2, 'noncoprime clock hostile')
    print('exact translated grid controls', count, 'noncoprime hostile retained')


def threshold_controls():
    for a in range(1, 7):
        b = 13-a
        threshold = Q*a//(7*(b+1))
        tree_bound = 355**(a-1)
        check((tree_bound <= threshold) == (a <= 5), 'automatic spanning-tree threshold boundary')
        check(F(2, 1)*(F(1, b+1)-F(1, 14)) == F(a, 7*(b+1)),
              'full normalized safe-arc length')
        print('split', a, b, 'threshold', threshold, 'tree-bound', tree_bound)
    check(6*Q//56 == 3*Q//28, 'balanced threshold')


def actual_rows():
    count = 0
    for a in range(1, 7):
        b = 13-a
        u = [3**j for j in range(b)]
        k = max(u)
        families = [[3**j for j in range(a)]]
        if a >= 2:
            families.append(sorted([2*3**j for j in range((a+1)//2)]
                                   + [3*3**j for j in range(a//2)]))
        for v in families:
            t, g = 2*Q*k+1, 1
            row = [t*x for x in v]+u
            check(len(row) == len(set(row)) == 13 and reduce(gcd, row) == 1
                  and sum(row) <= Q*Q, 'actual primitive physical-box row')
            check(component_sizes(row) == [a, b], 'actual decoder component sizes')
            check(k <= Q and min(v) <= Q*a//(7*(b+1)), 'unit-core and small-minimum gates')
            check(t > 2*Q*k, 'independent dominance proves decoder equality')
            check(reduce(gcd, v) == 1 and max(v) <= 355**(a-1), 'primitive small tree size')
            # eta=3/4 is V-safe; t=3 mod4 gives the literal lift x=1/4.
            eta = F(3, 4)
            j = (t-3)//4
            check(F(eta+j, t) == F(1, 4), 'physical lift with the retained small-core phase')
            check(min(distance(x*eta) for x in v) > F(1, 14), 'explicit chosen small-core phase is strictly safe')
            check(min(distance(F(x, 4)) for x in row) == F(1, 4), 'literal full-row strict safe phase')
            count += 1
    print('actual equality and safe-phase rows', count, 'including nonunit small cores')


def declared_producer_controls():
    all_u = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22, 24)
    small = [((1,), F(1, 2)), ((2, 3), F(1, 5)), ((2, 3, 6), F(1, 4)),
             ((2, 3, 4, 6), F(1, 5)), ((2, 3, 4, 6, 9), F(1, 5)),
             ((2, 3, 4, 6, 9, 12), F(1, 5))]
    centers = {7: F(4, 9), 8: F(1, 9), 9: F(5, 11), 10: F(1, 11),
               11: F(9, 19), 12: F(1, 13)}
    star = (1013861907, 1036929995, 1060507875, 1084612635, 1096868145)
    p = 179*181*183*185
    check(set(star) == {p} | {(356-q)*(p//q) for q in (179, 181, 183, 185)},
          'billion-minimum star independently reconstructed')
    check(min(star) > 10**9 and reduce(gcd, star) == 1, 'primitive nonunit stress')
    cases = small+[(star, F(1, 2))]
    for v, eta in cases:
        a, b = len(v), 13-len(v)
        u, zeta = all_u[:b], centers[b]
        k, t = max(u), 2*Q*max(u)+1
        row = [t*x for x in v]+list(u)
        check(len(row) == len(set(row)) == 13 and reduce(gcd, row) == 1
              and sum(row) <= Q*Q, 'producer control physical domain')
        check(component_sizes(row) == [a, b] and t > 2*Q*k, 'producer actual entry via graph and dominance')
        check(min(distance(x*eta) for x in v) >= F(1, a+1)
              and min(distance(x*zeta) for x in u) >= F(1, b+1), 'both literal supplier phases')
        # Choose the nearest grid point to the center, unlike the producer's left endpoint.
        target = t*zeta-eta+F(1, 2)
        j = (target.numerator//target.denominator) % t
        phase = (eta+j)/t
        check(min(distance(x*phase) for x in row) > F(1, 14), 'independent nearest-grid full physical witness')
        check(7*(b+1)*k*min(v) <= a*Q*(k+1), 'native delta=1 refinement')
        if v == star:
            check(sum(row) == 90168220083627413785737, 'billion-star physical sum')
            print('billion-star independently selected phase', phase, 'minV', min(v))
    print('all seven declared producer controls independently reconstructed')


def main():
    sys.stdout.reconfigure(newline='\n')
    print('Independent LRC two-component unit-core closure controls')
    positive_box_controls()
    grid_controls()
    threshold_controls()
    actual_rows()
    declared_producer_controls()
    print('PASS', gates, 'always-active gates')


if __name__ == '__main__':
    main()
