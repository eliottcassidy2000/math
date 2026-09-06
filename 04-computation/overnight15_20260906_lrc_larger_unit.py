"""Exact controls for actual two-component LRC entry with a larger unit shape.

The theorem is conditional on the stated actual equality entry and cited
lower-dimensional LRC. Finite controls do not establish universal LRC(14).
No imports from mathematical producers; all gates survive python -O.
"""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, prod
import sys

sys.stdout.reconfigure(newline="\n")
Q = 91**6
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(label)


def norm(x):
    return min(x % 1, (-x) % 1)


def ceildiv(x, y):
    return -((-x)//y)


def atlas_sum(s):
    if not 3 <= s <= 356:
        return False
    p = 2
    while p*p <= s:
        e = 0
        while s % p == 0:
            s //= p
            e += 1
        if e and (p % 3 != 2 or e > 2):
            return False
        p += 1
    return s == 1 or s % 3 == 2


def graph(speeds):
    adj = [set() for _ in speeds]
    for i, j in combinations(range(len(speeds)), 2):
        d = gcd(speeds[i], speeds[j])
        if atlas_sum((speeds[i]+speeds[j])//d):
            adj[i].add(j)
            adj[j].add(i)
    unseen = set(range(len(speeds)))
    components = []
    while unseen:
        seen = {min(unseen)}
        stack = list(seen)
        for i in stack:
            for j in sorted(adj[i]-seen):
                seen.add(j)
                stack.append(j)
        unseen -= seen
        components.append(sorted(seen))
    return sorted(components)


def crossing(A, B, Y, bound):
    """Exact signed-box membership after minimal distinguished coefficient."""
    D = gcd(A, B)
    a, b = sorted((A//D, B//D))
    need(1 <= a < b <= bound, 'internal pair height for full mixed-support test')
    delta = gcd(D, Y)
    c, x = D//delta, Y//delta
    r = x*pow(a, -1, b) % b
    s = (x-a*r)//b
    lower = max(ceildiv(-bound-r, b), ceildiv(s-bound, a))
    upper = min((bound-r)//b, (s+bound)//a)
    return c <= bound and lower <= upper


def phase_lift(V, U, t, g, eta, zeta):
    a, b, L = len(V), len(U), max(U)
    need(all(norm(v*eta) >= F(1, a+1) for v in V), 'literal smaller-core supplier phase')
    need(all(norm(u*zeta) >= F(1, b+1) for u in U), 'literal larger-core supplier phase')
    radius = F(a, 14*(b+1)*L)
    left = t*(zeta-radius)-g*eta
    k = left.numerator//left.denominator+1
    j = k*pow(g, -1, t) % t
    x = (eta+j)/t
    need(0 < k-left < 2*t*radius, 'selected grid point strictly inside safe arc')
    need(all(norm(t*v*x) > F(1, 14) for v in V), 'strict smaller physical clearances')
    need(all(norm(g*u*x) > F(1, 14) for u in U), 'strict larger physical clearances')
    return x


def actual_control(name, V, U, eta, zeta):
    V, U = tuple(sorted(V)), tuple(sorted(U))
    a, b = len(V), len(U)
    L, v = max(U), min(V)
    t, g = 2*Q*L+1, 1
    speeds = tuple(t*w for w in V)+U
    need(a+b == 13 and a <= b and U[0] == 1, 'native split and larger primitive unit')
    need(reduce(gcd, V) == reduce(gcd, U) == gcd(t, g) == 1, 'primitive shapes and coprime scales')
    need(len(set(speeds)) == 13 and sum(speeds) <= Q**2, 'distinct physical box control')
    need(graph(speeds) == [list(range(a)), list(range(a, 13))], 'actual atlas components exactly declared')
    need(max(V) <= 355**(a-1) and L <= Q, 'actual connected and internal height bounds')
    need(7*(b+1)*v <= a*Q, 'native minimum threshold')
    need(t > 2*Q*L, 'literal all-crossing dominance gives W=V_dec')
    checked = 0
    for P, O in ((speeds[:a], speeds[a:]), (speeds[a:], speeds[:a])):
        for A, B in combinations(P, 2):
            for Y in O:
                need(not crossing(A, B, Y, Q), 'every mixed support independently absent')
                checked += 1
    need(checked == 11*a*b//2, 'all mixed triple supports in both orientations')
    delta = gcd(g, v)
    need(t*v > delta*Q*(L+1), 'forced positive-box inequality in actual control')
    need(a*t > 7*(b+1)*L, 'strict cyclic-grid threshold')
    x = phase_lift(V, U, t, g, eta, zeta)
    print(f'ACTUAL {name} a={a} b={b} minV={v} maxV={max(V)} L={L} t={t} g={g} mixed={checked}')
    print(f'  V={V} U={U} sum={sum(speeds)} phase={x}')


def main():
    need(Q == 567869252041, 'native Q')
    print('Q', Q, 'BOX', Q**2)
    for a in range(1, 7):
        b = 13-a
        cap = a*Q//(7*(b+1))
        automatic = 7*(b+1)*355**(a-1) <= a*Q
        need(automatic == (a <= 5), 'exact automatic split range')
        print(f'THRESHOLD a={a} b={b} minV<={cap} connectedMax={355**(a-1)} automatic={automatic}')
    need(3*Q//28 == 60843134147, 'balanced native cutoff')
    need(355**4 == 15882300625, 'fifth connected product')

    boxes = points = 0
    for B in range(2, 22):
        for L in range(2, B+1):
            signed = {r+s*L for r in range(-B, B+1) for s in range(-B, B+1)}
            need(signed == set(range(-B*(L+1), B*(L+1)+1)), 'literal complete signed interval')
            boxes += 1
            for x in range(B*(L+1)+1):
                s = min(x//L, B)
                r = x-s*L
                need(0 <= r <= B and 0 <= s <= B and r+s*L == x, 'positive-box explicit witness')
                points += 1
    need(3 not in {r+6*s for r in range(-2, 3) for s in range(-2, 3)}, 'unit pair outside height scope has a hole')
    print('BOX_CONTROLS', boxes, points, 'HOSTILE B=2 L=6 missing=3')

    scale_cases = 0
    for g in range(1, 31):
        for t in range(1, 31):
            if gcd(t, g) != 1:
                continue
            for v in range(1, 16):
                delta = gcd(g, t*v)
                c = g//delta
                need(delta == gcd(g, v), 'coprime normalization of distinguished coordinate')
                need(c*t*v == g*(t*v//delta), 'literal physical coefficient normalization')
                scale_cases += 1
    print('SCALE_CONTROLS', scale_cases)

    grids = 0
    for t in range(1, 24):
        for g in range(1, 24):
            if gcd(t, g) != 1:
                continue
            eta = F(2, 7)
            need(len({(g*(eta+j)/t) % 1 for j in range(t)}) == t, 'complete shifted grid')
            for h in range(11):
                start, length = F(h, 11), F(8, 7*t)
                left = t*start-g*eta
                k = left.numerator//left.denominator+1
                j = k*pow(g, -1, t) % t
                need((g*j-k) % t == 0 and 0 < k-left < t*length, 'strict arc hit by exact modular inverse')
                grids += 1
    need(all(not F(1, 10) <= F(2*j, 4) % 1 <= F(2, 5) for j in range(4)), 'noncoprime grid hostile')
    print('GRID_CONTROLS', grids, 'HOSTILE t=4 g=2 misses[1/10,2/5]')

    Uall = (1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22, 24)
    small = {1: ((1,), F(1, 2)), 2: ((2, 3), F(1, 5)),
             3: ((2, 3, 6), F(1, 4)), 4: ((2, 3, 4, 6), F(1, 5)),
             5: ((2, 3, 4, 6, 9), F(1, 5)), 6: ((2, 3, 4, 6, 9, 12), F(1, 5))}
    phases = {7: F(4, 9), 8: F(1, 9), 9: F(5, 11), 10: F(1, 11), 11: F(9, 19), 12: F(1, 13)}
    for a, (V, eta) in small.items():
        b = 13-a
        actual_control(f'split_{a}_{b}', V, Uall[:b], eta, phases[b])
    qs = (179, 181, 183, 185)
    P = prod(qs)
    star = (P,)+tuple((356-q)*(P//q) for q in qs)
    need(reduce(gcd, star) == 1 and min(star) > 10**9, 'nonunit primitive fifth-star stress')
    actual_control('billion_minimum_star', star, Uall[:8], F(1, 2), phases[8])
    print('ACTUAL_CONTROLS 7; all six splits, larger scale g=1, singleton and unitless fifth-star retained')
    print('PASS', GATES, 'always-active gates')


if __name__ == '__main__':
    main()
