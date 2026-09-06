#!/usr/bin/env python3
"""Exact tree-cofactor and actual nonunit decoder-entry controls.

No mathematical producer imports; universal arguments are in the result note.
All gates remain active under python -O.
"""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations, product
from math import gcd, prod

Q = 91**6
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def tree(code, n):
    degrees = [1]*n
    for x in code:
        degrees[x] += 1
    edges = []
    for x in code:
        y = next(i for i in range(n) if degrees[i] == 1)
        edges.append((y, x))
        degrees[y] -= 1
        degrees[x] -= 1
    ends = [i for i in range(n) if degrees[i] == 1]
    return edges+[tuple(ends)]


def cuts(edges, n):
    answer = []
    for forbidden, (u, v) in enumerate(edges):
        seen, queue = {u}, [u]
        for x in queue:
            for k, (i, j) in enumerate(edges):
                if k == forbidden:
                    continue
                y = j if x == i else i if x == j else None
                if y is not None and y not in seen:
                    seen.add(y)
                    queue.append(y)
        answer.append(seen)
    return answer


def cofactors(edge_weights, sides, n):
    values = [1]*n
    for (a, b), side in zip(edge_weights, sides):
        for i in range(n):
            values[i] *= b if i in side else a
    return values


def tree_bank():
    count = 0
    for n, C in ((3, 8), (4, 6), (5, 4)):
        weights = ([(a, b) for a in range(1, C) for b in range(1, C-a+1)]
                   if n <= 4 else [(1,1),(1,2),(2,1),(1,3),(3,1)])
        for code in product(range(n), repeat=n-2):
            edges = tree(code, n)
            sides = cuts(edges, n)
            degrees = [sum(i in e for e in edges) for i in range(n)]
            ell = sum(d == 1 for d in degrees)
            for ws in product(weights, repeat=n-1):
                values = cofactors(ws, sides, n)
                for (u, v), (a, b) in zip(edges, ws):
                    need(a*values[u] == b*values[v], 'literal weighted tree kernel')
                m = min(values)
                need(m*(n-1)**(n-1) <= C**(n-1)*(n-2)**(n-2), 'uniform cofactor bound')
                need(m**ell*ell**(ell*(n-1)) <=
                     C**((n-1)*ell)*(ell-1)**((ell-1)*(n-1)), 'leaf refinement')
                content = reduce(gcd, values)
                need(min(v//content for v in values) <= m, 'primitive content direction')
                count += 1
    for n in range(3, 10):
        edges = [(0, j) for j in range(1, n)]
        values = cofactors([(1,n-2)]*(n-1), cuts(edges,n), n)
        need(min(values) == (n-2)**(n-2), 'sharp real/cofactor star')
        need(values[0] == (n-2)*min(values), 'star center versus minimum')
    need(27 > 2**4, 'half-sum cofactor guess hostile: n5 star')
    return count


def admissible_sum(s):
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


def graph_components(speeds):
    adj = [set() for _ in speeds]
    for i, j in combinations(range(len(speeds)), 2):
        if admissible_sum((speeds[i]+speeds[j])//gcd(speeds[i],speeds[j])):
            adj[i].add(j)
            adj[j].add(i)
    unseen, components = set(range(len(speeds))), []
    while unseen:
        seen, queue = {min(unseen)}, [min(unseen)]
        for x in queue:
            for y in sorted(adj[x]-seen):
                seen.add(y)
                queue.append(y)
        unseen -= seen
        components.append(sorted(seen))
    return components


def ceildiv(x,y):
    return -((-x)//y)


def crossing(P, S, Y):
    common = gcd(P,S)
    a,b = sorted((P//common,S//common))
    need(b <= Q, 'internal pair height')
    delta = gcd(common,Y)
    c,x = common//delta,Y//delta
    if c > Q:
        return False
    r = (x*pow(a,-1,b)) % b
    s = (x-a*r)//b
    low = max(ceildiv(-Q-r,b),ceildiv(s-Q,a))
    high = min((Q-r)//b,(s+Q)//a)
    return low <= high


def norm(x):
    return min(x % 1, (-x) % 1)


def actual_controls():
    larger = (6,10,15,20,24,30,60,18,2,3,40,45)
    small = {1:((1,),F(1,2)),2:((2,3),F(1,5)),3:((2,3,6),F(1,4)),
             4:((2,3,4,6),F(1,5)),5:((2,3,4,6,9),F(1,5)),
             6:((2,3,4,6,9,12),F(1,5))}
    for a,(V,eta) in small.items():
        b = 13-a
        U = tuple(sorted(larger[:b]))
        t,g = 120*Q+1,1
        L,v = max(U),min(V)
        speeds = tuple(t*w for w in V)+U
        need(1 not in U and reduce(gcd,U) == reduce(gcd,V) == 1,'primitive nonunit shapes')
        need(len(set(speeds)) == 13 and sum(speeds) <= Q*Q,'actual distinct physical box')
        need(graph_components(speeds) == [list(range(a)),list(range(a,13))], 'literal actual graph')
        need(t > 2*Q*L, 'dominance excludes every bounded mixed row')
        mixed = 0
        for inside,outside in ((speeds[:a],speeds[a:]),(speeds[a:],speeds[:a])):
            for P,S in combinations(inside,2):
                for Y in outside:
                    need(not crossing(P,S,Y),'every actual mixed support absent')
                    mixed += 1
        need(mixed == 11*a*b//2,'both mixed support orientations')
        u = min((w for w in U if w < L), key=lambda w:gcd(w,L))
        D = gcd(u,L)
        need(D > 1, 'no maximum-coprime shortcut')
        delta = gcd(g*D,t*v)
        A,B = u//D,L//D
        R = Q*(A+B)-(A-1)*(B-1)
        need(g*D//delta <= Q and a*delta*R >= 7*(b+1)*L*v,'native unit-free gate')
        need(all(norm(w*eta) >= F(1,a+1) for w in V),'smaller supplier phase')
        zeta = F(1,7)
        need(all(norm(w*zeta) >= F(1,b+1) for w in U),'larger supplier phase')
        radius = F(a,14*(b+1)*L)
        left = t*(zeta-radius)-eta
        j = left.numerator//left.denominator+1
        x = (eta+j)/t
        need(zeta-radius < x < zeta+radius,'actual selected grid point')
        clearance = min(norm(w*x) for w in speeds)
        need(clearance > F(1,14),'full physical consequence')
        print(f'ACTUAL split={a}+{b} D={D} minV={v} maxU={L} mixed={mixed} t={t}')
        print(f'  V={V} U={U} phase={x} clearance={clearance}')


def main():
    cases = tree_bank()
    print('LEAF_COFACTOR_AND_UNITFREE_ENTRY')
    print(f'tree_weight_cases={cases}; all labelled trees n3,4,5 in declared weight banks')
    expected = (6240321451,76388115,698294,4854,26,0)
    for a in range(1,7):
        b = 13-a
        bound = (F(1) if a == 1 else F(177) if a == 2 else
                 F(356**(a-1)*(a-2)**(a-2),(a-1)**(a-1)))
        integer = bound.numerator//bound.denominator
        scale = {12:1,11:2,10:4,9:9,8:30,7:90}[b]
        cutoff = min(Q//scale,a*Q//(7*(b+1)*integer))
        need(cutoff == expected[a-1],'exact endpoint cutoff')
        print(f'BOUND split={a}+{b} minV<={integer} raw={bound} gcd_cutoff={cutoff}')
    actual_controls()
    print('actual_nonunit_controls=6; all maxima lack a coprime partner')
    print('balanced automatic closure=OPEN; real-cofactor sharpness is not primitive sharpness')
    print(f'PASS gates={GATES}')


if __name__ == '__main__':
    main()
