#!/usr/bin/env python3
"""Exact signed-PPT bridge census and universal elliptic-extension controls.

Standard library only. No assertion is suppressed by Python -O.
Universe: primitive positive parity-oriented Pythagorean triples with c<=10000
and c*c>4*a*b. All four signed contacts between distinct petals are tested.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import gcd, isqrt, lcm

CHECKS = 0

def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def square_root(n):
    q = isqrt(n)
    return q if q*q == n else None


def rational_root(x):
    if x < 0:
        return None
    a, b = square_root(x.numerator), square_root(x.denominator)
    return None if a is None or b is None else F(a, b)


def spread(p):
    a, b, c = p[:3]
    return F(2*a*b, c*c)


def primitive_petals(bound):
    out = []
    for n in range(2, isqrt(bound)+1):
        for m in range(1, n):
            a, b, c = n*n-m*m, 2*n*m, n*n+m*m
            if c > bound or gcd(m, n) != 1 or (m+n) % 2 == 0:
                continue
            if c*c > 4*a*b:
                out.append((a, b, c, m, n))
    return sorted(out, key=lambda p: (p[2], p[0], p[1]))


def contact_edges(petals):
    edges = []
    for i, (a, b, c, *_) in enumerate(petals):
        for j in range(i+1, len(petals)):
            d, e, f, *_ = petals[j]
            for s in (-1, 1):
                for t in (-1, 1):
                    numerator = c*c*f*f + 2*s*a*b*f*f + 2*t*d*e*c*c
                    q = square_root(numerator)
                    if q is not None:
                        root = F(q, c*f)
                        check(root*root == 1+s*spread(petals[i])+t*spread(petals[j]), 'contact equation')
                        edges.append((i, j, s, t, root))
    return edges


def components_and_chains(petals, edges):
    adj = defaultdict(list)
    for i, j, s, t, _ in edges:
        adj[i].append((j, s, t))
        adj[j].append((i, t, s))
    seen, components = set(), []
    for u in range(len(petals)):
        if u in seen:
            continue
        seen.add(u)
        comp = [u]
        for v in comp:
            for w, _, _ in adj[v]:
                if w not in seen:
                    seen.add(w)
                    comp.append(w)
        components.append(comp)
    best = []
    def visit(path, used, required_port):
        nonlocal best
        if len(path) > len(best):
            best = list(path)
        u = path[-1][0]
        for v, s, t in adj[u]:
            if v in used or (required_port is not None and s != required_port):
                continue
            visit(path+[(v, t)], used | {v}, -t)
    for comp in components:
        if len(comp) > 1:
            for u in comp:
                visit([(u, None)], {u}, None)
    return components, best


def cycle_certificate(petals, ordered_ports):
    H = lcm(*(p[2] for p in petals))
    middle = 2*H*H
    ports = {}
    for i, (a, b, c) in enumerate(petals):
        g = H//c
        ports[i, -1] = 2*g*g*(c*c-4*a*b)
        ports[i, 1] = 2*g*g*(c*c+4*a*b)
        check(ports[i, -1]+ports[i, 1] == 2*middle, 'AP pairing')
        for x, y in [(ports[i, -1], middle), (middle, ports[i, 1]), (ports[i, -1], ports[i, 1])]:
            check(square_root(x+y) is not None, 'petal edge')
    cycle = [middle] + [ports[item] for item in ordered_ports]
    check(len(set(cycle)) == len(cycle) and min(cycle) > 0, 'cycle distinct positivity')
    roots = []
    for i in range(len(cycle)):
        q = square_root(cycle[i]+cycle[(i+1) % len(cycle)])
        check(q is not None, 'cycle edge')
        roots.append(q)
    return H, cycle, roots


def on_curve(P, A, B):
    if P is None:
        return True
    X, Y = P
    return Y*Y == (X-A*A)*(X*X-B**4)


def add(P, Q, A, B):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and y == -v:
        return None
    slope = ((3*x*x-2*A*A*x-B**4)/(2*y) if P == Q else (v-y)/(u-x))
    xx = slope*slope+A*A-x-u
    return xx, -y+slope*(x-xx)


def val2(x):
    if not x:
        raise ValueError('v2(0) not finite')
    n, d, out = abs(x.numerator), x.denominator, 0
    while n % 2 == 0:
        out += 1
        n //= 2
    while d % 2 == 0:
        out -= 1
        d //= 2
    return out


def quartic_to_curve(u, z, k):
    X = (k*z+k*k+2*u+k*k*u*u)/(2*u*u)
    Y = (X-k*k)*(X*u-1)/k
    return X, Y


def curve_to_quartic(X, Y, k):
    u = (k*Y+X-k*k)/(X*(X-k*k))
    z = (2*u*u*X-k*k-2*u-k*k*u*u)/k
    return u, z


def elliptic_checks(petals):
    for a, b, c, m, n in petals:
        t = F(n, m)
        for sign in (-1, 1):
            A, B = abs(a+sign*b), c
            k, u, z = F(A, B), sign*t, t*t+1
            check(A % 2 == 1 and B % 2 == 1 and gcd(A, B) == 1, 'odd reduced k')
            check(k*k != 1, 'nonsingular elliptic curve')
            check(z*z == k*k*(u*u+1)**2-4*u*(u*u-1), 'anchor quartic')
            X, Y = quartic_to_curve(u, z, k)
            check(Y*Y == (X-k*k)*(X*X-1), 'birational curve equation')
            check(X != 0 and X != k*k, 'anchor inverse chart')
            check(curve_to_quartic(X, Y, k) == (u, z), 'birational round trip')
            P = (F(0), F(A*B*B))
            check(on_curve(P, A, B), 'universal point')
            Q = add(P, P, A, B)
            check(Q[0] == A*A+F(B**4, 4*A*A), 'first doubling x')
            for j in range(1, 5):
                check(on_curve(Q, A, B), 'doubled point on curve')
                check(val2(Q[0]) == -2*j and val2(Q[1]) == -3*j, 'doubling valuation')
                Q = add(Q, Q, A, B)
    print('ELLIPTIC CONTROLS: both signs on all 528 petals; four doublings per sign')
    # A genuine distinct neighbor maps in both directions.
    k, u, z = F(241, 221), F(33, 4), F(941, 16)
    X, Y = quartic_to_curve(u, z, k)
    check(curve_to_quartic(X, Y, k) == (u, z), 'bridge round trip')
    check(Y*Y == (X-k*k)*(X*X-1), 'bridge curve point')
    print('ANCHOR 221 BRIDGE POINT (translated model x=X-k^2):', (X-k*k)*221**2, Y*221**3)



def rank_and_normalization_checks():
    # These test the ordering lemma, not the existence of square-sum contacts.
    for k in range(2, 51):
        d = [10*k+i for i in range(k)]
        M = 100*k
        path = [M]
        expected_rank = [k]
        for i, gap in enumerate(d, 1):
            path.extend([M-gap, M+gap])
            expected_rank.extend([k-i, k+i])
        rank = {v: i for i, v in enumerate(sorted(path))}
        ranks = [rank[v] for v in path]
        differences = [abs(y-x) for x,y in zip(path, path[1:])]
        rank_differences = [abs(y-x) for x,y in zip(ranks, ranks[1:])]
        check(ranks == expected_rank, 'Walecki rank sequence')
        check(rank_differences == list(range(1,2*k+1)), 'graceful rank differences')
        check(all(x<y for x,y in zip(differences, differences[1:])), 'strictly increasing original gaps')
        check(d[-1] not in differences, 'distinct closing original gap')
    # Odd/odd rational parameters are valid after removing the factor two.
    t = F(9,7)
    n,m = t.numerator,t.denominator
    raw = (n*n-m*m,2*n*m,n*n+m*m)
    delta = gcd(gcd(raw[0],raw[1]),raw[2])
    normalized = tuple(x//delta for x in raw)
    check(delta == 2 and normalized == (16,63,65), 'odd/odd primitive normalization')
    check(spread(normalized) == 4*t*(t*t-1)/(t*t+1)**2, 'normalization preserves spread')
    print('RANK CONTROLS: k=2..50; original gaps increase, rank gaps are 1..2k; odd/odd normalization checked')


def main():
    petals = primitive_petals(10000)
    check(len(petals) == 528, 'census size')
    check(len({spread(p) for p in petals}) == len(petals), 'no repeated primitive spread')
    edges = contact_edges(petals)
    check(len(edges) == 15, 'contact count')
    # Independent fraction-square predicate checks all pairs and all ports.
    direct = []
    for i in range(len(petals)):
        for j in range(i+1, len(petals)):
            for s in (-1, 1):
                for t in (-1, 1):
                    r = rational_root(1+s*spread(petals[i])+t*spread(petals[j]))
                    if r is not None:
                        direct.append((i, j, s, t, r))
    check(edges == direct, 'independent rational-square contact test')
    comps, best = components_and_chains(petals, edges)
    histogram = dict(sorted(Counter(map(len, comps)).items()))
    check(histogram == {1: 503, 2: 6, 3: 3, 4: 1}, 'component census')
    check(len(best) == 3, 'maximum compatible chain')
    print('UNIVERSE: primitive admissible PPTs, c<=10000; N=', len(petals))
    print('SIGNED CONTACTS:', len(edges), 'PORT COUNTS:', dict(sorted(Counter((e[2], e[3]) for e in edges).items())))
    print('COMPONENT SIZES:', histogram)
    print('MAXIMUM PORT-COMPATIBLE SIMPLE CHAIN:', len(best))
    for i, j, s, t, r in edges:
        print('EDGE', petals[i][:3], s, petals[j][:3], t, 'ROOT', r)
    print('BEST CHAIN:', [(petals[i][:3], sign) for i, sign in best])
    five = cycle_certificate([(21,220,221), (1073,264,1105)], [(0,-1),(0,1),(1,-1),(1,1)])
    seven = cycle_certificate([(1443,76,1445),(1419,8260,8381),(495,4888,4913)], [(0,1),(0,-1),(1,-1),(1,1),(2,-1),(2,1)])
    for name, (H, cycle, roots) in [('FIVE', five), ('SEVEN', seven)]:
        print(name, 'H=', H, 'CYCLE=', cycle, 'ROOTS=', roots)
    check(five[0] == 1105 and seven[0] == 712385, 'scales')
    # Hostile control: all contacts at the 1105 petal hit the same port.
    hub = next(i for i, p in enumerate(petals) if p[:3] == (1073,264,1105))
    signs = [s if i == hub else t for i,j,s,t,_ in edges if hub in (i,j)]
    check(signs == [-1,-1,-1], 'unsigned-star port obstruction')
    print('HOSTILE CONTROL: 1105 petal has three contacts, all at its minus port')
    elliptic_checks(petals)
    rank_and_normalization_checks()
    print('PASS:', CHECKS, 'explicit checks; rational arithmetic and integer square tests only')

if __name__ == '__main__':
    main()
