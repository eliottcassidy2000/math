#!/usr/bin/env python3
"""Side-based and factor-of-a-square checks, independent of the tree producer."""
from itertools import combinations
from math import gcd, isqrt
import json

gates = 0


def need(x, label):
    global gates
    gates += 1
    if not x:
        raise RuntimeError(label)


def prime(x):
    return x > 1 and all(x % d for d in range(2, isqrt(x) + 1))


U = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
A = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
D = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))


def mv(M, v):
    return tuple(sum(M[i][j] * v[j] for j in range(3)) for i in range(3))


def P(t):
    return t, (t*t-1)//2, (t*t+1)//2


def F(p):
    return p*(p+2), 2*(p+1), (p+1)**2+1


side_universe = []
prime_edges = []
for c in range(5, 2001, 2):
    for a in range(1, c, 2):
        b = isqrt(c*c-a*a)
        if b*b != c*c-a*a or b % 2 or gcd(a, b) != 1:
            continue
        T = (a, b, c)
        side_universe.append(T)
        x, y = isqrt(c-b), isqrt(c+b)
        need(x*x == c-b and y*y == c+b, 'square-root chart on direct sides')
        need(x*y == a and (y-x)**2 == 2*(c-a), 'factor and height chart')
        if prime(a):
            need(T == P(a), 'prime odd leg direct universe')
            kids = [mv(M, T) for M in (U, A, D)]
            need([v[0] for v in kids] == [a+2, a*(2*a+1), a*(2*a-1)],
                 'three odd-leg factors')
            for label, child in zip('UAD', kids):
                if prime(child[0]):
                    need(label == 'U' and child[0] == a+2, 'global prime edge filter')
                    prime_edges.append((a, child[0]))

# Independently solve c^2-b^2=a^2 through factor pairs of a^2.
semiprime_cases = 0
for a in range(3, 1000, 2):
    ps = [p for p in range(3, isqrt(a)+1, 2) if a % p == 0 and prime(p) and prime(a//p) and p != a//p]
    if not ps:
        continue
    p, q = ps[0], a//ps[0]
    semiprime_cases += 1
    actual = set()
    for d in range(1, a, 2):
        if a*a % d:
            continue
        e = a*a//d
        c, b = (e+d)//2, (e-d)//2
        if gcd(a, b) == 1:
            actual.add((a, b, c))
    expected = {P(a), (a, (q*q-p*p)//2, (q*q+p*p)//2)}
    need(actual == expected, 'two primitive nodes via factors of a-square')

twins = [p for p in range(3, 2000, 2) if prime(p) and prime(p+2)]
edges = [(p, p+2) for p in twins]
images = {p:F(p) for p in twins}
for p, q in edges:
    need(mv(U, P(p)) == P(q), 'actual U edge')
    need(mv(D, F(p)) == F(p+2), 'actual D continuation')
    # Repeated matrix application verifies the root offset without Euclid code.
    image = (3, 4, 5)
    for k in range((p-1)//2):
        image = mv(D, image)
    need(image == F(p), 'D word length and root shift')
sharing = []
for e, f in combinations(edges, 2):
    p, q = e[0], f[0]
    shared = bool(set(e) & set(f))
    adjacent = mv(D, images[p]) == images[q] or mv(D, images[q]) == images[p]
    need(shared == adjacent, 'adjacency preserved and reflected')
    if shared:
        sharing.append([p,q])
need(sharing == [[3,5]], 'unique edge in finite target')
need(prime(23) and not prime(21) and not prime(25), 'isolated H vertex disappears')
need(F(3)[2]+F(3)[1] == F(5)[2]-F(5)[1] == 25, 'oppositely flagged shared square')
need(sum(d for d in range(1,25) if 25%d==0) == 6,
     'shared square has perfect image')
need(sum(d for d in range(1,6) if 6%d==0) == 6,
     'perfect-six direct divisors')
need(F(3)[1]//2 == 4 and F(5)[1]//2 == 6, 'distinct centres retained')

print(json.dumps({'status':'PASS','gates':gates,'direct_hypotenuse_max':2000,
                  'direct_triples':len(side_universe),'factor_square_odd_leg_max':999,
                  'semiprime_fibres':semiprime_cases,'twin_lower_bound_exclusive':2000,
                  'twin_pairs':len(twins),'shared_edge_pairs':sharing,
                  'prime_edges_seen_in_direct_universe':prime_edges},indent=2,sort_keys=True))
