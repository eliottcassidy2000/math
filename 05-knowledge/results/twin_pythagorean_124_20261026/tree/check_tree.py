#!/usr/bin/env python3
"""Exact finite controls for RESULT.md; standard library, no producer imports."""
from collections import Counter
from itertools import combinations
from math import gcd, isqrt
import json

CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def prime(n):
    if n < 2:
        return False
    return all(n % d for d in range(2, isqrt(n) + 1))


def factors(n):
    out = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def phi(m, n):
    return m * m - n * n, 2 * m * n, m * m + n * n


def xi(x, y):
    return x * y, (y * y - x * x) // 2, (y * y + x * x) // 2


def ppt(t):
    a, b, c = t
    return (min(t) > 0 and a % 2 == 1 and b % 2 == 0
            and a * a + b * b == c * c and gcd(gcd(a, b), c) == 1)


def P(t):
    return t, (t * t - 1) // 2, (t * t + 1) // 2


def Q(k):
    m = 2 * k + 2
    return m * m - 1, 2 * m, m * m + 1


def F(p):
    return p * (p + 2), 2 * (p + 1), (p + 1) ** 2 + 1


MATRICES = {
    'U': ((1, -2, 2), (2, -1, 2), (2, -2, 3)),
    'A': ((1, 2, 2), (2, 1, 2), (2, 2, 3)),
    'D': ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3)),
}


def apply(label, t):
    return tuple(sum(v * w for v, w in zip(row, t)) for row in MATRICES[label])


def child(label, m, n):
    return {'U': (2 * m - n, m), 'A': (2 * m + n, m),
            'D': (m + 2 * n, n)}[label]


def word(m, n):
    letters = []
    while (m, n) != (2, 1):
        old = m + n
        if m < 2 * n:
            m, n = n, 2 * n - m
            letters.append('U')
        elif m < 3 * n:
            m, n = n, m - 2 * n
            letters.append('A')
        else:
            m, n = m - 2 * n, n
            letters.append('D')
        require(m > n > 0 and m + n < old, 'strict inverse descent')
    return ''.join(reversed(letters))


def aliquot(n):
    if n == 1:
        return 0
    total = 1
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            total += d
            if d * d != n:
                total += n // d
    return total


def main():
    root = (3, 4, 5)
    euclid = {}
    for m in range(2, 121):
        for n in range(1, m):
            if gcd(m, n) != 1 or (m + n) % 2 != 1:
                continue
            t = phi(m, n)
            a, b, c = t
            require(ppt(t), 'Euclid primitive gate')
            require(t not in euclid, 'Euclid uniqueness')
            euclid[t] = (m, n)
            x, y = isqrt(c - b), isqrt(c + b)
            require(x * x == c - b and y * y == c + b, 'two exact roots')
            require((x, y) == (m - n, m + n) and xi(x, y) == t, 'chart inverse')
            require(c - a == (y - x) ** 2 // 2, 'gap height')
            require((c - a == 2) == (n == 1), 'height-two ray criterion')
            if prime(a):
                require(m - n == 1 and t == P(a), 'global prime-leg rigidity')
            twin_image = c - a == 2 and len(factors(a)) == 2
            require(twin_image == (y - x == 2 and prime(x) and prime(y)),
                    'twin image iff')
            address = word(m, n)
            replay = root
            for label in address:
                replay = apply(label, replay)
            require(replay == t, 'descent address matrix replay')
            for label in MATRICES:
                u = apply(label, t)
                require(u == phi(*child(label, m, n)), 'symmetric-square branch')
                require(ppt(u) and u[2] > c, 'positive ascending branch')
                if prime(a) and prime(u[0]):
                    require(label == 'U' and u[0] == a + 2,
                            'all prime-to-prime branches')

    # Independent full ternary prefix, including every word of length <=8.
    tree = {'': (root, (2, 1))}
    level = ['']
    seen = {root}
    for depth in range(1, 9):
        nxt = []
        for w in level:
            t, (m, n) = tree[w]
            for label in MATRICES:
                ww = w + label
                tt = apply(label, t)
                mn = child(label, m, n)
                require(tt not in seen, 'tree prefix has no duplicate node')
                require(tt == phi(*mn), 'tree prefix against Euclid')
                require(word(*mn) == ww, 'tree prefix exact address inverse')
                seen.add(tt)
                tree[ww] = (tt, mn)
                nxt.append(ww)
        level = nxt
    require(len(tree) == (3 ** 9 - 1) // 2, 'full ternary universe size')
    for t, mn in tree.values():
        if mn[0] <= 120:
            require(euclid.get(t) == mn, 'tree/Euclid census intersection')

    chart_count = 0
    for y in range(3, 200, 2):
        for x in range(1, y, 2):
            if gcd(x, y) != 1:
                continue
            chart_count += 1
            t = xi(x, y)
            require(ppt(t) and t == phi((x + y) // 2, (y - x) // 2),
                    'odd-root independent universe')

    # Odd-leg fibres directly from factor pairs, checked against Euclid solves.
    semiprime_fibres = 0
    for a in range(3, 2000, 2):
        splits = [(x, a // x) for x in range(1, isqrt(a) + 1)
                  if a % x == 0 and x < a // x and gcd(x, a // x) == 1]
        triples = {xi(x, y) for x, y in splits}
        direct = set()
        for n in range(1, (a - 1) // 2 + 1):
            m = isqrt(a + n * n)
            if m * m == a + n * n and gcd(m, n) == 1 and (m + n) % 2:
                direct.add(phi(m, n))
        require(triples == direct, 'factor split vs independent difference-square solve')
        fs = factors(a)
        if len(fs) == 2 and fs[0] != fs[1]:
            semiprime_fibres += 1
            require(triples == {P(a), xi(*fs)}, 'exact two semiprime nodes')
            require((xi(*fs)[2] - a == 2) == (fs[1] - fs[0] == 2),
                    'semiprime minimal-height iff')

    # Prime sieve independent of the trial primality function.
    cap = 20001
    sieve = bytearray(b'\x01') * (cap + 1)
    sieve[:2] = b'\x00\x00'
    for d in range(2, isqrt(cap) + 1):
        if sieve[d]:
            sieve[d * d:cap + 1:d] = b'\x00' * ((cap - d * d) // d + 1)
    primes = [p for p in range(3, cap + 1, 2) if sieve[p]]
    for p in range(3, cap + 1, 2):
        require(prime(p) == bool(sieve[p]), 'trial primality vs sieve')
    twins = [p for p in primes if p + 2 <= cap and sieve[p + 2]]
    H_edges = [(p, p + 2) for p in twins]
    J_vertices = {Q((p - 1) // 2): p for p in twins}
    require(len(J_vertices) == len(H_edges), 'edge-to-node bijection cardinality')
    for p in twins:
        t = F(p)
        require(t == Q((p - 1) // 2) == xi(p, p + 2), 'three twin constructions')
        require(apply('U', P(p)) == P(p + 2), 'actual prime U edge')
        require(tuple(y - x for x, y in zip(P(p), P(p + 2)))
                == (2, t[1], t[1]), 'endpoint difference recovers even leg')
        require(apply('D', t) == F(p + 2), 'unfiltered edge-shift intertwining')
        if p < 500:
            a, b, c = t
            require(aliquot(c - b) == b // 2, 'lower square aliquot')
            require(aliquot(c + b) == b // 2 + 2, 'upper square aliquot')
            require(aliquot(a) == b + 1, 'semiprime aliquot')
    line_edges = {tuple(sorted((e[0], f[0]))) for e, f in combinations(H_edges, 2)
                  if set(e) & set(f)}
    d_edges = {tuple(sorted((p, J_vertices[apply('D', t)])))
               for t, p in J_vertices.items() if apply('D', t) in J_vertices}
    require(line_edges == d_edges == {(3, 5)}, 'actual graph isomorphism and unique edge')
    degrees = Counter()
    for p, q in H_edges:
        degrees[p] += 1
        degrees[q] += 1
    require([p for p in primes if degrees[p] == 2] == [5], 'only degree-two prime vertex')
    require(max(degrees.values()) == 2, 'prime induced degree bound')

    require(Q(0) == root and not (prime(1) and prime(3)), 'root excluded')
    require(Q(3) == (63, 16, 65), 'minimal distinct-factor hostile')
    require(len(set(factors(63))) == 2 and len(factors(63)) == 3,
            'omega is not Omega')
    require(xi(3, 7) == (21, 20, 29) and xi(3, 7)[2] - 21 == 8,
            'semiprime does not force twins')
    require(phi(3, 1) == (8, 6, 10) and gcd(gcd(8, 6), 10) == 2,
            'half-sized Euclid map parity hostile')
    require(apply('U', P(7)) == P(9) and not prime(9), 'primality not invariant')
    require(F(5) == (35, 12, 37) and aliquot(25) == aliquot(6) == 6,
            'unique-six positive anchor')
    require(F(3) == (15, 8, 17) and aliquot(9) == 4 and aliquot(4) == 3,
            'preceding center four hostile')
    require(F(3)[2] + F(3)[1] == F(5)[2] - F(5)[1] == 25,
            'unique edge shares perfect-image square with opposite flags')

    print(json.dumps({
        'status': 'PASS', 'scope': 'FINITE-EXACT; global conclusions use RESULT.md proofs',
        'checks': CHECKS,
        'universes': {'Euclid_m_at_most': 120, 'Euclid_nodes': len(euclid),
                      'ternary_depth_at_most': 8, 'ternary_nodes': len(tree),
                      'coprime_odd_root_y_at_most': 199, 'odd_root_pairs': chart_count,
                      'odd_leg_fibre_a_at_most': 1999, 'distinct_semiprime_fibres': semiprime_fibres,
                      'prime_cap': cap, 'prime_odd_vertices': len(primes),
                      'twin_edges': len(twins)},
        'restricted_graph': {'shared_edge_pairs': sorted(line_edges),
                             'D_graph_edges': [[list(F(3)), list(F(5))]],
                             'isolated_D_vertices_in_universe': len(twins) - 2},
        'first_twin_images': [{'pair': [p, p + 2], 'triple': F(p),
                              'U_edge_depths': [(p - 3) // 2, (p - 1) // 2],
                              'D_depth': (p - 1) // 2} for p in twins[:10]],
        'hostiles': {'root_not_twin': root, 'omega_not_Omega': Q(3),
                     'non_twin_semiprime': xi(3, 7), 'failed_prime_continuation': P(9),
                     'half_sized_nonprimitive': phi(3, 1)}
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
