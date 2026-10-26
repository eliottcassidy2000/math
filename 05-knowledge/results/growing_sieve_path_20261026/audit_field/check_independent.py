#!/usr/bin/env python3
"""Independent field-lane audit: chart tournaments and vectors modulo 105."""
from collections import Counter, deque
from itertools import permutations
from math import gcd, isqrt, prod
import json

GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


P = tuple(range(7)) + (None,)


def vec(t):
    return (1, 0) if t is None else (t, 1)


def det(v, w):
    return (v[0]*w[1]-v[1]*w[0]) % 7


def slope(v):
    return None if v[1] % 7 == 0 else v[0]*pow(v[1] % 7, 5, 7) % 7


def apply(mat, t):
    x, y = vec(t)
    a, b, c, d = mat
    return slope((a*x+b*y, c*x+d*y))


def chart(p, x):
    # Inverse of determinant-one matrix [[p,-1],[1,0]].
    if p is None:
        return x
    return 0 if x is None else pow((p-x) % 7, 5, 7)


def arcs(p, bit):
    vertices = [x for x in P if x != p]
    return frozenset((u, v) for u in vertices for v in vertices if u != v
                     and pow((chart(p, v)-chart(p, u)) % 7, 3, 7) == (1 if bit == 0 else 6))


def integer_word(w):
    m, n = 2, 1
    for letter in w:
        if letter == 'U':
            m, n = 2*m-n, m
        elif letter == 'A':
            m, n = 2*m+n, m
        elif letter == 'D':
            m, n = m+2*n, n
        else:
            raise RuntimeError('letter')
    return m, n, w.count('A') % 2


PRIMES = [p for p in range(2, 500) if all(p % d for d in range(2, isqrt(p)+1))]


def factor(n):
    out = {}
    for p in PRIMES:
        while n % p == 0:
            n //= p
            out[p] = out.get(p, 0)+1
        if p*p > n:
            break
    if n > 1:
        need(all(n % d for d in range(2, isqrt(n)+1)), 'final factor prime')
        out[n] = out.get(n, 0)+1
    return out


def record(w):
    m, n, bit = integer_word(w)
    a, b, c = m*m-n*n, 2*m*n, m*m+n*n
    fs = [factor(x) for x in (a, b, c)]
    support = sorted(set().union(*(f.keys() for f in fs)))
    for x, f in zip((a, b, c), fs):
        need(x == prod(p**e for p, e in f.items()), 'factor reconstruction')
    need(a*a+b*b == c*c and gcd(a, b) == gcd(b, c) == gcd(a, c) == 1, 'primitive abc')
    return {'pair': [m, n], 'triple': [a, b, c], 'bit': bit,
            'support': support, 'radical': prod(support), 'C': c*c,
            'hit': c*c > prod(support)}


def main():
    tournaments = {(p, bit): arcs(p, bit) for p in P for bit in (0, 1)}
    need(len(set(tournaments.values())) == 16, 'distinct actual marked tournaments')
    for (p, bit), relation in tournaments.items():
        need(len(relation) == 21, '21 arcs')
        for u, v in permutations([x for x in P if x != p], 2):
            K = det(vec(p), vec(u))*det(vec(u), vec(v))*det(vec(v), vec(p)) % 7
            need(((u, v) in relation) == (pow(K, 3, 7) == (1 if bit == 0 else 6)), 'chart equals bracket formula')
    # Reconstruct projectivities by their ordered images of infinity,0,1.
    maps = {}
    for aa, bb, cc in permutations(P, 3):
        va, vb, vc = vec(aa), vec(bb), vec(cc)
        alpha, beta = det(vc, vb), det(va, vc)
        matrix = (alpha*va[0] % 7, beta*vb[0] % 7,
                  alpha*va[1] % 7, beta*vb[1] % 7)
        determinant = (matrix[0]*matrix[3]-matrix[1]*matrix[2]) % 7
        need(determinant != 0 and tuple(apply(matrix, t) for t in (None, 0, 1)) == (aa, bb, cc), 'ordered triple reconstruction')
        action = tuple(apply(matrix, t) for t in P)
        maps[action] = int(pow(determinant, 3, 7) == 6)
    need(len(maps) == 336, 'sharp three-transitivity')
    states = Counter()
    for action, bit in maps.items():
        g = dict(zip(P, action))
        states[g[None], bit] += 1
        for (p, e), relation in tournaments.items():
            need(frozenset((g[u], g[v]) for u, v in relation) == tournaments[g[p], e ^ bit], 'chart arc transport')
    need(len(states) == 16 and set(states.values()) == {21}, 'homogeneous fibres')
    need(sum(action[-1] is None for action in maps) == 42, 'full pole stabilizer')
    root_gauge = (2, -1, 1, 0)
    need(apply(root_gauge, None) == 2 and 2*0-(-1)*1 == 1, 'square determinant root gauge')
    hcounts = []
    for (p, bit), relation in tournaments.items():
        hcounts.append(sum(all((w[i], w[i+1]) in relation for i in range(6))
                           for w in permutations([x for x in P if x != p])))
    need(set(hcounts) == {189}, 'literal Hamilton path enumeration')
    base = tournaments[None, 0]
    aut = [w for w in permutations(range(7))
           if all((w[u], w[v]) in base for u, v in base)]
    need(len(aut) == 21 and all(any(w == tuple((a*t+b) % 7 for t in range(7))
                                  for a in (1, 2, 4) for b in range(7)) for w in aut), 'literal full automorphism enumeration')
    need(189 // len(aut) == 9 and 189 % len(aut) == 0, 'normalized Hamilton count')

    # Use full vectors mod 105, rather than the producer's projective observer.
    q, seen = deque([(2, 1, 0)]), {(2, 1, 0)}
    while q:
        m, n, bit = q.popleft()
        for mm, nn, ee in ((2*m-n, m, bit), (2*m+n, m, bit ^ 1), (m+2*n, n, bit)):
            nxt = (mm % 105, nn % 105, ee)
            if nxt not in seen:
                seen.add(nxt)
                q.append(nxt)
    ambient = {(m, n, e) for m in range(105) for n in range(105) for e in (0, 1) if gcd(gcd(m, n), 105) == 1}
    need(seen == ambient and len(seen) == 18432, 'full vector reachability')
    fibres = Counter((slope((m, n)), bit, m % 15, n % 15) for m, n, bit in seen)
    need(len(fibres) == 3072 and set(fibres.values()) == {6}, 'six-to-one full-vector projection')
    primitive15 = {(m, n) for m in range(15) for n in range(15) if gcd(gcd(m, n), 15) == 1}
    need(len(primitive15) == 192 and set(fibres) == {(p, e, m, n) for p in P for e in (0, 1) for m, n in primitive15}, 'exact mod15 carrier')
    ray = [(m % 7, 0, m % 15, 1) for m in range(2, 212, 2)]
    survivors = [s for s in ray if all((s[2]+delta) % p != 0 for p in (3, 5) for delta in (-1, 1))]
    need(len(set(ray)) == 105 and len(survivors) == 21 and {s[2] for s in survivors} == {0, 3, 12}, 'literal D wheel')
    need((2+2*105) % 105 == 2 and all((2+2*k) % 105 != 2 for k in range(1, 105)), 'minimal D period')

    twins = [record('D'*k) for k in (5, 215)]
    for rec in twins:
        m, n = rec['pair']
        need(factor(m-1) == {m-1: 1} and factor(m+1) == {m+1: 1}, 'literal twin primality')
        need((m % 7, rec['bit'], m % 15, n % 15) == (5, 0, 12, 1), 'same refined state')
    need([r['hit'] for r in twins] == [False, True], 'opposite abc-hit predicates')
    need([r['radical'] for r in twins] == [124410, 8358844170], 'literal radicals')
    same_radical = [record(w) for w in ('UUUDUDDDU', 'DUUDDDDDUA')]
    need([r['pair'] for r in same_radical] == [[187, 100], [386, 153]], 'literal radical witness words')
    need([r['bit'] for r in same_radical] == [0, 1] and same_radical[0]['support'] == same_radical[1]['support'], 'same support opposite bits')
    need([r['radical'] for r in same_radical] == [2099693897070]*2, 'support product witness')
    controls = [record(w) for w in ('', 'UU', 'A', 'ADU')]
    need({(r['bit'], r['hit']) for r in controls} == {(0, False), (0, True), (1, False), (1, True)}, 'four sign-quality controls')
    need(integer_word('UA')[:2] == (8, 3) and integer_word('AU')[:2] == (8, 5), 'commutativization hostile')
    for pset in ((2,), (2, 3), (2, 3, 5), (2, 3, 7)):
        N = prod(pset)
        divisors = [d for d in range(1, N) if N % d == 0]
        need(len(divisors) == 2**len(pset)-1, 'augmentation counts proper divisors')
        need(sum((-1)**len(factor(d)) for d in divisors) == (-1)**(len(pset)+1), 'sign evaluation')
    print(json.dumps({'status': 'FINITE-EXACT independent audit controls', 'gates': GATES,
                      'projectivities_from_ordered_triples': len(maps),
                      'chart_tournaments': len(tournaments), 'hamilton_paths': sorted(set(hcounts)),
                      'full_automorphisms': len(aut), 'full_vector_states_mod105': len(seen),
                      'projected_states': len(fibres), 'projection_fibre_size': 6,
                      'D_period': len(ray), 'wheel_survivors': len(survivors),
                      'twin_quality_witness': twins, 'radical_bit_witness': same_radical}, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
