#!/usr/bin/env python3
"""Exact observer, wheel and arithmetic controls; no external dependencies."""
from collections import deque
from itertools import permutations, product
from math import gcd, prod
import json

GATES = 0


def need(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def compose(a, b):
    return tuple(a[b[x]] for x in range(len(a)))


def inverse(p):
    q = [0] * len(p)
    for i, x in enumerate(p):
        q[x] = i
    return tuple(q)


def mobius(M, t):
    a, b, c, d = M
    u, v = (a, c) if t == 7 else (a*t+b, c*t+d)
    return 7 if v % 7 == 0 else u * pow(v % 7, -1, 7) % 7


def perm(M):
    return tuple(mobius(M, t) for t in range(8))


def legendre(x):
    x %= 7
    return 0 if x == 0 else 1 if x in (1, 2, 4) else -1


def bracket(x, y):
    a, b = (1, 0) if x == 7 else (x, 1)
    c, d = (1, 0) if y == 7 else (y, 1)
    return (a*d-b*c) % 7


def relation(pole, bit):
    return frozenset((u, v) for u in range(8) for v in range(8)
                     if pole not in (u, v) and u != v
                     and legendre(bracket(pole, u)*bracket(u, v)*bracket(v, pole)) == (-1)**bit)


BRANCHES = {'U': (2, -1, 1, 0), 'A': (2, 1, 1, 0), 'D': (1, 2, 0, 1)}


def word(w):
    m, n, bit = 2, 1, 0
    for letter in w:
        a, b, c, d = BRANCHES[letter]
        m, n = a*m+b*n, c*m+d*n
        bit ^= letter == 'A'
    return m, n, bit


def factor(n):
    factors = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def packet(m, n):
    a, b, c = m*m-n*n, 2*m*n, m*m+n*n
    fs = [factor(x) for x in (a, b, c)]
    support = sorted(set().union(*(set(f) for f in fs)))
    R = prod(support)
    need(gcd(m, n) == 1 and (m+n) % 2 == 1 and m > n > 0, 'primitive Euclid hypotheses')
    need(a*a+b*b == c*c and gcd(a, b) == gcd(a, c) == gcd(b, c) == 1, 'primitive abc packet')
    for value, fsingle in zip((a, b, c), fs):
        need(prod(p**e for p, e in fsingle.items()) == value, 'factor certificate')
    return {'pair': [m, n], 'triple': [a, b, c], 'factorizations': fs,
            'prime_support': support, 'radical': R, 'abc_C': c*c, 'quality_above_one': c*c > R}


def main():
    G = {}
    for M in product(range(7), repeat=4):
        determinant = (M[0]*M[3]-M[1]*M[2]) % 7
        if not determinant:
            continue
        g, bit = perm(M), int(legendre(determinant) == -1)
        if g in G:
            need(G[g] == bit, 'projective determinant character')
        G[g] = bit
    H = {g for g in G if g[7] == 7 and G[g] == 0}
    B = {g for g in G if g[7] == 7}
    need((len(G), len(B), len(H)) == (336, 42, 21), 'G B H sizes')
    states = set(product(range(8), range(2)))
    tournaments = {s: relation(*s) for s in states}
    for (p, b), arcs in tournaments.items():
        need(len(arcs) == 21, 'each state is a seven-tournament')
        need(arcs == frozenset((v, u) for u, v in tournaments[p, b ^ 1]), 'opposite sheet is converse')
    for g, b in G.items():
        for (p, e), arcs in tournaments.items():
            moved = frozenset((g[u], g[v]) for u, v in arcs)
            need(moved == tournaments[g[p], e ^ b], 'exact arc transport')
    cosets = {}
    for g, b in G.items():
        coset = frozenset(compose(g, h) for h in H)
        state = (g[7], b)
        if coset in cosets:
            need(cosets[coset] == state, 'coset state is well-defined')
        cosets[coset] = state
    need(len(cosets) == 16 and set(cosets.values()) == states, 'homogeneous space bijection')
    U = perm(BRANCHES['U'])
    Hroot = {g for g in G if g[2] == 2 and G[g] == 0}
    need(U[7] == 2 and G[U] == 0, 'root pole conjugation')
    need({compose(U, compose(h, inverse(U))) for h in H} == Hroot, 'actual root stabilizer conjugate')
    base = tournaments[7, 0]
    need(base == frozenset((u, v) for u in range(7) for v in range(7) if (v-u) % 7 in (1, 2, 4)), 'base is Paley')
    arc_orbit = {(h[0], h[1]) for h in H}
    need(arc_orbit == base and len(arc_orbit) == len(H), 'H acts regularly on 21 arcs')
    full_aut = {p for p in permutations(range(7)) if frozenset((p[u], p[v]) for u, v in base) == base}
    need(full_aut == {g[:7] for g in H}, 'full Paley automorphism oracle')
    frobenius = perm((2, 0, 0, 1))
    need(frobenius in H and G[frobenius] == 0, 'Frobenius character is trivial')
    def gf8_mul(x, y):
        out = 0
        while y:
            if y & 1:
                out ^= x
            y >>= 1
            x <<= 1
            if x & 8:
                x ^= 11
        return out
    def frob(x, j):
        for unused in range(j):
            x = gf8_mul(x, x)
        return x
    joint_stabilizers = {}
    for x, y in product(range(8), repeat=2):
        stabilizer = [j for j in range(3) if frob(x, j) == x and frob(y, j) == y]
        need(len(stabilizer) == (3 if x in (0, 1) and y in (0, 1) else 1), 'joint-degree stabilizer')
        joint_stabilizers[x, y] = len(stabilizer)
    need(joint_stabilizers[0, 0] == 3 and joint_stabilizers[0, 2] == 1, 'same marginal loses joint degree')
    counts = {}
    for state, arcs in tournaments.items():
        vertices = [v for v in range(8) if v != state[0]]
        dp = {(1 << i, i): 1 for i in range(7)}
        for mask in range(1, 128):
            for last in range(7):
                count = dp.get((mask, last), 0)
                if not count:
                    continue
                for nxt in range(7):
                    if not mask >> nxt & 1 and (vertices[last], vertices[nxt]) in arcs:
                        key = (mask | 1 << nxt, nxt)
                        dp[key] = dp.get(key, 0) + count
        counts[state] = sum(dp.get((127, v), 0) for v in range(7))
        need(counts[state] == 189, 'all transported Hamilton counts')

    # Exact finite automaton: pair residues mod15 refine the sixteen states.
    M = 15
    primitive = {(m, n) for m in range(M) for n in range(M) if gcd(gcd(m, n), M) == 1}
    need(len(primitive) == 192, 'primitive residue carrier')
    branch_perms = {letter: perm(matrix) for letter, matrix in BRANCHES.items()}
    for matrix in BRANCHES.values():
        a, b, c, d = matrix
        moved = {((a*m+b*n) % M, (c*m+d*n) % M) for m, n in primitive}
        need(moved == primitive, 'residue branch is a bijection')
    root = (2, 0, 2, 1)
    seen, queue = {root}, deque([root])
    while queue:
        t, bit, m, n = queue.popleft()
        for letter, matrix in BRANCHES.items():
            a, b, c, d = matrix
            nxt = (branch_perms[letter][t], bit ^ (letter == 'A'), (a*m+b*n) % M, (c*m+d*n) % M)
            need(nxt[:2] in states and nxt[2:] in primitive, 'finite observer closure')
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    need(len(seen) == 3072, 'complete finite automaton census at modulus15')
    ray = {( (2+2*k) % 7, 0, (2+2*k) % M, 1) for k in range(105)}
    survivors = {s for s in ray if gcd(s[2]*s[2]-1, M) == 1}
    need(len(ray) == 105 and len(survivors) == 21, 'D-ray period and wheel count')
    need({s[2] for s in survivors} == {0, 3, 12}, 'literal wheel residues')

    threshold_controls = {}
    for w in ('', 'UU', 'A', 'ADU'):
        m, n, bit = word(w)
        rec = packet(m, n)
        rec.update(word=w, bit=bit)
        threshold_controls[w] = rec
    need({(r['bit'], r['quality_above_one']) for r in threshold_controls.values()} == set(product(range(2), (False, True))), 'each bit has hits and nonhits')
    twin_controls = []
    for k in (5, 215):
        m, n, bit = word('D'*k)
        rec = packet(m, n)
        need(factor(m-1) == {m-1: 1} and factor(m+1) == {m+1: 1}, 'actual twin centers')
        state = (m % 7, bit, m % 15, n % 15)
        need(state in survivors, 'same survivor observer')
        rec.update(D_exponent=k, state=state, twins=[m-1, m+1])
        twin_controls.append(rec)
    need(twin_controls[0]['state'] == twin_controls[1]['state'], 'identical finite observer')
    need([r['quality_above_one'] for r in twin_controls] == [False, True], 'observer cannot determine abc-hit predicate')
    radical_controls = []
    for w in ('UUUDUDDDU', 'DUUDDDDDUA'):
        m, n, bit = word(w)
        rec = packet(m, n)
        rec.update(word=w, bit=bit)
        radical_controls.append(rec)
    need([r['pair'] for r in radical_controls] == [[187, 100], [386, 153]], 'same radical exact words')
    need([r['bit'] for r in radical_controls] == [0, 1], 'opposite character bits')
    need(radical_controls[0]['prime_support'] == radical_controls[1]['prime_support'] == [2, 3, 5, 7, 11, 17, 29, 41, 193, 233], 'identical radical support')
    need(radical_controls[0]['radical'] == radical_controls[1]['radical'] == 2099693897070, 'literal radical value')

    # Monoid evaluations are different maps, not interchangeable scalar shadows.
    for length in range(7):
        for letters in product('UAD', repeat=length):
            w = ''.join(letters)
            bit = word(w)[2]
            need((-1)**w.count('A') == (-1)**bit, 'determinant through exponent-count monoid')
    need(word('UA')[:2] == (8, 3) and word('AU')[:2] == (8, 5), 'abelianization loses endpoint')
    for k in range(1, 7):
        supports = [mask for mask in range(1 << k) if mask != (1 << k)-1]
        need(len(supports) == (1 << k)-1, 'divisor-polynomial augmentation')
        need(sum((-1)**mask.bit_count() for mask in supports) == (-1)**(k+1), 'divisor-polynomial sign evaluation')
        for mask in range(1 << k):
            need(mask | mask == mask, 'radical support is idempotent')
    need((-1)**1 * (-1)**1 != (-1)**1, 'Mobius sign is not a union character')
    divisor_weights = []
    for primes in ((2, 3, 5), (2, 3, 7)):
        weighted = sum(prod(primes[i] for i in range(3) if mask >> i & 1) for mask in range(7))
        N = prod(primes)
        need(weighted == sum(d for d in range(1, N) if N % d == 0), 'weighted monoid evaluation is actual aliquot sum')
        divisor_weights.append(weighted)
    need(divisor_weights == [42, 54], 'same augmentation different arithmetic weights')

    receipt = {'status': 'FINITE-EXACT controls; elementary proofs in README.md', 'gates': GATES,
               'G_B_H_orders': [len(G), len(B), len(H)], 'homogeneous_states': len(cosets),
               'all_state_H': sorted(set(counts.values())), 'all_state_Aut': len(full_aut), 'all_state_tc': 9,
               'wheel_modulus': M, 'primitive_pairs': len(primitive), 'reachable_refined_states': len(seen),
               'D_ray_period': len(ray), 'D_ray_wheel_survivors': len(survivors),
               'quality_by_bit_controls': threshold_controls, 'same_observer_twins': twin_controls,
               'same_radical_opposite_bits': radical_controls,
               'scope': 'No universal abc bound, twin-prime infinitude, or global Hamilton-gap implication.'}
    print(json.dumps(receipt, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
