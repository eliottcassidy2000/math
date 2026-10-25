"""Literal checks for cyclic-prefix repair. Standard library; explicit gates.

Exhaustive small controls: all stated guards, n9/10. Deterministic larger
controls: n11..60, every legal h, two sources each. No class or SCC filter.
No inherited producer imports. The all-order conclusion rests on THEOREM.md.
"""
from itertools import combinations, product
from random import Random
from collections import Counter
from hashlib import sha256
import json


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def pair(u, v):
    return (min(u, v), max(u, v))


def edges(q):
    return {pair(u, v) for u, v in zip(q, q[1:])}


def directed(T, q):
    return all((u, v) in T for u, v in zip(q, q[1:]))


def blue(T, q):
    r = dict(zip(q, q[::-1]))
    return all((r[v], r[u]) in T for u, v in T)


def move(T, q):
    E = edges(q)
    return {(u, v) if pair(u, v) in E else (v, u) for u, v in T}


def rows(T, n):
    out = [0] * n
    for u, v in T:
        out[u] |= 1 << v
    return out


def row_move(A, q):
    n = len(A)
    out = [((1 << n) - 1) ^ row ^ (1 << u) for u, row in enumerate(A)]
    for u, v in zip(q, q[1:]):
        out[u] ^= 1 << v
        out[v] ^= 1 << u
    return out


def row_blue(A, q):
    n = len(A)
    r = [0] * n
    for i, v in enumerate(q):
        r[v] = q[n - 1 - i]
    return all(bool(A[u] >> v & 1) == bool(A[r[v]] >> r[u] & 1)
               for u in range(n) for v in range(n) if u != v)


def orbits(n):
    seen = set()
    fixed, free = [], []
    for u, v in combinations(range(1, n - 1), 2):
        if (u, v) in seen:
            continue
        orb = {(u, v), (n - 1 - v, n - 1 - u)}
        seen |= orb
        (fixed if v == u + 1 else free).append(orb)
    return fixed, free


def assemble(n, h, fixed, free, bits):
    T = set().union(*fixed)
    for orb, bit in zip(free, bits):
        T.update(orb if bit else {(v, u) for u, v in orb})
    I = {h, n - 3, n - 2}
    J = {n - 1 - v for v in I}
    for v in range(1, n - 1):
        T.add((v, 0) if v in I else (0, v))
        T.add((n - 1, v) if v in J else (v, n - 1))
    T.add((n - 1, 0))
    return T


def sample(n, h, rng, guards):
    fixed, free = orbits(n)
    bits = []
    for orb in free:
        forward = bool(orb & guards)
        backward = bool({(v, u) for u, v in orb} & guards)
        need(not (forward and backward), 'compatible reflected guards')
        bits.append(1 if forward else 0 if backward else rng.getrandbits(1))
    T = assemble(n, h, fixed, free, bits)
    need(guards <= T, 'all requested source guards')
    return T


def ledger(n, h):
    W = tuple(range(1, h))
    H = tuple(range(h, n - 3))
    a, y, c, z = 0, n - 3, n - 2, n - 1
    return ((z, a) + W + H + (y, c),
            (z, a) + H + W + (y, c),
            (a,) + H + W + (y, c, z))


def guards(n, h):
    return {(n - 2, 1), (3, n - 2), (n - h, 2)}


def replay(T, n, h, proof=True):
    P = tuple(range(n))
    I = {h, n - 3, n - 2}
    need(len(T) == n * (n - 1) // 2, 'tournament size')
    need(all((v, u) not in T for u, v in T), 'tournament pair orientation')
    need(directed(T, P) and blue(T, P), 'directed blue source')
    need({u for u, v in T if v == 0} == I | {n - 1}, 'exact incoming ports')
    need(guards(n, h) <= T, 'legal source guards')
    qs = ledger(n, h)
    S, B = move(T, qs[0]), move(T, P)
    now, A, xor, flags = T, rows(T, n), set(), []
    for q in qs:
        need(len(q) == n and set(q) == set(P), 'literal permutation word')
        need(directed(now, q), 'directed at actual intermediate tournament')
        flag = blue(now, q)
        need(flag == row_blue(A, q), 'independent blue predicate')
        flags.append(flag)
        now, A = move(now, q), row_move(A, q)
        need(rows(now, n) == A, 'independent K evaluation')
        xor ^= edges(q)
    need(xor == edges(P) and now == B, 'literal partner and XOR')
    need(not any(flags), 'three black words')
    if proof:
        need(n >= 11, 'claimed proof range')
        r = dict(zip(qs[1], qs[1][::-1]))
        e, y, c, z = n - 4, n - 3, n - 2, n - 1
        w = 2 * h - 1 if 2 * h <= n - 3 else 2 * h - n + 3
        need(r[0] == y and r[e] == w, 'Q2 positional formula')
        if not (n % 2 == 0 and h == (n - 2) // 2):
            need((e, y) in S and (0, w) not in S, 'Q2 retained-edge witness')
        elif n != 12:
            need(r[3] == n - 6 and r[c] == z, 'Q2 exceptional reflection')
            need((3, c) not in S and (z, n - 6) in S,
                 'Q2 exceptional chord witness')
        else:
            need(h == 5 and r[9] == 0 and r[1] == 8 and r[10] == 11 and r[2] == 7,
                 'n12 anti equations')
            need((0, 8) in T and ((2, 10) in T) == ((1, 9) in T),
                 'n12 original source identities')
            need((11, 7) in S, 'n12 final true side')
            need(((1, 9) in S) != ((0, 8) in S) or
                 ((2, 10) in S) != ((11, 7) in S), 'n12 two-equation contradiction')
        r = dict(zip(qs[2], qs[2][::-1]))
        g = lambda v: r[n - 1 - v]
        ports = {1, h, y, c}
        need(g(0) == 0 and g(z) == z and {v for u, v in B if u == 0} == ports | {z},
             'Q3 fixed endpoint and exact port set')
        need({g(v) for v in ports} != ports, 'Q3 port-set obstruction')
    return rows(T, n)


def main():
    count = Counter()
    digest = sha256()
    for n in (9, 10):
        fixed, free = orbits(n)
        for h in range(2, n - 3):
            for bits in product((0, 1), repeat=len(free)):
                T = assemble(n, h, fixed, free, bits)
                if not guards(n, h) <= T:
                    continue
                source = replay(T, n, h, proof=False)
                count['exhaustive_n' + str(n)] += 1
                digest.update(str((n, h, source)).encode())
    rng = Random(202610250341)
    for n in range(11, 61):
        for h in range(2, n - 3):
            for _ in range(2):
                T = sample(n, h, rng, guards(n, h))
                replay(T, n, h)
                count['sampled_n11_to60'] += 1
                if n == 12 and h == 5:
                    count['n12_exception'] += 1
    # A new h=3 source defeats all old and the sibling extended odd-prefix
    # guards; also rejects the growing-prefix c->4 guard.
    n, h = 19, 3
    G = guards(n, h) | {(n - 2, 2), (4, n - 2)}
    applicable = []
    for ell in range(1, n - 3, 2):
        if n < max(15, 2 * ell + 3):
            continue
        if h in {2, ell + 2} or (ell >= 3 and h in {4, ell - 1}):
            continue
        applicable.append(ell)
        G.add((2, ell + 3))
    T = sample(n, h, rng, G)
    new_source = replay(T, n, h)
    need((n - 3, 2) in T and (n - 2, 2) in T, 'old residual chord branch')
    need((n - 2, 4) not in T, 'reject inherited growing-prefix guard')
    need(all((ell + 3, 2) not in T for ell in applicable), 'reject every odd-prefix guard')
    hostiles = []
    for bad in ((3, n - 2), (n - h, 2)):
        Gbad = guards(n, h) - {bad}
        Gbad.add((bad[1], bad[0]))
        H = sample(n, h, rng, Gbad)
        qs = ledger(n, h)
        S = move(H, qs[0])
        need(directed(H, qs[0]) and not directed(S, qs[1]), 'omitted guard illegal Q2')
        hostiles.append({'reversed_guard': bad, 'rows': rows(H, n),
                         'failed_joins': [(u, v) for u, v in zip(qs[1], qs[1][1:]) if (u, v) not in S]})
    # Independent permutation-only control of port-set loss, with no source generator.
    for n in range(11, 201):
        for h in range(2, n - 3):
            q = ledger(n, h)[2]
            r = dict(zip(q, q[::-1]))
            ports = {1, h, n - 3, n - 2}
            need({r[n - 1 - v] for v in ports} != ports, 'port permutation control')
            count['permutation_pairs_n11_to200'] += 1
    print(json.dumps({'status': 'PASS; structural proof supplemented by bounded controls',
                      'counts': dict(sorted(count.items())), 'small_source_digest': digest.hexdigest(),
                      'newly_covered': {'n': 19, 'h': 3, 'rows': new_source,
                                        'blocked_odd_prefix_lengths': applicable},
                      'omitted_guard_hostiles': hostiles,
                      'filters': 'literal stated guards only; no mixedness, strongness, clone or isomorphism filter'},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
