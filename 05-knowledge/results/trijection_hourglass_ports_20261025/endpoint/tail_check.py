"""Exact tail-family checks; shares only local arc/row primitive evaluators.

All reflected h=n-4 sources at n9..11; no middle, isomorphism or SCC filter.
n9/10 are diagnostics. n11 is a complete finite control of the proof's head.
The unique exceptional source is reconstructed independently by orbit transport.
"""
from itertools import product
from random import Random
from collections import Counter
from hashlib import sha256
import json
from check import need, edges, directed, blue, move, rows, row_move, row_blue, orbits, assemble, sample

EXCEPTION = (126, 220, 472, 976, 1952, 1102, 1936, 1825, 1571, 1063, 15)
REPAIR = ((0, 1, 2, 3, 4, 7, 5, 6, 8, 9, 10),
          (0, 1, 2, 3, 4, 7, 5, 9, 10, 6, 8),
          (0, 1, 2, 3, 4, 5, 9, 10, 6, 7, 8))


def ledgers(n, x_forward):
    if x_forward:
        return (tuple([0, 1, 2] + list(range(4, n)) + [3]),
                tuple([1, 2] + list(range(4, n - 4)) + [0] + list(range(n - 4, n)) + [3]),
                tuple(list(range(1, n - 4)) + [0] + list(range(n - 4, n))))
    V = tuple(range(4, n - 1))
    return ((0, 2, 3) + V + (n - 1, 1),
            (3, 0, 2) + V + (n - 1, 1),
            (3, 0, 1, 2) + V + (n - 1,))


def replay_words(T, qs, require_black):
    n = len(qs[0])
    P = tuple(range(n))
    now, A, xor, flags, witnesses = T, rows(T, n), set(), [], []
    for q in qs:
        need(set(q) == set(P) and len(q) == n, 'literal word')
        need(directed(now, q), 'actual intermediate directedness')
        flag = blue(now, q)
        need(flag == row_blue(A, q), 'arc/row anti agreement')
        flags.append(flag)
        r = dict(zip(q, q[::-1]))
        witnesses.append(next(((u, v, r[v], r[u]) for u, v in sorted(now)
                               if (r[v], r[u]) not in now), None))
        now, A, xor = move(now, q), row_move(A, q), xor ^ edges(q)
        need(rows(now, n) == A, 'arc/row K agreement')
    need(xor == edges(P) and now == move(T, P), 'literal K_P target')
    if require_black:
        need(not any(flags), 'three black certificate')
    return flags, witnesses


def source_check(T, n):
    e, y, c, z = n - 4, n - 3, n - 2, n - 1
    need(directed(T, tuple(range(n))) and blue(T, tuple(range(n))), 'marked source')
    need({u for u, v in T if v == 0} == {e, y, c, z}, 'exact incoming tail')
    need({v for u, v in T if u == z} == {0, 1, 2, 3}, 'exact reflected ports')


def proof_checks(T, n, x_forward, flags):
    e, y, c, z = n - 4, n - 3, n - 2, n - 1
    qs = ledgers(n, x_forward)
    S, B = move(T, qs[0]), move(T, tuple(range(n)))
    if x_forward:
        r1 = dict(zip(qs[0], qs[0][::-1]))
        r2 = dict(zip(qs[1], qs[1][::-1]))
        need((r1[0], r1[c]) == (3, 2), 'Q1 literal anti images')
        need((c, 0) in T and (3, 2) not in T, 'Q1 black witness')
        need((r2[1], r2[z]) == (3, 2), 'Q2 literal anti images')
        need((z, 1) not in S and (3, 2) in S, 'Q2 black witness')
        r = dict(zip(qs[2], qs[2][::-1]))
        g = [r[n - 1 - v] for v in range(n)]
        expected = [None] * n
        expected[:4] = [1, 2, 3, 4]
        for v in range(4, n - 6):
            expected[v] = v + 2
        expected[n - 6:] = [0, n - 4, n - 3, n - 2, n - 1, 5]
        need(g == expected, 'Q3 explicit g formula')
        seen, cycles = set(), []
        for root in range(n):
            if root in seen:
                continue
            cycle, v = [], root
            while v not in seen:
                seen.add(v)
                cycle.append(v)
                v = g[v]
            cycles.append(cycle)
        ports = {1, e, y, c, z}
        need({v for u, v in B if u == 0} == ports, 'Q3 exact five-point row')
        if n % 2:
            need(len(cycles) == 1, 'odd transitive orbit')
        else:
            need([len(cycle) for cycle in cycles] == [n // 2, n // 2], 'even two orbits')
            need(set(cycles[0]) & ports == {1}, 'even a-orbit has one outport')
        if flags[2]:
            need(n == 11 and tuple(rows(T, n)) == EXCEPTION,
                 'only exact n11 exceptional source')
    else:
        need(sum(u == 0 for u, v in S) == 5, 'R2 a outdegree five')
        need(sum(v == z for u, v in S) == 4, 'R2 z indegree four')
        r2 = dict(zip(qs[1], qs[1][::-1]))
        need(r2[0] == z, 'R2 pairs a and z')
        r3 = dict(zip(qs[2], qs[2][::-1]))
        g = [r3[n - 1 - v] for v in range(n)]
        need(all(g[g[g[g[v]]]] == v for v in range(n)) and
             any(g[g[v]] != v for v in range(n)), 'R3 product has order four')


def reconstruct_exception():
    n = 11
    q = ledgers(n, True)[2]
    r = dict(zip(q, q[::-1]))
    g = [r[n - 1 - v] for v in range(n)]
    cycle, v = [], 0
    while v not in cycle:
        cycle.append(v)
        v = g[v]
    need(len(cycle) == n, 'exception transitive cycle')
    pos = {v: i for i, v in enumerate(cycle)}
    ports = {1, 7, 8, 9, 10}
    B = {(u, v) for u in range(n) for v in range(n)
         if u != v and cycle[(pos[v] - pos[u]) % n] in ports}
    need(len(B) == n * (n - 1) // 2 and all((v, u) not in B for u, v in B),
         'transported row defines a tournament')
    T = move(B, tuple(range(n)))
    source_check(T, n)
    need(tuple(rows(T, n)) == EXCEPTION, 'independent exception reconstruction')
    need(blue(B, q), 'exception is genuinely blue for generic Q3')
    flags, witnesses = replay_words(T, REPAIR, True)
    return cycle, witnesses


def main():
    cycle, exceptional_witnesses = reconstruct_exception()
    count, digest = Counter(), sha256()
    for n in (9, 10, 11):
        fixed, free = orbits(n)
        for bits in product((0, 1), repeat=len(free)):
            T = assemble(n, n - 4, fixed, free, bits)
            source_check(T, n)
            x_forward = (2, 4) in T
            flags, _ = replay_words(T, ledgers(n, x_forward), False)
            count[f'exhaustive_n{n}'] += 1
            count[f'n{n}_generic_' + ('black' if not any(flags) else 'blue_word')] += 1
            if n >= 11:
                proof_checks(T, n, x_forward, flags)
                if any(flags):
                    need(flags == [False, False, True] and tuple(rows(T, n)) == EXCEPTION,
                         'unique generic-head exception')
                    replay_words(T, REPAIR, True)
                    count['n11_repaired_exception'] += 1
            digest.update(bytes((n, x_forward, sum(int(f) << i for i, f in enumerate(flags)))))
            digest.update(str(rows(T, n)).encode())
    rng = Random(202610250342)
    for n in range(12, 81):
        for x_forward in (False, True):
            for _ in range(4):
                T = sample(n, n - 4, rng, {(2, 4) if x_forward else (4, 2)})
                source_check(T, n)
                flags, _ = replay_words(T, ledgers(n, x_forward), True)
                proof_checks(T, n, x_forward, flags)
                count['sampled_n12_to80'] += 1
    print(json.dumps({'status': 'PASS; structural proof with explicit unique-head certificate',
                      'counts': dict(sorted(count.items())),
                      'exhaustive_digest': digest.hexdigest(),
                      'exception': {'source_rows': EXCEPTION, 'cycle': cycle, 'repair': REPAIR,
                                    'anti_failure_witnesses': exceptional_witnesses},
                      'filters': 'exact tail ports only; no middle, strongness, mixedness, clone or isomorphism filter',
                      'independence': 'arc sets versus bit rows; exception independently rebuilt from cycle transport'},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
