"""Bounded all-k controls and direct exceptional-family reconstruction.

No exhaustive tournament census is claimed. All gates use explicit exceptions.
Shares only local arc/row primitive evaluators and the literal replay kernel.
"""
from random import Random
from collections import Counter
import json
from check import need, orbits, directed, blue, move, rows
from tail_check import replay_words


def generic(n, k, x_forward):
    e = n - k
    if x_forward:
        return (tuple([0, 1, 2] + list(range(4, n)) + [3]),
                tuple([1, 2] + list(range(4, e)) + [0] + list(range(e, n)) + [3]),
                tuple(list(range(1, e)) + [0] + list(range(e, n))))
    V = tuple(range(4, n - 1))
    return ((0, 2, 3) + V + (n - 1, 1),
            (3, 0, 2) + V + (n - 1, 1),
            (3, 0, 1, 2) + V + (n - 1,))


def random_source(n, k, x_forward, rng):
    fixed, free = orbits(n)
    T = set().union(*fixed)
    guard = (2, 4) if x_forward else (4, 2)
    for orb in free:
        reverse = {(v, u) for u, v in orb}
        bit = 1 if guard in orb else 0 if guard in reverse else rng.getrandbits(1)
        T.update(orb if bit else reverse)
    for v in range(1, n - 1):
        T.add((v, 0) if v >= n - k else (0, v))
        T.add((n - 1, v) if v < k else (v, n - 1))
    T.add((n - 1, 0))
    return T


def source_check(T, n, k):
    need(len(T) == n * (n - 1) // 2 and all((v, u) not in T for u, v in T), 'tournament')
    need(directed(T, tuple(range(n))) and blue(T, tuple(range(n))), 'source marked path')
    need({u for u, v in T if v == 0} == set(range(n - k, n)), 'exact incoming tail')
    need({v for u, v in T if u == n - 1} == set(range(k)), 'exact outgoing prefix')


def exceptional(k):
    n, e = 2 * k + 3, k + 3
    cycle = list(range(k + 1)) + list(range(k + 2, n)) + [k + 1]
    pos = {v: i for i, v in enumerate(cycle)}
    D = {1} | set(range(k + 2, 2 * k + 2))
    need(all((d in D) != ((n - d) in D) for d in range(1, n)), 'one of every signed pair')
    B = {(u, v) for u in range(n) for v in range(n)
         if u != v and (pos[v] - pos[u]) % n in D}
    T = move(B, tuple(range(n)))
    source_check(T, n, k)
    need((2, 4) in T, 'exception is x-forward')
    need(all(pos[n - 1 - v] == (2 * k + 1 - pos[v]) % n for v in range(n)),
         'P reflection in cyclic coordinates')
    q = generic(n, k, True)[2]
    need(blue(B, q), 'generic final word really blue at exceptional family')
    r = dict(zip(q, q[::-1]))
    need(all(r[n - 1 - cycle[i]] == cycle[(i + 1) % n] for i in range(n)),
         'generic reflection product is stated cycle')
    A, Z = tuple(range(k + 1)), tuple(range(k + 5, n))
    u, v, w, s = k + 1, k + 2, k + 3, k + 4
    repair = (A + (w, u, v, s) + Z, A + (w, u) + Z + (v, s), A + (u,) + Z + (v, w, s))
    replay_words(T, repair, True)
    S1 = move(T, repair[0])
    for now, q, p, target in ((T, repair[0], (0, v), (k - 1, n - 1)),
                               (S1, repair[1], (0, w), (w, s)),
                               (B, repair[2], (0, 3), (n - 1, s))):
        r = dict(zip(q, q[::-1]))
        need((r[p[1]], r[p[0]]) == target, 'uniform repair actual reflected pair')
        need((p in now) != (target in now), 'uniform repair literal black witness')
    return T, B, repair


def proof_checks(T, n, k, x_forward):
    qs = generic(n, k, x_forward)
    S, B = move(T, qs[0]), move(T, tuple(range(n)))
    if x_forward:
        need((n - 2, 0) in T and (3, 2) not in T, 'generic Q1 witness values')
        need((n - 1, 1) not in S and (3, 2) in S, 'generic Q2 witness values')
        r = dict(zip(qs[2], qs[2][::-1]))
        g = [r[n - 1 - v] for v in range(n)]
        e = n - k
        expected = list(range(1, k + 1)) + list(range(k + 2, e)) + [0, e] + list(range(e + 1, n)) + [k + 1]
        need(g == expected, 'all-k explicit g formula')
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
        ports = {1} | set(range(e, n))
        need({v for u, v in B if u == 0} == ports, 'all-k exact B endpoint row')
        if n % 2:
            need(len(cycles) == 1, 'all-k odd transitive orbit')
        else:
            need([len(c) for c in cycles] == [n // 2, n // 2], 'all-k even two orbits')
            need(set(cycles[0]) & ports == {1}, 'all-k even one outport')
    else:
        need(sum(u == 0 for u, v in S) == k + 1, 'opposite-chord a outdegree')
        need(sum(v == n - 1 for u, v in S) == k, 'opposite-chord z indegree')
        r = dict(zip(qs[1], qs[1][::-1]))
        need(r[0] == n - 1, 'second word pairs incompatible endpoints')
        r = dict(zip(qs[2], qs[2][::-1]))
        g = [r[n - 1 - v] for v in range(n)]
        need(all(g[g[g[g[v]]]] == v for v in range(n)) and any(g[g[v]] != v for v in range(n)),
             'all-k final product order four')


def main():
    rng, counts = Random(202610250343), Counter()
    for k in range(4, 21):
        exception_T, _, repair = exceptional(k)
        for n in range(2 * k + 3, 2 * k + 23):
            for x_forward in (False, True):
                for _ in range(2):
                    T = random_source(n, k, x_forward, rng)
                    source_check(T, n, k)
                    flags, _ = replay_words(T, generic(n, k, x_forward), False)
                    proof_checks(T, n, k, x_forward)
                    if any(flags):
                        need(n == 2 * k + 3 and T == exception_T and flags == [False, False, True],
                             'only uniquely forced generic boundary exception')
                        replay_words(T, repair, True)
                    counts['sampled_k4_to20_twenty_orders_both_chords'] += 1
    for k in range(4, 81):
        exceptional(k)
        counts['direct_exception_family_k4_to80'] += 1
    # Hostile: blindly translating the six last labels of k=4 loses the moving junction.
    k, n = 5, 13
    T, _, _ = exceptional(k)
    prefix = tuple(range(n - 6))
    wrong = prefix + (n - 4, n - 6, n - 5, n - 3, n - 2, n - 1)
    failed = [(u, v) for u, v in zip(wrong, wrong[1:]) if (u, v) not in T]
    need(bool(failed), 'translated-six-label hostile')
    print(json.dumps({'status': 'PASS; scoped all-k theorem, bounded exact controls',
                      'counts': dict(sorted(counts.items())),
                      'naive_translation_hostile': {'k': k, 'n': n, 'rows': rows(T, n), 'word': wrong,
                                                   'failed_joins': failed},
                      'universe': 'k4..20, n2k+3..2k+22, two deterministic sources per chord; exact exceptional carrier k4..80',
                      'filters': 'consecutive tail ports only; no middle or quotient filter'},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
