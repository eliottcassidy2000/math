#!/usr/bin/env python3
"""Self-contained exact checker. All conditions survive python -O."""
import itertools
import json
import math
import random
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def induced(rows, vertices):
    return [sum(1 << j for j, y in enumerate(vertices) if rows[x] >> y & 1)
            for x in vertices]


def check_tournament(rows):
    n = len(rows)
    require(all(0 <= x < (1 << n) and not (x >> i & 1)
                for i, x in enumerate(rows)), 'diagonal/range')
    require(all(((rows[i] >> j) ^ (rows[j] >> i)) & 1
                for i in range(n) for j in range(i)), 'tournament arcs')


def profiles(rows):
    """Terminal-removal DP; starting counts are DP on the converse."""
    check_tournament(rows)
    n = len(rows)
    N = 1 << n
    incoming = [sum(1 << i for i in range(n) if rows[i] >> j & 1)
                for j in range(n)]

    def endpoints(pred):
        values = [[0] * n for _ in range(N)]
        for S in range(1, N):
            scan = S
            while scan:
                bit = scan & -scan
                scan -= bit
                j = bit.bit_length() - 1
                W = S ^ bit
                if not W:
                    values[S][j] = 1
                else:
                    allowed = W & pred[j]
                    total = 0
                    while allowed:
                        last = allowed & -allowed
                        allowed -= last
                        total += values[W][last.bit_length() - 1]
                    values[S][j] = total
        return values

    end = endpoints(incoming)
    start = endpoints(rows)
    H = [1] + [sum(end[S]) for S in range(1, N)]
    full = N - 1
    result = []
    for r in range(n):
        R = sum(H[full ^ S] * start[S][r]
                for S in range(1, N) if S >> r & 1)
        result.append((H[full], end[full][r], R, H[full ^ (1 << r)]))
    return result, H


def literal_profile(rows, r):
    n = len(rows)
    h = e = R = 0
    for p in itertools.permutations(range(n)):
        edges = [bool(rows[p[i]] >> p[i + 1] & 1) for i in range(n - 1)]
        legal = all(edges)
        h += legal
        e += legal and p[-1] == r
        cut = p.index(r)
        R += all(edges[:max(0, cut - 1)]) and all(edges[cut:])
    S = induced(rows, [x for x in range(n) if x != r])
    d = sum(all(S[p[i]] >> p[i + 1] & 1 for i in range(len(S) - 1))
            for p in itertools.permutations(range(len(S))))
    return h, e, R, d


def components(rows):
    n = len(rows)
    reach = [rows[i] | (1 << i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    todo = set(range(n))
    result = []
    while todo:
        x = min(todo)
        C = sorted(y for y in todo if reach[x] >> y & 1 and reach[y] >> x & 1)
        result.append(C)
        todo.difference_update(C)
    snapshot = result[:]
    result.sort(key=lambda C: -sum(rows[C[0]] >> D[0] & 1
                                  for D in snapshot if D is not C))
    require(all(rows[x] >> y & 1 for i, C in enumerate(result)
                for D in result[i + 1:] for x in C for y in D), 'SCC order')
    return result


def normalize(rows, r):
    vertices = [v for v in range(len(rows)) if v != r]
    S = induced(rows, vertices)
    cs = components(S)
    O = {j for j, v in enumerate(vertices) if rows[r] >> v & 1}
    interior = 0
    if not O:
        keep = [x for C in cs if len(C) > 1 for x in C]
        t = 0
    else:
        special = [i for i, C in enumerate(cs) if O.intersection(C)]
        first, last = min(special), max(special)
        keep = []
        t = 0
        for i, C in enumerate(cs):
            if len(C) > 1 or C[0] in O or first < i < last:
                keep.extend(C)
                interior += len(C) == 1 and C[0] not in O
            elif i > last:
                t += 1
    retained = [vertices[x] for x in keep] + [r]
    return induced(rows, retained), len(retained) - 1, t, interior


def k_profile(z, t):
    if not t:
        return z
    h, e, R, d = z
    U, f = R - d, R + e - h - d
    return d - f + (1 << t) * U, d, d + (1 << t) * U, d


def k_graph(rows, r, t):
    for _ in range(t):
        n = len(rows)
        rows = rows[:] + [1 << r]
        for j in range(n):
            if j != r:
                rows[j] |= 1 << n
    return rows


def clone(rows, r, multiplicity):
    n = len(rows)
    exterior = [x for x in range(n) if x != r]
    labels = exterior + [r] * multiplicity
    out = [0] * len(labels)
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if labels[i] == labels[j]:
                out[i] |= 1 << j
            elif rows[labels[i]] >> labels[j] & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
    return out


def labelled(n):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        rows = [0] * n
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
        yield rows


def canonical(rows):
    return min(tuple(induced(rows, p)) for p in itertools.permutations(range(len(rows))))


def verify_normalization(rows, literal=False):
    allp, _ = profiles(rows)
    for r, z in enumerate(allp):
        h, e, R, d = z
        f, U = R + e - h - d, R - d
        k = (1 + math.isqrt(1 + 8 * f)) // 2
        ell = max(1, f).bit_length() - 1
        require(f >= 0 and f % 2 == 0, 'slack sign/parity')
        if U > 0:
            require(e <= h - e and h <= 2 * (U - f) - 1, 'endpoint insertion bound')
        require(math.comb(rows[r].bit_count(), 2) <= f, 'outdegree/slack bound')
        C, c, t, interior = normalize(rows, r)
        cp = profiles(C)[0][c]
        require(k_profile(cp, t) == z, ('normalization', rows, r, z, C, cp, t))
        require(len(C) <= 1 + d + k + ell, 'core size bound')
        require(interior <= ell, 'interior bound')
        if interior:
            require(1 << interior <= f, 'interior exponential certificate')
        vs = [v for v in range(len(rows)) if v != r]
        S = induced(rows, vs)
        sp = profiles(S)[0] if S else []
        require(U == sum(sp[j][2] for j, v in enumerate(vs) if rows[r] >> v & 1),
                'U transport identity')
        if literal:
            require(z == literal_profile(rows, r), 'independent permutation/cut profile')
    return len(allp)


def make_deletion(core, before, after):
    n = before + len(core) + after
    rows = [sum(1 << j for j in range(i + 1, n)) for i in range(n)]
    mask = ((1 << len(core)) - 1) << before
    for i, row in enumerate(core):
        rows[before + i] = (rows[before + i] & ~mask) | (row << before)
    return rows


def add_apex(S, out):
    n = len(S)
    rows = S[:] + [out]
    for i in range(n):
        if not out >> i & 1:
            rows[i] |= 1 << n
    return rows, n


def classify(core, d, f):
    k = (1 + math.isqrt(1 + 8 * f)) // 2
    limit = k + max(1, f).bit_length() - 1
    witnesses = {}
    examined = 0
    for singles in range(limit + 1):
        positions = range(singles + 1) if core else [0]
        for before in positions:
            S = make_deletion(core, before, singles - before)
            n = len(S)
            for out in range(1 << n):
                if out.bit_count() > k:
                    continue
                T, r = add_apex(S, out)
                C, c, t, interior = normalize(T, r)
                if len(C) != len(T):
                    continue
                examined += 1
                z = profiles(T)[0][r]
                h, e, R, dd = z
                require(dd == d, 'deletion core count')
                if R + e - h - dd == f and R > dd:
                    U = R - dd
                    witnesses.setdefault(U, {'rows': T, 'root': r, 'profile': list(z)})
    return examined, witnesses


def main():
    rooted = literal = 0
    strong_types = {3: set(), 5: set()}
    for n in range(1, 6):
        for rows in labelled(n):
            rooted += verify_normalization(rows, literal=n <= 4)
            if n <= 4:
                literal += n
            if n > 1 and len(components(rows)) == 1:
                h = profiles(rows)[0][0][0]
                if h in strong_types:
                    strong_types[h].add(canonical(rows))
    require(strong_types[3] == {canonical([4, 1, 2])}, 'unique strong H3 core')
    require(strong_types[5] == {canonical([8, 1, 3, 6])}, 'unique strong H5 core')

    rng = random.Random(197)
    larger_roots = 0
    for n, count in [(6, 160), (7, 80), (8, 20)]:
        for _ in range(count):
            rows = [0] * n
            for i, j in itertools.combinations(range(n), 2):
                a, b = (i, j) if rng.randrange(2) else (j, i)
                rows[a] |= 1 << b
            larger_roots += verify_normalization(rows)

    # Deliberately large neutral prefixes, interior gaps, and movable tails.
    padded = 0
    for core in [[4, 1, 2], [8, 1, 3, 6]]:
        for before, after in [(0, 5), (3, 3), (5, 0)]:
            S = make_deletion(core, before, after)
            for positions in [[before], [before, before + 1], [0, len(S) - 1]]:
                T, r = add_apex(S, sum(1 << x for x in positions))
                padded += verify_normalization(T)

    expected = {(1, 1): [1], (1, -1): [3, 5],
                (3, 1): [6], (3, -1): [9, 15],
                (5, 1): [11, 12, 16], (5, -1): [15, 39]}
    table = {}
    extension_checks = 0
    for core, d in [([], 1), ([4, 1, 2], 3), ([8, 1, 3, 6], 5)]:
        for beta in [1, -1]:
            examined, witnesses = classify(core, d, d - beta)
            require(sorted(witnesses) == expected[d, beta], ('table', d, beta, witnesses))
            for U, w in witnesses.items():
                z = tuple(w['profile'])
                for m in range(1, 4):
                    graph = k_graph(w['rows'], w['root'], m)
                    actual = profiles(graph)[0][w['root']]
                    require(actual == k_profile(z, m), 'K graph profile')
                    require(actual[0] == U * (1 << m) + beta, 'unit K family')
                    extension_checks += 1
            table[f'{d}:{beta}'] = {'examined_reduced_cores': examined,
                                    'U_base': sorted(witnesses), 'witnesses': witnesses}

    E = [2, 12, 9, 17, 7]
    h13 = [profiles(clone(E, 0, m))[0][0][0] for m in range(1, 7)]
    require(h13 == [6 * (1 << m) + 1 for m in range(1, 7)], 'H13 sign control')
    require(all(7 not in [u * (1 << j) for u in expected[d, beta] for j in range(4)]
                for beta in [-1, 1] for d in [1, 3, 5]), 'both amplitude 7 exclusions')

    out = {'status': 'PASS', 'labelled_root_checks_through_order_5': rooted,
           'independent_literal_profile_checks': literal,
           'larger_random_root_checks': larger_roots, 'padded_root_checks': padded,
           'unit_K_graph_checks': extension_checks,
           'strong_H3_core': sorted(strong_types[3]),
           'strong_H5_core': sorted(strong_types[5]),
           'unit_source_table': table, 'H13_clone_counts': h13,
           'scope': 'Finite base certificates plus universal proofs; full spectrum remains open.'}
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
