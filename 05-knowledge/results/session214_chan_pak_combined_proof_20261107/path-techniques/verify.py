#!/usr/bin/env python3
"""Independent finite controls for the sole-exit fixed-defect theorem.
No repository imports; default universe is self-contained.
"""
import argparse
import hashlib
import itertools
import json
import random
from collections import Counter
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def validate(a):
    n = len(a)
    need(all(not (a[i] >> i & 1) and a[i] < (1 << n) for i in range(n)), 'rows/loops')
    need(all(((a[i] >> j & 1) + (a[j] >> i & 1)) == 1
             for i in range(n) for j in range(i + 1, n)), 'tournament')


def sub(a, vertices):
    return [sum(((a[v] >> w) & 1) << j for j, w in enumerate(vertices)) for v in vertices]


def dual(a):
    return [sum(((a[j] >> i) & 1) << j for j in range(len(a))) for i in range(len(a))]


def endings(a):
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        rest = mask
        while rest:
            bit = rest & -rest
            rest -= bit
            v = bit.bit_length() - 1
            old = mask ^ bit
            if old:
                pred = old
                value = 0
                while pred:
                    b = pred & -pred
                    pred -= b
                    u = b.bit_length() - 1
                    if a[u] & bit:
                        value += dp[old][u]
                dp[mask][v] = value
    return dp


def tables(a):
    e = endings(a)
    starts = endings(dual(a))
    h = [sum(row) for row in e]
    h[0] = 1
    return h, e, starts


def covers(h, e, starts, umask, omask):
    L = K = 0
    q = umask
    while q:
        p = umask ^ q
        so = sum(starts[q][v] for v in range(len(starts[q])) if omask >> v & 1)
        eo = sum(e[p][v] for v in range(len(e[p])) if omask >> v & 1)
        L += h[p] * so
        K += eo * so
        q = (q - 1) & umask
    return L, K


def components(a):
    n = len(a)
    reach = [a[i] | (1 << i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    groups = []
    unseen = set(range(n))
    while unseen:
        v = min(unseen)
        c = [w for w in sorted(unseen) if reach[v] >> w & 1 and reach[w] >> v & 1]
        groups.append(c)
        unseen.difference_update(c)
    groups = sorted(groups, key=lambda c: -sum(a[c[0]] >> g[0] & 1 for g in groups if g != c))
    need(all(a[x] >> y & 1 for i, c in enumerate(groups) for d in groups[i+1:] for x in c for y in d), 'SCC order')
    return groups


def sole_exit_extension(u, omask):
    n = len(u)
    r, mark = n, n + 1
    a = [row | (1 << mark) for row in u] + [omask, 1 << r]
    for i in range(n):
        if not (omask >> i & 1):
            a[i] |= 1 << r
    return a, mark, r


def chain(n):
    return [sum(1 << j for j in range(i + 1, n)) for i in range(n)]


def triangle_chain(p, q):
    a = chain(p + q + 3)
    a[p] ^= 1 << (p + 2)
    a[p + 2] |= 1 << p
    return a


def power_two(n):
    return n > 0 and (n & (n - 1)) == 0


def in_defect3_spectrum(c):
    return (c == 1 or (c % 3 == 0 and power_two(c // 3))
            or ((c + 2) % 3 == 0 and (c + 2) // 3 >= 2 and power_two((c + 2) // 3))
            or ((c + 2) % 5 == 0 and (c + 2) // 5 >= 2 and power_two((c + 2) // 5)))


def classify_strong(a, mark, defect):
    r = a[mark].bit_length() - 1
    verts = [v for v in range(len(a)) if v not in (mark, r)]
    u = sub(a, verts)
    O = {i for i, v in enumerate(verts) if a[r] >> v & 1}
    cs = components(u)
    if not cs or not O.intersection(cs[0]):
        return False
    if defect == 1:
        return len(O) == 1 and all(len(c) == 1 for c in cs)
    need(defect == 3, 'classifier defect')
    if all(len(c) == 1 for c in cs):
        positions = sorted(i for i, c in enumerate(cs) if c[0] in O)
        return positions in ([0, 1], [0, 2])
    return (len(O) == 1 and [len(c) for c in cs if len(c) > 1] == [3])


def classify_defect3(a, mark):
    if a[mark].bit_count() != 1:
        return False
    cs = components(a)
    if mark not in cs[-1] or len(cs[-1]) < 3:
        return False
    pre_nontrivial = [len(c) for c in cs[:-1] if len(c) > 1]
    if pre_nontrivial == []:
        d = 3
    elif pre_nontrivial == [3]:
        d = 1
    else:
        return False
    last = cs[-1]
    return classify_strong(sub(a, last), last.index(mark), d)


def inspect(a, counts, exhaustive_compression=False):
    validate(a)
    n = len(a)
    h, e, starts = tables(a)
    total = h[-1]
    need(total % 2 == 1, 'oddness')
    cs = components(a)
    strong = len(cs) == 1 and n >= 3
    for mark in range(n):
        if a[mark].bit_count() != 1:
            continue
        counts['sole_exit_marks'] += 1
        r = a[mark].bit_length() - 1
        umask = ((1 << n) - 1) ^ (1 << mark) ^ (1 << r)
        omask = a[r] & umask
        L, K = covers(h, e, starts, umask, omask)
        c = e[-1][mark]
        d = total - 2 * c
        need(c == L - K and d == h[umask] + K, ('identity', a, mark, total, c, L, K))
        need(K >= omask.bit_count() * (omask.bit_count() - 1), 'capacity weak bound')
        pairfloor = 2 * sum(h[umask ^ (1 << x) ^ (1 << y)]
                            for x in range(n) if omask >> x & 1
                            for y in range(x + 1, n) if omask >> y & 1)
        need(K >= pairfloor, 'capacity weighted bound')
        need((c > 0 and d == 3) == classify_defect3(a, mark), ('defect3 structure', a, mark, c, d))
        if c > 0 and d == 3:
            counts['defect3_marks'] += 1
            need(in_defect3_spectrum(c), ('defect3 scalar', a, mark, c))
        if not strong or not exhaustive_compression:
            continue
        counts['compression_cases'] += 1
        verts = [v for v in range(n) if umask >> v & 1]
        u = sub(a, verts)
        O = {i for i, v in enumerate(verts) if omask >> v & 1}
        ucs = components(u)
        need(O.intersection(ucs[0]), 'strong entry')
        last = max(i for i, C in enumerate(ucs) if O.intersection(C))
        removed = {v for C in ucs[last+1:] if len(C) == 1 for v in C}
        keep = [v for v in range(len(u)) if v not in removed]
        u0 = sub(u, keep)
        o0 = sum(1 << i for i, v in enumerate(keep) if v in O)
        h0, e0, s0 = tables(u0)
        A, K0 = covers(h0, e0, s0, (1 << len(u0)) - 1, o0)
        need(h0[-1] == h[umask] and K0 == K and L == (1 << len(removed)) * A, 'compression')
        B = d + (max(1, d - 1).bit_length() - 1) + 4
        need(len(u0) + 2 <= B, ('skeleton bound', a, mark, d, len(u0) + 2, B))
        t0, m0, _ = sole_exit_extension(u0, o0)
        need(len(components(t0)) == 1, 'compressed strongness')
        counts['removed_singletons'] += len(removed)


def literal_covers(a, O):
    n = len(a)
    L = K = 0
    pairs = Counter()
    for word in itertools.permutations(range(n)):
        for cut in range(n):
            p, q = word[:cut], word[cut:]
            if not all(a[x] >> y & 1 for part in (p, q) for x, y in zip(part, part[1:])):
                continue
            if q[0] in O:
                L += 1
                if p and p[-1] in O:
                    K += 1
            if p:
                pairs[p[-1], q[0]] += 1
    return L, K, pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--class-files', nargs='*', default=[])
    args = parser.parse_args()
    report = {'status': 'PASS', 'claim': 'FINITE-EXACT controls; all-order proof is in PROOF.md', 'labeled': []}
    for n in range(1, 7):
        counts = Counter()
        edges = list(itertools.combinations(range(n), 2))
        for code in range(1 << len(edges)):
            a = [0] * n
            for bit, (i, j) in enumerate(edges):
                if code >> bit & 1:
                    a[i] |= 1 << j
                else:
                    a[j] |= 1 << i
            inspect(a, counts, exhaustive_compression=True)
            counts['tournaments'] += 1
            if n <= 4:
                h, e, starts = tables(a)
                for omask in range(1 << n):
                    O = {i for i in range(n) if omask >> i & 1}
                    literal_L, literal_K, _ = literal_covers(a, O)
                    L, K = covers(h, e, starts, (1 << n) - 1, omask)
                    need((L, K) == (literal_L, literal_K), 'literal covers')
                    counts['literal_cover_checks'] += 1
        report['labeled'].append({'order': n, **dict(sorted(counts.items()))})
    controls = Counter()
    for ell in range(1, 11):
        a = chain(ell)
        h, e, starts = tables(a)
        for i in range(ell):
            for j in range(i + 1, ell):
                _, K = covers(h, e, starts, (1 << ell) - 1, (1 << i) | (1 << j))
                expected = 2 if j == i + 1 else (1 << (j - i - 1))
                need(K == expected, 'chain capacity')
                controls['chain_pair_weights'] += 1
        for O, expected in [({0, 1}, 3 * (1 << (ell - 2)) - 2) if ell >= 2 else (set(), None),
                            ({0, 2}, 5 * (1 << (ell - 3)) - 2) if ell >= 3 else (set(), None)]:
            if not O:
                continue
            t, mark, _ = sole_exit_extension(a, sum(1 << i for i in O))
            ht, et, _ = tables(t)
            need(et[-1][mark] == expected and ht[-1] == 2 * expected + 3, 'chain ray')
            inspect(t, controls, exhaustive_compression=True)
            controls['chain_ray_constructions'] += 1
    for p in range(5):
        for q in range(5):
            a = triangle_chain(p, q)
            t, mark, _ = sole_exit_extension(a, 1)
            ht, et, _ = tables(t)
            expected = 3 * (1 << q) if p == 0 else 3 * (1 << (p + q + 1))
            need(et[-1][mark] == expected and ht[-1] == 2 * expected + 3, 'triangle ray')
            inspect(t, controls, exhaustive_compression=True)
            controls['triangle_ray_constructions'] += 1
    rng = random.Random(2142607)
    for order in range(2, 9):
        for _ in range(30):
            u = [0] * order
            for i in range(order):
                for j in range(i + 1, order):
                    if rng.randrange(2):
                        u[i] |= 1 << j
                    else:
                        u[j] |= 1 << i
            omask = rng.randrange(1, 1 << order)
            t, _, _ = sole_exit_extension(u, omask)
            inspect(t, controls, exhaustive_compression=True)
            controls['random_attachment_tournaments'] += 1
    # A real H=13 tournament; this unmarked count is not excluded by the theorem.
    hostile = [6, 28, 8, 17, 5]
    validate(hostile)
    h, e, _ = tables(hostile)
    need(h[-1] == 13 and not any(hostile[v].bit_count() == 1 and e[-1][v] == 5 for v in range(5)), 'H13 hostile')
    need(not in_defect3_spectrum(5), 'forbidden marked register c=5')
    report['controls'] = dict(sorted(controls.items()))
    report['unmarked_H13_hostile'] = {'rows': hostile, 'H': h[-1], 'endpoints': e[-1]}
    report['defect3_first_30_endpoint_values'] = [c for c in range(1, 100000) if in_defect3_spectrum(c)][:30]
    if args.class_files:
        report['additional_class_inputs'] = []
        for name in args.class_files:
            path = Path(name)
            counts = Counter()
            n = None
            for line in path.read_text().splitlines():
                if not line.strip():
                    continue
                a = list(map(int, line.split()))
                n = len(a)
                inspect(a, counts, exhaustive_compression=True)
                counts['tournaments'] += 1
            report['additional_class_inputs'].append({'filename': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                                                     'order': n, **dict(sorted(counts.items()))})
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
