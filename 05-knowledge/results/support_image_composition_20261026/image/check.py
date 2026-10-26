#!/usr/bin/env python3
"""Exact partial-support census, defect coefficients and joint line controls."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
import json

CHECKS = 0


def need(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(message)


def rows_from_code(n, code, radix):
    rows = [0] * n
    for x, y in combinations(range(n), 2):
        code, bit = divmod(code, radix)
        if radix == 2:
            if bit:
                rows[x] |= 1 << y
            else:
                rows[y] |= 1 << x
        elif bit == 1:
            rows[x] |= 1 << y
        elif bit == 2:
            rows[y] |= 1 << x
    return tuple(rows)


def closure(rows):
    reach = [row | (1 << x) for x, row in enumerate(rows)]
    for k in range(len(rows)):
        for x in range(len(rows)):
            if reach[x] >> k & 1:
                reach[x] |= reach[k]
    return reach


def decode(rows):
    n = len(rows)
    if any(row >> n or row >> x & 1 for x, row in enumerate(rows)):
        return None
    if any(rows[x] >> y & 1 and rows[y] >> x & 1 for x, y in combinations(range(n), 2)):
        return None
    reach = closure(rows)
    if any(not (reach[x] >> y & 1 or reach[y] >> x & 1)
           for x, y in combinations(range(n), 2)):
        return None
    remaining, blocks = set(range(n)), []
    while remaining:
        x = min(remaining)
        block = {y for y in remaining if reach[x] >> y & 1 and reach[y] >> x & 1}
        blocks.append(block)
        remaining -= block
    blocks.sort(key=lambda b: -reach[min(b)].bit_count())
    owner = {v: i for i, block in enumerate(blocks) for v in block}
    result = list(rows)
    for i, block in enumerate(blocks):
        for j in range(i + 1, len(blocks)):
            target = sum(1 << v for v in blocks[j])
            for u in block:
                if rows[u] & target != (target if j == i + 1 else 0):
                    return None
                result[u] |= target
    for x, y in combinations(range(n), 2):
        if owner[x] != owner[y] or rows[x] >> y & 1 or rows[y] >> x & 1:
            continue
        xy = sum(bool(rows[x] >> w & 1 and rows[w] >> y & 1)
                 for w in blocks[owner[x]])
        yx = sum(bool(rows[y] >> w & 1 and rows[w] >> x & 1)
                 for w in blocks[owner[x]])
        if (xy >= 2) == (yx >= 2):
            return None
        if xy >= 2:
            result[x] |= 1 << y
        else:
            result[y] |= 1 << x
    return tuple(result)


@lru_cache(maxsize=None)
def paths(rows):
    n = len(rows)
    result = []

    def visit(word, unseen):
        if not unseen:
            result.append(word)
            return
        choices = unseen if not word else unseen & rows[word[-1]]
        while choices:
            bit = choices & -choices
            choices ^= bit
            visit(word + (bit.bit_length() - 1,), unseen ^ bit)

    visit((), (1 << n) - 1)
    return tuple(result)


def support(rows):
    result = [0] * len(rows)
    for word in paths(rows):
        for x, y in zip(word, word[1:]):
            result[x] |= 1 << y
    return tuple(result)


def polynomial(T, D):
    n = len(T)
    if n == 0:
        return (1,)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = [1] + [0] * (n - 1)
    for mask in range(1, 1 << n):
        for v in range(n):
            a = dp.get((mask, v))
            if a is None:
                continue
            for w in range(n):
                if mask >> w & 1 or not T[v] >> w & 1:
                    continue
                shift = int(not D[v] >> w & 1)
                b = dp.setdefault((mask | (1 << w), w), [0] * n)
                for i in range(n - shift):
                    b[i + shift] += a[i]
    return tuple(sum(dp.get(((1 << n) - 1, v), [0] * n)[i]
                     for v in range(n)) for i in range(n))


def opposite(T):
    n = len(T)
    return tuple(((1 << n) - 1) ^ (1 << v) ^ row for v, row in enumerate(T))


def partner(T, P):
    result = list(opposite(T))
    for x, y in zip(P, P[1:]):
        result[x] |= 1 << y
        result[y] &= ~(1 << x)
    return tuple(result)


def blue(T, P):
    reverse = dict(zip(P, reversed(P)))
    return all(bool(T[reverse[x]] >> reverse[y] & 1) != bool(T[x] >> y & 1)
               for x, y in combinations(range(len(T)), 2))


def agreement_path(T, S):
    n = len(T)
    common = tuple(a & b for a, b in zip(T, S))
    if sum(row.bit_count() for row in common) != n - 1:
        return None
    indeg = [sum(row >> v & 1 for row in common) for v in range(n)]
    if any(d > 1 for d in indeg) or any(row.bit_count() > 1 for row in common):
        return None
    start = [v for v in range(n) if indeg[v] == 0]
    if len(start) != 1:
        return None
    word, seen = [], set()
    v = start[0]
    while v not in seen:
        word.append(v)
        seen.add(v)
        if common[v] == 0:
            break
        v = common[v].bit_length() - 1
    return tuple(word) if len(word) == n else None


def relabel(T, p):
    result = [0] * len(T)
    for x, y in combinations(range(len(T)), 2):
        if T[x] >> y & 1:
            result[p[x]] |= 1 << p[y]
        else:
            result[p[y]] |= 1 << p[x]
    return tuple(result)


@lru_cache(maxsize=None)
def canonical(T):
    return min(relabel(T, p) for p in permutations(range(len(T))))


def pad_join(T):
    n = len(T)
    return (((1 << (n + 2)) - 2),
            *(row << 1 | (1 << (n + 1)) for row in T), 0)


def sandwich(h, s, t):
    n = len(h)
    result = [1 << (s + 1)]
    for v, row in enumerate(h):
        result.append((row << 1) | (int(v != s)) | ((1 << (n + 1)) if v == t else 0))
    result.append(1 | sum(1 << (v + 1) for v in range(n) if v != t))
    return tuple(result)


def main():
    census = []
    first = None
    for n in range(6):
        images = {support(rows_from_code(n, code, 2)) for code in range(1 << (n * (n - 1) // 2))}
        counts = Counter(total=3 ** (n * (n - 1) // 2), images=len(images))
        for code in range(counts['total']):
            D = rows_from_code(n, code, 3)
            T = decode(D)
            if T is None:
                need(D not in images, 'decoder rejected actual image')
                continue
            counts['decoded'] += 1
            if support(D) != D:
                need(D not in images, 'unsaturated actual image')
                continue
            counts['saturated'] += 1
            coeffs = polynomial(T, D)
            need(coeffs[0] == len(paths(D)) and sum(coeffs) == len(paths(T)), 'DP endpoints')
            by_cover = Counter()
            for P in paths(T):
                gaps = [i for i, (x, y) in enumerate(zip(P, P[1:])) if not D[x] >> y & 1]
                boundaries = [-1] + gaps + [n - 1]
                blocks = [P[a + 1:b + 1] for a, b in zip(boundaries, boundaries[1:])]
                need(all(block and all(D[x] >> y & 1 for x, y in zip(block, block[1:]))
                         for block in blocks) if n else P == (), 'ordered path-cover decomposition')
                by_cover[len(gaps)] += 1
            need(tuple(by_cover[i] for i in range(max(1, n))) == coeffs, 'path-cover coefficients')
            accepted = sum(coeffs[1:]) == 0
            need(accepted == (D in images), 'exact image criterion')
            if accepted:
                counts['accepted'] += 1
            else:
                counts['false_local'] += 1
                if first is None:
                    first = {'n': n, 'code': code, 'D': D, 'completion': T, 'coefficients': coeffs}
        census.append({'n': n, **dict(counts)})
    need(first['n'] == 5 and first['code'] == 12295, 'first local false image')
    need(first['coefficients'] == (7, 2, 0, 0, 0), 'first defect polynomial')
    need(census[-1]['false_local'] == 480 and census[-1]['accepted'] == 1024, 'five-vertex image counts')

    cycle = (2, 4, 8, 1)
    completed = (6, 12, 8, 1)
    need(support(cycle) == cycle and polynomial(completed, cycle) == (4, 0, 1, 0), 'joint-only gap hostile')
    for x, y in [(0, 2), (1, 3)]:
        added = list(cycle)
        added[x] |= 1 << y
        need(paths(tuple(added)) == paths(cycle), 'single addition is inert')
    need(set(paths(completed)) - set(paths(cycle)) == {(1, 3, 0, 2)}, 'joint new word')
    transitive = (6, 4, 0)
    need(support(transitive) == (2, 4, 0), 'transitive triple')
    need(support((2, 0)) == (2, 0) and support((0, 0)) == (0, 0)
         and decode((0, 0)) is None, 'nonhereditary induced deletion')
    need(support(completed) == completed, 'insertion reveals old hidden arc')

    line_census = []
    for n in range(3, 6):
        all_T = [rows_from_code(n, k, 2) for k in range(1 << (n * (n - 1) // 2))]
        reps = sorted({canonical(T) for T in all_T})
        aut = {S: sum(relabel(S, p) == S for p in permutations(range(n))) for S in reps}
        incidence = {}
        hub_tests = 0
        for T in reps:
            literal = Counter()
            for P in paths(T):
                S = partner(T, P)
                need(agreement_path(T, S) == P, 'line agreement recovers word')
                need(blue(T, P) == blue(S, P), 'colour preserved by complement')
                literal[canonical(S), blue(T, P)] += 1
                padded = (0,) + tuple(v + 1 for v in P) + (n + 1,)
                hub = partner(pad_join(S), padded)
                need(hub == sandwich(T, P[0], P[-1]), 'same joint ports give same hub')
                need(blue(pad_join(S), padded) == blue(T, P), 'padding preserves colour')
                need(all(r == (1 << (n + 2)) - 1 for r in closure(hub)), 'hub strong')
                hub_tests += 1
            for S in reps:
                transport = Counter()
                for p in permutations(range(n)):
                    target = relabel(S, p)
                    P = agreement_path(T, target)
                    if P is not None:
                        need(partner(T, P) == target, 'agreement sufficient')
                        transport[blue(T, P)] += 1
                for colour in (False, True):
                    need(transport[colour] == aut[S] * literal[S, colour], 'transport multiplicity formula')
                incidence[T, S] = literal[S, False]
        for T in reps:
            for S in reps:
                need(incidence[T, S] * aut[S] == incidence[S, T] * aut[T], 'normalized black incidence symmetry')
        line_census.append({'n': n, 'ordinary_types': len(reps), 'hub_path_tests': hub_tests})
    print(json.dumps({'status': 'FINITE-EXACT PASS', 'checks': CHECKS,
                      'universe': 'All loopless oriented graphs n=0..5; all ordinary tournament-type pairs n=3..5 with all transports; every representative Hamilton word for hub tests.',
                      'partial_graph_census': census, 'first_false_image': first,
                      'joint_only_control': {'D': cycle, 'T': completed, 'polynomial': [4, 0, 1, 0], 'new_word': [1, 3, 0, 2]},
                      'line_census': line_census}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
