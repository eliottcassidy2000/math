#!/usr/bin/env python3
"""Second adversarial audit: no producer imports, explicit hostile control."""
import argparse
import itertools
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def act(s, letter):
    D, E, U, d = s
    if letter == 'O':
        return D, 0, 2*U, d
    if letter == 'Z':
        return D+U-E, d, 2*U-d, d
    return D+U-E, D-E, D+U, D


def literal_profile(a, r):
    n = len(a)
    full = (1 << n)-1
    counts = [0]*(1 << n)
    counts[0] = 1
    starts = [0]*(1 << n)
    endpoint = 0
    for size in range(1, n+1):
        for path in itertools.permutations(range(n), size):
            if all(a[u] >> v & 1 for u, v in zip(path, path[1:])):
                mask = sum(1 << v for v in path)
                counts[mask] += 1
                if path[0] == r:
                    starts[mask] += 1
                if size == n and path[-1] == r:
                    endpoint += 1
    return counts[full], endpoint, sum(starts[s]*counts[full^s] for s in range(1, full+1)), counts[full^(1 << r)]


def check_literal_letters():
    checks = 0
    for n in range(1, 5):
        edges = list(itertools.combinations(range(n), 2))
        for code in range(1 << len(edges)):
            a = [0]*n
            for k, (u, v) in enumerate(edges):
                if code >> k & 1:
                    a[u] |= 1 << v
                else:
                    a[v] |= 1 << u
            for r in range(n):
                source = literal_profile(a, r)
                for c in 'OZJ':
                    b = list(a)+[0 if c == 'O' else 1 << r]
                    for v in range(n):
                        if c == 'O' or v != r:
                            b[v] |= 1 << n
                    dest = literal_profile(b, n if c == 'J' else r)
                    need(act(source, c) == dest, 'independent literal physical letter')
                    checks += 1
    return checks


def sharp_core(k):
    if k == 1:
        return [2, 4, 1]
    n = k+1
    # A transitive tournament with its long edge reversed is strong.
    a = [sum(1 << v for v in range(u+1, n)) for u in range(n)]
    a[0] ^= 1 << (n-1)
    a[n-1] |= 1
    owners = [u for u in range(n) for _ in range(k-a[u].bit_count())]
    for u in owners:
        v = len(a)
        a.append(((1 << n)-1) ^ (1 << u))
        a[u] |= 1 << v
        for old in range(n, v):
            a[old] |= 1 << v
    return a


def peel(a, k, omitted):
    mask = ((1 << len(a))-1) ^ (1 << omitted)
    while mask:
        found = next((v for v in range(len(a)) if mask >> v & 1 and (a[v]&mask).bit_count() < k), None)
        if found is None:
            return False
        mask ^= 1 << found
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fail', action='store_true')
    args = ap.parse_args()
    need(not args.fail, 'intentional failure control')
    physical_checks = check_literal_letters()
    todo = [((1, 0, 0, 0), False, 0)]
    checked = 0
    while todo:
        s, has_j, depth = todo.pop()
        need((s == (1, 0, 0, 0)) == (not has_j), 'pure-count iff no J')
        if has_j:
            need(s[3] > 0, 'positive deletion cannot be erased')
            need(s[0] > 0 and 0 <= s[1] <= s[0] <= s[2], 'physical singleton cone')
        checked += 1
        if depth < 10:
            todo.extend((act(s, c), has_j or c == 'J', depth+1) for c in 'OZJ')
    rows = []
    for k in range(1, 21):
        a = sharp_core(k)
        need(len(a) == (k+1)*(k+2)//2, 'sharp order')
        need(min(x.bit_count() for x in a) >= k, 'minimum degree')
        need(all(peel(a, k, v) for v in range(len(a))), 'minimality by every vertex deletion')
        rows.append([k, len(a)])
    a = [218, 252, 57, 240, 224, 193, 132, 4]
    starts = [0]*8
    ends = [0]*8
    for p in itertools.permutations(range(8)):
        if all(a[u] >> v & 1 for u, v in zip(p, p[1:])):
            starts[p[0]] += 1
            ends[p[-1]] += 1
    need(starts == [23, 72, 5, 27, 17, 14, 10, 5], 'literal 173 start profile')
    need(sum(starts) == 173 and 1 not in starts, '173 is not a unit-start witness')
    need(ends[2] == 9, '173 marked endpoint')
    print(json.dumps({'status': 'PASS', 'pure_count_words_through_length_10': checked, 'literal_physical_letter_checks': physical_checks,
        'sharp_cores_k_order': rows, 'literal_173_masks': a, 'literal_173_starts': starts,
        'method': 'independent arithmetic, direct deletion peeling, literal full permutations; no producer imports'}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
