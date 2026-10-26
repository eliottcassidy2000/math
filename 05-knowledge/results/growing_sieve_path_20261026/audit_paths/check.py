#!/usr/bin/env python3
"""Independent subset-DP support and cut reconstruction; no producer imports."""
import itertools
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def vertices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def tournament(n, code):
    rows = [0] * n
    for k, (u, v) in enumerate(itertools.combinations(range(n), 2)):
        if (code >> k) & 1:
            rows[u] |= 1 << v
        else:
            rows[v] |= 1 << u
    return tuple(rows)


def reverse(rows):
    n = len(rows)
    return tuple(sum(1 << u for u in range(n) if rows[u] >> v & 1)
                 for v in range(n))


def ending_sets(rows):
    n = len(rows)
    incoming = reverse(rows)
    result = [0] * (1 << n)
    for s in range(1, 1 << n):
        if s & (s - 1) == 0:
            result[s] = s
        else:
            for v in vertices(s):
                if result[s ^ (1 << v)] & incoming[v]:
                    result[s] |= 1 << v
    return result


def used_arcs(rows):
    n = len(rows)
    full = (1 << n) - 1
    ends, starts = ending_sets(rows), ending_sets(reverse(rows))
    used = [0] * n
    for s in range(1, full):
        for u in vertices(ends[s]):
            used[u] |= rows[u] & starts[full ^ s]
    return tuple(used)


def literal_used(rows):
    n = len(rows)
    result = [0] * n
    count = 0
    for p in itertools.permutations(range(n)):
        if all(rows[u] >> v & 1 for u, v in zip(p, p[1:])):
            count += 1
            for u, v in zip(p, p[1:]):
                result[u] |= 1 << v
    return tuple(result), count


def reach(rows, allowed):
    result = [(row & allowed) | (1 << u) if allowed >> u & 1 else 0
              for u, row in enumerate(rows)]
    for k in vertices(allowed):
        for u in vertices(allowed):
            if result[u] >> k & 1:
                result[u] |= result[k]
    return tuple(result)


def components(rows, allowed):
    closure = reach(rows, allowed)
    todo, blocks = allowed, []
    while todo:
        u = (todo & -todo).bit_length() - 1
        block = sum(1 << v for v in vertices(todo)
                    if closure[u] >> v & 1 and closure[v] >> u & 1)
        blocks.append(block)
        todo ^= block
    return blocks, closure


def reconstruct(used):
    n = len(used)
    full = (1 << n) - 1
    blocks, closure = components(used, full)
    result = list(used)
    ownership = {v: block for block in blocks for v in vertices(block)}
    cut_reach = {block: [(z, reach(used, block ^ (1 << z)))
                         for z in vertices(block)]
                 for block in blocks if block & (block - 1)}
    witnesses = []
    for x, y in itertools.combinations(range(n), 2):
        if (used[x] >> y & 1) or (used[y] >> x & 1):
            continue
        signs = set()
        if ownership[x] != ownership[y]:
            xy, yx = (closure[x] >> y) & 1, (closure[y] >> x) & 1
            require(xy != yx, 'outer components not totally ordered')
            signs.add(1 if xy else -1)
            witness = 'outer'
        else:
            witness = None
            for z, cut in cut_reach[ownership[x]]:
                if z in (x, y):
                    continue
                xy, yx = (cut[x] >> y) & 1, (cut[y] >> x) & 1
                if xy != yx:
                    signs.add(1 if xy else -1)
                    witness = z
        require(len(signs) == 1, 'missing, inconsistent or ambiguous cut witness')
        sign = signs.pop()
        u, v = (x, y) if sign == 1 else (y, x)
        result[u] |= 1 << v
        witnesses.append([u, v, witness])
    return tuple(result), witnesses


def main():
    census = []
    first_strong_missing = None
    checks = 0
    for n in range(1, 7):
        total = 1 << (n * (n - 1) // 2)
        full = (1 << n) - 1
        signatures = set()
        strong_count = missing_count = missing_arcs = cut_checks = 0
        for code in range(total):
            rows = tournament(n, code)
            used = used_arcs(rows)
            signatures.add(used)
            recovered, witnesses = reconstruct(used)
            require(recovered == rows, ('reconstruction', n, code))
            checks += 1
            blocks, _ = components(rows, full)
            support_blocks, _ = components(used, full)
            require(set(blocks) == set(support_blocks), 'component mismatch')
            checks += 1
            if n <= 5:
                literal, count = literal_used(rows)
                require(literal == used and count > 0, 'DP versus permutation support')
                checks += 1
            if used != rows:
                missing_count += 1
                missing_arcs += len(witnesses)
            if len(blocks) == 1:
                strong_count += 1
                if used != rows and first_strong_missing is None:
                    first_strong_missing = {'n': n, 'code': code,
                        'rows': rows, 'used': used, 'witnesses': witnesses}
                for z in range(n):
                    allowed = full ^ (1 << z)
                    require(reach(rows, allowed) == reach(used, allowed),
                            ('strong deletion reachability', n, code, z))
                    checks += 1
                    cut_checks += 1
        require(len(signatures) == total, 'used support signatures collide')
        checks += 1
        census.append({'n': n, 'labelled_tournaments': total,
            'distinct_used_supports': len(signatures), 'strong': strong_count,
            'with_missing_arcs': missing_count, 'missing_arcs_total': missing_arcs,
            'strong_vertex_deletion_checks': cut_checks})
    transitive = (6, 4, 0)
    used = used_arcs(transitive)
    require(reach(transitive, 5) != reach(used, 5),
            'non-strong scope hostile absent')
    checks += 1
    require(first_strong_missing['n'] == 5, 'first strong unused arc boundary')
    checks += 1
    print(json.dumps({'status': 'PASS', 'checks': checks, 'census': census,
        'independent_method': 'subset endpoint DP versus literal permutations through5',
        'first_strong_unused_arc': first_strong_missing,
        'non_strong_hostile': {'rows': transitive, 'used': used,
                              'deleted_vertex': 1}}, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
