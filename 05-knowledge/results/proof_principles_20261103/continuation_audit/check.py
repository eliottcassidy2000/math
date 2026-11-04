"""Independent literal-permutation audit; imports no producer code."""
import hashlib
import json
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path

SOURCE = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


@lru_cache(None)
def h(rows):
    return sum(all((rows[u] >> v) & 1 for u, v in zip(word, word[1:]))
               for word in permutations(range(len(rows))))


def reverse(rows):
    return tuple(((1 << len(rows)) - 1) ^ (1 << u) ^ row for u, row in enumerate(rows))


def extend(rows, bits):
    n = len(rows)
    return tuple(row | (0 if bits & (1 << u) else 1 << n) for u, row in enumerate(rows)) + (bits,)


def response(rows):
    return tuple(h(extend(rows, mask)) for mask in range(1 << len(rows)))


def tournament(n, mask):
    rows = [0] * n
    for i, (u, v) in enumerate(combinations(range(n), 2)):
        a, b = (u, v) if mask & (1 << i) else (v, u)
        rows[a] |= 1 << b
    return tuple(rows)


def substitute(quotient, hole, module):
    blocks = []
    offset = 0
    for u in range(len(quotient)):
        block = list(range(offset, offset + (len(module) if u == hole else 1)))
        blocks.append(block)
        offset += len(block)
    out = [0] * offset
    for u, block in enumerate(blocks):
        for v, other in enumerate(blocks):
            if quotient[u] & (1 << v):
                for a in block:
                    for b in other:
                        out[a] |= 1 << b
    for i, a in enumerate(blocks[hole]):
        for j, b in enumerate(blocks[hole]):
            if module[i] & (1 << j):
                out[a] |= 1 << b
    return tuple(out)


def weights(resp, n):
    values = {}
    for u, v in combinations(range(n), 2):
        twice = resp[1 << u] + resp[1 << v] - resp[(1 << u) | (1 << v)] - resp[0]
        require(twice > 0 and twice % 2 == 0, 'positive integral capacity')
        values[u, v] = twice // 2
    return values


def main():
    module = (2, 4, 1)
    contexts = probes = 0
    for n in range(1, 4):
        for code in range(1 << (n * (n - 1) // 2)):
            q = tournament(n, code)
            for hole in range(n):
                a = substitute(q, hole, module)
                b = substitute(q, hole, reverse(module))
                ra, rb = response(a), response(b)
                require(ra == rb and h(a) == h(b), 'contextual module reversal')
                contexts += 1
                probes += len(ra)
    a, b = (10, 12, 9, 16, 7), (12, 9, 10, 16, 7)
    require(h(a) == h(b) == 15, 'strong base')
    require(response(a) == response(b), 'strong one-probe equality')
    ax, bx = extend(a, 1), extend(b, 1)
    require(h(ax) == h(bx) == 31, 'strong common intermediate')
    require((h(extend(ax, 2)), h(extend(bx, 2))) == (81, 87), 'two probes separate')
    m = (22, 88, 50, 69, 104, 11, 37)
    require(all(row.bit_count() == 3 for row in m), 'regularity')
    words = [word for word in permutations(range(7))
             if all((m[u] >> v) & 1 for u, v in zip(word, word[1:]))]
    require((len(words), sum(word[0] == 0 for word in words), sum(word[-1] == 0 for word in words)) == (171, 24, 25),
            'endpoint imbalance')
    require((h(extend(m, 1)), h(extend(reverse(m), 1))) == (403, 407), 'regularity is insufficient')
    collision_pairs = 0
    for n in range(1, 5):
        groups = defaultdict(list)
        for code in range(1 << (n * (n - 1) // 2)):
            t = tournament(n, code)
            groups[response(t)].append(t)
        for resp, group in groups.items():
            w = weights(resp, n)
            for t, u in combinations(group, 2):
                divergences = [0] * n
                changed = 0
                for (i, j), value in w.items():
                    if bool(t[i] & (1 << j)) != bool(u[i] & (1 << j)):
                        a, b = (i, j) if t[i] & (1 << j) else (j, i)
                        divergences[a] += value
                        divergences[b] -= value
                        changed += 1
                require(not any(divergences) and changed >= 3, 'collision circulation')
                collision_pairs += 1
    a, b = (24, 1, 3, 6, 14), (10, 16, 3, 6, 13)
    wa, wb = weights(response(a), 5), weights(response(b), 5)
    require([wa[e] for e in ((0, 1), (0, 4), (1, 4))] == [4, 4, 4], 'old triangle balanced')
    require((h(a), h(b)) == (11, 15) and wa != wb, 'endogenous weight hostile')
    require(hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == SOURCE, 'source changed during run')
    print(json.dumps({'status': 'PASS', 'source_sha256': SOURCE, 'checks': checks,
                      'contexts': contexts, 'one_probe_comparisons': probes,
                      'collision_pairs_through_order4': collision_pairs,
                      'strong_two_probe_counts': [81, 87], 'regular_hostile': [403, 407],
                      'balanced_old_weight_hostile': [11, 15]}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
