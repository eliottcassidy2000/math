"""Independent boundary counts via independent colour classes, not vertex search."""
from functools import lru_cache
from itertools import combinations, product, permutations
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def chromatic_evaluation(n, edges, q):
    rows = [0] * n
    for u, v in edges:
        rows[u] |= 1 << v
        rows[v] |= 1 << u
    independent = [True] * (1 << n)
    for mask in range(1, 1 << n):
        bit = mask & -mask
        rest = mask ^ bit
        independent[mask] = independent[rest] and not (rows[bit.bit_length() - 1] & rest)
    @lru_cache(None)
    def count(mask, colors):
        if not mask:
            return 1
        if not colors:
            return 0
        ans, sub = 0, mask
        while True:
            if independent[sub]:
                ans += count(mask ^ sub, colors - 1)
            if not sub:
                return ans
            sub = (sub - 1) & mask
    return count((1 << n) - 1, q)


def boundary(n, edges, ports, partition, q=6):
    blocks = max(partition) + 1
    labels = {v: partition[i] for i, v in enumerate(ports)}
    for v in range(n):
        if v not in labels:
            labels[v] = len(set(labels.values()))
    quotient = {tuple(sorted((labels[u], labels[v]))) for u, v in edges}
    if any(u == v for u, v in quotient):
        return 0
    quotient.update(combinations(range(blocks), 2))
    total = chromatic_evaluation(len(set(labels.values())), quotient, q)
    orbit = 1
    for i in range(blocks):
        orbit *= q - i
    need(total % orbit == 0, 'fixed-colour orbit divisibility')
    return total // orbit


def main():
    partitions = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (0, 1, 2)]
    base = list(combinations(range(3), 2)) + [(i, i + 3) for i in range(3)]
    base_weights = [boundary(6, base, (3, 4, 5), p) for p in partitions]
    need(base_weights == [60, 68, 68, 68, 71], 'unpinning base weights')
    k4_join = list(combinations(range(3, 7), 2)) + [(i, j) for i in range(3) for j in range(3, 7)]
    weights5 = [boundary(7, k4_join, (0, 1, 2), p) for p in partitions]
    need(weights5 == [120, 24, 24, 24, 0], 'sharp chromatic-five control')
    nae = list(combinations(range(3, 6), 2)) + [(i, i + 3) for i in range(3)]
    nae += list(combinations(range(6, 9), 2)) + [(i, j) for i in range(6) for j in range(6, 9)]
    weights6 = [boundary(9, nae, (0, 1, 2), p) for p in partitions]
    need(weights6 == [0, 48, 48, 48, 12], 'sharp chromatic-six NAE control')
    separate = list(combinations(range(3), 2)) + list(combinations(range(3, 6), 2))
    need(chromatic_evaluation(7, separate, 6) == 86400, 'actual six-colouring count')
    need([boundary(7, separate, (0, 3, 6), p) for p in partitions] == [400] * 5,
         'all full boundary weights, including coincidences')
    success = sum(len(set(q[:3])) > 1 and len(set(q[3:6])) > 1 for q in product(range(3), repeat=7))
    need(success == (27 - 3)**2 * 3 == 1728, 'all ternary pair maps')
    even = [0, 3, 5, 6]
    deterministic = []
    for perm in permutations(even):
        f = dict(zip(even, perm))
        # Equality pattern across two frames sees precisely x xor y.
        need(all(f[x] in even for x in even), 'state type')
        differences = {x ^ f[x] for x in even}
        if len(differences) == 1:
            deterministic.append(perm)
    need(len(deterministic) == 4, 'only four deterministic difference kernels')
    print(json.dumps({'status': 'PASS', 'method': 'independent-subset colour-class recurrence + quotient ports',
                      'base_boundary_weights': base_weights, 'chi5_boundary_weights': weights5,
                      'chi6_NAE_boundary_weights': weights6, 'ordinary_six_colorings': 86400,
                      'successful_ternary_maps': success, 'deterministic_four_state_kernels': 4},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
