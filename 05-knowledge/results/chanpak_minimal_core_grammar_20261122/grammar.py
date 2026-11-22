"""Target-capped canonical augmentation of strong tournament isomorphism classes."""
import argparse
from collections import Counter, deque
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / '05-knowledge/results/rank_geometry_metagraph_20261026/enumeration'))
import enumerate as canon


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def strongly_connected(rows, mask):
    if not mask:
        return False
    root = mask & -mask

    def reached(reverse):
        seen = 0
        frontier = root
        while frontier:
            take = frontier & -frontier
            frontier ^= take
            seen |= take
            i = take.bit_length() - 1
            out = ~rows[i] if reverse else rows[i]
            frontier |= mask & out & ~seen
        return seen == mask

    return reached(False) and reached(True)


def endpoint_cut(rows):
    n = len(rows)
    full = (1 << n) - 1
    dp = [[[0] * n for _ in range(n)] for _ in range(full + 1)]
    paths = [1] + [0] * full
    for t in range(n):
        dp[1 << t][t][t] = 1
    for s in range(1, full + 1):
        paths[s] = sum(sum(row) for row in dp[s])
        for a in bits(s):
            for b in bits(s):
                count = dp[s][a][b]
                if count:
                    for c in bits(rows[b] & ~s):
                        dp[s | (1 << c)][a][c] += count
    weights = [[0] * n for _ in range(n)]
    for s in range(1, full + 1):
        if s.bit_count() % 2 == 0:
            for i in bits(s):
                for j in bits(s):
                    weights[i][j] += paths[full ^ s] * dp[s][i][j]
    return paths[-1], weights


def cut(weights, mask):
    full = (1 << len(weights)) - 1
    return sum(weights[i][j] for i in bits(mask) for j in bits(full ^ mask))


def minimum_completion(weights, yes, no):
    """Directed integer min cut with each assigned side contracted to a terminal."""
    n = len(weights)
    assert yes and no and not (yes & no)
    vertices = [0 if (yes >> i) & 1 else 1 if (no >> i) & 1 else i + 2
                for i in range(n)]
    residual = [[0] * (n + 2) for _ in range(n + 2)]
    for i, row in enumerate(weights):
        for j, weight in enumerate(row):
            if vertices[i] != vertices[j]:
                residual[vertices[i]][vertices[j]] += weight
    flow = 0
    while True:
        parent = [-1] * (n + 2)
        parent[0] = 0
        queue = deque([0])
        while queue and parent[1] < 0:
            i = queue.popleft()
            for j, capacity in enumerate(residual[i]):
                if capacity and parent[j] < 0:
                    parent[j] = i
                    queue.append(j)
        if parent[1] < 0:
            return flow
        arcs = []
        j = 1
        while j:
            i = parent[j]
            arcs.append((i, j))
            j = i
        increment = min(residual[i][j] for i, j in arcs)
        flow += increment
        for i, j in arcs:
            residual[i][j] -= increment
            residual[j][i] += increment


def capped_masks(weights, budget):
    """Yield each nonconstant contact mask with cut weight at most budget."""
    n = len(weights)
    def visit(yes, no, next_bit):
        if next_bit == n:
            value = cut(weights, yes)
            if value <= budget:
                yield yes, value
            return
        if minimum_completion(weights, yes, no) > budget:
            return
        bit = 1 << next_bit
        yield from visit(yes | bit, no, next_bit + 1)
        yield from visit(yes, no | bit, next_bit + 1)

    for change in range(1, n):
        prefix = (1 << change) - 1
        bit = 1 << change
        yield from visit(prefix, bit, change + 1)
        yield from visit(bit, prefix, change + 1)


def busch_floor(n):
    if n < 3:
        raise ValueError('nontrivial strong core requires at least three vertices')
    k, digit = divmod(n - 3, 3)
    return 5 ** k * (3, 5, 9)[digit]


def strong_classes(cap, max_order=None):
    """Generate (order, canonical code, H). The same class is emitted once."""
    if cap < 3 or (max_order is not None and max_order < 3):
        return
    seed = canon.canonical(3, canon.encode((2, 4, 1)))
    level = {seed.code: (3, seed.canonical_automorphisms())}
    n = 3
    while level:
        for code, (h, _) in level.items():
            yield n, code, h
        if (max_order is not None and n >= max_order) or busch_floor(n + 1) > cap:
            break
        following = {}
        for code, (h, group) in level.items():
            parent_h, weights = endpoint_cut(canon.rows(n, code))
            assert parent_h == h
            for mask, gain in capped_masks(weights, (cap - h) // 2):
                if any(canon.subset_image(mask, pi) < mask for pi in group):
                    continue
                child = canon.extend(n, code, ((1 << n) - 1) ^ mask)
                cf = canon.canonical(n + 1, child)
                rows = canon.rows(n + 1, child)
                full = (1 << (n + 1)) - 1
                index = next(i for i, v in enumerate(cf.orders[0])
                             if strongly_connected(rows, full ^ (1 << v)))
                if not any(order[index] == n for order in cf.orders):
                    continue
                if cf.code in following:
                    raise ValueError('canonical extension accepted twice')
                following[cf.code] = (h + 2 * gain, cf.canonical_automorphisms())
        level = following
        n += 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('target', type=int)
    args = parser.parse_args()
    if args.target >= 3 and args.target % 2:
        result = Counter()
        for size, code, count in strong_classes(args.target):
            if count == args.target:
                result[size] += 1
                print('order', size, 'rows', canon.rows(size, code))
        print('strong isomorphism multiplicity by order', sorted(result.items()))
    else:
        print('strong isomorphism multiplicity by order [] (orders at least three)')
