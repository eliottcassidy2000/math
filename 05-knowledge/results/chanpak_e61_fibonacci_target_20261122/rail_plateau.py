"""Explore every count-preserving flip from the strong three-port lower witness, up to isomorphism."""

import argparse
from collections import Counter, deque
from itertools import combinations, permutations, product
from math import comb
from rail_arcs import LOWER, LPORTS, arc_changes, inflate, two_slot_kernel


def canonical(rows):
    """Refine by oriented neighbor counts and individualize any surviving color ties."""
    n = len(rows)
    partition = [tuple(range(n))]
    while True:
        masks = [sum(1 << i for i in cell) for cell in partition]
        refined = []
        for cell in partition:
            colors = {}
            for v in cell:
                key = tuple((rows[v] & mask).bit_count() for mask in masks)
                colors.setdefault(key, []).append(v)
            refined.extend(tuple(colors[key]) for key in sorted(colors))
        if len(refined) == len(partition):
            break
        partition = refined
    def code(order):
        return tuple(sum((rows[i] >> j & 1) << a for a, j in enumerate(order)) for i in order)
    return min(code(tuple(v for cell in cells for v in cell))
               for cells in product(*(permutations(cell) for cell in partition)))


def transitive_blocks(rows):
    """Collapse chains of vertices with the same orientation toward all other vertices."""
    n = len(rows)
    groups = [{i} for i in range(n)]
    for i, j in combinations(range(n), 2):
        if (rows[i] ^ rows[j]) & ~((1 << i) | (1 << j)) == 0:
            group = groups[i] | groups[j]
            for v in group:
                groups[v] = group
    blocks = []
    seen = set()
    for group in groups:
        if min(group) in seen:
            continue
        seen.update(group)
        mask = sum(1 << v for v in group)
        block = sorted(group, key=lambda v: -(rows[v] & mask).bit_count())
        assert [(rows[v] & mask).bit_count() for v in block] == list(reversed(range(len(block))))
        assert len({rows[v] & ~mask for v in group}) == 1
        blocks.append(block)
    quotient = tuple(sum((rows[b[0]] >> c[0] & 1) << j for j, c in enumerate(blocks)) for b in blocks)
    ports = tuple(i for i, b in enumerate(blocks) if len(b) > 1)
    return quotient, ports, tuple(len(blocks[i]) for i in ports), blocks


def component(exponent):
    widths = (exponent, (exponent-3)//2, (exponent-3)//2)
    initial = tuple(inflate(LOWER, LPORTS, widths))
    todo = deque([initial])
    seen = {canonical(initial): 0}
    summaries = []
    while todo:
        graph = todo.popleft()
        quotient, ports, sizes, blocks = transitive_blocks(graph)
        kernel = two_slot_kernel(quotient, ports)
        neutral, neighbors, positives, negatives = 0, Counter(), Counter(), Counter()
        for ((i, r), (j, s)), change in arc_changes(quotient, ports, sizes, kernel):
            if change > 0:
                positives[change] += 1
            elif change < 0:
                negatives[change] += 1
            else:
                neutral += 1
                if i == j and abs(r-s) == 1:
                    continue  # Transposing these consecutive twins is an explicit isomorphism.
                u, v = blocks[i][r], blocks[j][s]
                following = list(graph)
                following[u] ^= 1 << v
                following[v] ^= 1 << u
                key = canonical(following)
                if key not in seen:
                    seen[key] = len(seen)
                    todo.append(tuple(following))
                neighbors[seen[key]] += 1
        assert neutral + sum(positives.values()) + sum(negatives.values()) == comb(len(graph), 2)
        record = (neutral, sorted(neighbors.items()), min(positives), max(negatives), positives[2])
        print('type', len(summaries), 'neutral, other neutral arcs by type, min positive, max negative, +2 arcs:', record, flush=True)
        summaries.append(record)
    return summaries


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exponent', type=int, default=61)
    args = parser.parse_args()
    if args.exponent < 5 or args.exponent % 2 == 0:
        parser.error('choose an odd exponent at least five')
    component(args.exponent)
