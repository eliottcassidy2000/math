#!/usr/bin/env python3
"""Exhaust all Bell(8)=4140 candidate partitions of the n=4 cube."""
import json
from pathlib import Path
def require(ok, message):
    if not ok:
        raise RuntimeError(message)

N = 8
E = {(0, 1), (0, 4), (1, 2), (1, 5), (2, 4), (4, 5)}
merged_E = E | {(3, 6)}

def partitions(n):
    a = [0] * n

    def rec(i, m):
        if i == n:
            yield tuple(a)
            return
        for b in range(m + 2):
            a[i] = b
            yield from rec(i + 1, max(m, b))
    yield from rec(1, 0)

def passes(labels, edges):
    close = all((labels[x] == labels[y] for x, y in edges))
    neigh = [set() for _ in labels]
    for x, y in edges:
        neigh[x].add(y)
        neigh[y].add(x)
    connected = True
    for block in set(labels):
        vertices = {i for i, a in enumerate(labels) if a == block}
        seen = {min(vertices)}
        todo = list(seen)
        while todo:
            for y in neigh[todo.pop()] & vertices:
                if y not in seen:
                    seen.add(y)
                    todo.append(y)
        connected = connected and seen == vertices
    return (connected, close)
rows = []
for name, edges in [('ordinary', E), ('converse_merged', merged_E)]:
    counts = {f'{a}{b}': 0 for a in (0, 1) for b in (0, 1)}
    positive = []
    total = 0
    for p in partitions(N):
        a, b = passes(p, edges)
        counts[f'{int(a)}{int(b)}'] += 1
        total += 1
        if a and b:
            positive.append(p)
    require(total == 4140 and len(positive) == 1, 'verification failed')
    require(passes(tuple(range(N)), edges) == (True, False), 'verification failed')
    require(passes((0,) * N, edges) == (False, True), 'verification failed')
    rows.append({'mode': name, 'partitions_tested': total, 'condition_pair_counts_connected_closed': counts, 'unique_successful_partition': positive[0], 'split_hostile_singletons_connected_not_closed': True, 'merge_hostile_whole_cube_closed_not_connected': True})
Path('recognition_hostiles.json').write_text(json.dumps({'status': 'PASS', 'n': 4, 'bit_coordinates': [[0, 2], [0, 3], [1, 3]], 'rows': rows}, indent=2) + '\n')
print(Path('recognition_hostiles.json').read_text())
