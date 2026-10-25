#!/usr/bin/env python3
"""Explicit endpoint-four hostile; not an exhaustive order-seven layer."""
from itertools import permutations
import json

n = 7
target = {(i, j) for i in range(n) for j in range(i + 1, n)}
for a, b in ((0, 2), (4, 6)):
    target.remove((a, b))
    target.add((b, a))
word = tuple(range(n))
path_pairs = {frozenset((a, b)) for a, b in zip(word, word[1:])}
source = {(a, b) if frozenset((a, b)) in path_pairs else (b, a) for a, b in target}


def rows(arcs):
    return [sum(1 << b for a, b in arcs if a == v) for v in range(n)]


def blue(arcs, path):
    reverse = dict(zip(path, reversed(path)))
    return all((reverse[b], reverse[a]) in arcs for a, b in arcs)


paths = [p for p in permutations(range(n)) if all((a, b) in target for a, b in zip(p, p[1:]))]
scores = [sum(a == 0 for a, b in source), sum(a == 6 for a, b in source)]
if scores != [2, 4] or (6, 0) not in source or not blue(source, word):
    raise RuntimeError('source not in closed blue endpoint-four layer')
if not all((a, b) in source for a, b in zip(word, word[1:])):
    raise RuntimeError('source word not directed')
if len(paths) != 9 or any(not blue(target, p) for p in paths):
    raise RuntimeError('target not pure blue')
print(json.dumps({'source_rows': rows(source), 'target_rows': rows(target),
                  'source_endpoint_scores': scores, 'source_closed': True,
                  'source_blue': True, 'target_paths': len(paths), 'target_black_paths': 0,
                  'all_target_paths': paths}, indent=2))
