#!/usr/bin/env python3
"""Literal permutation verification of one reversal plus one birth."""
import itertools
import json

def need(b, s):
    if not b:
        raise RuntimeError(s)

graphs = [[8, 13, 17, 20, 35, 15],
          [0, 13, 17, 21, 35, 15],
          [64, 13, 81, 85, 99, 79, 2]]
expected = [(25, 2), (11, 0), (51, 3)]
result = []
for rows, counts in zip(graphs, expected):
    n = len(rows)
    for i, row in enumerate(rows):
        need(0 <= row < 1 << n and not row & (1 << i), 'bad row')
        for j in range(i):
            need(((row >> j) & 1) + ((rows[j] >> i) & 1) == 1, 'bad pair')
    h = sum(all((rows[a] >> b) & 1 for a, b in zip(p, p[1:]))
            for p in itertools.permutations(range(n)))
    c = sum(all((rows[a] >> b) & 1 for a, b in zip(p, p[1:] + p[:1]))
            for p in ((0,) + q for q in itertools.permutations(range(1, n))))
    need((h, c) == counts, 'unexpected path/cycle count')
    result.append({'rows': rows, 'H': h, 'C': c})
need([a ^ b for a, b in zip(graphs[0], graphs[1])] == [8, 0, 0, 1, 0, 0], 'wrong reversal')
need([r & 63 for r in graphs[2][:6]] == graphs[1], 'birth altered another old edge')
need(graphs[2][-1] == 2, 'wrong birth mask')
print(json.dumps({'status': 'PASS', 'graphs': result}, indent=2, sort_keys=True))
