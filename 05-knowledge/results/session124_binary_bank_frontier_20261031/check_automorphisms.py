"""Literal permutation check of the independent census automorphism counts.

All producer representatives through order 7 are used only as test inputs;
finite census coverage itself is independently generated and orbit-certified.
"""
import argparse
from collections import Counter
from itertools import permutations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

def need(condition, message):
    if not condition:
        raise ValueError(message)

def decode(code, n):
    g = [[False] * n for _ in range(n)]
    bit = 0
    for u in range(n):
        for v in range(u + 1, n):
            g[u][v] = bool(code & (1 << bit))
            g[v][u] = not g[u][v]
            bit += 1
    return g

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', type=Path, default=HERE / 'small_representatives.json')
    parser.add_argument('--census', type=Path, default=HERE / 'census.normal.json')
    args = parser.parse_args()
    inputs = json.loads(args.inputs.read_text())
    census = {row['n']: row for row in json.loads(args.census.read_text())['orders']}
    result = []
    total = 0
    for n in range(1, 8):
        hist = Counter()
        for code in inputs[str(n)]:
            graph = decode(code, n)
            degree = list(map(sum, graph))
            automorphisms = 0
            for perm in permutations(range(n)):
                total += 1
                if any(degree[u] != degree[perm[u]] for u in range(n)):
                    continue
                if all(graph[u][v] == graph[perm[u]][perm[v]] for u in range(n) for v in range(u + 1, n)):
                    automorphisms += 1
            hist[str(automorphisms)] += 1
        need(dict(hist) == census[n]['automorphism_histogram'], f'automorphism histogram mismatch n={n}')
        result.append({'n': n, 'classes': len(inputs[str(n)]), 'automorphism_histogram': dict(hist)})
    print(json.dumps({'method': 'literal_all_permutations', 'labelings_checked': total, 'orders': result}, sort_keys=True))

if __name__ == '__main__':
    main()
