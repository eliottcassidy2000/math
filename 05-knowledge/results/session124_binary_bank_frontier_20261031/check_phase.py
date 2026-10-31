"""Independent phase-bank audit by LITERAL permutations, not subset DP.

The certificate may be a positional argument. Every failure survives python -O.
"""
import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import permutations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

def need(condition, message):
    if not condition:
        raise ValueError(message)

@lru_cache(None)
def count_paths(graph):
    return sum(all(graph[u] & (1 << v) for u, v in zip(word, word[1:]))
               for word in permutations(range(len(graph))))

def verify(data):
    root = tuple(data['parent'])
    need(len(root) == 5, 'unexpected seed order')
    for u, row in enumerate(root):
        need(type(row) is int and 0 <= row < 1 << len(root) and not row & (1 << u), 'invalid adjacency')
        for v in range(u):
            need(((root[v] >> u) & 1) + ((row >> v) & 1) == 1, 'not a tournament')
    code, k = 0, 0
    for u in range(5):
        for v in range(u + 1, 5):
            code |= ((root[u] >> v) & 1) << k
            k += 1
    need(code == 4, 'wrong parent edgecode')
    need(data['offsets'] == [3, -3, 27], 'wrong phase offsets')
    levels, nodes, edges = defaultdict(list), [], []

    def walk(graph, proof, depth):
        h = count_paths(graph)
        levels[str(depth)].append(h)
        nodes.append({'order': len(graph), 'H': h})
        need(type(proof) is list, 'proof is not a list')
        if depth == 3:
            need(proof == [], 'nonempty leaf')
            return
        need(len(proof) == 2, 'missing binary branch')
        for digit, branch in enumerate(proof):
            need(type(branch) is list and len(branch) == 2, 'bad branch format')
            cut, tail = branch
            n = len(graph)
            need(type(cut) is int and 0 <= cut < (1 << n), 'invalid cut')
            # Build all adjacency entries directly; no imported extension helper.
            child = []
            for u in range(n + 1):
                row = 0
                for v in range(n + 1):
                    if u == v:
                        continue
                    arc = ((graph[u] >> v) & 1) if u < n and v < n else ((cut >> v) & 1) if u == n else 1 - ((cut >> u) & 1)
                    row |= arc << v
                child.append(row)
            child = tuple(child)
            target = 2 * h + data['offsets'][depth] + 2 * digit
            need(count_paths(child) == target, 'wrong child target')
            edges.append((h, cut, target))
            walk(child, tail, depth + 1)

    walk(root, data['proof'], 0)
    levels = dict(levels)
    need(len(nodes) == 15 and len(edges) == 14, 'wrong tree size')
    need(levels == data['direct_DP_levels'], 'claimed level mismatch')
    need(levels == {'0': [5], '1': [13, 15], '2': [23, 25, 27, 29], '3': list(range(73, 88, 2))}, 'wrong levels')
    return {'method': 'literal_vertex_permutations', 'nodes': len(nodes), 'edges': len(edges),
            'levels': levels, 'offsets': data['offsets'], 'seed_edgecode': code}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('certificate', type=Path, nargs='?', default=HERE / 'phase_bank.json')
    args = parser.parse_args()
    data = json.loads(args.certificate.read_text())
    result = verify(data)
    controls = []
    for name, mutate in [
        ('wrong_cut', lambda d: d['proof'][0].__setitem__(0, 0)),
        ('wrong_offset', lambda d: d['offsets'].__setitem__(2, 25)),
        ('missing_branch', lambda d: d['proof'].pop()),
        ('invalid_tournament', lambda d: d['parent'].__setitem__(0, 9)),
        ('nonempty_leaf', lambda d: d['proof'][0][1][0][1][0][1].append([0, []]))]:
        changed = json.loads(json.dumps(data))
        mutate(changed)
        try:
            verify(changed)
        except (ValueError, TypeError, IndexError, KeyError):
            controls.append(name)
        else:
            raise ValueError('hostile control accepted: ' + name)
    result['rejected_mutations'] = controls
    print(json.dumps(result, sort_keys=True))

if __name__ == '__main__':
    main()
