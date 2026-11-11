"""Hostile probe for a stronger, fibre-respecting domination injection.

The counts are exact, but no all-order theorem is inferred from this scout.
"""
from itertools import combinations, permutations
import argparse
import json
from pathlib import Path


def probe(rows, A):
    n = len(rows)
    full = (1 << n) - 1
    B = full ^ A
    paths = [[] for _ in range(1 << n)]
    paths[0] = [()]
    for size in range(1, n + 1):
        for vertices in combinations(range(n), size):
            mask = sum(1 << x for x in vertices)
            paths[mask] = [p for p in permutations(vertices)
                           if all(rows[x] >> y & 1 for x, y in zip(p, p[1:]))]
    old = {p: i for i, p in enumerate(paths[full])}
    dom = []
    for L in range(1 << n):
        R = full ^ L
        left_ok = all(rows[a] & L & B for a in range(n) if L & A & (1 << a))
        right_ok = all(any(rows[a] >> b & 1 for a in range(n) if R & A & (1 << a))
                       for b in range(n) if R & B & (1 << b))
        if left_ok and right_ok:
            dom.append(L)
    cuts = {L: len(old) + i for i, L in enumerate(dom)}
    parent = list(range(len(old) + len(cuts)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    edges, stubs, unused = [], [], 0
    total = 0
    for L in range(1 << n):
        R = full ^ L
        for p in paths[L]:
            if p and not (B >> p[-1] & 1):
                continue
            for q in paths[R]:
                if q and not (A >> q[0] & 1):
                    continue
                total += 1
                valid = not p or not q or (rows[p[-1]] >> q[0] & 1)
                pi = old[p + q] if valid else None
                ci = cuts.get(L)
                if pi is not None and ci is not None:
                    x, y = find(pi), find(ci)
                    parent[x] = y
                    edges.append((pi, ci))
                elif pi is not None or ci is not None:
                    stubs.append(pi if pi is not None else ci)
                else:
                    unused += 1
    components = {}
    for i in range(len(parent)):
        c = components.setdefault(find(i), {'nodes': 0, 'edges': 0, 'stubs': 0})
        c['nodes'] += 1
    for x, y in edges:
        components[find(x)]['edges'] += 1
    for x in stubs:
        components[find(x)]['stubs'] += 1
    trees = [v for v in components.values() if not v['stubs'] and v['edges'] == v['nodes'] - 1]
    return {'old_H': len(old), 'new_H': total, 'domination': len(dom),
            'slack': total - len(old) - len(dom) + 1,
            'closed_tree_components': len(trees), 'unused_paths': unused,
            'components': list(components.values())}


def tournament(n, code):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        if code >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    checked = 0
    max_trees = 0
    stages = []
    for n in range(1, 6):
        for code in range(1 << (n * (n - 1) // 2)):
            rows = tournament(n, code)
            for A in range(1 << n):
                result = probe(rows, A)
                checked += 1
                max_trees = max(max_trees, result['closed_tree_components'])
                if result['closed_tree_components'] > 1:
                    print(json.dumps({'status': 'HOSTILE to fibre-respecting matching',
                                      'n': n, 'rows': rows, 'A': A, 'checked': checked,
                                      **result}, indent=2))
                    return
        stage = {'finished_order': n, 'checked': checked,
                 'maximum_closed_tree_components': max_trees}
        stages.append(stage)
        print(json.dumps(stage), flush=True)
    hostile = probe([6,12,8,1], 10)
    receipt = {'status': 'FINITE-EXACT; all-order component bound remains OPEN',
               'universe': 'Every labelled old tournament of orders 1..5 and every newborn cut',
               'checked': checked, 'stages': stages,
               'nontransitive_tree_without_cut_B': {'rows': [6,12,8,1], 'A': 10, **hostile}}
    if args.write:
        Path(__file__).with_name('receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    main()
