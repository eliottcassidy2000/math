"""Reproduce the bounded closed-tree shape scout; no all-order inference."""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import hashlib
import json
import random

HERE = Path(__file__).resolve().parent


def tournament(n, code):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        if code >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return rows


def components(rows, A):
    n = len(rows)
    full = (1 << n) - 1
    B = full ^ A
    paths = {}
    for mask in range(1 << n):
        vertices = [v for v in range(n) if mask >> v & 1]
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
    adjacency = [[] for _ in range(len(old) + len(cuts))]
    stubs = Counter()
    total = unused = 0
    for L in range(1 << n):
        R = full ^ L
        for p in paths[L]:
            if p and not (B >> p[-1] & 1):
                continue
            for q in paths[R]:
                if q and not (A >> q[0] & 1):
                    continue
                total += 1
                pi, ci = old.get(p + q), cuts.get(L)
                if pi is not None and ci is not None:
                    adjacency[pi].append(ci)
                    adjacency[ci].append(pi)
                elif pi is not None or ci is not None:
                    stubs[pi if pi is not None else ci] += 1
                else:
                    unused += 1
    seen, out = set(), []
    for vertex in range(len(adjacency)):
        if vertex in seen:
            continue
        pending, members = [vertex], set()
        while pending:
            x = pending.pop()
            if x in members:
                continue
            members.add(x)
            pending.extend(adjacency[x])
        seen.update(members)
        component = dict(nodes=len(members),
                         edges=sum(len(adjacency[x]) for x in members) // 2,
                         stubs=sum(stubs[x] for x in members),
                         old_paths=[p for p, i in old.items() if i in members],
                         cuts=[L for L, i in cuts.items() if i in members])
        component['closed_tree'] = (not component['stubs'] and
                                    component['edges'] == component['nodes'] - 1)
        out.append(component)
    return dict(old_H=len(old), new_H=total, domination=len(dom),
                unused_paths=unused, components=out)


def main():
    witness_rows, witness_A = [6, 12, 8, 1], 10
    witness = components(witness_rows, witness_A)
    trees = [c for c in witness['components'] if c['closed_tree']]
    assert witness['old_H'] == 5
    assert len(trees) == 1
    assert trees[0]['old_paths'] == [(0, 1, 2, 3)]
    assert trees[0]['cuts'] == [1, 7]
    assert 5 not in trees[0]['cuts']
    rng, shapes, hostile = random.Random(125), Counter(), None
    checked = 0
    for n in range(3, 8):
        for trial in range(500):
            rows = tournament(n, rng.getrandbits(n * (n - 1) // 2))
            A = rng.randrange(1 << n)
            result = components(rows, A)
            checked += 1
            for component in result['components']:
                if not component['closed_tree']:
                    continue
                shape = (len(component['old_paths']), len(component['cuts']))
                shapes[shape] += 1
                if min(shape) > 1 and hostile is None:
                    hostile = dict(n=n, trial=trial, rows=rows, A=A, component=component)
    receipt = dict(status='FINITE-EXACT scout; star theorem OPEN',
                   universe=dict(seed=125, orders=list(range(3, 8)), cases_per_order=500,
                                 checked=checked, filters='none'),
                   nonstar_closed_tree=hostile,
                   closed_tree_shapes=[dict(old_paths=p, cuts=c, occurrences=count)
                                       for (p, c), count in sorted(shapes.items())],
                   hostile_to_canonical_cut=dict(rows=witness_rows, A=witness_A, **witness),
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE / 'exchange_shape_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(checked=checked, nonstar_closed_tree=hostile,
                          closed_tree_shapes=len(shapes))))


if __name__ == '__main__':
    main()
