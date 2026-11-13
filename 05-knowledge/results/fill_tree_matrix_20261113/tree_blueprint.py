"""Exact controls for the contracted arithmetic target tree, not a tournament compiler."""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
HOLES = {7, 21}

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def parent(h):
    need(h > 1 and h % 2 == 1, 'parent domain')
    return (h + 1) // 2 if h % 4 == 1 else (h - 1) // 2

def contracted_parent(h):
    p = parent(h)
    return parent(p) if p in HOLES else p

def odd_collatz(n):
    n = 3 * n + 1
    while n % 2 == 0:
        n //= 2
    return n

def main():
    limit = 131071
    counts = {}
    depth = {1: 0}
    splices = {}
    for h in range(3, limit + 1, 2):
        if h in HOLES:
            continue
        p = contracted_parent(h)
        need(p in depth and p < h, ('rooted descent', h, p))
        depth[h] = depth[p] + 1
        counts[p] = counts.get(p, 0) + 1
        if p != parent(h):
            splices[h] = p
        m = (h - 1) // 2
        need((parent(h) - 1) // 2 == m // 2, 'binary parent')
    need(splices == {13: 3, 15: 3, 41: 11, 43: 11}, 'exact splices')
    expected = {5: [5, 1], 17: [17, 13, 5, 1],
                7: [7, 11, 17, 13, 5, 1], 21: [21, 1]}
    for start, route in expected.items():
        need(all(odd_collatz(a) == b for a, b in zip(route, route[1:])), ('Collatz route', start))
    need(2 * (2 * 5 - 1) - 1 == 17 and 3 * 5 + 2 == 17 and 3 * 7 == 21, 'four-number relations')
    report = dict(status='PASS', scope='Arithmetic tree only; physical tournament lift OPEN',
                  universe=dict(odd_limit=limit, holes=sorted(HOLES)),
                  vertices=len(depth), max_depth=max(depth.values()),
                  splice_parent=splices, collatz_odd_routes=expected,
                  proof_sha256=hashlib.sha256((HERE / 'TREE_BLUEPRINT.md').read_bytes()).hexdigest())
    suffix = '.optimized' if sys.flags.optimize else ''
    (HERE / ('tree_blueprint' + suffix + '.json')).write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report))

if __name__ == '__main__':
    main()
