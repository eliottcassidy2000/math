"""Exact controls for NEXT_BOUNDARY.md; no assertion is disabled by -O."""
import importlib.util
import json
from pathlib import Path


spec = importlib.util.spec_from_file_location('base_ports_check', Path(__file__).with_name('check.py'))
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
need = base.need


def repair(n, k):
    a = tuple(range(k + 1))
    u, v, w, s = k + 1, k + 2, k + 3, k + 4
    z = tuple(range(k + 5, n))
    return [a + (w, u, v, s) + z,
            a + (w, u) + z + (v, s),
            a + (u,) + z + (v, w, s)]


def certify(src, words):
    n = src.n
    flips = set()
    universe = set(src.edges)
    for step, q in enumerate(words):
        need(sorted(q) == list(range(n)), ('word vertex set', n, step))
        for u, v in zip(q, q[1:]):
            need(src.forced(u, v, 1 ^ (base.pair(u, v) in flips)),
                 ('join', n, step, u, v))
        trial = src.clone()
        trial.anti(q, flips)
        need(not trial.ok, ('possibly blue', n, step))
        flips ^= universe ^ base.pathset(q)
    need(flips == universe ^ base.pathset(tuple(range(n))), 'literal target')
    for seed in (0, 917):
        base.literal(src.model(seed), words)


def cell(n, interior, forward):
    k = len(interior) + 1
    h = min(interior)
    words, guards = base.general(n, h, forward)
    src, report = base.symbolic(n, interior, words, guards, False)
    need(not report['illegal'] and report['target'], ('ledger', n, interior, forward))
    need(not any(report['blue_feasible'][:2]), ('first two moves', n, interior))
    if not forward:
        need(src.forced(2, 4, 0) and src.forced(n-3, n-1, 1), 'local cyclic witness')
    if not report['blue_feasible'][2]:
        certify(src, words)
        if forward and n == 2*k+3 and h == k+2:
            end = set(src.edges) ^ base.pathset(tuple(range(n)))
            need(src.forced(0, 1, 1 ^ (base.pair(0, 1) in end)), 'B(a,b)=1')
            need(src.forced(n-1, 0, 0 ^ (base.pair(n-1, 0) in end)), 'B(z,a)=0')
        return 0
    need(forward and n == 2*k+3 and h == k+3 and
         interior == set(range(k+3, 2*k+2)), ('unexpected exception', n, interior))
    flips = set(src.edges) ^ base.pathset(tuple(range(n)))
    src.anti(words[2], flips)
    need(src.ok, 'exception must exist')
    root, _ = src.root(src.const)
    need(all(src.root(i)[0] == root for i in range(src.const)), 'exception must be unique')
    certify(src, repair(n, k))
    return 1


def main():
    protected = exceptional = 0
    for n in range(9, 18):
        fixed = {n-4, n-2}
        free = sorted(set(range(5, n-1)) - fixed)
        for mask in range(1 << len(free)):
            interior = fixed | {v for i, v in enumerate(free) if (mask >> i) & 1}
            for forward in (False, True):
                exceptional += cell(n, interior, forward)
                protected += 1
    one_hole = 0
    for k in range(5, 25):
        for n in range(k+6, 2*k+8):
            interior = {n-k-1} | set(range(n-k+1, n-1))
            for forward in (False, True):
                need(cell(n, interior, forward) == 0, 'one-hole exception impossible')
                one_hole += 1
    for k in range(4, 41):
        n = 2*k+3
        interior = set(range(k+3, n-1))
        need(cell(n, interior, True) == 1, 'consecutive cyclic exception expected')
    words, guards = base.general(15, 4, True)
    hostile, report = base.symbolic(15, {4, 11, 13}, words, guards, False)
    flips = set()
    universe = set(hostile.edges)
    for q in words[:2]:
        flips ^= universe ^ base.pathset(q)
    need(hostile.forced(3, 0, 0 ^ (base.pair(3, 0) in flips)), 'h=4 false join')
    need([2, 3, 0] in report['illegal'], 'h=4 hostile retained')
    print(json.dumps({
        'status': 'PASS',
        'protected_affine_guard_cells': protected,
        'protected_order_range': [9, 17],
        'protected_conditional_exceptions': exceptional,
        'one_hole_affine_guard_cells': one_hole,
        'one_hole_k_range': [5, 24],
        'one_hole_n_range': 'k+6 .. 2k+7',
        'cyclic_exception_k_range': [4, 40],
        'literal_models_per_certified_cell': 2,
        'hostile': 'n=15,I={4,11,13},2->4 forces the final join 3->0 false',
        'scope': 'all protected port sets; unrestricted black-core connectivity remains open'
    }, indent=2))


if __name__ == '__main__':
    main()
