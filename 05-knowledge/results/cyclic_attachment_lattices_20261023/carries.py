#!/usr/bin/env python3
"""Exact section carries and both preferred-tree zero-slot switches."""
import argparse
import json
from itertools import product
from pathlib import Path
from nonabelian import mul, inv, ID

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def carry(a, b, m):
    return (a+b-(a+b) % m)//m


def run():
    for m in (3, 5, 9, 15):
        for a, b, c in product(range(m), repeat=3):
            need(carry(a, b, m)+carry((a+b) % m, c, m) ==
                 carry(b, c, m)+carry(a, (b+c) % m, m), 'representative cocycle')
    native = []
    for d in (3, 5, 9):
        size = d*d
        actual = {(a, b) for a in range(d) for b in range(size)
                  if all(((v-u) % d == 0) == ((v+b-u-a) % d == 0)
                         for u in range(d) for v in range(size))}
        reconstructed = {(a, (a+d*t) % size) for a, t in product(range(d), repeat=2)}
        need(actual == reconstructed and len(actual) == size, 'native quotient lifts complete')
        for a, b in product(range(d), repeat=2):
            need((a+b-(a+b) % d) % size == d*carry(a, b, d), 'native kernel carry')
        native.append(dict(cell_sizes=[d, size], compatible_frames=len(actual),
                           quotient_return_lift=[0, d]))

    # A deliberately nonhomomorphic section in the incoming nonabelian nine-cell group.
    # The extension itself splits; nonzero chosen-section carries do not refute that.
    def section(a):
        return (a, a, 0, 0)
    def kappa(a, b):
        return mul(mul(section(a), section(b)), inv(section((a+b) % 3)))
    untwisted_failures = 0
    for a, b, c in product(range(3), repeat=3):
        lhs = mul(kappa(a, b), kappa((a+b) % 3, c))
        transported = mul(mul(section(a), kappa(b, c)), inv(section(a)))
        rhs = mul(transported, kappa(a, (b+c) % 3))
        need(lhs == rhs, 'twisted kernel identity')
        need(kappa(a, b)[0] == 0, 'kernel membership')
        untwisted_failures += lhs != mul(kappa(b, c), kappa(a, (b+c) % 3))
    need(untwisted_failures > 0, 'untwisted hostile must fail')
    need(any(kappa(a, b) != ID for a, b in product(range(3), repeat=2)), 'chosen section carry')

    V = [(1, 0), (0, 1)]
    for k in range(2, 61):
        V.append(tuple(a+b for a, b in zip(V[-1], V[-2])))
    kinds = {'B-to-A': 0, 'C-to-B': 0}
    for k in range(4, 61):
        indices = list(range(k-1, 1, -2)) if k % 2 else list(range(k-1, 2, -2))+[0]
        changed = indices.copy()
        old, new, kind = (2, 0, 'B-to-A') if k % 2 else (3, 2, 'C-to-B')
        changed[changed.index(old)] = new
        for support, debt in ((indices, 1), (changed, 2)):
            need(len(set(support)) == len(support) and 1 not in support and
                 all(abs(a-b) > 1 for a, b in product(support, repeat=2) if a != b), 'proper support')
            total = tuple(sum(V[i][j] for i in support) for j in range(2))
            need(total == (V[k][0], V[k][1]-debt), 'full-branch zero-slot charge')
        kinds[kind] += 1
    return dict(status='FINITE-EXACT; proof in CARRIES.md',
                cyclic_moduli=[3, 5, 9, 15], native_rectangles=native,
                nonabelian_section_triples=27, untwisted_identity_failures=untwisted_failures,
                proper_split_indices=[4, 60], switched_node_types=kinds, checks=CHECKS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'carries.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    else:
        need(json.loads(path.read_text()) == result, 'carry receipt mismatch')
    print(json.dumps(result, indent=2))
