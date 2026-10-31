#!/usr/bin/env python3
"""Exact reachable/observable reduction; proofs use rational identities, not fits."""
from itertools import product, permutations
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent


def model():
    states = list(product(range(3), range(3), range(2)))
    indexes = {s: i for i, s in enumerate(states)}
    matrices = [sp.zeros(36) for _ in range(4)]
    start, finish = sp.zeros(1, 36), sp.zeros(36, 1)
    for block, p in enumerate(permutations(range(2))):
        offset = 18 * block
        start[offset] = 1
        for state in states:
            i = offset + indexes[state]
            if 2 not in state and (state[1] or p[0] == 0):
                finish[i] = 1
            for a in range(4):
                for gap in range(3):
                    if state[gap] == 0 and gap and not (a >> p[gap - 1] & 1):
                        continue
                    t = list(state)
                    t[gap] = int(gap != 2 and bool(a >> p[gap] & 1)) + 1
                    matrices[a][i, offset + indexes[tuple(t)]] += 1
    return matrices, start, finish


def build():
    matrices, u, v = model()
    rows, prefixes = [u], [()]
    i = 0
    while i < len(rows):
        for a, m in enumerate(matrices):
            r = rows[i] * m
            if sp.Matrix.vstack(*rows, r).rank() > len(rows):
                rows.append(r)
                prefixes.append(prefixes[i] + (a,))
        i += 1
    cols, suffixes = [v], [()]
    i = 0
    while i < len(cols):
        for a, m in enumerate(matrices):
            c = m * cols[i]
            if sp.Matrix.hstack(*cols, c).rank() > len(cols):
                cols.append(c)
                suffixes.append((a,) + suffixes[i])
        i += 1
    R, C = sp.Matrix.vstack(*rows), sp.Matrix.hstack(*cols)
    pivots = (R * C).rref()[1]
    O = C[:, pivots]
    row_pivots = (R * O).T.rref()[1]
    B = R[list(row_pivots), :]
    K = B * O
    if len(pivots) != 17 or K.det() == 0:
        raise RuntimeError('rank claim failed')
    reduced = [K.inv() * B * m * O for m in matrices]
    for m, a in zip(matrices, reduced):
        if R * m * O != R * O * a:
            raise RuntimeError('universal transfer identity failed')
    if sp.Matrix.vstack(R, *(R * m for m in matrices)).rank() != len(rows):
        raise RuntimeError('reachable closure failed')
    # O includes v as its first column, so H is literally coordinate zero.
    if O[:, 0] != v:
        raise RuntimeError('lost ordinary H coordinate')
    encode = lambda m: [[str(m[i, j]) for j in range(m.cols)] for i in range(m.rows)]
    result = {'status': 'passed', 'reachable_dimension': len(rows),
              'observable_dimension': len(cols), 'minimal_rational_dimension': 17,
              'hankel_minor_determinant': str(K.det()),
              'prefix_basis': [prefixes[i] for i in row_pivots],
              'suffix_basis': [suffixes[i] for i in pivots],
              'initial_observations': encode(u * O),
              'hankel_minor': encode(K), 'matrices': [encode(a) for a in reduced],
              'scope': 'Exact scalar series H for a transitive chain with two marked feedback vertices x->y.'}
    (HERE / 'minimal_feedback.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('matrices', 'hankel_minor')}, sort_keys=True, indent=2))


if __name__ == '__main__':
    build()
