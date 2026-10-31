#!/usr/bin/env python3
"""Finite-field reachable and observable closure with explicit Gaussian reduction."""
import json
from minimal_feedback import model


def basis(vectors, p):
    pivots, original = {}, []
    for v in vectors:
        z = [x % p for x in v]
        for j, row in sorted(pivots.items()):
            c = z[j]
            z = [(a - c * b) % p for a, b in zip(z, row)]
        j = next((i for i, x in enumerate(z) if x), None)
        if j is not None:
            inv = pow(z[j], -1, p)
            pivots[j] = [(x * inv) % p for x in z]
            original.append(v)
    return original


def main():
    matrices, u, v = model()
    matrices = [[[int(x) for x in r] for r in m.tolist()] for m in matrices]
    results = []
    for p in (2, 3, 5, 7, 11):
        R, C = [list(map(int, u))], [list(map(int, v))]
        while True:
            rows = R + [[sum(r[i] * m[i][j] for i in range(36)) % p
                         for j in range(36)] for r in R for m in matrices]
            new = basis(rows, p)
            if len(new) == len(R):
                break
            R = new
        while True:
            cols = C + [[sum(m[i][j] * c[j] for j in range(36)) % p
                         for i in range(36)] for c in C for m in matrices]
            new = basis(cols, p)
            if len(new) == len(C):
                break
            C = new
        H = [[sum(a * b for a, b in zip(r, c)) % p for c in C] for r in R]
        rank = len(basis(H, p))
        if rank != (1 if p == 2 else 17):
            raise RuntimeError(('unexpected closed Hankel rank', p, rank))
        results.append({'characteristic': p, 'reachable_dimension': len(R),
                        'observable_dimension': len(C), 'minimal_dimension': rank})
    print(json.dumps({'status': 'passed', 'closed_finite_field_models': results},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
