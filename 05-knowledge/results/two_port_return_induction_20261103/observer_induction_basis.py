#!/usr/bin/env python3
"""Exact adaptive/static 17-query change of basis on the universal 25 carrier."""
from fractions import Fraction as Q
from pathlib import Path
from math import lcm
import hashlib
import json

BASE = Path(__file__).resolve().parents[1] / 'adaptive_state_induction_20261103'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def transpose(A):
    return [list(row) for row in zip(*A)]


def product(A, B):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*B)] for row in A]


def reduce(A):
    A = [list(map(Q, row)) for row in A]
    pivots = []
    for j in range(len(A[0])):
        r = len(pivots)
        k = next((k for k in range(r, len(A)) if A[k][j]), None)
        if k is None:
            continue
        A[r], A[k] = A[k], A[r]
        scale = A[r][j]
        A[r] = [x/scale for x in A[r]]
        for k in range(len(A)):
            if k != r and A[k][j]:
                scale = A[k][j]
                A[k] = [x-scale*y for x, y in zip(A[k], A[r])]
        pivots.append(j)
        if len(pivots) == len(A):
            break
    return A, pivots


def inverse(A):
    n = len(A)
    R, pivots = reduce([list(row)+[int(i == j) for j in range(n)] for i, row in enumerate(A)])
    need(pivots == list(range(n)), 'invertible square query minor')
    return [row[n:] for row in R]


def determinant(A):
    A = [list(map(Q, row)) for row in A]
    result = Q(1)
    for j in range(len(A)):
        k = next((k for k in range(j, len(A)) if A[k][j]), None)
        if k is None:
            return Q(0)
        if k != j:
            A[j], A[k] = A[k], A[j]
            result = -result
        pivot = A[j][j]
        result *= pivot
        for k in range(j+1, len(A)):
            scale = A[k][j]/pivot
            A[k] = [x-scale*y for x, y in zip(A[k], A[j])]
    return result


def main():
    paths = {name: BASE/name for name in ('live_port.json', 'two_port_envelope.json')}
    source = {name: json.loads(path.read_text()) for name, path in paths.items()}
    old = source['live_port.json']
    common = source['two_port_envelope.json']
    matrices, observer = common['matrices'], common['observer']
    def column(word):
        v = list(observer)
        for a in reversed(word):
            w = [0]*len(v)
            for i, j, c in matrices[a]:
                w[i] += c*v[j]
            v = w
        return v
    old_words = [[3*(1 << (a//3))+a%3 for a in word] for word in old['rank']['suffixes']]
    static_words = [[common['rank']['observation_letter_map'][a] for a in word]
                    for word in common['rank']['suffixes']]
    adaptive = transpose([column(w) for w in old_words])
    static = transpose([column(w) for w in static_words])
    _, rows = reduce(transpose(adaptive))
    need(len(rows) == 17, 'adaptive query basis on full universal carrier')
    A, B = [adaptive[i] for i in rows], [static[i] for i in rows]
    forward = product(inverse(A), B)
    backward = inverse(forward)
    need(product(adaptive, forward) == static, 'all 25 coordinates: static=adaptive*forward')
    need(product(static, backward) == adaptive, 'all 25 coordinates: adaptive=static*backward')
    denominators = [lcm(*(x.denominator for row in M for x in row)) for M in (forward, backward)]
    print(json.dumps({
        'scope': 'exact rational query bases on all 25 universal carrier coordinates; no source-image inference',
        'source_sha256': {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()},
        'coordinate_relation': 'X_static = X_adaptive * adaptive_to_static',
        'selected_carrier_rows': rows,
        'adaptive_minor_det': str(determinant(A)), 'static_minor_det': str(determinant(B)),
        'change_det': str(determinant(forward)),
        'denominator_lcm_forward_backward': denominators,
        'both_dyadic': all(d > 0 and d & (d-1) == 0 for d in denominators),
        'adaptive_to_static': [[str(x) for x in row] for row in forward],
        'static_to_adaptive': [[str(x) for x in row] for row in backward],
        'full_carrier_identities_checked': 2*25*17,
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
