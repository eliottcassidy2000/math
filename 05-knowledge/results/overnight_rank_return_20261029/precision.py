"""Unimodular Smith certificate for the ten-by-fifteen guard incidence map."""
import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise ValueError(why)


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def product(A, B):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*B)] for row in A]


def smith(B):
    A = [r[:] for r in B]
    nr, nc = len(A), len(A[0])
    U, V, ops = identity(nr), identity(nc), []
    def rowadd(i, j, c):
        A[i] = [a+c*b for a, b in zip(A[i], A[j])]
        U[i] = [a+c*b for a, b in zip(U[i], U[j])]
        ops.append(['rowadd', i, j, c])
    def coladd(i, j, c):
        for M in (A, V):
            for row in M:
                row[i] += c*row[j]
        ops.append(['coladd', i, j, c])
    def rowswap(i, j):
        for M in (A, U):
            M[i], M[j] = M[j], M[i]
        ops.append(['rowswap', i, j])
    def colswap(i, j):
        for M in (A, V):
            for row in M:
                row[i], row[j] = row[j], row[i]
        ops.append(['colswap', i, j])
    for k in range(min(nr, nc)):
        entries = [(abs(A[i][j]), i, j) for i in range(k, nr)
                   for j in range(k, nc) if A[i][j]]
        if not entries:
            break
        _, i, j = min(entries)
        rowswap(k, i)
        colswap(k, j)
        while True:
            moved = False
            for i in range(k+1, nr):
                if A[i][k]:
                    rowadd(i, k, -(A[i][k] // A[k][k]))
                    if A[i][k]:
                        rowswap(i, k)
                    moved = True
                    break
            if moved:
                continue
            for j in range(k+1, nc):
                if A[k][j]:
                    coladd(j, k, -(A[k][j] // A[k][k]))
                    if A[k][j]:
                        colswap(j, k)
                    moved = True
                    break
            if moved:
                continue
            bad = next(((i, j) for i in range(k+1, nr) for j in range(k+1, nc)
                        if A[i][j] % A[k][k]), None)
            if bad is None:
                break
            rowadd(k, bad[0], 1)
        if A[k][k] < 0:
            A[k] = [-x for x in A[k]]
            U[k] = [-x for x in U[k]]
            ops.append(['rowneg', k])
    need(product(product(U, B), V) == A, 'dense transform verification')
    return A, U, V, ops


def replay(B, ops):
    M = [r[:] for r in B]
    for op in ops:
        kind, i, *rest = op
        if kind == 'rowneg':
            M[i] = [-x for x in M[i]]
        elif kind == 'rowswap':
            j, = rest
            M[i], M[j] = M[j], M[i]
        elif kind == 'colswap':
            j, = rest
            M = [[r[j] if k == i else r[i] if k == j else r[k]
                  for k in range(len(r))] for r in M]
        elif kind == 'rowadd':
            j, c = rest
            need(i != j, 'unimodular row addition')
            M[i] = [x+c*y for x, y in zip(M[i], M[j])]
        else:
            need(kind == 'coladd', 'known operation')
            j, c = rest
            need(i != j, 'unimodular column addition')
            M = [[r[k]+c*r[j] if k == i else r[k] for k in range(len(r))] for r in M]
    return M


def direct_image_size(B, modulus):
    image = {(0,)*len(B)}
    for column in zip(*B):
        image = {tuple((a+c*b) % modulus for a, b in zip(v, column))
                 for v in image for c in range(modulus)}
    return len(image)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    data = json.loads((HERE/'fifteen.json').read_text(encoding='utf-8'))
    B = [[row >> i & 1 for i in range(15)] for row in data['incidence_rows']]
    D, U, V, ops = smith(B)
    expected = [1]*5+[2]*4+[4]
    need(D == [[expected[i] if i == j else 0 for j in range(15)] for i in range(10)],
         'Smith diagonal')
    need(replay(B, ops) == D, 'elementary-operation replay')
    direct = {str(k): direct_image_size(B, 2**k) for k in (1, 2)}
    exponents = [5*k+4*max(k-1, 0)+max(k-2, 0) for k in range(11)]
    need(direct == {'1': 32, '2': 16384}, 'independent finite image enumeration')
    result = dict(smith_diagonal=expected, row_transform=U, column_transform=V,
                  unimodular_operations=ops, direct_image_sizes=direct,
                  log2_image_sizes_through_precision_10=exponents,
                  lattice_index=64, cokernel='(Z/2)^4 + Z/4')
    path = Path(__file__).with_suffix('.json')
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    else:
        need(json.loads(path.read_text(encoding='utf-8')) == result, 'stored precision receipt')
    print(json.dumps({k: result[k] for k in ('smith_diagonal', 'direct_image_sizes',
                     'log2_image_sizes_through_precision_10', 'lattice_index')}))


if __name__ == '__main__':
    main()
