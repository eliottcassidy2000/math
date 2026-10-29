"""Exact probes of native three-run observations after legal complements."""
from collections import Counter
from itertools import permutations
import json
from pathlib import Path
from time import perf_counter
from functools import reduce
from math import lcm, gcd
from random import Random

import numpy as np
from sympy import Matrix
from sympy.polys.matrices import DomainMatrix

from three_run import native, need

HERE = Path(__file__).resolve().parent


def rank_mod(M, p=1009):
    A = np.array([[int(x) % p for x in row] for row in M], dtype=np.int64)
    rank = 0
    for col in range(A.shape[1]):
        choices = np.flatnonzero(A[rank:, col])
        if not len(choices):
            continue
        pivot = rank + int(choices[0])
        A[[rank, pivot]] = A[[pivot, rank]]
        A[rank] = A[rank]*pow(int(A[rank, col]), -1, p) % p
        A[rank+1:] = (A[rank+1:] - A[rank+1:, col, None]*A[rank]) % p
        rank += 1
        if rank == min(A.shape):
            break
    return rank


def rational_kernel(M):
    K = DomainMatrix.from_Matrix(Matrix(M)).to_field().nullspace().to_Matrix()
    vectors = []
    for row in K.tolist():
        denominator = lcm(*(int(x.q) for x in row))
        vector = [int(x*denominator) for x in row]
        content = reduce(gcd, vector)
        vectors.append([x//content for x in vector])
    need(all(sum(x*y for x, y in zip(row, v)) == 0 for row in M for v in vectors),
         'integer null vector certificate')
    return vectors


def exact_rank_certificate(M):
    size = len(M[0])
    lower = rank_mod(M)
    if lower == size:
        return dict(rank=size, prime=1009, kernel=[])
    K = rational_kernel(M)
    upper = size-len(K)
    if lower == upper and rank_mod(K) == len(K):
        return dict(rank=upper, prime=1009, kernel=K)
    for p in (1013, 1019, 10007):
        if rank_mod(M, p) == upper and (not K or rank_mod(K, p) == len(K)):
            return dict(rank=upper, prime=p, kernel=K)
    raise ValueError('need another good prime')


def independent_rows(M, p, order=None):
    order = list(range(len(M))) if order is None else order
    A = np.array([[int(M[i][j]) % p for i in order] for j in range(len(M[0]))], dtype=np.int64)
    rank, chosen = 0, []
    for col in range(A.shape[1]):
        choices = np.flatnonzero(A[rank:, col])
        if not len(choices):
            continue
        pivot = rank+int(choices[0])
        A[[rank, pivot]] = A[[pivot, rank]]
        A[rank] = A[rank]*pow(int(A[rank, col]), -1, p) % p
        A[rank+1:] = (A[rank+1:] - A[rank+1:, col, None]*A[rank]) % p
        chosen.append(order[col])
        rank += 1
        if rank == A.shape[0]:
            return chosen
    return []


def unit_certificate(M):
    # Coprime full-size minors yield an integral left inverse by Bezout + adjugates.
    records, divisor = [], 0
    rng = Random(20261029)
    for p in (2, 3, 5, 7, 11, 1009, 1009, 1009):
        order = list(range(len(M)))
        if len(records) >= 3:
            rng.shuffle(order)
        chosen = independent_rows(M, p, order)
        if not chosen:
            return None
        det = int(DomainMatrix.from_Matrix(Matrix([M[i] for i in chosen])).det())
        need(det % p != 0, 'selected minor is a unit at its certifying prime')
        records.append(dict(prime=p, rows=chosen, determinant=det))
        divisor = gcd(divisor, det)
        if divisor == 1:
            return records
    return None


def unit_census():
    atlas = json.loads((HERE.parent / 'strong_atom_continuation_vitali_20261026' / 'strong-atlas.json').read_text())
    records = []
    for obj in atlas:
        n, rows = obj['n'], obj['source']
        if not 3 <= n <= 7:
            continue
        M = native(rows)[2]
        cert = unit_certificate(M)
        chosen_path = None
        if cert is None:
            for path in paths(rows):
                N = native(complement(rows, path))[2]
                cert = unit_certificate(M+N)
                if cert is not None:
                    chosen_path = path
                    break
        need(cert is not None, 'one selected complement supplies an integral left inverse')
        records.append(dict(n=n, source=rows, path=chosen_path, minors=cert))
        if len(records) % 10 == 0 or n <= 5:
            print('unit',len(records),'seconds',round(perf_counter()-START,1),flush=True)
        (HERE / 'complement_units.json').write_text(json.dumps(records,separators=(',',':'))+'\n')
    print('unit census',Counter(x['n'] for x in records),'depth0',Counter(x['n'] for x in records if x['path'] is None))


def complement(rows, path):
    need(sorted(path) == list(range(len(rows))), 'a Hamilton marking is a permutation of every vertex')
    full = (1 << len(rows)) - 1
    target = [full ^ row ^ (1 << i) for i, row in enumerate(rows)]
    for u, v in zip(path, path[1:]):
        need(rows[u] >> v & 1, 'legal Hamilton path')
        target[u] ^= 1 << v
        target[v] ^= 1 << u
    return target


def paths(rows):
    return [p for p in permutations(range(len(rows)))
            if all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))]


def census():
    atlas = json.loads((HERE.parent / 'strong_atom_continuation_vitali_20261026' / 'strong-atlas.json').read_text())
    records = []
    for obj in atlas:
        n = obj['n']
        if not 3 <= n <= 7:
            continue
        rows = obj['source']
        _, ports, M = native(rows)
        cert = exact_rank_certificate(M)
        item = dict(n=n, source=rows, prime_atom=obj['prime'], aut=obj['Aut'], source_certificate=cert,
                    size=len(ports), path_count=obj['H'])
        if cert['rank'] < len(ports):
            checks = []
            for path in paths(rows):
                target = complement(rows, path)
                need(complement(target, path) == rows, 'actual complement involution')
                _, _, N = native(target)
                if n <= 6:
                    target_cert = exact_rank_certificate(M+N)
                else:
                    rank = rank_mod(M+N)
                    if rank < len(ports):
                        continue
                    target_cert = dict(rank=rank, prime=1009, kernel=[])
                checks.append(dict(path=path, target=target, certificate=target_cert))
                if n == 7:
                    break
            item['complements'] = checks
            need(any(c['certificate']['rank'] == len(ports) for c in checks),
                 'at least one successful complement for each deficient strong source')
        records.append(item)
        if n <= 6 or len(records) % 20 == 0:
            print(n, len(records), 'rank',cert['rank'],'seconds',round(perf_counter()-START,1), flush=True)
    (HERE / 'complement_kernel.json').write_text(json.dumps(records, separators=(',', ':'))+'\n')
    print('census',Counter(x['n'] for x in records),'deficient',Counter(x['n'] for x in records if x['source_certificate']['rank']<x['size']))


if __name__ == '__main__':
    START = perf_counter()
    import sys
    if '--build-units' in sys.argv:
        unit_census()
    elif '--build-ranks' in sys.argv:
        census()
    else:
        from complement_verify import verify
        verify()
    print('seconds', perf_counter()-START)
