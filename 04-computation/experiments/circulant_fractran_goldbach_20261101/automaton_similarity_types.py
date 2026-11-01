#!/usr/bin/env python3
"""mac-mini-2026-11-01-lenses: does the GF(2) automaton x -> Ax classify the merged metagraph's nodes?

Exploratory observation (tournament_automaton.out, section A5): the coarse dynamical type (rank, nilpotent depth,
cycle multiset) of x -> Ax takes exactly 10 values on the 12 classes of 5-tournaments, one per merged node
(converse pair). Here the FINEST linear invariant is used: GF(2)-similarity of A, i.e. the GF(2)[x]-module type,
determined by the ranks of p(A)^k for all irreducible p over GF(2) and all k (elementary divisor theory). Since
similarity is preserved by isomorphism and A ~ A^T, the number of similarity types is at most the number of merged
nodes; the automaton classifies the merged nodes exactly when the two numbers are equal.

PRE-REGISTERED (before the first run of this file):
  S1  n = 3, 4, 5: the number of similarity types equals the number of merged nodes (2, 3, 10).
  S2  n = 6: the number of similarity types is smaller than the number of merged nodes (34), i.e. the automaton
      stops classifying merged nodes at n = 6 (the coarse invariant had 23 types there).
  S3  Report the colliding merged nodes at n = 6 (no prediction).
usage: python3 automaton_similarity_types.py
"""
import itertools
import sys

FAIL = []


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAIL.append(name)


def irreducibles(maxdeg):
    out = []
    for d in range(1, maxdeg + 1):
        for p in range(1 << d, 1 << (d + 1)):
            if all(pmod(p, q) != 0 for q in out if (q.bit_length() - 1) <= d // 2):
                out.append(p)
    return out


def pmod(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def matmul(A, B, n):
    C = []
    for r in A:
        acc = 0
        for j in range(n):
            if r >> j & 1:
                acc ^= B[j]
        C.append(acc)
    return C


def rank(rows, n):
    rows = list(rows)
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, n) if rows[i] >> c & 1), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(n):
            if i != r and rows[i] >> c & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def polyeval(p, A, n):
    I = [1 << i for i in range(n)]
    R = [0] * n
    for k in range(p.bit_length() - 1, -1, -1):
        R = matmul(R, A, n)
        if p >> k & 1:
            R = [R[i] ^ I[i] for i in range(n)]
    return R


def similarity_type(A, n, IRR):
    inv = []
    for p in IRR:
        d = p.bit_length() - 1
        if d > n:
            continue
        P = polyeval(p, A, n)
        ranks = []
        M = P
        for _ in range(n // d + 1):
            rk = rank(M, n)
            ranks.append(rk)
            if rk == n or (len(ranks) > 1 and ranks[-1] == ranks[-2]):
                break
            M = matmul(M, P, n)
        if ranks[0] != n:
            inv.append((p, tuple(ranks)))
    return tuple(inv)


def tilings(n):
    tiles = [(i, j) for i in range(n) for j in range(i + 2, n)]
    for bits in range(1 << len(tiles)):
        A = [[0] * n for _ in range(n)]
        for i in range(n - 1):
            A[i + 1][i] = 1
        for k, (i, j) in enumerate(tiles):
            if bits >> k & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def canon(A, n, perms):
    return min(tuple(A[p[i]][p[j]] for i in range(n) for j in range(n)) for p in perms)


def main():
    IRR = irreducibles(7)
    merged_expected = {3: 2, 4: 3, 5: 10, 6: 34}
    for n in (3, 4, 5, 6):
        perms = list(itertools.permutations(range(n)))
        cls = {}
        for A in tilings(n):
            key = canon(A, n, perms)
            if key not in cls:
                cls[key] = A
        # merged nodes: identify a class with its converse
        merged = {}
        for key, A in cls.items():
            At = [[A[j][i] for j in range(n)] for i in range(n)]
            k2 = canon(At, n, perms)
            merged[key] = min(key, k2)
        nodes = set(merged.values())
        types = {}
        for key, A in cls.items():
            rows = [sum(1 << v for v in range(n) if A[u][v]) for u in range(n)]
            t = similarity_type(rows, n, IRR)
            types.setdefault(t, set()).add(merged[key])
        ntypes = len(types)
        print(f'  n={n}: {len(cls)} classes, {len(nodes)} merged nodes, {ntypes} GF(2)-similarity types')
        if n <= 5:
            check(f'S1: n={n} similarity types = merged nodes', ntypes == len(nodes) == merged_expected[n])
        else:
            check('S2: n=6 similarity types < 34 merged nodes', ntypes < len(nodes) == 34, f'{ntypes} types')
            coll = [sorted(v) for v in types.values() if len(v) > 1]
            print(f'  S3: {len(coll)} similarity types hold more than one merged node; sizes '
                  f'{sorted(len(c) for c in coll)}')
            # describe one collision by score sequences
            for c in coll[:3]:
                desc = []
                for key in c:
                    A = [[key[i * n + j] for j in range(n)] for i in range(n)]
                    desc.append(sorted(sum(r) for r in A))
                print(f'      collision: score sequences {desc}')
    print()
    print('ALL PRE-REGISTERED CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
    return 0 if not FAIL else 1


if __name__ == '__main__':
    sys.exit(main())
