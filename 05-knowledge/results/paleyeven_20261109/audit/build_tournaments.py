#!/usr/bin/env python3
"""Audit P1 (independent): build the Part A tournaments from the definitions, without the session's builders.

Output, one line per tournament:  name n out_0 out_1 ... out_{n-1}
where out_i is the decimal bitmask of out-neighbours of vertex i.

Objects:
  QR_p      x -> y iff y - x is a nonzero square mod p (p = 3 mod 4)
  C9_1235   x -> x + s mod 9, s in {1, 2, 3, 5}
  D(B)      a_x = x, b_x = q + x; A and B' copies of B; a_x -> b_y (x != y); b_x -> a_x
Also checks: tournament axioms, double regularity of QR_p (every pair has (q-3)/4 common out-neighbours),
S S^T = qI - J, and the automorphism / anti-automorphism claims for the doubles.
"""
import sys


def qr_tour(p):
    sq = {(x * x) % p for x in range(1, p)}
    return [[((y - x) % p) in sq and x != y for y in range(p)] for x in range(p)]


def circ_tour(m, S):
    S = {s % m for s in S}
    return [[((y - x) % m) in S for y in range(m)] for x in range(m)]


def double(B):
    q = len(B)
    n = 2 * q
    A = [[False] * n for _ in range(n)]
    for x in range(q):
        for y in range(q):
            if B[x][y]:
                A[x][y] = True            # a_x -> a_y
                A[q + x][q + y] = True    # b_x -> b_y
            if x != y:
                A[x][q + y] = True        # a_x -> b_y
        A[q + x][x] = True                # b_x -> a_x
    return A


def is_tournament(A):
    n = len(A)
    for i in range(n):
        if A[i][i]:
            return False
        for j in range(i + 1, n):
            if A[i][j] == A[j][i]:
                return False
    return True


def masks(A):
    return [sum(1 << j for j in range(len(A)) if A[i][j]) for i in range(len(A))]


def doubly_regular(B):
    q = len(B)
    lam = (q - 3) // 4
    for u in range(q):
        for v in range(u + 1, q):
            c = sum(1 for w in range(q) if B[u][w] and B[v][w])
            if c != lam:
                return False
    return True


def sst_check(B):
    q = len(B)
    S = [[(1 if B[w][y] else (-1 if B[y][w] else 0)) for y in range(q)] for w in range(q)]
    for i in range(q):
        for j in range(q):
            e = sum(S[i][k] * S[j][k] for k in range(q))
            want = (q - 1) if i == j else -1
            if e != want:
                return False
    return True


def is_iso_map(A, Bm, f):
    n = len(A)
    return all(A[i][j] == Bm[f[i]][f[j]] for i in range(n) for j in range(n))


def main():
    objs = []
    for p in (7, 11):
        objs.append((f'QR{p}', qr_tour(p)))
    for p in (7, 11):
        objs.append((f'T{p}', double(qr_tour(p))))
    C = circ_tour(9, [1, 2, 3, 5])
    objs.append(('C9_1235', C))
    objs.append(('D_C9_1235', double(C)))
    objs.append(('QR3', qr_tour(3)))
    objs.append(('T3', double(qr_tour(3))))
    for name, A in objs:
        assert is_tournament(A), name
        print(name, len(A), *masks(A))
    # side checks to stderr
    for p in (7, 11, 19, 23):
        B = qr_tour(p)
        print(f'# QR{p}: tournament={is_tournament(B)} doubly_regular={doubly_regular(B)} SST=qI-J:{sst_check(B)}',
              file=sys.stderr)
    print(f'# C9_1235 doubly regular: {doubly_regular(C)}', file=sys.stderr)
    for p in (7, 11):
        B = qr_tour(p)
        D = double(B)
        q = p
        # translation x -> x + 1 acting diagonally
        f = [((x + 1) % q) for x in range(q)] + [q + ((x + 1) % q) for x in range(q)]
        # multiplier x -> 2x is an automorphism of QR_p iff 2 is a square mod p (p = 7: yes; p = 11: no -> use 3)
        mult = 2 if pow(2, (p - 1) // 2, p) == 1 else 3
        g = [((mult * x) % q) for x in range(q)] + [q + ((mult * x) % q) for x in range(q)]
        # anti-automorphism a_x -> b_{-x}, b_x -> a_{-x}: an isomorphism D -> D^op
        Dop = [[D[j][i] for j in range(2 * q)] for i in range(2 * q)]
        h = [q + ((-x) % q) for x in range(q)] + [((-x) % q) for x in range(q)]
        print(f'# T{p}: translation aut={is_iso_map(D, D, f)}, multiplier {mult} aut={is_iso_map(D, D, g)}, '
              f'a_x->b_-x anti-aut={is_iso_map(D, Dop, h)}', file=sys.stderr)
    # C9 double: translation and anti-automorphism x -> -x of C9_1235?
    D = double(C)
    q = 9
    f = [((x + 1) % q) for x in range(q)] + [q + ((x + 1) % q) for x in range(q)]
    Dop = [[D[j][i] for j in range(18)] for i in range(18)]
    Cop = [[C[j][i] for j in range(9)] for i in range(9)]
    anti = [x for x in range(9) if False]
    # search all affine maps x -> u x + t of Z_9 (u a unit) for an anti-automorphism of C
    found = None
    for u in (1, 2, 4, 5, 7, 8):
        for t in range(9):
            phi = [((u * x + t) % 9) for x in range(9)]
            if is_iso_map(C, Cop, phi):
                found = (u, t)
                break
        if found:
            break
    print(f'# D(C9_1235): translation aut={is_iso_map(D, D, f)}; affine anti-automorphism of C9_1235: {found}',
          file=sys.stderr)
    if found:
        u, t = found
        phi = [((u * x + t) % 9) for x in range(9)]
        h = [9 + phi[x] for x in range(9)] + [phi[x] for x in range(9)]
        print(f'#   induced a_x -> b_phi(x) is an anti-automorphism of D: {is_iso_map(D, Dop, h)}', file=sys.stderr)


if __name__ == '__main__':
    main()
