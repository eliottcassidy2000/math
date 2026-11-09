#!/usr/bin/env python3
# audit K1: explicit check of 275 = 5^2 * 11 at order 13 as B5 => C11 => B5^op, and of 121 = 11^2 at order 10 as
# B11 => B11^op (B5: the strong 4-tournament, H = 5; C11 = B11: the unique order-5 class with H = 11, SC and strong).
# Codes are THM-4922 census codes (bit b for the b-th pair i<j, set iff i -> j); everything else is recomputed here.
import sys, itertools
def decode(code, n):
    A = [[0] * n for _ in range(n)]; b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
            b += 1
    return A
def H(A):
    n = len(A); f = [[0] * n for _ in range(1 << n)]
    for v in range(n): f[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            x = f[S][v]
            if not x: continue
            for u in range(n):
                if not (S >> u & 1) and A[v][u]: f[S | 1 << u][u] += x
    return sum(f[(1 << n) - 1])
def converse(A): return [[A[j][i] for j in range(len(A))] for i in range(len(A))]
def join(*blocks):            # blocks[0] => blocks[1] => ...
    n = sum(len(B) for B in blocks); A = [[0] * n for _ in range(n)]; off = []; o = 0
    for B in blocks: off.append(o); o += len(B)
    for k, B in enumerate(blocks):
        for i in range(len(B)):
            for j in range(len(B)): A[off[k] + i][off[k] + j] = B[i][j]
            for l in range(k + 1, len(blocks)):
                for j in range(len(blocks[l])): A[off[k] + i][off[l] + j] = 1
    return A
def self_converse(A):          # backtracking anti-isomorphism search
    n = len(A); out = [sum(r) for r in A]; p = [None] * n; used = [False] * n
    def bt(k):
        if k == n: return True
        for w in range(n):
            if used[w] or out[w] != n - 1 - out[k]: continue
            if all(A[j][k] == A[w][p[j]] for j in range(k)):
                p[k] = w; used[w] = True
                if bt(k + 1): return True
                used[w] = False
        return False
    return bt(0)
B5 = decode(int(sys.argv[1]), 4); C11 = decode(int(sys.argv[2]), 5)
print("H(B5) =", H(B5), " H(C11) =", H(C11), " C11 self-converse:", self_converse(C11))
T13 = join(B5, C11, converse(B5)); T10 = join(C11, converse(C11))
print("order", len(T13), ": H(B5 => C11 => B5^op) =", H(T13), " self-converse:", self_converse(T13))
print("order", len(T10), ": H(B11 => B11^op) =", H(T10), " self-converse:", self_converse(T10))
