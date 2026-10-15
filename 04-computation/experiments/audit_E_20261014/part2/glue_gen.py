#!/usr/bin/env python3
"""Audit E part 2 (2026-10-15): own generator for the HYP-9279 'glued extremal blocks' family, written independently of
glue_blocks.py.  T = z + (B1 => B2): B1, B2 given as gentourng ascii strings (upper triangle row by row, 1 iff i -> j),
all arcs B1 -> B2, and z attached by EVERY pattern, keeping only strong T (z -> some vertex of B1, some vertex of B2 -> z).
Output: digraph6 lines (own encoder; vertex 0 = z) for strong_dfs_balance (DFS-based exact delta).
Usage: glue_gen.py B1 B2      or      glue_gen.py --ascii STRING [STRING ...]   (convert given ascii tournaments)"""
import sys

def dec(s):
    n = int(round((1 + (1 + 8 * len(s)) ** 0.5) / 2)); assert n * (n - 1) // 2 == len(s)
    A = [[0] * n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1': A[i][j] = 1
            else: A[j][i] = 1
            k += 1
    return A

def d6(A):
    n = len(A); bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    return '&' + chr(n + 63) + ''.join(chr(63 + int(''.join(map(str, bits[k:k + 6])), 2)) for k in range(0, len(bits), 6))

BLOCKS = {"S4": "100111", "M5": "1100111111", "M6a": "110001111111101", "M6b": "111001011111111",
          "M7": "111000101111111111101", "M8": "1100011111111111111100111111"}

if __name__ == '__main__':
    if sys.argv[1] == '--ascii':
        for s in sys.argv[2:]: print(d6(dec(s)))
        sys.exit()
    B1, B2 = dec(BLOCKS.get(sys.argv[1], sys.argv[1])), dec(BLOCKS.get(sys.argv[2], sys.argv[2]))
    n1, n2 = len(B1), len(B2); n = 1 + n1 + n2
    for pat in range(1 << (n1 + n2)):
        A = [[0] * n for _ in range(n)]
        for i in range(n1):
            for j in range(n1): A[1 + i][1 + j] = B1[i][j]
            for j in range(n2): A[1 + i][1 + n1 + j] = 1
        for i in range(n2):
            for j in range(n2): A[1 + n1 + i][1 + n1 + j] = B2[i][j]
        for x in range(1, n):
            if pat >> (x - 1) & 1: A[0][x] = 1
            else: A[x][0] = 1
        if not any(A[0][x] for x in range(1, 1 + n1)) or not any(A[x][0] for x in range(1 + n1, n)): continue
        print(d6(A))
