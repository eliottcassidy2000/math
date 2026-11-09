#!/usr/bin/env python3
"""d8_explore.py - doubles D(B) of near-regular 8-tournaments (order 16, minimum semidegree 4): random sample of B,
output "16 code tag" lines for mdefect/tprop.  usage: python3 d8_explore.py N seed > list"""
import random, sys
N = int(sys.argv[1]); random.seed(int(sys.argv[2]))
def rand_nearreg(q):
    while True:
        out = [[False]*q for _ in range(q)]
        for i in range(q):
            for j in range(i+1, q):
                if random.random() < 0.5: out[i][j] = True
                else: out[j][i] = True
        sc = sorted(sum(r) for r in out)
        if sc == [3]*4 + [4]*4: return out
def double(B):
    q = len(B); n = 2*q
    def arc(i, j):
        if j < q: return B[i][j]
        if i >= q: return B[i-q][j-q]
        return (j - q) != i          # a_x -> b_y unless y == x (back arc b_x -> a_x)
    return ''.join('1' if arc(i, j) else '0' for i in range(n) for j in range(i+1, n))
seen = set()
for k in range(N):
    B = rand_nearreg(8)
    c = double(B)
    if c in seen: continue
    seen.add(c)
    print(16, c, f"D8_{k}")
