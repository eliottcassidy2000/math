#!/usr/bin/env python3
"""d87_explore.py - order-15 'unequal doubles': A = near-regular 8-tournament, B' = QR_7, all arcs A -> B' except
back arcs b_y -> a_{f(y)} for an injection f: Z_7 -> A covering the four A-vertices of in-degree 3 (so delta0 >= 4).
usage: python3 d87_explore.py N seed > list"""
import random, sys
N = int(sys.argv[1]); random.seed(int(sys.argv[2]))
Q = {1, 2, 4}
def rand_nearreg(q):
    while True:
        out = [[False]*q for _ in range(q)]
        for i in range(q):
            for j in range(i+1, q):
                if random.random() < 0.5: out[i][j] = True
                else: out[j][i] = True
        if sorted(sum(r) for r in out) == [3]*4 + [4]*4: return out
seen = set()
for k in range(N):
    A = rand_nearreg(8)
    needy = [x for x in range(8) if sum(A[x]) == 4]        # score 4 -> in-degree 3
    others = [x for x in range(8) if x not in needy]
    tgt = needy + random.sample(others, 3)
    random.shuffle(tgt)                                      # f(y) = tgt[y]
    n = 15
    def arc(i, j):                                           # i < j; A = 0..7, B' = 8..14 (b_y = 8 + y)
        if j < 8: return A[i][j]
        if i >= 8: return ((j - i) % 7) in Q
        y = j - 8
        return tgt[y] != i                                   # a_i -> b_y unless back arc b_y -> a_i
    code = ''.join('1' if arc(i, j) else '0' for i in range(n) for j in range(i+1, n))
    if code in seen: continue
    seen.add(code); print(n, code, f"D87_{k}")
