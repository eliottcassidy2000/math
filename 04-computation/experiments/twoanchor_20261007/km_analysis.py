#!/usr/bin/env python3
"""Karlin-McGregor product law check from km_exponent_n20000.out: q_R(T) / q_1(T)^(R(R+1)/2)."""
import ast, math
rows = {}
for line in open('km_exponent_n20000.out'):
    if line.startswith('R='):
        R = int(line[2]); lst = ast.literal_eval(line[line.index('['):line.index(']')+1])
        rows[R] = dict(lst)
N = 20000
for R in (2, 3):
    e = R*(R+1)//2
    out = []
    for T in sorted(rows[1]):
        q1, qR = rows[1][T], rows[R][T]
        if qR*N >= 20: out.append((T, round(qR/q1**e, 2), int(qR*N)))
    print(f"R={R}: q_R / q_1^{e} at (T, ratio, events):", out)
for R in (1, 2, 3):
    T = sorted(rows[R]); sl = []
    for a, b in zip(T, T[1:]):
        if rows[R][b]*N >= 20: sl.append(round(-math.log(rows[R][b]/rows[R][a])/math.log(b/a), 3))
    print(f"R={R} local slopes:", sl)
