import json, re, sys
exec(open('identify_new.py').read().split('results = []')[0].replace('print(f"{len(qs)} primes', 'pass  # print(f"{len(qs)} primes'))
from cone_reduce import Y_H, tables
cheap = [2, 4, 5, 6, 9, 18]
Yv = {i: [] for i in cheap}
for q in qs:
    T = tables(q)
    for i in cheap: Yv[i].append(Y_H(q, G[i]['edges'], T=T)[0])
for i in cheap:
    diff = [a - b for a, b in zip([S[i][q] for q in qs], Yv[i])]
    print(f"H{i:02d}: S_H - Y_H in base span: {rank(base + [diff]) == rB}; Y_H new part independent: {rank(base + [Yv[i]]) > rB}")
r = rank(base + [Yv[i] for i in cheap]) - rB
print(f"rank of the six cheap Y_H modulo base: {r}; do they span the new space: {rank(base + [Yv[i] for i in cheap] + newcols) == rank(base + [Yv[i] for i in cheap])}")
json.dump({'qs': qs, 'Y': {str(i): Yv[i] for i in cheap}}, open('cheap_Y_values.json', 'w'))
