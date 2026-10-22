#!/usr/bin/env python3
"""AUDIT pass 2, item 3 (side check of THM-161 'Key structural fact 5: H-class ordering follows c_p'): compare the
ordering of circulant classes by H with the ordering by the Hamiltonian-cycle count c_n, n = 19, 23 (a12_all*_hk.out).
Output: a12c_order.out"""
import os, re
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'a12c_order.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
    for k, i in enumerate(o): r[i] = k
    return r
for n in (19, 23):
    rows = []
    for line in open(os.path.join(HERE, f'a12_all{n}_hk.out')):
        d = dict((k, int(v)) for k, v in re.findall(r'(H|c\d+)=(\d+)', line)); rows.append(d)
    H = [d['H'] for d in rows]; C = [d['c%d' % n] for d in rows]
    rH, rC = ranks(H), ranks(C); N = len(H); mH = sum(rH) / N; mC = sum(rC) / N
    sp = sum((a - mH) * (b - mC) for a, b in zip(rH, rC)) / (sum((a - mH) ** 2 for a in rH) * sum((b - mC) ** 2 for b in rC)) ** 0.5
    inv = sum(1 for i in range(N) for j in range(N) if H[i] > H[j] and C[i] < C[j])
    say(f"n={n}: {N} classes; Spearman(H, c_{n}) = {sp:.4f}; discordant pairs = {inv} of {N*(N-1)//2}; argmax H = argmax c_{n}: {H.index(max(H)) == C.index(max(C))}")
out.close()
