#!/usr/bin/env python3
"""AUDIT of THM-4751 (1): brute-force BFS distance between the classes of R_q and P_q in the circulant class graph
(vertices = connection sets up to multipliers, edges = one difference-class flip) for q = 7, 11, 19, 23, 31,
compared with min over multipliers of the set distance and with B_q.  No H values are needed.
Output: a1b_classgraph.out"""
import os
from math import gcd
from collections import deque
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'a1b_classgraph.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
for q in (7, 11, 19, 23, 31):
    m = (q - 1) // 2; U = [a for a in range(1, q) if gcd(a, q) == 1]
    # encode a set by bitmask b: bit k-1 set iff q-k in S
    def tomask(S): return sum(1 << (k - 1) for k in range(1, m + 1) if (q - k) in S)
    def toset(b): return {(q - k) if (b >> (k - 1)) & 1 else k for k in range(1, m + 1)}
    canon = {}
    for b in range(1 << m):
        if b in canon: continue
        S = toset(b); orb = {tomask({(a * s) % q for s in S}) for a in U}
        c = min(orb)
        for x in orb: canon[x] = c
    R = canon[0]; QR = {x * x % q for x in range(1, q)}; P = canon[tomask(QR)]
    dist = {R: 0}; dq = deque([R])
    while dq:
        x = dq.popleft()
        for k in range(m):
            y = canon[x ^ (1 << k)]
            if y not in dist: dist[y] = dist[x] + 1; dq.append(y)
    md = min(len({(a * s) % q for s in QR} - set(range(1, m + 1))) for a in U)
    A = len(QR & set(range(1, m + 1))); B = m - A
    say(f"q={q}: classes={len(set(canon.values()))}, BFS d([R],[P]) = {dist[P]}, min over multipliers = {md}, B_q = {B}, eccentricity of [R] = {max(dist.values())}")
out.close()
