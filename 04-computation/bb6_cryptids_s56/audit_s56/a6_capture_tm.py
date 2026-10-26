#!/usr/bin/env python3
"""Referee: run the A6 TM from the family configuration J(0,4,45) and time successive returns to J(0,v0,4v0+29),
v0 even: TM steps per period and tape length (translated cycler would have constant period)."""
import sys, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from tmcheck import parse_tm
from tm_rules_check import A6, a6_enc, a6_parse
tab = parse_tm(A6)
st, L, R = a6_enc((0, 4, 45))
N = 1 << 16
tape = bytearray(N); off = N // 2
for i, ch in enumerate(L): tape[off + i] = int(ch)
head = off + len(L)
for i, ch in enumerate(R): tape[head + i] = int(ch)
state = st; lo, hi = off, head + len(R)
last = 0; hits = []
for step in range(1, 30_000_000):
    tr = tab[(state, tape[head])]
    if tr is None:
        print('HALT', step); break
    w, mv, ns = tr; tape[head] = w; head += mv; state = ns
    lo = min(lo, head); hi = max(hi, head)
    if state == 0 and tape[head] == 1 and tape[head - 1] == 1:
        a = lo
        while tape[a] == 0: a += 1
        b = hi
        while tape[b] == 0: b -= 1
        left = ''.join(map(str, tape[a:head])); right = ''.join(map(str, tape[head:b + 1]))
        p = a6_parse(0, left, right)
        if p and p[0] == 0 and p[1] % 2 == 0 and p[2] == 4 * p[1] + 29:
            hits.append((step, step - last, p, b - a + 1)); last = step
            if len(hits) >= 12: break
for h in hits: print('step', h[0], 'period', h[1], 'config', h[2], 'tape length', h[3])
