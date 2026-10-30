#!/usr/bin/env python3
"""Direct integer content vs independently accumulated epoch envelope."""
from fractions import Fraction
import json
from pathlib import Path

def vp(n):
    if not n:
        return 1000000
    e = 0
    while n % 5 == 0:
        e += 1
        n //= 5
    return e

def words(A):
    if A == 0:
        yield ()
    for a in range(1, A + 1):
        for tail in words(A - a):
            yield ((a,) + tail)

def check(w, E):
    q = 5 ** E
    D = 8 * q
    N = 5 * q - 5
    P = Q = 1
    B = 0
    S = 0
    C = 0
    entered = False
    for a in w:
        P, Q, B = (3 * P, 2 ** a * Q, 3 * B + Q)
        active = pow(Q, -1, 5) * P % 5 in (1, 4)
        h = min(vp(B + Q - P), E)
        S += min(vp(P - Q), vp(B), E)
        if a > 1:
            entered = True
            C += (D * max(h, min(vp(2 ** a - 2), E)) if active else 0) - N * a
        elif entered:
            C += -(5 * q // 5 ** h - 5) if active else -N
    if not D * S - N * sum(w) <= C:
        raise RuntimeError((w, E, S, C))
    if not C <= 0:
        raise RuntimeError((w, E, S, C))
    return S
count = 0
for A in range(1, 15):
    for w in words(A):
        for E in (1, 2, 3, 6, 7, 8):
            check(w, E)
        count += 1
for E in range(1, 9):
    t = 2 * 5 ** (E - 1)
    S = sum((t // (2 * 5 ** (e - 1)) for e in range(1, E + 1)))
    if not 8 * 5 ** E * S == (5 * 5 ** E - 5) * t:
        raise RuntimeError('verification failed at original line 43')
for w in ((2,) + (1,) * 92 + (3,), (3,) + (1,) * 16, (2, 2), (1,) * 30 + (2,)):
    for E in range(1, 9):
        check(w, E)
r = {'exhaustive_words': count, 'maximum_clock': 14, 'precisions': [1, 2, 3, 6, 7, 8], 'equality_controls_through_E': 8, 'status': 'PASS; direct integer words are controls, not proof of all words'}
p = Path(__file__).with_name('word_controls.json')
p.write_text(json.dumps(r, indent=2) + '\n')
print(p.read_text())
