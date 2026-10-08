#!/usr/bin/env python3
"""Shortest absorbing coin sequences from (0, 1) for the px+1 pair chain, and the exact number of absorbing
sequences of that length (an exact lower bound q_p >= count * 2^-L).  p = 3 .. 31 odd."""
from fractions import Fraction as Fr
def par(e): return (e.numerator * pow(e.denominator, -1, 2)) % 2
def step(p, k, e, b):
    s = par(e)
    if s == 0:
        return (k, e / 2) if b == 0 else (k, (p * e + 1 - Fr(p) ** k) / 2)
    return (k + 1, (p * e + 1) / 2) if b == 0 else (k - 1, (e - Fr(p) ** (k - 1)) / 2)
def shortest(p, depth=26, cap=1500000):
    frontier = {(0, Fr(1)): 1}
    for L in range(1, depth + 1):
        new = {}; hits = 0
        for (k, e), mult in frontier.items():
            for b in (0, 1):
                k2, e2 = step(p, k, e, b)
                if k2 == 0 and e2 == 0: hits += mult; continue
                if abs(e2) > Fr(p) ** max(k2, 0) * 10**4: continue      # escape region; cannot be first-shortest
                new[(k2, e2)] = new.get((k2, e2), 0) + mult
        if hits: return L, hits
        frontier = new
        if len(frontier) > cap: return None, None
    return None, None
for p in range(3, 33, 2):
    L, h = shortest(p)
    print(f"p={p:2d}: shortest absorption from (0,1) at Terras time {L}, {h} absorbing sequences of that length -> q_p >= {h}/2^{L} = {h/2**L:.5f}" if L else f"p={p}: none within depth or state cap", flush=True)
