#!/usr/bin/env python3
"""Audit E (2026-10-14), HYP-9279(B) 'no poset analogue' sentence.  For every labelled poset on n <= 4 points (and the
two named examples), uniform random linear extension: which elements that are incomparable to at least one other element
have NO partner w with 1/3 <= P(v before w) <= 2/3?"""
from itertools import permutations, product
from fractions import Fraction as F

def lin_ext(n, rel):
    return [p for p in permutations(range(n)) if all(p.index(a) < p.index(b) for a, b in rel)]

def report(n, rel, name=None):
    L = lin_ext(n, rel); H = len(L)
    pos = [{v: k for k, v in enumerate(p)} for p in L]
    comparable = lambda a, b: all(q[a] < q[b] for q in pos) or all(q[a] > q[b] for q in pos)
    lonely = []
    for v in range(n):
        if all(comparable(v, w) for w in range(n) if w != v):
            continue
        ok = any(F(1, 3) <= F(sum(1 for q in pos if q[v] < q[w]), H) <= F(2, 3) for w in range(n) if w != v)
        if not ok:
            lonely.append(v)
    if name:
        probs = {(a, b): str(F(sum(1 for q in pos if q[a] < q[b]), H)) for a in range(n) for b in range(n) if a < b}
        print(f"{name}: {H} linear extensions; P(a before b): {probs}; non-comparable-to-all elements without a balanced partner: {lonely}")
    return lonely

# 2 + 1: elements 0 = a (isolated), 1 = b < 2 = c
report(3, [(1, 2)], "2+1 = {a} + {b<c}")
# 3 + 1: elements 0 = y (isolated), 1 < 2 < 3 chain
report(4, [(1, 2), (2, 3), (1, 3)], "3+1 = {y} + {x<u1<u2}")
# all labelled posets n <= 4
for n in range(2, 5):
    pairs = [(a, b) for a in range(n) for b in range(n) if a != b]
    cnt = bad = 0; first = None
    for bits in product((0, 1), repeat=len(pairs)):
        rel = [pairs[k] for k in range(len(pairs)) if bits[k]]
        S = set(rel)
        if any((b, a) in S for a, b in S): continue                       # antisymmetric
        if any((a, c) not in S for a, b in S for b2, c in S if b == b2 and a != c): continue   # transitive
        cnt += 1
        if lonely := report(n, rel):
            bad += 1
            if first is None: first = (rel, lonely)
    print(f"n={n}: labelled posets {cnt}; posets with an element (incomparable to something) lacking a balanced partner: {bad}; first: {first}")
