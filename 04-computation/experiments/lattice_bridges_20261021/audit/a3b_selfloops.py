#!/usr/bin/env python3
"""AUDIT: the classes that are strict local maxima of H against distinct neighbour classes but not in the
'hypercube' convention (a flip that lands in the same multiplier class counts as a tie): list the self-loop flip and
the multiplier that realises it, and the H values of the distinct neighbour classes.  Reads a3_classes_n.txt.
Output: a3b_selfloops.out"""
import json, os
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'a3b_selfloops.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
for n in range(7, 26, 2):
    f = os.path.join(HERE, f'a3_classes_{n}.txt')
    if not os.path.exists(f): continue
    rows = [json.loads(l) for l in open(f)]
    H = {tuple(r['S']): r['H'] for r in rows}
    m = (n - 1) // 2; U = [a for a in range(1, n) if gcd(a, n) == 1]
    for r in rows:
        S = tuple(r['S']); nb = [tuple(x) for x in r['nbrs']]
        others = [y for y in nb if y != S]
        if len(others) == len(nb) or not all(H[S] > H[y] for y in others): continue
        # self-loop flips
        loops = []
        for k, y in zip(range(1, m + 1), nb):
            if y == S:
                T = set(S)
                if k in T: T.discard(k); T.add(n - k)
                else: T.discard(n - k); T.add(k)
                a = next(a for a in U if {(a * s) % n for s in S} == T)
                loops.append((k, a))
        say(f"n={n}: S={list(S)} H={H[S]} dR={r['dR']}: self-loop flips (class {{k,n-k}}, multiplier a with aS = flipped S): {loops}; "
            f"distinct neighbour classes {len(set(others))}, max neighbour H = {max(H[y] for y in others)}")
out.close()
