#!/usr/bin/env python3
"""Verifier (lane digits), check V12: values of the one-mark words in a, b, c over a point (own DFS).
a: (h-e+B, g, 2B-g, g)   b: (h-e+B, h-e, B+h, h)   c: (h, 0, 2B, g);  start (1, 1, 1, 1).
These are the letters Z, J, O of the O/Z/J source language (THM-4904, Session 197/198), the point being the state
before the rooted edge (c or b applied to the point gives (1, 0, 2, 1))."""
import sys, time
import numpy as np
from numba import njit


@njit(cache=True)
def run(X):
    seen = np.zeros(X + 1, np.uint8)
    # explicit stack of (h, e, B, g, next_child)
    cap = 4096
    st = np.zeros((cap, 5), np.int64)
    st[0, 0] = 1; st[0, 1] = 1; st[0, 2] = 1; st[0, 3] = 1; st[0, 4] = 0
    top = 0
    nodes = 1
    seen[1] = 1
    while top >= 0:
        h = st[top, 0]; e = st[top, 1]; B = st[top, 2]; g = st[top, 3]; c = st[top, 4]
        if c == 3 or h - e + B > X:
            top -= 1
            continue
        st[top, 4] = c + 1
        if c == 0:
            nh = h - e + B; ne = g; nB = 2 * B - g; ng = g
        elif c == 1:
            nh = h - e + B; ne = h - e; nB = B + h; ng = h
        else:
            nh = h; ne = 0; nB = 2 * B; ng = g
            if nh + nB > X:          # nothing new below a sink step
                continue
        if c == 0 and nh == h and ne == e and nB == B and ng == g:
            continue                 # a fixes the state (mark is a sink): no new node
        seen[nh] = 1
        nodes += 1
        top += 1
        st[top, 0] = nh; st[top, 1] = ne; st[top, 2] = nB; st[top, 3] = ng; st[top, 4] = 0
    return seen, nodes


X = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
t0 = time.time()
seen, nodes = run(X)
miss = [n for n in range(1, X + 1, 2) if not seen[n]]
even = [n for n in range(2, X + 1, 2) if seen[n]]
print(f"values <= {X}: nodes {nodes}; odd numbers not attained: {len(miss)}: {miss[:60]}; even values attained: {len(even)}   [{time.time() - t0:.1f} s]")
ref = [7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 583, 637, 821, 1067, 1339, 1861, 2647]
print("equal to the 29 values of THM-4904 / depth_K7_1e6.txt:", miss == ref)
