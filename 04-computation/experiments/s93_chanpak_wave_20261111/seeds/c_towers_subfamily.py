#!/usr/bin/env python3
"""S93 lane seeds: collisions of slopes inside sub-families of the towers.
Sub-family F(imin): all levels have i >= imin (imin = 2: the E-coordinate never enters; reduced state (r, M, G),
G = h + r:  M+ = alpha 2^m G + gamma 2^i 3^m M,  G+ = 2^(i+delta+m) M).
Bases: fans (1, 2^g) with G = 2, pairs (-1, alpha 2^m) with G = 0.
Counts: distinct reduced states with M <= 2^K, distinct (sign, M), classes with several states."""
import sys
from collections import defaultdict

ALPHA = (0, 3, 5)
GAMMA = (0, 5, 13)


def run(K, imin):
    X = 1 << K
    states = {}
    stack = []
    g = 0
    while (1 << g) <= X:
        stack.append(((1, 1 << g, 2), ("fan", g)))
        g += 1
    for delta in (1, 2):
        m = 0
        while ALPHA[delta] << m <= X:
            stack.append(((-1, ALPHA[delta] << m, 0), ("pair", delta, m)))
            m += 1
    while stack:
        st, hist = stack.pop()
        if st in states:
            continue
        states[st] = hist
        r, M, G = st
        for delta in (1, 2):
            al, ga = ALPHA[delta], GAMMA[delta]
            i = imin
            while True:
                m = 0
                any_m = False
                while True:
                    if ga * (1 << i) * 3 ** m * M > X:
                        break
                    any_m = True
                    M2 = al * (1 << m) * G + ga * (1 << i) * 3 ** m * M
                    if M2 <= X:
                        stack.append(((-r, M2, (1 << (i + delta + m)) * M), hist + ((delta, i, m),)))
                    m += 1
                if not any_m:
                    break
                i += 1
    byslope = defaultdict(list)
    for st in states:
        byslope[(st[0], st[1])].append(st)
    coll = {k: v for k, v in byslope.items() if len(v) > 1}
    print(f"i >= {imin}, X = 2^{K}: distinct reduced states {len(states)}; distinct (sign, slope) {len(byslope)}; classes with several states {len(coll)}; largest class {max(len(v) for v in byslope.values())}")
    for (sg, M), v in sorted(coll.items(), key=lambda kv: kv[0][1])[:6]:
        print(f"    sign {sg:+d} slope {M}: " + " | ".join(f"G={s[2]} via {states[s]}" for s in v))
    return len(states), len(byslope)


if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    for imin in (2, 3):
        run(K, imin)
