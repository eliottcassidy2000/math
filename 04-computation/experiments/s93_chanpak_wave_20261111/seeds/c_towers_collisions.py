#!/usr/bin/env python3
"""S93 lane seeds: are the tower slopes distinct?  Pure python enumeration of the tower states (r, M, h, E) with
M <= 2^K (level map of c_towers.py), number of distinct states, of distinct slopes, and the collisions
(two different states with the same sign and slope)."""
import sys
from collections import defaultdict

ALPHA = (0, 3, 5)
GAMMA = (0, 5, 13)


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    X = 1 << K
    states = {}
    stack = []
    g = 0
    while (1 << g) <= X:
        stack.append(((1, 1 << g, 1, 1 if g == 0 else 0), ("fan", g)))
        g += 1
    for delta in (1, 2):
        m = 0
        while ALPHA[delta] << m <= X:
            stack.append(((-1, ALPHA[delta] << m, 1, 1 if m == 0 else 0), ("pair", delta, m)))
            m += 1
    nstates = 0
    while stack:
        st, hist = stack.pop()
        nstates += 1
        if st not in states:
            states[st] = hist
        else:
            continue                      # same state: same descendants
        r, M, h, E = st
        for delta in (1, 2):
            al, ga = ALPHA[delta], GAMMA[delta]
            i = 1
            while True:
                m = 0
                any_m = False
                while True:
                    lam = ga * (1 << i) * 3 ** m - (1 if i == 1 else 0)
                    if lam * M > X:
                        break
                    any_m = True
                    M2 = al * (1 << m) * (h + r) + 2 * ga * (1 << (i - 1)) * 3 ** m * M - (E if i == 1 else 0)
                    if M2 <= X:
                        a = i + delta + m
                        st2 = (-r, M2, r + (1 << a) * M, al * (1 << (i - 1)) * M + ((r - h) if m == 0 else 0))
                        stack.append((st2, hist + ((delta, i, m),)))
                    m += 1
                if not any_m:
                    break
                i += 1
    byslope = defaultdict(list)
    for st in states:
        byslope[(st[0], st[1])].append(st)
    coll = {k: v for k, v in byslope.items() if len(v) > 1}
    print(f"X = 2^{K}: distinct tower states {len(states)}; distinct (sign, slope) {len(byslope)}; (sign, slope) classes with more than one state: {len(coll)}; max class size {max(len(v) for v in byslope.values())}")
    ex = sorted(coll.items(), key=lambda kv: kv[0][1])[:8]
    for (sg, M), v in ex:
        print(f"  sign {sg:+d} slope {M}: " + " | ".join(f"h={s[2]} E={s[3]} via {states[s]}" for s in v))
    # states that differ only in E, or only in h
    sameMh = defaultdict(set)
    for st in states:
        sameMh[(st[0], st[1], st[2])].add(st[3])
    print(f"  distinct (sign, slope, h): {len(sameMh)}")


if __name__ == "__main__":
    main()
