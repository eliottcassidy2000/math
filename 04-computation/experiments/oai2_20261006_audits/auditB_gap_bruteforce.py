#!/usr/bin/env python3
"""Third, solver-free method for the fully gap-determined n = 2 game:
enumerate ALL subsets E of the lex-positive gap vectors of [t]^2 (bitmasks, numpy) and test
sum-freeness inside the box and hitting of every binary subgrid. t = 3 (12 gaps) and t = 4 (24 gaps).
Here the constraints are derived from the sum-free reformulation (d1, d2, d1+d2 all in the box,
d1 = d2 allowed) -- a different derivation from the point-triple encoding of auditB_n2_games.py."""
import itertools, numpy as np


def gaps(t):
    return [(0, c) for c in range(1, t)] + [(a, c) for a in range(1, t) for c in range(-(t - 1), t)]


def run(t):
    G = gaps(t)
    gi = {g: i for i, g in enumerate(G)}
    tri = set()
    for d1 in G:
        for d2 in G:
            s = (d1[0] + d2[0], d1[1] + d2[1])
            if s in gi:
                tri.add(frozenset((gi[d1], gi[d2], gi[s])))
    trimasks = sorted({sum(1 << i for i in f) for f in tri})
    hit = set()
    pairs = list(itertools.combinations(range(t), 2))
    for a, b in itertools.combinations(range(t), 2):
        for p in pairs:
            for q in pairs:
                pts = [(a, p[0]), (a, p[1]), (b, q[0]), (b, q[1])]
                m = 0
                for x, y in itertools.combinations(pts, 2):
                    if y < x:
                        x, y = y, x
                    m |= 1 << gi[(y[0] - x[0], y[1] - x[1])]
                hit.add(m)
    hitmasks = sorted(hit)
    nG = len(G)
    total = 0
    good_examples = []
    chunk = 1 << 22
    for start in range(0, 1 << nG, chunk):
        m = np.arange(start, min(start + chunk, 1 << nG), dtype=np.int64)
        ok = np.ones(len(m), dtype=bool)
        for c in trimasks:
            ok &= (m & c) != c
        for h in hitmasks:
            ok &= (m & h) != 0
        cnt = int(ok.sum())
        total += cnt
        if cnt and len(good_examples) < 3:
            for x in m[ok][:3]:
                good_examples.append([G[i] for i in range(nG) if (int(x) >> i) & 1])
    return nG, len(trimasks), len(hitmasks), total, good_examples


for t in (3, 4):
    nG, nt, nh, total, ex = run(t)
    print(f"t={t}: {nG} gap vectors, {nt} Schur-triple masks, {nh} distinct hitting masks -> "
          f"{total} winning sets E" + (f"; e.g. {ex[0]}" if ex else ""), flush=True)
