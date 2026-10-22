#!/usr/bin/env python3
"""(c) Paper 1: root-subsystem multiplicities by brute-force frame counting.
m_R(Y) = #(ordered tuples of roots with Gram = Cartan(Y)) / (|W(Y)| |Aut(Dynkin Y)|).
For irreducible R we fix v1 = root 0 and multiply by #roots (W(R) is transitive on roots)."""
import sys, time
sys.path.insert(0, '.')
from math import comb
import numpy as np
from rootsys import roots, gram, cartan, count_tuples

WA = {'A5': 720*2, 'D5': 1920*2, 'D4': 192*6, 'A4': 120*2}
STAR = np.ones((5, 5), dtype=np.int64) + np.eye(5, dtype=np.int64)   # Gram 2 on diag, 1 off diag

def frames(typ, k, G):
    R, sc = roots(typ, k)
    IP = gram(R, sc)
    return R.shape[0] * count_tuples(IP, G, first_candidates=[0])

rows = []
t0 = time.time()
systems = [('A', n) for n in range(5, 18)] + [('D', n) for n in range(4, 17)] + [('E', 6), ('E', 7), ('E', 8)]
expected = {}
for n in range(5, 18):
    expected[('A', n)] = {'D5': 0, 'A5': comb(n+1, 6), 'D4': 0}
for n in range(4, 17):
    expected[('D', n)] = {'D5': comb(n, 5), 'A5': 32*comb(n, 6), 'D4': comb(n, 4)}
expected[('E', 6)] = {'D5': 27, 'A5': 36, 'D4': 45}
expected[('E', 7)] = {'D5': 378, 'A5': 1344, 'D4': 315}
expected[('E', 8)] = {'D5': 7560, 'A5': 40320, 'D4': 3150}
allok = True
for typ, k in systems:
    res = {}
    for Y, G in [('D5', cartan('D', 5)), ('A5', cartan('A', 5)), ('D4', cartan('D', 4))]:
        f = frames(typ, k, G)
        assert f % WA[Y] == 0
        res[Y] = f // WA[Y]
    fs = frames(typ, k, STAR)
    res['A5(star)'] = fs / 1440
    ok = all(res[Y] == expected[(typ, k)][Y] for Y in ('D5', 'A5', 'D4')) and res['A5(star)'] == res['A5']
    allok &= ok
    rows.append(f"{typ}{k}: m(D5)={res['D5']}, m(A5)={res['A5']}, m(A5 via star frames)={res['A5(star)']:.0f}, m(D4)={res['D4']}   expected {expected[(typ,k)]}  {'OK' if ok else 'MISMATCH'}")
print("\n".join(rows))
print("ALL MULTIPLICITY FORMULAS OK:", allok, f"({time.time()-t0:.1f}s)")

# E7 / E6 / E8 orbit split of A5 (and D5) by the size of the orthogonal root system
def frame_complements(typ, k, G):
    R, sc = roots(typ, k)
    IP = gram(R, sc)
    N = R.shape[0]
    m = G.shape[0]
    hist = {}
    eqm = {v: (IP == v) for v in (-1, 0, 1)}
    def rec(ch):
        t = len(ch)
        if t == m:
            orth = np.ones(N, bool)
            for v in ch: orth &= eqm[0][v]
            c = int(orth.sum())
            hist[c] = hist.get(c, 0) + 1
            return
        cand = np.ones(N, bool) if t else np.zeros(N, bool)
        if t == 0: cand[0] = True
        for i, v in enumerate(ch):
            cand &= eqm[int(G[i, t])][v]
        for v in np.nonzero(cand)[0]:
            rec(ch + [int(v)])
    rec([])
    # frames with v1 fixed -> multiply by N, divide by frames per subsystem
    return {c: cnt * N for c, cnt in hist.items()}

for typ, k in [('E', 6), ('E', 7), ('E', 8)]:
    for Y, G, w in [('A5', cartan('A', 5), 1440), ('D5', cartan('D', 5), 3840)]:
        h = frame_complements(typ, k, G)
        print(f"{typ}{k} {Y}-subsystems split by #roots orthogonal to the subsystem: " +
              ", ".join(f"{c} orth. roots: {v//w}" for c, v in sorted(h.items())))
