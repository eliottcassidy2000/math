"""Audit (A4) and HYP-9428 (1): window-restricted catalytic debt on root lattices with roots as unit charges.
Independent enumeration of each lattice; Knuth-Dijkstra solver (knuth_debt.py).  g_W = 0 on a window is a valid
upper bound g <= g_W = 0, hence exact.  Also reducible controls A2+A1, A2+A2."""
import itertools, sys, time
import numpy as np
from knuth_debt import knuth_debt, INF

def run(name, pts, G, s, expect_roots=None):
    t0 = time.time()
    g, Q, _ = knuth_debt(pts, G, s)
    nroots = int((Q == 1).sum()) if Q[Q > 0].min() == 1 else int((Q == Q[Q > 0].min()).sum())
    unreach = int((g >= INF).sum())
    gm = int(g[g < INF].max())
    Wmax = int(Q.max())
    print(f"{name:10s} roots {nroots:4d}{'' if expect_roots is None else (' (expect %d)' % expect_roots)}  window Q<={Wmax}: "
          f"{len(pts):6d} points, unreachable {unreach}, max g {gm}   [{time.time()-t0:.1f}s]", flush=True)
    return g, Q

def box(d, R): return itertools.product(range(-R, R + 1), repeat=d)

def Dn(n, W):
    R = int((2 * W) ** 0.5)
    return [p for p in box(n, R) if sum(p) % 2 == 0 and sum(a*a for a in p) <= 2 * W]
def An(n, W):
    R = int((2 * W) ** 0.5)
    return [p for p in box(n + 1, R) if sum(p) == 0 and sum(a*a for a in p) <= 2 * W]
def E8(W):   # doubled coordinates: one parity, sum = 0 mod 4, Q = y.y/8
    R = int((8 * W) ** 0.5)
    out = []
    for par in (0, 1):
        vals = [a for a in range(-R, R + 1) if a % 2 == par]
        def rec(pre, budget):
            if len(pre) == 8:
                if sum(pre) % 4 == 0: out.append(tuple(pre))
                return
            for a in vals:
                if a*a <= budget: rec(pre + [a], budget - a*a)
        rec([], 8 * W)
    return out
def E7(W):   # doubled: one parity, sum 0, Q = y.y/8   (E7 = A7 u (A7 + (1/2^4, -1/2^4)))
    R = int((8 * W) ** 0.5)
    out = []
    for par in (0, 1):
        vals = [a for a in range(-R, R + 1) if a % 2 == par]
        def rec(pre, budget):
            if len(pre) == 8:
                if sum(pre) == 0: out.append(tuple(pre))
                return
            for a in vals:
                if a*a <= budget: rec(pre + [a], budget - a*a)
        rec([], 8 * W)
    return out
def E6(W):   # doubled (a1..a6; t): one parity, sum a = 0, Q = (a.a + 2 t^2)/8   (E6 = A5+A1 u glue)
    R = int((8 * W) ** 0.5)
    out = []
    for par in (0, 1):
        vals = [a for a in range(-R, R + 1) if a % 2 == par]
        for t in vals:
            if 2*t*t > 8*W: continue
            def rec(pre, budget):
                if len(pre) == 6:
                    if sum(pre) == 0: out.append(tuple(pre) + (t,))
                    return
                for a in vals:
                    if a*a <= budget: rec(pre + [a], budget - a*a)
            rec([], 8 * W - 2*t*t)
    return out

I = lambda d: np.eye(d, dtype=np.int64)
which = sys.argv[1:] or ['theorem', 'larger', 'new', 'reducible']
if 'theorem' in which:
    print("== THM-4824 (A4) windows ==")
    run("D3=A3", Dn(3, 12), I(3), 2, 12)
    run("D4", Dn(4, 7), I(4), 2, 24)
    run("A4", An(4, 6), I(5), 2, 20)
    run("A5", An(5, 4), I(6), 2, 30)
    run("E8", E8(2), I(8), 8, 240)
if 'larger' in which:
    print("== larger windows ==")
    run("D3=A3", Dn(3, 24), I(3), 2)
    run("D4", Dn(4, 10), I(4), 2)
    run("D5", Dn(5, 6), I(5), 2, 40)
    run("A4", An(4, 9), I(5), 2)
    run("A5", An(5, 6), I(6), 2)
if 'new' in which:
    print("== lattices HYP-9428 (1) calls OPEN ==")
    run("A6", An(6, 4), I(7), 2, 42)
    run("A7", An(7, 3), I(8), 2, 56)
    G6 = np.diag([1, 1, 1, 1, 1, 1, 2]).astype(np.int64)
    run("E6", E6(4), G6, 8, 72)
    run("E7", E7(3), I(8), 8, 126)
    run("E8", E8(3), I(8), 8, 240)
if 'reducible' in which:
    print("== reducible controls (A2 = Z[w] block [[2,-1],[-1,2]]/2) ==")
    A2 = np.array([[2, -1], [-1, 2]])
    G = np.zeros((3, 3), dtype=np.int64); G[:2, :2] = A2; G[2, 2] = 2
    W = 40
    pts = [p for p in box(3, 13) if p[0]**2 - p[0]*p[1] + p[1]**2 + p[2]**2 <= W]
    g, Q = run("A2+A1", pts, G, 2, 8)
    G = np.zeros((4, 4), dtype=np.int64); G[:2, :2] = A2; G[2:, 2:] = A2
    W = 16
    pts = [p for p in box(4, 9) if p[0]**2 - p[0]*p[1] + p[1]**2 + p[2]**2 - p[2]*p[3] + p[3]**2 <= W]
    g2, Q2 = run("A2+A2", pts, G, 2, 12)
