# -*- coding: utf-8 -*-
# klein-2026-07-11-S252: THE k=8 DEGREE-3 ROW, EXHAUSTIVE -- the moment ladder's last
# analytic gap verified over the FULL bounded-core box with exact rationals.
#
# CONTEXT: the LRC(14)-S3 ladder (opus-S220 two-moment residue; mac-mini THM-703/705/710)
# reduces to {deg-3 @ k=8} + {deg-2 @ k=9}, higher k inherited via THM-710's eigen-transfer.
# The deg-2 majorant q2 = 1 - (2/3)t + (1/6)t(t-1) FAILS at k=8 by ~0.033.  THIS SCRIPT:
#
# (1) derives the OPTIMAL degree-3 majorant by LP-vertex enumeration (all triples of tight
#     constraints among {phi(1) = 1/3} u {phi(n) = 0 : n = 2..6}; feasibility phi >= g on
#     {0..6}) -- the unique new feasible vertex touches {1,3,6}:
#         phi3(t) = 1 - (2/3) t + (17/90) t(t-1) - (1/45) t(t-1)(t-2),
#     values [1, 1/3, 2/45, 0, 1/15, 1/9, 0]; the deg-3 gain over deg-2 = (1/45)(m3 - m2).
# (2) verifies Phi(F) <= E[phi3(N)] pointwise-majorant style and computes the ROW
#         (2/3) m1 - (17/90) m2 + (1/45) m3  >=  1 - cap_9      (k = 8)
#     EXHAUSTIVELY over ALL 0-anchored 8-cores in [1..W8] (exact factorial moments via the
#     all-integer breakpoint sweep), reporting min margin + argmin.
# (3) same exhaustive treatment for the deg-2 row at k = 9 (mac-mini THM-705's linear
#     pair-correlation requirement): (2/3) m1 - (1/6) m2 >= 1 - cap_10 over [1..W9].
# (4) cross-checks: deg-2 row at k=8 fails (the known +0.033 gap) -- locating exactly
#     where deg-3 rescues it; consec values vs opus-S220/mac-mini printouts.
#
# Convention == kps/mac-mini/opus: sectors [c/7,(c+1)/7), N = #empty among inner {1..6}
# (0-anchor covers sector 0); factorial moments m_r = E[N(N-1)...(N-r+1)].

import sys
from math import gcd, comb
from fractions import Fraction as F
from itertools import combinations

CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91),
       12: F(6, 7), 13: F(1)}

def lcm(a, b):
    return a // gcd(a, b) * b

def moments_exact(E):
    """Exact (p0, p1, m1, m2, m3) for 0-anchored offset set E via all-integer sweep."""
    nz = [e for e in E if e != 0]
    L = 1
    for e in nz:
        L = lcm(L, e)
    D = 7 * L
    pts = set()
    for e in nz:
        step = L // e
        pts.update(range(0, D + 1, step))
    pts.add(0); pts.add(D)
    pts = sorted(pts)
    pn = [0] * 7                    # numerators of P(N = n), units 1/D
    for t1, t2 in zip(pts, pts[1:]):
        s = t1 + t2
        hit = 1
        for e in nz:
            c = (7 * e * s // (2 * D)) % 7
            hit |= 1 << c
        n = 7 - bin(hit).count("1")
        pn[n] += t2 - t1
    p = [F(x, D) for x in pn]
    m1 = sum(F(n) * p[n] for n in range(7))
    m2 = sum(F(n * (n - 1)) * p[n] for n in range(7))
    m3 = sum(F(n * (n - 1) * (n - 2)) * p[n] for n in range(7))
    return p, m1, m2, m3

def lp_vertices():
    """Enumerate all deg-3 majorant vertices: phi = 1 + c1 t + c2 t(t-1) + c3 t(t-1)(t-2),
    3 tight constraints among {phi(1)=1/3} u {phi(n)=0, n=2..6}; feasibility on {0..6}."""
    import itertools
    def phi_of(c1, c2, c3, t):
        return 1 + c1 * t + c2 * t * (t - 1) + c3 * t * (t - 1) * (t - 2)
    cons = [(1, F(1, 3))] + [(n, F(0)) for n in range(2, 7)]
    out = []
    for triple in itertools.combinations(cons, 3):
        # solve the 3x3 linear system for (c1, c2, c3)
        rows = []
        rhs = []
        for t, val in triple:
            rows.append((F(t), F(t * (t - 1)), F(t * (t - 1) * (t - 2))))
            rhs.append(val - 1)
        # Cramer
        def det3(M):
            return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                    - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                    + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
        M = [list(r) for r in rows]
        d = det3(M)
        if d == 0:
            continue
        cs = []
        for col in range(3):
            Mc = [list(r) for r in rows]
            for i in range(3):
                Mc[i][col] = rhs[i]
            cs.append(det3(Mc) / d)
        c1, c2, c3 = cs
        vals = [phi_of(c1, c2, c3, t) for t in range(7)]
        g = [F(1), F(1, 3)] + [F(0)] * 5
        if all(vals[t] >= g[t] for t in range(7)):
            out.append(((c1, c2, c3), tuple(t for t, _ in triple), vals))
    return out

def census_row(k, W, c1, c2, c3, capk1, label):
    """Exhaustively verify 1 + c1 m1 + c2 m2 + c3 m3 <= capk1 over 0-anchored k-cores in [1..W]."""
    worst = None; worst_set = None; consec_val = None
    n_eval = 0
    for combo in combinations(range(1, W + 1), k - 1):
        g0 = 0
        for c in combo:
            g0 = gcd(g0, c)
        if g0 > 1:
            continue
        n_eval += 1
        E = (0,) + combo
        p, m1, m2, m3 = moments_exact(E)
        bound = 1 + c1 * m1 + c2 * m2 + c3 * m3
        margin = capk1 - bound
        phi = p[0] + p[1] / 3
        assert phi <= bound + F(0), f"majorant violated at {E}"  # sanity: Phi <= E[phi(N)]
        if E == tuple(range(k)):
            consec_val = (bound, margin)
        if worst is None or margin < worst:
            worst = margin; worst_set = (E, bound)
    print(f"{label}: {n_eval} cores, box [1..{W}]")
    print(f"  min margin (cap - bound) = {worst} ~ {float(worst):+.6f}  at {worst_set[0]}")
    print(f"    (bound there = {worst_set[1]} ~ {float(worst_set[1]):.6f}; cap = {capk1} ~ {float(capk1):.6f})")
    if consec_val:
        print(f"  consec: bound = {consec_val[0]} ~ {float(consec_val[0]):.6f}, margin = {consec_val[1]} ~ {float(consec_val[1]):+.6f}")
    print(f"  VERDICT: {'ROW HOLDS on the whole box' if worst > 0 else '*** ROW FAILS ***'}")
    return worst, worst_set

def main():
    print("(1) DEG-3 MAJORANT LP-VERTEX ENUMERATION (phi >= g on {0..6}, g = [1,1/3,0,0,0,0,0]):")
    for (c1, c2, c3), touch, vals in lp_vertices():
        print(f"    touch {touch}: c = ({c1}, {c2}, {c3})  phi = {[str(v) for v in vals]}")
    print()
    c1, c2, c3 = F(-2, 3), F(17, 90), F(-1, 45)
    print(f"(2) THE k=8 DEG-3 ROW, EXHAUSTIVE: Phi <= 1 - (2/3)m1 + (17/90)m2 - (1/45)m3 <= cap_9")
    W8 = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    census_row(8, W8, c1, c2, c3, CAP[9], f"k=8 deg-3 row")
    print()
    print(f"(3) THE k=9 DEG-2 ROW, EXHAUSTIVE (THM-705's linear requirement): Phi <= 1 - (2/3)m1 + (1/6)m2 <= cap_10")
    census_row(9, 17, F(-2, 3), F(1, 6), F(0), CAP[10], "k=9 deg-2 row")
    print()
    print("(4) CROSS-CHECK: the k=8 DEG-2 row fails (the known gap):")
    census_row(8, 14, F(-2, 3), F(1, 6), F(0), CAP[9], "k=8 deg-2 row (expected to fail)")

if __name__ == "__main__":
    main()
