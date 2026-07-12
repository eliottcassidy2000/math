# -*- coding: utf-8 -*-
# klein-2026-07-11-S253: EXHAUSTIVE SERVICE on the two newest extremality objects.
#
# (A) mac-mini cont.36 (THM-712 cubic base form): the k=8 FREE-CORE deg-3 requirement
#     bound(E) = min over feasible deg-3 vertices of  a + b m1 + c m2 + d m3  <=  cap_9,
#     verified there ADVERSARIALLY (50 families, worst {1..8} at 0.4483).  HERE: EXHAUSTIVE
#     over ALL free 8-cores E in [0..W8] up to dilation (all-integer breakpoint sweep,
#     factorial moments over the FULL seven-sector empty count, support {0..7}).
# (B) kps cont.33 (HYP-6015): consec MAXIMIZES Var(N) -- their exhaustive box was [1..14].
#     HERE: wider boxes for k = 8, 9 (free cores, Var = m2 + m1 - m1^2).
#
# Conventions == mac-mini lrc14_cubic_base_macmini_S65cont36.py exactly:
# sectors [c/7,(c+1)/7); N = 7 - #occupied (ALL seven sectors; e = 0 occupies sector 0);
# g = [1, 1/3, 0, 0, 0, 0, 0, 0] on {0..7}; deg-3 vertices = 4 tight constraints, feasible.

import sys
from math import gcd
from fractions import Fraction as F
from itertools import combinations, permutations

CAP9 = F(1979, 4004)
G = [F(1), F(1, 3)] + [F(0)] * 6

def det4(M):
    s = F(0)
    for perm in permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        prod = F(1)
        for r in range(4):
            prod *= M[r][perm[r]]
        s += (-1) ** inv * prod
    return s

def basis(N):
    return [F(1), F(N), F(N * (N - 1)), F(N * (N - 1) * (N - 2))]

def vertices():
    VS = []
    for quad in combinations(range(8), 4):
        M = [basis(N) for N in quad]
        R = [G[N] for N in quad]
        D = det4(M)
        if D == 0:
            continue
        sol = []
        for c in range(4):
            Mc = [row[:] for row in M]
            for r in range(4):
                Mc[r][c] = R[r]
            sol.append(det4(Mc) / D)
        a, b, c3, d = sol
        if all(a + b * N + c3 * N * (N - 1) + d * N * (N - 1) * (N - 2) >= G[N]
               for N in range(8)):
            VS.append((a, b, c3, d, quad))
    return VS

def lcm(a, b):
    return a // gcd(a, b) * b

def moments_free(E):
    """(m0..m3, p) for FREE core E (0 allowed): N = 7 - #occupied among ALL 7 sectors."""
    nz = [e for e in E if e != 0]
    has0 = len(nz) < len(E)
    L = 1
    for e in nz:
        L = lcm(L, e)
    D = 7 * L
    pts = set([0, D])
    for e in nz:
        step = L // e
        pts.update(range(0, D + 1, step))
    pts = sorted(pts)
    pn = [0] * 8
    base_hit = 1 if has0 else 0
    for t1, t2 in zip(pts, pts[1:]):
        s = t1 + t2
        hit = base_hit
        for e in nz:
            c = (7 * e * s // (2 * D)) % 7
            hit |= 1 << c
        n = 7 - bin(hit).count("1")
        pn[n] += t2 - t1
    p = [F(x, D) for x in pn]
    m1 = sum(F(j) * p[j] for j in range(8))
    m2 = sum(F(j * (j - 1)) * p[j] for j in range(8))
    m3 = sum(F(j * (j - 1) * (j - 2)) * p[j] for j in range(8))
    return m1, m2, m3

def norm_dilation(E):
    """dilation-normalize: divide by gcd of nonzero elements (0 kept)."""
    nz = [e for e in E if e != 0]
    g0 = 0
    for e in nz:
        g0 = gcd(g0, e)
    if g0 <= 1:
        return tuple(E)
    return tuple(e // g0 for e in E)

def service_A(W8):
    VS = vertices()
    print(f"(A) k=8 FREE-CORE deg-3 requirement, exhaustive [0..{W8}] "
          f"({len(VS)} feasible vertices):")
    seen = set()
    worst = None; worst_set = None
    n_eval = 0
    for combo in combinations(range(0, W8 + 1), 8):
        key = norm_dilation(combo)
        if key in seen:
            continue
        seen.add(key)
        n_eval += 1
        m1, m2, m3 = moments_free(key)
        bnd = min(a + b * m1 + c3 * m2 + d * m3 for a, b, c3, d, _ in VS)
        margin = CAP9 - bnd
        if worst is None or margin < worst:
            worst = margin; worst_set = (key, bnd)
    print(f"  {n_eval} dilation-normalized free cores")
    print(f"  min margin = {worst} ~ {float(worst):+.6f} at {worst_set[0]}")
    print(f"  (bound = {worst_set[1]} ~ {float(worst_set[1]):.6f}; cap_9 ~ {float(CAP9):.6f})")
    print(f"  VERDICT: {'REQUIREMENT HOLDS on the whole box' if worst > 0 else '*** FAILS ***'}")

def service_B(k, W):
    print(f"(B) Var(N)-max over free {k}-cores in [1..{W}] (kps cont.33 extension):")
    best = None; best_set = None
    consec = tuple(range(1, k + 1))
    consec_var = None
    n_eval = 0
    seen = set()
    for combo in combinations(range(1, W + 1), k):
        key = norm_dilation(combo)
        if key in seen:
            continue
        seen.add(key)
        n_eval += 1
        m1, m2, m3 = moments_free(key)
        var = m2 + m1 - m1 * m1
        if key == consec:
            consec_var = var
        if best is None or var > best:
            best = var; best_set = key
    is_consec = (best_set == consec)
    print(f"  {n_eval} normalized cores; argmax Var = {best_set} "
          f"{'== CONSEC' if is_consec else '*** NON-CONSEC ***'}")
    print(f"  max Var = {best} ~ {float(best):.6f}"
          + (f"; consec Var = {consec_var} ~ {float(consec_var):.6f}" if consec_var is not None and not is_consec else ""))

if __name__ == "__main__":
    W8 = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    service_A(W8)
    print()
    service_B(9, 16)
    print()
    service_B(8, 17)
