#!/usr/bin/env python3
"""
klein-2026-07-01-S88 -- HYP-3842: the ARC x RADIUS TWO-VARIABLE LP

Adversary relaxation: fractional speed profile x_v in [0,1], v <= W, sum x_v = 13,
covering pins sum_{q | v} x_v >= 1 for q = 2..14. Speed v's danger geometry at radius r
is DETERMINISTIC: m[v][c][j] = |danger(v, r_j) cap arc_c| (exact, precomputed).
Coverage variables: cov_{c,j} <= arc length, cov_{c,j} <= sum_v x_v m[v][c][j].
Adversary MAXIMIZES total coverage at the target row; floor = 1 - max.

Three instruments compared:
  L1: single-radius (target row only)                        -- baseline (HYP-3824-style)
  L2: arc x radius (all rows share x; target objective)      -- the two-variable coupling
  L3: L2 + pairwise-overlap cuts from the exact small-speed overlap law
      |D_p cap D_q| at r: for p+q <= 1/r-band the overlap is forced (kps S28 Farey law
      F(p,q)=1/(7 max) iff p+q<=14 at r=1/14); McCormick-linearized via y_{pq} >= x_p+x_q-1.
      Overlap forces WASTED mass: effective coverage <= sum x m - sum_pq y_pq OV(p,q,j).

Honest goal: measure whether radius-coupling and overlap cuts move the relaxation value
off zero at sub-critical target rows. A zero is an informative negative (the relaxation
still smears the atom); a positive value = a UNIVERSAL floor for all covering sets with
speeds <= W (modulo the fractional relaxation, which only HELPS the adversary).
"""
from fractions import Fraction as F
import itertools
import numpy as np
from scipy.optimize import linprog

W = 200                      # speed cap for the relaxation
Q = 28                       # arcs
ROWS = [F(1, 16), F(1, 15), F(1, 14)]
TARGETS = [0, 1, 2]          # certify each row in turn

def danger_arc_mass(v, r, Q):
    """Exact |danger(v,r) cap arc_c| for all c; returns list of floats (exact fractions -> float)."""
    out = [F(0)] * Q
    rv = r / v
    for a in range(v + 1):
        lo, hi = F(a, v) - rv, F(a, v) + rv
        lo, hi = max(lo, F(0)), min(hi, F(1))
        if hi <= lo:
            continue
        c0, c1 = int(lo * Q), min(int(hi * Q), Q - 1)
        for c in range(c0, c1 + 1):
            alo, ahi = F(c, Q), F(c + 1, Q)
            seg = min(hi, ahi) - max(lo, alo)
            if seg > 0:
                out[c] += seg
    return [float(x) for x in out]

def pair_overlap(p, q, r):
    """Exact |D_p cap D_q| at radius r (p != q)."""
    ivs = []
    for v in (p, q):
        rv = r / v
        for a in range(v + 1):
            lo, hi = F(a, v) - rv, F(a, v) + rv
            lo, hi = max(lo, F(0)), min(hi, F(1))
            if hi > lo:
                ivs.append((v, lo, hi))
    tot = F(0)
    P = [(lo, hi) for v, lo, hi in ivs if v == p]
    Qi = [(lo, hi) for v, lo, hi in ivs if v == q]
    for lo1, hi1 in P:
        for lo2, hi2 in Qi:
            seg = min(hi1, hi2) - max(lo1, lo2)
            if seg > 0:
                tot += seg
    return float(tot)

print("precomputing danger-arc masses (W=%d, Q=%d, %d rows)..." % (W, Q, len(ROWS)))
M = {}  # M[j][v] = list over arcs
for j, r in enumerate(ROWS):
    M[j] = {v: danger_arc_mass(v, r, Q) for v in range(1, W + 1)}

# overlap cuts: small-speed pairs (both <= 14, p+q <= 20) -- the structured, forced overlaps
PAIRS = [(p, q) for p, q in itertools.combinations(range(1, 15), 2)]
OV = {}
for j, r in enumerate(ROWS):
    OV[j] = {(p, q): pair_overlap(p, q, r) for p, q in PAIRS}

def build_lp(target, rows_used, use_overlap):
    """Vars: x_v (W), cov_{c,j} for j in rows_used (Q per row), y_pq (if overlap).
    Maximize sum_c cov_{c,target}  ==  linprog minimizes -that."""
    nx = W
    row_ix = {j: i for i, j in enumerate(rows_used)}
    ncov = Q * len(rows_used)
    pairs = PAIRS if use_overlap else []
    ny = len(pairs)
    N = nx + ncov + ny
    def xi(v): return v - 1
    def ci(c, j): return nx + row_ix[j] * Q + c
    def yi(k): return nx + ncov + k

    c_obj = np.zeros(N)
    for c in range(Q):
        c_obj[ci(c, target)] = -1.0

    A_ub, b_ub = [], []
    # cov_{c,j} <= sum_v x_v m[v][c][j]  (+ overlap waste subtracted at row j)
    for j in rows_used:
        for c in range(Q):
            row = np.zeros(N)
            row[ci(c, j)] = 1.0
            for v in range(1, W + 1):
                row[xi(v)] -= M[j][v][c]
            A_ub.append(row); b_ub.append(0.0)
    # overlap waste: total coverage at row j <= sum_v x_v * 2r - sum_pairs y_pq OV
    if use_overlap:
        for j in rows_used:
            row = np.zeros(N)
            for c in range(Q):
                row[ci(c, j)] = 1.0
            for v in range(1, W + 1):
                row[xi(v)] -= 2 * float(ROWS[j])
            for k, (p, q) in enumerate(pairs):
                row[yi(k)] += OV[j][(p, q)]
            A_ub.append(row); b_ub.append(0.0)
        # McCormick: y_pq >= x_p + x_q - 1  ->  x_p + x_q - y_pq <= 1
        for k, (p, q) in enumerate(pairs):
            row = np.zeros(N)
            row[xi(p)] += 1.0; row[xi(q)] += 1.0; row[yi(k)] -= 1.0
            A_ub.append(row); b_ub.append(1.0)

    A_eq, b_eq = [], []
    row = np.zeros(N)
    for v in range(1, W + 1):
        row[xi(v)] = 1.0
    A_eq.append(row); b_eq.append(13.0)
    # covering pins
    for qq in range(2, 15):
        row = np.zeros(N)
        for v in range(qq, W + 1, qq):
            row[xi(v)] = -1.0
        A_ub.append(row); b_ub.append(-1.0)

    bounds = [(0, 1)] * nx + [(0, 1.0 / Q)] * ncov + [(0, 1)] * ny
    res = linprog(c_obj, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds,
                  method="highs")
    return res

print("\n%-46s %14s %14s" % ("instrument", "max coverage", "floor 1-max"))
for target in TARGETS:
    rt = ROWS[target]
    for label, rows_used, ov in [
        ("L1 single-radius", [target], False),
        ("L2 arc x radius (3 rows)", [0, 1, 2], False),
        ("L3 + overlap cuts (small-speed pairs)", [0, 1, 2], True),
    ]:
        res = build_lp(target, rows_used, ov)
        if res.status == 0:
            mx = -res.fun
            print("%-46s %14.6f %14.6f   [target r=%s]" % (label, mx, 1 - mx, rt))
        else:
            print("%-46s FAILED %s" % (label, res.message))
    print()

print("interpretation: floor > 0 at a sub-critical target = a universal covering floor")
print("(for all covering sets, speeds <= %d, modulo fractional relaxation)." % W)
print("floor = 0 = the relaxation smears the atom; identifies the missing constraint layer.")
print("DONE.")
