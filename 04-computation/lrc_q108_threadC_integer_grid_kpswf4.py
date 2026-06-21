#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadC_integer_grid_kpswf4.py  (kind-pasteur 2026-06-21, THREAD C)

THREAD C: prove the SHARP integer-grid discrepancy constant for the (q,p) torus geodesic
cell occupancy at apex prime P=7.

The (q,p) torus geodesic v -> (qv, pv) mod 1, with the 7x7 sector grid, has an EXACT
integer description on the common 7pq-grid: subdivide [0,1) into 7pq equal subintervals.
Each subinterval lies entirely in one (i,j) sector cell. Let c_ij = number of subintervals
in cell (i,j). Then mu(i,j) = c_ij / (7pq) and

    D_{p,q} = sum_ij |mu(i,j) - 1/49| = (1/(7pq)) sum_ij |c_ij - pq/7|.

GOAL: characterize the integer matrix c = (c_ij), prove
    - row sums = col sums = pq  (doubly balanced)   [Lemma 1, marginals]
    - each row is a CYCLIC SHIFT of row 0 by slope s = p q^{-1} mod 7  [shift structure]
    - hence D = (1/(7pq)) * 7 * sum_j |r_j - pq/7|  where r = row 0
    - then a SHARP 1D bound on sum_j |r_j - pq/7|.

This file establishes the EXACT integer model and the shift structure, and computes the
1D row-0 vector r_j exactly, to feed the sharp combinatorial bound.
"""
from fractions import Fraction as Fr
from math import gcd
import sys

P = 7

def cmatrix(p, q):
    """Exact integer 7x7 cell-count matrix on the 7pq grid.

    Subinterval k = [k/(7pq), (k+1)/(7pq)), k=0..7pq-1.
    Its midpoint v=(k+0.5)/(7pq); cell = (floor(7 frac(qv)), floor(7 frac(pv))).
    Using midpoints avoids boundary ambiguity (boundaries are exactly grid lines).
    """
    N = 7 * p * q
    c = [[0] * P for _ in range(P)]
    for k in range(N):
        # v = (2k+1)/(2N); compute sectors exactly with integers
        # frac(q v): q*(2k+1)/(2N) ; sector = floor(7 * frac)
        num_q = q * (2 * k + 1)          # /(2N)
        # 7*frac(qv) = 7 * (num_q mod 2N)/(2N); floor =>
        i = (P * (num_q % (2 * N))) // (2 * N)
        num_p = p * (2 * k + 1)
        j = (P * (num_p % (2 * N))) // (2 * N)
        c[i][j] += 1
    return c

def D_from_c(p, q, c):
    N = 7 * p * q
    # D = (1/N) sum |c_ij - pq/7|, but pq/7 may be non-integer; use exact Fraction
    target = Fr(p * q, P)
    return sum(abs(Fr(c[i][j]) - target) for i in range(P) for j in range(P)) / N

def mu_full(p, q):
    """Reference exact occupancy via breakpoint subdivision (independent check)."""
    bp = {Fr(0), Fr(1)}
    for f in (p, q):
        for t in range(0, P * f):
            bp.add(Fr(t, P * f))
    vs = sorted(bp)
    cell = {}
    for a, b in zip(vs, vs[1:]):
        mid = (a + b) / 2
        key = (int(P * ((q * mid) % 1)), int(P * ((p * mid) % 1)))
        cell[key] = cell.get(key, Fr(0)) + (b - a)
    return cell

def main():
    print("THREAD C: exact integer-grid model c_ij on the 7pq grid")
    print("=" * 70)

    # ---- 1. Verify c_ij/(7pq) == mu(i,j), and D matches ----------------------
    ok_c = True
    ok_marg = True
    ok_shift = True
    n_checked = 0
    worst = []
    for q in range(1, 26):
        for p in range(q + 1, int(Fr(43, 20) * q) + 1):
            if gcd(p, q) != 1:
                continue
            if not (Fr(1) < Fr(p, q) <= Fr(43, 20)):
                continue
            c = cmatrix(p, q)
            N = 7 * p * q
            # cross-check against mu_full
            cell = mu_full(p, q)
            for i in range(P):
                for j in range(P):
                    if Fr(c[i][j], N) != cell.get((i, j), Fr(0)):
                        ok_c = False
            # marginals = pq
            for i in range(P):
                if sum(c[i]) != p * q:
                    ok_marg = False
            for j in range(P):
                if sum(c[i][j] for i in range(P)) != p * q:
                    ok_marg = False
            # shift structure: row i = row 0 shifted by s*i, s = p q^{-1} mod 7
            if q % P != 0:
                qinv = pow(q % P, -1, P)
                s = (p * qinv) % P
                for i in range(P):
                    for j in range(P):
                        if c[i][j] != c[0][(j - s * i) % P]:
                            ok_shift = False
            D = D_from_c(p, q, c)
            worst.append((D * q, D, p, q))
            n_checked += 1
    print(f"checked {n_checked} window ratios (q<=25)")
    print(f"  c_ij/(7pq) == mu(i,j) exact:        {'YES' if ok_c else 'NO'}")
    print(f"  row sums = col sums = pq:           {'YES' if ok_marg else 'NO'}")
    print(f"  row_i = cyclic shift of row_0 by s*i: {'YES' if ok_shift else 'NO'}  (s=p q^{{-1}} mod 7)")

    worst.sort(reverse=True)
    print("\nTop D*q (sharp constant target sup D*q = 12/7):")
    for dq, D, p, q in worst[:8]:
        print(f"  p/q={p}/{q:<3d} D*q={dq} ({float(dq):.6f})  D={D} ({float(D):.6f})")

    # ---- 2. Row-0 integer vector r_j and its 1D discrepancy ------------------
    print("\n" + "=" * 70)
    print("ROW-0 integer vector r_j (c[0][:]) and 1D L1 discrepancy")
    print("Goal: D = (7/(7pq)) * sum_j |r_j - pq/7| = (1/(pq)) sum_j |r_j - pq/7|")
    print("=" * 70)
    ok_d1d = True
    worst1d = []
    for q in range(1, 26):
        for p in range(q + 1, int(Fr(43, 20) * q) + 1):
            if gcd(p, q) != 1 or not (Fr(1) < Fr(p, q) <= Fr(43, 20)):
                continue
            c = cmatrix(p, q)
            r = c[0]
            target = Fr(p * q, P)
            D1 = sum(abs(Fr(rj) - target) for rj in r)  # sum_j |r_j - pq/7|
            D = D_from_c(p, q, c)
            # claim D == (1/(pq)) * D1  (since 7 rows each contribute D1, /(7pq))
            if D != D1 / (p * q):
                ok_d1d = False
            worst1d.append((Fr(D1) / q, D1, p, q, r, sum(r)))
    print(f"  D == (1/pq) * sum_j|r_j - pq/7|:  {'YES' if ok_d1d else 'NO'}")
    worst1d.sort(reverse=True, key=lambda t: t[0])
    print("\nTop (sum_j|r_j-pq/7|)/q  [== D*q]:")
    for d1q, D1, p, q, r, rs in worst1d[:8]:
        print(f"  p/q={p}/{q:<3d} rowsum={rs}=pq={p*q}  r={r}  sum|r-pq/7|={D1}  /q={d1q} ({float(d1q):.5f})")

    print("\nDONE.")

if __name__ == "__main__":
    main()
