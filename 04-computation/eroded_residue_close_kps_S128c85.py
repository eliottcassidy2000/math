#!/usr/bin/env python3
"""eroded_residue_close_kps_S128c85.py -- kind-pasteur-2026-07-19-S128c85

eroded_start_complex_kps_S128c85.py reduced supplier (3) to a FINITE residue:
the measure bound |E_w(P)| > 2/21 holds for every core once k1 > K* = 203.87,
while the clustered regime only guarantees k1 > 13*max(P).  This script closes
the residue by DIRECT EXACT COUNTING instead of measure.

For a core P and a killer k1, the k1-gaps sit at FIXED positions
    gap_j = [ j/k1 + 1/(14 k1),  (j+1)/k1 - 1/(14 k1) ],   j = 0,...,k1-1,
each of width w = 6/(7 k1).  Define the exact count

    G(P,k1) = #{ j : gap_j is contained in S(P) }.

This is the discrete, finite-k1 truth that |E_w(P)| only approximates: no
equidistribution, no error term.  By codex's PROVED ceiling (THM-1203) the bad
phases occupy measure <= 2/21, so at most about k1*(2/21) of the k1 gaps can be
bad.  Supplier (3) therefore holds at (P,k1) as soon as

    G(P,k1)  >  floor(k1 * 2/21)                                 (*)

(with the gap between "measure <= 2/21" and "count <= floor(k1*2/21)" being
exactly codex's supplier (2), flagged and NOT claimed here).

WHAT IS REPORTED
 1. The PER-CORE residue window ( 13*max(P), K*(P) ] -- narrower than the global
    (156, 203.87], because a core with a small maximum enters the clustered
    regime earlier.
 2. G(P,k1) over every residue pair (P,k1), and whether (*) holds throughout.
 3. The tightest pair, i.e. the smallest surplus G - floor(k1*2/21).
 4. A sanity check that G(P,k1)/k1 -> |E_w(P)| as k1 grows, which validates the
    measure computation against the exact count.
"""
import sys
from fractions import Fraction as F
from itertools import combinations

BAD = F(2, 21)
LAM = F(1, 14)


def safe_components(P):
    bps = {F(0), F(1)}
    for p in P:
        for j in range(p + 1):
            for s in (F(1, 14 * p), -F(1, 14 * p)):
                v = F(j, p) + s
                if 0 <= v <= 1:
                    bps.add(v)
    B = sorted(bps)
    out = []
    for i in range(len(B) - 1):
        a, b = B[i], B[i + 1]
        if b <= a:
            continue
        mid = (a + b) / 2
        if all(min((p * mid) % 1, 1 - (p * mid) % 1) >= LAM for p in P):
            out.append((a, b))
    merged = []
    for a, b in out:
        if merged and a == merged[-1][1]:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return merged


def critical_width(lengths, target=BAD):
    Ls = sorted(lengths, reverse=True)
    if sum(Ls, F(0)) <= target:
        return F(0)
    pref = F(0)
    best = F(0)
    for m in range(1, len(Ls) + 1):
        pref += Ls[m - 1]
        if pref <= target:
            continue
        w = (pref - target) / m
        if w < Ls[m - 1] and (True if m == len(Ls) else w >= Ls[m]):
            if w > best:
                best = w
    return best


def G(comps, k1):
    """Exact count of complete k1-gaps contained in S(P)."""
    n = 0
    for j in range(k1):
        lo = F(j, k1) + F(1, 14 * k1)
        hi = F(j + 1, k1) - F(1, 14 * k1)
        for a, b in comps:
            if a <= lo and hi <= b:
                n += 1
                break
    return n


CORES = [tuple(c) for c in combinations(range(1, 13), 8)]
print("=" * 78)
print("PER-CORE RESIDUE WINDOWS  ( 13*max(P) , K*(P) ]")
print("=" * 78)
info = {}
tot_pairs = 0
worst_window = None
for P in CORES:
    comps = safe_components(P)
    Ls = [b - a for a, b in comps]
    w = critical_width(Ls)
    K = F(6, 7) / w if w else None
    floor_k = 13 * max(P)
    hi = int(K) if K else 10 ** 9
    npairs = max(0, hi - floor_k)
    tot_pairs += npairs
    info[P] = (comps, Ls, w, K, floor_k, hi, npairs)
    if worst_window is None or npairs > worst_window[0]:
        worst_window = (npairs, P, floor_k, hi)
print("  cores with a NON-EMPTY residue window : %d of %d"
      % (sum(1 for P in CORES if info[P][6] > 0), len(CORES)))
print("  total (core, k1) residue pairs        : %d" % tot_pairs)
print("  widest single window                  : %d values, at P = %s, k1 in (%d, %d]"
      % (worst_window[0], worst_window[1], worst_window[2], worst_window[3]))

print()
print("=" * 78)
print("EXACT COUNTING OVER THE RESIDUE:  is  G(P,k1) > floor(k1 * 2/21) ?")
print("=" * 78)
viol = []
tightest = None
checked = 0
for P in CORES:
    comps, Ls, w, K, floor_k, hi, npairs = info[P]
    if npairs <= 0:
        continue
    for k1 in range(floor_k + 1, hi + 1):
        g = G(comps, k1)
        need = (k1 * 2) // 21
        checked += 1
        surplus = g - need
        if surplus <= 0:
            viol.append((P, k1, g, need))
        if tightest is None or surplus < tightest[0]:
            tightest = (surplus, P, k1, g, need)
print("  (core, k1) pairs checked exactly : %d" % checked)
print("  pairs FAILING  G > floor(k1*2/21) : %d  -> %s"
      % (len(viol), "RESIDUE CLOSES" if not viol else "RESIDUE DOES NOT CLOSE"))
if viol:
    print("  first 10 failures (P, k1, G, needed):")
    for row in viol[:10]:
        print("     %s" % (row,))
if tightest:
    s, P, k1, g, need = tightest
    print("  TIGHTEST pair : P = %s, k1 = %d" % (P, k1))
    print("     G(P,k1) = %d   floor(k1*2/21) = %d   surplus = %d" % (g, need, s))
    print("     |S(P)| = %.6f   N = %d   K*(P) = %.2f"
          % (float(sum(info[P][1], F(0))), len(info[P][1]), float(info[P][3])))

print()
print("=" * 78)
print("VALIDATION:  G(P,k1)/k1  ->  |E_w(P)|  as k1 grows (measure vs exact count)")
print("=" * 78)
Pv = (1, 2, 3, 5, 7, 8, 9, 11)
comps = safe_components(Pv)
Ls = [b - a for a, b in comps]
print("  bottleneck core P* = %s" % (Pv,))
for k1 in (160, 200, 300, 600, 1200, 2400):
    w = F(6, 7 * k1)
    Em = sum((L - w for L in Ls if L > w), F(0))
    g = G(comps, k1)
    print("     k1=%-5d  G=%-5d  G/k1 = %.6f   |E_w| = %.6f   ratio %.4f   (2/21 = %.6f)"
          % (k1, g, g / k1, float(Em), (g / k1) / float(Em) if Em else 0, 2 / 21))
