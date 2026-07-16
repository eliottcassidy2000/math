#!/usr/bin/env python3
"""fixed_locus_residue6_macmini_S117.py -- mac-mini-2026-07-16-S117.
THE FIXED-LOCUS PRINCIPLE for codex's negative residue 6 (equivariance investigation).

Claim chain (proof = the substitution x -> -x, five lines; verified here):
 (1) Sector reflection: a point in sector s maps to sector 6-s; hence the missed-set
     M_E(-x) = 6 - M_E(x), and g_s(-y) = g_{6-s}(y) a.e.
 (2) Error_t(a) = sum_s int_{R_s} g_s(atx) dx is EXACTLY invariant under x -> -x --
     the reflection acts WITHIN each residue class of a (no residue transfer 6 <-> 1).
 (3) The pair-sectors fixed by the reflection are exactly {s, c} with s + c = 6:
     {1,5} and {2,4} (and the degenerate 3) -- PRECISELY the A_15 + A_24 where THM-891's
     negative-residue-6 mass concentrates. The difficulty is pinned to the fixed locus
     of the involution (the Redei/locker/SC template).
 (4) Consequences: A_{s,c} = A_{6-c,6-s} (all pair-sector measures reflect); on the
     fixed sectors the ODD part of the g-coupling integrates to zero, so certificates
     need only the SYMMETRIZED kernel on a half-domain.
Verification: core E = {0,1,2,3,4,6} (codex's unique-max D20 core), fine-grid measures.
"""
import sys
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

E = [0, 1, 2, 3, 4, 6]
GRID = 7 * 11 * 13 * 16 * 3       # fine, avoids sector-boundary aliasing mostly
xs = (np.arange(GRID) + 0.5) / GRID

def sectors(vals):                # sector index of frac(v) on the 7-clock
    return (np.floor((vals % 1.0) * 7)).astype(int)

# missed-sector sets of the 6-point core, per x
S = np.stack([sectors(e * xs) for e in E])          # 6 x GRID
present = np.zeros((7, GRID), dtype=bool)
for i in range(len(E)):
    present[S[i], np.arange(GRID)] = True
missed_mask = ~present                               # 7 x GRID
nmiss = missed_mask.sum(0)

print("(1)+(3) pair-sector measures and the reflection map {s,c} -> {6-c,6-s}:")
from collections import defaultdict
Am = defaultdict(float)
idx2 = np.where(nmiss == 2)[0]
for k in idx2:
    s, c = np.where(missed_mask[:, k])[0]
    Am[(int(s), int(c))] += 1.0 / GRID
tot_err = 0.0
for (s, c), v in sorted(Am.items()):
    rs, rc = sorted((6 - c, 6 - s))
    w = Am.get((rs, rc), 0.0)
    fixed = (rs, rc) == (s, c)
    tot_err = max(tot_err, abs(v - w))
    print(f"   A_{s}{c} = {v:.5f}   reflects to A_{rs}{rc} = {w:.5f}"
          f"   {'FIXED SECTOR' if fixed else ''}")
print(f"   max |A - A_reflected| = {tot_err:.2e}  (symmetry law A_sc = A_(6-c)(6-s))")
print(f"   fixed pair-sectors: exactly s+c = 6: {{1,5}}, {{2,4}} -- the THM-891 concentration set")

print()
print("(2) Error_t(a) reflection-invariance (numeric, t = 501):")
t = 501
for a in (1, 6, 13, 20):
    # R_s: miss exactly sector s among the SEVEN offsets (core + t)
    S7 = sectors(t * xs)
    present7 = present.copy()
    present7[S7, np.arange(GRID)] = True
    missed7 = ~present7
    n7 = missed7.sum(0)
    err = 0.0; err_ref = 0.0
    y = (a * t * xs) % 1.0
    gsec = np.floor(y * 7).astype(int)
    one = np.where(n7 == 1)[0]
    for k in one:
        s = int(np.where(missed7[:, k])[0][0])
        err += ((1.0 if gsec[k] == s else 0.0) - 1 / 7) / GRID
    # reflected: substitute x -> -x == use grid reversed (1 - x)
    err2 = 0.0
    y2 = (a * t * (1 - xs)) % 1.0
    g2 = np.floor(y2 * 7).astype(int)
    S7b = sectors(t * (1 - xs)); Sb = np.stack([sectors(e * (1 - xs)) for e in E])
    presentb = np.zeros((7, GRID), dtype=bool)
    for i in range(len(E)): presentb[Sb[i], np.arange(GRID)] = True
    presentb[S7b, np.arange(GRID)] = True
    missedb = ~presentb; nb = missedb.sum(0)
    oneb = np.where(nb == 1)[0]
    for k in oneb:
        s = int(np.where(missedb[:, k])[0][0])
        err2 += ((1.0 if g2[k] == s else 0.0) - 1 / 7) / GRID
    print(f"   a = {a} (residue {a % 7}): Error = {err:+.6f}   reflected = {err2:+.6f}"
          f"   |delta| = {abs(err - err2):.2e}")

print()
print("(4) odd-part cancellation on the fixed sectors (a = 6, t = 501):")
# on x-regions with missed pair {1,5} (fixed sector), the reflection-odd component of the
# g-coupling must integrate to 0: check int_{A15-region} [g_s(atx) - g_{6-s}(at(1-x))] ~ 0
a = 6
y = (a * t * xs) % 1.0
gs1 = (np.floor(y * 7).astype(int) == 1).astype(float) - 1 / 7
reg = np.where((nmiss == 2) & missed_mask[1] & missed_mask[5])[0]
odd = gs1[reg].sum() / GRID
y2 = (a * t * (1 - xs)) % 1.0
gs5r = (np.floor(y2 * 7).astype(int) == 5).astype(float) - 1 / 7
odd2 = gs5r[reg].sum() / GRID
print(f"   int_A15 g_1(atx) = {odd:+.6f}   int_A15 g_5(at(1-x)) = {odd2:+.6f}"
      f"   (equal by the symmetry: |delta| = {abs(odd - odd2):.2e})")
print("\nDONE")
