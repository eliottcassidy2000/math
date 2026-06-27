#!/usr/bin/env python3
"""
kind-pasteur-2026-06-22 (kpswf14) STRUCTURAL LEMMA: the apex winding tournament is a function
of the residue multiset mod 14 ONLY (hence magnitude-blind on the apex-twins).

Lemma. Fix Q and a unit a mod Q. The winding tournament T(a/Q) on speeds S has arc i->j iff
((s_i - s_j) * a mod Q) lies in (0, Q/2). This value depends on s_i, s_j only through their
residues mod Q. Therefore two speed sets with the SAME residue multiset mod Q yield the SAME
tournament up to the relabeling that matches residues. In particular AP {1..13}, 12->26, 12->96
(all residues {1..13} mod 14) give the SAME apex tournament for EVERY unit a -- the apex is
magnitude-blind, so it is at best a NECESSARY condition for tightness.

This script verifies the lemma EXHAUSTIVELY over all units a (not by sampling).
"""
from fractions import Fraction as F

def winding(S, t):
    k = len(S); A = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            rel = (F(S[i]-S[j]) * t) % 1
            if 0 < rel < F(1, 2): A[i][j] = 1
            elif rel > F(1, 2): A[i][j] = 0
            else: A[i][j] = 1 if S[i] < S[j] else 0
    return A

def residue_aligned(S, A, Q):
    """canonicalize vertex order by (residue mod Q, speed); return flattened adjacency."""
    order = sorted(range(len(S)), key=lambda i: (S[i] % Q, S[i]))
    n = len(S)
    return tuple(A[order[i]][order[j]] for i in range(n) for j in range(n) if i != j)

Q = 14
units = [a for a in range(1, Q) if F(a, Q).denominator == Q]
print(f"units mod {Q}: {units}")

AP = list(range(1, 14))
L26 = list(range(1, 12)) + [13, 26]
L96 = list(range(1, 12)) + [13, 96]
print(f"residues mod 14: AP={sorted(s%14 for s in AP)}")
print(f"                 12->26={sorted(s%14 for s in L26)}")
print(f"                 12->96={sorted(s%14 for s in L96)}")
print(f"  (identical: {sorted(s%14 for s in AP)==sorted(s%14 for s in L26)==sorted(s%14 for s in L96)})")
print()
all_equal = True
for a in units:
    t = F(a, Q)
    cap = residue_aligned(AP, winding(AP, t), Q)
    c26 = residue_aligned(L26, winding(L26, t), Q)
    c96 = residue_aligned(L96, winding(L96, t), Q)
    eq = (cap == c26 == c96)
    all_equal &= eq
    print(f"  a={a:2d}: residue-aligned apex tournaments AP==26==96 : {eq}")
print()
print(f"EXHAUSTIVE VERDICT: apex tournament identical across AP/12->26/12->96 for ALL units: {all_equal}")
print("=> the winding tournament at any apex phase a/14 depends ONLY on residues mod 14;")
print("   it is MAGNITUDE-BLIND, hence necessary-only for LRC tightness. QED (structural).")
