#!/usr/bin/env python3
"""THM-786 exact replay (opus-S304, scope-corrected by codex-S10).

This file exactly checks the two extreme-ratio wall-count certificates and
their extents.  It does *not* implement the session-reported generic /
extreme-ratio / balanced-pair / near-multiple census or annealing run; the
legacy 0.589 line is printed only as an explicitly unreplayed report.

Tournament Analysis is intentionally not recomputed here: the wall chronology
tournament is transitive (one Hamiltonian path) and forgets both metric extent
and which companion serves which g-wall.  The theorem-facing carrier is the
labelled incidence relation (g-wall, f-period, companion); the exact tournament
guardrail is audited in lrc14_r8_raw_wall_refuter_codex_S10.py.
"""
import random
from fractions import Fraction as F
def inv7(w): return pow(w % 7, 5, 7)
def walls_window(W8, lo, hi):
    ev = []
    for w in W8:
        for tm in range(int(2*w*lo), int(2*w*hi)+2):
            if tm % 2:
                x = F(tm, 2*w)
                if lo < x < hi: ev.append((x, w, (tm-1)//2))
    ev.sort(); return ev
def token_at(v, x):
    q = v*x; n,d = q.numerator, q.denominator
    if d == 2 and n % 2 == 1: return None
    return (-inv7(v) * ((2*n+d)//(2*d))) % 7
def wall_ok(W8, x, w):
    toks = []
    for v in W8:
        if v == w: continue
        t = token_at(v, x)
        if t is None: return False
        toks.append(t)
    return len(set(toks)) == 7
def run_stats(W8, lo, hi):
    ev = walls_window(W8, lo, hi)
    bk, bext, cur = 0, F(0), []
    for x, w, m in ev:
        if wall_ok(W8, x, w): cur.append(x)
        else:
            bk = max(bk, len(cur))
            if cur: bext = max(bext, cur[-1]-cur[0])
            cur = []
    bk = max(bk, len(cur))
    if cur: bext = max(bext, cur[-1]-cur[0])
    return bk, bext
if __name__ == '__main__':
    # (1) the artifact certificates
    for W8 in ([10,12,17,18,22,32,39,2445], [8,10,18,24,32,34,39,3887]):
        ws = sorted(W8); bound = F(1, ws[-2]) + F(2, ws[-1])
        k, ext = run_stats(W8, F(1,100), F(1,100)+F(1,4))
        print(f"certificate {W8}: run = {k} WALLS (K0=6 refuted), "
              f"extent = {float(ext):.5f} < bound {float(bound):.5f}: {ext < bound}")
    # (2) Legacy session report; no census generator was stored in this file.
    random.seed(3041)
    print("extent census NOT REPLAYED: S304 reported peak ratio 0.589 (annealed).")
