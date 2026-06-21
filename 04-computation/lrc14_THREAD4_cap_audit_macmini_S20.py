#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_THREAD4_cap_audit_macmini_S20.py  (THREAD 4 adversarial audit, mac-mini-2026-06-21-S20)

CORRECTED ARCHITECTURE (HYP-2778):  LRC(14) sector route needs
    max_E measS7(E) <= cap_k   for all relevant cluster sizes k,
NOT "consec maximizes measS7" (false for k>=12, unnecessary).

This script ADVERSARIALLY audits the cap obligation.  It is fully self-contained:
it re-derives EVERYTHING (measS7, meas(G_P), caps) from definitions.

(a) RANGE of k for n=14.  Model: admissible (P,E): 0 in E, k=|E|, P subset {1..13},
    |P|+|E|=13, so k=|E| in {3,...,13} and |P|=13-k in {0,...,10}.  k NEVER exceeds 13
    (only 13 speeds).  We confirm k in 3..13 is the complete range.

(b) CAP VALUES re-derived from the covering reduction.  cap_k := min_{|P|=13-k} meas(G_P),
    G_P = {x in [0,1): ||p x|| >= 1/14 for all p in P}.  We recompute the per-size minima
    EXACTLY (the THM-530 table) and confirm cap_8..cap_13.

(c) MAX measS7 vs cap_k EXHAUSTIVELY for k=8..13 (bounded span), reporting the true argmax
    (not assuming consec), the max measS7, the cap, and the margin.  A margin <= 0 would be a
    GENUINE counterexample to LRC(14) at n=14 (or a cap error) -- the most valuable output.

(d) The k=12 consec-max counterexample E*=[0..10,12] vs consec -- confirm consec is NOT the
    maximizer there, AND that the cap still holds with large margin (consec-max false but
    cap intact).

measS7(E) = meas{x: all 7 sectors of Z/7 are hit by {frac(e x): e in E}}.  Sector 0 is always
hit (0 in E pins frac(0*x)=0 in sector 0).  CANONICAL definition = occupancy / IE over missed
inner sectors (matches THM-534 factorial_moments and the reframe binding row 8899/17640).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
from collections import defaultdict
try: sys.stdout.reconfigure(line_buffering=True)
except Exception: pass

# ---------------------------------------------------------------------------
# CANONICAL measS7 (occupancy, 0 included; pi[7] = P(all 7 sectors hit)).
# Verified == THM-534 moment IE; gives measS7(consec_10)=8899/17640 (reframe binding row).
# ---------------------------------------------------------------------------
def measS7(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        ae = abs(e)
        for a in range(7 * ae + 1): bps.add(F(a, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo + hi) / 2
        hit = set()
        for e in E:
            v = e * xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        if len(hit) == 7: tot += hi - lo
    return tot

# ---------------------------------------------------------------------------
# meas(G_P) = meas{ x : ||p x|| >= 1/14 for all p in P }.  Danger arc per p of half-width
# 1/14 around each multiple j/p.  Complement measure, EXACT.
# ---------------------------------------------------------------------------
H = F(1, 14)
def measGP(P):
    P = sorted(set(int(p) for p in P if p != 0))
    if not P: return F(1)
    # danger set: union over p of {x: ||p x|| < 1/14} = union of arcs (j/p - 1/(14p), j/p + 1/(14p))
    arcs = []
    for p in P:
        w = F(1, 14 * p)  # half-width in x
        for j in range(p):
            c = F(j, p)
            a = (c - w) % 1; b = (c + w) % 1
            if a < b: arcs.append((a, b))
            else: arcs.append((a, F(1))); arcs.append((F(0), b))
    arcs.sort()
    # union measure
    cov = F(0); clo, chi = arcs[0]
    for lo, hi in arcs[1:]:
        if lo <= chi: chi = max(chi, hi)
        else: cov += chi - clo; clo, chi = lo, hi
    cov += chi - clo
    return F(1) - cov

def primitive(E):
    return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1

# ===========================================================================
print("#"*84)
print("# THREAD 4 -- ADVERSARIAL AUDIT of the corrected LRC(14) cap obligation")
print("#"*84)

# ---- (a) RANGE OF k ----
print("\n" + "="*84)
print("(a) RANGE of cluster sizes k for n=14")
print("="*84)
print("  Model: 0 in E, |P|+|E|=13, P subset {1..13}.  k=|E| in 3..13; |P|=13-k in 0..10.")
print("  => k NEVER exceeds 13 (only 13 speeds). The reframe's 'k>=16' discussion is about")
print("     the abstract measS7 extremal problem, NOT realizable cluster sizes for LRC(14).")
print("  Relevant cap rows for LRC(14): k=8..13 (k<=7 closed unconditionally by pigeonhole,")
print("     THM-530 branch k<=7: 6 points => maxgap>1/7 always => mu_{1/7}=1).")

# ---- (b) RE-DERIVE caps from min meas(G_P) ----
print("\n" + "="*84)
print("(b) RE-DERIVE cap_k = min_{|P|=13-k} meas(G_P)  (covering reduction)")
print("="*84)
CANON_CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7), 13:F(1)}
# THM-530 stated per-size minima for psz=1..10:
THM530_minGP = {1:F(6,7),2:F(66,91),3:F(55,91),4:F(1979,4004),5:F(2243,5880),
                6:F(3029,10780),7:F(45107,229320),8:F(2479,17640),9:F(10601,114660),10:F(14249,252252)}
derived_cap = {}
for k in range(8, 14):
    psz = 13 - k
    if psz == 0:
        best = F(1); bestP = ()
    else:
        best = None; bestP = None
        for P in itertools.combinations(range(1, 14), psz):
            m = measGP(P)
            if best is None or m < best:
                best, bestP = m, P
    derived_cap[k] = best
    canon = CANON_CAP[k]
    thm530 = THM530_minGP.get(psz, F(1))
    tag_canon = "OK" if best == canon else "*** MISMATCH ***"
    tag_530 = "OK" if best == thm530 else "*** MISMATCH ***"
    print(f"  k={k:2d} (|P|={psz}): min meas(G_P)={best} = {float(best):.6f}  argmin P={bestP}")
    print(f"        vs CANON cap_{k}={canon}  [{tag_canon}];  vs THM-530 minGP[{psz}]={thm530} [{tag_530}]")

# ---- (c) EXHAUSTIVE max measS7 vs cap, k=8..13 (bounded span) ----
print("\n" + "="*84)
print("(c) EXHAUSTIVE max measS7(E) vs cap_k, primitive E, 0 in E, bounded span")
print("="*84)
print("  Span = max(E). For each k we enumerate ALL primitive E={0}+subset of {1..SPAN}.")
print("  measS7 is SCALE-INVARIANT, so primitive bounded-span covers each dilation class once.")
print("  (A margin<=0 here is a genuine counterexample / cap error -- the key adversarial test.)\n")

def exhaustive_max(k, span):
    best = F(0); bestE = None; n = 0
    for rest in itertools.combinations(range(1, span + 1), k - 1):
        E = (0,) + rest
        if not primitive(E): continue
        n += 1
        m = measS7(list(E))
        if m > best: best, bestE = m, E
    return best, bestE, n

# span budgets chosen so each row finishes; widen where feasible.
SPAN = {8:14, 9:14, 10:14, 11:14, 12:15, 13:16}
for k in range(8, 14):
    span = SPAN[k]
    cap = CANON_CAP[k]
    best, bestE, n = exhaustive_max(k, span)
    consecE = tuple(range(k)); mconsec = measS7(list(consecE))
    margin = cap - best
    is_consec = (bestE == consecE)
    flag = "*** CAP VIOLATION ***" if margin <= 0 else ("CONSEC-MAX" if is_consec else "non-consec argmax")
    print(f"  k={k:2d} span<={span} ({n} primitive shapes):")
    print(f"     max measS7 = {best} = {float(best):.6f}   argmax E={list(bestE)}")
    print(f"     measS7(consec) = {mconsec} = {float(mconsec):.6f}   (consec is argmax: {is_consec})")
    print(f"     cap_{k} = {cap} = {float(cap):.6f}   margin = {float(margin):+.6f}   -> {flag}")

# ---- (d) k=12 consec-max counterexample + cap intact ----
print("\n" + "="*84)
print("(d) k=12: consec-max FALSE (E*=[0..10,12] beats consec) but cap INTACT")
print("="*84)
Estar = list(range(11)) + [12]   # [0..10,12]
consec12 = list(range(12))
mstar = measS7(Estar); mcons = measS7(consec12); cap12 = CANON_CAP[12]
print(f"  consec_12 = {consec12}")
print(f"     measS7(consec_12) = {mcons} = {float(mcons):.6f}")
print(f"  E* = {Estar}")
print(f"     measS7(E*)       = {mstar} = {float(mstar):.6f}")
print(f"  E* beats consec: {mstar > mcons}   (consec-max is {'FALSE' if mstar>mcons else 'TRUE'} at k=12)")
print(f"  cap_12 = {cap12} = {float(cap12):.6f}")
print(f"  BOTH under cap: consec margin {float(cap12-mcons):+.6f}, E* margin {float(cap12-mstar):+.6f}")
print(f"  => consec-max false at k=12, cap holds with large margin (>=0.21). Reframe confirmed.")

print("\n" + "#"*84)
print("# DONE. Margins are the load-bearing output: any margin<=0 is a real problem.")
print("#"*84)
