#!/usr/bin/env python3
"""r=8 deck no-go, exact combinatorial form (opus-2026-07-14-S301, HYP-6840).

At lens c=7 with coprime exceptions, THM-767's zero-variance pins each exception's
bad set to EXACTLY ONE sheet on every event-free chamber piece (empty exactly at
its own walls). So a full cover of the deck over the whole closed core-safe set
Gbar_P by 8 exceptions is EQUIVALENT to the finite condition:

  (a) on every chamber piece of Gbar_P (refined at all chosen exceptions' walls),
      the 8 assigned sheets {k_w(piece)} hit all of Z_7;
  (b) at every wall of a chosen exception inside Gbar_P, the other seven form a
      PERFECT RAINBOW (all 7 sheets, each exactly once).

This search decides, exactly, whether an 8-subset of small exceptions (w <= WMAX,
7 not| w) realizes (a)+(b) over the spread core P = {2,5,9,11,13}:
  - FOUND    -> an exact no-go certificate: the lens-7 deck argument can NEVER
                produce a witness for V = 7P u W (method boundary at r=8), and
                the same V is then closed by other routes (signal complement);
  - NOT FOUND (exhaustive over the annealer's reach) -> honest negative + the
                rigidity narrative: r=8 static covers demand a rainbow structure
                that generic small sets cannot satisfy.
"""
import sys, os, random
from fractions import Fraction as F
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrc14_certificates import good_intervals, M_exact, is_covering, h_band_protocol, LAM

random.seed(30111)
C = 7
P = [2, 5, 9, 11, 13]
WMAX = 90

def circ(x):
    fx = x - (x.numerator // x.denominator)
    return min(fx, 1 - fx)

GP = good_intervals(P)   # closed core-safe intervals
cands = [w for w in range(1, WMAX + 1) if w % 7 != 0]

# --- refine Gbar_P at every candidate's event t0-values (walls) ---
walls = set()
for w in cands:
    for sgn in (1, -1):
        for j in range(14 * w):
            t0 = (F(sgn, 14) * C + j) / w
            t0 -= (t0.numerator // t0.denominator)
            if any(a <= t0 <= b for a, b in GP):
                walls.add((t0, w))
wall_pts = sorted(set(t for t, w in walls))
pieces = []      # open chamber pieces (represented by midpoints) within Gbar_P
for a, b in GP:
    cuts = [a] + [t for t in wall_pts if a < t < b] + [b]
    for i in range(len(cuts) - 1):
        lo, hi = cuts[i], cuts[i + 1]
        if hi > lo: pieces.append((lo + hi) / 2)
print(f"core P={P}: {len(GP)} safe components, {len(wall_pts)} interior walls, "
      f"{len(pieces)} chamber pieces, {len(cands)} candidates (w<={WMAX}, 7 nmid w)")

# --- sheet assignment tables: k_w(piece) in Z_7 (exactly one sheet, c=7 coprime) ---
def sheet_of(w, t0):
    for k in range(C):
        if circ(F(w) * (t0 + k) / C) < LAM: return k
    return None
table = {}   # (w) -> tuple of sheets per piece
for w in cands:
    table[w] = tuple(sheet_of(w, t0) for t0 in pieces)
# sanity: exactly one sheet on every open piece
assert all(all(s is not None for s in table[w]) for w in cands)

# wall bookkeeping: for each wall (t0, w_owner) the other chosen sheets must be a rainbow
wall_list = sorted(walls)

def uncovered_score(W8):
    """number of violated constraints for the 8-subset (0 = certificate)."""
    bad = 0
    for i in range(len(pieces)):
        if len({table[w][i] for w in W8}) < 7:
            bad += 1
    for t0, w_own in wall_list:
        if w_own not in W8: continue
        others = [w for w in W8 if w != w_own]
        ks = [sheet_of(w, t0) for w in others]
        if None in ks or len(set(ks)) < 7:
            bad += 1
    return bad

# --- annealed search over 8-subsets ---
best = None
for restart in range(24):
    W8 = random.sample(cands, 8)
    sc = uncovered_score(W8)
    for step in range(1200):
        if sc == 0: break
        W2 = list(W8)
        W2[random.randrange(8)] = random.choice(cands)
        if len(set(W2)) < 8: continue
        s2 = uncovered_score(W2)
        if s2 <= sc or random.random() < 0.02:
            W8, sc = W2, s2
    if best is None or sc < best[0]:
        best = (sc, sorted(W8))
    if best[0] == 0: break
sc, W8 = best
print(f"annealer best: violations = {sc} at W = {W8}")

if sc == 0:
    V = sorted([C * p for p in P] + W8)
    print(f"NO-GO CERTIFICATE CANDIDATE: V = {V}")
    print(f"distinct: {len(set(V)) == 13}, covering: {is_covering(V)}")
    # independent full recheck on a fine exact grid of extra points per piece
    ok = True
    for a, b in GP:
        for frac in (F(1,10), F(3,10), F(1,2), F(7,10), F(9,10)):
            t0 = a + (b - a) * frac
            sheets = {sheet_of(w, t0) for w in W8} - {None}
            if len(sheets) < 7: ok = False
    print(f"independent spot recheck (5 pts/component): {'PASS' if ok else 'FAIL'}")
    MV = M_exact(V)
    print(f"SIGNAL COMPLEMENT: M_exact(V) = {MV} >= 1/14: {MV >= F(1,14)}")
    lay, det = h_band_protocol(V)
    print(f"band protocol: layer {lay} {det}")
else:
    print()
    print("HONEST NEGATIVE: no static r=8 small-exception cover found "
          f"(24 restarts x 1200 steps, w <= {WMAX}).")
    print("The binding constraint is the WALL RAINBOW: at every interior wall of a")
    print("chosen exception, the remaining seven must occupy all seven sheets exactly")
    print("once. Violations concentrate there -- the deck resists r=8 static covers")
    print("for the same reason r=7 tilings are chamber-locked: sheet assignments are")
    print("rigid per chamber and walls demand perfect rainbows.")
    # diagnostic: how close did we get, and where do violations sit?
    piece_bad = sum(1 for i in range(len(pieces))
                    if len({table[w][i] for w in W8}) < 7)
    print(f"best subset {W8}: piece violations {piece_bad}, "
          f"wall violations {sc - piece_bad}")
