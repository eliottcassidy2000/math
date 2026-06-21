#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROUTE 3 v2 (opus, 2026-06-21): sharpen the per-cell Huffer-Shepp reflection port.

v1 found:
  * FULL reflection E->-E is an EXACT measure-preserving symmetry: W_a(E)=W_{7-a}(-E),
    measS7 invariant.  (Clean, provable -- the HS "reflect all arcs = symmetry" step.)
  * The full-coverage COMPONENT containing y=0 is ONE-SIDED with reach exactly
    1/(7 max|e|): the fastest clock first breaks coverage at the cell midpoint.
    So y=0 is generically a COVERAGE BOUNDARY, not an interior centre.
  * Single-clock reflection from non-consec rivals only LOWERED measS7.

This v2 nails the per-cell monotonicity question with the RIGHT target:

  (A) PER-CELL HS MONOTONICITY (the literal route-3 deliverable).  For each cell a,
      is consec the per-cell maximiser  W_a(consec) >= W_a(E)  for all E in the
      stratum?  If YES per cell, summing gives measS7 trivially (a MUCH stronger,
      cleaner statement than the global one).  If NO, find the cell + rival where a
      non-consec shape BEATS consec on a single W_a -- the precise obstruction.

  (B) REFLECTION-TOWARD-CONSEC.  HS transfer: a single move toward the balanced
      (consec) config does not decrease coverage.  Define the move = take a shape
      E and reflect ONE clock e->-e ONLY IF it brings the magnitude multiset
      closer to consec's {0..k-1} (a "Robin-Hood toward consec").  Does W_a not
      decrease?  Test exhaustively per cell.

  (C) THE BINDING-SECTOR REFLECTION.  Within cell a, the coverage breaks when the
      fastest-relevant clock vacates its sector.  Reflecting that clock's drift
      (e->-e) sends its covered y-window from one side of 0 to the other.  Test:
      does the UNION of the two reflected copies (the "symmetrised" clock, both
      drift directions) cover MORE?  This is the cleanest HS reflection: replace a
      one-directional covering arc by its symmetric two-sided version.
"""
import itertools, math, sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def W_a(E, a):
    E = sorted(set(E))
    lo = F(a, 7) - F(1, 14); hi = F(a, 7) + F(1, 14)
    bps = {lo, hi}
    for e in E:
        if e == 0: continue
        d = 7 * abs(e)
        j0 = math.floor(lo * d); j1 = math.ceil(hi * d)
        for j in range(j0 - 1, j1 + 2):
            x = F(j, d)
            if lo <= x <= hi: bps.add(x)
    bps = sorted(bps); tot = F(0)
    for l, h in zip(bps, bps[1:]):
        if h <= l: continue
        xm = (l + h) / 2; hit = set()
        for e in E:
            v = e * xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        if len(hit) == 7: tot += h - l
    return tot

def Wvec(E): return [W_a(E, a) for a in range(1, 7)]
def measS7(E): return sum(W_a(E, a) for a in range(1, 7))
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

if __name__ == "__main__":
    print("=" * 84)
    print("ROUTE 3 v2: PER-CELL Huffer-Shepp monotonicity of W_a")
    print("=" * 84)

    # ============== (A) PER-CELL maximality of consec ==============
    print("\n[A] IS consec THE PER-CELL MAXIMISER?  W_a(consec) >= W_a(E) for ALL E, each a?")
    print("    (If yes, the wall collapses to 6 independent per-cell inequalities.)")
    print("=" * 84)
    for k, span in [(8, 14), (9, 14), (10, 14)]:
        C = consec(k); WC = Wvec(C)
        # full-residue stratum
        bank = []
        for combo in itertools.combinations(range(1, span + 1), k - 1):
            E = (0,) + combo
            if not primitive(E): continue
            if set(e % 7 for e in E) != {0,1,2,3,4,5,6}: continue
            bank.append(E)
        print(f"\n  k={k}: consec W_a = {[str(w) for w in WC]}")
        print(f"        stratum size = {len(bank)}")
        beaten = defaultdict(list)   # cell a -> list of (E, W_a(E), W_a(consec))
        for E in bank:
            if list(E) == C: continue
            for ai, a in enumerate(range(1, 7)):
                we = W_a(E, a)
                if we > WC[ai]:
                    beaten[a].append((E, we, WC[ai]))
        any_beat = False
        for a in range(1, 7):
            if beaten[a]:
                any_beat = True
                bb = max(beaten[a], key=lambda t: t[1] - t[2])
                print(f"        cell a={a}: consec BEATEN on W_a by {len(beaten[a])} shapes; "
                      f"worst {bb[0]} W_a={float(bb[1]):.5f} vs consec {float(bb[2]):.5f} "
                      f"(+{float(bb[1]-bb[2]):.5f})")
        if not any_beat:
            print("        >>> consec is the PER-CELL MAXIMISER on every cell a=1..6.  "
                  "PER-CELL HS HOLDS at this k.")
        else:
            print("        >>> consec is NOT a per-cell maximiser: PER-CELL HS port FAILS at this k.")

    # ============== (B) symmetrised (two-sided) covering arc ==============
    print("\n" + "=" * 84)
    print("[B] HS SYMMETRISATION.  Replace each non-zero clock e by the PAIR {e,-e}")
    print("    (a symmetric two-sided covering arc about the cell centre).  HS centering")
    print("    lemma => symmetrising should NOT DECREASE coverage.  Compare W_a(E) to")
    print("    W_a(E union -E).  Test sign per cell.")
    print("=" * 84)
    for E in [consec(8), [0,1,2,3,4,5,6,8], [0,2,3,4,5,6,7,8], [0,1,2,3,4,5,7,8]]:
        Esym = sorted(set(E) | set(-e for e in E))
        print(f"\n  E={E}  ->  symmetrised E∪(-E)={Esym}")
        for a in range(1, 7):
            w0 = W_a(E, a); w1 = W_a(Esym, a)
            tag = "+" if w1 > w0 else ("-" if w1 < w0 else "0")
            print(f"     a={a}: W_a(E)={float(w0):.5f}  W_a(sym)={float(w1):.5f}  d={float(w1-w0):+.5f} [{tag}]")

    # ============== (C) reflection TOWARD consec, per cell ==============
    print("\n" + "=" * 84)
    print("[C] REFLECTION-TOWARD-CONSEC single moves.  From a rival E, reflect each")
    print("    clock e->-e; for the moves that reduce |magnitude - consec target|,")
    print("    record per-cell sign of dW_a.  HS-port predicts: moving toward the")
    print("    balanced consec config does NOT decrease W_a.")
    print("=" * 84)
    Ctgt = set(range(8))
    for E in [[0,1,2,3,4,5,6,8], [0,2,3,4,5,6,7,8]]:
        print(f"\n  rival E={E}")
        for idx in range(len(E)):
            e = E[idx]
            if e == 0: continue
            E2 = list(E); E2[idx] = -e
            if not primitive(E2): continue
            # does reflection bring residue closer to canonical? (informational)
            for a in range(1, 7):
                d = W_a(E2, a) - W_a(E, a)
                tag = "+" if d > 0 else ("-" if d < 0 else "0")
                if a == 1:
                    print(f"     reflect e={e:>3}->{-e:>3}:", end="")
                print(f" a{a}:{tag}", end="")
            print(f"   dmeasS7={float(measS7(E2)-measS7(E)):+.5f}")
