#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Working the team questions (mac-mini-S35 broadcast): does the 2-adic DESCENT land the WORST covering on
the minimal-gap (DOUBLET) core, and what is the tournament-side image of the binding doublet C_7?
mac-mini-2026-06-30-S36.

Q2 (descent lands worst on doublet): the descent (THM-580: S=O u E, S'=E/2, recurse; cores = O_j mod 7)
produces a CHAIN of Z_7-cores. The per-level apex floor of each core is THM-590's gap; the 5 values are
{0, 0.198=4cos^2(3pi/7), 0.308, 1, 2}. The total floor ~ product of gaps over the chain. So the WORST
(minimal per-level) factor a covering can have is the DOUBLET gap 0.198 -- a covering achieves it IFF its
descent chain contains a 2-residue (doublet) core. Verify: the binding R={1..13}\{7} hits a doublet; sweep
which coverings do; confirm THM-590 (no positive gap < 0.198).

Q1 (the matching): the doublet's autocorrelation = 2I + A(C_7); the binding atom is the 7-CYCLE C_7, an
EVEN graph = the cusp of the even-graph dual E_7 -- NOT a single strongly-connected tournament. The
tournament side has the 3-cycle (minimal cyclicity, Z_3); the apex binding object is its EVEN-GRAPH DUAL at
length 7 (HYP-3590). So the 'tournament-side image' is a category clarification: the atom is dual/cycle-space.
"""
from __future__ import annotations
import functools, math, cmath
print = functools.partial(print, flush=True)
W = cmath.exp(2j * math.pi / 7)


def gap(core):
    O = set(x % 7 for x in core)
    if not O: return None                       # empty core: no apex factor at this level
    if O == {0, 1, 2, 3, 4, 5, 6}: return 0.0   # full Z_7 = the gap-0 cusp (measure-0, existence carries)
    return min(abs(sum(W**(k * x) for x in O))**2 for k in range(1, 7))


def descend(S):
    """2-adic descent: return the chain of odd-cores (mod 7)."""
    cores = []
    S = sorted(set(S))
    while S:
        O = [x for x in S if x % 2 == 1]
        E = [x for x in S if x % 2 == 0]
        if O: cores.append(sorted(set(o % 7 for o in O)))
        if not E: break
        S = sorted(set(e // 2 for e in E))
    return cores


def chain_report(S):
    cores = descend(S)
    gaps = [(c, gap(c)) for c in cores]
    pos = [g for _, g in gaps if g is not None and g > 1e-9]
    min_gap = min(pos) if pos else 0.0
    hit_doublet = any(len(c) == 2 for c in cores)
    hit_full = any(g == 0.0 for _, g in gaps)
    return cores, min_gap, hit_doublet, hit_full


def main():
    DOUBLET = 4 * math.cos(3 * math.pi / 7)**2
    print("=" * 84)
    print(f"DOES THE DESCENT LAND THE WORST COVERING ON THE DOUBLET? (4cos^2(3pi/7)={DOUBLET:.4f})")
    print("=" * 84)

    # ---- (1) the binding covering and the drop-one family ----
    print("\n[1] R = {1..13}\\{x}, drop each x -- the descent chain's min per-level gap:")
    print(f"    {'drop x':>6} {'depth':>6} {'min per-level gap':>18} {'doublet?':>9} {'full-Z7(gap0)?':>14}")
    best = (None, 9.9)
    for x in range(1, 14):
        R = [v for v in range(1, 14) if v != x]
        cores, mg, hd, hf = chain_report(R)
        if 0 < mg < best[1]: best = (x, mg)
        print(f"    {x:>6} {len(cores):>6} {mg:>18.4f} {('YES' if hd else 'no'):>9} {('YES' if hf else 'no'):>14}")
    print(f"    => the worst BOUNDED (positive) per-level floor over drop-one = {best[1]:.4f} at x={best[0]} "
          f"({'=4cos^2(3pi/7), the DOUBLET' if abs(best[1]-DOUBLET)<1e-3 else ''}).")
    print("    (drop x=7 = the apex prime = klein-S8's binding R; its chain hits the doublet {1,3}.)")

    # ---- (2) broader sweep: consec prefixes, dense, and the tightest coverings ----
    print("\n[2] broader sweep -- min per-level gap and whether the chain hits a doublet:")
    fams = {
        "consec {1..13}": list(range(1, 14)),
        "consec {1..12}": list(range(1, 13)),
        "tightest {1..12,182}": list(range(1, 13)) + [182],
        "skip-12 {1..11,13,84}": [v for v in range(1, 14) if v != 12] + [84],
        "even-heavy {2,4,..,26}": list(range(2, 27, 2)),
        "{1..13}\\{7} (binding)": [v for v in range(1, 14) if v != 7],
    }
    for name, S in fams.items():
        cores, mg, hd, hf = chain_report(S)
        coresz = [len(c) for c in cores]
        print(f"    {name:>24}: core sizes {coresz}, min gap {mg:.4f}, doublet={'Y' if hd else 'N'}, "
              f"full-Z7={'Y' if hf else 'N'}")

    # ---- (3) THM-590 check: no positive gap below the doublet ----
    print("\n[3] THM-590 (per-level): over ALL 127 nonempty Z_7-cores, the positive gaps and their MIN:")
    import itertools
    vals = set()
    for sz in range(1, 8):
        for O in itertools.combinations(range(7), sz):
            g = gap(O)
            if g is not None: vals.add(round(g, 4))
    posmin = min(v for v in vals if v > 1e-9)
    print(f"    distinct gap values: {sorted(vals)}; MIN POSITIVE = {posmin} = 4cos^2(3pi/7) "
          f"(at the doublets); gap=0 only at full Z_7.")
    print("    => no covering's chain can have a positive per-level factor below the doublet 0.198 (THM-590).")

    # ---- (4) Q1: the tournament-side image of the doublet ----
    print("\n[4] Q1 -- the matching (tournament-side image of the binding doublet C_7):")
    print("    The doublet {a,b} depends only on d=b-a; its autocorrelation = 2I + A(C_7), C_7=Cay(Z_7,{+-d}).")
    print("    So the binding atom = the 7-CYCLE C_7 = an EVEN graph = the CUSP of the even-graph dual E_7,")
    print("    NOT a strongly-connected TOURNAMENT (a tournament connection set has size 3, the doublet is 2).")
    print("    Tournament side has the 3-CYCLE (minimal cyclicity, Z_3); the apex binding object is its")
    print("    EVEN-GRAPH DUAL at length 7 (HYP-3590). The 'image' is a CATEGORY clarification: the binding")
    print("    atom lives on the cycle-space/dual side; matching it to one SC tournament is a category error.")

    print("\n" + "=" * 84)
    print("ANSWERS: Q2 -- the WORST bounded covering descends THROUGH a DOUBLET (min per-level gap = 0.198 =")
    print("4cos^2(3pi/7), THM-590, attained at drop-x=7 = klein-S8's binding R); deeper/full-Z_7 chains are")
    print("the gap-0 cusp (measure-0, existence carries, klein-S16). The total floor is the PRODUCT over the")
    print("chain (-> 0 for deep chains); the per-level binding atom is the doublet. Q1 -- the binding atom is")
    print("the EVEN-GRAPH 7-cycle C_7 (cusp of E_7), not a tournament; tournament side has the 3-cycle, its")
    print("length-7 dual. So the descent reduces the floor to ONE irreducible atom: the doublet/C_7.")
    print("=" * 84)


if __name__ == "__main__":
    main()
