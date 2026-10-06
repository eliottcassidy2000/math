#!/usr/bin/env python3
"""Task D (chessboard weave, LRC reading), 2026-10-06: the parity ("scaffold") weave.
Exact arithmetic throughout.

v = O u E (odd / even speeds), E/2 = {e/2 : e in E}.
D1  exact values for the 8x8 board: O = {1,3,5,7}, E = {2,4,6,8}, v = {1..8}.
D2  weave identity (proved below) and the tested candidate inequalities
      (I1) delta(v) <= min(delta(O), delta(E/2))                       [trivial]
      (I2) delta(v) >= delta(E/2)/2                                      [candidate]
      (I3) 1/delta(v) <= |O| + 1/delta(E/2)   (1/delta(empty) := 1)      [candidate;
           iterating it down the 2-adic tower would imply LRC]
      (I4) 1/delta(v) <= 1/delta(O) + 1/delta(E/2)                       [candidate]
    over the census universes (primitive n-subsets of {1..B}):
      n=2 B=200, n=3 B=60, n=4 B=30, n=5 B=18
D3  the same quantities on every tight set found in the A-script (n <= 7).
Reproduce: python3 chessboard_weave_20261006_lrc_D_parity.py > chessboard_weave_20261006_lrc_D_parity.out
"""
import sys
import time
from fractions import Fraction
from functools import lru_cache

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from chessboard_weave_20261006_lrc_core import (delta, delta_argmax, safe_components,
                                                primitive_sets, check)


@lru_cache(maxsize=None)
def dl(v):
    return delta(v)


def split(v):
    O = tuple(x for x in v if x % 2)
    E2 = tuple(x // 2 for x in v if x % 2 == 0)
    return O, E2


def inv_delta(v):
    """1/delta, with the convention 1/delta(empty) = 1 (LRC with 0 speeds)."""
    return Fraction(1) if not v else 1 / dl(v)


def is_tight(v):
    comps, D = safe_components(v)
    return bool(comps) and all(lo == hi for lo, hi in comps)


T0 = time.time()
print("=" * 78)
print("D1. The board's two scaffolds and their weave (exact).")
for name, v in [("odd lengths  O", (1, 3, 5, 7)), ("even lengths E", (2, 4, 6, 8)),
                ("E/2          ", (1, 2, 3, 4)), ("union {1..8}  ", tuple(range(1, 9)))]:
    d, pts = delta_argmax(v)
    n = len(v)
    print(f"  {name} = {v}: delta = {d}; 1/(n+1) = {Fraction(1, n + 1)}; "
          f"tight: {d == Fraction(1, n + 1) and is_tight(v)}; argmax in (0,1/2]: {[str(p) for p in pts]}")
check(dl((1, 3, 5, 7)) == Fraction(1, 2), "odd")
check(dl((2, 4, 6, 8)) == Fraction(1, 5) and is_tight((2, 4, 6, 8)), "even")
check(dl(tuple(range(1, 9))) == Fraction(1, 9) and is_tight(tuple(range(1, 9))), "union")
print("  verified: delta(O) = 1/2, delta(E) = delta(1,2,3,4) = 1/5 (tight, n=4), "
      "delta({1..8}) = 1/9 (tight, n=8).")

print()
print("=" * 78)
print("D2. Weave identity (PROVED, one line).  For u in [0,1) the two preimages of u under")
print("    t -> 2t are u/2 and u/2 + 1/2; even speeds 2e give ||e u|| at both, and for odd o")
print("    ||o(u/2 + 1/2)|| = 1/2 - ||o u/2||.  Hence")
print("      delta(v) = max_u min{ f_{E/2}(u),  max( min_o ||o u/2||,  1/2 - max_o ||o u/2|| ) }.")
print("    Consequences: delta(O) = 1/2 for every nonempty set O of odd speeds (t = 1/2, i.e. u = 0")
print("    on the shifted branch),")
print("    so 'delta of the odd part' carries no information, and (I1) reads delta(v) <= delta(E/2).")
print()
print("    Census test of (I1)-(I4):")
U = {2: 200, 3: 60, 4: 30, 5: 18}
summary = {}
for n, B in U.items():
    t1 = time.time()
    num = 0
    viol = {"I1": [], "I2": [], "I3": [], "I4": []}
    eq_I1 = 0
    worst_I2 = None      # min delta(v)/delta(E/2)
    worst_I3 = None      # max 1/delta(v) - |O| - 1/delta(E/2)
    noE = 0
    for v in primitive_sets(n, B):
        num += 1
        d = dl(v)
        O, E2 = split(v)
        if not E2:
            noE += 1
            check(d == Fraction(1, 2), v)
        else:
            dE = dl(E2)
            # (I1)
            if d > dE or (O and d > Fraction(1, 2)):
                viol["I1"].append(v)
            if d == dE:
                eq_I1 += 1
            # (I2)
            r = d / dE
            if worst_I2 is None or r < worst_I2[0]:
                worst_I2 = (r, [v])
            elif r == worst_I2[0] and len(worst_I2[1]) < 12:
                worst_I2[1].append(v)
            if 2 * d < dE:
                viol["I2"].append(v)
            # (I4)
            if O and 1 / d > 2 + 1 / dE:
                viol["I4"].append(v)
        # (I3)
        X = 1 / d - len(O) - inv_delta(E2)
        if worst_I3 is None or X > worst_I3[0]:
            worst_I3 = (X, [v])
        elif X == worst_I3[0] and len(worst_I3[1]) < 12:
            worst_I3[1].append(v)
        if X > 0:
            viol["I3"].append(v)
    summary[n] = (num, viol, worst_I2, worst_I3)
    print(f"  n={n}, B={B}: {num} primitive sets ({time.time() - t1:.1f}s); "
          f"{noE} all-odd (delta = 1/2 for all of them)")
    print(f"    (I1) violations: {len(viol['I1'])}; equality delta(v) = delta(E/2) in {eq_I1} sets")
    print(f"    (I2) min delta(v)/delta(E/2) = {worst_I2[0]} at {worst_I2[1][:8]}; "
          f"violations (ratio < 1/2): {len(viol['I2'])} {viol['I2'][:8]}")
    print(f"    (I3) max [1/delta(v) - |O| - 1/delta(E/2)] = {worst_I3[0]} at {worst_I3[1][:8]}; "
          f"violations (> 0): {len(viol['I3'])} {viol['I3'][:8]}")
    print(f"    (I4) violations of 1/delta(v) <= 2 + 1/delta(E/2): {len(viol['I4'])}; "
          f"first: {viol['I4'][:8]}")

print()
print("=" * 78)
print("D3. The same quantities on the tight sets of the A-script (delta = 1/(n+1)).")
TIGHT = [(1, 2), (1, 2, 3), (1, 2, 3, 4), (1, 3, 4, 7), (1, 2, 3, 4, 5), (1, 3, 4, 5, 9),
         (1, 2, 3, 4, 5, 6), (1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 7, 12),
         (1, 4, 5, 6, 7, 11, 13), tuple(range(1, 9))]
print(f"  {'v':<26} {'O':<18} {'E/2':<14} {'delta(v)':>8} {'delta(E/2)':>10} "
      f"{'I2 ratio':>8} {'I3 excess':>9}")
for v in TIGHT:
    d = dl(v)
    check(d == Fraction(1, len(v) + 1) and is_tight(v), v)
    O, E2 = split(v)
    dE = dl(E2) if E2 else None
    r = d / dE if E2 else None
    X = 1 / d - len(O) - inv_delta(E2)
    print(f"  {str(v):<26} {str(O):<18} {str(E2):<14} {str(d):>8} {str(dE):>10} "
          f"{str(r):>8} {str(X):>9}   E/2 tight: {is_tight(E2) if E2 else '-'}")
print("  (I2 needs ratio >= 1/2; I3 needs excess <= 0.)")
print()
print(f"Total time {time.time() - T0:.1f}s")
