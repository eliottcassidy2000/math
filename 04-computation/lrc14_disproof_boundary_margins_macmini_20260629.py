#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""How thin is the line between proof and disproof of LRC(14)?
(mac-mini-2026-06-29-S15)

LRC(14): for 13 nonzero speeds S, M(S) = sup_t min_i ||v_i t|| >= 1/14 (a lonely time exists).
A DISPROOF = a config with M(S) < 1/14 (danger combs cover [0,1), lonely set EMPTY).

This script measures the RAZOR-THIN LINE:
 (1) the extremal {1,...,13}: M=1/14 EXACTLY, lonely set = isolated points (measure 0, NONEMPTY).
     -> the proof/disproof boundary IS this measure-zero-but-nonempty locus.
 (2) margins of structured covering families: how close to 1/14 from above?
 (3) the sector/residue pattern at the tight time t*=1/14 (the apex-7 / mod-14 structure).
 (4) rigidity: perturbing the extremal -- does M ever dip below 1/14? (it must not, if LRC holds.)
"""
from __future__ import annotations
import functools, math, itertools
print = functools.partial(print, flush=True)


def frac(x): return x - math.floor(x)
def dist(x):  # ||x|| = distance to nearest integer
    f = frac(x); return min(f, 1 - f)


def M_of(speeds, grid=2_000_004, refine=True):
    """M(S) = max_t min_i ||v_i t||, t in (0,1). Fine grid + local refine at the argmax."""
    best_t, best_v = 0.0, -1.0
    inv = 1.0 / grid
    for k in range(1, grid):
        t = k * inv
        mn = 1.0
        for v in speeds:
            d = dist(v * t)
            if d < mn:
                mn = d
                if mn <= best_v: break
        if mn > best_v:
            best_v, best_t = mn, t
    if refine:
        lo, hi = best_t - inv, best_t + inv
        for _ in range(60):
            m1 = lo + (hi - lo) / 3; m2 = hi - (hi - lo) / 3
            f1 = min(dist(v * m1) for v in speeds)
            f2 = min(dist(v * m2) for v in speeds)
            if f1 < f2: lo = m1
            else: hi = m2
        tm = (lo + hi) / 2
        vm = min(dist(v * tm) for v in speeds)
        if vm > best_v: best_v, best_t = vm, tm
    return best_v, best_t


def main():
    thr = 1.0 / 14
    print("=" * 78)
    print(f"LRC(14): threshold 1/14 = {thr:.6f}. DISPROOF needs M(S) < 1/14.")
    print("=" * 78)

    # (1) the extremal {1,...,13}
    S0 = list(range(1, 14))
    M0, t0 = M_of(S0)
    print(f"\n(1) EXTREMAL S={{1,...,13}}: M = {M0:.6f}, argmax t* = {t0:.6f} "
          f"(= {round(t0*14)}/14 = {round(t0*14)/14:.6f})")
    print(f"    M - 1/14 = {M0 - thr:+.6e}  => TIGHT (equality). LRC holds at EXACTLY 1/14.")
    # find all near-lonely times (min dist within eps of 1/14) -- the lonely 'set'
    eps = 2e-4
    lonely_pts = []
    g = 14 * 6  # sample a/14 and refined
    for a in range(1, 14):
        t = a / 14
        md = min(dist(v * t) for v in S0)
        if md >= thr - 1e-9:
            who = [i + 1 for i, v in enumerate(S0) if abs(dist(v * t) - thr) < 1e-9]
            lonely_pts.append((f"{a}/14", round(md, 6), who))
    print(f"    lonely points (min dist >= 1/14) at a/14: {lonely_pts}")
    print(f"    => the lonely SET is isolated rational points (measure 0, NONEMPTY). The boundary.")

    # (2) margins of structured covering families
    print(f"\n(2) MARGINS of structured 13-speed families (M and M-1/14):")
    fams = {
        "{1..13}":            list(range(1, 14)),
        "{1..12,14}":         list(range(1, 13)) + [14],
        "{2,3,...,14} (1+i)": list(range(2, 15)),
        "odds {1,3,..,25}":   list(range(1, 26, 2)),
        "{1..6,8..14} (no7)": [v for v in range(1, 15) if v != 7],
        "all even-ish 2*{1..13}":[2 * v for v in range(1, 14)],   # = 2*{1..13}, scales: M same as {1..13}
        "{1,2,3,5,6,9,11,...}": [1,2,3,5,6,9,11,13,15,17,19,21,23],
        "{1..7,9..14} (no8)": [v for v in range(1, 15) if v != 8],
    }
    for name, S in fams.items():
        if len(set(S)) != 13:
            print(f"    {name:24s}: (not 13 distinct, skip)"); continue
        M, t = M_of(S, grid=1_000_002)
        flag = "  <-- BELOW 1/14 (DISPROOF?!)" if M < thr - 1e-6 else ("  = tight" if abs(M-thr)<1e-5 else "")
        print(f"    {name:24s}: M={M:.6f}  M-1/14={M-thr:+.5f}{flag}")

    # (3) the residue/sector pattern at t*=1/14 (apex-7 / mod-14)
    print(f"\n(3) SECTOR pattern at t*=1/14 for {{1..13}}: v_i*14*t = i mod 14")
    res = [(i, i % 14, "DANGER(0)" if (i % 14) == 0 else
            ("boundary 1/14" if dist(i/14) == thr else "safe")) for i in range(1, 14)]
    odd = [i for i in range(1, 14) if i % 2 == 1]; even = [i for i in range(1, 14) if i % 2 == 0]
    print(f"    odd residues  (apex-7 face): {odd}  -> incl boundary runners 1,13 (the touching pair)")
    print(f"    even residues (2*apex face): {even}  = 2*{{1,2,3,4,5,6}} (the 2-adic-descended core)")
    print(f"    the EXACT-TOUCHING runners at 1/14 are i=1 and i=13 (= +-1 mod 14, the units).")
    print(f"    => a disproof must KILL the t=1/14 lonely point = cover where runners 1,13 touch.")

    # (4) rigidity: integer perturbations of one speed -- does M ever drop below 1/14?
    print(f"\n(4) RIGIDITY: replace one speed of {{1..13}} by a larger integer; track M vs 1/14:")
    below = 0
    for j in range(13):
        for newv in [14, 15, 20, 27, 41, 99]:
            S = list(range(1, 14)); S[j] = newv
            if len(set(S)) != 13: continue
            M, _ = M_of(S, grid=600_002)
            if M < thr - 1e-6:
                below += 1
                print(f"    {{1..13}} with v_{j+1}->{newv}: M={M:.6f} < 1/14  <-- candidate!")
    print(f"    perturbations dipping below 1/14: {below}  (0 expected if LRC holds)")
    print(f"    (note: grid resolution ~1e-6; a true sub-1/14 needs exact certification, not grid.)")

    print("\n" + "=" * 78)
    print("THE LINE: M(S)=1/14 is achieved with EQUALITY and an ISOLATED (measure-0) lonely set")
    print("at the extremal {1..13}, touching at runners +-1 mod 14 (the units). A disproof must")
    print("push M strictly below 1/14 = destroy that isolated touching point. The proof target is")
    print("the STRICT inequality at the measure-zero boundary, not a bulk measure bound.")
    print("=" * 78)


if __name__ == "__main__":
    main()
