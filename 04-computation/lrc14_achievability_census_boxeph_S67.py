#!/usr/bin/env python3
"""
THE (phi, o) ACHIEVABILITY CENSUS (boxeph-2026-07-17-S67)
Finishing sweep, part 2: LEM-036's named open -- which survivor sets
S_r subset Z_7 are achievable by the affine six-line arrangements
v_f(j) = phi_f j + o_f (survivor condition: 0 not in v, {1..5} <= v)?

Strata enumerated EXHAUSTIVELY (all 7^6 = 117,649 offset vectors each):
  P6:  phi = (1,2,3,4,5,6)      -- the permutation stratum (family clusters)
  C1:  phi = (1,2,3,4,5,0)      -- one constant coordinate (two-large shape)
  R1:  phi = (1,1,2,3,4,5)      -- a repeated slope
  C2:  phi = (1,2,3,4,0,0)      -- two constants (balanced-21 shape)
Plus the all-constant stratum phi = (0,...,0): v = o constant in j, so
S = Z_7 iff o covers ({1..5} <= o, 0 not in o) else S = empty -- the ONLY
way |S| = 7 is achieved (no moving coordinate), per LEM-036(B).

Questions answered:
  (1) the |S|-distribution per stratum; re-verification of the full-column
      theorem by exhaustion (|S| = 6 only in P6, exactly at o = t phi;
      |S| = 7 never with a moving coordinate);
  (2) WHICH subsets are achievable: counts of distinct S, and orbit count
      under the affine action j -> aj + b (order 42) which acts within the
      achievable family;
  (3) achievability extremes: are all singletons/pairs/... achievable?
      what is the largest non-full achievable |S|< 6 per stratum?
"""

import sys
from itertools import product

FIVE = frozenset(range(1, 6))


def census(phi, label):
    ach = {}
    full = []
    for o in product(range(7), repeat=6):
        S = []
        for j in range(7):
            v = [(phi[i] * j + o[i]) % 7 for i in range(6)]
            if 0 not in v and FIVE <= set(v):
                S.append(j)
        fs = frozenset(S)
        ach[fs] = ach.get(fs, 0) + 1
        if len(S) >= 6:
            full.append((o, fs))
    hist = {}
    for fs, cnt in ach.items():
        hist[len(fs)] = hist.get(len(fs), 0) + cnt
    distinct = {}
    for fs in ach:
        distinct.setdefault(len(fs), set()).add(fs)
    print(f"  [{label}] phi = {phi}:")
    print(f"      configs by |S|: " + ", ".join(
        f"{k}: {hist.get(k, 0)}" for k in range(8) if hist.get(k)))
    print(f"      distinct S by |S|: " + ", ".join(
        f"{k}: {len(distinct.get(k, ()))}" for k in range(8)
        if distinct.get(k)))
    return ach, full, distinct


def affine_orbits(sets):
    """orbit count of a family of subsets of Z_7 under j -> a j + b."""
    seen = set()
    orbits = 0
    for S in sets:
        if S in seen:
            continue
        orbits += 1
        for a in range(1, 7):
            for b in range(7):
                seen.add(frozenset((a * j + b) % 7 for j in S))
    return orbits


if __name__ == "__main__":
    print("THE (phi, o) ACHIEVABILITY CENSUS -- LEM-040 referee (boxeph S67)")
    print("=" * 78)
    strata = [((1, 2, 3, 4, 5, 6), "P6 permutation"),
              ((1, 2, 3, 4, 5, 0), "C1 one-constant"),
              ((1, 1, 2, 3, 4, 5), "R1 repeat"),
              ((1, 2, 3, 4, 0, 0), "C2 two-constant")]
    all_ach = {}
    for phi, label in strata:
        ach, full, distinct = census(phi, label)
        moving = any(p != 0 for p in phi)
        for o, fs in full:
            assert len(fs) < 7, (label, "seven-column with moving coord!")
        if label.startswith("P6"):
            # full-column theorem by exhaustion: |S| = 6 iff o = t*phi
            sixes = [(o, fs) for o, fs in full if len(fs) == 6]
            pred = {tuple((t * p) % 7 for p in (1, 2, 3, 4, 5, 6))
                    for t in range(7)}
            assert {o for o, _ in sixes} == pred, "full-column mismatch"
            print(f"      FULL-COLUMN THEOREM by exhaustion: the {len(sixes)}"
                  f" configs with |S| = 6 are exactly o = t phi, t = 0..6")
        else:
            assert not any(len(fs) >= 6 for _, fs in full), (label,)
        for fs in ach:
            all_ach.setdefault(fs, 0)
        print()
    # all-constant stratum (the only |S| = 7 route)
    n7 = sum(1 for o in product(range(7), repeat=6)
             if 0 not in o and FIVE <= set(o))
    print(f"  [C6 all-constant] phi = 0: S = Z_7 for the {n7} covering "
          f"offsets, else empty -- the ONLY |S| = 7 route (LEM-036(B))")
    print()
    nonempty = [S for S in all_ach if S]
    print(f"  ACROSS MOVING STRATA: {len(all_ach)} distinct achievable S "
          f"({len(nonempty)} nonempty); affine orbits: "
          f"{affine_orbits(nonempty)}")
    total = 2 ** 7
    print(f"  (subset space: {total}; nonempty achievable fraction "
          f"{len(nonempty) / (total - 1) * 100:.0f}% across these strata)")
    print("=" * 78)
    print("done")
