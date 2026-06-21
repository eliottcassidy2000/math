#!/usr/bin/env python3
"""ANGLE F (II) -- decisive test of CA dwell decomposition.

measS7(E) = (number of all-ones intervals) x (mean dwell per interval)
          = n_runs(E) * mean_dwell(E).

HYPOTHESES to test against ALL bounded k-subsets (the THM-535 box):
  (H-RC)  consec MAXIMIZES n_runs(E)  (the "return count" / # visits to all-ones)
  (H-MD)  mean_dwell(E) is (nearly) shape-INVARIANT
  (H-LD)  consec MAXIMIZES longest single dwell
  (H-AT)  consec MAXIMIZES #all-ones atoms (pre-merge)

If (H-RC) holds with zero violators it would be a NEW clean CA characterization
of consec extremality (more local than measS7 itself).  If it fails, record where.

We enumerate primitive, scale/reflection-deduped k-subsets with 0 in E and
max(E) <= span cap, exactly as the verified extremality box.  For speed we use
the same far-box bounds (B2 gives span <= N*(k): 7,8,10 for k=8,9,10) but we scan
the FULL box max(E)<=13 to match the prompt's "0 violators over the full box".

stdlib only; exact Fractions.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import gcd


def breakpoints(E):
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for j in range(0, 7 * e + 1):
            bps.add(F(j, 7 * e))
    return sorted(bps)


def color_of(e, x):
    y = (F(e) * x) % 1
    return (y.numerator * 7) // y.denominator


def dwell_profile(E):
    """Return (measS7, n_runs, longest, n_allones_atoms, n_atoms)."""
    bps = breakpoints(E)
    all_ones = 7
    total = F(0)
    n_runs = 0
    longest = F(0)
    cur = F(0)
    in_run = False
    n_ao_atoms = 0
    n_atoms = 0
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi:
            continue
        n_atoms += 1
        mid = (lo + hi) / 2
        Ncov = len({color_of(e, mid) for e in E})
        w = hi - lo
        if Ncov == all_ones:
            total += w
            cur += w
            n_ao_atoms += 1
            if not in_run:
                n_runs += 1
                in_run = True
        else:
            if in_run and cur > longest:
                longest = cur
            cur = F(0)
            in_run = False
    if in_run and cur > longest:
        longest = cur
    return total, n_runs, longest, n_ao_atoms, n_atoms


def primitive(E):
    g = 0
    for e in E:
        g = gcd(g, e)
    return g == 1


def canon(E):
    """Scale+reflect canonical form for dedup. E has 0, sorted."""
    E = tuple(sorted(E))
    mx = E[-1]
    refl = tuple(sorted(mx - e for e in E))
    # scale-reduce both
    def reduce(t):
        g = 0
        for e in t:
            g = gcd(g, e)
        if g <= 1:
            return t
        return tuple(e // g for e in t)
    a = reduce(E)
    b = reduce(refl)
    return min(a, b)


def enumerate_box(k, span_max):
    """All k-subsets of {0,...,span_max} with 0 and max=given, primitive,
    scale/reflect-deduped."""
    seen = set()
    out = []
    inner = range(1, span_max + 1)
    for rest in combinations(inner, k - 1):
        E = (0,) + rest
        if E[-1] != max(E):  # ensure max present; combinations sorted so ok
            pass
        if not primitive(E):
            continue
        c = canon(E)
        if c in seen:
            continue
        seen.add(c)
        out.append(E)
    return out


def main():
    print("=" * 78)
    print("ANGLE F (II) -- CA return-count / dwell decomposition extremality")
    print("=" * 78)

    SPAN = 13  # full box per prompt
    for k in (8, 9, 10):
        print(f"\n=== k={k}, full box max(E)<=13 ===")
        E0 = tuple(range(k))
        m0, r0, l0, a0, na0 = dwell_profile(E0)
        print(f"consec_{k}: measS7={float(m0):.6f}  n_runs={r0}  "
              f"longest={float(l0):.6f}  ao_atoms={a0}  atoms={na0}")

        shapes = enumerate_box(k, SPAN)
        print(f"  scanning {len(shapes)} deduped primitive shapes...")

        viol_measS7 = []      # measS7 > consec (should be 0 -- known)
        viol_nruns = []       # n_runs > consec
        viol_longest = []     # longest > consec
        viol_ao_atoms = []    # ao_atoms > consec
        # also collect: among the measS7-maximizers, do they share consec's stats?
        best_meas = m0
        meandwell_vals = set()

        for E in shapes:
            if E == E0:
                continue
            m, r, l, a, na = dwell_profile(E)
            md = m / r if r else F(0)
            meandwell_vals.add(md)
            if m > best_meas:
                best_meas = m
            if m > m0:
                viol_measS7.append((E, float(m)))
            if r > r0:
                viol_nruns.append((E, r, float(m)))
            if l > l0:
                viol_longest.append((E, float(l), float(m)))
            if a > a0:
                viol_ao_atoms.append((E, a, float(m)))

        print(f"  measS7 violators (>consec):     {len(viol_measS7)}  "
              f"[expect 0]")
        print(f"  n_runs violators (>consec):     {len(viol_nruns)}")
        print(f"  longest-dwell violators:        {len(viol_longest)}")
        print(f"  ao_atoms violators:             {len(viol_ao_atoms)}")
        print(f"  distinct mean-dwell values:     {len(meandwell_vals)}  "
              f"(H-MD invariance => 1)")
        if viol_nruns[:5]:
            print(f"    n_runs viol examples: {viol_nruns[:5]}")
        if viol_longest[:5]:
            print(f"    longest viol examples: {viol_longest[:5]}")


if __name__ == "__main__":
    main()
