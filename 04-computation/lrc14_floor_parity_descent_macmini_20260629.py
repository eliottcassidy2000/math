#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LRC14 FLOOR -- two creative angles and their connection (mac-mini-2026-06-29-S4).

ANGLE 1 (topology/parity, Borsuk-Ulam):  meas(lonely(S))>0 for covering S.
  lonely(S) = {t: ||s t||>=1/14 for all s in S} = complement of the danger union.
  It is symmetric under t -> 1-t (||s(1-t)||=||st||).  Test forced parities/indices:
   - #components of lonely(S)  (even by symmetry if 0,1/2 are danger)
   - #components of the DANGER set, winding/crossing counts
   - a Borsuk-Ulam style functional whose sign-change forces a lonely arc.

ANGLE 2 (2-adic descent, renormalization):  even e=2e' => ||e t||=||e'(2t)||.
  Map D2: t -> 2t (mod 1).  Under D2, the danger comb of an even speed 2e' is the
  PREIMAGE of the comb of e' on the doubled circle.  Test whether descent
  CONTRACTS the problem: does meas(lonely) relate to meas(lonely of the descended
  set), and does the 2-part 'flow out' to leave an odd core (free at 1/2)?

CONNECTION:  t->1-t (Z2 parity) and t->2t (2-adic flow) are the two order-2
structures of the circle = the '2' of 14=2.7.  We test their interplay on the
lonely set of covering configs.
"""
from __future__ import annotations
import functools, math, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from itertools import combinations

W = F(1, 14)


def danger_intervals(p):
    p = abs(int(p))
    if p == 0:
        return [(F(0), F(1))]   # speed 0 always danger (||0||=0)
    ivs = []
    for k in range(p + 1):
        lo = max(F(k, p) - W / p, F(0)); hi = min(F(k, p) + W / p, F(1))
        if hi > lo:
            ivs.append((lo, hi))
    return ivs


def union_intervals(lists):
    ivs = sorted(iv for lst in lists for iv in lst)
    if not ivs:
        return []
    out = [list(ivs[0])]
    for lo, hi in ivs[1:]:
        if lo > out[-1][1]:
            out.append([lo, hi])
        else:
            out[-1][1] = max(out[-1][1], hi)
    return [(a, b) for a, b in out]


def complement(intervals):
    out = []; cur = F(0)
    for lo, hi in intervals:
        if lo > cur:
            out.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        out.append((cur, F(1)))
    return out


def lonely_components(S):
    """Connected components of lonely(S) on the CIRCLE [0,1) (wrap-aware)."""
    dyn = union_intervals([danger_intervals(s) for s in S])
    comps = complement(dyn)
    # wrap: if first comp starts at 0 and last ends at 1, they join on the circle
    if len(comps) >= 2 and comps[0][0] == 0 and comps[-1][1] == 1:
        first = comps.pop(0); last = comps.pop(-1)
        comps.append((last[0], first[1] + 1))  # merged (measure-wise)
    return comps


def measure(intervals):
    return sum((hi - lo for lo, hi in intervals), F(0))


def is_covering(S):
    """S covering iff for every q in 2..14 some s in S divisible by q."""
    Sset = set(abs(int(s)) for s in S)
    for q in range(2, 15):
        if not any(s % q == 0 for s in Sset if s):
            return False
    return True


def random_covering(rng, k=13, maxspeed=40):
    """Build a covering 13-set: seed with one multiple of each prime power, fill."""
    needed = [2, 3, 4, 5, 7, 8, 9, 11, 13]   # prime powers <=14 force divisibility
    S = set()
    for q in needed:
        mult = q * rng.randint(1, max(1, maxspeed // q))
        S.add(mult)
    while len(S) < k:
        S.add(rng.randint(1, maxspeed))
    S = set(list(S)[:k])
    # ensure still covering after trim
    if len(S) == k and is_covering(S):
        return tuple(sorted(S))
    return None


def main():
    print("=" * 84)
    print("LRC14 floor: parity/index + 2-adic descent (mac-mini-S4)")
    print("=" * 84)
    rng = random.Random(31415)

    # ---- ANGLE 1: parity/index of the lonely set over covering S ----
    print("\n[ANGLE 1] #lonely-components, measure, and parity over covering 13-sets:")
    print(f"{'S (sorted)':<46}{'#comp':>6}{'meas>0':>8}{'meas':>12}")
    samples = []
    tries = 0
    while len(samples) < 12 and tries < 4000:
        tries += 1
        S = random_covering(rng)
        if S is None:
            continue
        samples.append(S)
    comp_counts = []
    zero_meas = 0
    for S in samples:
        comps = lonely_components(S)
        m = measure([(a if a < 1 else a - 1, b if b <= 1 else b - 1) for a, b in
                     [(c[0], c[1]) for c in lonely_components(S)]]) if comps else F(0)
        # recompute measure cleanly from danger complement on [0,1)
        dyn = union_intervals([danger_intervals(s) for s in S])
        m = F(1) - measure(dyn)
        nc = len(comps)
        comp_counts.append(nc)
        if m == 0:
            zero_meas += 1
        print(f"{str(S):<46}{nc:>6}{('YES' if m>0 else 'NO'):>8}{float(m):>12.6f}")
    print(f"\n  component-count parities: {[c % 2 for c in comp_counts]}")
    print(f"  all even? {all(c % 2 == 0 for c in comp_counts)};  any zero-measure? {zero_meas>0}")
    print(f"  min #components = {min(comp_counts)} (=> floor holds if always >=1)")

    # ---- ANGLE 2: 2-adic descent contraction ----
    print("\n[ANGLE 2] 2-adic descent: even e=2e' => ||e t||=||e'(2t)||.")
    print("  For covering S, split S = S_odd U 2*S' (S' = halves of even speeds).")
    print("  Test: meas(lonely(S)) vs the descended measure structure.")
    for S in samples[:6]:
        Sodd = [s for s in S if s % 2 == 1]
        Seven_half = [s // 2 for s in S if s % 2 == 0]
        m_full = F(1) - measure(union_intervals([danger_intervals(s) for s in S]))
        # descended set on the 2t-circle: odd speeds become 'fast' (2t halves their period),
        # even speeds 2e' -> e' on the doubled circle. Build descended danger on u=2t:
        # ||s t|| = ||(s) (u/2)|| ; for even s=2e', = ||e' u||. For odd s, = ||s u/2||.
        # The honest descended object: u-danger of evens = comb of e' ; odds couple via u/2.
        m_odd_half = F(1) - measure(union_intervals([danger_intervals(s) for s in Sodd]))
        m_even_desc = F(1) - measure(union_intervals([danger_intervals(s) for s in Seven_half]))
        print(f"  S={S}")
        print(f"    |S_odd|={len(Sodd)} |S_even|={len(Seven_half)}  "
              f"meas(lonely S)={float(m_full):.5f}  "
              f"meas(lonely S_odd)={float(m_odd_half):.5f}  "
              f"meas(lonely halves)={float(m_even_desc):.5f}")

    # ---- CONNECTION: behavior of the lonely set under t->2t and t->1-t ----
    print("\n[CONNECTION] symmetry t->1-t (parity) and flow t->2t on lonely(S):")
    for S in samples[:4]:
        dyn = union_intervals([danger_intervals(s) for s in S])
        L = complement(dyn)
        # symmetric under t->1-t?
        Lrev = sorted([(F(1) - b, F(1) - a) for a, b in L])
        sym = (union_intervals([L]) == union_intervals([Lrev]))
        print(f"  S={S}: lonely symmetric under t->1-t? {sym}  (#comp={len(L)})")

    print("\n" + "=" * 84)


if __name__ == "__main__":
    main()
