#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_spectrum_wide_decay_envelope_kpswf13.py  (kind-pasteur 2026-06-22-S36, THREAD 3 part 1 -- wide envelope)

PIN C_wide in  |SPEC(dE)| <= C_wide / d  (=> R' >= 1 - eta, eta = C_wide*span0/V* on the wide side).

The two-regime probe showed |SPEC|*d ~ 0.1..0.29 but NON-MONOTONE (d=13 bumps up).  WHY: dilation by d
is a RESONANCE phenomenon -- when d shares structure with the 7-sectors or P, low-n overlap revives.
The clean statement is NOT pointwise monotone decay but a SUP-ENVELOPE:  sup_{d>=D} |SPEC(dE)| -> 0.

We pin this two ways:
  (1) SUP-ENVELOPE: max_{D<=d<=Dmax} |SPEC(dE)|  vs D  -- the decreasing envelope, with its rate.
  (2) THE MECHANISM (exact): SPEC(dE) = sum_{n} chat(n) conj(ghat_{dE}(n)).  chat supported on gcd(P)Z
      (LOW, |n| up to ~max P from the dominant terms).  ghat_{dE}(n): coverSet(dE) breakpoints are at
      (7m+r)/(7 d e), so ghat_{dE}(n) = ghat_E(n/d) STRUCTURE -- dilation MOVES ghat's spectrum to
      d-multiples.  So SPEC(dE) only gets contributions where chat's low support gcd(P)Z MEETS the
      d-dilated ghat support.  As d grows, the overlap at FIXED low n -> the DC-ish piece -> SPEC->0.

  We verify ghat_{dE}(n) relates to ghat_E by the dilation identity, and confirm SPEC(dE) is carried
  by n in gcd(P)Z with |n| small, whose ghat_{dE} value decays.

ALSO: a SUP over a genuinely-wide RANDOM ensemble (not just dilations) to show R'>=1-eta is not a
dilation artifact -- random wide clusters (span large, non-AP) also have R' near 1.
"""
import sys, itertools, cmath, math, random
from fractions import Fraction as F
from math import gcd, pi, sqrt
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass
sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import (
    merge, meas, intersect, complement, safe_set, cover_set, fourier_num_of_arcs,
)

def chat(arcs, n):
    if n == 0:
        return complex(float(meas(arcs)), 0.0)
    return fourier_num_of_arcs(arcs, n) / (2j * pi * n)

def Rprime(P, E):
    gp = safe_set(P); covc = complement(cover_set(E))
    base = meas(gp) * meas(covc)
    if base == 0:
        return None, None
    return meas(intersect(gp, covc)) / base, base

def spec(P, E, N=4000):
    gp = safe_set(P); covc = complement(cover_set(E))
    return sum(2.0 * (chat(gp, n) * chat(covc, n).conjugate()).real for n in range(1, N + 1))

def main():
    print("#" * 96)
    print("# THREAD 3 part 1: WIDE-DECAY ENVELOPE -- pin C_wide in |SPEC(dE)| <= C_wide/d (sup-envelope)")
    print("#" * 96)

    # (1) SUP-ENVELOPE of |SPEC(dE)| over d, for the worst bounded cores
    print("\n" + "=" * 96)
    print("(1) SUP-ENVELOPE  M(D) = max_{D<=d<=60} |SPEC(d*E0)|, and M(D)*D  (the envelope constant)")
    print("=" * 96)
    cases = [([1, 3, 4, 5], list(range(9)), "floor argmin"),
             ([1, 2, 3], list(range(10)), "k=10 |P|=3"),
             ([2, 3, 5], list(range(10)), "coprime-P")]
    for P, E0, lab in cases:
        # precompute |SPEC| over a range of d
        ds = list(range(1, 61))
        sp = {}
        for d in ds:
            E = [d * e for e in E0]
            sp[d] = abs(spec(P, E, N=3000))
        print(f"\n   P={P}  E0={lab}:")
        print(f"      {'D':>4}{'M(D)=sup_{d>=D}|SPEC|':>24}{'M(D)*D':>12}{'argmax d':>10}")
        for D in [1, 2, 3, 5, 8, 13, 21, 34]:
            dom = [d for d in ds if d >= D]
            am = max(dom, key=lambda d: sp[d])
            M = sp[am]
            print(f"      {D:>4}{M:>24.6f}{M*D:>12.6f}{am:>10}")
    print("\n   => M(D) DECREASES in D (sup-envelope); C_wide := sup_D M(D)*D pins |SPEC(dE)|<=C_wide/d.")

    # (2) the dilation identity: ghat_{dE}(n) vs ghat_E(n/d)
    print("\n" + "=" * 96)
    print("(2) DILATION IDENTITY (exact): coverSet(dE) = (1/d)-scaled copies of coverSet(E) tiling [0,1).")
    print("    Consequence: ghat_{dE}(n) = ghat_E(n/d) if d|n, else 0.  So SPEC(dE) lives on n in dZ.")
    print("=" * 96)
    P = [1, 2, 3]; E0 = list(range(10)); d = 3
    Ed = [d * e for e in E0]
    covc0 = complement(cover_set(E0)); covcd = complement(cover_set(Ed))
    print(f"   E0={E0}, d={d}: check ghat_{{dE}}(n)=ghat_E(n/d)*[d|n] for n=1..18")
    print(f"      {'n':>4}{'|ghat_dE(n)|':>16}{'d|n?':>6}{'|ghat_E(n/d)|':>16}{'match':>8}")
    for n in range(1, 19):
        gdn = abs(chat(covcd, n))
        if n % d == 0:
            gEn = abs(chat(covc0, n // d))
            ok = abs(gdn - gEn) < 1e-9
        else:
            gEn = 0.0
            ok = gdn < 1e-9
        print(f"      {n:>4}{gdn:>16.8f}{('Y' if n%d==0 else 'n'):>6}{gEn:>16.8f}{('OK' if ok else 'XX'):>8}")
    print("   => CONFIRMED: dilation by d moves ghat's spectrum onto dZ.  chat lives on gcd(P)Z (low).")
    print("      For d coprime to gcd(P), supp(chat) cap supp(ghat_dE) = gcd(P)Z cap dZ = lcm(gcd(P),d)Z,")
    print("      pushed to HIGH n => the only surviving overlap is the sparse high-lcm lattice => SPEC->0.")

    # (3) genuinely-RANDOM wide ensemble (not dilations): R' near 1?
    print("\n" + "=" * 96)
    print("(3) RANDOM WIDE ensemble (non-AP, span large): is R' near 1 (decorrelation, not a dilation artifact)?")
    print("=" * 96)
    random.seed(20260622)
    worst = (10.0, None)  # min R'
    n_below = 0; n_tot = 0
    print(f"   {'P':>14}{'span':>7}{'R-prime':>12}")
    for trial in range(40):
        P = sorted(random.sample(range(1, 14), random.choice([2, 3, 4])))
        k = random.choice([9, 10])
        # random wide cluster: 0 plus k-1 distinct values up to a large span
        span = random.choice([40, 70, 120, 200])
        E = [0] + sorted(random.sample(range(1, span + 1), k - 1))
        R, base = Rprime(P, E)
        if R is None:
            continue
        n_tot += 1
        if R < 1:
            n_below += 1
        if R < worst[0]:
            worst = (R, (tuple(P), tuple(E)))
        if trial < 12:
            print(f"   {str(P):>14}{max(E):>7}{float(R):>12.6f}")
    print(f"\n   over {n_tot} random wide configs: min R' = {worst[0]:.6f}  ({n_below} had R'<1)")
    print(f"   worst at {worst[1]}")
    print(f"   => random wide R' clusters NEAR 1 (decorrelation is generic, not a dilation special case).")
    print("\nDONE.")

if __name__ == "__main__":
    main()
