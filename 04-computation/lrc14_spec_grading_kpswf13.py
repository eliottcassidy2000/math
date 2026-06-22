#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_spec_grading_kpswf13.py   (kind-pasteur 2026-06-22-S36, THREAD 2 part 3 -- the grading)

ON THE TRUE GOOD (GOOD_true = {maxgap{frac(e_i x)}>1/7}, ghat REAL, verified in
lrc14_true_maxgap_good_kpswf13.py), decompose
        SPEC_true = sum_{n!=0} chat(n) ghat_true(n)
and look for a Z_2 grading n -> +/- that splits SPEC into a CONTROLLED part + a (near) MEAN-ZERO
part, per THREAD-2 question (2).  We test, exactly/at high precision:

  G1  parity of n          : SPEC = S[even n] + S[odd n]
  G2  n mod 7              : SPEC = S[7|n] + S[7 nmid n]   (apex prime)
  G3  n mod 2 within 7nmid : finer cross-grading
  G4  the CLUSTER even/odd channel: build GOOD from the complement-even SYMMETRIZATION of the
      phase structure and see if the odd-n SPEC mass tracks the cluster's complement-asymmetry.

OBSERVED in the true-maxgap run: for the symmetric small-P cases the ODD-n SPEC is ~1e-5 (machine
zero up to the 4000-cutoff tail), while even-n carries all of SPEC.  We test whether this is a
THEOREM (odd-n SPEC = 0 for some structural reason) or a near-cancellation, by:
  (i)  computing odd-n SPEC to high cutoff and watching convergence,
  (ii) checking the structural cause: chat(n) for n odd vs the ghat(n) odd-n pattern,
  (iii) the half-integer translate test: SPEC_odd = <c, g> restricted to odd freqs
        = <c, g_odd> where g_odd(x)=(g(x)-g(x+1/2))/2; we compute meas-space
        <1_{G_P}, (1_{GOOD}-1_{GOOD+1/2})/2> EXACTLY (Fraction) and confirm it equals S[odd n].

The EXACT identity used:
  sum_{n odd} hat_c(n) hat_g(n) e(...)  <-> the bilinear form with the ODD-projector
  Podd f (x) = (f(x) - f(x+1/2))/2.   Then  S[odd n] = <1_{G_P}, Podd 1_{GOOD}> in L^2
  = integral 1_{G_P}(x) * (1_{GOOD}(x) - 1_{GOOD}(x+1/2))/2 dx
  = ( meas(G_P cap GOOD) - meas(G_P cap (GOOD - 1/2)) ) / 2     [EXACT Fraction].
  Likewise S[even n] via Peven = (f(x)+f(x+1/2))/2.
This gives an EXACT real-space value for the parity-graded SPEC, no truncation.
"""
import sys
from fractions import Fraction as F
from math import pi
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import (
    merge, meas, intersect, complement, safe_set, fourier_num_of_arcs,
)
from lrc14_true_maxgap_good_kpswf13 import good_true, chat_c

def shift(arcs, s):
    """translate arc-set by +s on the circle (s a Fraction in [0,1))."""
    out = []
    for a, b in arcs:
        a2 = a + s; b2 = b + s
        # reduce mod 1, splitting if it crosses 1
        if b2 <= 1:
            out.append((a2 % 1 if a2 >= 1 else a2, b2))
        else:
            if a2 < 1:
                out.append((a2, F(1)))
                out.append((F(0), b2 - 1))
            else:
                out.append((a2 - 1, b2 - 1))
    return merge([(a, b) for a, b in out if a < b])

def parity_graded_spec_exact(P, E):
    """EXACT parity-graded SPEC via the half-shift projector.
       S_odd = (meas(G_P cap GOOD) - meas(G_P cap (GOOD+1/2)))/2
       S_even = SPEC - S_odd  (and independently (meas(G_P cap GOOD)+meas(G_P cap (GOOD+1/2)))/2
                               - meas(G_P)meas(GOOD), as a check)."""
    gp = safe_set(P); good = good_true(E)
    mGP = meas(gp); mG = meas(good); base = mGP * mG
    if base == 0:
        return None
    good_half = shift(good, F(1, 2))
    i0 = meas(intersect(gp, good))
    ih = meas(intersect(gp, good_half))
    SPEC = i0 - base
    S_odd = (i0 - ih) / 2
    # even-n piece: <1_GP, Peven 1_GOOD> - base where Peven f=(f(x)+f(x+1/2))/2
    #   = (i0 + ih)/2 - base
    S_even = (i0 + ih) / 2 - base
    assert S_even + S_odd == SPEC, (S_even, S_odd, SPEC)
    return dict(SPEC=SPEC, S_even=S_even, S_odd=S_odd, base=base, Rp=i0 / base,
                mGP=mGP, mG=mG, i0=i0, ih=ih, gp=gp, good=good)

def sevens_graded_spec_exact(P, E):
    """EXACT 7-graded SPEC via the 1/7-shift averaging projector.
       P_{7|n} f (x) = (1/7) sum_{j=0}^{6} f(x + j/7).  Then
       S[7|n] = <1_GP, P_{7|n}1_GOOD> - base = (1/7) sum_j meas(GP cap (GOOD+j/7)) - base.
       S[7 nmid] = SPEC - S[7|n]."""
    gp = safe_set(P); good = good_true(E)
    mGP = meas(gp); mG = meas(good); base = mGP * mG
    if base == 0:
        return None
    acc = F(0)
    for j in range(7):
        gj = shift(good, F(j, 7))
        acc += meas(intersect(gp, gj))
    i0 = meas(intersect(gp, good))
    S_sev = acc / 7 - base
    SPEC = i0 - base
    S_nonsev = SPEC - S_sev
    return dict(SPEC=SPEC, S_sev=S_sev, S_nonsev=S_nonsev, base=base)

def main():
    print("#" * 96)
    print("# THREAD 2 part 3: Z_2 / Z_7 GRADING of SPEC_true via EXACT shift-projectors")
    print("#   S_odd = (meas(GP cap GOOD) - meas(GP cap GOOD+1/2))/2   [EXACT Fraction]")
    print("#" * 96)

    bank = [
        ([1, 2, 3, 12, 13], list(range(8)),                       "k=8 consec"),
        ([1, 2, 3, 4, 5],   list(range(8)),                       "k=8 consec smallP"),
        ([1, 2, 3],         list(range(10)),                      "k=10 consec"),
        ([2, 3, 4, 5, 6],   list(range(8)),                       "killer P"),
        ([5, 7, 11],        list(range(10)),                      "coprime P"),
        ([1, 2, 6],         [0, 4, 6, 8, 10, 12, 14, 15, 16, 17], "wide"),
        ([1, 3, 4, 5],      list(range(9)),                       "k=9 consec floor-ish"),
        ([1, 2, 3, 12],     list(range(9)),                       "k=9 consec"),
    ]

    print("\n" + "=" * 96)
    print("(A) PARITY-GRADED SPEC (EXACT via half-shift projector).  S_even + S_odd = SPEC.")
    print("=" * 96)
    print(f"  {'case':<22}{'SPEC':>11}{'S_even':>11}{'S_odd':>11}{'R_true':>9}{'S_odd/SPEC':>11}")
    rows = []
    for P, E, lab in bank:
        r = parity_graded_spec_exact(P, sorted(set(E)))
        if not r:
            continue
        rows.append((lab, P, sorted(set(E)), r))
        ratio = float(r['S_odd'] / r['SPEC']) if r['SPEC'] != 0 else 0.0
        print(f"  {lab:<22}{float(r['SPEC']):>11.6f}{float(r['S_even']):>11.6f}"
              f"{float(r['S_odd']):>11.6f}{float(r['Rp']):>9.5f}{ratio:>11.3f}")
    print("\n  EXACT S_odd values (Fraction):")
    for lab, P, E, r in rows:
        print(f"    {lab:<22} S_odd = {r['S_odd']}   S_even = {r['S_even']}")

    print("\n" + "=" * 96)
    print("(B) 7-GRADED SPEC (EXACT via 1/7-shift averaging projector).  S[7|n] + S[7nmid] = SPEC.")
    print("=" * 96)
    print(f"  {'case':<22}{'SPEC':>11}{'S[7|n]':>11}{'S[7nmid]':>11}{'S7/SPEC':>10}")
    for P, E, lab in bank:
        r = sevens_graded_spec_exact(P, sorted(set(E)))
        if not r:
            continue
        ratio = float(r['S_sev'] / r['SPEC']) if r['SPEC'] != 0 else 0.0
        print(f"  {lab:<22}{float(r['SPEC']):>11.6f}{float(r['S_sev']):>11.6f}"
              f"{float(r['S_nonsev']):>11.6f}{ratio:>10.3f}")

    print("\n" + "=" * 96)
    print("VERDICT on the parity grading:")
    print("=" * 96)
    print("""  S_odd is the EXACT half-shift correlation defect (meas(GP cap GOOD) - meas(GP cap GOOD+1/2))/2.
  If GOOD were 1/2-periodic, S_odd=0.  It is NOT periodic, but the table shows whether |S_odd| is
  uniformly small (a near-mean-zero channel) vs S_even (the controlled bulk).  The SIGN of SPEC is
  set by S_even when |S_odd|<<|S_even|.""")
    print("DONE.")

if __name__ == "__main__":
    main()
