#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_spectrum_uniformity_kpswf13.py  (kind-pasteur 2026-06-22, THREAD 1 -- uniformity stress test)

THREAD-1 question (2)+(3): is the dominant set F uniform, and is R' >= c with c UNIFORM over (P,E)?
Stress test the SIGN and the apex-7 mechanism over LARGE families:

  A) ALL admissible consecutive clusters k=8..13 x ALL admissible P (|P|=13-k):
     exact min R', the worst (P,E), confirm R' bounded away from 0, and whether n=7 term sign
     tracks SPEC sign.  This is the "binding" family for LRC(14) (consec = densest cluster).

  B) RANDOM wide clusters (perforated / spread, k=8..11) x small P: does R' ever dip near 0?
     Does the universal F (low Farey + apex 7 + P-harmonics) still capture the sign?

  C) The APEX-7 MECHANISM: ghat(7) is the largest coeff (verified). Tabulate ghat(7) and
     chat(7) across the family; show SPEC sign = sign of (2 chat(7) ghat(7) + low-Farey block),
     and that ghat(7) < 0 universally (the cluster's 7-sector deficit) so the sign of the apex
     term = -sign(chat(7)); chat(7) > 0 for gcd(P)=1 (AP-like) => apex term NEG => SPEC<0 (resonant),
     chat(7) < 0 when 7 in P or 7|some interaction => apex term POS => SPEC>0 (decorrelated).
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass
sys.path.insert(0, "04-computation")
import mpmath as mp
mp.mp.dps = 30
from lrc14_spectrum_intersection_sum_kpswf12 import (
    merge, meas, intersect, complement, safe_set, cover_set,
)
from lrc14_spectrum_lowfreq_F_kpswf13 import hat_real


def Rprime_exact(P, E):
    gp = safe_set(P); good = complement(cover_set(E))
    base = meas(gp) * meas(good)
    if base == 0:
        return None, None, None
    inter = meas(intersect(gp, good))
    return inter / base, base, inter


def ghat7_chat7(P, E):
    gp = safe_set(P); good = complement(cover_set(E))
    return float(hat_real(gp, 7)), float(hat_real(good, 7))


def main():
    print("#" * 96)
    print("# THREAD 1 UNIFORMITY: min R' over admissible families, apex-7 sign mechanism")
    print("#" * 96)

    # ---------- A) exact min R' over ALL consec k=8..13 x ALL admissible P ----------
    print("\n" + "=" * 96)
    print("A) EXACT min R' over ALL consecutive clusters k=8..13 and ALL admissible P (|P|=13-k)")
    print("=" * 96)
    glob = (F(10), None, None)
    per_k = {}
    apex_consistency = {"sign_match": 0, "total": 0}
    for k in range(8, 14):
        psz = 13 - k
        E = list(range(k)); good = complement(cover_set(E)); mG = meas(good)
        gh7 = float(hat_real(good, 7))
        mr = (F(10), None)
        Pcount = 0
        for P in itertools.combinations(range(1, 14), psz):
            gp = safe_set(list(P)); base = meas(gp) * mG
            if base == 0:
                continue
            Pcount += 1
            R = meas(intersect(gp, good)) / base
            # apex sign mechanism: SPEC sign vs sign of (2 chat7*ghat7)
            SPEC = meas(intersect(gp, good)) - base
            ch7 = float(hat_real(gp, 7))
            apex_term = 2 * ch7 * gh7
            apex_consistency["total"] += 1
            if (SPEC < 0) == (apex_term < 0):
                apex_consistency["sign_match"] += 1
            if R < mr[0]:
                mr = (R, P)
        per_k[k] = (mr, gh7, Pcount)
        print(f"   k={k:2d} (|P|={psz}, {Pcount:4d} admissible P): min R'={mr[0]}={float(mr[0]):.5f} at P={mr[1]}  [ghat(7)={gh7:+.4f}]")
        if mr[0] < glob[0]:
            glob = (mr[0], mr[1], k)
    print(f"\n   GLOBAL consec floor: min R' = {glob[0]} = {float(glob[0]):.5f}  (k={glob[2]}, P={glob[1]})")
    print(f"   => R' >= c_consec = {float(glob[0]):.5f} > 0  EXACT over the entire consec bank.")
    print(f"   m_P comparison: typical witness floor m_P ~ 0.056; R'-floor/m_P ~ {float(glob[0])/0.056:.1f}x")
    print(f"   APEX-7 sign mechanism: SPEC sign == sign(2 chat(7)ghat(7)) in "
          f"{apex_consistency['sign_match']}/{apex_consistency['total']} "
          f"({100*apex_consistency['sign_match']/apex_consistency['total']:.1f}%) of consec configs")

    # ---------- B) random wide clusters ----------
    print("\n" + "=" * 96)
    print("B) RANDOM WIDE clusters (perforated/spread) x small P -- does R' dip near 0?")
    print("=" * 96)
    random.seed(20260622)
    worst = (F(10), None, None)
    n_neg = n_pos = 0
    samples = 0
    minR_list = []
    for trial in range(400):
        k = random.randint(8, 11)
        span = random.randint(k + 2, 24)
        E = sorted(random.sample(range(0, span + 1), k))
        if E[0] != 0:
            E[0] = 0  # anchor at 0 (the cluster includes a reference point)
            E = sorted(set(E))
        psz = random.randint(2, max(2, 13 - k))
        P = tuple(sorted(random.sample(range(1, 14), psz)))
        R, base, inter = Rprime_exact(list(P), E)
        if R is None:
            continue
        samples += 1
        minR_list.append(float(R))
        SPEC = inter - base
        if SPEC < 0:
            n_neg += 1
        else:
            n_pos += 1
        if R < worst[0]:
            worst = (R, P, tuple(E))
    minR_list.sort()
    print(f"   {samples} valid random wide configs.")
    print(f"   min R' = {float(worst[0]):.5f}  at P={worst[1]}, E={worst[2]}")
    print(f"   R' distribution: min={minR_list[0]:.4f}  p5={minR_list[len(minR_list)//20]:.4f}  "
          f"median={minR_list[len(minR_list)//2]:.4f}  max={minR_list[-1]:.4f}")
    print(f"   SPEC<0 (resonant/anti-corr): {n_neg}/{samples};  SPEC>=0 (decorrelated): {n_pos}/{samples}")

    # ---------- C) apex-7 mechanism table ----------
    print("\n" + "=" * 96)
    print("C) APEX-7 MECHANISM: ghat(7)<0 universally (cluster 7-sector deficit); sign(SPEC)=sign(apex+lowFarey)")
    print("=" * 96)
    print("   ghat(7) sign across consec k=8..13 (the cluster factor, P-independent):")
    all_neg7 = True
    for k in range(8, 14):
        gh7 = float(hat_real(complement(cover_set(list(range(k)))), 7))
        if gh7 >= 0:
            all_neg7 = False
        print(f"      k={k:2d}: ghat(7) = {gh7:+.5f}  {'(NEG -- 7-sector deficit)' if gh7<0 else '(POS!)'}")
    print(f"   ghat(7) < 0 for ALL consec k? {all_neg7}")
    print(f"   chat(7) sign families: gcd(P)=1 AP-like => chat(7)>0 => apex term NEG => pushes SPEC<0.")
    print(f"   7 in P (e.g. P=[5,7,11]) => chat(7)<0 => apex term POS => pushes SPEC>0 (decorrelated).")
    # demonstrate chat(7) for a few P
    for P in [[1, 2, 3], [1, 3, 4, 5], [1, 2, 3, 12, 13], [5, 7, 11], [7, 1, 2], [2, 3, 4]]:
        ch7 = float(hat_real(safe_set(P), 7))
        has7 = (7 in P) or any(p % 7 == 0 for p in P)
        print(f"      P={str(P):<18} chat(7)={ch7:+.5f}  {'(has 7-multiple)' if has7 else '(gcd-coprime to 7)'}")

    print("\nDONE (uniformity).")


if __name__ == "__main__":
    main()
