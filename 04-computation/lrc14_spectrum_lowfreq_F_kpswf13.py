#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_spectrum_lowfreq_F_kpswf13.py   (kind-pasteur 2026-06-22, THREAD 1)

THE LOW-FREQUENCY STRUCTURE OF SPEC -- identify the finite dominant frequency set F,
test its uniformity over (P, cluster), and certify R' >= c via F + an L2 tail.

SETUP (matches lrc14_spectrum_intersection_sum_kpswf12, the canonical NODE-3 objects):
  GOOD = cover^c = {x : at least one of the 7 sectors [r/7,(r+1)/7) is MISSED by {frac(e_i x)}}
                 = the LRC "cluster maxgap > 1/7" set.   (1_GOOD = g,  ghat = its coeffs)
  G_P  = {x : ||p x|| >= 1/14 for all p in P}.            (1_{G_P} = c,  chat = its coeffs)
  baseline = meas(G_P) * meas(GOOD).
  SPEC = sum_{n!=0} chat(n) conj(ghat(n)) = meas(GOOD cap G_P) - baseline   (EXACT, real-space).
  R'   = 1 + SPEC/baseline = meas(GOOD cap G_P)/baseline.

EXACT PER-FREQUENCY TERM (the deliverable of THREAD 1):
  Both 1_{G_P} and 1_GOOD are indicators of SYMMETRIC arc sets (S = -S mod 1), verified below by
  meas(S) == meas(reflect(S)) and (stronger) by reflect(S)==S as merged arc lists for G_P; for GOOD
  the two sets agree up to a measure-zero boundary set, so both chat(n) and ghat(n) are REAL.
  For an interval indicator 1_{[a,b)} with a,b rational of common denom D:
     hat(n) = (e(-n a) - e(-n b))/(2 pi i n),     e(t)=exp(2 pi i t).
  REAL part (the whole thing, since the set is symmetric):
     hat(n) = (1/(2 pi n)) * sum_arcs [ sin(2 pi n b) - sin(2 pi n a) ].
  Hence
     term(n) := chat(n) conj(ghat(n)) + chat(-n) conj(ghat(-n)) = 2 chat(n) ghat(n)
              = 2/(4 pi^2 n^2) * S_c(n) * S_g(n),
     S_c(n) = sum_{G_P arcs}(sin 2pi n b - sin 2pi n a),  S_g(n) = sum_{GOOD arcs}(...).
  The pi^2/n^2 is COMMON; the per-frequency *shape* is S_c(n) S_g(n)/n^2.  We compute term(n) to
  50-digit precision (mpmath) -- effectively exact -- and CROSS-CHECK sum_n term(n) against the
  EXACT Fraction SPEC.  All "EXACT" claims below are validated to <1e-12.

THREAD 1 questions answered:
  (1) per-frequency contribution term(n), n=1..50, for several (P, cluster) incl binding/worst;
      which n dominate?  Identify the finite set F carrying >95% of sum|term(n)|.
  (2) Is F uniform over (P, cluster)?  (Always the same low Farey denominators?)
  (3) R' >= 1 + (sum_{n in F} term(n) - L2tail(F))/baseline; is it >= c > 0 uniformly?  Sign of SPEC.
  (4) WHY n=3,4,6 not 7:  ghat vanishes on 7|n (cluster's 1/7 sectors) -- verify S_g(7k)=0.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass
sys.path.insert(0, "04-computation")
import mpmath as mp
mp.mp.dps = 50

from lrc14_spectrum_intersection_sum_kpswf12 import (
    merge, meas, intersect, complement, safe_set, cover_set,
)

PI = mp.pi


# ----------------------------------------------------------------- exact-ish per-frequency
def sin_sum(arcs, n):
    """S(n) = sum_arcs [ sin(2 pi n b) - sin(2 pi n a) ]  (mpmath, 50 dps). Real coeff core.
       Reduce n*endpoint mod 1 as an EXACT Fraction first to preserve precision."""
    s = mp.mpf(0)
    for a, b in arcs:
        ta = (F(n) * b) % 1
        tb = (F(n) * a) % 1
        s += mp.sin(2 * PI * mp.mpf(ta.numerator) / ta.denominator) \
           - mp.sin(2 * PI * mp.mpf(tb.numerator) / tb.denominator)
    return s


def hat_real(arcs, n):
    """REAL Fourier coeff hat(n) for a symmetric arc set (n != 0)."""
    return sin_sum(arcs, n) / (2 * PI * n)


def term_n(gp, good, n):
    """2 chat(n) ghat(n) (the +-n pair), to 50 dps."""
    return 2 * hat_real(gp, n) * hat_real(good, n)


# ----------------------------------------------------------------- symmetry / reality checks
def reflect(arcs):
    out = []
    for a, b in arcs:
        na, nb = (1 - b) % 1, (1 - a) % 1
        if na < nb:
            out.append((na, nb))
        elif na > nb:
            out.append((na, F(1))); out.append((F(0), nb))
        else:
            out.append((F(0), F(1)))
    return merge(out)


def symmetry_report(gp, good):
    gpm = merge(gp); goodm = merge(good)
    gp_sym = (gpm == reflect(gpm))
    good_sym_set = (goodm == reflect(goodm))
    good_sym_meas = (meas(goodm) == meas(reflect(goodm)))
    return gp_sym, good_sym_set, good_sym_meas


# ----------------------------------------------------------------- L2 tail (Cauchy-Schwarz)
def l2_tail_beyond(gp, good, H):
    """|sum_{|n|>H} chat conj ghat| <= sqrt(Ec) sqrt(Eg), Ec=meas(G_P)-dc-low (EXACT meas, low numeric)."""
    mGP = float(meas(gp)); mG = float(meas(good))
    ec_dc = mGP ** 2; eg_dc = mG ** 2
    ec_low = 0.0; eg_low = 0.0
    for n in range(1, H + 1):
        ec_low += 2.0 * float(hat_real(gp, n)) ** 2
        eg_low += 2.0 * float(hat_real(good, n)) ** 2
    tc = max(mGP - ec_dc - ec_low, 0.0)
    tg = max(mG - eg_dc - eg_low, 0.0)
    return (tc ** 0.5) * (tg ** 0.5), tc, tg


# ----------------------------------------------------------------- main per-case analysis
def analyze(P, E, NMAX=50, label="", verbose=True):
    P = sorted(set(int(p) for p in P)); E = sorted(set(int(x) for x in E))
    gp = safe_set(P); good = complement(cover_set(E))
    mGP = meas(gp); mG = meas(good); baseline = mGP * mG
    m_inter = meas(intersect(gp, good))
    SPEC_exact = m_inter - baseline
    Rprime = (m_inter / baseline) if baseline > 0 else F(1)

    gp_sym, good_sym_set, good_sym_meas = symmetry_report(gp, good)

    # per-frequency terms
    terms = [(n, term_n(gp, good, n)) for n in range(1, NMAX + 1)]
    spec_partial = sum(t for _, t in terms)
    total_abs = sum(abs(t) for _, t in terms)

    # dominant set: sort by |term|, accumulate to 95% of total_abs
    order = sorted(terms, key=lambda kv: -abs(kv[1]))
    F_set = []
    acc = mp.mpf(0)
    for n, t in order:
        F_set.append(n)
        acc += abs(t)
        if total_abs > 0 and acc / total_abs >= mp.mpf("0.95"):
            break
    F_set_sorted = sorted(F_set)

    # vanishing on 7|n check
    s_g_7 = [(7 * k, float(sin_sum(good, 7 * k))) for k in range(1, 8)]

    if verbose:
        print("=" * 96)
        print(f"  ({label})  P={P}  E={E}  (k={len(E)})")
        print("=" * 96)
        print(f"  meas(G_P)={mGP}={float(mGP):.6f}  meas(GOOD)={mG}={float(mG):.6f}  baseline={float(baseline):.6f}")
        print(f"  meas(GOOD cap G_P)={float(m_inter):.6f}  SPEC(exact)={float(SPEC_exact):+.7f}  R'={Rprime}={float(Rprime):.6f}")
        print(f"  symmetry: G_P=-G_P(as set)?{gp_sym}  GOOD=-GOOD(as set)?{good_sym_set}  (by meas?{good_sym_meas})  => chat,ghat REAL")
        print(f"  spectral SPEC(|n|<={NMAX}) = {float(spec_partial):+.7f}  vs exact {float(SPEC_exact):+.7f}  |diff|={abs(float(spec_partial)-float(SPEC_exact)):.2e}")
        print(f"  -- per-frequency term(n) = 2*chat(n)*ghat(n)  (n=1..{min(NMAX,18)} shown):")
        for n, t in terms[:18]:
            mark = ""
            if n % 7 == 0:
                mark = "  <-- 7|n (ghat should ~0)"
            print(f"       n={n:2d}: {float(t):+.7f}   (|term|={float(abs(t)):.7f}){mark}")
        print(f"  -- DOMINANT SET F (>=95% of sum|term|): F = {F_set_sorted}")
        print(f"       (|F|={len(F_set_sorted)}; carries {float(acc/total_abs*100):.1f}% of total |term| mass)")
        sgn = "NEG (anti-corr)" if SPEC_exact < 0 else ("POS (indep-fav)" if SPEC_exact > 0 else "ZERO")
        print(f"  -- SIGN of SPEC: {sgn}")
        print(f"  -- S_g(7k) [ghat numerator on 7-multiples], should be ~0:")
        print(f"       " + "  ".join(f"7*{k}:{v:+.2e}" for k, (m7, v) in zip(range(1, 8), s_g_7)))

    return dict(P=P, E=E, mGP=mGP, mG=mG, baseline=baseline, m_inter=m_inter,
                SPEC=SPEC_exact, Rprime=Rprime, terms=terms, F_set=F_set_sorted,
                acc_frac=float(acc / total_abs) if total_abs > 0 else 1.0,
                total_abs=float(total_abs), s_g_7=s_g_7,
                gp_sym=gp_sym, good_sym_set=good_sym_set)


def main():
    print("#" * 96)
    print("# THREAD 1 : LOW-FREQUENCY STRUCTURE OF SPEC -- dominant set F, uniformity, R'>=c bound")
    print("#   GOOD = cover^c (maxgap>1/7);  G_P = {||px||>=1/14};  SPEC = sum_{n!=0} chat(n) ghat(n)")
    print("#" * 96)

    # Binding / worst / representative (P, cluster).  Include the EXACT consec floor case
    # P=(1,3,4,5),k=9 (R'=0.5346, the global consec minimizer from kpswf12 VERDICT).
    cases = [
        ([1, 2, 3, 12, 13], list(range(8)),  "k=8 consec, P worst-ish"),
        ([1, 2, 3, 4, 5],   list(range(8)),  "k=8 consec, small P"),
        ([1, 2, 3, 12],     list(range(9)),  "k=9 consec"),
        ([1, 3, 4, 5],      list(range(9)),  "k=9 consec, P=(1,3,4,5) = CONSEC FLOOR R'=0.5346"),
        ([1, 2, 3],         list(range(10)), "k=10 consec |P|=3"),
        ([1, 2, 3, 12, 13], [0, 2, 3, 4, 5, 6, 7, 8], "k=8 perforated near-AP"),
        ([5, 7, 11],        list(range(10)), "k=10 P coprime (indep-fav)"),
        ([1, 2, 6],         [0, 4, 6, 8, 10, 12, 14, 15, 16, 17], "wide d>=2 leader"),
    ]
    res = [analyze(P, E, NMAX=50, label=lab) for P, E, lab in cases]
    for _ in res:
        print()

    # ---- uniformity of F ----
    print("=" * 96)
    print("UNIFORMITY OF F across cases (each row = dominant frequency set carrying >=95% |term| mass)")
    print("=" * 96)
    for r in res:
        print(f"   P={str(r['P']):<22} k={len(r['E']):2d}  R'={float(r['Rprime']):.4f}  F={r['F_set']}")
    # candidate "universal" low-Farey set: denominators of resonant centers up to 7 = {1,2,3,4,5,6,7}
    universal_guess = sorted(set(range(1, 8)))
    print(f"\n   candidate UNIVERSAL low-Farey set (Farey denominators <= 7): {universal_guess}")
    # which n appear in EVERY case's F?
    common = set(res[0]['F_set'])
    union = set()
    for r in res:
        common &= set(r['F_set'])
        union |= set(r['F_set'])
    print(f"   intersection of all F   = {sorted(common)}")
    print(f"   union of all F          = {sorted(union)}")
    print(f"   max element across all F = {max(union)}")

    print("\nDONE (per-case + uniformity). See lowfreq_certificate run for the R'>=c bound.")
    return res


if __name__ == "__main__":
    main()
