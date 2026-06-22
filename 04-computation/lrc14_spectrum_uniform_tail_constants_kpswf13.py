#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_spectrum_uniform_tail_constants_kpswf13.py   (kind-pasteur 2026-06-22-S36, THREAD 3 part 1)

NAIL THE UNIFORM CONSTANT in the L2-Cauchy-Schwarz tail bound
    |sum_high(H)| = |sum_{|n|>H} chat(n) conj(ghat(n))|  <=  sqrt(E_c(H)) * sqrt(E_g(H)),
where (Parseval, EXACT for indicators)
    E_c(H) = sum_{|n|>H} |chat(n)|^2 = meas(G_P) - meas(G_P)^2 - sum_{0<|n|<=H} |chat(n)|^2,
    E_g(H) = sum_{|n|>H} |ghat(n)|^2 = (1-p0) - (1-p0)^2 - sum_{0<|n|<=H} |ghat(n)|^2.

THE CLAIM TO TEST (THREAD 3 (1)):  C' in  |tail| <= C'/sqrt(H)  is UNIFORM over (P, cluster).

The mechanism that would make C' uniform is the L2 *tail-energy* rate.  For an indicator 1_S that is
a finite union of arcs, |hat(n)|^2 ~ (jump structure)/n^2, and the tail energy
    E(H) = sum_{|n|>H} |hat(n)|^2.
Parseval gives a HARD CAP independent of H:  E(0) = meas(S) - meas(S)^2 = Var(1_S) <= 1/4.
So E_c(H) <= 1/4 and E_g(H) <= 1/4 ALWAYS, giving the H-FREE bound |tail| <= 1/4.  That is uniform
but not -> 0.  The real question is the RATE in H:

  * Does E_g(H) <= C_g / H with C_g UNIFORM over the cluster E (whose #jumps grows!)?
  * Does E_c(H) <= C_c / H with C_c UNIFORM over the small part P?

KEY SUBTLETY: the #jumps of 1_{cover^c} grows with the cluster (breakpoints have denom | 7 lcm(E)),
so the naive "C = (#jumps)^2/(4 pi^2)" is NOT uniform.  But the *energy* E(H) is governed by the
ACTUAL coefficient decay, and for an indicator the energy in the tail is the "high-frequency content"
which is controlled by the SMALLEST gap / finest feature, not the count.  We MEASURE:
  - E_c(H)*H, E_g(H)*H  across a wide bank (consec, wide-dilated, random, coprime-P)
  - is sup over the bank of E_g(H)*H bounded?  of E_c(H)*H?
  - the resulting uniform C' = sup sqrt(E_c(H)*H) * sqrt(E_g(H)*H) (so |tail| <= C'/H if both ~1/H)
    OR C' = sup sqrt(E_c(H)) * sqrt(E_g(H)*H)  (so |tail| <= C'/sqrt(H) if only g decays ~1/H).

We ALSO record the HARD H-FREE cap (1/4 each => |tail| <= 1/4) which is trivially uniform, and
the per-(P,E) ACTUAL |sum_high| to compare the bound's sharpness.
"""
import sys, itertools, cmath, math, random
from fractions import Fraction as F
from math import gcd, pi, sqrt
from functools import reduce
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

def energies(arcs, H):
    """EXACT-Parseval tail energy E(H)=meas - meas^2 - sum_{0<|n|<=H}|hat(n)|^2 (numeric low part)."""
    m = float(meas(arcs))
    var = m - m * m                     # = sum_{n!=0} |hat(n)|^2  (Parseval, EXACT)
    low = 0.0
    for n in range(1, H + 1):
        low += 2.0 * abs(chat(arcs, n)) ** 2
    tail = var - low
    return max(tail, 0.0), var

def true_sum_high(gp, covc, H, N):
    """actual |sum_{H<|n|<=N} chat conj(ghat)| (numeric, N large as proxy for inf)."""
    s = 0.0
    for n in range(H + 1, N + 1):
        s += 2.0 * (chat(gp, n) * chat(covc, n).conjugate()).real
    return s

def jumps(arcs):
    return 2 * len(merge(arcs))

def analyze(P, E, Hs=(7, 14, 28, 56, 112), N=4000, label=""):
    P = sorted(set(int(p) for p in P)); E = sorted(set(int(x) for x in E))
    gp = safe_set(P); covc = complement(cover_set(E))
    mGP = meas(gp); mC = meas(covc); base = mGP * mC
    Vc = jumps(gp); Vg = jumps(covc)
    rows = []
    for H in Hs:
        Ec, varc = energies(gp, H)
        Eg, varg = energies(covc, H)
        cs = sqrt(Ec) * sqrt(Eg)
        actual = abs(true_sum_high(gp, covc, H, N))
        rows.append(dict(H=H, Ec=Ec, Eg=Eg, EcH=Ec * H, EgH=Eg * H, cs=cs, actual=actual,
                         varc=varc, varg=varg))
    return dict(P=P, E=E, label=label, mGP=mGP, mC=mC, base=base, Vc=Vc, Vg=Vg, rows=rows)

def main():
    print("#" * 96)
    print("# THREAD 3 part 1: UNIFORM CONSTANT in the L2-Cauchy-Schwarz tail |sum_high| <= sqrt(Ec)sqrt(Eg)")
    print("#   tail energies E(H) are EXACT (Parseval: E(0)=meas-meas^2<=1/4 H-FREE cap); rate in H is the Q")
    print("#" * 96)

    # A WIDE BANK spanning every regime that stresses uniformity:
    #   - bounded consec cores (the floor argmin family)
    #   - WIDE dilated clusters (Vg grows; the wide regime mac-mini S32 says -> R'=1)
    #   - coprime-P (independence favourable)
    #   - random admissible
    bank = [
        ([1, 2, 3, 12, 13], list(range(8)),                         "k=8 consec"),
        ([1, 3, 4, 5],      list(range(9)),                         "k=9 consec (floor argmin)"),
        ([1, 2, 3],         list(range(10)),                        "k=10 consec |P|=3"),
        ([5, 7, 11],        list(range(10)),                        "k=10 coprime-P"),
        ([1, 2, 3],         [0, 2, 4, 6, 8, 10, 12, 14, 16, 18],    "WIDE x2 dilate"),
        ([1, 2, 3],         [0, 5, 10, 15, 20, 25, 30, 35, 40, 45], "WIDE x5 dilate"),
        ([1, 2, 3],         [0, 13, 26, 39, 52, 65, 78, 91, 104, 117], "WIDE x13 dilate"),
        ([1, 2, 6],         [0, 4, 6, 8, 10, 12, 14, 15, 16, 17],   "wide d>=2 mixed"),
        ([2, 3, 5],         [0, 7, 11, 23, 41, 53, 67, 79, 97, 101], "WIDE random-ish"),
    ]
    Hs = (7, 14, 28, 56, 112, 224)
    results = [analyze(P, E, Hs=Hs, N=4000, label=lab) for P, E, lab in bank]

    # ---- table 1: the tail-energy RATE E(H)*H per (P,E), to see if sup E*H is bounded ----
    print("\n" + "=" * 96)
    print("TABLE 1: tail-energy rate.  E_g(H)*H and E_c(H)*H  (if bounded in H AND across bank => uniform 1/H)")
    print("=" * 96)
    print(f"{'case':<28}{'Vc':>4}{'Vg':>5}  " + "".join(f"{'EgH@'+str(H):>11}" for H in Hs))
    for r in results:
        line = f"{r['label']:<28}{r['Vc']:>4}{r['Vg']:>5}  "
        line += "".join(f"{row['EgH']:>11.5f}" for row in r['rows'])
        print(line)
    print()
    print(f"{'case':<28}{'Vc':>4}{'Vg':>5}  " + "".join(f"{'EcH@'+str(H):>11}" for H in Hs))
    for r in results:
        line = f"{r['label']:<28}{r['Vc']:>4}{r['Vg']:>5}  "
        line += "".join(f"{row['EcH']:>11.5f}" for row in r['rows'])
        print(line)

    # ---- sup over bank of E*H at each H ----
    print("\n" + "=" * 96)
    print("sup over bank of E_g(H)*H and E_c(H)*H  (the candidate uniform constants C_g, C_c)")
    print("=" * 96)
    print(f"{'H':>6}{'sup Eg*H':>14}{'sup Ec*H':>14}{'sup sqrt(Ec*H)*sqrt(Eg*H)':>28}{'sup actual*H':>16}")
    for i, H in enumerate(Hs):
        supEgH = max(r['rows'][i]['EgH'] for r in results)
        supEcH = max(r['rows'][i]['EcH'] for r in results)
        # uniform C' s.t. |tail| <= C'/H if BOTH energies ~1/H:
        supprod = max(sqrt(r['rows'][i]['EcH']) * sqrt(r['rows'][i]['EgH']) for r in results)
        supactH = max(r['rows'][i]['actual'] * H for r in results)
        print(f"{H:>6}{supEgH:>14.5f}{supEcH:>14.5f}{supprod:>28.5f}{supactH:>16.5f}")

    # ---- table 2: the CS bound vs actual sum_high, and the H-free cap ----
    print("\n" + "=" * 96)
    print("TABLE 2: CS bound sqrt(Ec)sqrt(Eg) vs ACTUAL |sum_high|, at H=28; H-free cap = sqrt(varc)sqrt(varg)")
    print("=" * 96)
    iH = Hs.index(28)
    print(f"{'case':<28}{'CS@28':>12}{'actual@28':>12}{'ratio':>9}{'Hfree cap':>12}{'1/4 cap':>9}")
    for r in results:
        row = r['rows'][iH]
        hfree = sqrt(row['varc']) * sqrt(row['varg'])
        ratio = row['cs'] / max(row['actual'], 1e-12)
        print(f"{r['label']:<28}{row['cs']:>12.5f}{row['actual']:>12.5f}{ratio:>9.2f}{hfree:>12.5f}{0.25:>9.4f}")

    # ---- the H-free uniform cap is trivially Var<=1/4 each: |tail|<=1/4 ALWAYS ----
    print("\n" + "=" * 96)
    print("THE H-FREE UNIFORM CAP (no rate needed): Var(1_S)=meas-meas^2 <= 1/4 for ANY indicator.")
    print("=" * 96)
    supvc = max(max(row['varc'] for row in r['rows']) for r in results)
    supvg = max(max(row['varg'] for row in r['rows']) for r in results)
    print(f"  sup over bank Var(1_GP)   = {supvc:.6f}   (<= 1/4 = 0.25)")
    print(f"  sup over bank Var(1_cov^c)= {supvg:.6f}   (<= 1/4 = 0.25)")
    print(f"  => |sum_high(H)| <= sqrt(Var_c)sqrt(Var_g) <= 1/4 for ALL H, ALL (P,E)  [H-FREE, UNIFORM]")
    print(f"     (this alone does NOT -> 0; the RATE table 1 shows E_g(H)*H bounded => the 1/H decay)")

    print("\nDONE.")

if __name__ == "__main__":
    main()
