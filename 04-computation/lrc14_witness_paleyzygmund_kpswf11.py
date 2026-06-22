#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_witness_paleyzygmund_kpswf11.py   (kind-pasteur 2026-06-22, THREAD 3 probe)

CREATIVE ANGLE PROBE for the LRC(14) witness floor (OPEN-Q-108 / THM-530 / HYP-2602).

THE OPEN GAP (honest):  the k>=8 branch of the global-witness floor needs
    mu_{1/7}(E) := meas{x in [0,1): maxgap{frac(e_i x): e in E} > 1/7}  >=  thr_k.
The PROVED route (EWLB) bounds mu >= meas(union_a W_a)  (W_a = empty 1/7-window),
but the SYMBOLIC closure "consec minimizes EWLB" is blocked because the windows
are POSITIVELY correlated, so a Bonferroni / union LOWER bound on the union is
not available term-by-term.  (Documented: THM-530 sec E, "cross-window
correlation is positive so Bonferroni fails".)

NEW IDEA TESTED HERE:  POSITIVE correlation is exactly the regime where the
PALEY-ZYGMUND second-moment inequality gives a LOWER bound on P(Z>0):

    P(Z >= 1)  >=  (E Z)^2 / E[Z^2].

Take Z(x) = N(x) = # of the SIX inner sectors {1..6} of Z/7 left EMPTY by the
orbit {frac(e_i x)} (sectors S_j = [ (j-1)/7 ... ] -- we use the canon S7 sectors,
half-open [j/7,(j+1)/7), j=0..6; "inner" excludes the sector containing 0).
An empty inner sector => an orbit gap > 1/7 (canon lemma, THM-530/codex
mu17_sector_cover).  So  {N >= 1} subset {maxgap > 1/7}  and
    mu_{1/7}(E)  >=  P(N >= 1)  >=  (E N)^2 / E[N^2].
E N and E[N^2] are FIRST and SECOND sector-emptiness moments -- computable
EXACTLY and (key hope) of LOW combinatorial complexity (pairwise), so a uniform
bound on them may be reachable where the full measure is not.

WHAT WE MEASURE (exact rational, all):
  1. q_t-distribution of N (t=0..6) for consec_k and a bank of E, exact.
  2. E N = sum_j P(sector j empty);  E[N^2] = sum_{j,j'} P(both j,j' empty).
  3. The Paley-Zygmund bound PZ(E) = (E N)^2 / E[N^2]  and how it compares to
     thr_k and to m_P, and whether consec is the PZ-minimizer.
  4. Does PZ(E) >= m_P hold over the bank?  (m_P = 14249/252252.)
  5. Cantelli one-sided variant and the simpler "E N >= c" first-moment floor.

This is a 5-min scout: rate the angle DEAD/WEAK/PROMISING from the numbers.
"""
import itertools
from fractions import Fraction as Fr

m_P = Fr(14249, 252252)

# thr_k = 1 - min_{|P|=13-k} meas(G_P)  (THM-530); consec mu_{1/7} anchors:
THR = {8: Fr(3637,5880), 9: Fr(2025,4004), 10: Fr(36,91), 11: Fr(25,91), 12: Fr(1,7), 13: Fr(0)}
MU_CONSEC = {8: Fr(691,735), 9: Fr(247,294), 10: Fr(38,49), 11: Fr(1381,2205),
             12: Fr(13823,24255), 13: Fr(477,1078)}

# ---- exact sector-emptiness engine -----------------------------------------
# Sector j = [j/7, (j+1)/7),  j=0..6.  Orbit point e*x mod 1 is in sector j iff
# floor(7*frac(e*x)) == j.  Breakpoints of the whole sector-incidence pattern:
# frac(e*x) crosses a multiple of 1/7  <=>  7 e x in Z  <=>  x = t/(7 e).
def sector_breakpoints(E):
    pts = set()
    for e in E:
        if e == 0:
            continue
        for t in range(0, 7*e + 1):
            pts.add(Fr(t, 7*e))
    pts.add(Fr(0)); pts.add(Fr(1))
    return sorted(p for p in pts if Fr(0) <= p <= Fr(1))

def empty_inner_sectors_at(E, x):
    """Return the SET of inner sectors (j in 1..6) left empty by {frac(e x)}.
       Sector 0 = [0,1/7) always contains e=0's point (0) -> never 'inner empty'."""
    hit = set()
    for e in E:
        j = int((Fr(e)*x % 1) * 7)   # floor(7 frac)
        if j == 7:
            j = 6
        hit.add(j)
    return set(range(1,7)) - hit

def N_distribution_and_moments(E):
    """Exact: P(N=t) for t=0..6, plus E N, E[N^2], and pairwise-empty matrix."""
    bps = sector_breakpoints(E)
    qt = {t: Fr(0) for t in range(7)}
    EN = Fr(0); EN2 = Fr(0)
    # pairwise sector-empty probabilities P(sector j empty AND sector j' empty)
    pair = {(a,b): Fr(0) for a in range(1,7) for b in range(1,7)}
    single = {a: Fr(0) for a in range(1,7)}
    for a, b in zip(bps, bps[1:]):
        if b <= a:
            continue
        mid = (a + b) / 2
        L = b - a
        empt = empty_inner_sectors_at(E, mid)
        n = len(empt)
        qt[n] += L
        EN += n * L
        EN2 += n*n * L
        for s in empt:
            single[s] += L
        for s in empt:
            for s2 in empt:
                pair[(s,s2)] += L
    return qt, EN, EN2, single, pair

def paley_zygmund(E):
    qt, EN, EN2, single, pair = N_distribution_and_moments(E)
    PN1 = sum(qt[t] for t in range(1,7))   # P(N>=1) exact
    PZ = (EN*EN)/EN2 if EN2 != 0 else Fr(0)
    return PN1, PZ, EN, EN2, qt

def consec(k):
    return list(range(k))

def banks(k, span_cap):
    """primitive E, 0 in E, |E|=k, max<=span_cap."""
    out = []
    top = span_cap
    for rest in itertools.combinations(range(1, top+1), k-1):
        E = (0,) + rest
        from math import gcd
        g = 0
        for e in E: g = gcd(g, e)
        if g != 1:
            continue
        out.append(E)
    return out

def fr(x):
    return f"{x}={float(x):.6f}"

print("="*86)
print("LRC(14) WITNESS FLOOR -- PALEY-ZYGMUND (2nd-moment) PROBE  [kpswf11, THREAD 3]")
print("="*86)
print(f"m_P = {fr(m_P)}    (the binding admissible floor, THM-530)")
print()
print("Object: N(x)=# inner sectors {1..6} empty.  {N>=1} subset {maxgap>1/7}.")
print("So mu_{1/7}(E) >= P(N>=1) >= (E N)^2 / E[N^2]  (Paley-Zygmund).")
print()

print("-"*86)
print("[1] CONSEC: exact P(N>=1), Paley-Zygmund bound, vs thr_k and m_P")
print("-"*86)
print(f"{'k':>3} {'P(N>=1)':>12} {'PZ=(EN)^2/EN2':>16} {'thr_k':>10} {'PZ>=thr?':>9} {'PZ>=m_P?':>9} {'E N':>10}")
for k in range(8, 14):
    E = consec(k)
    PN1, PZ, EN, EN2, qt = paley_zygmund(E)
    thr = THR[k]
    print(f"{k:>3} {float(PN1):>12.6f} {float(PZ):>16.6f} {float(thr):>10.6f} "
          f"{str(PZ>=thr):>9} {str(PZ>=m_P):>9} {float(EN):>10.6f}")
    # sanity: P(N>=1) should be <= mu_consec (since {N>=1} subset {maxgap>1/7})
    assert PN1 <= MU_CONSEC[k] + Fr(1,10**6), (k, PN1, MU_CONSEC[k])

print()
print("-"*86)
print("[2] Is PALEY-ZYGMUND BOUND a UNIFORM FLOOR?  min PZ over primitive E bank")
print("-"*86)
print(f"{'k':>3} {'#E':>6} {'min P(N>=1)':>13} {'argmin shape':>26} {'min PZ':>12} {'minPZ>=m_P?':>11}")
for k in range(8, 13):
    span_cap = 14 if k <= 11 else 14
    bank = banks(k, span_cap)
    best_pn1 = None; arg_pn1 = None
    best_pz = None
    for E in bank:
        PN1, PZ, EN, EN2, qt = paley_zygmund(E)
        if best_pn1 is None or PN1 < best_pn1:
            best_pn1 = PN1; arg_pn1 = E
        if best_pz is None or PZ < best_pz:
            best_pz = PZ
    print(f"{k:>3} {len(bank):>6} {float(best_pn1):>13.6f} {str(arg_pn1):>26} "
          f"{float(best_pz):>12.6f} {str(best_pz>=m_P):>11}")

print()
print("-"*86)
print("[3] FIRST-MOMENT floor:  is E N bounded below uniformly?  (E N >= ? )")
print("    E N = sum_{j=1..6} P(sector j empty).  Lowest-complexity functional.")
print("-"*86)
print(f"{'k':>3} {'min E N':>12} {'argmin':>26} {'EN(consec)':>12}")
for k in range(8, 13):
    bank = banks(k, 14)
    best = None; arg = None
    for E in bank:
        qt, EN, EN2, single, pair = N_distribution_and_moments(E)
        if best is None or EN < best:
            best = EN; arg = E
    qtc, ENc, _, _, _ = N_distribution_and_moments(consec(k))
    print(f"{k:>3} {float(best):>12.6f} {str(arg):>26} {float(ENc):>12.6f}")

print()
print("-"*86)
print("[4] WIDE STRESS: does PZ survive far elements (where Bonferroni-EWLB is")
print("    suspect)?  consec_{k-1} + one far element f.")
print("-"*86)
for k in (9, 10, 11):
    base = list(range(k-1))
    print(f"  k={k}:")
    for f in [15, 21, 40, 100]:
        E = base + [f]
        PN1, PZ, EN, EN2, qt = paley_zygmund(E)
        print(f"    f={f:>4}: P(N>=1)={float(PN1):.6f}  PZ={float(PZ):.6f}  "
              f"E N={float(EN):.6f}  PZ>=m_P:{PZ>=m_P}")

print()
print("="*86)
print("VERDICT NOTES printed above; interpretation in session message.")
print("="*86)
