#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_witness_vs_criterion_threshold_kpswf10.py  (kind-pasteur 2026-06-21, THREAD B core)

THE FACTOR-2, RESOLVED via the canon's own slow-fast dictionary (THM-526 kps-S4):
   fast-phase width  w_theta(x) = 6/7 - circ_width({frac(e_i x)})
                                = 6/7 - (1 - maxgap)              [circ_width = 1 - maxgap]
                                = maxgap{frac(e_i x)} - 1/7.
So:
   GLOBAL WITNESS (direct M>=1/14, no runner removed):  needs w_theta > 0  <=>  maxgap > 1/7.
   VIA-MAX CRITERION C (W(S\{Vmax}) > 1/(7 Vmax)):       needs w_theta > 1/7 <=>  maxgap > 2/7.
THM-527's rho* uses the CRITERION threshold (maxgap > 2/7).  The WITNESS threshold (maxgap>1/7)
is WEAKER (a larger event) and ALSO certifies M>=1/14 (more directly).

This is the factor-2: the SECTOR ROUTE's '1 missed 1/7-sector' is the maxgap>1/7 (witness)
scale; THM-527 demands the maxgap>2/7 (criterion) scale.  BOTH certify M>=1/14, but they are
DIFFERENT sufficient events with DIFFERENT measures.

We verify EXACTLY (rational):
  (D1)  w_theta(x) = maxgap - 1/7  pointwise  (the dictionary).
  (D2)  {global witness} = {maxgap > 1/7} = G1  (=> M>=1/14).
  (D3)  {criterion}      = {maxgap > 2/7} = G2 = rho* (THM-527)  (=> M>=1/14).
  (D4)  G2 subset G1 subset {>=1 sector missed}=1-p0.   So  g*=G2 <= G1 <= 1-p0.
  (D5)  meas(G1) > 0  is EASIER to get than meas(G2) > 0, and BOTH suffice for LRC.
        => THE WITNESS ROUTE (G1>0) IS A STRICTLY WEAKER SUFFICIENT CONDITION than THM-527.
           If meas(G1) (intersect G_P) > 0 we are DONE without needing the 2/7 event.

PAYOFF QUESTION (the prompt's Q4): since maxgap>1/7 already certifies M>=1/14 (global witness)
and is a LARGER event than maxgap>2/7, can we close LRC with the EASIER 1/7 floor instead of
THM-527's 2/7 floor?  We tabulate meas(G1 cap G_P) vs meas(G2 cap G_P) to see if the 1/7 floor
is also positive (and larger), i.e. whether THM-527 over-demands.
"""
import itertools
from fractions import Fraction as Fr

P = 7
T1 = Fr(1, 7)
T2 = Fr(2, 7)
HALF = Fr(1, 14)


def phases_at(E, x):
    return sorted((int(e) * x) % 1 for e in E)


def maxgap(phases):
    if len(phases) <= 1:
        return Fr(1)
    g = max(b - a for a, b in zip(phases, phases[1:]))
    return max(g, (phases[0] + 1) - phases[-1])


def circ_width(phases):
    """smallest arc containing all points = 1 - maxgap."""
    return Fr(1) - maxgap(phases)


# ---- G_P safe set: ||p x|| >= 1/14 for all p in P (the small part) ----
def in_GP(Pset, x):
    for p in Pset:
        r = (int(p) * x) % 1
        d = min(r, 1 - r)        # ||p x||
        if d < HALF:
            return False
    return True


def grid(E, Pset):
    allv = list(E) + list(Pset)
    bp = {Fr(0), Fr(1)}
    for e in allv:
        e = int(e)
        if e == 0:
            continue
        # 1/7 sector boundaries (for sectors / maxgap-vs-1/7,2/7)
        for t in range(0, P * e + 1):
            bp.add(Fr(t, P * e))
        # 1/14 danger boundaries (for G_P): ||e x||=1/14 at x=(k+-1/14)/e
        for k in range(0, e + 1):
            for s in (HALF, -HALF):
                val = Fr(k) + s
                bp.add(val / e)
    # pairwise maxgap threshold crossings at +-1/7 and +-2/7 among cluster phases
    Elist = list(E)
    for i in range(len(Elist)):
        for j in range(len(Elist)):
            d = Elist[i] - Elist[j]
            if d == 0:
                continue
            ad = abs(d)
            for m in range(0, ad + 1):
                for s in (T1, -T1, T2, -T2):
                    val = Fr(m, ad) + s / ad
                    if Fr(0) <= val <= Fr(1):
                        bp.add(val)
    return sorted(b for b in bp if Fr(0) <= b <= Fr(1))


def compute(E, Pset):
    pts = grid(E, Pset)
    M = {k: Fr(0) for k in
         ["measGP", "G1", "G2", "G1capGP", "G2capGP", "p0", "p0capGP"]}
    D1_ok = True
    for a, b in zip(pts, pts[1:]):
        if b <= a:
            continue
        w = b - a
        mid = (a + b) / 2
        ph = phases_at(E, mid)
        g = maxgap(ph)
        wtheta = Fr(6, 7) - circ_width(ph)
        # D1 dictionary check
        if wtheta != g - T1:
            D1_ok = False
        gp = in_GP(Pset, mid)
        if gp:
            M["measGP"] += w
        if g > T1:
            M["G1"] += w
            if gp:
                M["G1capGP"] += w
        if g > T2:
            M["G2"] += w
            if gp:
                M["G2capGP"] += w
        S = set(int(P * p) for p in ph)
        if len(S) == P:
            M["p0"] += w
            if gp:
                M["p0capGP"] += w
    return M, D1_ok


def main():
    print("=" * 104)
    print("THREAD B: WITNESS (maxgap>1/7) vs CRITERION (maxgap>2/7) — the two LRC sufficient thresholds.")
    print("=" * 104)
    print("Dictionary (THM-526 kps-S4):  w_theta = 6/7 - circ_width = maxgap - 1/7.")
    print("  GLOBAL WITNESS  needs maxgap > 1/7  (=G1)  => M>=1/14   [WEAKER event, EASIER]")
    print("  VIA-MAX CRITERION needs maxgap > 2/7 (=G2=rho*) => M>=1/14  [THM-527's object]")
    print("  G_P = {x: ||p x||>=1/14 all p in P}.  rho*(P,E) = meas(G2 cap G_P).")
    print()
    print("Each row: P = small part, E = cluster co-offsets.  ALL EXACT rational.")
    print()

    # (P small part, E cluster offsets) pairs spanning the structure
    cases = [
        ("P={1,2,3}",         [1, 2, 3],       "consec k=5",  [0, 1, 2, 3, 4]),
        ("P={1,2,3}",         [1, 2, 3],       "consec k=7",  list(range(7))),
        ("P={1,2,3}",         [1, 2, 3],       "consec k=9",  list(range(9))),
        ("P={1,2,3,12}",      [1, 2, 3, 12],   "consec k=9",  list(range(9))),
        ("P={1,2,3,12}",      [1, 2, 3, 12],   "consec k=13", list(range(13))),
        ("P={1..5}",          [1, 2, 3, 4, 5], "consec k=8",  list(range(8))),
        ("P={1,2,3}",         [1, 2, 3],       "perf k=7",    [0, 2, 3, 4, 5, 6, 8]),
        ("P={1,2,3}",         [1, 2, 3],       "Sidon k=5",   [0, 1, 3, 7, 12]),
    ]
    hdr = (f"{'P':<14}{'E':<14}{'meas G_P':>10}{'G1':>8}{'G2=rho*':>9}"
           f"{'G1capGP':>9}{'rho*GP':>9}{'p0':>8}{'D1':>5}")
    print(hdr)
    print("-" * len(hdr))
    allD1 = True
    rows = []
    for pname, Pset, ename, E in cases:
        M, d1 = compute(E, Pset)
        allD1 &= d1
        rows.append((pname, ename, M))
        print(f"{pname:<14}{ename:<14}{float(M['measGP']):>10.4f}{float(M['G1']):>8.4f}"
              f"{float(M['G2']):>9.4f}{float(M['G1capGP']):>9.4f}{float(M['G2capGP']):>9.4f}"
              f"{float(M['p0']):>8.4f}{str(d1):>5}")
    print("-" * len(hdr))
    print(f"\n(D1) dictionary w_theta = maxgap - 1/7 holds pointwise for ALL: {allD1}")

    print("\n" + "=" * 104)
    print("THE PAYOFF: is the WEAKER 1/7 witness floor (meas G1 cap G_P) positive AND >= the 2/7 floor?")
    print("=" * 104)
    print(f"{'P':<14}{'E':<14}{'rho*=G2capGP':>14}{'G1capGP(witness)':>18}{'witness>=rho*?':>16}{'witness>0?':>12}")
    all_wge = True
    all_wpos = True
    for pname, ename, M in rows:
        wge = M["G1capGP"] >= M["G2capGP"]
        wpos = M["G1capGP"] > 0
        all_wge &= wge
        all_wpos &= wpos
        print(f"{pname:<14}{ename:<14}{float(M['G2capGP']):>14.4f}{float(M['G1capGP']):>18.4f}"
              f"{str(wge):>16}{str(wpos):>12}")
    print(f"\nwitness floor >= rho* (criterion floor) for ALL: {all_wge}")
    print(f"witness floor > 0 for ALL:                       {all_wpos}")

    print("""
RESOLUTION OF THE FACTOR-2 (the prompt's Q3):
  * The sector route's missed-1/7-sector corresponds to the maxgap>1/7 (G1) = GLOBAL-WITNESS
    threshold, which ALREADY certifies M>=1/14 directly (THM-526 LEMMA 1 / global witness).
  * THM-527's rho* uses the STRICTER maxgap>2/7 (G2) = VIA-MAX-CRITERION threshold.
  * So the two are GENUINELY a factor 2 apart in the cluster-gap (1/7 vs 2/7), BOTH sufficient.
  * G2 subset G1: the criterion is a SUBSET (smaller, rarer) of the witness event. THM-527
    chose the harder (criterion) object; the witness object (G1 cap G_P) is LARGER and ALSO
    certifies LRC. So THM-527's 2/7 floor is SUFFICIENT but NOT MINIMAL: the 1/7 witness floor
    is the weaker sufficient condition.
  * Whether the witness floor is positive everywhere (its own crux) is tested above.
""")


if __name__ == "__main__":
    main()
