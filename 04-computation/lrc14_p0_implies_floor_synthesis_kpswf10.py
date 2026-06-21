#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_p0_implies_floor_synthesis_kpswf10.py  (kind-pasteur 2026-06-21, THREAD B synthesis)

FINAL THREAD-B SYNTHESIS.  Does the CLOSED sector route p0<=cap_k give EITHER LRC floor
(witness G1capGP or criterion G2capGP=rho*) > 0 ?  i.e. is the sector-route CLOSURE the
same as THM-527's crux, or strictly parallel?

Three exact-rational quantities per config (P, E), with |P|+|E|=13 (admissible):
  p0          = meas{all 7 sectors hit by {frac(e_i x)}}   (the L7 sector object; cap-bounded)
  G1capGP     = witness floor = meas{x in G_P: maxgap{frac(e_i x)} > 1/7}
  G2capGP     = rho* (criterion floor) = meas{x in G_P: maxgap > 2/7}

EXACT CONTAINMENTS (verify): on x in G_P,
   {maxgap>2/7} subset {maxgap>1/7} subset {>=1 sector missed} subset {NOT all hit}.
   So  rho* <= witness <= meas(G_P) - p0capGP <= 1 - p0.

THE DIRECTION TEST (the prompt's Q2):
   p0<=cap  =>  1-p0 >= 1-cap  =>  meas{some sector missed} >= 1-cap.
   But the floors are UPPER-bounded by 1-p0, not lower-bounded.  So p0<=cap does NOT
   force witness>0 or rho*>0.  We CONFIRM by exhibiting (if any) admissible configs with
   p0 <= cap_k yet small/zero floor, OR prove p0<=cap and floor>0 are independent.

THE REVERSE (does a floor>0 follow from p0 being BOUNDED AWAY from 1)?  Test whether
   1 - p0 - witness  (the 'sector-but-not-1/7-gap' slack) is bounded, i.e. whether
   controlling p0 from above ALSO controls the slack.  If witness >= (1-p0) - slack with
   slack small, then p0<=cap with cap<1 would give witness >= 1-cap-slack -- a usable LOWER
   bound.  We measure the slack to see if this transfer is viable.

EXACT rational throughout.
"""
import itertools
from fractions import Fraction as Fr

P7 = 7
T1 = Fr(1, 7)
T2 = Fr(2, 7)
HALF = Fr(1, 14)
CAP = {8: Fr(2243, 5880), 9: Fr(1979, 4004), 10: Fr(55, 91),
       11: Fr(66, 91), 12: Fr(6, 7), 13: Fr(1, 1),
       3: Fr(1), 4: Fr(1), 5: Fr(1), 6: Fr(1), 7: Fr(1)}  # small k: cap=1 (no constraint)


def phases_at(E, x):
    return sorted((int(e) * x) % 1 for e in E)


def maxgap(ph):
    if len(ph) <= 1:
        return Fr(1)
    g = max(b - a for a, b in zip(ph, ph[1:]))
    return max(g, (ph[0] + 1) - ph[-1])


def in_GP(Pset, x):
    return all(min((int(p) * x) % 1, 1 - (int(p) * x) % 1) >= HALF for p in Pset)


def grid(E, Pset):
    bp = {Fr(0), Fr(1)}
    for e in list(E) + list(Pset):
        e = int(e)
        if e == 0:
            continue
        for t in range(0, P7 * e + 1):
            bp.add(Fr(t, P7 * e))
        for k in range(0, e + 1):
            for s in (HALF, -HALF):
                bp.add((Fr(k) + s) / e)
    El = list(E)
    for i in range(len(El)):
        for j in range(len(El)):
            d = El[i] - El[j]
            if d == 0:
                continue
            ad = abs(d)
            for m in range(0, ad + 1):
                for s in (T1, -T1, T2, -T2):
                    v = Fr(m, ad) + s / ad
                    if Fr(0) <= v <= Fr(1):
                        bp.add(v)
    return sorted(b for b in bp if Fr(0) <= b <= Fr(1))


def measures(E, Pset):
    pts = grid(E, Pset)
    p0 = p0capGP = G1capGP = G2capGP = measGP = Fr(0)
    for a, b in zip(pts, pts[1:]):
        if b <= a:
            continue
        w = b - a
        mid = (a + b) / 2
        ph = phases_at(E, mid)
        S = set(int(P7 * p) for p in ph)
        allhit = (len(S) == P7)
        if allhit:
            p0 += w
        gp = in_GP(Pset, mid)
        if gp:
            measGP += w
            if allhit:
                p0capGP += w
            g = maxgap(ph)
            if g > T1:
                G1capGP += w
            if g > T2:
                G2capGP += w
    return dict(p0=p0, p0capGP=p0capGP, measGP=measGP, G1capGP=G1capGP, G2capGP=G2capGP)


def main():
    print("=" * 104)
    print("THREAD B SYNTHESIS: does p0<=cap_k force a positive LRC floor? (sector closure vs THM-527 crux)")
    print("=" * 104)
    print("Admissible consec configs |P|+k=13.  For each k, the worst (min-rho*) P from the prior probe.")
    print()
    worst_P = {
        3: [1, 2, 3, 5, 7, 8, 9, 11, 12, 13],
        4: [1, 2, 3, 5, 7, 8, 9, 11, 13],
        5: [1, 2, 3, 6, 7, 9, 11, 13],
        6: [1, 2, 3, 7, 9, 11, 13],
        7: [1, 2, 3, 9, 11, 13],
        8: [1, 2, 3, 11, 13],
        9: [1, 2, 3, 12],
        10: [1, 2, 3],
        11: [1, 6],
        12: [1],
        13: [],
    }
    hdr = (f"{'k':>3}{'P (worst)':<24}{'cap_k':>9}{'p0':>9}{'p0<=cap?':>10}"
           f"{'rho*':>9}{'witness':>9}{'1-p0':>8}{'slack=1-p0-wit':>16}")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for k in range(3, 14):
        Pset = worst_P[k]
        E = list(range(k))
        m = measures(E, Pset)
        cap = CAP[k]
        p0 = m["p0"]
        slack = (Fr(1) - p0) - m["G1capGP"]
        rows.append((k, Pset, cap, m, slack))
        print(f"{k:>3}{str(Pset):<24}{float(cap):>9.4f}{float(p0):>9.4f}{str(p0 <= cap):>10}"
              f"{float(m['G2capGP']):>9.4f}{float(m['G1capGP']):>9.4f}{float(Fr(1)-p0):>8.4f}"
              f"{float(slack):>16.4f}")
    print("-" * len(hdr))

    print("""
DIRECTION ANALYSIS (the prompt's Q2):
  * p0 <= cap_k holds (the sector route's closure).  But the floors rho*, witness are
    UPPER-bounded by 1-p0 (since {gap>thresh} subset {sector missed}).  p0<=cap gives
    1-p0 >= 1-cap (a LOWER bound on 1-p0), which does NOT lower-bound the floors.
  * => The closed sector route p0<=cap does NOT, by itself, give rho*>0 or witness>0.
    The sector route certifies survival THROUGH p0<=cap by a DIFFERENT mechanism (the
    measS7 union-cover argument), not by exhibiting a gap.  THM-527's crux (a positive gap
    floor) is a SEPARATE statement.  They are PARALLEL.
""")

    # ----- the reverse: is the slack (1-p0 - witness) small enough to transfer? -----
    print("=" * 104)
    print("REVERSE Q: can a p0 UPPER bound give a witness LOWER bound via small slack? (1-p0 - witness)")
    print("=" * 104)
    maxslack = max(s for (_, _, _, _, s) in rows)
    print(f"  max slack (1-p0 - witness) over worst configs = {float(maxslack):.4f}")
    print("  If slack were < 1-cap, then witness >= (1-p0)-slack >= (1-cap)-slack > 0 would TRANSFER.")
    for k, Pset, cap, m, slack in rows:
        oneminuscap = Fr(1) - cap
        transfers = slack < oneminuscap and (Fr(1) - cap) - slack > 0
        print(f"   k={k:>2}: 1-cap={float(oneminuscap):.4f}, slack={float(slack):.4f}, "
              f"(1-cap)-slack={float(oneminuscap-slack):+.4f}  transfer-bound>0: {transfers}")
    print("""
  READING: where (1-cap)-slack > 0, the chain  witness >= (1-p0)-slack >= (1-cap)-slack > 0
  would give a POSITIVE witness floor FROM the sector cap p0<=cap PLUS a slack bound.
  This is the ONLY way the sector route could feed the gap route: not directly, but if one
  ALSO bounds the slack (= the measure of x where a sector is missed but no >1/7 gap forms,
  i.e. the missed sectors are SCATTERED not contiguous).  That slack bound is a NEW lemma,
  not implied by p0<=cap alone.
""")

    print("=" * 104)
    print("BOTTOM LINE (THREAD B):")
    print("=" * 104)
    print("""  p0 (sector, 1/7 cover) and rho* (THM-527, 2/7 gap) are the SAME phi-phase object at
  DIFFERENT thresholds.  Exact dictionary: w_theta = maxgap - 1/7, so
      sector-missed (1/7 cover gap)  <->  maxgap > 1/7  = WITNESS floor (global witness)
      THM-527 criterion              <->  maxgap > 2/7  = rho* (via-max criterion).
  Containment: rho* <= witness <= 1-p0 (all on G_P).  The inequality points the WRONG way
  for the sector CAP (p0<=cap upper-bounds p0, lower-bounds 1-p0, but the floors sit BELOW
  1-p0).  So the closed sector route does NOT close THM-527; they are PARALLEL.
  HOWEVER: the WITNESS floor (maxgap>1/7) is itself LRC-sufficient (global witness, finite-
  Vmax-verified) and is STRICTLY LARGER than rho* (floor 0.0565 vs 1/84).  So THM-527's 2/7
  object is sufficient-but-not-minimal; the easier 1/7 witness floor is the weaker target.""")


if __name__ == "__main__":
    main()
