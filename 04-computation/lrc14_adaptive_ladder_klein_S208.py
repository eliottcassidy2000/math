#!/usr/bin/env python3
"""
lrc14_adaptive_ladder_klein_S208.py

The ADAPTIVE-SPLIT LADDER for the LRC(14) covering case (klein-S208):
verification of the three legs and the end-to-end composite.

THE LADDER. Given a covering 13-set S with the THM-527 split S = P ∪ L
(P = S∩[1,13], L the cluster, Vmax = max L), re-split ADAPTIVELY: move the
slow cluster members u < (9/14)·Vmax into the small part:
    P~ = P ∪ {u ∈ L : u < (9/14)Vmax},   L~ = {u ∈ L : u ≥ (9/14)Vmax},
    k~ = |L~|, |P~| = 13 − k~,  E~ = {Vmax − u : u ∈ L~} ⊆ [0, (5/14)Vmax],
so spread(E~) ≤ (5/14)Vmax, i.e. split-ratio r~ = Vmax/spread(E~) ≥ 14/5 = 2.8:
the residual cluster is ALWAYS in the a-priori (aliasing+embed) zone.

MEASURE-FLOOR LEGS (this script verifies both numerically/exactly):
 (leg A, k~ ≤ 7)  μ_{1/7}(E~) = 1 a.e. (THM-530 pigeonhole), so
    ρ~* ≥ meas(G_P~).  And meas(G_P~) > 0 ALWAYS, by LRC(≤13) + LIPSCHITZ
    FATTENING: |P~| ≤ 12 distinct speeds + observer = an LRC(≤13) instance,
    so ∃τ0 with min_p ||p τ0|| ≥ 1/13 (SETTLED, owner directive 2026-07-02);
    each ||p·|| is p-Lipschitz, so the ball |τ−τ0| ≤ (1/13 − 1/14)/max(P~)
    stays ≥ 1/14:  meas(G_P~) ≥ 2·(1/182)/max(P~) = 1/(91·max P~) > 0.
    ** The apex-7 / Fraenkel-tiling wall CANNOT occur for P~: exact tiling
    of the danger sets would need M(P~) < 1/13, contradicting LRC(≤13). **
 (leg B, k~ ≥ 8)  |P~| = 13 − k~ ≤ 5 and Bonferroni
    meas(G_P~) ≥ 1 − |P~|/7, with the THM-661/663 shape-universal floors
    μ_{1/7}(E~) ≥ floor_{k~}: ρ~* ≥ (1 − |P~|/7) + floor_{k~} − 1 > 0,
    margins: k~=8: 2/7−0.2389=+0.047; 9: 3/7−0.3554=+0.073; 10: +0.125;
             11: +0.119; 12: +0.213; 13: +0.309.

REALIZATION (the remaining node, localized): the snap/embed realizes the
witness from a good x for (P~, E~) only if the P~-members can hold safety at
the realized τ. Members p ≤ 13 ride any τ-move (slack p/V tiny). Members
t ∈ (13, V/14] ride the SNAP INTERVAL (width ≤ (6/7)/V, phase move ≤ 6t/(7V)
≤ 3/49). Members in the MID-BAND (V/14, 9V/14) can neither ride nor cluster:
*** the entire remaining difficulty of LRC(14) is the mid-band × realization ***

This script:
 [1] verifies leg A exactly on slow-heavy census instances: exact arc-union
     meas(G_P~) vs the Lipschitz floor 1/(91·max P~) (both positive);
 [2] verifies leg B arithmetic against the THM-663 floors;
 [3] builds MID-BAND-FREE covering instances and runs the END-TO-END
     composite chain: exact meas floor > 0, aliasing existence on E~,
     Lean-faithful drift-embed with P~-safety at the realized τ, and exact
     ground truth M(S) ≥ 1/14;
 [4] measures how the mid-band population statistics look over random
     covering sets (how big is the residual family?).
"""

from fractions import Fraction
from math import gcd
import random

random.seed(20260709)

QS = list(range(2, 15))
F = Fraction


def is_covering(S):
    return all(any(s % q == 0 for s in S) for q in QS)


# ---------------------------------------------------------------- exact G_T measure
def G_measure_exact(T, c=F(1, 14), slack=None):
    """exact meas{x in [0,1): forall t in T, ||t x|| >= c + slack_t} via
    arc-union. slack: dict t->Fraction or None."""
    arcs = []
    for t in T:
        st = (slack or {}).get(t, F(0))
        w = c + st
        if w >= F(1, 2):
            return F(0)
        for m in range(t + 1):
            lo, hi = F(m, t) - w / t, F(m, t) + w / t
            arcs.append((max(lo, F(0)), min(hi, F(1))))
    arcs.sort()
    covered = F(0)
    cur_lo, cur_hi = None, None
    for (a, b) in arcs:
        if b <= a:
            continue
        if cur_lo is None:
            cur_lo, cur_hi = a, b
        elif a <= cur_hi:
            cur_hi = max(cur_hi, b)
        else:
            covered += cur_hi - cur_lo
            cur_lo, cur_hi = a, b
    if cur_lo is not None:
        covered += cur_hi - cur_lo
    return F(1) - covered


# ---------------------------------------------------------------- M(S) exact
def dist_num(v, num, den):
    r = (v * num) % den
    return F(min(r, den - r), den)


def M_exact(S):
    best = F(0)
    dens = set()
    for i in range(len(S)):
        dens.add(2 * S[i])
        for j in range(i, len(S)):
            dens.add(S[i] + S[j])
    for D in sorted(dens):
        for m in range(1, D // 2 + 1):
            if gcd(m, D) != 1:
                continue
            val = min(dist_num(v, m, D) for v in S)
            if val > best:
                best = val
    return best


# ---------------------------------------------------------------- W, intW, TV (exact)
THETA = F(1, 7)


def exact_intW_and_TV(E):
    dens = set()
    E = sorted(E)
    for i in range(len(E)):
        if E[i] > 0:
            dens.add(E[i])
        for j in range(i + 1, len(E)):
            if E[j] - E[i] > 0:
                dens.add(E[j] - E[i])
    pts = set()
    for d in dens:
        for m in range(d):
            pts.add(F(m, d))
    bps = sorted(pts)
    n = len(bps)
    intW = F(0)
    slopes_seq = []
    for idx in range(n):
        a = bps[idx]
        b = bps[idx + 1] if idx + 1 < n else F(1)
        if b <= a:
            continue
        mid = (a + b) / 2
        ph = sorted((F(e) * mid - int(F(e) * mid), e) for e in E)
        k = len(ph)
        gmid = [(ph[t + 1][0] - ph[t][0], ph[t + 1][1] - ph[t][1]) for t in range(k - 1)]
        gmid.append((1 - ph[-1][0] + ph[0][0], ph[0][1] - ph[-1][1]))
        kinks = {a, b}
        for (gm, s) in gmid:
            if s != 0:
                xc = mid + (THETA - gm) / s
                if a < xc < b:
                    kinks.add(xc)
        pp = sorted(kinks)
        for t in range(len(pp) - 1):
            u, v = pp[t], pp[t + 1]
            m2 = (u + v) / 2
            Wm, Sm = F(0), 0
            for (gm, s) in gmid:
                gv = gm + s * (m2 - mid)
                if gv > THETA:
                    Wm += gv - THETA
                    Sm += s
            intW += Wm * (v - u)
            slopes_seq.append(Sm)
    tv = 0
    m = len(slopes_seq)
    for i in range(m):
        tv += abs(slopes_seq[(i + 1) % m] - slopes_seq[i])
    return intW, tv


# ---------------------------------------------------------------- drift embed (Lean-faithful)
def embed_with_Ptilde(E, Pt, V, spread):
    """scan grid points; at each, all circular gaps; fire iff the EXACT
    minReach_ge_of_driftGap hypothesis holds AND P~ safe at realized tau."""
    for j in range(1, V):
        x = F(j, V)
        ok = True
        for p in Pt:
            d = F(p) * x
            d -= int(d)
            d = min(d, 1 - d)
            if d < F(1, 14) + F(p, V):   # ride slack for the full period
                ok = False
                break
        if not ok:
            continue
        ph = sorted(F(e) * x - int(F(e) * x) for e in E)
        k = len(ph)
        cands = [(ph[t], ph[t + 1] - ph[t]) for t in range(k - 1)]
        cands.append((ph[-1], 1 - ph[-1] + ph[0]))
        for (a, g) in cands:
            if F(1, 7) + 2 * (F(spread) * (a + g / 2) / V) < g:
                tau = (F(j) + a + g / 2) / V
                safe = True
                for p in Pt:
                    d = F(p) * tau
                    d -= int(d)
                    d = min(d, 1 - d)
                    if d < F(1, 14):
                        safe = False
                        break
                if safe:
                    return j, a, g, tau
    return None


FLOORS = {8: F(761132, 10 ** 6), 9: F(644603, 10 ** 6), 10: F(553111, 10 ** 6),
          11: F(404751, 10 ** 6), 12: F(355876, 10 ** 6), 13: F(308844, 10 ** 6)}


def main():
    print("=" * 78)
    print("ADAPTIVE-SPLIT LADDER verification (klein-S208)")
    print("=" * 78)

    # ---- [1] leg A: LRC(<=13) Lipschitz fattening vs exact meas(G_P~) ----
    print("\n[1] Leg A (k~<=7): exact meas(G_P~) vs Lipschitz floor 1/(91 max P~)")
    print("    (P~ = small part + slow cluster members; |P~| up to 12)")
    legA_cases = [
        # (P, slow L-members) -- slow-heavy shapes from the census families
        ((8, 9, 10, 12), (22, 26, 35, 49, 60)),          # |P~|=9
        ((7, 9, 10, 11, 12), (22, 26, 39, 52, 65)),      # |P~|=10, structured 13k
        ((11, 12, 13), (14, 22, 26, 30, 42, 55)),        # |P~|=9
        ((8, 9, 10, 12), (22, 26, 33, 39, 45, 51, 57)),  # |P~|=11
        ((1, 2, 3, 4, 5), (22, 26, 33, 39, 45, 51, 57)), # |P~|=12 worst-ish
    ]
    for (P, slow) in legA_cases:
        Pt = sorted(set(P) | set(slow))
        m = G_measure_exact(Pt)
        lip = F(1, 91 * max(Pt))
        print(f"   P~={Pt}")
        print(f"     exact meas(G_P~) = {m} = {float(m):.6f}   "
              f"Lipschitz floor = {lip} = {float(lip):.6f}   "
              f"positive: {m > 0}   floor<=meas: {lip <= m}")

    # ---- [2] leg B arithmetic --------------------------------------------
    print("\n[2] Leg B (k~>=8): Bonferroni(|P~|=13-k~) + shape-universal floors")
    for kt in range(8, 14):
        pt = 13 - kt
        bonf = 1 - F(pt, 7)
        need = 1 - FLOORS[kt]
        print(f"   k~={kt}: |P~|<={pt}, meas(G_P~)>={bonf}={float(bonf):.4f}, "
              f"need >{float(need):.4f}, margin {float(bonf - need):+.4f}  "
              f"{'OK' if bonf > need else 'FAIL'}")

    # ---- [3] end-to-end composite on MID-BAND-FREE instances -------------
    print("\n[3] End-to-end composite on mid-band-free covering instances")
    print("    zones: P=[1,13], T=(13,V/14], MID=(V/14,9V/14) EMPTY, cluster=[9V/14,V]")
    composites = []
    # V=280: T-zone=(13,20], cluster zone=[180,280]
    #   P={8,9,10,12} covers q in {2,3,4,5,6,8,9,10,12}; missed {7,11,13,14}
    #   cluster covering duty: 7|252, 11|220, 13|260, 14|252 -- all in [180,280]
    composites.append(((8, 9, 10, 12), (14,), (252, 220, 260, 189, 200, 230, 270, 280)))
    #   with T-member 14 (covers 7,14 too); k~=8, |P~|=5
    # V=420: T=(13,30], cluster=[270,420]
    composites.append(((8, 9, 10, 11), (26,), (280, 308, 300, 330, 350, 390, 405, 420)))
    #   P covers {2,3,4,5,8,9,10,11,6? 6|no: 8,9,10,11 -> 2,4,8|3,9|2,5,10|11; 6 missed!
    #   6|300? yes 300=6*50 in cluster; 7|280; 13|26 (T); 12|300? 300=12*25 yes; 14|280.
    for (P, T, Lc) in composites:
        V = max(Lc)
        S = sorted(set(P) | set(T) | set(Lc))
        cov = is_covering(S)
        if len(S) != 13:
            print(f"   SKIP (|S|={len(S)}): {S}")
            continue
        Ltilde = sorted(Lc)
        E = sorted(V - u for u in Ltilde)
        spread = E[-1]
        kt = len(Ltilde)
        Pt = sorted(set(P) | set(T))
        # measure leg
        if kt <= 7:
            gm = G_measure_exact(Pt)
            floor_txt = f"legA meas(G_P~)={float(gm):.5f}>0"
            floor_ok = gm > 0
        else:
            pt = len(Pt)
            gm = G_measure_exact(Pt)
            need = 1 - FLOORS[kt]
            floor_ok = gm > need
            floor_txt = f"legB meas(G_P~)={float(gm):.5f} > 1-floor_{kt}={float(need):.5f}: {floor_ok}"
        # aliasing existence on E~
        intW, tv = exact_intW_and_TV(E)
        v0sq = F(tv, 12) / intW if intW > 0 else None
        alias = intW > 0 and F(V) ** 2 > v0sq
        # embed with P~ safety
        emb = embed_with_Ptilde(E, Pt, V, spread)
        # ground truth
        M = M_exact(S) if V <= 460 else None
        print(f"   S={S} (covering={cov})")
        print(f"     V={V} r~={V / spread:.2f} k~={kt} |P~|={len(Pt)}")
        print(f"     measure: {floor_txt}")
        print(f"     aliasing: intW={float(intW):.5f} TV={tv} V0={float(v0sq) ** 0.5:.1f} fires={alias}")
        print(f"     embed(+P~ at tau): {'FIRES j=%d gap=%.4f' % (emb[0], float(emb[2])) if emb else 'no'}")
        if M is not None:
            print(f"     ground truth M(S) = {M} = {float(M):.5f} >= 1/14: {M >= F(1, 14)}")

    # ---- [4] mid-band population over random covering sets ----------------
    print("\n[4] Mid-band population over random k>=8 covering sets (how big is the residual?)")
    tot, empty_mid = 0, 0
    strata = {}
    for _ in range(4000):
        V = random.choice([120, 200, 300, 500, 800])
        P = random.sample(range(1, 14), random.randint(3, 5))
        L = {V}
        missed = [q for q in QS if not any(p % q == 0 for p in P)]
        okgen = True
        for q in missed:
            if any(u % q == 0 for u in L):
                continue
            lo = (14 + q - 1) // q
            hi = V // q
            if lo > hi:
                okgen = False
                break
            L.add(q * random.randint(lo, hi))
        if not okgen:
            continue
        while len(L) < 13 - len(P):
            L.add(random.randint(14, V))
        S = sorted(set(P) | L)
        if len(S) != 13 or not is_covering(S):
            continue
        k = len(L)
        if k < 8:
            continue
        tot += 1
        mid = [u for u in L if F(V, 14) < u < F(9, 14) * V]
        nm = len(mid)
        strata[nm] = strata.get(nm, 0) + 1
        if nm == 0:
            empty_mid += 1
    print(f"   random k>=8 covering instances: {tot}")
    for nm in sorted(strata):
        print(f"     #mid-band members = {nm}: {strata[nm]} ({100 * strata[nm] / tot:.1f}%)")
    print(f"   mid-band-FREE (composite-closable): {empty_mid} ({100 * empty_mid / tot:.1f}%)")
    print("   NOTE: random generation is uniform-ish; the MEASURE of the residual family")
    print("   is not the point -- the point is it is nonempty and characterized:")
    print("   *** the open content of LRC(14) = covering sets with a speed in (V/14, 9V/14)")
    print("       x the realization (Kronecker) node; everything else is closed. ***")

    print("\nDONE.")


if __name__ == '__main__':
    main()
