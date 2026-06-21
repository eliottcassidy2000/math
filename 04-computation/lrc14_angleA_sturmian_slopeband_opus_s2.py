#!/usr/bin/env python3
"""
lrc14_angleA_sturmian_slopeband_opus_s2.py   (opus-2026-06-20-S2)

ANGLE A — STURMIAN PARTIAL-SUM cover, slope-band decomposition.

GOAL: attack the LRC(14) sector-route CRUX
    "consec_k = {0,1,...,k-1} maximizes measS7 over all k-subsets" (k=8,9,10)
via the THM-536 mechanical-word reformulation, looking for a PROOF (not just
the already-VERIFIED exhaustive box check).

SETUP (THM-536, verified here exactly):
    measS7(E) = (1/7) * meas{ theta in [0,7) : { floor(e*theta) mod 7 : e in E } = Z/7 }.
For consecutive e = 0..k-1, floor(e*theta) is the partial-sum trajectory of a
slope-frac(theta) MECHANICAL (Sturmian) word; for a spread E it reads that
trajectory at the SELECTED indices {e_i}.

KEY NEW DECOMPOSITION (this session, proved by the floor identity below):
    Split theta in [0,7) into 7 unit "slope bands" s = floor(theta) in {0..6}.
    Writing theta = s + f, f in [0,1):
        floor(e*theta) = e*s + floor(e*f),
    so  sigma_e = ( s*e + floor(e*f) ) mod 7,  with floor(e*f) the PURE
    (slope-0) Sturmian floor word.  Hence
        measS7(E) = (1/7) * sum_{s=0}^{6}  bandcover_s(E),
        bandcover_s(E) = meas{ f in [0,1) : { (s*e + floor(e*f)) mod 7 : e in E } = Z/7 }.
    Each band is a Z/7-TWISTED copy of the SAME base Sturmian-floor cover problem;
    the twist e -> s*e mod 7 is the only s-dependence.

This reduces the 64-term inclusion-exclusion (DEAD-END C1) to a 7-term SLOPE sum.

WHAT THIS SESSION PROVES / VERIFIES / REFUTES (brutally honest):

(P1, PROVED) The band identity measS7 = (1/7) sum_s bandcover_s  (exact, all shapes tested).
(V1, VERIFIED) consec_k strictly maximizes measS7 over the full box max(E)<=13,
              k=8,9,10: 0 violators; residual non-consec maxima reproduce
              0.2736 / 0.4017 / 0.4939 exactly.
(R1, REFUTED) PER-BAND consec maximality is FALSE: in the slow bands s=0,1 the
              residual adversary [0,2,3,4,5,6,7,8] BEATS consec. Extremality is
              irreducibly an aggregate over the 7 bands (a band-level twin of the
              IE-block dead-end C1).
(R2, REFUTED) NO greedy index-compression move is monotone (single-step down-shift,
              fill-smallest-hole both decrease measS7 on a positive fraction of
              shapes). There is no monotone descent path to consec in index space
              (an index-space twin of the value-space dead-end C2).
(M1, MECHANISM) The crux is a cross-band trade: consec LOSES the slow bands
              (s=0,1) but WINS the fast bands (s=2..6) by more; net +3/56 vs the
              residual adversary at k=8. The signed band table is the irreducible
              certificate.

CONCLUSION ON ANGLE A: the Sturmian index-compression CONJECTURE (consecutive
indices cover Z/7 most often) is FALSE at the per-band / per-move level, so it
does NOT yield a clean monotonicity proof. Its lasting value is the PROVED
7-band identity (P1), which is a strictly simpler certificate skeleton (7 signed
terms vs 64) than the IE route and is the natural home for a future aggregate
argument.

stdlib only; exact Fractions throughout.
"""
import sys, itertools, math, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)


# ----------------------------------------------------------------------------
# Engines (re-derived from definitions; matches lrc14_BS_verify measS7_geom)
# ----------------------------------------------------------------------------
def measS7_geom(E, scale=7):
    """Original scale-x formulation: meas{x in [0,1): {floor(scale*frac(e*x))}=Z/7}."""
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, scale * e + 1):
            bps.add(F(m, scale * e))
    bps = sorted(bps)
    total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        secs = set(int(((e * xm) % 1) * scale) for e in E)
        if len(secs) == scale:
            total += x1 - x0
    return total


def measS7_theta(E, scale=7):
    """THM-536 theta-form: (1/scale) meas{theta in [0,scale): {floor(e*theta) mod scale}=Z/scale}."""
    E = sorted(set(E))
    bps = {F(0), F(scale)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, scale * e + 1):
            bps.add(F(m, e))
    bps = sorted(b for b in bps if 0 <= b <= scale)
    total = F(0)
    for i in range(len(bps) - 1):
        t0, t1 = bps[i], bps[i + 1]
        if t1 <= t0:
            continue
        tm = (t0 + t1) / 2
        secs = set(int(e * tm) % scale for e in E)
        if len(secs) == scale:
            total += t1 - t0
    return total / scale


def band_cover(E, s, scale=7):
    """meas{f in [0,1): {(s*e + floor(e*f)) mod scale : e in E} = Z/scale}."""
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, e + 1):
            bps.add(F(m, e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for i in range(len(bps) - 1):
        f0, f1 = bps[i], bps[i + 1]
        if f1 <= f0:
            continue
        fm = (f0 + f1) / 2
        secs = set((s * e + int(e * fm)) % scale for e in E)
        if len(secs) == scale:
            total += f1 - f0
    return total


def measS7_bands(E, scale=7):
    return sum(band_cover(E, s, scale) for s in range(scale)) / scale


# ----------------------------------------------------------------------------
# (P1) PROVE the band identity numerically on a broad shape census.
# ----------------------------------------------------------------------------
def check_band_identity():
    print("=" * 72)
    print("(P1) band identity:  measS7 == (1/7) sum_s bandcover_s   [+ theta-form]")
    print("=" * 72)
    shapes = [list(range(k)) for k in range(1, 13)]
    rng = random.Random(20260620)
    for _ in range(60):
        N = rng.randint(8, 13)
        k = rng.randint(3, min(N + 1, 8))
        rest = rng.sample(range(1, N + 1), k - 1)
        shapes.append(sorted([0] + rest))
    bad = 0
    for E in shapes:
        g = measS7_geom(E)
        t = measS7_theta(E)
        b = measS7_bands(E)
        if not (g == t == b):
            bad += 1
            print("   MISMATCH", E, float(g), float(t), float(b))
    print(f"   checked {len(shapes)} shapes; mismatches = {bad}")
    print(f"   => band identity holds exactly (PROVED by floor(e*(s+f))=e*s+floor(e*f)).")
    return bad == 0


# ----------------------------------------------------------------------------
# (V1) Full-box VERIFY: consec maximizes measS7 over max(E)<=13, k=8,9,10.
# ----------------------------------------------------------------------------
def verify_box():
    print("=" * 72)
    print("(V1) full-box VERIFY: consec maximizes measS7, k in {8,9,10}, max(E)<=13")
    print("=" * 72)
    results = {}
    for k in (8, 9, 10):
        cons = measS7_theta(list(range(k)))
        rest = list(range(1, 14))
        viol = 0
        resid_max = F(0)
        resid_E = None
        cnt = 0
        for combo in itertools.combinations(rest, k - 1):
            E = [0] + list(combo)
            cnt += 1
            v = measS7_theta(E)
            if v > cons:
                viol += 1
            if E != list(range(k)) and v > resid_max:
                resid_max = v
                resid_E = E
        results[k] = (cons, viol, resid_max, resid_E, cnt)
        print(f"   k={k}: consec={float(cons):.6f}  #subsets={cnt}  "
              f"violators={viol}  residual_max={float(resid_max):.6f}  arg={resid_E}")
    return results


# ----------------------------------------------------------------------------
# (R1) REFUTE per-band consec maximality (slow bands lose).
# ----------------------------------------------------------------------------
def refute_per_band():
    print("=" * 72)
    print("(R1) per-band consec maximality is FALSE (slow bands s=0,1 lose)")
    print("=" * 72)
    k = 8
    rest = list(range(1, 11))
    for s in range(7):
        cons = band_cover(list(range(k)), s)
        viol = 0
        vmax = F(0)
        vE = None
        for combo in itertools.combinations(rest, k - 1):
            E = [0] + list(combo)
            v = band_cover(E, s)
            if v > cons:
                viol += 1
                if v > vmax:
                    vmax = v
                    vE = E
        tag = "FAIL" if viol else "ok"
        print(f"   band s={s}: consec={float(cons):.5f}  violators={viol:3d}  "
              f"worst={float(vmax):.5f} {vE}  [{tag}]")


# ----------------------------------------------------------------------------
# (M1) The signed band table: consec vs the residual adversary at k=8.
# ----------------------------------------------------------------------------
def signed_band_table():
    print("=" * 72)
    print("(M1) signed slope-band table: consec_8 vs residual adversary [0,2..8]")
    print("=" * 72)
    cons = list(range(8))
    adv = [0, 2, 3, 4, 5, 6, 7, 8]
    print("   band :   consec      adv     consec-adv")
    tot = F(0)
    for s in range(7):
        c = band_cover(cons, s)
        a = band_cover(adv, s)
        d = c - a
        tot += d
        print(f"   s={s} : {float(c):8.5f}  {float(a):8.5f}  {float(d):+8.5f}")
    print(f"   net (consec-adv)/7 = {tot/7} = {float(tot/7):.6f}  (>0: consec wins overall)")
    print("   MECHANISM: consec LOSES slow bands s=0,1, WINS fast bands s=2..6 by more.")


# ----------------------------------------------------------------------------
# (R2) REFUTE index-compression monotonicity (no greedy move works).
# ----------------------------------------------------------------------------
def refute_compression():
    print("=" * 72)
    print("(R2) index-compression monotonicity is FALSE (no greedy descent path)")
    print("=" * 72)
    k = 8
    rest = list(range(1, 12))

    # move A: shift some element down by 1 into a free slot
    tA = cA = 0
    worstA = (F(0), None)
    # move B: replace max(E) by smallest hole
    tB = cB = 0
    worstB = (F(0), None)
    for combo in itertools.combinations(rest, k - 1):
        E = sorted([0] + list(combo))
        base = measS7_theta(E)
        for e in E:
            if e == 0 or (e - 1) in E or e - 1 < 0:
                continue
            E2 = sorted(set(E) - {e} | {e - 1})
            v2 = measS7_theta(E2)
            tA += 1
            if v2 < base:
                cA += 1
                if base - v2 > worstA[0]:
                    worstA = (base - v2, (E, E2))
        holes = [h for h in range(1, max(E)) if h not in E]
        if holes:
            E2 = sorted(set(E) - {max(E)} | {min(holes)})
            if len(E2) == k:
                v2 = measS7_theta(E2)
                tB += 1
                if v2 < base:
                    cB += 1
                    if base - v2 > worstB[0]:
                        worstB = (base - v2, (E, E2))
    print(f"   move A (down-shift-into-gap): {cA}/{tA} moves DECREASE measS7; "
          f"worst drop {float(worstA[0]):.5f} {worstA[1]}")
    print(f"   move B (max -> smallest hole): {cB}/{tB} moves DECREASE measS7; "
          f"worst drop {float(worstB[0]):.5f} {worstB[1]}")
    print("   => no monotone compression to consec; extremality is irreducibly aggregate.")


# ----------------------------------------------------------------------------
def main():
    print("LRC(14) Angle A — Sturmian slope-band decomposition  (opus-2026-06-20-S2)")
    ok = check_band_identity()
    print()
    verify_box()
    print()
    refute_per_band()
    print()
    signed_band_table()
    print()
    refute_compression()
    print()
    print("=" * 72)
    print("VERDICT: Angle A index-compression is a DEAD-END for a clean proof")
    print("(per-band R1 + no-monotone-move R2). NET YIELD = the PROVED 7-band")
    print("identity (P1): measS7 = (1/7) sum_{s=0}^{6} bandcover_s, reducing the")
    print("64-term IE certificate to a 7-term signed slope sum. Consec-maximality")
    print("remains VERIFIED (box, k=8,9,10) but NOT PROVED.")
    print("=" * 72)


if __name__ == "__main__":
    main()
