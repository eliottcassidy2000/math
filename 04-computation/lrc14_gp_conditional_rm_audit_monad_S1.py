#!/usr/bin/env python3
"""
monad-explorer-2026-07-07-S1 -- G_P-conditional reverse-Markov audit of the E[maxgap] program.

CONTEXT (synthesizing kps-S57/S58, opus-S133, death-star-S1, THM-530, LRCFourteenSkeleton):
The fleet reduced "Route-1 density floor" to inf_E E[maxgap] > 1/7 by reverse-Markov and
celebrates a ~+0.06 margin.  This script audits three seams that the reduction does NOT cover:

 SEAM 1 (quantitative threshold).  The skeleton's hlarge obligation is G2 >= m_P
   (m_P = 14249/252252, THM-530), not G2 > 0.  Via reverse-Markov the true k=13 target is
      inf_E E[maxgap] >= 1/7 + (6/7)*m_P  =  T*  (exact rational, ~0.1913),
   so against death-star's exact minima (~0.201) the margin is ~0.010, NOT 0.06.

 SEAM 2 (G_P intersection, k=8..12).  G2(P,E) = meas(G_P /\ {maxgap_E > 1/7}).  THM-530's
   k>=8 route is the UNION bound rho* >= meas(G_P) + mu_{1/7}(E) - 1, which needs
   mu_{1/7}(E) > 1 - min_{|P|=13-k} meas(G_P)  (~0.62 at k=8 down to 1/7 at k=12).
   Reverse-Markov can never deliver mu beyond (7/6)(E[maxgap]-1/7) ~ 0.07-0.18: the
   E[maxgap] reduction is STRUCTURALLY VACUOUS for k=8..12 under the union bound.
   It only serves the k=13 / P=empty leg.  (The two "density floors" -- four-leg leg-3
   vs skeleton hlarge -- were conflated across architectures.)

 SEAM 3 (proposed fix, NEW object).  G_P-CONDITIONAL reverse-Markov:
      rho*_{1/7}(P,E) >= meas(G_P) * (7/6) * ( E[maxgap | x in G_P] - 1/7 ),
   so the corrected reduced target for k=8..12 is the CONDITIONAL mean
      E[maxgap_E | G_P] >= 1/7 + (6/7) * m_P / meas(G_P).
   This restores a mean-route for ALL k iff the conditional means clear the bar --
   the threat is exactly THM-530's anti-correlation pathology (the 2/7 zeros).
   Nobody has computed this object; we compute it adversarially here.

 Plus: crux-class-constrained descent (saturated+primitive+single-scale) on E[maxgap],
 exact-verified with death-star's corrected integrator, against the SEAM-1 threshold T*.

Tournament Analysis declaration:
  vertices: proof-obligation carriers (skeleton nodes hsmall/hlarge/hpartA and their
            (P,E)-shape instances), not runners;
  pairwise observable: does leg A's floor survive leg B's consumption threshold
            (positivity vs m_P vs union-bound requirement);
  switch/gauge: orient toward the leg with the smaller certified margin;
  tie Hamiltonian path: exact rational ledger -> conditional-mean probe -> constrained
            descent -> arc-count growth shadow.
"""
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import random

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# ---------- import death-star's exact integrator (cite, don't re-derive) ----------
_spec = spec_from_file_location(
    "ds", ROOT / "04-computation" / "lrc_maxgap_ap_minimality_check_deathstar_S1.py")
_ds = module_from_spec(_spec); _spec.loader.exec_module(_ds)
Emaxgap_exact = _ds.Emaxgap_exact

M_P = F(14249, 252252)          # THM-530 admissible small-part floor
THR = F(1, 7)

# ---------- Part 0: exact interval machinery for G_P ----------
def gp_intervals(P):
    """Exact maximal intervals of G_P = {x in [0,1): ||p x|| >= 1/14 for all p in P}."""
    bad = []
    for p in P:
        w = F(1, 14 * p)
        for j in range(p + 1):
            lo, hi = F(j, p) - w, F(j, p) + w
            bad.append((max(lo, F(0)), min(hi, F(1))))
    bad.sort()
    merged = []
    for lo, hi in bad:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    good, cur = [], F(0)
    for lo, hi in merged:
        if lo > cur:
            good.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        good.append((cur, F(1)))
    return good

def gp_measure(P):
    return sum(hi - lo for lo, hi in gp_intervals(P))

def min_gp_over_size(psz, pool=range(1, 14)):
    best, bestP = F(2), None
    for P in combinations(pool, psz):
        m = gp_measure(P)
        if m < best:
            best, bestP = m, P
    return best, bestP

# ---------- numeric maxgap over a grid restricted to G_P ----------
def cond_stats(P, E, n=24000, thr=1.0 / 7.0):
    """Returns (meas(G_P) exact, E[mg|G_P], P(mg>thr | G_P), min mg on G_P) numerically.
    E are co-offsets, tooth at 0 included iff 0 in E (skeleton convention: 0 in E)."""
    iv = gp_intervals(P) if P else [(F(0), F(1))]
    meas = sum(hi - lo for lo, hi in iv)
    xs = []
    for lo, hi in iv:
        m = max(8, int(round(n * float(hi - lo))))
        xs.append(np.linspace(float(lo), float(hi), m + 2)[1:-1])
    xs = np.concatenate(xs)
    ph = np.mod(np.outer(xs, np.array(E, dtype=np.float64)), 1.0)
    ph.sort(axis=1)
    gaps = np.empty_like(ph)
    gaps[:, :-1] = np.diff(ph, axis=1)
    gaps[:, -1] = ph[:, 0] + 1.0 - ph[:, -1]
    mg = gaps.max(axis=1)
    return meas, float(mg.mean()), float((mg > thr).mean()), float(mg.min())

def plain_mu_and_mean(E, n=24000, thr=1.0 / 7.0):
    _, mean, mu, _ = cond_stats((), E, n, thr)
    return mu, mean

# ---------- family predicates (crux class) ----------
def saturated(v):
    return all(any(x % q == 0 for x in v) for q in range(2, 15))

def primitive(v):
    from math import gcd
    g = 0
    for x in v:
        g = gcd(g, x)
    return g == 1

def single_scale(v):
    return max(v) <= 13 * min(v)

def in_crux(v):
    return saturated(v) and primitive(v) and single_scale(v)

def report(s):
    print(s, flush=True)

if __name__ == "__main__":
    random.seed(20260707)

    report("=" * 78)
    report("PART 0 -- verify THM-530 exact G_P machinery (audit trail)")
    report("=" * 78)
    per_size = {}
    for psz in range(0, 6):           # |P| = 13-k for k=8..13
        if psz == 0:
            per_size[0] = (F(1), ())
            report("  |P|=0: meas=1")
            continue
        m, Pmin = min_gp_over_size(psz)
        per_size[psz] = (m, Pmin)
        report(f"  |P|={psz}: min meas(G_P) = {m} = {float(m):.5f} at P={Pmin}")
    # spot-verify m_P at |P|=10 (THM-530 claims 14249/252252 at {1,2,3,5,7,8,9,11,12,13})
    mp_claim = gp_measure((1, 2, 3, 5, 7, 8, 9, 11, 12, 13))
    report(f"  THM-530 m_P shape {(1,2,3,5,7,8,9,11,12,13)}: meas = {mp_claim} "
           f"(claimed 14249/252252: {'MATCH' if mp_claim == M_P else 'MISMATCH'})")

    report("")
    report("=" * 78)
    report("PART 1 -- SEAM 1: the quantitative k=13 threshold T* (exact)")
    report("=" * 78)
    Tstar = THR + F(6, 7) * M_P
    report(f"  hlarge needs G2 >= m_P = {M_P} ~ {float(M_P):.5f}")
    report(f"  reverse-Markov k=13 target: E[maxgap] >= T* = 1/7 + (6/7)m_P = "
           f"{Tstar} ~ {float(Tstar):.6f}")
    for nm, val in [("AP 93/440", F(93, 440)),
                    ("death-star prim-sat 2*{1..12}+13", F(145091, 720720)),
                    ("S57 adversarial (exact, death-star)", F(168192271811, 829905476400))]:
        report(f"    {nm}: E[mg] = {float(val):.6f}, margin over T* = {float(val - Tstar):+.6f}"
               f"  (over 1/7: {float(val - THR):+.6f})")

    report("")
    report("=" * 78)
    report("PART 2 -- SEAM 2: union-bound requirement vs reverse-Markov ceiling, per k")
    report("=" * 78)
    report("  k  |P|  min meas(G_P)      union needs mu >   (+mP quant)   RM ceiling*")
    for k in range(8, 14):
        psz = 13 - k
        m = per_size[psz][0]
        need_pos = 1 - m
        need_qnt = 1 - m + M_P
        # RM ceiling: even a generous E[maxgap]<=0.30 for k-point families
        rm_ceil = F(7, 6) * (F(30, 100) - THR)
        verdict = "VACUOUS" if rm_ceil < need_pos else "ok"
        report(f"  {k:2d}  {psz}   {float(m):.5f}            {float(need_pos):.5f}"
               f"            {float(need_qnt):.5f}        {float(rm_ceil):.5f}  {verdict}")
    report("  (*) RM ceiling = (7/6)(E[mg]-1/7) at a generous E[mg]=0.30; actual minimizing")
    report("      families have E[mg]~0.20 => ceiling ~0.07.  => k=8..11 unreachable by RM;")
    report("      k=12 needs E[mg|.]>0.2653 -- above every observed minimizing family.")

    report("")
    report("=" * 78)
    report("PART 3 -- SEAM 3: the NEW G_P-conditional mean E[maxgap | G_P] (k=8..12)")
    report("=" * 78)
    report("  conditional-RM: rho* >= meas(G_P)*(7/6)*(E[mg|G_P]-1/7); need >= m_P")
    report("  i.e. E[mg|G_P] >= 1/7 + (6/7)*m_P/meas(G_P)  [bar depends on shape]")
    report("")
    # (a) THM-530's 2/7-pathology shapes (the anti-correlation stress test), at 1/7:
    pathology = [
        ((1, 2, 3, 12, 13), tuple([0] + list(range(2, 9)))),          # k=8
        ((1, 2, 3, 13),     tuple([0] + list(range(2, 10)))),         # k=9
        ((1, 2, 3),         tuple([0] + list(range(2, 11)))),         # k=10
    ]
    report("  (a) THM-530 via-max-zero pathology shapes (worst known anti-correlation):")
    for P, E in pathology:
        meas, cmean, cmu, cmin = cond_stats(P, E)
        bar = float(THR + F(6, 7) * M_P / meas)
        rho_lb = float(meas) * (7.0 / 6) * (cmean - 1.0 / 7)
        report(f"    P={P} E={E}: meas(G_P)={float(meas):.4f}  E[mg|G_P]={cmean:.4f}  "
               f"bar={bar:.4f}  margin={cmean - bar:+.4f}  condRM rho*>={rho_lb:.4f} "
               f"{'>=m_P OK' if rho_lb >= float(M_P) else '< m_P FAIL'}")

    # (b) worst-P consecutive shapes:
    report("  (b) consecutive E with the per-size WORST P (min meas(G_P)):")
    for k in range(8, 13):
        psz = 13 - k
        P = per_size[psz][1]
        E = tuple(range(k))
        meas, cmean, cmu, cmin = cond_stats(P, E)
        bar = float(THR + F(6, 7) * M_P / meas)
        rho_lb = float(meas) * (7.0 / 6) * (cmean - 1.0 / 7)
        report(f"    k={k} P={P}: meas={float(meas):.4f}  E[mg|G_P]={cmean:.4f}  "
               f"bar={bar:.4f}  margin={cmean - bar:+.4f}  condRM rho*>={rho_lb:.4f} "
               f"{'OK' if rho_lb >= float(M_P) else 'FAIL'}")

    # (c) joint adversarial descent on (P, E) minimizing the conditional-RM bound:
    report("  (c) joint adversarial descent minimizing condRM rho* bound (per k):")
    worst = {}
    for k in range(8, 13):
        psz = 13 - k
        best_obj, best_PE, best_row = 1e9, None, None
        for trial in range(6):
            P = tuple(sorted(random.sample(range(1, 14), psz))) if psz else ()
            spread = random.choice([k + 1, k + 3, k + 6, 2 * k, 30])
            E = tuple([0] + sorted(random.sample(range(1, spread + 1), k - 1)))
            meas, cmean, cmu, cmin = cond_stats(P, E, n=12000)
            obj = float(meas) * (7.0 / 6) * (cmean - 1.0 / 7)
            for step in range(60):
                move = random.random()
                P2, E2 = P, E
                if move < 0.35 and psz:
                    P2 = list(P); i = random.randrange(psz)
                    cand = random.randrange(1, 14)
                    if cand not in P2:
                        P2[i] = cand; P2 = tuple(sorted(P2))
                    else:
                        continue
                else:
                    E2 = list(E); i = random.randrange(1, k)  # keep 0
                    cand = random.randrange(1, 41)
                    if cand not in E2:
                        E2[i] = cand; E2 = tuple([0] + sorted(E2[1:]))
                    else:
                        continue
                meas2, cmean2, cmu2, cmin2 = cond_stats(P2, E2, n=12000)
                obj2 = float(meas2) * (7.0 / 6) * (cmean2 - 1.0 / 7)
                if obj2 < obj - 1e-6:
                    P, E, obj = P2, E2, obj2
                    meas, cmean, cmu, cmin = meas2, cmean2, cmu2, cmin2
            if obj < best_obj:
                best_obj, best_PE = obj, (P, E)
                best_row = (meas, cmean, cmu, cmin)
        meas, cmean, cmu, cmin = best_row
        # re-measure the winner at high resolution
        meas, cmean, cmu, cmin = cond_stats(best_PE[0], best_PE[1], n=40000)
        best_obj = float(meas) * (7.0 / 6) * (cmean - 1.0 / 7)
        worst[k] = (best_PE, best_obj, cmean, float(meas), cmu)
        report(f"    k={k}: worst (P,E)={best_PE}  meas={float(meas):.4f}  "
               f"E[mg|G_P]={cmean:.4f}  mu|G_P={cmu:.4f}  condRM rho*>={best_obj:.4f} "
               f"{'OK' if best_obj >= float(M_P) else '< m_P FAIL'}")
    report("  NOTE: condRM uses maxgap>=1/7 only via the mean; direct mu|G_P (col above)")
    report("        times meas(G_P) is the true G2 -- compare both against m_P.")

    report("")
    report("=" * 78)
    report("PART 4 -- crux-class-constrained descent on plain E[maxgap] (vs T*)")
    report("=" * 78)
    seed_fams = [
        [2, 4, 6, 8, 10, 12, 13, 14, 16, 18, 20, 22, 24],     # death-star refuter
        list(range(14, 27)),                                   # consecutive saturated block
        [2, 6, 8, 9, 10, 12, 14, 15, 16, 18, 22, 26, 39],      # S57-style, made saturated
    ]
    best, bestv = 10.0, None
    for seed in seed_fams:
        v = sorted(seed); assert len(set(v)) == 13
        if not in_crux(v):
            report(f"    seed {v} not in crux class "
                   f"(sat={saturated(v)} prim={primitive(v)} ss={single_scale(v)}) -- skip")
            continue
        cur = plain_mu_and_mean(v, n=16000)[1]
        for step in range(400):
            i = random.randrange(13)
            lo, hi = 1, 13 * min(v)  # keep single-scale-able
            cand = random.randrange(max(1, min(v) // 2), min(2 * max(v), 3000))
            w = sorted(set(v[:i] + v[i + 1:] + [cand]))
            if len(w) != 13 or not in_crux(w):
                continue
            c = plain_mu_and_mean(w, n=8000)[1]
            if c < cur - 1e-5:
                v, cur = w, c
        cur = plain_mu_and_mean(v, n=40000)[1]
        report(f"    descended to {v}: E[mg]~{cur:.5f}")
        if cur < best:
            best, bestv = cur, v
    ex = Emaxgap_exact(bestv)
    report(f"  best crux-class family: {bestv}")
    report(f"  EXACT E[maxgap] = {ex} = {float(ex):.6f}")
    report(f"  vs 1/7:  {float(ex - THR):+.6f}    vs T* (m_P quantitative): "
           f"{float(ex - Tstar):+.6f}")

    report("")
    report("=" * 78)
    report("PART 5 -- Part-A finite-bridge shadow: arc-count growth of the good set")
    report("=" * 78)
    report("  #maximal arcs of {x: maxgap_E(x) > 1/7} for E={0..k-2, S} as spread S grows")
    report("  (Part A's correction is O(#arcs/Vmax); #arcs ~ k^2*spread makes it vacuous")
    report("   when cluster spread ~ Vmax)")
    k = 10
    for S in [12, 24, 48, 96, 192, 384]:
        E = list(range(k - 1)) + [S]
        xs = np.linspace(0, 1, 300000, endpoint=False) + 0.5 / 300000
        ph = np.mod(np.outer(xs, np.array(E, dtype=np.float64)), 1.0)
        ph.sort(axis=1)
        gaps = np.empty_like(ph)
        gaps[:, :-1] = np.diff(ph, axis=1)
        gaps[:, -1] = ph[:, 0] + 1.0 - ph[:, -1]
        good = (gaps.max(axis=1) > 1.0 / 7)
        arcs = int(np.sum(good[1:] & ~good[:-1])) + int(good[0])
        report(f"    spread S={S:4d}: #arcs ~ {arcs:5d}   arcs/S ~ {arcs / S:.2f}")
    report("")
    report("DONE.")
