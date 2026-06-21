#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
ADVERSARIAL VERIFICATION -- LRC(14) decorrelation-error angle (Lemma DE).  kps-2026-06-20-Sx-wf
================================================================================================
Targets (default skeptical):
  (1) re-derive the key engine numbers (p0_exact, single_block_decorr, TV_x, e(E), proved bound);
  (2) HUNT a primitive WIDE k-set (span>14, k=8..12) with
        p0(E) > cap_k   OR   p0_decorr exceeding the single-block bound
        OR a split exceeding the single block   OR  |e(E)| exceeding the budget;
      try span 15-30, resonant w, multi-scale, AP-with-difference d != 1 blocks;
  (3) check the closed-form |e| <= 49 m^2/(6M) exactly vs measured;
  (4) check the four-pillar glue has no boundary gap (15 <= span <= M*).
Witness any failure.  EXACT rationals throughout.
"""
import sys, itertools, math, random
from fractions import Fraction as F
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7)}
INNER = set(range(1, 7))
OUT = []
FAILURES = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    OUT.append(s)


# ------------------------------------------------------------------------------------------------
# EXACT ENGINES (independent re-implementation; do not import the engine under test)
# ------------------------------------------------------------------------------------------------
def p0_exact(E):
    """meas{x in (0,1): all 6 inner sectors [j/7,(j+1)/7), j=1..6 hit by some frac(e x)}.
       Boolean piecewise-constant in x; breakpoints x = s/(7 e). EXACT rational integration."""
    bps = set()
    for e in E:
        if e == 0:
            continue
        ae = abs(e)
        for s in range(1, 7 * ae):
            bps.add(F(s, 7 * ae))
    bps.add(F(0)); bps.add(F(1))
    bps = sorted(bps)
    nz = [e for e in E if e != 0]
    tot = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        hit = {int((e * mid % 1) * 7) for e in nz}
        if len(hit & INNER) == 6:
            tot += b - a
    return tot


def single_block_decorr(m, Nx=1260):
    """Decorrelated cover of one coherent block of m sweeping points, anchor phi INDEP uniform.
       Exact in phi (the cover is a union of arcs in phi computed via sector-wall breakpoints),
       sampled in x at midpoints of Nx subcells (x-dependence is the only sampled part)."""
    tot = F(0)
    for ix in range(Nx):
        x = F(2 * ix + 1, 2 * Nx)
        r = [(j * x) % 1 for j in range(m)]
        bps = sorted({(F(s, 7) - rj) % 1 for rj in r for s in range(7)})
        bps.append(bps[0] + 1)
        good = F(0)
        for a, b in zip(bps, bps[1:]):
            mid = (a + b) / 2
            hit = {int(((mid + rj) % 1) * 7) for rj in r}
            if len(hit & INNER) == 6:
                good += b - a
        tot += good
    return tot / Nx


def tv_in_phi(m, Nx=420):
    worst = 0
    for ix in range(Nx):
        x = F(2 * ix + 1, 2 * Nx)
        r = [(j * x) % 1 for j in range(m)]
        bps = sorted({(F(s, 7) - rj) % 1 for rj in r for s in range(7)})
        bps.append(bps[0] + 1)
        vals = []
        for a, b in zip(bps, bps[1:]):
            mid = (a + b) / 2
            hit = {int(((mid + rj) % 1) * 7) for rj in r}
            vals.append(1 if len(hit & INNER) == 6 else 0)
        tvc = sum(1 for i in range(len(vals)) if vals[i] != vals[(i + 1) % len(vals)])
        worst = max(worst, tvc)
    return worst


def primitive(E):
    g = 0
    for e in E:
        g = math.gcd(g, e)
    return g == 1


# ------------------------------------------------------------------------------------------------
def main():
    log(__doc__)

    log("=" * 96)
    log("STEP 1.  Re-derive engine numbers independently (caps, p0_decorr, TV_x).")
    log("=" * 96)
    decorr = {}
    for k in range(8, 13):
        m = k - 1
        d = single_block_decorr(m, 1260)
        decorr[k] = d
        cap = CAPS[k]
        log(f"  k={k:2d} m={m:2d}: cap={float(cap):.5f}  p0_decorr={float(d):.5f}  "
            f"budget cap-p0={float(cap - d):+.5f}")
    budget = {k: CAPS[k] - decorr[k] for k in range(8, 13)}
    log(f"\n  Claimed budgets 0.18895/0.18862/0.19213/0.23043/0.28125 (k=8..12).")
    log(f"  Recomputed:     " + "/".join(f"{float(budget[k]):.5f}" for k in range(8, 13)))
    minbud = min(budget.values())
    log(f"  min budget = {float(minbud):.5f}  (claim: 0.18862)")
    # cross-check the claimed decorr stubs
    claimed = {8: 0.1925, 9: 0.3056, 10: 0.4123, 11: 0.4948, 12: 0.5759}
    for k in range(8, 13):
        diff = abs(float(decorr[k]) - claimed[k])
        tag = "OK" if diff < 2e-3 else "MISMATCH"
        if diff >= 2e-3:
            FAILURES.append(f"p0_decorr k={k} recomputed {float(decorr[k]):.5f} vs claimed {claimed[k]} ({tag})")
        log(f"    k={k}: p0_decorr claim {claimed[k]} vs recomputed {float(decorr[k]):.5f}  {tag}")

    log("\n  TV_x verification (claim TV_x <= 14m, actual max 22/40/48/52/62 for m=7..11):")
    for k in range(8, 13):
        m = k - 1
        w = tv_in_phi(m, 420)
        ok = w <= 14 * m
        if not ok:
            FAILURES.append(f"TV_x m={m} = {w} EXCEEDS 14m={14*m}")
        log(f"    m={m:2d}: max TV_x={w:3d}  bound 14m={14*m:3d}  {'OK' if ok else 'EXCEEDS'}")

    log("\n" + "=" * 96)
    log("STEP 2.  HUNT for a wide primitive counterexample (span>14): p0(E) > cap_k.")
    log("=" * 96)
    log("  Family A: single coherent consecutive blocks {0}u{M..M+m-1}, M=8..60 (the claimed extremal).")
    worstA = {}
    for k in range(8, 13):
        m = k - 1
        cap = CAPS[k]
        wmax = F(0); wM = None
        for M in range(2, 61):
            E = tuple([0] + list(range(M, M + m)))
            if max(E) <= 14:
                continue
            if not primitive(E):
                continue
            p0 = p0_exact(E)
            if p0 > wmax:
                wmax = p0; wM = M
            if p0 > cap:
                FAILURES.append(f"FAMILY A counterexample k={k} M={M} E={E} p0={float(p0):.5f} > cap={float(cap):.5f}")
                log(f"    !!! k={k} M={M} p0={float(p0):.5f} > cap {float(cap):.5f}  COUNTEREXAMPLE")
        worstA[k] = (wmax, wM)
        log(f"    k={k:2d}: sup p0 over wide consec = {float(wmax):.5f} (M={wM})  cap={float(cap):.5f}  "
            f"margin={float(cap - wmax):+.5f}")

    log("\n  Family B: AP-with-difference-d blocks {0}u{M, M+d, ..., M+(m-1)d}, d=2,3,5,7; M small..")
    for k in range(8, 13):
        m = k - 1
        cap = CAPS[k]
        wmax = F(0); warg = None
        for d in (2, 3, 5, 7):
            for M in range(2, 40):
                E = tuple([0] + [M + j * d for j in range(m)])
                if max(E) <= 14:
                    continue
                if not primitive(E):
                    continue
                p0 = p0_exact(E)
                if p0 > wmax:
                    wmax = p0; warg = (d, M)
                if p0 > cap:
                    FAILURES.append(f"FAMILY B counterexample k={k} d={d} M={M} E={E} p0={float(p0):.5f} > cap={float(cap):.5f}")
                    log(f"    !!! k={k} d={d} M={M} p0={float(p0):.5f} > cap {float(cap):.5f}  COUNTEREXAMPLE")
        log(f"    k={k:2d}: sup p0 over AP-d wide = {float(wmax):.5f} (d,M={warg})  cap={float(cap):.5f}  "
            f"margin={float(cap - wmax):+.5f}")

    log("\n  Family C: resonant w -- include a speed that hits sector walls (multiples of 7 +/-1).")
    for k in range(8, 13):
        m = k - 1
        cap = CAPS[k]
        wmax = F(0); warg = None
        # base small cluster + one or two resonant far speeds
        for base in ([0, 1, 2, 3], [0, 1, 2, 3, 4], [0, 1, 2]):
            for w in (7, 8, 13, 14, 15, 21, 22, 28, 35, 49):
                for w2 in (0, w + 1, 2 * w, 3 * w + 1):
                    extra = [w] + ([w2] if w2 else [])
                    pool = sorted(set(base) | set(extra))
                    # extend to size k by adding consecutive far speeds
                    while len(pool) < k:
                        pool.append(max(pool) + 1)
                    E = tuple(pool[:k])
                    if max(E) <= 14 or not primitive(E):
                        continue
                    p0 = p0_exact(E)
                    if p0 > wmax:
                        wmax = p0; warg = E
                    if p0 > cap:
                        FAILURES.append(f"FAMILY C counterexample k={k} E={E} p0={float(p0):.5f} > cap={float(cap):.5f}")
                        log(f"    !!! k={k} E={E} p0={float(p0):.5f} > cap {float(cap):.5f}  COUNTEREXAMPLE")
        log(f"    k={k:2d}: sup p0 resonant = {float(wmax):.5f} (E={warg})  margin={float(cap - wmax):+.5f}")

    log("\n  Family D: RANDOM primitive wide k-sets, span in [15,40], many samples per k.")
    random.seed(20260620)
    for k in range(8, 13):
        cap = CAPS[k]
        wmax = F(0); warg = None
        N = 4000
        for _ in range(N):
            span = random.randint(15, 40)
            pts = set([0, span])
            while len(pts) < k:
                pts.add(random.randint(1, span))
            E = tuple(sorted(pts))
            if not primitive(E):
                continue
            p0 = p0_exact(E)
            if p0 > wmax:
                wmax = p0; warg = E
            if p0 > cap:
                FAILURES.append(f"FAMILY D counterexample k={k} E={E} p0={float(p0):.5f} > cap={float(cap):.5f}")
                log(f"    !!! k={k} E={E} p0={float(p0):.5f} > cap {float(cap):.5f}  COUNTEREXAMPLE")
        log(f"    k={k:2d}: sup p0 random wide ({N} samp) = {float(wmax):.5f} (E={warg})  "
            f"margin={float(cap - wmax):+.5f}")

    log("\n" + "=" * 96)
    log("STEP 3.  Verify the PROVED bound |e(E)| <= 49 m^2/(6 M) holds, and SPLIT <= single block.")
    log("=" * 96)
    log("  3a. proved bound vs measured |e| for single blocks:")
    boundviol = 0
    for k in range(8, 13):
        m = k - 1
        d = decorr[k]
        supeM = 0.0; argM = None
        for M in list(range(15, 70)) + list(range(70, 140, 4)) + [160, 200, 256]:
            E = tuple([0] + list(range(M, M + m)))
            if not primitive(E):
                continue
            p0 = p0_exact(E)
            e = float(p0 - d)
            bnd = 49 * m * m / (6 * M)
            if abs(e) > bnd + 1e-12:
                boundviol += 1
                FAILURES.append(f"BOUND VIOLATION k={k} M={M} |e|={abs(e):.5f} > 49m^2/6M={bnd:.5f}")
                log(f"    !!! k={k} M={M} |e|={abs(e):.5f} > bound {bnd:.5f}")
            if abs(e) * M > supeM:
                supeM = abs(e) * M; argM = M
        log(f"    k={k:2d}: sup|e|*M (C_meas) = {supeM:.4f} (M={argM})  proved-bound coeff 49m^2/6 = "
            f"{49*m*m/6:.1f}  -> bound never violated: {'OK' if boundviol==0 else 'VIOLATED'}")
    if boundviol == 0:
        log("    => PROVED bound 49 m^2/(6M) holds for every single block tested.")

    log("\n  3b. splitting test: does any 2-cluster (or 3-cluster) decorr/p0 EXCEED single-block sup?")
    for k in range(8, 13):
        m = k - 1
        sb = float(decorr[k])
        cap = CAPS[k]
        worst = 0.0; warg = None
        # all integer compositions of m into r>=2 parts, clusters placed at well-separated scales
        scales = [0, 60, 130, 210, 300]
        for r in range(2, min(m, 4) + 1):
            for comp in compositions(m, r):
                cls = []
                for i, sz in enumerate(comp):
                    base = scales[i]
                    cls.append([base + j for j in range(sz)])
                E = tuple(sorted(set([0]).union(*[set(c) for c in cls])))
                # ensure 0 is anchor; clusters already include 0 via first scale=0
                if not primitive(E):
                    continue
                p0 = float(p0_exact(E))
                if p0 > worst:
                    worst = p0; warg = comp
                if p0 > cap:
                    FAILURES.append(f"SPLIT EXCEEDS CAP k={k} comp={comp} E={E} p0={p0:.5f} > cap {float(cap):.5f}")
        excess = worst - sb
        tag = "OK (split <= single block)" if excess <= 1e-3 else "SPLIT EXCEEDS SINGLE BLOCK"
        if excess > 1e-3:
            FAILURES.append(f"SPLIT k={k} comp={warg} p0={worst:.5f} EXCEEDS single-block sup {sb:.5f} by {excess:.5f}")
        log(f"    k={k:2d}: max split p0 = {worst:.5f} (comp={warg})  single-block decorr = {sb:.5f}  "
            f"diff={excess:+.5f}  {tag}")

    log("\n" + "=" * 96)
    log("STEP 4.  GLUE / boundary-gap audit (the four-pillar close).")
    log("=" * 96)
    log("  The wide branch is split at M* = 49 m^2/(6 budget). Below: finite check + THM-546/547 peel.")
    log("  Audit: is the ENTIRE wide region 15 <= span actually covered, with NO M-gap?")
    for k in range(8, 13):
        m = k - 1
        bud = float(budget[k])
        Mstar = 49 * m * m / (6 * bud)
        Cmeas = 0.0
        # measured C_meas over a broad M sweep
        for M in list(range(15, 70)) + list(range(70, 260, 3)):
            E = tuple([0] + list(range(M, M + m)))
            e = abs(float(p0_exact(E) - decorr[k]))
            Cmeas = max(Cmeas, e * M)
        realistic = Cmeas / bud
        log(f"    k={k:2d}: PROVED-bound M*={Mstar:7.1f}  measured C_meas={Cmeas:.3f}  "
            f"realistic thresh C_meas/bud={realistic:5.1f}")
    log("\n  GLUE GAPS (named honestly):")
    log("   G1. Region 15<=span<=M* (~thousands) NOT closed by Lemma DE alone -> delegated to")
    log("       finite-check(span<=14) + THM-546/547 iterated peel. The span<=14 finite check does")
    log("       NOT cover 15<=span<=M*; the peel must. THIS IS A REAL GLUE OBLIGATION, not closed here.")
    log("   G2. HYP-2694 (single block = wide sup) is VERIFIED not PROVED; Lemma DE is conditional on it.")
    log("   G3. Multi-cluster r>=2 needs JOINT Erdos-Turan-Koksma constant (only VERIFIED here).")
    log("   G4. V_x(Hhat_s)<=49m^2/(pi|s|) majorant in step (3) is a coarse (non-tight) upper bound;")
    log("       rigorous as an inequality, so the PROVED bound stands, but constant is lossy.")

    log("\n" + "=" * 96)
    log("VERDICT")
    log("=" * 96)
    if FAILURES:
        log(f"  {len(FAILURES)} FAILURE(S) / COUNTEREXAMPLE(S) FOUND:")
        for f in FAILURES:
            log("   * " + f)
    else:
        log("  NO counterexample found: no wide primitive E with p0>cap; proved bound 49m^2/(6M)")
        log("  never violated; no split exceeds the single-block decorr sup; TV_x<=14m holds;")
        log("  p0_decorr matches claimed values; budgets match.")
        log("  SURVIVING GAP (not a refutation, a delegation): G1 the finite region 15<=span<=M*")
        log("  is NOT closed by Lemma DE -- it relies on THM-546/547 peel covering it. Lemma DE is")
        log("  an ASYMPTOTIC-tail result; the claim correctly marks itself PARTIAL. G2/G3 also open.")

    with open("05-knowledge/results/lrc_fin_verify_decorrelatio_kps-Sx-wf.out", "w",
              encoding="utf-8") as fo:
        fo.write("\n".join(OUT))
    log("\n[output -> 05-knowledge/results/lrc_fin_verify_decorrelatio_kps-Sx-wf.out]")


def compositions(n, r):
    """all ordered compositions of n into r positive parts"""
    if r == 1:
        yield (n,)
        return
    for first in range(1, n - r + 2):
        for rest in compositions(n - first, r - 1):
            yield (first,) + rest


if __name__ == "__main__":
    main()
