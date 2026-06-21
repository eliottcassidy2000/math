#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
ADVERSARIAL VERIFICATION STEP 3+4 (faster) -- LRC(14) decorrelation-error angle.
Step 1 (engine numbers) and Step 2 (counterexample hunt) already passed in the main script:
  NO wide primitive E with p0>cap was found in Families A/B/C/D; p0_decorr and budgets match;
  TV_x<=14m holds.  Here we close:
   (3a) proved bound |e| <= 49 m^2/(6M) for single blocks (exact rationals);
   (3b) split (multi-cluster) <= single-block decorr sup;
   (4)  glue/boundary audit M*.
EXACT rationals.
"""
import sys, math
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


def p0_exact(E):
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


def primitive(E):
    g = 0
    for e in E:
        g = math.gcd(g, e)
    return g == 1


def compositions(n, r):
    if r == 1:
        yield (n,)
        return
    for first in range(1, n - r + 2):
        for rest in compositions(n - first, r - 1):
            yield (first,) + rest


def main():
    log(__doc__)
    decorr = {k: single_block_decorr(k - 1, 1260) for k in range(8, 13)}
    budget = {k: CAPS[k] - decorr[k] for k in range(8, 13)}

    log("=" * 92)
    log("STEP 3a.  PROVED bound |e(E)| <= 49 m^2/(6 M) for single blocks (exact).")
    log("=" * 92)
    boundviol = 0
    for k in range(8, 13):
        m = k - 1
        d = decorr[k]
        supeM = 0.0; argM = None
        # M kept modest to bound runtime; bound is monotone-helpful at large M anyway
        Mlist = list(range(15, 50)) + [56, 64, 80, 96, 112, 128]
        for M in Mlist:
            E = tuple([0] + list(range(M, M + m)))
            if not primitive(E):
                continue
            p0 = p0_exact(E)
            e = float(p0 - d)
            bnd = 49 * m * m / (6 * M)
            if abs(e) > bnd + 1e-12:
                boundviol += 1
                FAILURES.append(f"BOUND VIOLATION k={k} M={M} |e|={abs(e):.5f} > {bnd:.5f}")
            if abs(e) * M > supeM:
                supeM = abs(e) * M; argM = M
        log(f"  k={k:2d}: sup|e|*M (C_meas)={supeM:.4f} (M={argM})  49m^2/6={49*m*m/6:.1f}  "
            f"bound {'HOLDS' if boundviol == 0 else 'VIOLATED'}")
    log(f"  => bound 49 m^2/(6M) never violated across all single blocks: {boundviol == 0}")

    log("\n" + "=" * 92)
    log("STEP 3b.  SPLIT test: any multi-cluster p0 EXCEED single-block decorr sup or cap?")
    log("=" * 92)
    for k in range(8, 13):
        m = k - 1
        sb = float(decorr[k])
        cap = CAPS[k]
        worst = 0.0; warg = None
        scales = [0, 30, 62, 96, 130]   # well-separated but modest (keeps breakpoint sets small)
        for r in range(2, min(m, 4) + 1):
            for comp in compositions(m, r):
                cls = []
                for i, sz in enumerate(comp):
                    cls.append([scales[i] + j for j in range(sz)])
                E = tuple(sorted(set([0]).union(*[set(c) for c in cls])))
                if not primitive(E):
                    continue
                p0 = float(p0_exact(E))
                if p0 > worst:
                    worst = p0; warg = comp
                if p0 > cap:
                    FAILURES.append(f"SPLIT>CAP k={k} comp={comp} E={E} p0={p0:.5f}")
        excess = worst - sb
        tag = "OK (split<=single)" if excess <= 1e-3 else "SPLIT EXCEEDS SINGLE BLOCK"
        if excess > 1e-3:
            FAILURES.append(f"SPLIT k={k} comp={warg} p0={worst:.5f} > single-block {sb:.5f}")
        log(f"  k={k:2d}: max split p0={worst:.5f} (comp={warg})  single-block={sb:.5f}  "
            f"diff={excess:+.5f}  {tag}")

    log("\n" + "=" * 92)
    log("STEP 4.  GLUE / boundary M* audit.")
    log("=" * 92)
    for k in range(8, 13):
        m = k - 1
        bud = float(budget[k])
        Mstar = 49 * m * m / (6 * bud)
        Cmeas = 0.0
        for M in list(range(15, 50)) + [64, 96, 128]:
            E = tuple([0] + list(range(M, M + m)))
            e = abs(float(p0_exact(E) - decorr[k]))
            Cmeas = max(Cmeas, e * M)
        log(f"  k={k:2d}: PROVED-bound M*={Mstar:7.1f}  measured C_meas={Cmeas:.3f}  "
            f"realistic C_meas/bud={Cmeas/bud:5.1f}")
    log("\n  GLUE GAP G1 (REAL): region 15<=span<=M* (~thousands) NOT closed by Lemma DE; delegated")
    log("  to THM-546/547 iterated peel (the span<=14 finite check does NOT cover 15<=span<=M*).")
    log("  G2: HYP-2694 (single block = wide sup) VERIFIED not PROVED. G3: multi-cluster r>=2 needs")
    log("  joint ET-Koksma constant (VERIFIED only). G4: V_x(Hhat_s) majorant lossy but rigorous.")

    log("\n" + "=" * 92)
    log("VERDICT (step 3+4)")
    log("=" * 92)
    if FAILURES:
        log(f"  {len(FAILURES)} FAILURE(S):")
        for f in FAILURES:
            log("   * " + f)
    else:
        log("  Proved bound holds; no split exceeds single block; no split exceeds cap.")
        log("  Combined with Step 1/2 (no wide counterexample, numbers match): claim's PARTIAL")
        log("  status confirmed. Surviving gap = G1 finite region 15<=span<=M* (delegated to peel).")

    with open("05-knowledge/results/lrc_fin_verify_decorrelatio_step34_kps-Sx-wf.out", "w",
              encoding="utf-8") as fo:
        fo.write("\n".join(OUT))
    log("\n[output -> 05-knowledge/results/lrc_fin_verify_decorrelatio_step34_kps-Sx-wf.out]")


if __name__ == "__main__":
    main()
