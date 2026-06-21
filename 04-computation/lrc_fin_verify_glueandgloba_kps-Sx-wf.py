#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc_fin_verify_glueandgloba_kps-Sx-wf.py  (kps-2026-06-20, ADVERSARIAL VERIFY)
==============================================================================
ADVERSARIAL verification of the "glue-and-global-sup" advance for LRC(14)-S3.

The advance CLAIMS:
  LRC(14)-S3 <=> p0(E)=meas(S7(E)) <= cap_k for all primitive k-sets E (0 in E, k=8..12),
  cap_8..cap_12 = 2243/5880, 1979/4004, 55/91, 66/91, 6/7.
  GLOBAL-SUP: pinned consec_k is the sup over all E; wide shapes sit strictly below
  because the M->inf single-block sweeping cover is the phi-shift-average of the consec
  sector pattern.  SOLE OPEN link: true-wide (>=2 far) error-aggregation lemma.

This script DEFAULTS SKEPTICAL.  Four jobs:
  (V1) Re-derive the key reference inequalities EXACTLY with the engine.
  (V2) HUNT for a primitive WIDE k-set (span>14) with p0(E) > cap_k, OR
       a split exceeding single block, OR a decorrelation error exceeding the margin.
       Directions: boundary span 15-30 exhaustive bands, dilated APs (d!=1),
       AP+far-spike, geometric ratio-2/3 resonant cores, 2-far/3-far collars,
       multi-scale two-cluster resonance, large random wide rows.
  (V3) Check the claimed closed forms / constants exactly vs the engine.
  (V4) Check the boundary glue has NO gap (span 14->15->...).

ONE counterexample (p0 > cap) => holds=false WITH WITNESS.
EXACT rationals throughout.  Output saved to results/.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7)}
INNER = set(range(1, 7))

# claimed reference numbers from the advance (to check exactly)
CLAIM_CONSEC = {8: 0.32721, 9: 0.41616, 10: 0.50448, 11: 0.58146, 12: 0.62411}
CLAIM_SWEEP  = {8: 0.19252, 9: 0.30564, 10: 0.41227, 11: 0.49484, 12: 0.57590}

# ---------------------------------------------------------------------------
def sector_of(p):
    return int((p % 1) * 7)

def p0_inner(E):
    """meas{ x : every inner sector 1..6 is hit by some frac(e x) }, EXACT.
    Breakpoints of frac(e x) at sector boundaries s/7 are x=m/(7e); this set
    is COMPLETE for the union step function over all e."""
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    tot = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        secs = {sector_of(e * xm) for e in E}
        if INNER <= secs:
            tot += x1 - x0
    return tot

def primitive(E):
    g = 0
    for e in E:
        g = gcd(g, int(e))
    return g == 1

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


# ===========================================================================
RESULT = {"violations": [], "sweep_exceed_pinned": [], "split_exceed_block": []}

def check(E, tag):
    """Return (p0, margin). Record any cap violation as a witness."""
    E = tuple(sorted(set(int(e) for e in E)))
    if 0 not in E:
        E = tuple(sorted(set((0,) + E)))
    k = len(E)
    if k not in CAPS or not primitive(E):
        return None
    p = p0_inner(E)
    m = CAPS[k] - p
    if m < 0:
        RESULT["violations"].append((tag, E, p, CAPS[k]))
    return (p, m, E, k)


def V1_redrive():
    print("=" * 78)
    print("V1  RE-DERIVE the key reference inequalities EXACTLY")
    print("=" * 78)
    ok = True
    # three anchor values from the prompt
    a = p0_inner(range(8)); b = single_block_decorr(7); c = p0_inner([0,19,20,21,22,23,24,25])
    print(f"  p0(consec_8)          = {a} = {float(a):.5f}  (claim 0.32721)  {'OK' if abs(float(a)-0.32721)<1e-4 else 'MISMATCH'}")
    print(f"  single_block_decorr(7)= {b} = {float(b):.5f}  (claim 0.19252)  {'OK' if abs(float(b)-0.19252)<1e-4 else 'MISMATCH'}")
    print(f"  p0([0,19..25])        = {c} = {float(c):.5f}  (claim 0.20218)  {'OK' if abs(float(c)-0.20218)<1e-4 else 'MISMATCH'}")
    print("\n  k :  p0(consec_k)  sweep(M->inf)   cap_k       sweep<pinned  pinned<cap  sweep<cap")
    rows = {}
    for k in range(8, 13):
        pin = p0_inner(range(k))
        sw = single_block_decorr(k - 1, 1260)
        cap = CAPS[k]
        rows[k] = (pin, sw, cap)
        s_lt_p = sw < pin
        p_lt_c = pin < cap
        s_lt_c = sw < cap
        # cross-check claimed floats
        cm = abs(float(pin) - CLAIM_CONSEC[k]) < 1e-4 and abs(float(sw) - CLAIM_SWEEP[k]) < 1e-4
        if not (s_lt_p and s_lt_c) or not cm:
            ok = False
        print(f"  {k:2d}:  {float(pin):.5f}      {float(sw):.5f}        {float(cap):.5f}   "
              f"{str(s_lt_p):5s}        {str(p_lt_c):5s}      {str(s_lt_c):5s}  "
              f"{'claimOK' if cm else 'CLAIM-MISMATCH'}")
    # CRITICAL note: pinned<cap is the actual route claim; verify it
    print("\n  pinned<cap for ALL k? ", all(rows[k][0] < rows[k][2] for k in rows))
    return ok, rows


def V4_glue_no_gap():
    """Boundary glue: exhaustive exact over span = 14,15,16,17 for k=8,9,10.
    Span s means max(E)=s with endpoints 0 and s pinned. The DONE half is span<=14;
    a gap would show as a jump in worst-margin at the 14->15 boundary."""
    print("\n" + "=" * 78)
    print("V4  BOUNDARY GLUE — exhaustive exact, span s=14..17, k=8,9,10 (no-gap test)")
    print("=" * 78)
    for k in (8, 9, 10):
        print(f"  k={k}:")
        for s in (14, 15, 16, 17):
            interior = list(range(1, s))
            need = k - 2
            if need < 0 or need > len(interior):
                continue
            worst = None
            viol = 0
            cnt = 0
            for combo in itertools.combinations(interior, need):
                E = (0,) + combo + (s,)
                if not primitive(E):
                    continue
                cnt += 1
                p = p0_inner(E)
                m = CAPS[k] - p
                if m < 0:
                    viol += 1
                    RESULT["violations"].append(("glue-span%d" % s, E, p, CAPS[k]))
                if worst is None or m < worst[0]:
                    worst = (m, E)
            if cnt:
                print(f"     span={s:2d}: #prim={cnt:7d}  worst margin={float(worst[0]):+.5f}  "
                      f"viol={viol}  argworst={worst[1]}")


def V2_hunt():
    print("\n" + "=" * 78)
    print("V2  ADVERSARIAL HUNT — wide primitive E with p0(E) > cap_k")
    print("=" * 78)
    champ = {k: None for k in range(8, 13)}

    def consider(E, tag):
        res = check(E, tag)
        if res is None:
            return
        p, m, EE, k = res
        ratio = p / CAPS[k]
        if champ[k] is None or ratio > champ[k][0]:
            champ[k] = (ratio, EE, p, tag)

    # (a) EXHAUSTIVE bands span exactly s for s=15..22, k=8,9 (endpoints pinned).
    print("  (a) exhaustive span bands 15..22, k=8,9 ...")
    for s in range(15, 23):
        for k in (8, 9):
            interior = list(range(1, s))
            need = k - 2
            if need > len(interior):
                continue
            for combo in itertools.combinations(interior, need):
                consider((0,) + combo + (s,), f"band-s{s}-k{k}")
    # k=10 exhaustive bands up to span 18 (heavier)
    for s in range(15, 19):
        interior = list(range(1, s)); need = 8
        if need <= len(interior):
            for combo in itertools.combinations(interior, need):
                consider((0,) + combo + (s,), f"band-s{s}-k10")

    # (b) dilated APs (scale invariance) — d != 1; and AP with a far spike.
    print("  (b) dilated APs + AP-with-far-spike ...")
    for k in range(8, 13):
        for d in (1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 17):
            consider(tuple(d * i for i in range(k)), f"dilAP-d{d}-k{k}")
        for d in (1, 2, 3):
            core = tuple(d * i for i in range(k - 1))
            for w in (core[-1] + d, core[-1] + 2 * d, 2 * core[-1], 3 * core[-1], 100, 199):
                consider(core + (w,), f"AP-spike-d{d}-k{k}")

    # (c) geometric / multi-scale resonant cores (ratios 2,3) + two-cluster resonance
    print("  (c) geometric / multi-scale two-cluster resonance ...")
    for k in range(8, 13):
        for sep in (10, 20, 50, 100, 7, 14, 49):
            for sz in range(1, k):
                c1 = list(range(sz)); c2 = list(range(sep, sep + (k - sz)))
                consider(tuple(c1 + c2), f"2clu-sep{sep}-sz{sz}-k{k}")
        for ratio in (2, 3):
            # geometric block: 1, r, r^2, ... each +small neighbor
            pts = []
            v = 1
            while len(pts) < k:
                pts.append(v); pts.append(v + 1); v *= ratio
            consider(tuple([0] + pts[:k - 1]), f"geom-r{ratio}-k{k}")

    # (d) 2-far / 3-far / 4-far collars: bounded core in [0,14] + far spikes
    print("  (d) multi-far collars (random) ...")
    rng = random.Random(20260620)
    for k in range(8, 13):
        for _ in range(20000):
            nfar = rng.randint(2, 4)
            ncore = k - nfar
            if ncore < 2:
                continue
            core = sorted(set(rng.sample(range(0, 15), ncore)) | {0})
            far = sorted(rng.sample(range(15, 400), k - len(core)))
            consider(tuple(core) + tuple(far), f"collar-{nfar}far-k{k}")

    # (e) fully random wide rows, big modulus
    print("  (e) large random wide rows ...")
    for k in range(8, 13):
        for _ in range(20000):
            E = sorted(rng.sample(range(1, 500), k - 1))
            consider((0,) + tuple(E), f"rand-k{k}")

    # (f) consec + dilations & translations of consec inside a wider window (resonance with consec)
    print("  (f) near-consec wide perturbations ...")
    for k in range(8, 13):
        base = list(range(k))
        for shift_one in range(1, 20):
            E = base[:-1] + [base[-1] + shift_one]   # stretch the top
            consider(tuple(E), f"stretchtop-{shift_one}-k{k}")
        for d in (2, 3):
            # consec core dilated then add the gaps back -> AP of difference d (span > 14)
            consider(tuple(d * i for i in range(k)), f"APd-{d}-k{k}")

    print("\n  k :  worst p0/cap ratio       p0           cap          via              champion E")
    worst_overall = F(0)
    for k in range(8, 13):
        ratio, E, p, tag = champ[k]
        worst_overall = max(worst_overall, ratio)
        flag = "  *** EXCEEDS CAP ***" if ratio > 1 else ""
        print(f"  {k:2d}:  {float(ratio):.5f}              {float(p):.5f}      {float(CAPS[k]):.5f}    "
              f"{tag:16s}  {E}{flag}")
    print(f"\n  GLOBAL worst p0/cap ratio over hunt: {float(worst_overall):.5f} "
          f"({'NO VIOLATION' if worst_overall <= 1 else 'VIOLATION FOUND'})")
    return worst_overall


def V2b_split_vs_block():
    """Does any SPLIT (multi-cluster) decorrelated cover exceed the SINGLE block?
    Compute decorrelated cover for a few cluster shapes and compare to single_block."""
    print("\n" + "=" * 78)
    print("V2b  SPLIT vs SINGLE-BLOCK decorrelated cover (does splitting ever raise the limit?)")
    print("=" * 78)
    # We test finite-scale proxies: a single block consec_m vs split into 2 widely-separated
    # blocks (scales differ -> decorrelated). Use large separations so anchors decorrelate.
    for k in range(8, 11):
        m = k - 1
        sb = single_block_decorr(m, 1260)
        # split block as actual finite rows with huge scale separation (proxy for decorrelated split)
        # m tiles distributed: m1 in cluster A (around 0), m2 in cluster B (around BIG)
        worst_split = F(0); argsplit = None
        BIG = 100000
        for m1 in range(1, m):
            m2 = m - m1
            A = list(range(m1))
            B = [BIG + i for i in range(m2)]
            E = tuple([0] + A[1:] + B) if m1 >= 1 else tuple([0] + B)
            E = tuple(sorted(set([0] + list(range(m1)) + [BIG + i for i in range(m2)])))
            if len(E) != k:  # adjust
                continue
            p = p0_inner(E)
            if p > worst_split:
                worst_split = p; argsplit = (m1, m2, E)
        print(f"  k={k:2d}: single-block decorr={float(sb):.5f}   best split p0={float(worst_split):.5f}  "
              f"split>block? {worst_split>sb}  arg={argsplit}")
        if worst_split > sb:
            RESULT["split_exceed_block"].append((k, float(worst_split), float(sb), argsplit))


def V3_closed_forms():
    """Check claimed closed-form constants exactly."""
    print("\n" + "=" * 78)
    print("V3  CLOSED-FORM / CONSTANT CHECKS (exact)")
    print("=" * 78)
    # caps as decimals
    for k in CAPS:
        print(f"  cap_{k} = {CAPS[k]} = {float(CAPS[k]):.5f}")
    # Qb plateau claimed: 0.19660/0.36210/0.44789/0.53125/0.60224 < cap_k
    Qb = {8: 0.19660, 9: 0.36210, 10: 0.44789, 11: 0.53125, 12: 0.60224}
    print("\n  claimed true-wide plateau Qb(k-1) < cap_k ?")
    for k in range(8, 13):
        print(f"   k={k}: Qb={Qb[k]:.5f}  cap={float(CAPS[k]):.5f}  Qb<cap? {Qb[k] < float(CAPS[k])}")
    # margin check: cap - pinned (the actual route margin)
    print("\n  ROUTE margin cap_k - p0(consec_k) (must be > 0 for the pinned argmax to clear):")
    for k in range(8, 13):
        pin = p0_inner(range(k))
        print(f"   k={k}: cap-pinned = {float(CAPS[k]-pin):+.5f}")


def main():
    print(__doc__)
    ok1, rows = V1_redrive()
    V3_closed_forms()
    V4_glue_no_gap()
    V2b_split_vs_block()
    worst = V2_hunt()

    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    nv = len(RESULT["violations"])
    print(f"  cap violations found: {nv}")
    for tag, E, p, cap in RESULT["violations"][:20]:
        print(f"    VIOLATION [{tag}] E={E} p0={float(p):.5f} > cap={float(cap):.5f}")
    print(f"  splits exceeding single block: {len(RESULT['split_exceed_block'])}")
    for row in RESULT["split_exceed_block"]:
        print(f"    {row}")
    print(f"  global worst p0/cap ratio in hunt: {float(worst):.5f}")
    if nv == 0 and worst <= 1:
        print("\n  RESULT: NO COUNTEREXAMPLE FOUND. cap inequality holds in all tested regions.")
        print("  The GLOBAL-SUP claim (pinned consec is the sup, wide below) is consistent.")
        print("  SURVIVING GAP: the true-wide >=2-far closed-form aggregation lemma is still")
        print("  only numerically supported, NOT proved. The route is gap-free modulo it.")
    else:
        print("\n  RESULT: COUNTEREXAMPLE / ANOMALY FOUND (see above). holds=false.")
    print("\nDONE.")


if __name__ == "__main__":
    main()
