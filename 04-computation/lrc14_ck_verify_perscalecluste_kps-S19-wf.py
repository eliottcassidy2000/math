#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_ck_verify_perscalecluste_kps-S19-wf.py
========================================================================
ADVERSARIAL VERIFICATION of the claimed THEOREM:

  C(k) := sup_{E', w} w*|Delta_w(E')|  is BOUNDED,
          C(k) <= c*(#clusters(E')-1) <= c*(k-2),
  with explicit per-cluster slope c ~ 2.72 resting on a single-cluster
  saturation constant K_1 ~ 1.11.

EXACT (fractions.Fraction).  KEY EFFICIENCY: for a fixed core E', the
runs (a,b,s) are computed ONCE; then  w*Delta_w  is just a sum over those
runs of  G0(w*b - s/7) - G0(w*a - s/7)  -- only w varies.  This lets us
scan w over a huge resonant range cheaply, which is exactly where the
claim is most fragile.

FRONTS OF ATTACK:
  (1) re-derive engine exactly; cross-check task engine; telescoping (A).
  (2) HUNT for slope violations -- multi-scale resonant w (w = lcm of
      cluster scales & multiples), large #scales.
  (3) STRESS K_1 saturation: single block {0..L-1}, DEEP resonant w-scan.
  (4) STRESS linearity-in-r and report worst per-cluster slope.

A single violation -> holds=false with witness.

Run:
  python3 04-computation/lrc14_ck_verify_perscalecluste_kps-S19-wf.py 2>&1 | tee \
     05-knowledge/results/lrc14_ck_verify_perscalecluste_kps-S19-wf.out
"""
import sys
import itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def pr(*a):
    print(*a, flush=True)


# ---------- G_0 ----------
def G0(y):
    y = y - int(y)
    if y < 0:
        y += 1
    if y < F(1, 7):
        return y * F(6, 7)
    return F(6, 49) - (y - F(1, 7)) * F(1, 7)


def _f(y):
    y = y - int(y)
    if y < 0:
        y += 1
    return y


# ---------- breakpoints / cells / runs (computed ONCE per core) ----------
def breakpoints(Ep):
    Ep = sorted(set(e for e in Ep if e != 0))
    bp = {F(0), F(1)}
    for e in Ep:
        for j in range(7):
            c = F(j, 7)
            m = 0
            while True:
                xv = (c + m) / e
                if xv >= 1:
                    break
                if xv >= 0:
                    bp.add(xv)
                m += 1
    return sorted(b for b in bp if 0 <= b < 1)


def cells(Ep):
    Ep = sorted(set(e for e in Ep if e != 0))
    bp = breakpoints(Ep)
    out = []
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set(int(_f(e * mid) * 7) for e in Ep)
        miss = set(range(1, 7)) - hit
        out.append((lo, hi, next(iter(miss)) if len(miss) == 1 else None))
    return out


def runs(Ep):
    cm = cells(Ep)
    out = []
    i = 0
    n = len(cm)
    while i < n:
        lo, hi, s = cm[i]
        if s is None:
            i += 1
            continue
        a = lo
        j = i
        while j + 1 < n and cm[j + 1][2] == s and cm[j + 1][0] == cm[j][1]:
            j += 1
        b = cm[j][1]
        out.append((a, b, s))
        i = j + 1
    return out


def wDelta_signed_from_cells(cm, w):
    """signed w*Delta_w from precomputed cells (the raw engine)."""
    D = F(0)
    for lo, hi, s in cm:
        if s is None:
            continue
        D += G0(w * hi - F(s, 7)) - G0(w * lo - F(s, 7))
    return D


def wDelta_signed_from_runs(rn, w):
    D = F(0)
    for a, b, s in rn:
        D += G0(w * b - F(s, 7)) - G0(w * a - F(s, 7))
    return D


# task-supplied engine, EXACT copy (cross-check on a few cases only -- slow)
def G0_task(y):
    y = y - int(y); y = y + 1 if y < 0 else y
    return y * F(6, 7) if y < F(1, 7) else F(6, 49) - (y - F(1, 7)) * F(1, 7)


def wDelta_task(Ep, w):
    Ep = sorted(set(Ep)); bp = set([F(0), F(1)])
    for e in Ep:
        if e == 0:
            continue
        for j in range(7):
            c = F(j, 7); m = 0
            while True:
                xv = (c + m) / e
                if xv >= 1:
                    break
                if xv >= 0:
                    bp.add(xv)
                m += 1
    bp = sorted(b for b in bp if 0 <= b < 1); D = F(0)
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set(int((e * mid) % 1 * 7) for e in Ep)
        miss = set(range(1, 7)) - hit
        if len(miss) == 1:
            s = next(iter(miss)); D += G0_task(w * hi - F(s, 7)) - G0_task(w * lo - F(s, 7))
    return abs(float(D))


# ---------- helpers ----------
def prim(E):
    nz = [e for e in E if e != 0]
    return len(nz) > 0 and reduce(gcd, nz) == 1


def lcm(a, b):
    return a * b // gcd(a, b)


def lcml(xs):
    return reduce(lcm, xs, 1)


def clusters(Ep, ratio=3):
    es = sorted(e for e in Ep if e > 0)
    if not es:
        return []
    cl = [[es[0]]]
    for a, b in zip(es, es[1:]):
        if b > a * ratio:
            cl.append([b])
        else:
            cl[-1].append(b)
    return cl


def nclu(Ep):
    return len(clusters(Ep))


def scan_w(rn, wlist):
    worst, ww = 0.0, 0
    for w in wlist:
        v = abs(float(wDelta_signed_from_runs(rn, w)))
        if v > worst:
            worst, ww = v, w
    return worst, ww


def banner(t):
    pr("\n" + "=" * 72)
    pr(t)
    pr("=" * 72)


# ======================================================================
# (1)
# ======================================================================
def step1():
    banner("(1) ENGINE CROSS-CHECK + TELESCOPING IDENTITY (claim A)")
    cases = [
        ([0, 1, 2, 3, 4, 5, 6, 7], 384),
        ([0, 1, 2, 3, 40, 41, 42, 43], 320),
        ([0, 1, 2, 30, 31, 32, 60, 61], 480),
        ([0, 1, 2, 40, 41, 42, 80, 81, 82], 82),
        ([0, 1, 6, 7, 36, 37], 165),
    ]
    ok = True
    for core, w in cases:
        cm = cells(core)
        rn = runs(core)
        d_cell = wDelta_signed_from_cells(cm, w)
        d_run = wDelta_signed_from_runs(rn, w)
        d_task = wDelta_task(core, w)
        eng = (abs(float(d_cell)) == d_task)
        tel = (d_cell == d_run)
        ok = ok and eng and tel
        pr(f"  core={str(core):<30} w={w:4d}: cell={float(d_cell):+.6f} "
           f"task={d_task:.6f} eng={eng}  run={float(d_run):+.6f} tel={tel}")
    pr(f"  => engine+telescope exact: {ok}")
    return ok


# ======================================================================
# (3) K_1 single-cluster saturation -- DEEP w-scan (the load-bearing lemma)
# ======================================================================
def step3():
    banner("(3) K_1 single-cluster saturation: block {0..L-1}, DEEP w-scan")
    pr("  Claim: stays O(1) ~ 1.11 independent of L.  We scan w densely up to")
    pr("  a large bound and ALSO at lcm-resonant multiples.")
    pr(f"  {'L':>4}{'worst w|D|':>12}{'w*':>8}{'scan_wmax':>11}")
    K1 = 0.0
    wit = None
    for L in (6, 7, 8, 9, 10, 11, 12, 13):
        core = list(range(L))
        rn = runs(core)
        nz = [e for e in core if e > 0]
        Lall = lcml(nz)
        wmax = 5000
        cand = set(range(2, wmax))
        for m in range(1, 80):
            v = Lall * m
            if v < 2_000_000:
                cand.add(v)
        worst, ww = scan_w(rn, sorted(cand))
        if worst > K1:
            K1, wit = worst, (core, ww)
        pr(f"  {L:>4}{worst:>12.4f}{ww:>8}{wmax:>11}")
    pr(f"  => K_1 (deep scan) = {K1:.4f} at {wit}   (claim ~1.11)")
    return K1, wit


# ======================================================================
# (2)+(4) HUNT multi-scale resonant; report worst & slope
# ======================================================================
def step24():
    banner("(2)+(4) HUNT multi-scale resonant cores; worst w|D| and slope")
    pr(f"  {'core':<46}{'k':>3}{'r':>3}{'worst':>9}{'/(k-2)':>8}{'/(r-1)':>8}  w*")
    worstg = 0.0
    gw = None
    worst_slope = 0.0
    # designs: geometric scales 1, S, S^2 ... with small blocks
    designs = []
    for S in (10, 12, 20, 30, 50):
        for sz in (2, 3, 4):
            for r in (1, 2, 3, 4):
                core = set()
                for i in range(r):
                    base = 0 if i == 0 else S ** i
                    core |= set(range(base, base + sz))
                core = sorted(core)
                if 0 in core and len(core) <= 13 and prim(core):
                    designs.append(core)
    # plus a few hand designs with deeper resonance
    designs += [
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        [0, 1, 2, 14, 15, 16, 196, 197, 198],
        [0, 1, 2, 3, 21, 22, 23, 24, 441, 442],
    ]
    seen = set()
    for core in designs:
        core = sorted(set(core))
        key = tuple(core)
        if key in seen or 0 not in core or not prim(core) or len(core) > 13:
            continue
        seen.add(key)
        rn = runs(core)
        nz = [e for e in core if e > 0]
        cl = clusters(core)
        reps = [c[0] for c in cl if c[0] > 0]
        Lrep = lcml(reps)
        Lall = lcml(nz)
        cand = set(range(2, 1500))
        for base in (Lrep, Lall):
            for m in range(1, 60):
                v = base * m
                if v < 3_000_000:
                    cand.add(v)
        worst, ww = scan_w(rn, sorted(cand))
        k = len(core)
        r = nclu(core)
        slope = worst / max(r - 1, 1)
        worst_slope = max(worst_slope, slope)
        if worst > worstg:
            worstg, gw = worst, (core, ww)
        pr(f"  {str(core):<46}{k:>3}{r:>3}{worst:>9.4f}{worst/max(k-2,1):>8.3f}{slope:>8.3f}  {ww}")
    pr(f"\n  => worst w|D| = {worstg:.4f} at {gw}")
    pr(f"  => worst per-cluster slope (worst/(r-1)) = {worst_slope:.4f}  (claim c~2.72)")
    return worstg, gw, worst_slope


# ======================================================================
# (5) exhaustive small (independent re-check of claim E)
# ======================================================================
def step5():
    banner("(5) EXHAUSTIVE small primitive cores (re-check claim E worst=1.000)")
    worst = 0.0
    wc = None
    cnt = 0
    for k in (5, 6):
        for combo in itertools.combinations(range(0, 11), k):
            if 0 not in combo or not prim(combo):
                continue
            cnt += 1
            rn = runs(combo)
            for w in range(2, 56):
                v = abs(float(wDelta_signed_from_runs(rn, w)))
                if v > worst:
                    worst, wc = v, (combo, w)
    pr(f"  scanned {cnt} cores (k=5,6, entries<=10), w<56")
    pr(f"  worst w|D| = {worst:.4f} at {wc}  (claim 1.000)")
    # deeper w on the same pool to see if 1.000 is scan-limited:
    worst2 = 0.0
    wc2 = None
    for k in (5, 6):
        for combo in itertools.combinations(range(0, 11), k):
            if 0 not in combo or not prim(combo):
                continue
            rn = runs(combo)
            nz = [e for e in combo if e > 0]
            La = lcml(nz)
            cand = set(range(2, 400))
            for m in range(1, 30):
                if La * m < 100000:
                    cand.add(La * m)
            for w in cand:
                v = abs(float(wDelta_signed_from_runs(rn, w)))
                if v > worst2:
                    worst2, wc2 = v, (combo, w)
    pr(f"  DEEPER w-scan (to lcm multiples): worst = {worst2:.4f} at {wc2}")
    return worst, wc, worst2, wc2


def main():
    pr(__doc__)
    ok = step1()
    K1, K1it = step3()
    wg, gw, slope = step24()
    we, weit, we2, we2it = step5()

    banner("VERDICT")
    pr(f"  (1) engine + telescoping            : {'PASS' if ok else 'FAIL'}")
    pr(f"  (3) K_1 deep-scan saturation        : {K1:.4f} at {K1it}  (claim 1.11)")
    pr(f"  (2/4) worst multi-scale w|D|        : {wg:.4f} at {gw}")
    pr(f"        worst per-cluster slope       : {slope:.4f}  (claim c~2.72)")
    pr(f"  (5) exhaustive small worst (w<56)   : {we:.4f} at {weit}  (claim 1.000)")
    pr(f"      deeper-w small worst            : {we2:.4f} at {we2it}")
    overall = max(K1, wg, we, we2)
    pr(f"\n  OVERALL worst w|D| observed         : {overall:.4f}")
    pr("\n  ASSESSMENT:")
    pr(f"   - K_1 constant 1.11 certified?  {'NO -> ' if K1 > 1.20 else 'plausibly'} K_1={K1:.3f}")
    pr(f"   - slope c=2.72 a true sup?      {'NO -> ' if slope > 2.80 else 'survives'} slope={slope:.3f}")
    pr("   - STRUCTURAL claim (each cluster O(1) => C(k)=O(k)) is the dovetail need;")
    pr("     report whether worst w|D| stays bounded & ~linear in r across the hunt.")


if __name__ == "__main__":
    main()
