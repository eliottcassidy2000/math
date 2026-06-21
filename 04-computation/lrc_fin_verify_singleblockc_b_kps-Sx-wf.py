#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Follow-up adversarial checks for the single-block-cover-bound angle.
 (A) TASK 2c: does any decorrelated SPLIT exceed the single block?
 (B) The decorrelation-error budget: the prompt claims error ~0.01; measure it exactly.
 (C) Critical: is the single-block decorr p0_decorr even an UPPER BOUND for finite wide p0(E)?
     If finite p0 > p0_decorr by a margin comparable to / exceeding the cap-margin, the
     whole "wide cover <= single-block-decorr < cap" logic is BROKEN (it is not a majorant).
 (D) Wider/denser counterexample hunt: span up to ~60, many AP differences, two-scale, random.
"""
import sys, itertools, math, random
from fractions import Fraction as F
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7)}
INNER = set(range(1, 7))


def direct_p0(E):
    E = sorted(set(E))
    diffs = set(E)
    for a, b in itertools.combinations(E, 2):
        diffs.add(abs(a - b))
    diffs.discard(0)
    xb = set()
    for d in diffs:
        for a in range(0, 7 * d + 1):
            xb.add(F(a, 7 * d))
    xb = sorted(b for b in xb if 0 <= b <= 1)
    if not xb or xb[0] != F(0):
        xb = [F(0)] + xb
    if xb[-1] != F(1):
        xb = xb + [F(1)]
    tot = F(0)
    for a, b in zip(xb, xb[1:]):
        if b == a:
            continue
        xmid = (a + b) / 2
        hit = {int((F(e) * xmid % 1) * 7) for e in E}
        if len(hit & INNER) == 6:
            tot += b - a
    return tot


def is_primitive(E):
    g = 0
    for e in E:
        g = math.gcd(g, e)
    return g == 1


def multiblock_decorr(blocks, Nx, Nphi):
    tot = 0.0
    grid = [(p + 0.5) / Nphi for p in range(Nphi)]
    for ix in range(Nx):
        x = (ix + 0.5) / Nx
        bres = [[(d * x) % 1 for d in B] for B in blocks]
        cnt = 0
        N = Nphi ** len(blocks)
        for phis in itertools.product(grid, repeat=len(blocks)):
            hit = set()
            for i, B in enumerate(bres):
                ph = phis[i]
                for rj in B:
                    hit.add(int(((ph + rj) % 1) * 7))
            if len(hit & INNER) == 6:
                cnt += 1
        tot += cnt / N
    return tot / Nx


# exact single block decorr (x-exact)
def measure_union_arcs(intervals):
    segs = []
    for a, b in intervals:
        a %= 1; b %= 1
        if a < b: segs.append((a, b))
        elif a > b: segs.append((a, F(1))); segs.append((F(0), b))
    segs.sort()
    tot = F(0); cur_a = cur_b = None
    for a, b in segs:
        if cur_a is None: cur_a, cur_b = a, b
        elif a <= cur_b:
            if b > cur_b: cur_b = b
        else: tot += cur_b - cur_a; cur_a, cur_b = a, b
    if cur_a is not None: tot += cur_b - cur_a
    return tot


def phi_cover_at_x(m, x):
    r = [(j * x) % 1 for j in range(m)]
    bps = sorted({(F(s, 7) - rj) % 1 for rj in r for s in range(7)}); bps.append(bps[0] + 1)
    good = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        hit = {int(((mid + rj) % 1) * 7) for rj in r}
        if len(hit & INNER) == 6: good += b - a
    return good


def exact_single_block(m):
    xb = set()
    for d in range(1, m):
        for a in range(0, 7 * d + 1):
            xb.add(F(a, 7 * d))
    xb = sorted(b for b in xb if 0 <= b <= 1)
    if xb[0] != F(0): xb = [F(0)] + xb
    if xb[-1] != F(1): xb = xb + [F(1)]
    tot = F(0)
    for a, b in zip(xb, xb[1:]):
        if b == a: continue
        tot += phi_cover_at_x(m, (a + b) / 2) * (b - a)
    return tot


def main():
    print(__doc__)
    EXACT = {m: exact_single_block(m) for m in range(7, 12)}

    print("=" * 78)
    print("(A) TASK 2c: decorrelated SPLITS vs single block (grid). Does any split exceed?")
    print("=" * 78)
    for k in (8, 9, 10, 11, 12):
        m = k - 1
        sb = float(EXACT[m])
        print(f"  k={k} m={m}: single-block(exact)={sb:.5f}")
        parts = []
        def gen(rem, mx, cur):
            if rem == 0: parts.append(list(cur)); return
            for p in range(min(rem, mx), 0, -1): gen(rem - p, p, cur + [p])
        gen(m, m, [])
        maxsplit = 0.0; argmax = None
        for P in parts:
            if len(P) == 1: continue
            if len(P) > 4: continue  # cost control; many-block splits are far lower anyway
            blocks = [tuple(range(sz)) for sz in P]
            v = multiblock_decorr(blocks, 84, 14)
            if v > maxsplit: maxsplit, argmax = v, P
            if v > sb + 0.015:
                print(f"      *** split {P}: {v:.5f}  EXCEEDS single block ***")
        print(f"      max split (<=4 blocks) = {maxsplit:.5f} at {argmax}  "
              f"{'<= single' if maxsplit <= sb + 0.015 else 'EXCEEDS'}")

    print()
    print("=" * 78)
    print("(C) CRITICAL: is single-block-decorr an UPPER BOUND for finite wide p0(E)?")
    print("    Compare exact finite p0 of wide AP-d=1 (consec) blocks vs decorr limit.")
    print("=" * 78)
    print("    If finite p0 >> decorr, then 'p0(E) <= p0_decorr < cap' is NOT the proof;")
    print("    the real safety is p0(E) < cap DIRECTLY (decorr is a red herring as a majorant).")
    for k in range(8, 13):
        m = k - 1
        decorr = float(EXACT[m])
        worst_finite = F(0); wE = None
        # scan consecutive blocks for offsets 0..40 (wide)
        for M in range(15, 41):
            E = tuple(range(M, M + k))
            p0 = direct_p0(E)
            if p0 > worst_finite: worst_finite, wE = p0, E
        wf = float(worst_finite)
        print(f"  k={k}: max finite consec p0(wide)={wf:.5f}  decorr_limit={decorr:.5f}  "
              f"finite-MINUS-decorr={wf-decorr:+.5f}  cap={float(CAPS[k]):.5f}  "
              f"{'decorr is NOT an upper bound (finite>decorr)' if wf>decorr+1e-9 else 'decorr>=finite'}")

    print()
    print("=" * 78)
    print("(D) EXTENDED counterexample hunt: span up to ~60, AP-d up to 9, two-scale, random.")
    print("=" * 78)
    random.seed(7)
    worst = {k: (F(0), None) for k in range(8, 13)}
    tested = 0
    found_ce = False
    for k in range(8, 13):
        cands = []
        # AP with difference d, various offsets
        for d in range(1, 10):
            for off in range(0, 7):
                cands.append(tuple(off + d * i for i in range(k)))
        # consecutive at offsets up to 50 (wide)
        for start in range(15, 51):
            cands.append(tuple(range(start, start + k)))
        # two-scale: r small + far block
        for far in (15, 17, 19, 23, 29, 31, 50, 70, 100):
            for split in range(1, k):
                cands.append(tuple(range(split)) + tuple(far + i for i in range(k - split)))
        # near-resonant with 7
        for base in (6, 7, 8, 13, 14, 15):
            cands.append(tuple(base * i + 1 for i in range(k)))
        # random wide primitive
        for _ in range(2000):
            E = tuple(sorted(random.sample(range(0, 60), k)))
            cands.append(E)
        for E in cands:
            E = tuple(sorted(set(E)))
            if len(E) != k: continue
            if not is_primitive(E): continue
            if max(E) - min(E) <= 14: continue
            p0 = direct_p0(E)
            tested += 1
            if p0 > worst[k][0]: worst[k] = (p0, E)
            if p0 > CAPS[k]:
                found_ce = True
                print(f"  !!! COUNTEREXAMPLE k={k}: E={E} span={max(E)-min(E)} "
                      f"p0={float(p0):.6f} > cap={float(CAPS[k]):.6f}")
    print(f"  tested {tested} wide primitive sets; counterexample found = {found_ce}")
    for k in range(8, 13):
        p0, E = worst[k]
        print(f"  k={k}: max wide p0={float(p0):.6f}  cap={float(CAPS[k]):.6f}  "
              f"margin={float(CAPS[k]-p0):.6f}  E={E}")


if __name__ == "__main__":
    main()
