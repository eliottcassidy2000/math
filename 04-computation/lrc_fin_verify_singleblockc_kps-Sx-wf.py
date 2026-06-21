#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL VERIFICATION of the "single-block-cover-bound" angle for LRC(14).
============================================================================
Claim under test (angle "single-block-cover-bound", RESULT PROVED):
  The single-block decorrelated cover p0_decorr(m) equals an exact finite
  inclusion-exclusion, and p0_decorr(m) < cap_{m+1} for all needed k=m+1 in {8..12},
  with claimed exact rationals
     p0_decorr = 283/1470, 629/2058, 16969/41160, 30551/61740, 71111/123480  (m=7..11)
  against caps 2243/5880, 1979/4004, 55/91, 66/91, 6/7, minimum margin 111019/588588 > 0.18.
  Single-arc closed form A_single(w,m)=(7-w)^2/(49(m-1)) for w in {4,5,6}.
  Pinned-consec dominance: p0_decorr(m) = E_phi[coverage_phi(consec_m)] <= sup_phi.

Tasks:
 (1) re-derive p0_decorr with the prompt engine AND an independent exact engine; compare to claimed rationals.
 (2) HUNT for a wide primitive k-set (span>14) with p0(E) > cap_k, or p0_decorr exceeding single-block,
     or split exceeding single block, or decorrelation error exceeding margin.
 (3) check the single-arc closed form exactly vs the engine.
 (4) check the glue / boundary.
"""
import sys, itertools, math
from fractions import Fraction as F
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7)}
INNER = set(range(1, 7))
CLAIMED = {7: F(283, 1470), 8: F(629, 2058), 9: F(16969, 41160),
           10: F(30551, 61740), 11: F(71111, 123480)}  # keyed by m


# ---------------------------------------------------------------------------
# Engine A: prompt's single_block_decorr (midpoint sampling at given Nx)
# ---------------------------------------------------------------------------
def single_block_decorr(m, Nx):
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


# ---------------------------------------------------------------------------
# Engine B: EXACT single-block cover via exact x-integration.
# For a fixed rational x the phi-cover is exact (piecewise const in phi).
# The cover-as-function-of-x is piecewise constant with breakpoints only where
# some frac((j-j')*x) or frac(j x) crosses a multiple of 1/7, i.e. x = a/(7d),
# d = |j - j'| in 1..m-1 (and d=j in 0..m-1). We integrate EXACTLY by splitting
# [0,1) at ALL such rational breakpoints and evaluating phi-cover at an interior
# rational point of each x-cell. This is fully exact (no sampling assumption).
# ---------------------------------------------------------------------------
def phi_cover_at_x(m, x):
    """Exact measure over phi in [0,1) that union_{j} {phi+frac(jx)} hits all 6 inner sectors."""
    r = [(j * x) % 1 for j in range(m)]
    bps = sorted({(F(s, 7) - rj) % 1 for rj in r for s in range(7)})
    bps.append(bps[0] + 1)
    good = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        hit = {int(((mid + rj) % 1) * 7) for rj in r}
        if len(hit & INNER) == 6:
            good += b - a
    return good


def exact_single_block(m):
    """Fully exact x-integral of phi_cover_at_x over x in [0,1)."""
    # x-breakpoints: all a/(7d) for d in 1..m-1, a in 0..7d.  (d=0 trivial.)
    xb = set()
    for d in range(1, m):
        for a in range(0, 7 * d + 1):
            xb.add(F(a, 7 * d))
    xb = sorted(b for b in xb if 0 <= b <= 1)
    if xb[0] != F(0):
        xb = [F(0)] + xb
    if xb[-1] != F(1):
        xb = xb + [F(1)]
    tot = F(0)
    for a, b in zip(xb, xb[1:]):
        if b == a:
            continue
        xmid = (a + b) / 2
        tot += phi_cover_at_x(m, xmid) * (b - a)
    return tot


# ---------------------------------------------------------------------------
# Engine C: EXACT inclusion-exclusion over missed inner sectors S subset {1..6}.
# p0_decorr = sum_{S} (-1)^{|S|} A_S(m),  A_S(m) = mean_x [ 1 - meas_phi( union_j (D_S - frac(jx)) ) ]
# where D_S = union of the |S| arcs [s/7,(s+1)/7) for s in S. We compute A_S exactly
# by exact x-integration too. This is an INDEPENDENT route (different algorithm).
# ---------------------------------------------------------------------------
def measure_union_arcs(intervals):
    """Exact Lebesgue measure of a union of [a,b) arcs on the circle [0,1), inputs in [0,1)."""
    # normalize, split wrap-around
    segs = []
    for a, b in intervals:
        a %= 1
        b %= 1
        if a < b:
            segs.append((a, b))
        elif a > b:
            segs.append((a, F(1)))
            segs.append((F(0), b))
        # a==b means full or empty; our arcs are length 1/7 so never equal
    segs.sort()
    tot = F(0)
    cur_a = None
    cur_b = None
    for a, b in segs:
        if cur_a is None:
            cur_a, cur_b = a, b
        elif a <= cur_b:
            if b > cur_b:
                cur_b = b
        else:
            tot += cur_b - cur_a
            cur_a, cur_b = a, b
    if cur_a is not None:
        tot += cur_b - cur_a
    return tot


def A_S_at_x(S, m, x):
    """meas_phi that ALL m points avoid every sector in S = 1 - meas(union_j (D_S - frac(jx)))."""
    r = [(j * x) % 1 for j in range(m)]
    intervals = []
    for s in S:
        lo = F(s, 7)
        hi = F(s + 1, 7)
        for rj in r:
            intervals.append(((lo - rj) % 1, (hi - rj) % 1))
    cov = measure_union_arcs(intervals)
    return F(1) - cov


def exact_A_S(S, m):
    xb = set()
    for d in range(1, m):
        for a in range(0, 7 * d + 1):
            xb.add(F(a, 7 * d))
    xb = sorted(b for b in xb if 0 <= b <= 1)
    if xb[0] != F(0):
        xb = [F(0)] + xb
    if xb[-1] != F(1):
        xb = xb + [F(1)]
    tot = F(0)
    for a, b in zip(xb, xb[1:]):
        if b == a:
            continue
        xmid = (a + b) / 2
        tot += A_S_at_x(S, m, xmid) * (b - a)
    return tot


def exact_IE(m):
    tot = F(0)
    for sz in range(1, 7):
        for S in itertools.combinations(range(1, 7), sz):
            tot += (-1) ** sz * exact_A_S(set(S), m)
    return tot


# ---------------------------------------------------------------------------
# DIRECT cover p0(E) for a FINITE speed set E (no decorrelation): the real LRC quantity.
# S7(E) = { x in [0,1): all 6 inner sectors hit by some frac(e_i x) }.
# Exact via x-breakpoints at a/(7*e_i) and a/(7*(e_i-e_j)).
# ---------------------------------------------------------------------------
def direct_p0(E, exact=True, Nx=200000):
    E = sorted(set(E))
    if exact:
        diffs = set()
        for e in E:
            diffs.add(e)
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
    else:
        cnt = 0
        for ix in range(Nx):
            x = (ix + 0.5) / Nx
            hit = {int((e * x % 1) * 7) for e in E}
            if len(hit & INNER) == 6:
                cnt += 1
        return F(cnt, Nx)


def is_primitive(E):
    g = 0
    for e in E:
        g = math.gcd(g, e)
    return g == 1


def span(E):
    return max(E) - min(E)


def main():
    print(__doc__)
    Nx_exact = {m: 7 * math.lcm(*range(1, m)) if m >= 2 else 7 for m in range(7, 12)}

    print("=" * 78)
    print("TASK 1: re-derive p0_decorr exactly; compare prompt-engine, exact-x, exact-IE, claim")
    print("=" * 78)
    print(f"{'m':>3} {'k':>3} {'prompt(1260)':>14} {'exact-x':>22} {'exact-IE':>22} {'claim':>16} {'all-agree':>10}")
    exactvals = {}
    for m in range(7, 12):
        k = m + 1
        pe = single_block_decorr(m, 1260)
        ex = exact_single_block(m)
        ie = exact_IE(m)
        cl = CLAIMED[m]
        exactvals[m] = ex
        agree = (ex == ie == cl)
        print(f"{m:>3} {k:>3} {str(pe):>14} {str(ex):>22} {str(ie):>22} {str(cl):>16} {str(agree):>10}")
    print()
    print("  -> prompt-engine(1260) vs exact:")
    for m in range(7, 12):
        pe = single_block_decorr(m, 1260)
        ex = exactvals[m]
        print(f"     m={m}: prompt={float(pe):.8f}  exact={float(ex):.8f}  "
              f"equal={pe==ex}  prompt-was-{'EXACT' if pe==ex else 'APPROX'}")

    print()
    print("=" * 78)
    print("TASK 1b: margins p0_decorr_exact vs cap_k")
    print("=" * 78)
    margins = []
    for m in range(7, 12):
        k = m + 1
        ex = exactvals[m]
        cap = CAPS[k]
        mar = cap - ex
        margins.append((mar, k))
        print(f"  k={k}: p0_decorr={str(ex):>20} = {float(ex):.6f}  cap={float(cap):.6f}  "
              f"margin={str(mar):>20} = {float(mar):.6f}  {'OK' if ex < cap else 'EXCEEDS!!'}")
    mmin = min(margins)
    print(f"  MIN margin = {str(mmin[0])} = {float(mmin[0]):.6f}  at k={mmin[1]}  "
          f"claimed 111019/588588={float(F(111019,588588)):.6f}  match={mmin[0]==F(111019,588588)}")

    print()
    print("=" * 78)
    print("TASK 3: single-arc closed form A_single(w,m) = (7-w)^2/(49(m-1)) for w in {4,5,6}")
    print("=" * 78)
    # A_single(w,m) = mean_x meas_phi(all m points avoid a single arc-set of w sectors).
    # Use exact_A_S with S = a contiguous block of w sectors (any block; necklace-invariant).
    for m in range(7, 12):
        for w in (4, 5, 6):
            S = set(range(1, 1 + w)) if w <= 6 else set(range(1, 7))
            # ensure within {1..6}
            S = set(list(range(1, 7))[:w])
            av = exact_A_S(S, m)
            cf = F((7 - w) ** 2, 49 * (m - 1))
            print(f"  m={m} w={w}: exact A_S={str(av):>14}={float(av):.6f}  "
                  f"closed=(7-w)^2/(49(m-1))={str(cf):>14}={float(cf):.6f}  match={av==cf}")
    print("  (note: closed form should be necklace-invariant; testing the contiguous block.)")
    # also test a NON-contiguous w-arc to verify the closed form is only for SINGLE arc (contiguous)
    print("  non-contiguous w=2 (sectors {1,3}) vs contiguous {1,2}: (should differ -> not single arc)")
    for m in (7, 8):
        cont = exact_A_S({1, 2}, m)
        noncont = exact_A_S({1, 3}, m)
        print(f"     m={m}: contiguous{{1,2}}={float(cont):.6f}  noncontiguous{{1,3}}={float(noncont):.6f}  "
              f"differ={cont!=noncont}")

    print()
    print("=" * 78)
    print("TASK 2: HUNT for a wide primitive k-set with p0(E) > cap_k (direct cover)")
    print("=" * 78)
    worst = {k: (F(0), None) for k in range(8, 13)}
    tested = 0
    families = []
    # consecutive blocks at various offsets (the pinned shape) but WIDE via scaling won't change p0 (scale-inv).
    # Build genuinely WIDE primitive sets:
    for k in range(8, 13):
        cands = []
        # AP with difference d (block of k terms): {0, d, 2d, ..., (k-1)d}+offset, primitive
        for d in range(1, 6):
            for off in range(0, 4):
                E = tuple(off + d * i for i in range(k))
                cands.append(E)
        # consecutive starting away from 0 (still primitive if contains coprime)
        for start in range(0, 25):
            E = tuple(start + i for i in range(k))
            cands.append(E)
        # two-scale clusters: small block + far block
        for far in (15, 19, 23, 30, 50, 100):
            for split in range(2, k - 1):
                E = tuple(range(split)) + tuple(far + i for i in range(k - split))
                cands.append(E)
        # resonant: multiples of 7-ish, and near-7 differences
        for E in [tuple(range(k)), tuple(7 * i + (i % 2) for i in range(k))]:
            cands.append(E)
        # boundary span 15-30 consecutive-ish with a hole
        for s in range(15, 31):
            E = tuple([0] + [i for i in range(1, k)] ) # placeholder
        # spread sets with one far speed
        for far in range(15, 31):
            E = tuple(range(k - 1)) + (far,)
            cands.append(E)
        for E in cands:
            E = tuple(sorted(set(E)))
            if len(E) != k:
                continue
            if not is_primitive(E):
                continue
            if span(E) <= 14:
                continue  # wide only
            p0 = direct_p0(E, exact=True)
            tested += 1
            if p0 > worst[k][0]:
                worst[k] = (p0, E)
            if p0 > CAPS[k]:
                print(f"  !!! COUNTEREXAMPLE k={k}: E={E} span={span(E)} p0={float(p0):.6f} > cap={float(CAPS[k]):.6f}")
    print(f"  tested {tested} wide primitive sets")
    for k in range(8, 13):
        p0, E = worst[k]
        cap = CAPS[k]
        print(f"  k={k}: max wide p0={float(p0):.6f} (cap={float(cap):.6f}, margin {float(cap-p0):.6f}) at E={E} span={span(E) if E else '-'}")

    print()
    print("=" * 78)
    print("TASK 2b: decorrelation error = direct_p0(wide consec) vs single-block-decorr limit")
    print("=" * 78)
    for k in range(8, 13):
        m = k - 1
        ex = exactvals[m]
        # finite wide consecutive block far out
        for M in (19, 50, 200):
            E = tuple(range(M, M + k))
            if not is_primitive(E):
                # make primitive by ensuring gcd 1 (consecutive integers always gcd 1 if k>=2)
                pass
            p0 = direct_p0(E, exact=True)
            err = p0 - ex
            print(f"  k={k} M={M}: direct p0(consec)={float(p0):.6f}  decorr_limit={float(ex):.6f}  "
                  f"error={float(err):+.6f}  cap={float(CAPS[k]):.6f}  "
                  f"{'STILL OK' if p0 < CAPS[k] else 'EXCEEDS CAP!!'}")

    print()
    print("=" * 78)
    print("TASK 2c: split vs single-block (decorrelated) -- does any split exceed single block?")
    print("=" * 78)
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
    for k in (9, 10, 11):
        m = k - 1
        sb = float(exactvals[m])
        print(f"  k={k}: single-block (exact) = {sb:.5f}")
        # enumerate partitions into blocks of sizes summing to m
        parts = []
        # generate integer partitions of m with parts>=1
        def gen(rem, mx, cur):
            if rem == 0:
                parts.append(list(cur))
                return
            for p in range(min(rem, mx), 0, -1):
                gen(rem - p, p, cur + [p])
        gen(m, m, [])
        any_exceed = False
        for P in parts:
            if len(P) == 1:
                continue
            blocks = [tuple(range(sz)) for sz in P]
            v = multiblock_decorr(blocks, 120, 18)
            if v > sb + 0.01:
                any_exceed = True
                print(f"      split {P}: {v:.5f}  >>> EXCEEDS single block !!!")
            elif v > sb - 0.005:
                print(f"      split {P}: {v:.5f}  (near single block)")
        if not any_exceed:
            print(f"      -> no split exceeded single block (within grid resolution)")

    print()
    print("=" * 78)
    print("VERDICT SUMMARY printed above. See exact comparisons in TASK 1.")
    print("=" * 78)


if __name__ == "__main__":
    main()
