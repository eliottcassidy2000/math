# -*- coding: utf-8 -*-
# klein-2026-07-11-S252: THE {1..11,13,24} EDGE CASE -- exact B5 / clean-ruler scan.
#
# opus-S226: 71/72 two-element near-AP residuals have CLEAN rulers (live >= 1, maxBand <= 5)
# via pair-sums; the lone exception {1..11,13,24} has NO clean ruler q <= 400.  kps THM-707:
# EXACT B5 = liveCount - Penalty, Penalty = sum_p C(bandCount(p)-1, 5); B5 > 0 is the actual
# hB5 obligation (clean is the transparent sufficient case).  THIS SCRIPT decides the edge
# case exactly:
#   (a) first LIVE ruler (liveCount >= 1) and its structure;
#   (b) first B5 > 0 ruler (the exact-identity discharge -- what kps's Lean actually needs);
#   (c) first CLEAN ruler if any (extended range);
#   (d) the pair-sum diagnosis (why v_a + v_b fails here);
#   (e) the maxBand landscape (is maxBand <= 5 ever achievable?).
# Conventions == opus-S226 == kps THM-707: safe band = [ceil(q/14), floor(13q/14)] CLOSED;
# bandCount(p) = #{i : (v_i p) mod q outside}; live = bandCount 0.

from math import comb

V = list(range(1, 12)) + [13, 24]
assert len(V) == 13

def scan_q(v, q):
    lo = -(-q // 14); hi = (13 * q) // 14
    live = 0; maxband = 0; penalty = 0
    for p in range(1, q):
        bc = 0
        for vi in v:
            r = (vi * p) % q
            if not (lo <= r <= hi):
                bc += 1
        if bc == 0:
            live += 1
        elif bc >= 6:
            penalty += comb(bc - 1, 5)
        if bc > maxband:
            maxband = bc
    return live, maxband, penalty

def main():
    print(f"family v = {V}  (the opus-S226 edge case; max = 24 >= 23: residual class)")
    print("=" * 78)
    first_live = None; first_b5 = None; first_clean = None
    best = []  # (B5, q, live, maxband, penalty) top rulers
    QMAX = 2000
    for q in range(8, QMAX + 1):
        live, mb, pen = scan_q(V, q)
        b5 = live - pen
        if live >= 1 and first_live is None:
            first_live = (q, live, mb, pen, b5)
        if b5 > 0 and first_b5 is None:
            first_b5 = (q, live, mb, pen, b5)
        if live >= 1 and mb <= 5 and first_clean is None:
            first_clean = (q, live, mb, pen, b5)
        best.append((b5, q, live, mb, pen))
    best.sort(reverse=True)
    print(f"(a) first LIVE ruler:   {first_live}   (q, live, maxBand, penalty, B5)")
    print(f"(b) first B5>0 ruler:   {first_b5}")
    print(f"(c) first CLEAN ruler (q <= {QMAX}): {first_clean}")
    print(f"    top-5 by B5: {best[:5]}")
    mbs = {}
    for b5, q, live, mb, pen in best:
        mbs[mb] = mbs.get(mb, 0) + 1
    print(f"(e) maxBand histogram over q in [8, {QMAX}]: {dict(sorted(mbs.items()))}")
    print()
    print("(d) pair-sum diagnosis (all pair-sums q = v_a + v_b, deduped):")
    ps = sorted(set(a + b for i, a in enumerate(V) for b in V[i + 1:]))
    for q in ps:
        live, mb, pen = scan_q(V, q)
        b5 = live - pen
        tag = "CLEAN" if (live >= 1 and mb <= 5) else ("B5>0" if b5 > 0 else "")
        print(f"    q={q:>3}: live={live:>3} maxBand={mb:>2} penalty={pen:>5} B5={b5:>5}  {tag}")

if __name__ == "__main__":
    main()
