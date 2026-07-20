#!/usr/bin/env python3
"""
lrc14_q55_slice_census_boxeph_S133b.py  (HYP-7930 part B; owner-directed)

THE q=55 SLICE (4/55 analogue of the q=38 census).  Full exhaustion is
C(46,11) ~ 6.6e9 PER SCALE x 20 scales — out of reach here; delivered instead:

 (i) SAMPLED census, all 20 scales (deterministic seed): 13-sets
     {s1,s2} u (11 from A_m), A_m = {c in [1,54]: c*m mod 55 in [4,51]},
     s1 = 4*m^{-1} mod 55, s2 = 55-s1; covering{2..13}; rungs {17,19,23,25}
     (13 auto-blocked by the covering duty); survivor rate, blocker census,
     min sampled M (exact, prefiltered).
(ii) CAPPED EXHAUSTIVE WINDOW STATEMENT (strengthens S131's B-sweep semantics
     from "no 4/55 realizer" to "no WINDOW member"): enumerate ALL stack-passing
     13-sets with Vmax <= CAP over every scale; a leaf with some hole > 3/41
     has M > 3/41 (not in window); otherwise exact M decides: report any family
     with M in (1/14, 3/41) — expected NONE.

boxeph-2026-07-19-S133.  Pure Python, exact integers.
"""

import sys, random
from math import gcd

DUTIES = [7, 8, 9, 10, 11, 12, 13]
RUNGS = [17, 19, 23, 25]

def upairs(p):
    return [j for j in range(1, (p + 1) // 2) if gcd(j, p) == 1 and gcd(p - j, p) == 1]

def covering13(sp):
    return all(any(v % q == 0 for v in sp) for q in range(2, 14))

def rung_status(sp, p):
    if any(v % p == 0 for v in sp):
        return 'blk'
    hit = {min(v % p, p - v % p) for v in sp if gcd(v % p, p) == 1}
    return 'sprd' if all(j in hit for j in upairs(p)) else 'KILL'

def M_int(sp, cap2):
    bn, bd, qstar = 0, 1, None
    for q in range(2, cap2 + 1):
        for b in range(1, q // 2 + 1):
            md = q
            for c in sp:
                r = (c * b) % q
                d = r if r <= q - r else q - r
                if d < md:
                    md = d
                    if md * bd <= bn * q:
                        break
            if md * bd > bn * q:
                bn, bd, qstar = md, q, q
    g = gcd(bn, bd)
    return bn // g, bd // g, qstar

SCALES = [m for m in range(1, 28) if gcd(m, 55) == 1]

def sampled_census(per_scale):
    rng = random.Random(5533)
    print("=" * 96)
    print("(i) SAMPLED q=55 slice census: %d samples/scale x %d scales" % (per_scale, len(SCALES)))
    print("=" * 96)
    grand = {'cov': 0, 'surv': 0}
    blocker_hist = {}
    best = (1, 1, None, None, None)
    for m in SCALES:
        Am = [c for c in range(1, 55) if 4 <= min((c * m) % 55, 55 - (c * m) % 55)]
        s1 = (4 * pow(m, -1, 55)) % 55
        s2 = 55 - s1
        pool = [c for c in Am if c not in (s1, s2)]
        cov = surv = 0
        minm_scale = (1, 1)
        for _ in range(per_scale):
            C = sorted([s1, s2] + rng.sample(pool, 11))
            if not covering13(C):
                continue
            cov += 1
            sts = [rung_status(C, p) for p in RUNGS]
            if 'KILL' in sts:
                continue
            surv += 1
            blk = tuple(p for p, s in zip(RUNGS, sts) if s == 'blk')
            blocker_hist[blk] = blocker_hist.get(blk, 0) + 1
            # min-M hunt with cheap prefilter: skip if a small-q hole already
            # exceeds current global best (then M > best, useless for the min)
            bn, bd = best[0], best[1]
            deep = True
            for q in range(14, 34):
                mx = 0
                for b in range(1, q // 2 + 1):
                    md = min(min((c * b) % q, q - (c * b) % q) for c in C)
                    if md > mx: mx = md
                if mx * bd > bn * q:
                    deep = False; break
            if deep:
                nu, de, qs = M_int(list(C), 108)
                if nu * best[1] < best[0] * de:
                    best = (nu, de, qs, tuple(C), m)
                if nu * minm_scale[1] < minm_scale[0] * de:
                    minm_scale = (nu, de)
        grand['cov'] += cov; grand['surv'] += surv
        print("m=%2d straddle (%2d,%2d): covering=%5d survivors=%5d (%.1f%% of covering)"
              % (m, s1, s2, cov, surv, 100 * surv / max(1, cov)))
    print("\nTOTALS: covering=%d survivors=%d (%.2f%%)"
          % (grand['cov'], grand['surv'], 100 * grand['surv'] / max(1, grand['cov'])))
    print("blocker census (top 8):")
    for blk, c in sorted(blocker_hist.items(), key=lambda x: -x[1])[:8]:
        print("  %-22s: %d" % (str(blk) if blk else "(none: all-spread)", c))
    print("MIN SAMPLED M = %d/%d = %.4f  (scale m=%s, maximizer q=%s)"
          % (best[0], best[1], best[0] / best[1], best[4], best[2]))
    print("  family: %s" % (best[3],))
    print("  [window check: 1/14=%.4f < M? %s ; M < 3/41=%.4f? %s]"
          % (1 / 14, best[0] / best[1] > 1 / 14, 3 / 41, best[0] / best[1] < 3 / 41))

def capped_window_sweep(cap, leaf_limit=9_000_000):
    print("=" * 96)
    print("(ii) CAPPED EXHAUSTIVE WINDOW SWEEP: Vmax <= %d, all scales; target M in (1/14, 3/41)" % cap)
    print("=" * 96)
    reqbit = {q: 1 << k for k, q in enumerate(DUTIES)}
    BP, BM = 1 << len(DUTIES), 2 << len(DUTIES)
    REQ = (4 << len(DUTIES)) - 1
    off, width = {}, 0
    for p in RUNGS:
        off[p] = width; width += len(upairs(p)) + 1
    SEG = {p: (1 << len(upairs(p))) - 1 for p in RUNGS}
    leaves = crtpass = unbeaten = 0
    window_members, tights = [], []
    stack = []
    try:
        for m in SCALES:
            band = {r for r in range(55) if min(r, 55 - r) >= 4}
            A = [c for c in range(1, cap + 1) if (c * m) % 55 in band]
            if len(A) < 13: continue
            if any(all(c % q for c in A) for q in DUTIES): continue
            nA = len(A)
            met = [0] * nA; sbit = [0] * nA
            for i, c in enumerate(A):
                mm = 0
                for q in DUTIES:
                    if c % q == 0: mm |= reqbit[q]
                r = (c * m) % 55
                if r == 4: mm |= BP
                if r == 51: mm |= BM
                met[i] = mm
                sb = 0
                for p in RUNGS:
                    if c % p == 0:
                        sb |= 1 << (off[p] + len(upairs(p)))
                    elif gcd(c % p, p) == 1:
                        sb |= 1 << (off[p] + upairs(p).index(min(c % p, p - c % p)))
                sbit[i] = sb
            sufmet = [0] * (nA + 1); sufsb = [0] * (nA + 1)
            for i in range(nA - 1, -1, -1):
                sufmet[i] = sufmet[i + 1] | met[i]
                sufsb[i] = sufsb[i + 1] | sbit[i]

            def spread_ok(sb):
                for p in RUNGS:
                    seg = sb >> off[p]
                    L = len(upairs(p))
                    if not ((seg >> L) & 1 or (seg & SEG[p]) == SEG[p]):
                        return False
                return True

            def rec(start, got, sb):
                nonlocal leaves, crtpass, unbeaten
                slots = 13 - len(stack)
                if slots == 0:
                    if REQ & ~got: return
                    leaves += 1
                    if leaves > leaf_limit: raise KeyboardInterrupt
                    fam = stack
                    g0 = 0
                    for c in fam: g0 = gcd(g0, c)
                    if g0 != 1: return
                    if not spread_ok(sb): return
                    crtpass += 1
                    # beaten at the 3/41 level? (hole deeper than 3/41 => M > 3/41)
                    for q in range(2, 2 * cap + 1):
                        beat = False
                        for b in range(1, q // 2 + 1):
                            md = q
                            for c in fam:
                                r = (c * b) % q
                                d = r if r <= q - r else q - r
                                if d < md:
                                    md = d
                                    if md * 41 <= 3 * q: break
                            if md * 41 > 3 * q:
                                beat = True; break
                        if beat: break
                    else:
                        unbeaten += 1
                        nu, de, qs = M_int(list(fam), 2 * cap)
                        if nu * 14 > de:            # M > 1/14
                            window_members.append((nu, de, tuple(fam), m))
                        else:
                            tights.append((nu, de, tuple(fam), m))
                    return
                if (REQ & ~(got | sufmet[start])): return
                if not spread_ok(sb | sufsb[start]): return
                if nA - start < slots: return
                for i in range(start, nA):
                    stack.append(A[i])
                    rec(i + 1, got | met[i], sb | sbit[i])
                    stack.pop()

            rec(0, 0, 0)
        complete = True
    except KeyboardInterrupt:
        complete = False
    print("leaf visits=%d  stack-passing=%d  unbeaten-at-3/41=%d  complete=%s"
          % (leaves, crtpass, unbeaten, complete))
    print("WINDOW MEMBERS (M in (1/14,3/41)): %d  %s"
          % (len(set(window_members)), sorted(set(window_members))[:5]))
    print("M <= 1/14 leaves (tight-or-below controls): %d  %s"
          % (len(set(tights)), sorted(set(tights))[:5]))
    if complete and not window_members:
        print("=> NO family in the q=55 slice with Vmax <= %d has M in (1/14, 3/41)." % cap)

if __name__ == "__main__":
    if sys.argv[1] == "sample":
        sampled_census(int(sys.argv[2]))
    elif sys.argv[1] == "cap":
        capped_window_sweep(int(sys.argv[2]))
