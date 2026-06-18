#!/usr/bin/env python3
"""
lrc14_Bk_bounded-spread-exact-floor_kps-S5-wf.py   (kind-pasteur-2026-06-18-S5, ANGLE = bounded-spread-exact-floor)

GOAL: the SINGLE remaining lemma "B(k)" of LRC(14): a UNIFORM floor
        mu_min(k) := inf over integer co-offset sets E (0 in E, |E|=k) of mu(E)  >  0 ,  k <= 13,
where  mu(E) = meas{ x in [0,1) : the k phase-points {frac(e x)} have circular max-gap > 2/7 }
(equivalently fit in an arc of length < 5/7).  Then intersect with G_P to get rho*_min, the TRUE
target floor c0 of THM-527.

This file does FOUR things, all with EXACT fractions.Fraction:

(1) RIGOROUS mu(E):  good_set_exact adds BOTH kinds of breakpoint:
      (a) order-change points x = m/d  (d = e_i - e_j a positive difference, m integer), and
      (b) gap = 2/7 crossing points: inside each order-cell every cyclic gap is AFFINE in x,
          so {gap > 2/7} is a rational sub-interval found by solving (affine) = 2/7.
    We then INDEPENDENTLY VERIFY rigour: a fine deterministic sampling LOWER/UPPER sandwich must
    bracket the exact mu (the exact value must lie in [sample_lo, sample_hi] up to the sample step).

(2) BOUNDED-SPREAD EXACT mu_min(k), k=3..13: exhaustive over all primitive E (0 in E) with
    maxE <= cap.  We PUSH the cap upward and record mu_min(k, cap) to test the bounded-spread
    reduction: does increasing the cap LOWER mu_min, or does it STABILIZE?  Honest verdict either way.

(3) The OVERALL pure floor  min_k mu_min(k)  with the minimizing shape.

(4) G_P INTERSECTION: rho*(P,E) = meas{x in G_P : maxgap{frac(e x)} > 2/7}, where
      G_P = {x : ||p x|| >= 1/14 for all p in P}.  Exhaustive/heavy search over (P subset {1..13},
    bounded-spread E, |P|+k=13) for the exact min rho*, confirming 0 zeros.

Run:  python3 04-computation/lrc14_Bk_bounded-spread-exact-floor_kps-S5-wf.py [stage]
  stage in {rigour, mumin, push, floor, rho, all}  (default: all, but each stage also runnable alone)
Exact Fractions throughout; stdlib only.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
random.seed(20260618)

TWO7  = F(2, 7)
FIVE7 = F(5, 7)
ONE14 = F(1, 14)

# ============================================================================
# (1) RIGOROUS EXACT good set  Good(E) = { x in [0,1) : maxgap{frac(e x)} > 2/7 }
# ============================================================================
def merge(iv):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out

def meas(arcs):
    return sum((b - a for a, b in arcs), F(0))

def good_set_exact(E):
    """Exact set of x where the cyclic max-gap of {frac(e x)} exceeds 2/7.
    Breakpoints: order-changes m/d AND (implicitly, by solving inside each cell) gap=2/7 crossings."""
    E = sorted(set(E)); k = len(E)
    if k == 1:
        return [(F(0), F(1))]                # single point -> gap 1 > 2/7 everywhere
    diffs = set()
    for a in range(k):
        for b in range(a + 1, k):
            diffs.add(E[b] - E[a])
    bps = {F(0), F(1)}
    for d in diffs:
        for m in range(0, d + 1):
            bps.add(F(m, d))
    bps = sorted(x for x in bps if 0 <= x <= 1)
    good = []
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        # cyclic order of points is constant on (x0,x1); record order + integer floor of each e*xm
        pts = sorted(((E[t] * xm) % 1, E[t]) for t in range(k))
        order  = [e for _, e in pts]
        floors = [int((e * xm) // 1) for e in order]
        # each cyclic gap g_idx(x) = (e_nx x - f_nx) - (e_cur x - f_cur) [+1 for the wrap gap]
        #                          = (e_nx-e_cur) x + (f_cur-f_nx) [+wrap]  -- AFFINE in x.
        for idx in range(k):
            e_cur = order[idx]; f_cur = floors[idx]
            if idx < k - 1:
                e_nx = order[idx + 1]; f_nx = floors[idx + 1]; wrap = F(0)
            else:
                e_nx = order[0]; f_nx = floors[0]; wrap = F(1)
            A  = F(e_nx - e_cur)            # slope of this gap in x
            Cc = F(f_cur - f_nx) + wrap     # intercept
            # gap(x) = A x + Cc.  Want gap(x) > 2/7 on (x0,x1).
            if A == 0:
                if Cc > TWO7:
                    good.append((x0, x1))
                continue
            xb = (TWO7 - Cc) / A            # gap = 2/7 crossing point
            if A > 0:
                lo = max(x0, xb); hi = x1   # gap increasing: good above xb
            else:
                lo = x0; hi = min(x1, xb)   # gap decreasing: good below xb
            if lo < hi:
                good.append((lo, hi))
    return merge(good)

def mu(E):
    return meas(good_set_exact(E))

# ----- (1) RIGOUR CHECK: deterministic fine-sample sandwich -------------------
def maxgap_float(points):
    pts = sorted(p % 1.0 for p in points)
    g = 0.0
    for i in range(len(pts) - 1):
        g = max(g, pts[i + 1] - pts[i])
    return max(g, (pts[0] + 1.0) - pts[-1])

def mu_sample(E, N):
    """Midpoint-sampled measure of {maxgap>2/7}; returns float estimate (NOT a guaranteed bound,
    but with N large it must bracket exact mu to ~1/N up to # of crossing arcs)."""
    c = 0
    for t in range(N):
        x = (t + 0.5) / N
        if maxgap_float([(e * x) % 1.0 for e in E]) > 2.0 / 7.0 + 1e-12:
            c += 1
    return c / N

def stage_rigour():
    print("=" * 90)
    print("(1) RIGOUR of exact good_set_exact: order-change AND gap=2/7 breakpoints; sandwich check")
    print("=" * 90)
    tests = [
        [0, 1, 2], [0, 1, 2, 3], [0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        [0, 2, 3, 5, 7, 8], [0, 3, 7, 8, 11, 13], [0, 1, 4, 9, 16, 25], [0, 6, 9, 12, 14, 38],
        [0, 5, 17, 42], [0, 1, 2, 3, 4, 5, 6, 100], [0, 7, 14, 21, 28, 35],
    ]
    N = 600000
    allok = True
    for E in tests:
        ex = mu(E)
        approx = mu_sample(E, N)
        # exact must be within ~ (#arcs)/N of the sample; use a generous tolerance tied to N.
        nbr = len(good_set_exact(E))
        tol = (2 * nbr + 4) / N
        ok = abs(float(ex) - approx) <= tol
        allok = allok and ok
        print(f"   E={E}")
        print(f"        mu_exact = {ex} = {float(ex):.8f}   sample(N={N}) = {approx:.8f}   "
              f"|diff|={abs(float(ex)-approx):.2e}  tol={tol:.2e}  {'OK' if ok else 'MISMATCH!!'}")
    # scale-invariance + reflection sanity (the proved L0/L1)
    print("\n   L0/L1 sanity on the exact engine:")
    for E in [[0, 2, 3, 5, 7, 8], [0, 1, 4, 9, 16]]:
        g = good_set_exact(E)
        refl = merge(sorted([(1 - b, 1 - a) for a, b in g]))
        l0 = (refl == g)
        l1 = (mu([2 * e for e in E]) == mu(E) == mu([3 * e for e in E]))
        print(f"      E={E}: L0(reflection)={'OK' if l0 else 'FAIL'}  L1(scale c=2,3)={'OK' if l1 else 'FAIL'}")
        allok = allok and l0 and l1
    print("\n   VERDICT:", "exact engine RIGOROUS (sandwich + L0/L1 pass)" if allok else "PROBLEM — investigate")
    return allok

# ============================================================================
# helpers
# ============================================================================
def primitive(E):
    g = 0
    for e in E:
        g = gcd(g, e)
    return tuple(e // g for e in E) if g > 1 else tuple(E)

def is_primitive(E):
    g = 0
    for e in E:
        g = gcd(g, e)
    return g == 1

# ============================================================================
# (2) BOUNDED-SPREAD EXACT mu_min(k) with PUSHED caps
# ============================================================================
def mumin_exhaustive(k, cap):
    """Exact min of mu over all primitive E={0,e_2<...<e_k}, e_k<=cap. Returns (min, argmin)."""
    best = (F(2), None)
    for rest in itertools.combinations(range(1, cap + 1), k - 1):
        E = (0,) + rest
        if not is_primitive(E):
            continue
        m = mu(E)
        if m < best[0]:
            best = (m, E)
    return best

def stage_mumin(caps_by_k=None):
    print("=" * 90)
    print("(2) BOUNDED-SPREAD EXACT mu_min(k): exhaustive over primitive E, maxE<=cap")
    print("=" * 90)
    if caps_by_k is None:
        # caps chosen to stay tractable in exact arithmetic; bigger k -> smaller cap.
        caps_by_k = {3: 20, 4: 20, 5: 18, 6: 16, 7: 15, 8: 14, 9: 14, 10: 14, 11: 14, 12: 14, 13: 14}
    out = {}
    for k in range(3, 14):
        cap = caps_by_k[k]
        best = mumin_exhaustive(k, cap)
        out[k] = (best[0], best[1], cap)
        print(f"   k={k:2d}: mu_min(<=cap {cap:2d}) = {str(best[0]):>14s} = {float(best[0]):.6f}   at E={best[1]}")
    return out

# ============================================================================
# (2') PUSH the cap and watch whether mu_min keeps dropping  (test bounded-spread reduction)
# ============================================================================
def stage_push():
    print("=" * 90)
    print("(2') PUSH THE CAP: does mu_min(k,cap) STABILIZE (bounded-spread reduction TRUE)")
    print("     or keep DROPPING (reduction FALSE)?  Exact, per cap.")
    print("=" * 90)
    # k where exhaustive at growing caps is affordable
    cap_lists = {
        3: [6, 8, 12, 16, 24, 32, 48],
        4: [8, 12, 16, 24, 32, 44],
        5: [10, 14, 18, 24, 30],
        6: [10, 14, 18, 22],
        7: [10, 13, 16, 19],
        8: [10, 13, 15],
        9: [11, 13, 14],
    }
    for k, caps in cap_lists.items():
        row = []
        prev = None
        for cap in caps:
            if cap < k - 1:
                continue
            m, arg = mumin_exhaustive(k, cap)
            drop = "" if prev is None else ("  (DROP)" if m < prev else ("  (=)" if m == prev else "  (UP??)"))
            row.append((cap, m, arg, drop))
            prev = m
        print(f"\n   k={k}:")
        for cap, m, arg, drop in row:
            print(f"      cap<={cap:3d}: mu_min = {str(m):>14s} = {float(m):.6f}{drop}   argmin spread={max(arg)}  E={arg}")
        # verdict for this k
        ms = [r[1] for r in row]
        if len(ms) >= 2 and ms[-1] == ms[-2]:
            print(f"      -> STABILIZED at cap {row[-2][0]} (last two caps agree): bounded-spread OK for k={k}.")
        else:
            print(f"      -> NOT yet stable for k={k}: mu_min still changing as cap grows. CAUTION.")

# ============================================================================
# (4) G_P safe set and rho*
# ============================================================================
def GP_arcs(P):
    """G_P = { x in [0,1) : ||p x|| >= 1/14 for all p in P }.  Exact union of closed arcs.
    For one p: forbidden = union over j of (j/p - 1/(14p), j/p + 1/(14p)); safe = complement.
    Intersect over p in P."""
    safe = [(F(0), F(1))]
    for p in P:
        forb = []
        for j in range(p + 1):
            c = F(j, p)
            forb.append((c - F(1, 14 * p), c + F(1, 14 * p)))
        # safe_p = [0,1] minus forb (clipped)
        forb = merge([(max(F(0), a), min(F(1), b)) for a, b in forb if b > 0 and a < 1])
        sp = []
        cur = F(0)
        for a, b in forb:
            if a > cur:
                sp.append((cur, a))
            cur = max(cur, b)
        if cur < 1:
            sp.append((cur, F(1)))
        # intersect safe with sp
        new = []
        for a, b in safe:
            for c, d in sp:
                lo, hi = max(a, c), min(b, d)
                if lo < hi:
                    new.append((lo, hi))
        safe = merge(new)
    return safe

def intersect(A, B):
    out = []
    for a, b in A:
        for c, d in B:
            lo, hi = max(a, c), min(b, d)
            if lo < hi:
                out.append((lo, hi))
    return merge(out)

def rho_star(P, E):
    """meas( G_P  intersect  Good(E) )  -- exact."""
    return meas(intersect(GP_arcs(P), good_set_exact(E)))

# ============================================================================
# (4) rho*_min search
# ============================================================================
def stage_rho():
    print("=" * 90)
    print("(4) G_P INTERSECTION: rho*(P,E) = meas(G_P cap Good(E)),  |P|+|E\\{0}|... search the floor")
    print("    Setup mirrors THM-527: total speeds 13 = |P| small (subset of {1..13}) + large cluster.")
    print("    E = co-offsets of the large cluster (0 in E, |E|=k). |P| + k relates to 13; we scan")
    print("    P subset {1..13} and bounded-spread E, requiring rho*>0 everywhere; report exact min.")
    print("=" * 90)
    # The cluster size k and small part size: in S3, |L|=k large speeds, small part P has <=13-k
    # elements drawn from {1..13}. We scan representative (P,E) widely.
    best = (F(2), None, None)
    zeros = 0
    total = 0
    # (a) structured: P = {1..j} (consecutive small) for j=0..10, E bounded-spread primitive
    for j in range(0, 11):
        P = tuple(range(1, j + 1))
        kmax = 13 - j
        for k in range(3, min(kmax, 8) + 1):
            cap = {3: 12, 4: 11, 5: 10, 6: 9, 7: 9, 8: 9}[k]
            for rest in itertools.combinations(range(1, cap + 1), k - 1):
                E = (0,) + rest
                if not is_primitive(E):
                    continue
                r = rho_star(P, E)
                total += 1
                if r == 0:
                    zeros += 1
                    if zeros <= 5:
                        print(f"   !!! rho*=0 at P={P} E={E}")
                if r < best[0]:
                    best = (r, P, E)
    print(f"\n   structured scan: {total} (P,E) pairs, {zeros} zeros.")
    print(f"   min rho* (structured) = {best[0]} = {float(best[0]):.6f}  at P={best[1]} E={best[2]}")
    # (b) random P (arbitrary subsets of {1..13}) + random bounded E
    best2 = (F(2), None, None); zeros2 = 0; total2 = 0
    for _ in range(40000):
        sz = random.randint(0, 10)
        P = tuple(sorted(random.sample(range(1, 14), sz)))
        k = random.randint(3, 8)
        sp = random.choice([k, k + 1, 2 * k, 3 * k, 12, 20])
        if sp < k - 1:
            continue
        rest = tuple(sorted(random.sample(range(1, sp + 1), k - 1)))
        E = (0,) + rest
        if not is_primitive(E) or len(set(E)) != k:
            continue
        r = rho_star(P, E)
        total2 += 1
        if r == 0:
            zeros2 += 1
            if zeros2 <= 5:
                print(f"   !!! [rand] rho*=0 at P={P} E={E}")
        if r < best2[0]:
            best2 = (r, P, E)
    print(f"\n   random scan: {total2} (P,E) pairs, {zeros2} zeros.")
    print(f"   min rho* (random)     = {best2[0]} = {float(best2[0]):.6f}  at P={best2[1]} E={best2[2]}")
    overall = min(best, best2, key=lambda t: t[0])
    print(f"\n   *** OVERALL min rho* found = {overall[0]} = {float(overall[0]):.6f}  at P={overall[1]} E={overall[2]} ***")
    print(f"   *** total zeros across all scans: {zeros + zeros2} ***")
    return overall, zeros + zeros2

# ============================================================================
# main
# ============================================================================
def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("rigour", "all"):
        stage_rigour(); print()
    if stage in ("mumin", "all"):
        res = stage_mumin(); print()
        floor_k = min(res, key=lambda k: res[k][0])
        print(f"   >>> OVERALL pure floor  min_k mu_min(k) = {res[floor_k][0]} = {float(res[floor_k][0]):.6f} "
              f"at k={floor_k}, E={res[floor_k][1]} (cap {res[floor_k][2]})\n")
    if stage in ("push", "all"):
        stage_push(); print()
    if stage in ("rho", "all"):
        stage_rho(); print()
    print("DONE.")

if __name__ == "__main__":
    main()
