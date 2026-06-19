#!/usr/bin/env python3
"""
ADVERSARIAL VERIFICATION of the LRC(14) "k8-direct-net-bound" claimed advance.
kind-pasteur-2026-06-18-S6-wf.

Goal: verify the claims around the 1/7-spread bound (the last lemma of LRC(14)):
  For every integer co-offset set E (0 in E, |E|=k, 8<=k<=12),
  prove  mu_{1/7}(E) >= thr_k,  equiv  meas(N(E)) <= cap_k = 1 - thr_k.

We:
  (A) re-implement the EXACT mu_theta engine and SANITY-check it,
  (B) recompute thr_k and the consec_k mu values claimed,
  (C) HUNT for a counterexample E with mu_{1/7}(E) < thr_k (small + large spread),
  (D) check "consecutive minimizes mu_{1/7}" claims, perforated/structured E,
  (E) sanity-check the B_7 minorant & its iid floor formula,
  (F) check the SCALING / REFLECTION theorems, sum-g^2 relaxation breakage at k=10.

Everything in EXACT rationals (fractions.Fraction).
"""
from fractions import Fraction as F
from itertools import combinations
import random
import sys

# ---------------------------------------------------------------------------
# EXACT mu_theta engine (copied verbatim from the prompt spec).
# ---------------------------------------------------------------------------
def mu_theta(E, theta):
    E = sorted(set(E)); n = len(E); bp = set([F(0), F(1)])
    for i in range(n):
        for j in range(i+1, n):
            d = E[j]-E[i]
            for m in range(0, d+1): bp.add(F(m, d))
    bp = sorted(b for b in bp if 0 <= b <= 1); total = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        mid = (a+b)/2; order = sorted(range(n), key=lambda i: (E[i]*mid) % 1)
        ks = [(E[order[t]]*mid).__floor__() for t in range(n)]; subs = []
        for t in range(n):
            o1 = order[t]; o2 = order[(t+1) % n]; k1 = ks[t]; k2 = ks[(t+1) % n]
            wrap = 1 if t == n-1 else 0
            s = E[o2]-E[o1]; c = F(k1-k2+wrap)
            if s == 0:
                if c > theta: subs.append((a, b))
            elif s > 0:
                lo = max(a, (theta-c)/s); subs.append((lo, b)) if lo < b else None
            else:
                hi = min(b, (theta-c)/s); subs.append((a, hi)) if a < hi else None
        subs.sort(); cur = cb = None
        for lo, hi in subs:
            if cur is None: cur, cb = lo, hi
            elif lo <= cb: cb = max(cb, hi)
            else: total += cb-cur; cur, cb = lo, hi
        if cur is not None: total += cb-cur
    return total


# ---------------------------------------------------------------------------
# INDEPENDENT brute-force mu_theta via fine sampling at exact rationals,
# to cross-check the engine on small examples.
# meas{x in [0,1): maxgap of {frac(e_i x)} > theta}.
# We compute the exact gap function and integrate over a refined grid of
# breakpoints, evaluating maxgap at midpoints. Use the SAME breakpoint set as
# engine but compute maxgap directly (no LP). This is an independent code path.
# ---------------------------------------------------------------------------
def maxgap_at(E, x):
    """Exact maxgap of the cyclic point set {frac(e_i x)} on the circle."""
    pts = sorted(set((e * x) % 1 for e in E))
    if len(pts) == 1:
        return F(1)  # single point: gap is full circle
    gaps = []
    for i in range(len(pts)):
        nxt = pts[(i+1) % len(pts)]
        if i+1 < len(pts):
            gaps.append(nxt - pts[i])
        else:
            gaps.append((1 - pts[i]) + pts[0])
    return max(gaps)

def mu_theta_bruteforce(E, theta):
    """Independent computation of meas{maxgap > theta} using the same
    breakpoint partition but a direct maxgap evaluation at midpoints.
    Within each cell the *order* and integer parts are constant, but maxgap
    is piecewise-LINEAR (not constant), so a midpoint evaluation only tells us
    whether maxgap>theta in MOST of the cell -- it can be wrong near a cell
    boundary where maxgap crosses theta. So this is only a SANITY proxy that
    should agree to within the (small) number of crossing cells. We refine by
    subdividing each cell and summing midpoint indicators; with enough
    subdivisions it converges to the true measure. We use rational midpoints."""
    E = sorted(set(E)); n = len(E); bp = set([F(0), F(1)])
    for i in range(n):
        for j in range(i+1, n):
            d = E[j]-E[i]
            for m in range(0, d+1): bp.add(F(m, d))
    bp = sorted(b for b in bp if 0 <= b <= 1)
    total = F(0)
    SUB = 8  # subdivisions per cell for the proxy
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        w = (b - a)
        for s in range(SUB):
            lo = a + w * s // 1 / SUB if False else a + w * F(s, SUB)
            hi = a + w * F(s+1, SUB)
            mid = (lo + hi) / 2
            if maxgap_at(E, mid) > theta:
                total += (hi - lo)
    return total


# ---------------------------------------------------------------------------
# (A) SANITY CHECK the engine against the brute-force proxy on small sets.
# ---------------------------------------------------------------------------
def section_A():
    print("=== (A) Engine sanity vs brute-force proxy ===")
    tests = [
        ([0,1,2,3,4,5,6,7], F(1,7)),
        ([0,1,2,3,4,5,6], F(1,7)),
        ([0,1,2,3,4,5,6,7,8], F(1,7)),
        ([0,2,3,5,7,11,13,14], F(1,7)),
        ([0,1,2,3,4,5,6,7], F(2,7)),
        ([0,1,3,4,5,9], F(1,7)),
    ]
    for E, th in tests:
        m_eng = mu_theta(E, th)
        m_bf  = mu_theta_bruteforce(E, th)
        diff = abs(m_eng - m_bf)
        flag = "OK" if diff <= F(1,50) else "**MISMATCH**"
        print(f"  E={E} th={th}: engine={m_eng} ({float(m_eng):.5f}) "
              f"bf~={float(m_bf):.5f} diff={float(diff):.5f} {flag}")
    # The key verified anchors from the prompt:
    print("  ANCHOR consec_8 mu_1/7 =", mu_theta(list(range(8)), F(1,7)),
          "expect 691/735 =", F(691,735))
    print("  ANCHOR consec_7 mu_1/7 =", mu_theta(list(range(7)), F(1,7)),
          "expect 1")
    print()


# ---------------------------------------------------------------------------
# (B) Recompute thr_k and consec_k mu values.
# ---------------------------------------------------------------------------
# G_P measure machinery. P is a set of small integers (the "small" part of S).
# G_P = {x : ||p x|| >= 1/14 for all p in P}.  meas(G_P) is exact.
def meas_GP(P):
    """meas{x in [0,1): frac dist of p*x to nearest integer >= 1/14 for all p in P}."""
    P = sorted(set(p for p in P if p != 0))
    if not P:
        return F(1)
    # ||p x|| >= 1/14  <=> there exists no integer within 1/14 of p x.
    # frac(p x) in [1/14, 13/14] mod 1 ... more precisely ||t|| = min(frac(t),1-frac(t)).
    # ||p x|| >= 1/14  <=>  frac(p x) in [1/14, 13/14].
    # Build breakpoints from all p.
    bp = set([F(0), F(1)])
    for p in P:
        for m in range(0, p+1):
            bp.add(F(m, p))
            # also the 1/14 and 13/14 offsets: frac(px)=1/14 => x=(m+1/14)/p
            bp.add(F(14*m + 1, 14*p))
            bp.add(F(14*m + 13, 14*p))
    bp = sorted(b for b in bp if 0 <= b <= 1)
    total = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        mid = (a+b)/2
        ok = True
        for p in P:
            t = (p*mid) % 1
            nd = min(t, 1-t)
            if nd < F(1,14):
                ok = False; break
        if ok:
            total += (b - a)
    return total

def section_B():
    print("=== (B) thr_k and consec_k mu values ===")
    # thr_k = 1 - min_{|P|=13-k} meas(G_P).  We need the min over all primitive
    # covering small-sets P with |P| = 13-k. The prompt GIVES the thr_k values
    # and the consec mu values; we recompute the consec mu values directly and
    # check thr_k claims where feasible.
    consec_expected = {
        8: F(691,735), 9: F(247,294), 10: F(38,49),
        11: F(1381,2205), 12: F(13823,24255), 13: F(477,1078),
    }
    thr_expected = {
        8: F(3637,5880), 9: None, 10: None, 11: None,
        12: F(1,7), 13: F(0),
    }
    for k in range(8, 14):
        mc = mu_theta(list(range(k)), F(1,7))
        exp = consec_expected[k]
        flag = "OK" if mc == exp else "**MISMATCH**"
        print(f"  k={k}: mu_1/7(consec_{k}) = {mc} ({float(mc):.5f}) "
              f"expect {exp} ({float(exp):.5f}) {flag}")
        te = thr_expected[k]
        if te is not None:
            margin = mc - te
            print(f"        thr_{k}={te} ({float(te):.5f})  "
                  f"consec - thr = {margin} ({float(margin):.5f}) "
                  f"{'>=0 OK' if margin>=0 else '**BELOW thr**'}")
    # cap_8 = 1 - thr_8
    cap8 = 1 - thr_expected[8]
    print(f"  cap_8 = 1-thr_8 = {cap8} ({float(cap8):.5f}) expect 2243/5880={F(2243,5880)}")
    print()
    return thr_expected, consec_expected


# ---------------------------------------------------------------------------
# (C) HUNT for a counterexample: integer E, 0 in E, 8<=k<=12,
#     with mu_{1/7}(E) < thr_k.  thr_k for k>=9 not given exactly; we use the
#     SUFFICIENT target form: it suffices that consecutive minimizes mu_1/7.
#     So we hunt for ANY E with mu_1/7(E) < mu_1/7(consec_k). That would break
#     the "consecutive minimizes" sufficient condition. We ALSO separately test
#     against the hard thr_8 for k=8 (the binding, exactly-known case).
# ---------------------------------------------------------------------------
def section_C(thr_expected, consec_expected):
    print("=== (C) Counterexample hunt: mu_1/7(E) < consec / thr ===")
    th = F(1,7)
    violations = []

    # ---- C1: exhaustive small spread, k=8, over subsets of {0..MAX} with 0,MAX in E.
    print("  -- C1: exhaustive k=8, primitive-ish, spread up to 14 --")
    thr8 = thr_expected[8]
    consec8 = consec_expected[8]
    best = (None, consec8)
    count = 0
    for MAX in range(7, 15):  # max element (spread). consec_8 has MAX=7.
        # choose 6 interior elements from 1..MAX-1, plus 0 and MAX
        interior = list(range(1, MAX))
        for combo in combinations(interior, 6):
            E = (0,) + combo + (MAX,)
            count += 1
            m = mu_theta(list(E), th)
            if m < best[1]:
                best = (E, m)
            if m < thr8:
                violations.append(("k8-hardthr", E, m))
            if m < consec8:
                violations.append(("k8-below-consec", E, m))
    print(f"     tested {count} sets; min mu found = {float(best[1]):.5f} at E={best[0]}")
    print(f"     consec_8 = {float(consec8):.5f}, thr_8 = {float(thr8):.5f}")
    print(f"     min >= consec_8 ? {best[1] >= consec8}")

    # ---- C2: larger spread random + structured descent for k=8 (the mu_2/7 crusher style)
    print("  -- C2: aggressive large-spread descent, k=8 --")
    random.seed(1234)
    best2 = (best[0], best[1])
    def try_E(E):
        nonlocal best2
        Es = tuple(sorted(set(E)))
        if len(Es) != 8 or Es[0] != 0: return None
        m = mu_theta(list(Es), th)
        if m < best2[1]:
            best2 = (Es, m)
        if m < consec8:
            violations.append(("k8-largespread-below-consec", Es, m))
        if m < thr8:
            violations.append(("k8-largespread-below-thr", Es, m))
        return m
    # random restarts with local descent (add/replace one coordinate to reduce mu)
    for restart in range(40):
        W = random.choice([14, 20, 30, 50, 80, 150, 300])
        E = sorted(random.sample(range(1, W+1), 7))
        E = [0] + E
        E[-1] = W  # ensure spread
        E = sorted(set(E))
        if len(E) < 8:
            continue
        cur = mu_theta(E, th)
        improved = True
        steps = 0
        while improved and steps < 60:
            improved = False
            steps += 1
            # try moving each interior coordinate to a nearby value
            for idx in range(1, len(E)-1):
                base = E[idx]
                for delta in (-3,-2,-1,1,2,3, -7,7, -W//7, W//7):
                    cand = base + delta
                    if cand <= E[0] or cand >= E[-1]: continue
                    newE = sorted(set(E[:idx] + [cand] + E[idx+1:]))
                    if len(newE) != 8: continue
                    mv = mu_theta(newE, th)
                    try_E(newE)
                    if mv < cur:
                        cur = mv; E = newE; improved = True; break
                if improved: break
    print(f"     best (lowest mu) found in descent: {float(best2[1]):.5f} at {best2[0]}")
    print(f"     still >= consec_8 ({float(consec8):.5f})? {best2[1] >= consec8}")

    # ---- C3: structured families known to stress nets: APs with gaps, near-AP, geometric
    print("  -- C3: structured families (perforated, geometric, Sidon-ish), k=8..12 --")
    fams = []
    # perforated AP: {0,1,..,k} minus some, scaled
    for k in range(8, 13):
        fams.append(("consec", list(range(k))))
        fams.append(("perf-skip2", [2*i for i in range(k)]))
        fams.append(("perf-odd", [0]+[2*i+1 for i in range(k-1)]))
        fams.append(("geom", sorted(set([0]+[ (3**i) % (5*k) for i in range(1,k)]))))
        # near-consecutive with one big outlier
        fams.append(("consec+outlier", list(range(k-1))+[5*k]))
        # Sidon-like (Mian-Chowla prefix)
        mc = [0,1,3,7,12,20,30,44,65,80,96,122,147,181,203,251,289]
        fams.append(("sidon", mc[:k]))
    for name, E in fams:
        E = sorted(set(E))
        if len(E) < 8 or E[0] != 0:
            continue
        k = len(E)
        m = mu_theta(E, th)
        cc = mu_theta(list(range(k)), th)
        rel = "consec is min" if m >= cc else "**BELOW consec**"
        if m < cc:
            violations.append((f"struct-{name}-k{k}", tuple(E), m))
        print(f"     {name:18s} k={k} E={E[:6]}{'...' if len(E)>6 else ''}: "
              f"mu={float(m):.5f} consec_k={float(cc):.5f} {rel}")
    print()
    return violations, best2


# ---------------------------------------------------------------------------
# (D) Test "consecutive minimizes mu_1/7" more exhaustively for k=8 and k=9.
#     Also test the FALSE-claims: stretch-monotonicity, single-difference.
# ---------------------------------------------------------------------------
def section_D():
    print("=== (D) 'consecutive minimizes mu_1/7' exhaustive-ish + monotonicity ===")
    th = F(1,7)
    # k=8 exhaustive over spread<=10 already done in C; here push k=8 to spread<=18
    # but that is huge; instead do a targeted exhaustive over spread 7..12 fully
    # and report the global min and whether consec_8 is it.
    consec8 = mu_theta(list(range(8)), th)
    gmin = (tuple(range(8)), consec8)
    cnt = 0
    for MAX in range(7, 13):
        for combo in combinations(range(1, MAX), 6):
            E = (0,) + combo + (MAX,)
            cnt += 1
            m = mu_theta(list(E), th)
            if m < gmin[1]:
                gmin = (E, m)
    print(f"  k=8 spread<=12: tested {cnt}, global min mu = {float(gmin[1]):.6f} "
          f"at {gmin[0]} (consec={float(consec8):.6f})")
    print(f"  consecutive IS the minimizer over this range? {gmin[1] >= consec8}")

    # k=9 spread<=11 exhaustive
    consec9 = mu_theta(list(range(9)), th)
    gmin9 = (tuple(range(9)), consec9)
    cnt9 = 0
    for MAX in range(8, 12):
        for combo in combinations(range(1, MAX), 7):
            E = (0,) + combo + (MAX,)
            cnt9 += 1
            m = mu_theta(list(E), th)
            if m < gmin9[1]:
                gmin9 = (E, m)
    print(f"  k=9 spread<=11: tested {cnt9}, global min mu = {float(gmin9[1]):.6f} "
          f"at {gmin9[0]} (consec={float(consec9):.6f})")
    print(f"  consecutive IS the minimizer over this range? {gmin9[1] >= consec9}")

    # stretch-monotonicity claimed FALSE: find E and a single "stretch" that DECREASES mu
    print("  -- stretch-monotonicity check (prompt claims FALSE) --")
    found_decrease = False
    random.seed(7)
    for _ in range(2000):
        k = random.randint(8, 10)
        W = random.randint(k, 25)
        E = sorted(random.sample(range(1, W+1), k-1)); E = [0]+E
        E = sorted(set(E))
        if len(E) != k: continue
        m0 = mu_theta(E, th)
        # stretch the top element up by 1
        E2 = E[:-1] + [E[-1]+1]
        E2 = sorted(set(E2))
        if len(E2) != k: continue
        m1 = mu_theta(E2, th)
        if m1 < m0:
            found_decrease = True
            print(f"     stretch DECREASES mu: E={E} mu={float(m0):.5f} -> "
                  f"E2={E2} mu={float(m1):.5f}  (confirms NON-monotone)")
            break
    if not found_decrease:
        print("     (no stretch-decrease found in 2000 random trials; inconclusive)")
    print()


# ---------------------------------------------------------------------------
# (E) B_7 minorant + iid floor formula sanity.
#     B_7(E) = meas{ at least one of 7 fixed arcs [(2i+1)/14,(2i+3)/14) empty of frac(e x)}.
#     Claim: B_7(E) <= mu_1/7(E) (rigorous lower bound on mu). And iid floor:
#     B_7^iid(k) = 1 - sum_{j=0}^{7} (-1)^j C(7,j) (1-j/7)^k.
# ---------------------------------------------------------------------------
from math import comb
def B7_exact(E):
    """Exact meas of the union over the 7 arcs A_i=[(2i+1)/14,(2i+3)/14) (i=0..6)
    of {x : no frac(e x) lands in A_i}.  We compute via the order-cell breakpoints
    (the same partition makes each indicator piecewise-constant? No: membership of
    frac(e x) in a FIXED arc changes at x where frac(e x) hits arc endpoints).
    So add breakpoints x = (m + a)/e for arc endpoints a in {(2i+1)/14}."""
    E = sorted(set(e for e in E))
    arcs = [(F(2*i+1,14), F(2*i+3,14)) for i in range(7)]
    bp = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        for m in range(0, e+1):
            for (lo,hi) in arcs:
                bp.add(F(m + lo*1, 1)/e if False else (F(m)+lo)/e)
                bp.add((F(m)+hi)/e)
            bp.add(F(m, e))
    bp = sorted(b for b in bp if 0 <= b <= 1)
    total = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        mid = (a+b)/2
        pts = [ (e*mid) % 1 for e in E ]
        # is there an arc with NO point inside?
        empty = False
        for (lo,hi) in arcs:
            has = any(lo <= p < hi for p in pts)
            if not has:
                empty = True; break
        if empty:
            total += (b - a)
    return total

def B7_iid(k):
    s = F(0)
    for j in range(0, 8):
        s += (-1)**j * comb(7, j) * F(7-j, 7)**k
    return 1 - s

def section_E():
    print("=== (E) B_7 minorant & iid floor ===")
    th = F(1,7)
    # check B_7(E) <= mu_1/7(E) on a batch
    print("  -- B_7(E) <= mu_1/7(E) sanity (rigorous-lower-bound claim) --")
    random.seed(99)
    ok = True
    samples = [list(range(8)), list(range(9)), [0,2,3,5,7,11,13,17],
               [0,1,2,3,5,8,13,21]]
    for _ in range(15):
        k = random.randint(8, 10)
        W = random.randint(k, 25)
        E = sorted(set([0]+random.sample(range(1,W+1), k-1)))
        if len(E) >= 8:
            samples.append(E)
    for E in samples:
        b7 = B7_exact(E)
        mu = mu_theta(E, th)
        flag = "OK" if b7 <= mu else "**B7 > mu (LOWER-BOUND VIOLATED)**"
        if b7 > mu: ok = False
        print(f"     E={E[:6]}{'...' if len(E)>6 else ''} k={len(E)}: "
              f"B7={float(b7):.5f} mu={float(mu):.5f} {flag}")
    print(f"  B_7 <= mu held on all samples? {ok}")
    # iid floor values
    print("  -- B_7^iid(k) floor values --")
    exp = {8:0.9755,9:0.9423,10:0.8951,11:0.8369,12:0.7715,13:0.7027}
    for k in range(8,14):
        v = B7_iid(k)
        print(f"     k={k}: B7_iid = {v} = {float(v):.5f}  (prompt {exp[k]})")
    print()


# ---------------------------------------------------------------------------
# (F) SCALING & REFLECTION theorems; sum-g^2 relaxation breakage at k=10.
# ---------------------------------------------------------------------------
def section_F():
    print("=== (F) Scaling, Reflection, sum-g^2 relaxation ===")
    # SCALING: mu_theta(cE) == mu_theta(E)
    print("  -- SCALING: mu_theta(cE)=mu_theta(E) --")
    bad = 0; tot = 0
    random.seed(5)
    for _ in range(40):
        k = random.randint(4, 9)
        W = random.randint(k, 20)
        E = sorted(set([0]+random.sample(range(1,W+1), k-1)))
        for th in (F(1,7), F(2,7), F(1,5), F(3,11)):
            base = mu_theta(E, th)
            for c in (2,3,5,7):
                cE = [c*e for e in E]
                v = mu_theta(cE, th)
                tot += 1
                if v != base:
                    bad += 1
                    if bad <= 3:
                        print(f"     **SCALING FAIL** E={E} c={c} th={th}: "
                              f"{v} != {base}")
    print(f"     scaling checks: {tot} total, {bad} violations")

    # REFLECTION: meas(N(E)) = meas(N(maxE - E)), where N = {maxgap<=1/7}.
    # meas(N(E)) = 1 - mu_1/7(E) (since mu = meas{maxgap>1/7}, and maxgap=1/7
    # exactly on a null set). So reflection <=> mu_1/7(E) = mu_1/7(maxE-E).
    print("  -- REFLECTION: mu_1/7(E) = mu_1/7(maxE - E) --")
    bad = 0; tot = 0
    th = F(1,7)
    random.seed(11)
    for _ in range(40):
        k = random.randint(4, 9)
        W = random.randint(k, 18)
        E = sorted(set([0]+random.sample(range(1,W+1), k-1)))
        mE = max(E)
        Er = sorted(set(mE - e for e in E))
        a = mu_theta(E, th); b = mu_theta(Er, th)
        tot += 1
        if a != b:
            bad += 1
            if bad <= 3:
                print(f"     **REFLECTION FAIL** E={E} -> {Er}: {a} != {b}")
    print(f"     reflection checks: {tot} total, {bad} violations")

    # sum-g^2 relaxation: N(E) subset {sum_t g_t^2 <= 1/7}? and meas of RHS.
    # g_t are the gaps of {frac(e x)} (cyclic). sum g_t = 1. maxgap<=1/7 => each
    # g<=1/7 => sum g^2 <= (1/7)*sum g = 1/7. So N(E) subset {sum g^2 <= 1/7}. valid.
    # Claim: meas{sum g^2 <= 1/7} at k=10 reaches 0.6197 > cap_10=0.6044 while
    # true net is only 0.0599. Verify the consec_10 numbers.
    print("  -- sum-g^2 relaxation: consec_10 net vs {sum g^2<=1/7} measure --")
    def net_consec(k):
        # net = meas{maxgap<=1/7} = 1 - mu_1/7(consec_k)
        return 1 - mu_theta(list(range(k)), F(1,7))
    def measSumSq(E, bound):
        """meas{x: sum of squared gaps of {frac(e x)} <= bound}. Gaps are
        piecewise-linear in x on order-cells, so sum g^2 is piecewise-quadratic.
        We approximate by fine rational sampling on the order-cell partition."""
        E = sorted(set(E)); n = len(E)
        bp = set([F(0), F(1)])
        for i in range(n):
            for j in range(i+1, n):
                d = E[j]-E[i]
                for m in range(0, d+1): bp.add(F(m, d))
        bp = sorted(b for b in bp if 0 <= b <= 1)
        total = F(0); SUB = 6
        for a, b in zip(bp, bp[1:]):
            if b <= a: continue
            w = b - a
            for s in range(SUB):
                lo = a + w*F(s,SUB); hi = a + w*F(s+1,SUB); mid=(lo+hi)/2
                pts = sorted(set((e*mid)%1 for e in E))
                if len(pts)==1:
                    ss = F(1)
                else:
                    gg = []
                    for i2 in range(len(pts)):
                        nxt = pts[(i2+1)%len(pts)]
                        gg.append(nxt-pts[i2] if i2+1<len(pts) else (1-pts[i2])+pts[0])
                    ss = sum(g*g for g in gg)
                if ss <= bound:
                    total += (hi-lo)
        return total
    for k in (8, 10):
        nc = net_consec(k)
        ms = measSumSq(list(range(k)), F(1,7))
        cap = {8: F(2243,5880), 10: None}[k]
        print(f"     k={k}: net(consec)={float(nc):.5f}  "
              f"meas{{sum g^2<=1/7}}~={float(ms):.5f}"
              + (f"  cap_{k}={float(cap):.5f}" if cap else ""))
    print()


if __name__ == "__main__":
    section_A()
    thr_expected, consec_expected = section_B()
    violations, best2 = section_C(thr_expected, consec_expected)
    section_D()
    section_E()
    section_F()
    print("=================================================================")
    if violations:
        print(f"!!! {len(violations)} VIOLATIONS FOUND:")
        for v in violations[:20]:
            print("   ", v[0], "E=", v[1], "mu=", float(v[2]))
    else:
        print("No counterexamples found in tested ranges "
              "(consec remained the minimizer; thr_8 never breached).")
