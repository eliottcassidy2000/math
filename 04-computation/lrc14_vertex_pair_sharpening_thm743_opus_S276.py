#!/usr/bin/env python3
"""
lrc14_vertex_pair_sharpening_thm743_opus_S276.py
================================================
opus-2026-07-13-S276.  THM-743: the VERTEX-TERM PAIR-SECTOR SHARPENING of THM-742.

THM-742's second-order constant C2 = 6 Sum j^2 + 4 max(J) V_R is dominated by the vertex term
(shape 1: 4*11*129 = 5676 of 8712).  That term charged EVERY breakpoint the WORST pair cost.
The perspective frame (the T(n-2) pair sector) says a vertex is a PAIR EVENT:

  * an arrangement vertex u* is a coincidence of the endpoints of arcs j and j' -- it exists
    only at delta*u* == 0 or +-1/7 (mod 1), delta = |j-j'|: the DIFFERENCE-RUNNER's grid
    (k +- 1/7)/delta -- the three-distance structure of J's difference set (klein-S294's
    windowed-overlap Farey resonances are the same pair-sector geometry);
  * the fiber interval that collapses there closes at the RELATIVE SPEED delta -- so the
    capping triangle in the event strip costs <= delta/(2W) on the strip-average (g) side;
  * the strand (f) side is SIGN-FAVORABLE for the lower bound: truth >= linear model
    (max(0,l) >= l), so annihilations only help L >= Area - bound;
  * events BURIED inside a third arc change nothing (the union is locally that third arc):
    only EXPOSED events cost.

THM-743.  In THM-742's bound the vertex part of C2 improves to
    two-sided:  Sum over exposed events of delta_v   (stated; derived: including the f-side)
    lower-bound-only (the LRC-relevant direction):  (1/2) Sum_v delta_v  (g-side only),
so C2 = 3 Sum_{j in J} j^2 + Sum_v^exposed delta_v  (stated)  replaces 6 Sum j^2 + 4 max(J) V_R.

This script: enumerates the actual events per breakpoint (which pairs coincide; exposure
against third arcs), reports the delta-histogram (perspective prediction: small-delta pairs --
adjacent cluster members -- are numerous but CHEAP), recomputes C2 and W0, and re-verifies the
full bound |L - Area| <= C1/W + C2'/W^2 over the W-spread in exact rationals (zero-violation
check), both THM-739 shapes.
"""
from fractions import Fraction as F

LAM = F(1, 14)

def normalize(ivs):
    ivs = sorted((a, b) for a, b in ivs if b > a)
    out = []
    for a, b in ivs:
        if out and a <= out[-1][1]:
            if b > out[-1][1]: out[-1] = (out[-1][0], b)
        else: out.append((a, b))
    return out

def intersect(A, B):
    out = []; i = j = 0
    while i < len(A) and j < len(B):
        a1, b1 = A[i]; a2, b2 = B[j]
        lo, hi = max(a1, a2), min(b1, b2)
        if hi > lo: out.append((lo, hi))
        if b1 <= b2: i += 1
        else: j += 1
    return out

def measure(A): return sum(b - a for a, b in A)

def safe_set(w):
    return [((k + LAM) / w, (k + 1 - LAM) / w) for k in range(w)]

def good_intervals(speeds):
    cur = [(F(0), F(1))]
    for w in sorted(speeds): cur = intersect(cur, safe_set(w))
    return cur

def A_J(u, J):
    ivs = []
    for j in J:
        c = (-j * u) % 1
        a, b = c - LAM, c + LAM
        if a < 0: ivs += [(a + 1, F(1)), (F(0), b)]
        elif b > 1: ivs += [(a, F(1)), (F(0), b - 1)]
        else: ivs.append((a, b))
    return 1 - measure(normalize(ivs))

def breakpoints(J, extra):
    pts = set(extra)
    diffs = sorted({abs(a - b) for a in J for b in J if a != b})
    for d in diffs:
        for k in range(d + 1):
            for off in (F(0), F(1, 7), -F(1, 7)):
                u = (F(k) + off) / d
                if 0 <= u <= 1: pts.add(u)
    pts.add(F(0)); pts.add(F(1))
    return sorted(pts)

def integrate_AJ(J, dom):
    pts = breakpoints(J, [e for iv in dom for e in iv])
    tot = F(0)
    def integ(a, b, fa, fb, d=0):
        m = (a + b) / 2; fm = A_J(m, J)
        assert fm * 2 == fa + fb or d <= 12
        if fm * 2 == fa + fb:
            return (fa + fb) * (b - a) / 2
        return integ(a, m, fa, fm, d + 1) + integ(m, b, fm, fb, d + 1)
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        for da, db in dom:
            lo, hi = max(a, da), min(b, db)
            if hi > lo: tot += integ(lo, hi, A_J(lo, J), A_J(hi, J))
    return tot, pts

def vertex_census(J, GB, pts):
    """enumerate actual pair events at breakpoints interior to G_B; exposure vs third arcs.
    returns (sum_delta_all, sum_delta_exposed, histogram{delta: (count_all, count_exposed)})"""
    interior = [p for p in pts if any(da < p < db for da, db in GB)]
    Jl = sorted(J)
    hist = {}
    s_all = 0; s_exp = 0
    for u in interior:
        for a in range(len(Jl)):
            for b in range(a + 1, len(Jl)):
                j1, j2 = Jl[a], Jl[b]
                d = j2 - j1
                r = (d * u) % 1
                if r == 0 or r == F(1, 7) or r == 1 - F(1, 7):
                    # meeting point(s): endpoints of arc j1 that coincide with arc j2's
                    if r == 0:
                        xs = [(-j1 * u + F(1, 14)) % 1, (-j1 * u - F(1, 14)) % 1]
                    elif r == F(1, 7):
                        xs = [(-j1 * u - F(1, 14)) % 1]   # left of j1 meets right of j2
                    else:
                        xs = [(-j1 * u + F(1, 14)) % 1]   # right of j1 meets left of j2
                    exposed = False
                    for x in xs:
                        buried = False
                        for j3 in Jl:
                            if j3 == j1 or j3 == j2: continue
                            dd = (x - (-j3 * u)) % 1
                            if (dd < LAM or dd > 1 - LAM) and dd != LAM and dd != 1 - LAM:
                                buried = True; break
                        if not buried:
                            exposed = True; break
                    cnt_a, cnt_e = hist.get(d, (0, 0))
                    hist[d] = (cnt_a + 1, cnt_e + (1 if exposed else 0))
                    s_all += d
                    if exposed: s_exp += d
    return s_all, s_exp, hist

# constants from THM-742 (S275, verified)
S275 = {
    "shape 1 {1}u{W..W+11}": dict(B=[1], J=list(range(12)), C1=F(1449, 100), oldC2=8712, oldW0=452, old739=1948),
    "shape 2 {1,2,3}u{W..W+9}": dict(B=[1, 2, 3], J=list(range(10)), C1=F(1914, 100), oldC2=4086, oldW0=584, old739=2676),
}

print("=" * 104)
print("THM-743: vertex-term pair-sector sharpening.  C2 vertex part: 4 max(J) V_R -> Sum_v^exposed delta_v")
print("=" * 104)

for name, cfg in S275.items():
    Bset, J = cfg["B"], cfg["J"]
    GB = good_intervals(Bset)
    area, pts = integrate_AJ(J, GB)
    s_all, s_exp, hist = vertex_census(J, GB, pts)
    sumj2 = sum(j * j for j in J); Jm = max(J)
    C2_stated = 3 * sumj2 + s_exp          # two-sided stated
    C2_derived = F(3, 2) * sumj2 + F(s_exp, 2)  # lower-bound-only derived
    print("\n%s:  Area = %.6f" % (name, float(area)))
    print("   events: Sum delta (all) = %d ; Sum delta (EXPOSED) = %d ; THM-742 charged 4*max*V_R = %d"
          % (s_all, s_exp, 4 * Jm * len([p for p in pts if any(da < p < db for da, db in GB)])))
    print("   delta-histogram (delta: all/exposed): %s"
          % ", ".join("%d: %d/%d" % (d, a, e) for d, (a, e) in sorted(hist.items())))
    print("   C2: %d (THM-742 stated) -> %d (THM-743 stated) -> %.1f (derived, lower-bound side)"
          % (cfg["oldC2"], C2_stated, float(C2_derived)))
    # re-verify the full bound with the new C2 over the W-spread
    C1 = cfg["C1"]
    print("   %6s %12s %14s %8s" % ("W", "|L-Area|", "new bound", "ok"))
    viol = 0
    for W in [10, 20, 30, 50, 90, 150, 250, 400, 800]:
        body = list(Bset) + [W + j for j in J]
        if len(set(body)) < 13: continue
        Lx = measure(good_intervals(body))
        err = abs(Lx - area)
        bnd = C1 / W + F(C2_stated, W * W)
        ok = err <= bnd
        viol += 0 if ok else 1
        print("   %6d %12.6f %14.6f %8s" % (W, float(err), float(bnd), ok))
    lo, hi = 1, 10 ** 7
    while lo < hi:
        mid = (lo + hi) // 2
        if C1 / mid + F(C2_stated, mid * mid) < area: hi = mid
        else: lo = mid + 1
    C1floor = 1
    lo2, hi2 = 1, 10 ** 7
    while lo2 < hi2:
        mid = (lo2 + hi2) // 2
        if C1 / mid < area: hi2 = mid
        else: lo2 = mid + 1
    print("   NEW W0 = %d  (THM-742: %d; crude THM-739: %d)  [C1-only floor: %d -- C1 now binds]"
          % (lo, cfg["oldW0"], cfg["old739"], lo2))
    print("   violations: %d" % viol)

print("\ndone.")
