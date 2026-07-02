#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S93 -- HYP-3840: the RADIUS-DERIVATIVE frame for the lonely measure.

Objects (n runners frame; here N=14 means 13 moving runners + observer, threshold 1/14):
  danger set   D_v(r) = { t in R/Z : ||v t|| < r }  (open; v intervals of length 2r/v)
  lonely set   L_S(r) = complement of union of D_v(r)  (closed)
  measure      m_S(r) = |L_S(r)|
  M(S)         = max_t min_v ||v t||   (covering-min of the set; LRC(N): M >= 1/N)

CLAIMS TESTED (all exact rational arithmetic):
 (1) m_S(r) is piecewise linear in r; ALL breakpoints lie on the Farey-type grid
     { d/(v+w) } u { d/(w-v) } for v,w in S (incl. v=w for d/(2v)), d in Z>0.
 (2) On each linear cell, slope(m_S) = - sum over components of L(r) of (1/v_left + 1/v_right)
     where v_left/v_right are the speeds owning the component's two boundary endpoints.
 (3) CONVEXITY: kinks at r=d/(v+w) (merge events) are convex; kinks at r=d/(w-v)
     (overtaking events) are concave.  If no overtaking event lies in (0, 1/N), m_S is
     convex there.  Overtaking below radius rho exists iff some pair has
     (w-v) > gcd(v,w)/rho, i.e. w-v > N*gcd(v,w) at rho=1/N.
 (4) AP closed form: for S={1..N-1}, near r=1/N:
        m_S(r) = C_AP(N) * (1 - N r),   C_AP(N) = (2/N) * sum_{a in (Z/N)*} 1/a.
     At N=14: C_AP = 2*(1/13+1/45+1/33) = 0.258891... (the S92 empirical "0.26").
 (5) UNIT-RESIDUE LEMMA: if M(S)=1/q exactly with a witness at a/q (lowest terms),
     then the residue set {a v mod q} contains both +1 and -1; hence every tight set
     (binding modulus q=N) represents EVERY unit residue; hence the dropped residue in
     the duplication+drop classification (HYP-3750) is a NON-UNIT; hence for q prime
     there are NO non-difference-closed tight sets (explains klein-S48 census zeros at
     n=4,6 i.e. q=5,7).
 (6) SLOPE (NON-)RIGIDITY: tight sets that duplicate a UNIT residue put a fast runner
     on a binding slot => slope contribution 1/(v w) with large v,w => total slope
     BELOW C_AP.  Test klein's {1,4,5,6,7,11,13} (q=8) and enumerate all dup+drop
     tight 13-sets mod 14 to find the TRUE inf slope at N=14 (S92 claimed AP=0.26 min).
"""
from fractions import Fraction as F
from math import gcd
import sys, itertools

def danger_intervals(S, r):
    """All open danger intervals (start,end,owner v) on [0,1), exact fractions."""
    iv = []
    for v in S:
        h = F(r, 1) if not isinstance(r, F) else r
        half = h / v
        for k in range(v):
            c = F(k, v)
            a, b = c - half, c + half
            iv.append((a, b, v))
    # wrap into [0,1)
    out = []
    for a, b, v in iv:
        if a < 0:
            out.append((a + 1, F(1), v)); out.append((F(0), b, v))
        elif b > 1:
            out.append((a, F(1), v)); out.append((F(0), b - 1, v))
        else:
            out.append((a, b, v))
    return out

def lonely_components(S, r):
    """Components of L(r) as list of (start, end, left_owner, right_owner). Exact."""
    iv = sorted(danger_intervals(S, r))
    # merge union, tracking owner of the current right frontier
    merged = []  # (a, b, owner_of_left_endpoint(a), owner_of_right_endpoint(b))
    for a, b, v in iv:
        if merged and a <= merged[-1][1]:
            pa, pb, la, ra = merged[-1]
            if b > pb:
                merged[-1] = (pa, b, la, v)
        else:
            merged.append((a, b, v, v))
    if not merged:
        return [(F(0), F(1), None, None)]
    comps = []
    for i in range(len(merged)):
        a, b, la, ra = merged[i]
        na, nb, nla, nra = merged[(i + 1) % len(merged)]
        gap_start, gap_end = b, (na if i + 1 < len(merged) else na + 1)
        if gap_end > gap_start:
            # component boundary owners: left endpoint of gap owned by ra (right end of
            # danger block), right endpoint owned by nla (left end of next block)
            comps.append((gap_start, gap_end, ra, nla))
    return comps

def lonely_measure(S, r):
    return sum(b - a for a, b, _, _ in lonely_components(S, r))

def slope_from_components(S, r):
    """Predicted d m / d r = - sum over components (1/v_left + 1/v_right)."""
    s = F(0)
    for a, b, lo, ro in lonely_components(S, r):
        if lo is None:
            return F(0)
        s += F(1, lo) + F(1, ro)
    return -s

def breakpoint_grid(S, rmax):
    """All candidate breakpoints d/(v+w), d/(w-v) in (0, rmax]."""
    pts = set()
    Sl = sorted(set(S))
    for i, v in enumerate(Sl):
        for w in Sl[i:]:
            for den in ([v + w] + ([w - v] if w > v else [])):
                if den <= 0:
                    continue
                d = 1
                while F(d, den) <= rmax:
                    pts.add(F(d, den)); d += 1
    return sorted(pts)

def M_exact(S, qcap_extra=0):
    """Exact covering-min via binding-pair candidate witnesses t=m/(v+w), m/(w-v), m/(2v)."""
    cands = set()
    Sl = sorted(set(S))
    dens = set()
    for i, v in enumerate(Sl):
        for w in Sl[i:]:
            dens.add(v + w)
            if w > v: dens.add(w - v)
    for den in dens:
        for m in range(1, den):
            cands.add(F(m, den))
    best = F(0); wit = None
    for t in cands:
        mn = min(dist_circle(v * t) for v in Sl)
        if mn > best:
            best, wit = mn, t
    return best, wit

def dist_circle(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def witnesses(S, M):
    """All binding-pair candidates achieving min ||vt|| == M exactly."""
    Sl = sorted(set(S)); dens = set()
    for i, v in enumerate(Sl):
        for w in Sl[i:]:
            dens.add(v + w)
            if w > v: dens.add(w - v)
    out = []
    for den in dens:
        for m in range(1, den):
            t = F(m, den)
            if min(dist_circle(v * t) for v in Sl) == M:
                out.append(t)
    return sorted(set(out))

def critical_slope(S, N=14):
    """Normalized slope c_S: m_S(r) = c_S*(1 - N r) on the last linear cell below 1/N.
    Returns (c_S, last_breakpoint_below)."""
    rmax = F(1, N)
    grid = [g for g in breakpoint_grid(S, rmax) if g < rmax]
    last = grid[-1] if grid else F(0)
    r1 = last + (rmax - last) * F(1, 3)
    r2 = last + (rmax - last) * F(2, 3)
    m1, m2 = lonely_measure(S, r1), lonely_measure(S, r2)
    # linear in cell: m = c*(1-N r) => c from two points (also cross-checks linearity)
    c1 = m1 / (1 - N * r1); c2 = m2 / (1 - N * r2)
    return c1, c2, last

def check_piecewise_linear_and_kinks(S, rmax, label):
    """Verify slopes constant per cell, kink signs classify by grid type."""
    grid = [g for g in breakpoint_grid(S, rmax) if 0 < g <= rmax]
    pts = [F(0)] + grid
    print(f"\n--- {label}: piecewise structure on (0, {rmax}] ---")
    print(f"S = {sorted(S)}; #grid breakpoint candidates = {len(grid)}")
    prev_slope = None; n_convex = n_concave = n_flat = 0; ok_slope_formula = True
    kink_report = []
    for i in range(len(pts) - 1):
        a, b = pts[i], pts[i + 1]
        r1 = a + (b - a) * F(1, 4); r2 = a + (b - a) * F(3, 4)
        m1, m2 = lonely_measure(S, r1), lonely_measure(S, r2)
        sl = (m2 - m1) / (r2 - r1)
        pred = slope_from_components(S, r1)
        if pred != sl:
            ok_slope_formula = False
        if prev_slope is not None:
            if sl > prev_slope: n_convex += 1
            elif sl < prev_slope:
                n_concave += 1
                kink_report.append((a, prev_slope, sl))
            else: n_flat += 1
        prev_slope = sl
    print(f"slope formula -sum(1/vL+1/vR) matches cell slope in ALL cells: {ok_slope_formula}")
    print(f"kinks: convex={n_convex}, CONCAVE={n_concave}, flat={n_flat}")
    if kink_report:
        print("concave kinks (r, slope_before, slope_after):")
        for r, s1, s2 in kink_report[:8]:
            print(f"   r={r} = {float(r):.6f}   {float(s1):.4f} -> {float(s2):.4f}")
    # overtaking prediction
    over = []
    Sl = sorted(set(S))
    for i, v in enumerate(Sl):
        for w in Sl[i + 1:]:
            g = gcd(v, w)
            if F(g, w - v) < rmax:
                over.append((v, w, F(g, w - v)))
    print(f"overtaking predicted below {rmax}: {len(over)} pair(s) "
          f"{[(v,w,str(r)) for v,w,r in over[:6]]}")
    print(f"=> convex on (0,{rmax}) predicted: {len(over)==0}; observed: {n_concave==0}")
    return n_concave == 0

# ---------------------------------------------------------------- experiments
def main():
    N = 14
    print("=" * 78)
    print("(1)+(2)+(3): piecewise linearity, slope formula, convexity criterion")
    print("=" * 78)
    AP = list(range(1, 14))
    check_piecewise_linear_and_kinks(AP, F(1, 14), "AP {1..13}")
    DW = list(range(1, 13)) + [182]
    check_piecewise_linear_and_kinks(DW, F(1, 14), "deep well {1..12,182}")
    # small mixed set with a spread pair to force overtaking
    SP = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 40]
    check_piecewise_linear_and_kinks(SP, F(1, 14), "spread {1..12,40}")

    print()
    print("=" * 78)
    print("(4): AP closed form  C_AP(N) = (2/N) sum_{a unit} 1/a")
    print("=" * 78)
    for NN, APn in [(14, list(range(1, 14))), (8, list(range(1, 8))), (6, list(range(1, 6)))]:
        cap = F(2, NN) * sum(F(1, a) for a in range(1, NN) if gcd(a, NN) == 1)
        c1, c2, last = critical_slope(APn, NN)
        print(f"N={NN}: predicted C_AP={cap} = {float(cap):.6f}; "
              f"measured c={c1}={float(c1):.6f} (2nd pt {float(c2):.6f}); last bp={last}")
        assert c1 == c2 == cap, f"AP closed form FAILS at N={NN}"
    print("AP closed form EXACT at N=6,8,14  [assert passed]")

    print()
    print("=" * 78)
    print("(5): unit-residue lemma sanity + (6): slope of known tight sets")
    print("=" * 78)
    tight_known = {
        "AP{1..5} q=6": (list(range(1, 6)), 6),
        "cross {1,3,4,5,9} q=6": ([1, 3, 4, 5, 9], 6),
        "AP{1..7} q=8": (list(range(1, 8)), 8),
        "GW {1,2,3,4,5,7,12} q=8": ([1, 2, 3, 4, 5, 7, 12], 8),
        "cross {1,4,5,6,7,11,13} q=8": ([1, 4, 5, 6, 7, 11, 13], 8),
    }
    for name, (S, q) in tight_known.items():
        M, wit = M_exact(S)
        cap = F(2, q) * sum(F(1, a) for a in range(1, q) if gcd(a, q) == 1)
        if M != F(1, q):
            print(f"{name}: NOT TIGHT, M={M}={float(M):.5f}")
            continue
        c1, c2, last = critical_slope(S, q)
        wits = witnesses(S, M)
        res = sorted(set(v % q for v in S))
        units = [a for a in range(1, q) if gcd(a, q) == 1]
        has_all_units = all(u in res for u in units)
        print(f"{name}: TIGHT M=1/{q}; slope c={c1}={float(c1):.6f} vs C_AP={float(cap):.6f} "
              f"({'RIGID' if c1 == cap else 'BEATS AP' if c1 < cap else 'ABOVE AP'}); "
              f"#witnesses={len(wits)}; residues={res}; all units present: {has_all_units}")
        assert c1 == c2

    print()
    print("=" * 78)
    print("(6): ENUMERATE dup+drop tight 13-sets mod 14 and find TRUE inf slope at N=14")
    print("=" * 78)
    q = 14
    units = [a for a in range(1, q) if gcd(a, q) == 1]
    nonunits = [a for a in range(1, q) if gcd(a, q) > 1]
    base = set(range(1, 14))
    found = []
    # dup+drop family: drop v, add lifted s+14j (residue s duplicated), j=1..3
    for v in range(1, 14):
        for s in range(1, 14):
            if s == v: continue
            for j in (1, 2, 3):
                lift = s + 14 * j
                S = sorted((base - {v}) | {lift})
                if len(S) != 13: continue
                M, wit = M_exact(S)
                if M == F(1, 14):
                    c1, c2, last = critical_slope(S, 14)
                    assert c1 == c2
                    found.append((float(c1), c1, v, s, j, S))
    cap14 = F(2, 14) * sum(F(1, a) for a in units)
    print(f"C_AP(14) = {cap14} = {float(cap14):.6f}")
    print(f"tight dup+drop sets found: {len(found)}")
    found.sort()
    for fc, c, v, s, j, S in found:
        vu = "UNIT" if gcd(v, 14) == 1 else "nonunit"
        su = "UNIT" if gcd(s, 14) == 1 else "nonunit"
        mark = "  <-- BEATS AP SLOPE" if c < cap14 else ("  (=AP)" if c == cap14 else "")
        print(f"  drop {v}({vu}), dup {s}({su}) via {s+14*j}: slope={c}={fc:.6f}{mark}")
        print(f"      S={S}")
    if found:
        cmin = min(c for _, c, *_ in [(f[0], f[1]) + f[2:] for f in found]) if False else min(f[1] for f in found)
        overall = min(cmin, cap14)
        print(f"\nTRUE inf slope over (AP + dup+drop tight locus) at N=14: {overall} = {float(overall):.6f}")
        if overall < cap14:
            print(">>> S92 'AP is the minimizer / floor 0.26(1-14r)' is CORRECTED: "
                  f"the floor slope is {float(overall):.6f}(1-14r).")
        # unit-residue lemma check: no tight set drops a unit
        dropped_units = [f for f in found if gcd(f[2], 14) == 1]
        print(f"tight sets dropping a UNIT residue: {len(dropped_units)} (lemma predicts 0)")

    print()
    print("=" * 78)
    print("(bonus) tangent-ladder demo on the deep well: certify m>0 up to 1/14")
    print("=" * 78)
    S = DW
    for r0 in [F(1, 28), F(1, 20), F(1, 16), F(13, 183), F(1, 14) - F(1, 500)]:
        m0 = lonely_measure(S, r0); sl = slope_from_components(S, r0)
        reach = r0 + (m0 / (-sl)) if sl < 0 else None
        print(f"r0={r0}={float(r0):.5f}: m={float(m0):.5f}, slope={float(sl):.4f}, "
              f"tangent zero at r={float(reach):.5f}" + ("  >= 1/14 ✓" if reach and reach >= F(1,14) else ""))

if __name__ == "__main__":
    main()
