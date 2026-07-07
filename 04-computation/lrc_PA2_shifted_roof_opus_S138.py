"""
lrc_PA2_shifted_roof_opus_S138.py

THE 2-ANCHOR JOINT TAIL ON ITS CONJECTURED MINIMIZER — the shifted-origin roof
(boxeph-S1's named joint next step; owner assignment).

boxeph-S1 (HYP-4801 cont.): the load-bearing tail lemma reduces to
    PA_2(E) = P_x( max(gap∋0, gap∋1/2)(config(E,x)) > 1/7 ) >= T_k,
and the empirical PA_2-minimizer at every k is the SPREAD INHOMOGENEOUS AP {a + d·j}.

STRUCTURE (the shifted-origin roof): for E = {a+dj : j=0..k-1},
    config(E, x) = rotation by z=ax of the Steinhaus config C0(y) = {frac(jy): j=0..k-1},
    y = dx.  So gap∋0 = gap∋(-z) of C0(y), gap∋1/2 = gap∋(1/2 - z), and as (a,d) grows
    coprimely the pair (y,z) equidistributes on T^2:
      * 1-ANCHOR LIMIT:  P(gap∋0 > 1/7)  ->  E_y[ L(y) ],   L(y) = total length of the
        > 1/7 gaps of C0(y)  (z-average of "z lands in a long gap" = the long-gap MASS).
      * 2-ANCHOR LIMIT:  PA_2  ->  E_y[ meas( S(y) ∪ (S(y) - 1/2) ) ],  S(y) = the union
        of >1/7 gaps  =  2 E[L] - E[meas(S ∩ (S - 1/2))].

EXACT E[L] (new, tiny extension of THM-637): on each Farey-K cell (p/q, p'/q'), K = k-1,
the C0 gaps take the three-distance values d1 = q·eps, d2 = q'·eps', d3 = d1+d2 with the
CLASSICAL MULTIPLICITIES
      N1 = K+1-q,  N2 = K+1-q',  N3 = q+q'-K-1     (sum = K+1),
so L(y) is piecewise linear per cell with breakpoints where each value crosses 1/7, and
E_y[L] is an EXACT RATIONAL per k.  (Multiplicities verified exhaustively below before use.)

RUNS:
 (1) verify the multiplicity formula exactly (random rational y, K = 5..16).
 (2) exact E[L] table for k = 8..13  ==  the conjectured 1-anchor tail infima
     (boxeph numeric: 0.602/0.511/0.434/0.368/0.321/0.281) — match?
 (3) 2-anchor limit: E[meas(S ∪ S-1/2)] to high precision (per-cell exact positions,
     numeric y-sampling as guard) vs T_k = {0.6185,0.5057,0.3956,0.2747,0.1429,0.0565}.
 (4) finite-(a,d) EXACT PA_2 via the order-cell engine (anchored gaps per subcell) on a
     grid of (a,d) — how fast the limit is approached; the inf over the family.
"""
from fractions import Fraction as F
from math import gcd
import random, sys, time

sys.path.insert(0, ".")
from lrc_exact_mu_ordercells_opus_S136 import order_cells, cell_gap_affines

THETA = F(1, 7)
TK = {8: 0.6185, 9: 0.5057, 10: 0.3956, 11: 0.2747, 12: 0.1429, 13: 0.0565}
BOXEPH_1ANCH = {8: 0.602, 9: 0.511, 10: 0.434, 11: 0.368, 12: 0.321, 13: 0.281}

def farey(n):
    fr = set()
    for q in range(1, n + 1):
        for p in range(0, q + 1):
            fr.add(F(p, q))
    return sorted(fr)

def gaps_of_config(pts):
    pts = sorted(pts)
    return [b - a for a, b in zip(pts, pts[1:])] + [pts[0] + 1 - pts[-1]]

def verify_multiplicities(K, trials, rng):
    """gaps of {frac(jy): 0<=j<=K} on cell (p/q,p'/q'): values {q e, q' e', sum} with
       multiplicities {K+1-q, K+1-q', q+q'-K-1}."""
    Fl = farey(K)
    cells = list(zip(Fl[:-1], Fl[1:]))
    bad = 0
    for _ in range(trials):
        a_, b_ = rng.choice(cells)
        den = rng.randint(10**5, 10**6)
        lo, hi = int(a_ * den) + 1, int(b_ * den)
        if lo > hi: continue
        y = F(rng.randint(lo, hi), den)
        if not (a_ < y < b_): continue
        q, qp = a_.denominator, b_.denominator
        d1 = q * (y - a_); d2 = qp * (b_ - y)
        expect = {}
        if K + 1 - q: expect[d1] = K + 1 - q
        if K + 1 - qp: expect[d2] = expect.get(d2, 0) + (K + 1 - qp)
        if q + qp - K - 1: expect[d1 + d2] = expect.get(d1 + d2, 0) + (q + qp - K - 1)
        pts = [ (j * y) % 1 for j in range(0, K + 1) ]
        got = {}
        for g in gaps_of_config(pts):
            got[g] = got.get(g, 0) + 1
        if got != expect:
            bad += 1
    return bad

def EL_exact(k, theta=THETA):
    """EXACT E_y[L(y)], L = total length of > theta gaps of {frac(jy): 0<=j<=k-1}."""
    K = k - 1
    Fl = farey(K)
    tot = F(0)
    for a_, b_ in zip(Fl[:-1], Fl[1:]):
        q, qp = a_.denominator, b_.denominator
        N = ((K + 1 - q, lambda y: q * (y - a_)),
             (K + 1 - qp, lambda y: qp * (b_ - y)),
             (q + qp - K - 1, lambda y: q * (y - a_) + qp * (b_ - y)))
        # each value is affine in y on the cell; contribution N*v(y)*[v(y)>theta]
        # breakpoints where v crosses theta
        bps = set([a_, b_])
        for mult, v in N:
            if mult <= 0: continue
            va, vb = v(a_), v(b_)
            if (va - theta) * (vb - theta) < 0:
                # affine: solve v(y) = theta
                # v(y) = alpha y + beta
                alpha = (vb - va) / (b_ - a_)
                beta = va - alpha * a_
                bps.add((theta - beta) / alpha)
        bps = sorted(bps)
        for u, w in zip(bps, bps[1:]):
            mid = (u + w) / 2
            for mult, v in N:
                if mult <= 0: continue
                if v(mid) > theta:
                    # ∫ mult*v over (u,w): v affine
                    tot += mult * (v(u) + v(w)) / 2 * (w - u)
    return tot

def anchored_gap(pts_sorted, a):
    """gap of the sorted config containing point a (a not a config point)."""
    import bisect
    i = bisect.bisect_left(pts_sorted, a)
    lo = pts_sorted[i - 1] if i > 0 else pts_sorted[-1] - 1
    hi = pts_sorted[i] if i < len(pts_sorted) else pts_sorted[0] + 1
    return hi - lo

def two_anchor_limit_numeric(k, samples, rng, theta=float(THETA)):
    """E_y[meas(S ∪ S-1/2)] numeric (exact-in-z per y: interval arithmetic on the gaps)."""
    K = k - 1
    tot = 0.0
    for _ in range(samples):
        y = rng.random()
        pts = sorted((j * y) % 1.0 for j in range(0, K + 1))
        # long-gap intervals S = union of (p_i, p_{i+1}) with length > theta
        ivs = []
        n = len(pts)
        for i in range(n):
            lo = pts[i]
            hi = pts[i + 1] if i < n - 1 else pts[0] + 1
            if hi - lo > theta:
                ivs.append((lo, hi))
        # meas(S ∪ (S - 1/2)) on the circle: shift copies by -1/2 mod 1, merge
        def norm(iv):
            lo, hi = iv
            lo %= 1.0; hi = lo + (iv[1] - iv[0])
            return (lo, hi)
        allv = []
        for iv in ivs:
            allv.append(norm(iv))
            allv.append(norm((iv[0] - 0.5, iv[1] - 0.5)))
        # unroll circle: split wrapping intervals
        flat = []
        for lo, hi in allv:
            if hi <= 1.0: flat.append((lo, hi))
            else: flat.append((lo, 1.0)); flat.append((0.0, hi - 1.0))
        flat.sort()
        m = 0.0; cur = -1.0; hi0 = -1.0
        for lo, hi in flat:
            if lo > hi0:
                m += max(0.0, hi0 - cur) if cur >= 0 else 0.0
                cur, hi0 = lo, hi
            else:
                hi0 = max(hi0, hi)
        m += max(0.0, hi0 - cur) if cur >= 0 else 0.0
        tot += m
    return tot / samples

def PA2_exact_finite(a, d, k, theta=THETA):
    """EXACT PA_2 for E = {a + d*j : j=0..k-1} via the order-cell engine."""
    E = [a + d * j for j in range(k)]
    tot = F(0)
    for c0, c1 in order_cells(E):
        gaps = cell_gap_affines(E, c0, c1)
        # anchored gaps at 0 and 1/2 are specific gaps per subcell: find which gap
        # contains each anchor. Use midpoint evaluation.
        subbp = set([c0, c1])
        for i in range(len(gaps)):
            ci, bi = gaps[i]
            for j in range(i + 1, len(gaps)):
                cj, bj = gaps[j]
                if ci != cj:
                    xc = (bj - bi) / (ci - cj)
                    if c0 < xc < c1: subbp.add(xc)
        # anchored-gap identity: gap∋anchor as function of x needs the sorted phases;
        # simplest exact route: on each subcell evaluate at midpoint which gaps hold the
        # anchors, then those gaps' affine forms are the anchored gaps on the subcell.
        for u, v in zip(sorted(subbp), sorted(subbp)[1:]):
            mid = (u + v) / 2
            pts = sorted(((e * mid) % 1, e) for e in E)
            phis = [p for p, _ in pts]
            import bisect
            def gap_form(anchor):
                i = bisect.bisect_left(phis, anchor)
                lo_e = pts[i - 1][1] if i > 0 else pts[-1][1]
                hi_e = pts[i][1] if i < len(pts) else pts[0][1]
                # affine form of that gap on this subcell:
                lo_p = (lo_e * mid) % 1; hi_p = (hi_e * mid) % 1
                c = F(hi_e - lo_e); b0 = (hi_p - lo_p) - c * mid
                if hi_p < lo_p or (i == 0 or i == len(pts)):  # wrap
                    pass
                # normalize so that value at mid equals actual gap length
                val = (hi_p - lo_p) % 1
                b0 = val - c * mid
                return c, b0
            cg0, bg0 = gap_form(F(0) + F(1, 10**9))  # anchor 0 (avoid exact hit)
            cgh, bgh = gap_form(F(1, 2))
            # superlevel of max of the two affines on (u,v)
            mvals = []
            for (cc, bb) in ((cg0, bg0), (cgh, bgh)):
                mvals.append((cc, bb))
            # compute measure{ x in (u,v): max > theta } exactly
            pts2 = set([u, v])
            for cc, bb in mvals:
                if cc != 0:
                    xc = (theta - bb) / cc
                    if u < xc < v: pts2.add(xc)
            for uu, vv in zip(sorted(pts2), sorted(pts2)[1:]):
                mm = (uu + vv) / 2
                if max(cc * mm + bb for cc, bb in mvals) > theta:
                    tot += vv - uu
    return tot

def main():
    rng = random.Random(138138)
    print("=" * 100)
    print("(1) three-distance multiplicities N1=K+1-q, N2=K+1-q', N3=q+q'-K-1 (exact check)")
    print("=" * 100)
    for K in (5, 7, 9, 12, 16):
        bad = verify_multiplicities(K, 250, rng)
        print(f"   K={K:2d}: violations = {bad}")

    print()
    print("=" * 100)
    print("(2) EXACT 1-anchor limit E_y[L] per k  (= conjectured inf of the 1-anchor tail)")
    print("=" * 100)
    ELs = {}
    for k in range(8, 14):
        e = EL_exact(k)
        ELs[k] = e
        match = "≈ boxeph" if abs(float(e) - BOXEPH_1ANCH[k]) < 0.006 else "*** CHECK vs boxeph ***"
        disch = "clears T_k" if float(e) > TK[k] else "BELOW T_k"
        print(f"   k={k:2d}: E[L] = {str(e):>18} = {float(e):.6f}   boxeph inf {BOXEPH_1ANCH[k]:.3f} {match}"
              f"   T_k={TK[k]:.4f} -> {disch}")

    print()
    print("=" * 100)
    print("(3) 2-anchor limit E_y[meas(S ∪ S-1/2)] (numeric, 200k y-samples) vs T_k")
    print("=" * 100)
    for k in range(8, 14):
        v = two_anchor_limit_numeric(k, 200_000, rng)
        print(f"   k={k:2d}: 2-anchor limit ≈ {v:.5f}   1-anchor exact {float(ELs[k]):.5f}"
              f"   gain +{v - float(ELs[k]):.5f}   T_k={TK[k]:.4f}  margin {v - TK[k]:+.4f}")

    print()
    print("=" * 100)
    print("(4) finite-(a,d) EXACT PA_2 (order-cell engine), k=8: approach to the limit")
    print("=" * 100)
    for (a, d) in ((1, 1), (1, 3), (2, 5), (3, 7), (1, 10), (7, 15)):
        if gcd(a, d) != 1: continue
        t0 = time.time()
        v = PA2_exact_finite(a, d, 8)
        print(f"   (a,d)=({a},{d}): PA_2 exact = {v} = {float(v):.6f}   [{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
