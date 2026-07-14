#!/usr/bin/env python3
"""
lrc14_two_cluster_hierarchical_thm740_opus_S274.py
==================================================
opus-2026-07-13-S274.  THM-740: TWO-CLUSTER HIERARCHICAL CLOCKS -- proof verification, exact.

V = B u (W1 + J1) u (W2 + J2),  |B| + |J1| + |J2| = 13,  W2 >> W1 >> 1.

OUTER RIDE (= THM-739 instantiated, base B' = B u (W1+J1)):
    (*)  L(V) >= Int_{G1} A_{J2}(u) du - C_out(W1)/W2,
         G1 = good set of B',  C_out = V(A_{J2}) + 2#comp(G1) + 2 Sum_{B'} b + |J2| max(J2).
INNER RIDE (NEW -- the W1-clock on the weighted area):
    (**) Int_{G1} A_{J2}(u) du >= Area2 - C_in/W1,
         Area2 = Int_{G_B} A_{J1}(u) A_{J2}(u) du    (the PRODUCT-AREA),
         C_in  = V(A_{J1} A_{J2}) + V(A_{J2}) + 2#comp(G_B) + 2 Sum_B b + |J1| max(J1).
    Proof of (**): branch u = (n+sigma)/W1; base freezes (cost 2 Sum_B b / W1 measure +
    Riemann (V(A1A2)+2#comp)/W1); cluster-1 conditions become the 1/14-fattened AP in sigma
    (measure >= A_{J1}(x_n) - |J1|max(J1)/W1); the weight A_{J2}((n+sigma)/W1) >= A_{J2}(x_n)
    - osc_n, and Sum_n osc_n <= V(A_{J2}) since the branches partition [0,1).  Using
    (a-x)(b-y) >= ab - x - y for a,b <= 1, x,y >= 0 (checked also when b-y < 0), sum and
    Riemann-sum the product A_{J1} A_{J2} over G_B^fat.  QED.

COMBINED:  L(V) >= Area2 - C_in/W1 - C_out(W1)/W2,  with C_out(W1) <= alpha + 4|J1| (W1+max J1)
(the linear-in-W1 form via #comp(G1) <= 1 + Sum_{B'} w).  Corollary: Area2 > 0 closes the
separated cone W1 >= 3 C_in/Area2, W2 >= 3 C_out(W1)/Area2 (each error < Area2/3); every fixed
W1 reduces to a THM-739 instance (base B') for the rest.

PERSPECTIVE READING: within-cluster pairs = coherent (A_J); cross-cluster pairs = DECOUPLED
(the product structure -- separated clocks factorize); base-cluster pairs = frozen fan.  Every
pair-sector role of the (n-1)^2 decomposition now has a theorem.

EXACT VERIFICATION BELOW (all Fractions):
  A. Area2 by piecewise-QUADRATIC Simpson (exact for quadratics; two-half self-check per piece)
     for shape A: B={1}, J1=J2={0..5}, and shape B: B={1,2}, J1={0..4}, J2={0..5}.
     Compare with the S273 taste test (L ~ 0.097-0.101 at W1=30..60, W2=400..1600).
  B. V(A_{J2}), V(A1 A2) exact (interior-vertex handling for the quadratic pieces).
  C. Inner bound (**): exact Int_{G1} A_{J2} at W1 in {30, 60, 120} vs Area2 - C_in/W1.
  D. Combined bound: exact L(V) on a (W1, W2) grid vs Area2 - C_in/W1 - C_out(W1)/W2.
"""
from fractions import Fraction as F

LAM = F(1, 14)

# ---------------- interval engine ----------------

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

# ---------------- exact integration ----------------

def simpson(a, b, fa, fm, fb):
    return (b - a) * (fa + 4 * fm + fb) / 6

def integ_quad(f, a, b, depth=0):
    """exact integral of piecewise-quadratic f on [a,b] via Simpson + two-half self-check"""
    m = (a + b) / 2
    fa, fm, fb = f(a), f(m), f(b)
    I = simpson(a, b, fa, fm, fb)
    m1, m2 = (a + m) / 2, (m + b) / 2
    I2 = simpson(a, m, fa, f(m1), fm) + simpson(m, b, fm, f(m2), fb)
    if I == I2: return I
    assert depth < 14, "non-quadratic piece; missing breakpoint"
    return integ_quad(f, a, m, depth + 1) + integ_quad(f, m, b, depth + 1)

def var_quad(f, pts):
    """total variation of piecewise-quadratic f with pieces between consecutive pts"""
    V = F(0)
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        m = (a + b) / 2
        fa, fm, fb = f(a), f(m), f(b)
        # Lagrange coefficients: f = c2 u^2 + c1 u + c0
        c2 = (fa * (m - b) + fm * (b - a) + fb * (a - m)) / ((a - m) * (a - b) * (m - b))
        if c2 != 0:
            c1 = (fm - fa) / (m - a) - c2 * (m + a)
            us = -c1 / (2 * c2)
            if a < us < b:
                fus = fa + (us - a) * ((fm - fa) / (m - a) + c2 * (us - m))  # eval via form
                # safer: direct quadratic eval
                c0 = fa - c2 * a * a - c1 * a
                fus = c2 * us * us + c1 * us + c0
                V += abs(fus - fa) + abs(fb - fus)
                continue
        V += abs(fb - fa)
    return V

def area_and_consts(Bset, J1, J2):
    GB = good_intervals(Bset)
    extra = [e for iv in GB for e in iv]
    pts = breakpoints(sorted(set(J1) | set(J2)), extra)
    fprod = lambda u: A_J(u, J1) * A_J(u, J2)
    f2 = lambda u: A_J(u, J2)
    area2 = F(0)
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        for da, db in GB:
            lo, hi = max(a, da), min(b, db)
            if hi > lo: area2 += integ_quad(fprod, lo, hi)
    V2 = var_quad(f2, pts)      # A_{J2} is piecewise linear (quadratic machinery still exact)
    V12 = var_quad(fprod, pts)
    C_in = V12 + V2 + 2 * len(GB) + 2 * sum(Bset) + len(J1) * max(J1)
    return GB, area2, V2, V12, C_in

def C_out_exact(Bset, J1, J2, W1, V2):
    Bp = list(Bset) + [W1 + j for j in J1]
    G1 = good_intervals(Bp)
    return G1, V2 + 2 * len(G1) + 2 * sum(Bp) + len(J2) * max(J2)

def int_G1_A2(G1, J2):
    """exact integral of piecewise-LINEAR A_{J2} over interval set G1 (trapezoid + check)"""
    extra = [e for iv in G1 for e in iv]
    pts = breakpoints(J2, extra)
    f2 = lambda u: A_J(u, J2)
    tot = F(0)
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        for da, db in G1:
            lo, hi = max(a, da), min(b, db)
            if hi > lo: tot += integ_quad(f2, lo, hi)
    return tot

# ---------------- run ----------------

print("=" * 104)
print("THM-740 (two-cluster hierarchical clocks) -- exact verification.")
print("=" * 104)

SHAPES = [
    ("shape A", [1], list(range(6)), list(range(6))),       # {1} u {W1..W1+5} u {W2..W2+5}
    ("shape B", [1, 2], list(range(5)), list(range(6))),    # {1,2} u {W1..W1+4} u {W2..W2+5}
]

for name, Bset, J1, J2 in SHAPES:
    GB, area2, V2, V12, C_in = area_and_consts(Bset, J1, J2)
    print("\n%s: B=%s, J1={0..%d}, J2={0..%d}   (%d + %d + %d = 13 runners)"
          % (name, Bset, max(J1), max(J2), len(Bset), len(J1), len(J2)))
    print("   |G_B| = %s ;  PRODUCT-AREA Area2 = %s = %.6f" % (measure(GB), area2, float(area2)))
    print("   V(A_J2) = %s = %.4f ;  V(A_J1*A_J2) = %s = %.4f ;  C_in = %s = %.2f"
          % (V2, float(V2), V12, float(V12), C_in, float(C_in)))

    print("   INNER BOUND (**) check:  Int_{G1} A_J2  >=  Area2 - C_in/W1 :")
    for W1 in [30, 60, 120]:
        G1, C_out = C_out_exact(Bset, J1, J2, W1, V2)
        I1 = int_G1_A2(G1, J2)
        rhs = area2 - C_in / W1
        print("     W1=%4d : Int = %.6f   Area2 - C_in/W1 = %+.6f   valid=%s   (C_out(W1) = %s)"
              % (W1, float(I1), float(rhs), I1 >= rhs, C_out))

    print("   COMBINED BOUND check:  L(V) >= Area2 - C_in/W1 - C_out(W1)/W2 :")
    for W1 in [30, 120]:
        G1, C_out = C_out_exact(Bset, J1, J2, W1, V2)
        for W2 in [3200, 12800]:
            body = list(Bset) + [W1 + j for j in J1] + [W2 + j for j in J2]
            Lx = measure(good_intervals(body))
            bnd = area2 - C_in / W1 - F(C_out) / W2
            print("     W1=%4d W2=%6d : L = %.6f   bound = %+.6f   valid=%s"
                  % (W1, W2, float(Lx), float(bnd), Lx >= bnd))

print("\ndone.")
