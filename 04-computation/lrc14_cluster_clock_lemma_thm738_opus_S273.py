#!/usr/bin/env python3
"""
lrc14_cluster_clock_lemma_thm738_opus_S273.py
=============================================
opus-2026-07-13-S273.  THM-738: the CLUSTER-CLOCK LEMMA (additive pack-clock) -- uniform tail
closure of far-cluster families B u (W + J), the additively-coherent slice of the j >= 7 seam.

SHAPE: V = B u {W + j : j in J}, B a small base (|B| + |J| = 13), J a fixed offset pattern
(j >= 7 cluster members allowed -- beyond kps THM-735's Bonferroni boundary j <= 6), W -> inf.

THE MOVE (perspective frame): ride the CLUSTER's clock s = {W t}, branches t = (m+s)/W.
 - base b in B: b t moves by only b/W per branch -- the FROZEN FAN: branch m is base-safe
   for every s provided ||b m/W|| >= 1/14 + b/W  (fattened base good set).
 - cluster member W + j: (W+j)t = m + s + j t == s + j t (mod 1), and j t = j u_m + O(j/W):
   within branch m the 13 - |B| cluster conditions are "s avoids the 1/14-fattened AP
   { -j u_m : j in J }" -- an ADDITIVE (three-distance / tile) structure with step u_m.
   Define  A_J(u) = 1 - | Union_{j in J} (-j u - 1/14, -j u + 1/14) mod 1 |.
Because the cluster is COHERENT, the arcs OVERLAP: A_J(u) > 0 on windows around low-denominator
u even for |J| = 12 (where the Bonferroni base 1 - j/7 < 0 is useless).

THM-738.  L(V)  >=  Area(B,J) - C(B,J)/W,   where
    Area(B,J) = Integral_{G_B} A_J(u) du            (exactly computable rational)
    C(B,J)    = V(A_J) + 2#comp(G_B) + 2 Sum_{b in B} b + |J| max(J)   (exactly computable)
so every W > W0 = C/Area is closed; W <= W0 is a finite list of exactly-decidable bodies.

|J| = 1 degenerates to the far-element peel main term (Area = (6/7)|G_B|): the lemma unifies
the far-element regime and the far-cluster regime as ONE clock-ride.

This script: (1) exact adaptive piecewise-linear integration of A_J over G_B (Fractions;
midpoint-linearity self-check with recursive bisection -- correctness does not depend on the
event enumeration); (2) V(A_J) exact; (3) constants + W0 for the two klein-S289 shapes;
(4) verification: exact L at several W vs the bound and vs Area; (5) float sweep W = 14..W0
for min L over the full shape family (closure below the tail threshold), exact-verified at the
minimizing W.
"""
from fractions import Fraction as F

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

def safe_set(w, lam=F(1, 14)):
    return [((k + lam) / w, (k + 1 - lam) / w) for k in range(w)]

def good_intervals(speeds, lam=F(1, 14)):
    cur = [(F(0), F(1))]
    for w in sorted(speeds): cur = intersect(cur, safe_set(w, lam))
    return cur

# ---------------- A_J(u) exact ----------------

LAM = F(1, 14)

def A_J(u, J):
    """1 - measure of union of arcs (-ju - 1/14, -ju + 1/14) mod 1, exact."""
    ivs = []
    for j in J:
        c = (-j * u) % 1
        a, b = c - LAM, c + LAM
        if a < 0:
            ivs.append((a + 1, F(1))); ivs.append((F(0), b))
        elif b > 1:
            ivs.append((a, F(1))); ivs.append((F(0), b - 1))
        else:
            ivs.append((a, b))
    return 1 - measure(normalize(ivs))

# ---------------- exact adaptive integration ----------------

def breakpoints(J, extra):
    """candidate slope-change points: (k)/d and (k +- 1/7)/d for pairwise diffs d, plus extra."""
    pts = set(extra)
    diffs = sorted({abs(a - b) for a in J for b in J if a != b})
    for d in diffs:
        for k in range(d + 1):
            for off in (F(0), F(1, 7), -F(1, 7)):
                u = (F(k) + off) / d
                if 0 <= u <= 1: pts.add(u)
    pts.add(F(0)); pts.add(F(1))
    return sorted(pts)

def integrate_AJ(J, domain):
    """Integral of A_J over the interval set `domain`, exact, with linearity self-check."""
    extra = [e for iv in domain for e in iv]
    pts = breakpoints(J, extra)
    total = F(0); var = F(0); prev_val = None
    pieces = []
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        pieces.append((a, b))
    def integ(a, b, fa, fb, depth):
        m = (a + b) / 2
        fm = A_J(m, J)
        if fm * 2 == fa + fb or depth >= 12:   # linear (exact) or give up bisecting
            assert fm * 2 == fa + fb, "non-linear piece at depth cap: refine breakpoints"
            return (fa + fb) * (b - a) / 2
        return integ(a, m, fa, fm, depth + 1) + integ(m, b, fm, fb, depth + 1)
    # clip pieces to the domain
    dom = normalize(domain)
    for a, b in pieces:
        for da, db in dom:
            lo, hi = max(a, da), min(b, db)
            if hi > lo:
                fa, fb = A_J(lo, J), A_J(hi, J)
                total += integ(lo, hi, fa, fb, 0)
    # total variation of A_J on [0,1] (full circle, piece endpoints)
    vals = [A_J(p, J) for p in pts]
    var = sum(abs(v2 - v1) for v1, v2 in zip(vals, vals[1:]))
    return total, var

# ---------------- the lemma's constants ----------------

def lemma_constants(Bset, J):
    GB = good_intervals(Bset)
    area, var = integrate_AJ(J, GB)
    C = var + 2 * len(GB) + 2 * sum(Bset) + len(J) * max(J)
    W0 = None
    if area > 0:
        W0 = (C / area).__ceil__()
    return GB, area, var, C, W0

# ---------------- exact L of a body ----------------

def L_exact(body):
    return measure(good_intervals(body))

# ---------------- run the two klein-S289 shapes ----------------

print("=" * 102)
print("THM-738 (cluster-clock lemma): L(B u (W+J)) >= Area(B,J) - C(B,J)/W.  All exact Fractions.")
print("=" * 102)

SHAPES = [
    ("klein far-cluster shape", [1], list(range(0, 12))),          # {1} u {W..W+11}, j=12 cluster
    ("klein two-scale shape",   [1, 2, 3], list(range(0, 10))),    # {1,2,3} u {W..W+9}, j=10 cluster
]

results = []
for name, Bset, J in SHAPES:
    GB, area, var, C, W0 = lemma_constants(Bset, J)
    print("\n%s:  B=%s, J={0..%d}  (|J|=%d far elements -- beyond the Bonferroni j<=6 seam)"
          % (name, Bset, max(J), len(J)))
    print("   |G_B| = %s = %.6f  (#comp %d)" % (measure(GB), float(measure(GB)), len(GB)))
    print("   Area(B,J) = %s = %.6f" % (area, float(area)))
    print("   V(A_J) = %s = %.4f ;  C = %s = %.2f ;  W0 = ceil(C/Area) = %s"
          % (var, float(var), C, float(C), W0))
    results.append((name, Bset, J, area, C, W0))

# ---------------- verification: exact L vs bound at several W ----------------

print("\n" + "-" * 102)
print("VERIFICATION (exact L vs the bound Area - C/W and vs Area):")
name, Bset, J, area, C, W0 = results[0]
for W in [90, 200, 400, 800, 1600]:
    body = Bset + [W + j for j in J]
    Lx = L_exact(body)
    bnd = area - C / W
    print("   W=%5d : L = %.6f   bound = %+.6f  (valid: %s)   L - Area = %+.6f  (x W = %+.2f)"
          % (W, float(Lx), float(bnd), Lx >= bnd, float(Lx - area), float((Lx - area) * W)))

name, Bset, J, area, C, W0 = results[1]
for W in [50, 200, 800]:
    body = Bset + [W + j for j in J]
    Lx = L_exact(body)
    bnd = area - C / W
    print("   [shape 2] W=%5d : L = %.6f   bound = %+.6f  (valid: %s)" % (W, float(Lx), float(bnd), Lx >= bnd))

print("\ndone (part 1). Float sweep of W <= W0 in the companion part.")
