#!/usr/bin/env python3
"""
lrc14_segment_bijection_lemma_opus_S279.py
==========================================
opus-2026-07-13-S279.  THE SEGMENT-BIJECTION LEMMA -- completing THM-745's pairing theorem.

THE PROOF (three lemmas; verified exactly below).

LEMMA B (mirror segment bijection).  The mirror mu(u) = 1-u maps the exposed set of the
(j, +1/14)-line bijectively onto the exposed set of the (j, -1/14)-line, sending each maximal
segment (u1, u2) to (1-u2, 1-u1).  Proof: (M1) G_B is mu-invariant (||b(1-u)|| = ||bu||);
(M2) sigma in A_{j'}(u) <=> -sigma in A_{j'}(-u) (compute both arcs); (M3) the image of the
right endpoint r_j(u) = -ju + 1/14 at u is  1 - r_j(u) = ju - 1/14 = ell_j(1-u) mod 1 (using
j in Z), the LEFT endpoint at 1-u; (M4) exposure transports: r_j(u) interior to A_{j'}(u)
<=> ell_j(1-u) interior to A_{j'}(1-u), by M2+M3.  QED.

LEMMA C (grid-count matching).  Paired segments have EQUAL full-inside crossing counts:
K+1 = floor(u2 W) - ceil(u1 W)  vs  floor((1-u1)W) - ceil((1-u2)W), equal by
floor(W - x) = W - ceil(x) and ceil(W - x) = W - floor(x) (all cases incl. integral x).  QED.

LEMMA A (no-wrap).  If 0 in J and W > 14 max(J): every exposed crossing height lies in
[1/14, 13/14] (heights in (0, 1/14) or (13/14, 1) are strictly interior to the static j=0 arc,
hence buried), so h_m >= 1/14 > j/W = alpha: NO wrap terms; hence rho_seg = -(K+1) alpha^2 / 2
EXACTLY (purely deterministic).  QED.

PAIRING THEOREM (proved for W > 14 max(J)):  rho(j,+) = Sum_segs -(K+1)alpha^2/2 =
(Lemmas B+C, termwise) = rho(j,-).  Hence the orient-weighted residual total vanishes and the
first-order wedge content equals the exact dF_ext telescoping sum.  Below the threshold the
defect is Sum_wraps (2h - alpha) per pair (empirically absent at all tested W).

VERIFICATION: (1) Lemma B computationally: the (j,-) segment list EQUALS the mirrored (j,+)
list as exact Fractions, every j, both shapes; (2) Lemma C: K-counts match per pair at several
W; (3) Lemma A: no wraps at any exposed crossing for W in {160, 250, 800} (> 154) AND the
tested sub-threshold W in {90, 97, 150}; (4) rho termwise = -(K+1)alpha^2/2 in all no-wrap
cases; (5) the pairing totals (re-check).
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

def segments_with_heights(J, GB, pts):
    pieces = []
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        for da, db in GB:
            lo, hi = max(a, da), min(b, db)
            if hi > lo: pieces.append((lo, hi))
    pieces.sort()
    segs = []
    for j in J:
        if j == 0: continue
        for sign in (1, -1):
            run = None; prev_hi = None
            for lo, hi in pieces:
                mid = (lo + hi) / 2
                x = (-j * mid + F(sign, 14)) % 1
                exposed = True
                for j2 in J:
                    if j2 == j: continue
                    d = (x - (-j2 * mid)) % 1
                    if (d < LAM or d > 1 - LAM) and d != LAM and d != 1 - LAM:
                        exposed = False; break
                if exposed:
                    if run is not None and prev_hi == lo: run = (run[0], hi)
                    else:
                        if run is not None: segs.append((j, sign, run[0], run[1]))
                        run = (lo, hi)
                    prev_hi = hi
                else:
                    if run is not None: segs.append((j, sign, run[0], run[1])); run = None
            if run is not None: segs.append((j, sign, run[0], run[1]))
    return segs

CFG = [("shape 1 {1}u{W..W+11}", [1], list(range(12))),
       ("shape 2 {1,2,3}u{W..W+9}", [1, 2, 3], list(range(10)))]

print("=" * 100)
print("THE SEGMENT-BIJECTION LEMMA -- exact verification (completes THM-745's pairing theorem)")
print("=" * 100)

for name, Bset, J in CFG:
    GB = good_intervals(Bset)
    pts = breakpoints(J, [e for iv in GB for e in iv])
    segs = segments_with_heights(J, GB, pts)
    print("\n%s: %d segments total" % (name, len(segs)))

    # (1) Lemma B: mirrored (j,+) list == (j,-) list, exactly
    okB = True
    for j in sorted(set(j for (j, s, a, b) in segs)):
        plus = sorted([(a, b) for (jj, s, a, b) in segs if jj == j and s == 1])
        minus = sorted([(a, b) for (jj, s, a, b) in segs if jj == j and s == -1])
        mirrored = sorted([(1 - b, 1 - a) for (a, b) in plus])
        if mirrored != minus:
            okB = False
            print("   LEMMA B FAILS at j=%d" % j)
    print("   Lemma B (segment bijection, exact Fractions): %s" % ("VERIFIED -- (j,-) = mirror of (j,+) for every j" if okB else "**FAILED**"))

    # (2)+(3)+(4): per W checks
    for W in [90, 97, 150, 160, 250, 800]:
        okC = True; wraps_found = 0; okrho = True; okpair = True
        from collections import defaultdict
        tot = defaultdict(lambda: F(0))
        for (j, sign, u1, u2) in segs:
            m1 = (u1 * W).__ceil__(); m2 = (u2 * W).__floor__() - 1
            K1 = m2 - m1 + 1
            # partner count
            p1 = ((1 - u2) * W).__ceil__(); p2 = ((1 - u1) * W).__floor__() - 1
            if (p2 - p1 + 1) != K1: okC = False
            if K1 <= 0: continue
            alpha = F(j, W)
            h = (F(sign, 14) - F(j * m1, W)) % 1
            rho = F(0)
            for m in range(m1, m2 + 1):
                rho -= alpha * alpha / 2
                if h < alpha:
                    wraps_found += 1
                    rho += alpha - h
                    h = h - alpha + 1
                else:
                    h = h - alpha
            if rho != -K1 * alpha * alpha / 2 and wraps_found == 0: okrho = False
            tot[(j, sign)] += rho
        for j in sorted(set(j for (j, s) in tot)):
            if tot.get((j, 1), F(0)) != tot.get((j, -1), F(0)): okpair = False
        thr = "  (W %s 14*maxJ = %d)" % (">" if W > 14 * max(J) else "<=", 14 * max(J))
        print("   W=%4d: LemC K-match %s ; wraps at exposed crossings = %d ; rho = -(K+1)a^2/2 %s ; pairing %s%s"
              % (W, "OK" if okC else "**FAIL**", wraps_found,
                 "OK" if okrho else "**FAIL**", "EXACT" if okpair else "**FAIL**", thr))

print("\ndone.")
