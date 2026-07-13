#!/usr/bin/env python3
"""
lrc14_peel_survey_perspectives_opus_S270.py
===========================================
opus-2026-07-13-S270.  THE PEEL FAMILY AS A PERSPECTIVE FAMILY (owner's frame).

A 14-runner system = origin + 13 movers.  Each mover i is a PERSPECTIVE: its n-1 = 13 arcs are
its relative speeds {v_j - v_i} u {-v_i}.  Among the (n-1)^2 = 169 perspective-arc slots,
T(n-2) = C(13,2) = 78 pair-arcs are double-counted (each |v_i - v_j| seen from both ends) and
13 spoke-arcs are single-counted -- the cut (+) cycle split of K_14 (13 = dim cut, 78 = dim cycle
= the n=14 staircase tile count m = C(n-1,2)).

The peeling identity  L = (6/7)|G'_{~v}| - eps_v  holds for EVERY peel v: 13 identities = the 13
perspectives.  klein-S287/288 (THM-731/732) uses the FAR peel (fine v-grid); the runner-1 lemma
(opus-S265, kps cont.70) is the CORE peel.  This script maps the WHOLE family:

For each covering body and EVERY peel v:
  r      = #intervals of the leave-one-out good set G'_{~v}         (exact, Fractions)
  |G'|   = its measure                                              (exact)
  M732   = 6|G'| - sqrt(2)*r/v     crude certificate margin (THM-732: L >= M732/7; certifies iff >0)
  disc_v = v-grid discrepancy of the autocorrelation of G'_{~v}     (FFT float, klein's method)
  M731   = 6|G'| - sqrt(6*disc_v)  true-disc certificate margin (THM-731; certifies iff >0)
plus the body's true L (exact) and the PERSPECTIVE CENSUS:
  - pair-arc collapse: #distinct values among the 78 pairwise |v_i-v_j| and 13 spokes (91 arcs);
  - frame collapse: per frame i, #distinct relative speeds (<=12 distinct => LRC(<=13) makes
    runner i lonely at 1/13 -- the near-AP frames degenerate).

STRESS TEST: the compressed near-dilate 2*{1..12} u {13} (MISTAKE-141's true bounded-DC floor,
M = 1/13) has NO far element -- does ANY peel certify?  This probes the mid-band residual that
klein-S288 flags (small-far-element families need the shared density-Q_s cancellation).

ALSO: independent verification of the exact Bernoulli form (kps-S128's THM-732 reservation):
  disc_v = (1/(2 v^2)) * Sum_{p,q} sigma_p sigma_q B2({v(x_p - x_q)}),   B2(x) = x^2 - x + 1/6,
derived here via  chat_l = [Sum_p sigma_p e(-l x_p)]/(2 pi i l)  and  Sum_{l!=0} e(l theta)/l^2
= 2 pi^2 B2({theta}) -- checked in EXACT rational arithmetic against the FFT on showcase peels.
"""
from fractions import Fraction as F
import math
import numpy as np

THR = F(1, 14)
NG = 1 << 21

# ---------------- exact interval engine ----------------

def normalize(ivs):
    ivs = sorted((a, b) for a, b in ivs if b > a)
    out = []
    for a, b in ivs:
        if out and a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
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

def measure(A):
    return sum(b - a for a, b in A)

def safe_set(w):
    return [((k + THR) / w, (k + 1 - THR) / w) for k in range(w)]

def good_intervals(speeds):
    cur = [(F(0), F(1))]
    for w in sorted(speeds):
        cur = intersect(cur, safe_set(w))
    return cur

# ---------------- float FFT disc (klein's method) ----------------

def good_indicator_float(S):
    t = np.arange(NG, dtype=np.float64) / NG
    g = np.ones(NG, dtype=np.float64)
    for w in S:
        frac = (w * t) % 1.0
        g *= (np.minimum(frac, 1.0 - frac) >= 1.0 / 14.0)
    return g

def disc_fft(S, v):
    g = good_indicator_float(S)
    G = np.fft.rfft(g)
    A = np.fft.irfft(G * np.conj(G), n=NG) / NG
    idx = np.round((np.arange(v) / v) * NG).astype(np.int64) % NG
    return float(A[idx].mean() - A.mean())

# ---------------- exact Bernoulli disc (independent derivation) ----------------

def disc_bernoulli_exact(ivs, v):
    """disc_v = (1/(2 v^2)) Sum_{p,q} s_p s_q B2({v (x_p - x_q)}), exact Fractions."""
    jumps = []
    for a, b in ivs:
        jumps.append((a, 1)); jumps.append((b, -1))
    tot = F(0)
    B2 = lambda x: x * x - x + F(1, 6)
    for xp, sp in jumps:
        for xq, sq in jumps:
            d = (v * (xp - xq)) % 1
            tot += sp * sq * B2(d)
    return tot / (2 * v * v)

# ---------------- perspective census ----------------

def census(body):
    pts = [0] + sorted(body)
    diffs = set()
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            diffs.add(abs(pts[j] - pts[i]))
    frames = []
    for i, vi in enumerate(pts):
        rel = {abs(vj - vi) for j, vj in enumerate(pts) if j != i}
        frames.append((vi, len(rel)))
    return len(diffs), frames

# ---------------- survey ----------------

BODIES = [
    ("deep well {1..12,182}",            [1,2,3,4,5,6,7,8,9,10,11,12,182]),
    ("residue {1..11,13,84}",            [1,2,3,4,5,6,7,8,9,10,11,13,84]),
    ("consecutive {2..14}",              [2,3,4,5,6,7,8,9,10,11,12,13,14]),
    ("compressed 2{1..12}u{13} MID-BAND",[2,4,6,8,10,12,14,16,18,20,22,24,13]),
    ("klein limitation {1,3..14}",       [1,3,4,5,6,7,8,9,10,11,12,13,14]),
]

print("=" * 100)
print("PEEL SURVEY -- all 13 perspectives per body.  M732 = 6|G'|-sqrt2*r/v (crude, PROVED);")
print("M731 = 6|G'|-sqrt(6*disc_v) (true-disc).  Certifies L>0 iff margin>0 (L >= margin/7).")
print("=" * 100)

for name, body in BODIES:
    Lexact = measure(good_intervals(body))
    ndiff, frames = census(body)
    print("\n%s   true L = %.6f (exact %s)" % (name, float(Lexact), measure(good_intervals(body))))
    print("   perspective census: distinct arcs among 91 = %d (collapse %.0f%%);  frames with <=12"
          % (ndiff, 100.0 * (1 - ndiff / 91.0)))
    degen = [(vi, c) for vi, c in frames if c <= 12]
    print("   degenerate frames (LRC(<=13) applies -> that runner IS 1/13-lonely): %s"
          % (", ".join("v=%d(%d)" % (vi, c) for vi, c in degen) if degen else "none"))
    print("   %5s %4s %10s %12s %12s %12s %12s  %s" %
          ("peel", "r", "|G'_~v|", "M732", "disc_v", "M731", "L_cert731", "verdict"))
    best = None
    for v in sorted(body):
        rest = [w for w in body if w != v]
        ivs = good_intervals(rest)
        r = len(ivs)
        Gp = measure(ivs)
        M732 = 6 * float(Gp) - math.sqrt(2) * r / v
        dv = disc_fft(rest, v)
        M731 = 6 * float(Gp) - math.sqrt(max(6 * dv, 0.0))
        verdict = ("CERT-732" if M732 > 0 else ("cert-731" if M731 > 0 else "--"))
        if best is None or M731 > best[1]: best = (v, M731, M732)
        print("   %5d %4d %10.6f %12.6f %12.3e %12.6f %12.6f  %s" %
              (v, r, float(Gp), M732, dv, M731, M731 / 7, verdict))
    print("   -> best peel v=%d: M731=%.6f (L>=%.6f); crude-732 margin there %.6f" %
          (best[0], best[1], best[1] / 7, best[2]))

# ---------------- exact Bernoulli-vs-FFT verification ----------------

print("\n" + "=" * 100)
print("EXACT BERNOULLI FORM CHECK (independent support for kps-S128's THM-732 reservation)")
print("disc_v = (1/2v^2) Sum_{p,q} s_p s_q B2({v(x_p-x_q)})  -- exact Fractions vs FFT float:")
for name, body, v in [("deep well", [1,2,3,4,5,6,7,8,9,10,11,12,182], 182),
                      ("residue",   [1,2,3,4,5,6,7,8,9,10,11,13,84],  84),
                      ("compressed",[2,4,6,8,10,12,14,16,18,20,22,24,13], 13)]:
    rest = [w for w in body if w != v]
    ivs = good_intervals(rest)
    de = disc_bernoulli_exact(ivs, v)
    df = disc_fft(rest, v)
    print("   %-11s peel %3d:  exact = %-24s = %.6e   FFT = %.6e   ratio %.4f"
          % (name, v, str(de), float(de), df, float(de) / df if df else float('nan')))

print("\ndone.")
