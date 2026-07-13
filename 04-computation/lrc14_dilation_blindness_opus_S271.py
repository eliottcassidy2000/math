#!/usr/bin/env python3
"""
lrc14_dilation_blindness_opus_S271.py
=====================================
opus-2026-07-13-S271.  THE DILATION-BLINDNESS THEOREM (exact verification).

THEOREM (elementary, proof in the docstring's sketch):
  Let R be a speed set, c >= 1, and consider the dilated rest c*R.  For any peel modulus v with
  g = gcd(c, v):
     (i)   |G'_{cR}| = |G'_R|                       (danger sets rescale measure-preservingly)
     (ii)  A_{cR}(tau) = A_R(c*tau mod 1)            (autocorrelation transport)
     (iii) disc_v(cR) = disc_{v/g}(R)                ({c j/v mod 1} = the (v/g)-grid, g-fold)
  Hence the THM-731 certificate at peel v against a PURE-DILATE rest cR equals the certificate
  at modulus v/g against R.  In particular (g=1): the peel-v perspective is DILATION-BLIND --
  on runner v's clock the pack cR is indistinguishable from R.

COROLLARY (AP-shadow).  For the compressed near-dilate 2{1..12} u {13} (the bounded-DC floor,
MISTAKE-141), the peel-13 certificate EQUALS the certificate of peel 13 against {1..12}, i.e. of
the AP {1..13} -- the tight extremal with L = 0.  Soundness (the certificate is a rigorous lower
bound on L) then FORCES disc_13({1..12}) >= 6|G'_{1..12}|^2: the peel-13 perspective cannot
certify ANY body c*{1..12} u {13} with gcd(c,13)=1, for EVERY c -- an infinite family of blind
perspectives, explained (not just observed).  The mid-band certificate must come from a
dilation-BREAKING peel (peel an even element: rest = 2({1..12} minus {k}) u {13} is not a dilate).

Proof sketch: (i) ||cw t|| < 1/14 iff ||w s|| < 1/14 with s = ct mod 1; t -> ct mod 1 is
measure-preserving on the circle.  (ii) same substitution in |G cap (G - tau)|.  (iii)
{c j / v mod 1 : j mod v} hits each point of the (v/g)-grid exactly g times.  QED.

This script verifies (i)+(iii) EXACTLY (Fractions) on:
  disc_13(2{1..12}) = disc_13({1..12})          [g=1: blind peel]
  disc_16(2{1..12}) = disc_8({1..12})           [g=2: descent]
  disc_24(2{1..12}) = disc_12({1..12})          [g=2: descent]
  disc_13(3{1..12}) = disc_13({1..12})          [c=3, g=1: same blindness, new c]
and the FORCED AP-shadow inequality disc_13({1..12}) >= 6 |G'_{1..12}|^2 (exact).
"""
from fractions import Fraction as F

THR = F(1, 14)

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
    return [((k + THR) / w, (k + 1 - THR) / w) for k in range(w)]

def good_intervals(speeds):
    cur = [(F(0), F(1))]
    for w in sorted(speeds): cur = intersect(cur, safe_set(w))
    return cur

def disc_exact(ivs, v):
    """disc at grid modulus v via the Bernoulli jump-pair form (exact)."""
    jumps = []
    for a, b in ivs:
        jumps.append((a, 1)); jumps.append((b, -1))
    B2 = lambda x: x * x - x + F(1, 6)
    tot = F(0)
    for xp, sp in jumps:
        for xq, sq in jumps:
            tot += sp * sq * B2((v * (xp - xq)) % 1)
    return tot / (2 * v * v)

R = list(range(1, 13))
R2 = [2 * w for w in R]
R3 = [3 * w for w in R]

ivR = good_intervals(R)
ivR2 = good_intervals(R2)
ivR3 = good_intervals(R3)

print("=" * 100)
print("DILATION-BLINDNESS THEOREM -- exact verification (all Fractions)")
print("=" * 100)

gR, gR2, gR3 = measure(ivR), measure(ivR2), measure(ivR3)
print("\n(i) measure invariance: |G'_{1..12}| = %s ; |G'_{2{1..12}}| = %s ; |G'_{3{1..12}}| = %s"
      % (gR, gR2, gR3))
print("    equal: %s / %s" % (gR == gR2, gR == gR3))
print("    interval counts: r(R)=%d, r(2R)=%d, r(3R)=%d  (r scales with c; measure does not)"
      % (len(ivR), len(ivR2), len(ivR3)))

print("\n(iii) disc transport  disc_v(cR) = disc_{v/gcd(c,v)}(R):")
tests = [
    ("disc_13(2R) = disc_13(R)   [c=2,v=13,g=1: BLIND]", disc_exact(ivR2, 13), disc_exact(ivR, 13)),
    ("disc_16(2R) = disc_8(R)    [c=2,v=16,g=2: descent]", disc_exact(ivR2, 16), disc_exact(ivR, 8)),
    ("disc_24(2R) = disc_12(R)   [c=2,v=24,g=2: descent]", disc_exact(ivR2, 24), disc_exact(ivR, 12)),
    ("disc_13(3R) = disc_13(R)   [c=3,v=13,g=1: BLIND]", disc_exact(ivR3, 13), disc_exact(ivR, 13)),
]
for name, lhs, rhs in tests:
    print("   %-52s  %s   (lhs=%s=%.6e)" % (name, "EQUAL" if lhs == rhs else "** DIFFER **", lhs, float(lhs)))

print("\nCOROLLARY (AP-shadow, forced by soundness since L({1..13}) = 0):")
d13 = disc_exact(ivR, 13)
rhs = 6 * gR * gR
print("   disc_13({1..12}) = %s = %.6e" % (d13, float(d13)))
print("   6 |G'_{1..12}|^2 = %s = %.6e" % (rhs, float(rhs)))
print("   forced inequality disc_13 >= 6|G'|^2 : %s  (ratio %.3f)"
      % ("HOLDS" if d13 >= rhs else "** VIOLATED (would contradict soundness!) **", float(d13 / rhs)))
print("\n   => the peel-13 certificate is <= 0 for EVERY body c*{1..12} u {13}, gcd(c,13)=1 --")
print("      an infinite family of provably blind perspectives; mid-band certificates must use")
print("      dilation-BREAKING peels.  (Matches S270 survey: peel 13 fails, peels 8..24 certify.)")
print("done.")
