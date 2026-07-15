#!/usr/bin/env python3
"""
thm826_farey_profile_referee_kps_S128c10.py
===========================================
kind-pasteur-2026-07-15-S128 (cont.10). REFEREE for THM-826 (the Farey profile theorem)
+ the cont.10 FORMULA HARVEST (small exact statements, each checked in Q).
"""
import sys
from fractions import Fraction as F
from math import gcd, comb
sys.stdout.reconfigure(line_buffering=True)

def good_measure(k, lam):
    pieces = []
    for w in range(1, k + 1):
        r = lam / w
        for a in range(w):
            c = F(a, w); lo, hi = c - r, c + r
            if lo < 0: pieces.append((F(0), hi)); pieces.append((lo + 1, F(1)))
            elif hi > 1: pieces.append((lo, F(1))); pieces.append((F(0), hi - 1))
            else: pieces.append((lo, hi))
    pieces.sort(); tot = F(0); cur = F(0)
    for lo, hi in pieces:
        if lo > cur: tot += lo - cur
        cur = max(cur, hi)
    if cur < 1: tot += 1 - cur
    return tot

def farey_pairs(k):
    """consecutive pairs of F_k on the circle, returned as denominator pairs (i, j)."""
    fr = sorted(set(F(a, i) for i in range(1, k + 1) for a in range(i)))
    out = []
    for x, y in zip(fr, fr[1:] + [F(1)]):
        out.append((x.denominator, y.denominator if y != 1 else 1))
    return out

def profile(k, lam):
    return sum(max(F(0), (1 - lam * (i + j))) / (i * j) for i, j in farey_pairs(k))

print("=" * 90)
print("THM-826 REFEREE: profile formula vs direct measure, k=2..12")
print("=" * 90)
allok = True
for k in range(2, 13):
    lams = set(F(c, 720) for c in range(0, 720 // (k + 1) + 1)) | \
           set(F(1, s) for s in range(k + 1, 2 * k)) | {F(1, k + 2), F(0)}
    bad = 0
    for lam in lams:
        if lam > F(1, k + 1):
            continue
        if good_measure(k, lam) != profile(k, lam):
            bad += 1
            if bad < 3:
                print("  k=%d MISMATCH at lam=%s: direct=%s formula=%s"
                      % (k, lam, good_measure(k, lam), profile(k, lam)))
    pairs = farey_pairs(k)
    A = sum(F(1, i * j) for i, j in pairs)
    B = sum(F(i + j, i * j) for i, j in pairs)
    Btot = 2 * sum(F(sum(1 for a in range(1, l + 1) if gcd(a, l) == 1), l) for l in range(1, k + 1))
    # phi(1)=1: count a in 1..l coprime -> for l=1 that's a=1: phi(1)=1 ok
    ss = sorted(set(i + j for i, j in pairs))
    ok = (bad == 0) and (A == 1) and (B == Btot) and ss[0] == k + 1 and ss[-1] <= 2 * k - 1
    # THM-817 segment check
    q = k + 1
    Hprim = sum(F(1, u) for u in range(1, k + 1) if gcd(u, q) == 1)
    t817 = profile(k, F(1, k + 2)) == 2 * Hprim / (q * (k + 2))
    ok &= t817
    allok &= ok
    print("k=%2d gaps=%3d  formula==direct on %2d lams: %s  A=1:%s  B=2*sum(phi/l):%s  i+j range %d..%d  THM-817 seg:%s"
          % (k, len(pairs), len([l for l in lams if l <= F(1, k + 1)]), bad == 0,
             A == 1, B == Btot, ss[0], ss[-1], t817))
print("PROFILE THEOREM: %s" % ("ALL EXACT" if allok else "FAILURES"))

print()
print("=" * 90)
print("FORMULA HARVEST (cont.10): small exact statements, each verified here")
print("=" * 90)
# (H1) area under profile = sum 1/(2ij(i+j))  [Franel-type]
for k in [3, 5, 8]:
    pairs = farey_pairs(k)
    area_pred = sum(F(1, 2 * i * j * (i + j)) for i, j in pairs)
    # integrate the piecewise-linear profile exactly: breakpoints
    bps = sorted(set([F(0)] + [F(1, i + j) for i, j in pairs] + [F(1, k + 1)]))
    bps = [b for b in bps if b <= F(1, k + 1)]
    area = F(0)
    for x, y in zip(bps, bps[1:]):
        area += (profile(k, x) + profile(k, y)) / 2 * (y - x)
    print("(H1) k=%d: area under profile = %s == Franel sum %s : %s" % (k, area, area_pred, area == area_pred))
# (H2) the number of Farey gaps = |F_k| = 1 + sum phi(l) ; gaps with i+j = k+1 <-> THM-817 witnesses phi(k+1)
for k in [4, 6, 10]:
    pairs = farey_pairs(k)
    n_gaps = len(pairs)
    tot = 1 + sum(sum(1 for a in range(1, l + 1) if gcd(a, l) == 1) for l in range(1, k + 1)) - 1
    first_seg = sum(1 for i, j in pairs if i + j == k + 1)
    phi_k1 = sum(1 for a in range(1, k + 1) if gcd(a, k + 1) == 1)
    print("(H2) k=%d: #gaps=%d==|F_k|=%d:%s ; #(i+j=k+1 gaps)=%d==phi(k+1)=%d:%s"
          % (k, n_gaps, tot, n_gaps == tot, first_seg, phi_k1, first_seg == phi_k1))
# (H3) profile at lambda = 1/(2k): m = 1 - (1/2k)*B_partial... spot exact values as dictionary rows
for k in [5, 12]:
    lam = F(1, 2 * k)
    print("(H3) k=%d: m({1..k}; 1/(2k)) = %s (exact rational, from the formula: %s)"
          % (k, good_measure(k, lam), profile(k, lam)))
# (H4) tournament side: kappa(N_n) = (n-2)! = # spanning trees = # one-tile-per-row transversals;
#      the profile's first-segment witness count phi(k+1) at k=n-2 = phi(n-1) = # primitive witnesses
#      => at n with n-1 prime: witnesses n-2 = rows of the staircase (one witness per row width!)
for n in [8, 14]:
    k = n - 2
    phi = sum(1 for a in range(1, k + 1) if gcd(a, k + 1) == 1)
    print("(H4) n=%d: staircase rows = %d ; profile witnesses phi(n-1) = %d ; equal iff n-1 prime: %s"
          % (n, k, phi, phi == k))
print("\nDONE")
