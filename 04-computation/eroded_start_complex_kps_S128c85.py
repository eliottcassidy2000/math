#!/usr/bin/env python3
"""eroded_start_complex_kps_S128c85.py -- kind-pasteur-2026-07-19-S128c85

SUPPLIER (3) of THM-1162 s3 / THM-1203 s12 -- the last open piece of uniform r=5
that is not codex's finite-to-continuum error term:

  "an ERODED core-atlas bound ensuring that a safe phase outside the bad set
   contains a COMPLETE k1-gap, rather than merely a safe point."

FORMALISATION.  For a core P (8 speeds from {1,...,12}) the core-safe set is
    S(P) = { t : ||p t|| >= 1/14  for all p in P }
        = intersection over p in P of  union_j [ j/p + 1/(14p), (j+1)/p - 1/(14p) ],
a finite union of closed intervals (never wrapping, since t=0 is always killed).
A complete k1-gap is an interval of width
    w = 6/(7 k1)
(the run between consecutive teeth of the k1 comb, each tooth of width 1/(7k1)).
The set of START points of complete k1-gaps lying inside S(P) is the MORPHOLOGICAL
EROSION of S(P) by a window of width w:

    E_w(P) = { t : [t, t+w] subset S(P) },      |E_w(P)| = sum_j max(0, L_j - w)

over the component lengths L_j of S(P).  This is the "eroded start complex".

WHY MEASURE IS THE RIGHT CURRENCY.  The k1 gaps sit at fixed positions j/k1, and
#{j : gap_j inside S(P)} ~ k1 * |E_w(P)| by equidistribution, while the bad ones
number ~ k1 * (2/21) by codex's now-PROVED ceiling (THM-1203).  So a usable gap
exists once
    |E_w(P)|  >  2/21 = 0.095238...,
and the error terms in those two counts are exactly codex's supplier (2).  The raw
safe-set bound |S(P)| >= 0.164 is NOT enough, because erosion strips w off every
component and deletes outright every component shorter than w.

WHAT THIS SCRIPT COMPUTES, EXACTLY, IN RATIONALS.
 (a) S(P) for all 495 cores: components, lengths, count N(P), min |S(P)|.
 (b) |E_w(P)| as an exact piecewise-linear function of w.
 (c) The CRITICAL WIDTH w*(P): the largest w with |E_w(P)| > 2/21, solved exactly;
     and w* = min over the 495 cores, with the bottleneck core identified.
 (d) The resulting threshold K* = 6/(7 w*): for k1 > K* supplier (3) holds by
     measure alone, uniformly over all cores.
 (e) Comparison with the clustered-regime floor k1 > 13*max(P), to see how much of
     the range is left as a finite residue.
 (f) The crude bound |E| >= |S| - N*w, to show how much the exact computation gains.
"""
import sys
from fractions import Fraction as F
from itertools import combinations

BAD = F(2, 21)
LAM = F(1, 14)


def safe_components(P):
    """Exact components of S(P) = {t in [0,1] : ||p t|| >= 1/14 for all p in P}."""
    # candidate breakpoints: tooth edges of every comb
    bps = {F(0), F(1)}
    for p in P:
        for j in range(p + 1):
            for s in (F(1, 14 * p), -F(1, 14 * p)):
                v = F(j, p) + s
                if 0 <= v <= 1:
                    bps.add(v)
    B = sorted(bps)
    out = []
    for i in range(len(B) - 1):
        a, b = B[i], B[i + 1]
        if b <= a:
            continue
        mid = (a + b) / 2
        if all(min((p * mid) % 1, 1 - (p * mid) % 1) >= LAM for p in P):
            out.append((a, b))
    # merge touching components
    merged = []
    for a, b in out:
        if merged and a == merged[-1][1]:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return merged


def eroded_measure(lengths, w):
    return sum((L - w for L in lengths if L > w), F(0))


def critical_width(lengths, target=BAD):
    """Largest w with sum_j max(0, L_j - w) > target, solved exactly.
    On the interval where exactly the top m components survive,
    |E_w| = (sum of top m) - m*w, so w = ((sum of top m) - target)/m."""
    Ls = sorted(lengths, reverse=True)
    if sum(Ls, F(0)) <= target:
        return F(0)
    pref = F(0)
    best = F(0)
    for m in range(1, len(Ls) + 1):
        pref += Ls[m - 1]
        if pref <= target:
            continue
        w = (pref - target) / m
        # valid only if this w really leaves exactly the top m alive
        if w < Ls[m - 1] and (m == len(Ls) or w >= Ls[m] if m < len(Ls) else True):
            if w > best:
                best = w
    return best


CORES = [tuple(c) for c in combinations(range(1, 13), 8)]

print("=" * 78)
print("(a) THE CORE ATLAS: S(P) for all %d eight-speed cores P in {1,...,12}" % len(CORES))
print("=" * 78)
data = {}
for P in CORES:
    comps = safe_components(P)
    Ls = [b - a for a, b in comps]
    data[P] = Ls
minS = min((sum(Ls, F(0)), P) for P, Ls in data.items())
maxN = max((len(Ls), P) for P, Ls in data.items())
minN = min((len(Ls), P) for P, Ls in data.items())
maxLong = max((max(Ls), P) for P, Ls in data.items())
minLong = min((max(Ls), P) for P, Ls in data.items())
print("  min |S(P)| over all cores : %s = %.6f   at P = %s"
      % (minS[0], float(minS[0]), minS[1]))
print("     (the constant quoted in THM-1162 as '0.164 <= |S(P)|')")
print("  max |S(P)|                : %s = %.6f"
      % (max(sum(Ls, F(0)) for Ls in data.values()),
         float(max(sum(Ls, F(0)) for Ls in data.values()))))
print("  component count N(P)      : min %d (at %s), max %d (at %s)"
      % (minN[0], minN[1], maxN[0], maxN[1]))
print("  LONGEST component         : min over cores %s = %.6f at %s"
      % (minLong[0], float(minLong[0]), minLong[1]))
print("                              max over cores %s = %.6f"
      % (maxLong[0], float(maxLong[0])))
print("  -> if the longest component is <= w, E_w(P) is EMPTY: hard obstruction at")
print("     w >= %s, i.e. k1 <= 6/(7w) = %.2f"
      % (minLong[0], float(F(6, 7) / minLong[0])))

print()
print("=" * 78)
print("(c)(d) CRITICAL WIDTH  w*(P) = largest w with |E_w(P)| > 2/21, and K* = 6/(7w*)")
print("=" * 78)
crit = []
for P, Ls in data.items():
    w = critical_width(Ls)
    crit.append((w, P))
crit.sort()
wstar, Pstar = crit[0]
print("  BOTTLENECK core P* = %s" % (Pstar,))
print("     |S(P*)|            = %s = %.6f" % (sum(data[Pstar], F(0)), float(sum(data[Pstar], F(0)))))
print("     N(P*)              = %d components" % len(data[Pstar]))
print("     w*(P*)             = %s = %.8f" % (wstar, float(wstar)))
print("     K* = 6/(7 w*)      = %.3f" % (float(F(6, 7) / wstar) if wstar else float('inf')))
print()
print("  the five tightest cores:")
for w, P in crit[:5]:
    K = float(F(6, 7) / w) if w else float('inf')
    print("     P = %-32s w* = %-14s K* = %8.2f  N = %d  |S| = %.5f"
          % (str(P), str(w), K, len(data[P]), float(sum(data[P], F(0)))))
print()
print("  the five loosest cores:")
for w, P in crit[-5:]:
    K = float(F(6, 7) / w) if w else float('inf')
    print("     P = %-32s w* = %-14s K* = %8.2f  N = %d  |S| = %.5f"
          % (str(P), str(w), K, len(data[P]), float(sum(data[P], F(0)))))

KSTAR = F(6, 7) / wstar
print()
print("=" * 78)
print("(e) HOW MUCH OF THE RANGE DOES THIS CLOSE?")
print("=" * 78)
print("  UNIFORM MEASURE BOUND: for k1 > K* = %.3f, every core P has" % float(KSTAR))
print("     |E_w(P)| > 2/21, so a complete k1-gap survives inside S(P) outside the")
print("     bad set -- supplier (3) holds, uniformly over all 495 cores.")
print()
print("  clustered-regime floor  : k1 > 13*max(P) >= 13*12 = 156")
print("  so the FINITE RESIDUE is : 156 < k1 <= %d   (%d values of k1)"
      % (int(KSTAR), max(0, int(KSTAR) - 156)))
print("  per-core residue is smaller: the floor is 13*max(P), so a core with")
print("  max(P) = 12 needs k1 > 156, but a core with smaller max needs less.")

print()
print("=" * 78)
print("(f) HOW MUCH DOES THE EXACT EROSION GAIN OVER THE CRUDE |S| - N*w BOUND?")
print("=" * 78)
Ls = data[Pstar]
S = sum(Ls, F(0))
N = len(Ls)
wcrude = (S - BAD) / N
print("  at the bottleneck core P* = %s:" % (Pstar,))
print("     crude  w = (|S| - 2/21)/N = %s = %.8f  -> K_crude = %.2f"
      % (wcrude, float(wcrude), float(F(6, 7) / wcrude)))
print("     exact  w*                = %s = %.8f  -> K*      = %.2f"
      % (wstar, float(wstar), float(KSTAR)))
print("     exact/crude width ratio  = %.4f  (K* is %.2fx %s than crude)"
      % (float(wstar / wcrude), float(wcrude / wstar) if wstar > wcrude else float(wstar / wcrude),
         "SMALLER" if wstar > wcrude else "LARGER"))
print()
print("  component-length profile of the bottleneck core (top 12, exact):")
for L in sorted(Ls, reverse=True)[:12]:
    print("     %-16s = %.8f" % (str(L), float(L)))
print("     ... %d components total; %d of them shorter than w*"
      % (N, sum(1 for L in Ls if L <= wstar)))
