#!/usr/bin/env python3
"""rigidity_Sj_stress_kps_S128c106.py -- kind-pasteur-2026-07-20-S128c106

STRESS TEST OF THM-488 THEOREM A ON ITS OWN NEAR-BOUNDARY FAMILY.

THM-488 Theorem A: for a sign sequence eps differing from Euler's (-1)^{k+1} in
finitely many places S, F_eps is zero-free in the open unit disk IFF S is empty.
The IVT half handles B(S) = sum_{k in S} (-1)^{k+1} < 0.  The HARD half, B(S) >= 0,
is CERTIFIED by the argument principle on |S| <= 6, k <= 12 (1585/1585).

THE GAP.  THM-488 itself names its four near-boundary sets -- {1,3,5,7,9},
{1,3,5,7,9,10}, {1,3,5,7,9,11}, {1,3,5,7,9,12} -- and notes that flipping signs on
the smallest odd k "pushes the interior zero toward |x| = 1 but never onto or
outside it".  Every one of those is an initial segment of

    S_j = {1, 3, 5, ..., 2j-1},     B(S_j) = j > 0   (so squarely in the hard half)

and the certification stops at |S| <= 6, i.e. j <= 6 -- exactly where its own
diagnostic says the danger starts.  If |z_min(S_j)| -> 1 as j grows, Theorem A is
FALSE and the 1585/1585 certificate was a small-|S| artifact.

METHOD.  F_S(x) = 1 - sum_{k<=Kmax} eps_k (x^{g_k} + x^{gbar_k}), g_k = k(3k-1)/2,
gbar_k = k(3k+1)/2, eps_k = -(-1)^{k+1} for k in S and (-1)^{k+1} otherwise.  Only
2*Kmax terms, so evaluation is cheap at any Kmax.  Interior zeros are counted by
the ARGUMENT PRINCIPLE: the winding number of F_S(r e^{i theta}) about 0 equals the
number of zeros inside |x| = r.  Argument increments are accumulated with adaptive
refinement so no step exceeds pi/4 (no 2*pi aliasing).

TRUNCATION DISCIPLINE (THM-488's own, reused).  The truncated series approximates
Pi(1-x^n), whose zeros sit ON |x| = 1 and leak inside as r -> 1.  So:
  * the SAFE ZONE is r <= 0.95, where the truncated Euler series must have winding
    exactly 0 -- this is asserted and CHECKED here, not assumed;
  * above it we report the winding DIFFERENCE D(r) = w(F_S,r) - w(F_Euler,r), which
    cancels the shared artifact zeros since F_S and F_Euler differ only in finitely
    many low coefficients.
A winding >= 1 at r <= 0.95 certifies a genuine interior zero.  A winding of 0 at
every safe r means the zero (if any) has been pushed into the boundary annulus, and
D(r) is then the only honest instrument.
"""
import sys
import math
import cmath

JMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
KBUF = int(sys.argv[2]) if len(sys.argv) > 2 else 12


def pent(k):
    return k * (3 * k - 1) // 2, k * (3 * k + 1) // 2


def make_F(S, Kmax):
    """Return (exponent, coefficient) list for F_S truncated at index Kmax."""
    terms = [(0, 1.0)]
    for k in range(1, Kmax + 1):
        e = (-1) ** (k + 1)
        if k in S:
            e = -e
        g, gb = pent(k)
        terms.append((g, -float(e)))
        terms.append((gb, -float(e)))
    return terms


def evalF(terms, z):
    s = 0j
    for e, c in terms:
        s += c * (z ** e)
    return s


def winding(terms, r, base=2048, maxdepth=22):
    """Exact integer winding number of F(r e^{i t}) about 0, by accumulated
    argument with adaptive bisection keeping every step below pi/4."""
    def val(t):
        return evalF(terms, cmath.rect(r, t))

    ts = [2 * math.pi * i / base for i in range(base + 1)]
    vs = [val(t) for t in ts]
    if any(abs(v) == 0.0 for v in vs):
        return None
    total = 0.0
    for i in range(len(ts) - 1):
        a, b = ts[i], ts[i + 1]
        va, vb = vs[i], vs[i + 1]
        stack = [(a, b, va, vb, 0)]
        while stack:
            a2, b2, va2, vb2, d = stack.pop()
            d_arg = cmath.phase(vb2 / va2)
            if abs(d_arg) < math.pi / 4 or d >= maxdepth:
                total += d_arg
                continue
            m = 0.5 * (a2 + b2)
            vm = val(m)
            if vm == 0:
                return None
            stack.append((m, b2, vm, vb2, d + 1))
            stack.append((a2, m, va2, vm, d + 1))
    return int(round(total / (2 * math.pi)))


print("=" * 78)
print("STEP 0 -- truncation safety: winding of the TRUNCATED EULER series must be 0")
print("          on the safe zone r <= 0.95 (THM-488's condition, re-checked)")
print("=" * 78)
for Kmax in (20, 40, 80):
    terms = make_F(set(), Kmax)
    row = []
    for r in (0.80, 0.90, 0.95):
        row.append((r, winding(terms, r)))
    print("  Kmax = %-4d " % Kmax + "  ".join("w(%.2f) = %s" % (r, w) for r, w in row))
print("  -> safe zone confirmed where all these read 0.")

print()
print("=" * 78)
print("STEP 1 -- THE FAMILY  S_j = {1,3,...,2j-1},  j = 1..%d" % JMAX)
print("=" * 78)
print("  %-4s %-6s %-8s %-10s %-10s %-10s %-12s"
      % ("j", "|S|", "B(S)", "w(0.80)", "w(0.90)", "w(0.95)", "verdict"))
first_fail = None
rows = []
for j in range(1, JMAX + 1):
    S = set(range(1, 2 * j, 2))
    Kmax = max(2 * j - 1 + KBUF, 20)
    terms = make_F(S, Kmax)
    w80 = winding(terms, 0.80)
    w90 = winding(terms, 0.90)
    w95 = winding(terms, 0.95)
    safe_zero = max(x for x in (w80, w90, w95) if x is not None)
    if safe_zero >= 1:
        verdict = "interior zero"
    else:
        verdict = "NO SAFE ZERO"
        if first_fail is None:
            first_fail = j
    rows.append((j, len(S), j, w80, w90, w95, verdict, Kmax))
    print("  %-4d %-6d %-8d %-10s %-10s %-10s %-12s"
          % (j, len(S), j, w80, w90, w95, verdict))

print()
print("=" * 78)
print("STEP 2 -- THE BOUNDARY ANNULUS for any j with no certified safe-zone zero")
print("=" * 78)
if first_fail is None:
    print("  none: every j = 1..%d carries a certified interior zero at r <= 0.95," % JMAX)
    print("  so Theorem A SURVIVES on its own near-boundary family, far past the")
    print("  |S| <= 6 certification.  The 'pushes the zero toward |x| = 1' diagnostic")
    print("  does not become a counterexample in this direction.")
else:
    print("  first j with no safe-zone zero : %d" % first_fail)
    print("  winding DIFFERENCE D(r) = w(F_S,r) - w(F_Euler,r) in the annulus:")
    for j in range(first_fail, min(first_fail + 6, JMAX + 1)):
        S = set(range(1, 2 * j, 2))
        Kmax = max(2 * j - 1 + KBUF, 20)
        tS = make_F(S, Kmax)
        tE = make_F(set(), Kmax)
        ds = []
        for r in (0.96, 0.97, 0.98, 0.99):
            wS, wE = winding(tS, r), winding(tE, r)
            ds.append((r, None if (wS is None or wE is None) else wS - wE))
        print("     j = %-3d " % j + "  ".join("D(%.2f) = %s" % (r, d) for r, d in ds))
    print()
    print("  D(r) >= 1 anywhere  => the zero merely moved into the annulus, Theorem A")
    print("                         survives (the safe zone was just too small).")
    print("  D(r) = 0 for all r  => a genuine COUNTEREXAMPLE to Theorem A's hard half:")
    print("                         F_S zero-free in the disk with S non-empty.")

print()
print("=" * 78)
print("STEP 3 -- where is the zero?  smallest safe radius carrying a zero")
print("=" * 78)
for j in list(range(1, min(13, JMAX + 1))) + [x for x in (20, 30, 40) if x <= JMAX]:
    S = set(range(1, 2 * j, 2))
    Kmax = max(2 * j - 1 + KBUF, 20)
    terms = make_F(S, Kmax)
    lo, hi = 0.02, 0.95
    if (winding(terms, hi) or 0) < 1:
        print("  j = %-3d : no zero inside 0.95" % j)
        continue
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if (winding(terms, mid) or 0) >= 1:
            hi = mid
        else:
            lo = mid
    print("  j = %-3d : |z_min| ~ %.6f   (B = %d, Kmax = %d)" % (j, hi, j, Kmax))
print()
print("  If |z_min| is INCREASING toward 1 with j, Theorem A is under real strain")
print("  and the next j is worth computing.  If it stabilises or decreases, the")
print("  near-boundary diagnostic in THM-488 was a small-|S| effect.")
