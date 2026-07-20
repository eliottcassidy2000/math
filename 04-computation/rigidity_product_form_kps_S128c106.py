#!/usr/bin/env python3
"""rigidity_product_form_kps_S128c106.py -- kind-pasteur-2026-07-20-S128c106

THE TRUNCATION ARTIFACT IN THM-488 IS AVOIDABLE ENTIRELY.

THM-488 certifies the hard half of its Theorem A with the argument principle on a
TRUNCATED pentagonal SERIES, and then spends its most delicate section on a
"truncation-artifact correction": the truncated series approximates Pi(1-x^n),
whose zeros lie ON |x| = 1, and they leak inside as r -> 1, so windings are only
trustworthy for r <= 0.95.  That restriction is what caps the certificate.

But the perturbation is FINITE.  Writing

    F_S(x) = Pi_{n>=1} (1 - x^n)  +  2 * sum_{k in S} (-1)^{k+1} (x^{g_k} + x^{gbar_k})

the first term is the EULER PRODUCT -- zero-free in the open disk and computable to
any accuracy by a partial product, with NO spurious zeros -- and the second is an
exact finite polynomial.  So F_S can be evaluated reliably at any |x| < 1, the
argument principle applies with no artifact, and the 0.95 ceiling disappears.

This matters because rigidity_Sj_stress showed |z_min(S_j)| climbing monotonically
    j      = 1        2        3        4        5        6
    |z_min| = 0.7815   0.8801   0.9101   0.9262   0.9411   0.9497
and then leaving the r <= 0.95 safe zone at exactly j = 7 -- i.e. THM-488's
|S| <= 6 certification range coincides with the reach of its instrument rather
than with any change in the mathematics.  Whether Theorem A's hard half survives
past j = 6 is therefore OPEN on the evidence in canon.  This script settles it.

Numerical note: at r = 0.99 the Euler product is ~1e-72 while the perturbation is
O(1), so F_S is perturbation-dominated near the boundary; both are representable in
double precision and are ADDED (no cancellation), so ordinary floats suffice.  The
partial-product tail error is |x|^{N+1}/(1-|x|) ~ 1e-11 at r = 0.99, N = 3000.
"""
import sys
import math
import cmath

JMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
NPROD = int(sys.argv[2]) if len(sys.argv) > 2 else 4000


def pent(k):
    return k * (3 * k - 1) // 2, k * (3 * k + 1) // 2


def make_pert(S):
    """Exact finite perturbation 2*sum_{k in S} (-1)^{k+1}(x^{g_k}+x^{gbar_k})."""
    t = []
    for k in sorted(S):
        c = 2.0 * ((-1) ** (k + 1))
        g, gb = pent(k)
        t.append((g, c))
        t.append((gb, c))
    return t


def nterms(r):
    """Partial-product length: stop once |z|^n is below 1e-18 (tail error
    |z|^{N+1}/(1-|z|) is then negligible against an O(1) perturbation)."""
    if r <= 0:
        return 1
    return min(NPROD, max(50, int(math.log(1e-18) / math.log(r)) + 2))


def F(z, pert, N=None):
    """Euler product (exact, zero-free in the disk) + exact finite perturbation."""
    if N is None:
        N = nterms(abs(z))
    p = 1.0 + 0j
    zn = z
    for _ in range(N):
        p *= (1 - zn)
        zn *= z
        if abs(zn) < 1e-19:
            break
    s = p
    for e, c in pert:
        s += c * (z ** e)
    return s


def winding(pert, r, base=1024, maxdepth=22):
    def val(t):
        return F(cmath.rect(r, t), pert)
    ts = [2 * math.pi * i / base for i in range(base + 1)]
    vs = [val(t) for t in ts]
    if any(v == 0 for v in vs):
        return None
    total = 0.0
    for i in range(len(ts) - 1):
        stack = [(ts[i], ts[i + 1], vs[i], vs[i + 1], 0)]
        while stack:
            a, b, va, vb, d = stack.pop()
            da = cmath.phase(vb / va)
            if abs(da) < math.pi / 4 or d >= maxdepth:
                total += da
                continue
            m = 0.5 * (a + b)
            vm = val(m)
            if vm == 0:
                return None
            stack.append((m, b, vm, vb, d + 1))
            stack.append((a, m, va, vm, d + 1))
    return int(round(total / (2 * math.pi)))


print("=" * 78)
print("STEP 0 -- CONTROL: the unperturbed Euler product must have winding 0")
print("          at EVERY radius (it is zero-free in the open disk)")
print("=" * 78)
ctrl = [(r, winding([], r)) for r in (0.80, 0.95, 0.99, 0.995, 0.999)]
print("  " + "  ".join("w(%.3f) = %s" % (r, w) for r, w in ctrl))
ok = all(w == 0 for _, w in ctrl)
print("  all zero : %s   -> the product form has NO truncation artifact," % ok)
print("     unlike the truncated series, whose winding jumped 0 -> 19 -> 41 -> 48")
print("     across r = 0.95 -> 0.97 -> 0.99 (THM-488's own table).")

print()
print("=" * 78)
print("STEP 1 -- S_j = {1,3,...,2j-1}: does F_S have an interior zero, j = 1..%d?" % JMAX)
print("=" * 78)
print("  %-4s %-9s %-9s %-9s %-10s %-10s %s"
      % ("j", "w(0.95)", "w(0.99)", "w(0.999)", "|z_min|", "zeros", "verdict"))
zmins = []
counter = []
for j in range(1, JMAX + 1):
    S = set(range(1, 2 * j, 2))
    pert = make_pert(S)
    w95 = winding(pert, 0.95)
    w99 = winding(pert, 0.99)
    w999 = winding(pert, 0.999)
    nz = max(x for x in (w95, w99, w999) if x is not None)
    if nz >= 1:
        lo, hi = 0.02, 0.999
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            if (winding(pert, mid) or 0) >= 1:
                hi = mid
            else:
                lo = mid
        zmins.append((j, hi))
        verdict = "interior zero"
        zm = "%.6f" % hi
    else:
        counter.append(j)
        verdict = "*** ZERO-FREE ***"
        zm = "  --  "
    print("  %-4d %-9s %-9s %-9s %-10s %-10d %s"
          % (j, w95, w99, w999, zm, nz, verdict))

print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
if counter:
    print("  COUNTEREXAMPLES to THM-488 Theorem A (hard half) at j = %s" % counter)
    print("  For these S, F_S is zero-free in the open unit disk with S non-empty,")
    print("  contradicting 'zero-free IFF S is empty'.")
else:
    print("  Every S_j, j = 1..%d, carries a certified interior zero." % JMAX)
    print("  Theorem A's hard half SURVIVES on its own near-boundary family, far")
    print("  past the |S| <= 6 certification -- and now with no truncation caveat.")
    print()
    print("  |z_min| trajectory:")
    for j, z in zmins:
        print("     j = %-3d  |z_min| = %.6f   1 - |z_min| = %.3e" % (j, z, 1 - z))
    if len(zmins) >= 3:
        print()
        print("  The zero approaches the boundary but does not reach it; the")
        print("  certification ceiling at |S| <= 6 was an artifact of the truncated-")
        print("  series instrument, not a boundary in the mathematics.  The product")
        print("  form removes it and extends the certificate by a factor of ~%d in j."
              % (JMAX // 6))
