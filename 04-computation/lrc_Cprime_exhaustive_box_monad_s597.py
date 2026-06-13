#!/usr/bin/env python3
"""EXHAUSTIVE completeness certificate for C' inside a bounded box
(monad-compute-2026-06-03-S597).

THM-398 / HYP-2102 reduce LRC(n) to:

    C' : if a PRIMITIVE speed set S (gcd = 1, the n-1 moving runners) contains a
         multiple of n, then S is LOOSE, i.e. M(S) = max_t min_{v in S} ||v t|| > 1/n.

Prior monad-compute work (S595/S596) and S571 verified C' by SAMPLING the
multiple-of-n class through n=20 (random + hardest small-companion slices). That
is strong evidence but not a completeness statement: it cannot certify that *no*
tight-with-multiple config exists inside any region.

This script gives a genuinely EXHAUSTIVE certificate. For each n it enumerates
*every* primitive (n-1)-subset of {1,...,B} that contains a multiple of n
(gcd of the whole set = 1, since M is scale-invariant so primitive reps suffice
and avoid scaled duplicates), and tests each one EXACTLY for looseness. If none
is tight, C' holds with NO EXCEPTIONS inside the box [1,B]^{n-1} -- a stronger
statement than any sample.

WHY primitive: M(cS) = M(S) for any scalar c (substitute t -> t/c), so tightness
is scale-invariant; the LRC is stated for gcd-1 speed sets. The class "primitive
+ contains a multiple of n" is exactly the C' hypothesis class (matches the
prim()+any(%n==0) filter used by S595/S596). Note nS (all multiples of n) has
prim form S, which need NOT contain a multiple of n -- so this class is a genuine
proper subclass, not all of LRC.

LOOSENESS TEST (exact, rational, identical method to S595/S596):
  S is loose  <=>  the OPEN safe set  A = { t in [0,1) : ||v t|| > 1/n  for all v }
                   has positive measure.
The indicator of A is piecewise constant with breakpoints exactly at
t = (k*n +- 1)/(n*v), v in S. Each elementary interval is decided by its midpoint
via a STRICT inequality (open set). All arithmetic is fractions.Fraction -- no
floats anywhere (cf. MISTAKE-019 on overflow / float pitfalls).

We compute the FULL exact measure (not just positivity) so we can also report the
minimum looseness margin over the whole box -- the smallest mu(A) > 0 -- which is
the empirical "closest approach to tight" inside the box. A tight config would
have mu(A) = 0 exactly and would be printed loudly as a C' counterexample.

No proof; pure computation.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools
import sys


def dist_strict_gt(v, mid, n):
    """Return True iff ||v*mid|| > 1/n, exactly (mid a Fraction in [0,1))."""
    x = (v * mid) % 1
    d = x if x < 1 - x else 1 - x
    return n * d > 1  # ||v mid|| > 1/n  <=>  n*d > 1


def _breakpoints(V, n):
    eps = {F(0)}
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                eps.add(F(k * n + s, n * v) % 1)
    return sorted(eps)


def open_safe_measure(V, n):
    """Exact measure of { t in [0,1) : min_i ||v_i t|| > 1/n } (STRICT / open)."""
    pts = _breakpoints(V, n)
    meas = F(0)
    L = len(pts)
    for i in range(L):
        a = pts[i]
        b = pts[i + 1] if i + 1 < L else pts[0] + 1
        ln = b - a
        mid = (a + ln / 2) % 1
        if all(dist_strict_gt(v, mid, n) for v in V):
            meas += ln
    return meas


def is_loose(V, n):
    """True iff S is LOOSE (open safe set nonempty), with EARLY EXIT at the first
    safe interval. Exact; rigorously equivalent to open_safe_measure(V,n) > 0
    because both scan the same elementary intervals with the same strict midpoint
    test -- this just stops at the first one that passes."""
    pts = _breakpoints(V, n)
    L = len(pts)
    for i in range(L):
        a = pts[i]
        b = pts[i + 1] if i + 1 < L else pts[0] + 1
        mid = (a + (b - a) / 2) % 1
        if all(dist_strict_gt(v, mid, n) for v in V):
            return True
    return False


def set_gcd(V):
    return reduce(gcd, V)


def enumerate_box(n, B):
    """Yield every primitive (n-1)-subset of {1..B} containing a multiple of n.

    'Primitive' = gcd of the whole set is 1 (scale-invariance of M). Yielding only
    primitive reps avoids counting cS as a distinct config from S.
    """
    m = n - 1
    multiples = [x for x in range(n, B + 1, n)]  # n, 2n, ... <= B
    mset = set(multiples)
    # Enumerate all m-subsets that contain >=1 multiple of n, primitive, no dups.
    # Strategy: for each subset, require it meets `mset` and gcd==1.
    for combo in itertools.combinations(range(1, B + 1), m):
        if not (mset & set(combo)):
            continue
        if set_gcd(combo) != 1:
            continue
        yield combo


def check_n(n, B, exact_margin, report_every=50000):
    """Exhaustively certify every primitive multiple-of-n config in the box is
    loose. `exact_margin`: if True, also compute the full open-safe measure of
    every config to report the minimum looseness margin (only affordable for
    small boxes). Otherwise use early-exit positivity (is_loose) -- still a fully
    rigorous certificate, just without the margin statistic."""
    total = 0
    loose = 0
    tight = 0
    min_meas = None
    min_cfg = None
    tights = []
    for V in enumerate_box(n, B):
        total += 1
        if exact_margin:
            mu = open_safe_measure(V, n)
            ok = mu > 0
            if ok and (min_meas is None or mu < min_meas):
                min_meas, min_cfg = mu, V
        else:
            ok = is_loose(V, n)
        if ok:
            loose += 1
        else:
            tight += 1
            tights.append(V)
        if total % report_every == 0:
            print(f"    ... n={n} B={B}: scanned {total} configs "
                  f"({loose} loose, {tight} tight) so far", flush=True)
    return total, loose, tight, min_meas, min_cfg, tights


# Box sizes B = K*n. Chosen for exhaustive feasibility with exact arithmetic.
# K is the largest companion-to-n ratio covered (the C' residual lives at small
# multiples / companions, so a few multiples of n in each direction is the
# interesting band).
BOXES = [
    (4, 40),   # K=10  C(40,3) base
    (5, 40),   # K=8
    (6, 36),   # K=6
    (7, 28),   # K=4
    (8, 24),   # K=3
]


def main():
    print("=" * 72)
    print("C' EXHAUSTIVE BOX CERTIFICATE (monad-compute-S597)")
    print("Claim C': primitive speed set with a multiple of n  =>  LOOSE (M>1/n).")
    print("Method: enumerate EVERY primitive (n-1)-subset of {1..B} containing a")
    print("multiple of n; test looseness EXACTLY via open-safe-set measure > 0.")
    print("A single tight config inside the box would REFUTE C'.")
    print("=" * 72)
    print()
    print(f"{'n':>3} {'B':>4} {'K=B/n':>6} {'#configs':>10} {'loose':>9} "
          f"{'tight':>6} {'cert':>6}   min_margin (smallest mu>0)")
    grand_total = grand_loose = grand_tight = 0
    all_tights = []
    # exact min-margin only affordable on the small boxes (full measure on all
    # configs); larger boxes use rigorous early-exit positivity certificate.
    EXACT_MARGIN_N = {4, 5}
    for n, B in BOXES:
        total, loose, tight, min_meas, min_cfg, tights = check_n(
            n, B, exact_margin=(n in EXACT_MARGIN_N))
        cert = "PASS" if tight == 0 else "FAIL"
        if min_meas is not None:
            mm = f"{float(min_meas):.6f}  (={min_meas})  at {min_cfg}"
        else:
            mm = "(positivity-only certificate)"
        print(f"{n:>3} {B:>4} {B // n:>6} {total:>10} {loose:>9} "
              f"{tight:>6} {cert:>6}   {mm}", flush=True)
        grand_total += total
        grand_loose += loose
        grand_tight += tight
        all_tights += [(n, B, t) for t in tights]

    print()
    print("-" * 72)
    print(f"GRAND TOTAL: {grand_total} primitive multiple-of-n configs enumerated "
          f"exhaustively across the boxes above.")
    print(f"             {grand_loose} loose, {grand_tight} tight.")
    if grand_tight == 0:
        print()
        print("CERTIFICATE: C' holds with ZERO EXCEPTIONS over every primitive")
        print("multiple-of-n speed set inside the listed boxes. This is an")
        print("EXHAUSTIVE (not sampled) confirmation of THM-398/HYP-2102's C'")
        print("hypothesis within these finite regions.")
    else:
        print()
        print("*** C' COUNTEREXAMPLE(S) FOUND -- would REFUTE C' / LRC reduction ***")
        for n, B, t in all_tights[:50]:
            print(f"    n={n} B={B}  TIGHT (M=1/n) config: {t}")
    sys.stdout.flush()


if __name__ == '__main__':
    main()
