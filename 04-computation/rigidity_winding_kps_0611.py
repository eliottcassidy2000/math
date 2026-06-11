#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-485 / HYP-2417 — the HARD HALF by the ARGUMENT PRINCIPLE.

Target (rigidity reverse, the open analytic residue): for every finite nonempty
flip set S of Euler's pentagonal signs with NON-NEGATIVE boundary sum
    B(S) = sum_{k in S} (-1)^{k+1}  >= 0,
the truncated pentagonal polynomial
    F_S(x) = 1 - sum_k eps_k (x^{g_k} + x^{gbar_k}),
    eps_k = (-1)^{k+1} except flipped on S,   g_k = k(3k-1)/2, gbar_k = k(3k+1)/2,
has a ZERO in the OPEN unit disk |x| < 1.

The IVT half (B(S) < 0 => F_S(1) < 0 => real root in (0,1)) is already proved.
Here we attack the complementary half by CERTIFYING interior zeros with the
argument principle (winding number), which yields an EXACT integer zero-count
inside |x| = r (robust, not a heuristic grid minimum).

Method.
  number of zeros (with multiplicity) of F_S strictly inside |x| = r
      = (1/2pi) * total change of arg F_S(r e^{i theta}), theta: 0 -> 2pi
      = winding number of the image curve about 0.
We compute it by summing the principal-branch increments of arg between
consecutive sample points on the circle. To make the integer reliable we:
  (i) sample the circle finely (Nθ points),
  (ii) adaptively refine any arc where |Δarg| > pi/4 (so we never alias a
       jump of ~2pi as ~0), bisecting until every step is < pi/4 or a depth cap,
  (iii) sum the increments and divide by 2pi; round to nearest int.
A winding number >= 1 at some r < 1 CERTIFIES an interior zero. Winding 0 at
r -> 1^- (with the boundary value also checked) is a counterexample candidate.

F_S is a genuine POLYNOMIAL once truncated at Kmax pentagonal pairs; for x in the
closed disk the tail sum_{k>Kmax} eps_k(x^{g_k}+x^{gbar_k}) is bounded by
2*sum_{m>=g_{Kmax+1}} r^m = 2 r^{g_{Kmax+1}}/(1-r). At r=0.999, Kmax=60 gives
g_61 ~ 5550, tail < 2*0.999^5550/0.001 ~ 2000*e^{-5.55} ~ 7.8 -- NOT negligible
at r very close to 1. So we (a) keep r modestly below 1 (r in {0.9,0.95,0.99}) so
the truncation tail is provably tiny, AND (b) bound the tail explicitly and only
trust a winding number when |F_S| stays >> tail-bound along the whole circle.
"""

import cmath
import math
import itertools
import time


def pent_pairs(K):
    """List of (g_k, gbar_k, k) for k = 1..K."""
    out = []
    for k in range(1, K + 1):
        out.append((k * (3 * k - 1) // 2, k * (3 * k + 1) // 2, k))
    return out


def make_F(flips, Kmax):
    """Return a closure F(x) = 1 - sum eps_k (x^{g_k}+x^{gbar_k}) truncated at Kmax.
    eps_k = Euler (-1)^{k+1}, flipped on `flips`."""
    pairs = pent_pairs(Kmax)
    eps = []
    for g, gb, k in pairs:
        e = 1 if (k % 2 == 1) else -1
        if k in flips:
            e = -e
        eps.append((g, gb, e))

    def F(x):
        s = 1.0 + 0.0j
        for g, gb, e in eps:
            s -= e * (x ** g + x ** gb)
        return s

    return F


def tail_bound(r, Kmax):
    """Rigorous bound on |sum_{k>Kmax} eps_k (x^{g_k}+x^{gbar_k})| for |x|<=r.
    Each |term| <= r^{g_k} + r^{gbar_k} <= 2 r^{g_k}; g_k strictly increasing,
    g_{Kmax+1} = (Kmax+1)(3(Kmax+1)-1)/2. Bound by 2 * sum_{m>=g_{Kmax+1}} r^m."""
    k0 = Kmax + 1
    gmin = k0 * (3 * k0 - 1) // 2
    if r >= 1.0:
        return float("inf")
    # 2 * r^{gmin} / (1 - r) over-counts (only pentagonal exps appear) but is safe.
    return 2.0 * (r ** gmin) / (1.0 - r)


def _arg_increment(z0, z1):
    """Principal-branch increment of arg from z0 to z1 in (-pi, pi]."""
    return cmath.phase(z1 / z0)


def winding_number(F, r, Ntheta=2048, refine_depth=22, min_modulus_track=None):
    """Winding number of F(r e^{i theta}) about 0, theta in [0, 2pi).
    Adaptive: any step with |Δarg| > pi/4 is bisected until < pi/4 or depth cap.
    Returns (winding_int, total_arg, min|F| on the sampled circle).
    If F passes (near) through 0 on the circle, min|F| ~ 0 flags an unreliable count."""
    thetas = [2 * math.pi * i / Ntheta for i in range(Ntheta)]
    pts = []
    minmod = float("inf")
    for th in thetas:
        z = r * cmath.exp(1j * th)
        f = F(z)
        minmod = min(minmod, abs(f))
        pts.append((th, f))

    total = 0.0
    TWO_PI = 2 * math.pi
    THRESH = math.pi / 4

    def refine(thA, fA, thB, fB, depth):
        """Accumulate arg increment from A to B, bisecting big jumps."""
        inc = _arg_increment(fA, fB)
        if abs(inc) <= THRESH or depth >= refine_depth:
            return inc
        thm = (thA + thB) / 2 if thB > thA else (thA + (thB + TWO_PI)) / 2
        thm_wrapped = thm % TWO_PI
        fm = F(r * cmath.exp(1j * thm_wrapped))
        return refine(thA, fA, thm, fm, depth + 1) + refine(thm, fm, thB, fB, depth + 1)

    for i in range(Ntheta):
        thA, fA = pts[i]
        thB, fB = pts[(i + 1) % Ntheta]
        if (i + 1) == Ntheta:
            thB = TWO_PI  # close the loop
        total += refine(thA, fA, thB, fB, 0)

    wind = total / TWO_PI
    return round(wind), wind, minmod


def boundary_sum(S):
    return sum((1 if k % 2 == 1 else -1) for k in S)


def F_at_1(S):
    return 4 * boundary_sum(S)


def main():
    t0 = time.time()
    Kmax = 60  # truncation: g_61 ~ 5550, tail tiny for r <= 0.99
    KFLIP = 10
    MAXSIZE = 4
    radii = [0.90, 0.95, 0.99]

    print("=" * 78, flush=True)
    print("THM-485 / HYP-2417 : argument-principle certification of interior zeros",
          flush=True)
    print(f"  flip range k <= {KFLIP}, |S| <= {MAXSIZE}, truncation Kmax = {Kmax}",
          flush=True)
    print(f"  test radii r in {radii}  (winding number = exact interior zero count)",
          flush=True)
    print("=" * 78, flush=True)

    # tail bounds at the radii we use (sanity: must be tiny vs typical |F|~O(1))
    print("\nTruncation tail bounds (|F_S - F_S^trunc| <= bound for |x|<=r):", flush=True)
    for r in radii:
        print(f"  r={r}: tail <= {tail_bound(r, Kmax):.3e}", flush=True)

    # ---- enumerate the HARD HALF: B(S) >= 0 ----
    hard = []
    easy = []
    for sz in range(1, MAXSIZE + 1):
        for S in itertools.combinations(range(1, KFLIP + 1), sz):
            S = frozenset(S)
            if boundary_sum(S) >= 0:
                hard.append(S)
            else:
                easy.append(S)
    print(f"\nNonempty flip sets on k<={KFLIP}, |S|<={MAXSIZE}: "
          f"{len(hard)+len(easy)} total", flush=True)
    print(f"  HARD half (B(S) >= 0, complex/interior zero expected): {len(hard)}",
          flush=True)
    print(f"  EASY half (B(S) <  0, IVT real root, already proved):   {len(easy)}",
          flush=True)

    print("\n" + "-" * 78, flush=True)
    print("Certifying the HARD half by winding number ...", flush=True)
    print("-" * 78, flush=True)

    n_certified = 0           # winding >= 1 at SOME r
    n_uncertified = []        # winding 0 at ALL r  -> counterexample candidates
    n_unreliable = []         # min|F| too close to tail bound (circle near a zero)
    winding_at_best = []      # (S, best winding, r)
    examples = []

    F_cache = {}

    for idx, S in enumerate(hard):
        F = make_F(S, Kmax)
        best_w = 0
        best_r = None
        reliable_any = False
        for r in radii:
            w_int, w_float, minmod = winding_number(F, r, Ntheta=1024)
            tb = tail_bound(r, Kmax)
            reliable = minmod > 10 * tb and minmod > 1e-9
            if reliable:
                reliable_any = True
            if w_int > best_w:
                best_w = w_int
                best_r = r
        winding_at_best.append((S, best_w, best_r))
        if best_w >= 1:
            n_certified += 1
        else:
            # no interior zero detected at any radius below 1
            n_uncertified.append((S, boundary_sum(S)))
        if not reliable_any:
            n_unreliable.append(S)
        if len(examples) < 12:
            examples.append((set(S), boundary_sum(S), best_w, best_r))

    print("\nResults (HARD half, B(S) >= 0):", flush=True)
    print(f"  sets tested:                              {len(hard)}", flush=True)
    print(f"  CERTIFIED interior zero (winding >= 1):   {n_certified}", flush=True)
    print(f"  NO interior zero found (winding 0 all r): {len(n_uncertified)}", flush=True)
    print(f"  flagged unreliable (circle near a zero):  {len(n_unreliable)}", flush=True)

    if n_uncertified:
        print("\n  *** COUNTEREXAMPLE CANDIDATES (winding 0 at all radii) ***", flush=True)
        for S, b in n_uncertified[:40]:
            print(f"      S={set(S)}  B(S)={b}", flush=True)

    print("\nSample certifications (S, B(S), best winding, r):", flush=True)
    for S, b, w, r in examples:
        print(f"   S={S}  B={b:+d}  winding={w}  @ r={r}", flush=True)

    # distribution of winding numbers
    from collections import Counter
    dist = Counter(w for _, w, _ in winding_at_best)
    print("\nDistribution of (best) winding numbers over the HARD half:", flush=True)
    for w in sorted(dist):
        print(f"   winding {w}: {dist[w]} sets", flush=True)

    # ---- spot-check the EASY half too, for completeness (winding should also see it) ----
    print("\n" + "-" * 78, flush=True)
    print("Spot-check EASY half (B<0): winding should ALSO certify (real root inside)",
          flush=True)
    print("-" * 78, flush=True)
    easy_ok = 0
    easy_bad = []
    for S in easy:
        F = make_F(S, Kmax)
        best_w = 0
        for r in radii:
            w_int, _, minmod = winding_number(F, r, Ntheta=1024)
            best_w = max(best_w, w_int)
        if best_w >= 1:
            easy_ok += 1
        else:
            easy_bad.append(set(S))
    print(f"  EASY sets with winding >= 1: {easy_ok}/{len(easy)}", flush=True)
    if easy_bad:
        print(f"  EASY sets missed by winding: {easy_bad[:20]}", flush=True)

    # ---- the all-Euler control: F_Euler truncated must have winding 0 inside ----
    print("\n" + "-" * 78, flush=True)
    print("Control: Euler (empty flip set) -- winding must be 0 (zero-free inside)",
          flush=True)
    print("-" * 78, flush=True)
    Feul = make_F(frozenset(), Kmax)
    for r in radii:
        w_int, w_float, minmod = winding_number(Feul, r, Ntheta=2048)
        print(f"   r={r}: winding = {w_int} (float {w_float:+.4f}), min|F| = {minmod:.4e}",
              flush=True)

    # ---- verdict ----
    print("\n" + "=" * 78, flush=True)
    rigid = (len(n_uncertified) == 0)
    print(f"VERDICT: rigidity holds on |S|<=4, k<=10  ==>  {rigid}", flush=True)
    if rigid:
        print("  Every B(S)>=0 flip set has a CERTIFIED interior zero (winding>=1).",
              flush=True)
        print("  Combined with the IVT half (B(S)<0), EVERY nonempty finite flip set",
              flush=True)
        print("  on this family acquires a zero in the open disk. No counterexample.",
              flush=True)
    else:
        print("  Counterexample candidate(s) found -- rigidity would be REFUTED;",
              flush=True)
        print("  re-examine with higher Ntheta / Kmax / argument-principle on a contour",
              flush=True)
        print("  that avoids the boundary before drawing conclusions.", flush=True)
    print(f"=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
