#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_resonant_stranger_certificate_macmini_0615s6.py  (mac-mini-2026-06-15-S6)

ANGLE 5 — RESONANT-STRANGER ENUMERATION + UNIFORM PER-CORE BOUND.
Converts the INFINITE m-direction of the interior-drop extremizer family into a
FINITE, bounded check, with EXACT rational arithmetic.

THE FAMILY.  The numerically-observed extremizers of LRC(14) are the interior-drop
cores  A_j ∪ {14m},  A_j = {1,...,13}\{j}  (j=1..13),  stranger 14m  (m=1,2,3,...).
L(S) = meas{ tau in [0,1) : ||v tau|| > 1/14  for all v in S }.

EXACT L.  For the 12-speed core A_j (all speeds <= 13), Lonely(A_j) is a finite
union of r_j maximal open "safe" intervals computed exactly with Fractions.  Adding
the stranger w=14m removes from each safe interval the points with ||w tau|| <= 1/14
(the union of "danger bands" of width 2/(14w)=1/(7w), period 1/w).  Exact interval
subtraction gives L(A_j ∪ {14m}) as an exact rational.  (Cross-checked against the
brute breakpoint-sweep: identical rationals.)

THE TAIL BOUND (PROVED, the key to finiteness).  On a safe interval I of length ell,
the stranger's danger set meets at most floor(w*ell)+1 period-cells, each of measure
2/(14w)=1/(7w), so it removes at most (w*ell+1)/(7w) = ell/7 + 1/(7w).  Summing over
the r_j intervals:
        L(A_j ∪ {14m})  >=  (6/7)*meas(Lonely(A_j))  -  r_j/(7w)
                         =   limit_j  -  r_j/(98 m).
Hence the dip below the limit is AT MOST r_j/(98 m): it -> 0 as m -> inf (Weyl limit
recovered), and L(A_j∪{14m}) > c_j is GUARANTEED for every m > M_j := r_j/(98(limit_j-c_j)),
where c_j is the min over the finite window.  So min over ALL m is attained at finite m
and we verify it exhaustively up to max(window, M_j).

RESULT (exact):
  inf over the whole interior-drop family = min_j c_j = 1543/294294 ~= 0.0052431 > 0,
  uniquely at j=6, m=7, i.e. S = {1,2,3,4,5,7,8,9,10,11,12,13,98}  (98 = 2*7^2).
  Every per-core floor c_j is a strictly positive rational; the second-deepest is
  c_12 = 563/105105 ~= 0.0053565 at m=6 (84 = 2^2*3*7).
This is exactly the previously-observed numerical inf, now pinned as an EXACT rational
and certified positive over the entire infinite family by a finite + tail argument.
"""
from fractions import Fraction as Fr
from math import floor, ceil
import sys, time
sys.stdout.reconfigure(line_buffering=True)

n = 14
thr = Fr(1, n)          # danger half-threshold 1/14

def core_safe_intervals(core):
    """Maximal open safe intervals of Lonely(core) in [0,1); core speeds all small."""
    bps = set([Fr(0), Fr(1)])
    for v in core:
        for k in range(0, v + 1):
            for s in (-thr, thr):
                x = Fr(k) + s
                if 0 <= x <= v:
                    bps.add(x / v)
    bps = sorted(b for b in bps if 0 <= b <= 1)
    ivs = []
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        ok = True
        for v in core:
            x = v * mid
            f = x - (x.numerator // x.denominator)
            d = f if f <= Fr(1, 2) else 1 - f
            if d <= thr:
                ok = False
                break
        if ok:
            ivs.append((lo, hi))
    merged = []
    for lo, hi in ivs:
        if merged and merged[-1][1] == lo:
            merged[-1] = (merged[-1][0], hi)
        else:
            merged.append((lo, hi))
    return merged

def danger_overlap(w, a, b):
    """Exact meas{tau in [a,b]: ||w tau|| <= 1/14}; bands [(k-1/14)/w,(k+1/14)/w]."""
    half = thr / w
    klo = floor(w * (a - half)); khi = ceil(w * (b + half))
    tot = Fr(0)
    for k in range(klo, khi + 1):
        lo = Fr(k) / w - half; hi = Fr(k) / w + half
        L = max(lo, a); R = min(hi, b)
        if R > L:
            tot += R - L
    return tot

def L_with_stranger(core_intervals, w):
    return sum(((b - a) - danger_overlap(w, a, b) for a, b in core_intervals), Fr(0))

def L_brute(S):
    """Independent breakpoint-sweep, used only to cross-validate the fast method."""
    bps = set([Fr(0), Fr(1)])
    for v in S:
        for k in range(0, v + 1):
            for s in (-thr, thr):
                x = Fr(k) + s
                if 0 <= x <= v:
                    bps.add(x / v)
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = Fr(0)
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        ok = True
        for v in S:
            x = v * mid
            f = x - (x.numerator // x.denominator)
            d = f if f <= Fr(1, 2) else 1 - f
            if d <= thr:
                ok = False
                break
        if ok:
            tot += hi - lo
    return tot

def main():
    t0 = time.time()
    WINDOW = 200   # finite enumeration window for m

    print("=" * 88)
    print("CROSS-VALIDATION: fast core-interval method == brute breakpoint sweep (exact rationals)")
    print("=" * 88)
    core6 = sorted(x for x in range(1, 14) if x != 6)
    ci6 = core_safe_intervals(core6)
    for m in (1, 4, 7, 60):
        w = 14 * m
        fast = L_with_stranger(ci6, w)
        brute = L_brute(sorted(set(core6 + [w])))
        print(f"  j=6 m={m:2d}: fast={fast}  brute={brute}  MATCH={fast == brute}")

    print()
    print("=" * 88)
    print("PER-CORE CERTIFICATE  c_j = min_m L(A_j u {14m}) > 0")
    print("Tail bound (PROVED): L(A_j u {14m}) >= limit_j - r_j/(98 m),  limit_j=(6/7)*meas(Lonely(A_j))")
    print("  => L > c_j guaranteed for all m > M_j := r_j/(98(limit_j - c_j)); verify m<=max(WINDOW,M_j).")
    print("=" * 88)
    print(f"{'j':>2} {'c_j (exact)':>16} {'c_j~':>10} {'@m':>4} {'limit~':>9} {'r':>3} {'M_thr':>7} {'verified':>12}")

    glob = None
    per_core = {}
    for j in range(1, 14):
        core = sorted(x for x in range(1, 14) if x != j)
        ci = core_safe_intervals(core)
        meas = sum((b - a for a, b in ci), Fr(0))
        limit = Fr(6, 7) * meas
        r = len(ci)
        # candidate c_j over the finite window
        best = None
        for m in range(1, WINDOW + 1):
            L = L_with_stranger(ci, 14 * m)
            if best is None or L < best[1]:
                best = (m, L)
        mj, cj = best
        gap = limit - cj
        if gap <= 0:
            Mthr = 0.0            # no dip below limit anywhere; limit itself is the floor
            Mver = WINDOW
        else:
            Mthr = float(r) / (98.0 * float(gap))
            Mver = max(WINDOW, ceil(Mthr) + 2)
        # exhaustive verification up to Mver
        cmin, cm = cj, mj
        for m in range(1, Mver + 1):
            L = L_with_stranger(ci, 14 * m)
            if L < cmin:
                cmin, cm = L, m
        assert cmin > 0, f"NONPOSITIVE floor at j={j}"
        per_core[j] = (cmin, cm, limit, r, Mthr, Mver)
        print(f"{j:>2} {str(cmin):>16} {float(cmin):>10.7f} {cm:>4} {float(limit):>9.6f} "
              f"{r:>3} {Mthr:>7.1f} {('1..'+str(Mver)):>12}")
        if glob is None or cmin < glob[1]:
            glob = (j, cmin, cm)

    print("-" * 88)
    jg, cg, mg = glob
    Sg = sorted(set([x for x in range(1, 14) if x != jg] + [14 * mg]))
    print(f"inf over ALL interior-drop cores = min_j c_j = {cg} ~= {float(cg):.10f}")
    print(f"  achieved UNIQUELY at j={jg}, m={mg}:  S = {Sg}   (14m = {14*mg} = 2*7^2)")
    print(f"  STRICTLY POSITIVE  =>  the infinite m-family is reduced to a finite, positive floor.")

    # rank the cores by their floor
    print()
    print("Per-core floors c_j ranked (deepest first):")
    for j, (cmin, cm, limit, r, Mthr, Mver) in sorted(per_core.items(), key=lambda kv: kv[1][0]):
        print(f"  j={j:2d}: c_j={str(cmin):>16} ~={float(cmin):.7f} at m={cm}")
    print(f"\nDone in {round(time.time()-t0,1)}s")
    print("\nHONEST SCOPE: this certifies inf L>0 over the interior-drop FAMILY only")
    print("(the numerical extremizers). It is NOT a proof of inf L>0 over ALL primitive")
    print("13-speed sets containing a multiple of 14. The general inf still requires")
    print("ruling out every other configuration; this angle pins the believed-worst family.")

if __name__ == "__main__":
    main()
