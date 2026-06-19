#!/usr/bin/env python3
"""
LRC(14) S3 residual -- DIRECT COVERING-SYSTEM BOUND, refined to the EXCESS-GAP functional.
kind-pasteur-2026-06-18-S5-wf.

NEW RIGOROUS LOWER BOUND (proved here):
  For x in [0,1), let the points be {frac(e x): e in E} (k points, since 0 in E one is 0).
  Let g_1(x),...,g_r(x) be the circular GAPS (r = #distinct points). Define the
  EXCESS-GAP functional
        Psi(E) := integral_0^1  sum_j ( g_j(x) - 2/7 )_+  dx.
  CLAIM:  Psi(E) <= mu(E).
  PROOF: pointwise, sum_j (g_j - 2/7)_+ is 0 unless some gap > 2/7 (i.e. unless g(x)=1),
  and when nonzero it is <= sum_j g_j = 1. So the integrand is <= 1[maxgap>2/7] = the
  integrand of mu. Integrate. QED (rigorous, no approximation).

  AVERAGING INTERPRETATION (the 'covering-system' content): with
     Phi_theta(x) = 1[ every frac(e x) avoids the closed 2/7-arc centered at theta ],
  one has  integral_theta Phi_theta(x) dtheta = sum_j (g_j(x)-2/7)_+  (an arc of length
  2/7 fits in gap g_j iff its center lies in a sub-arc of length (g_j-2/7)_+). Hence
     Psi(E) = integral_x integral_theta Phi_theta(x) dtheta dx
            = average over arc-position theta of the 0-anchored-type strip-avoidance meas.
  So Psi is the ARC-AVERAGED covering measure -- it restores the bulk that a SINGLE fixed
  arc (mu_0) loses. THIS is the correct covering-system functional.

GOAL: is inf_E Psi(E) > 0 uniformly (k<=13)?  If yes -> mu(E) >= Psi(E) >= c0 -> B(k) closed.

EXACT COMPUTATION of Psi(E):
  On each order-cell (between collision breakpoints), the sorted points are linear in x,
  so each gap g_j(x) is linear (a_j + b_j x). Then (g_j - 2/7)_+ is piecewise linear with
  one extra breakpoint where g_j = 2/7 (rational). We integrate exactly by splitting at
  all collision breakpoints AND all gap=2/7 crossings, then on each atomic interval the
  integrand sum_j (g_j-2/7)_+ is linear -> integrate via trapezoid of exact endpoints
  (exact for linear). RIGOROUS.
"""
from fractions import Fraction as F
from math import gcd

G0 = F(2, 7)
def _frac(q): return q - q.__floor__()

def _collision_bps(E):
    E = sorted(set(E)); k = len(E); bp = set([F(0), F(1)])
    for i in range(k):
        for j in range(i+1, k):
            d = E[j]-E[i]
            if d == 0: continue
            for m in range(0, d+1): bp.add(F(m, d))
    return bp

def _gap_eq_bps(E):
    # gap between two consecutive points equals 2/7: difference of two frac() = 2/7 mod 1.
    # frac(e_b x) - frac(e_a x) = (e_b-e_a) x mod 1. Set = 2/7 or =5/7 (=-2/7). Also the WRAP
    # gap = 1 - sum, but a single gap=2/7 between SOME ordered pair is captured by pair diffs.
    E = sorted(set(E)); bp = set()
    diffs = set()
    for i in range(len(E)):
        for j in range(len(E)):
            if i != j: diffs.add(abs(E[j]-E[i]))
    for D in diffs:
        if D == 0: continue
        n = 0
        while True:
            cand = [(F(n)+F(2,7))/D, (F(n)+F(5,7))/D]
            anyin = False
            for x in cand:
                if F(0) <= x < F(1): bp.add(x); anyin = True
            if min(cand) >= F(1): break
            n += 1
            if n > D+2: break
    return bp

def gaps_at(E, x):
    pts = sorted(set(_frac(e*x) for e in E))
    if len(pts) == 1: return [F(1)]
    gaps = [pts[t+1]-pts[t] for t in range(len(pts)-1)]
    gaps.append(pts[0]+1-pts[-1])
    return gaps

def excess_at(E, x):
    return sum((g - G0) for g in gaps_at(E, x) if g > G0)

def psi_exact(E):
    bp = sorted((_collision_bps(E) | _gap_eq_bps(E)))
    bp = [b for b in bp if F(0) <= b <= F(1)]
    if bp[0] != F(0): bp = [F(0)] + bp
    if bp[-1] != F(1): bp = bp + [F(1)]
    tot = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        # integrand sum_j (g_j-2/7)_+ is LINEAR on the open cell (a,b) (constant point-order,
        # constant active-gap set since we split at every collision AND every gap=2/7 crossing).
        # Exact integral of a linear function = value at the midpoint times length. RIGOROUS, EXACT.
        mid = (a + b) / 2
        tot += excess_at(E, mid) * (b - a)
    return tot

def mu_exact(E):
    bp = sorted((_collision_bps(E) | _gap_eq_bps(E)))
    bp = [b for b in bp if F(0) <= b < F(1)]
    tot = F(0)
    pts_bp = bp + [F(1)]
    for a, b in zip(bp, pts_bp[1:]):
        if b <= a: continue
        mid = (a+b)/2
        gs = gaps_at(E, mid)
        if max(gs) > G0: tot += (b-a)
    return tot

def normalize(E):
    E = sorted(set(E)); g = 0
    for e in E: g = gcd(g, e)
    return [e//g for e in E] if g else E

def header(t):
    print("\n"+"="*74); print(t); print("="*74)

if __name__ == "__main__":
    import random
    header("(A) VALIDATE Psi(E) <= mu(E) and exactness of psi_exact")
    tests = [[0,1,2],[0,1,2,3],[0,1,2,3,4,5,6],list(range(13)),
             [0,2,3,5,7,8],[0,1,4,9,16,25],[0,7,14,21,28,35],[0,1,500,501]]
    for E in tests:
        E = normalize(E)
        ps = psi_exact(E); m = mu_exact(E)
        ok = "OK" if ps <= m else "VIOLATION!!"
        print(f"  E={E if len(E)<=8 else ('0..'+str(E[-1])+f' (k={len(E)})')}: "
              f"Psi={float(ps):.5f}={ps}  mu={float(m):.5f}={m}  Psi<=mu {ok}")

    header("(B) Psi for CONSECUTIVE E={0..k-1}, k=3..13 (the canonical near-extremal)")
    for k in range(3, 14):
        E = list(range(k))
        ps = psi_exact(E); m = mu_exact(E)
        print(f"  k={k:2d}: Psi={float(ps):.6f}={ps}   mu={float(m):.6f}   ratio Psi/mu={float(ps/m):.4f}")

    header("(C) INFIMUM SEARCH: worst-case Psi over bounded-spread shapes (the floor)")
    # consecutive, perforated near-APs (the known mu-minimizers), random bounded
    random.seed(202606185)
    worst = None
    # structured: all k-subsets of {0..k+5} containing 0, for small k; random for large k
    from itertools import combinations
    for k in range(3, 11):
        cap = k + 6
        for sub in combinations(range(1, cap+1), k-1):
            E = normalize([0]+list(sub))
            if len(E) < k: continue
            ps = psi_exact(E)
            if worst is None or ps < worst[1]:
                worst = (E, ps, mu_exact(E))
    print(f"  structured (k<=10, subsets of 0..k+6): worst Psi={float(worst[1]):.6f}={worst[1]}")
    print(f"     at E={worst[0]}  (mu there={float(worst[2]):.6f})")

    header("(D) Psi floor at k=13 (the binding case): structured + random bounded")
    worst13 = None
    # consecutive and single/double perforations of {0..12+pad}
    base_caps = [13, 14, 15, 16, 18]
    cand = []
    for cap in base_caps:
        # choose 13 from {0..cap} containing 0 -- enumerate perforations (remove cap+1-13 elts)
        from itertools import combinations as C
        removable = list(range(1, cap+1))
        nrem = (cap+1) - 13
        if nrem < 0: continue
        if nrem == 0:
            cand.append(normalize(list(range(cap+1))))
        else:
            for rem in C(removable, nrem):
                E = normalize([e for e in range(cap+1) if e not in set(rem)])
                if len(E) == 13: cand.append(E)
    seen = set()
    for E in cand:
        key = tuple(E)
        if key in seen: continue
        seen.add(key)
        ps = psi_exact(E)
        if worst13 is None or ps < worst13[1]:
            worst13 = (E, ps)
    print(f"  k=13 structured perforations of 0..cap (caps {base_caps}): {len(seen)} shapes")
    print(f"     worst Psi={float(worst13[1]):.6f}={worst13[1]}  at E={worst13[0]}")
    print(f"     mu there = {float(mu_exact(worst13[0])):.6f}")

    header("(E) does Psi -> 0 with SPREAD at fixed k? (decides uniformity of THIS bound)")
    print("  k=4, E={0,1,N,N+1} (one short relation), push N:")
    for N in (2,5,20,100,500,2000):
        E = normalize([0,1,N,N+1])
        print(f"    N={N:5d}: Psi={float(psi_exact(E)):.6f}  mu={float(mu_exact(E)):.6f}")
    print("  k=4, relation-free large spread {0,1,N,N^2}:")
    for N in (3,7,20,50):
        E = normalize([0,1,N,N*N])
        print(f"    N={N:4d}: Psi={float(psi_exact(E)):.6f}  mu={float(mu_exact(E)):.6f}")
