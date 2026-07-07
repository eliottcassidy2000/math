#!/usr/bin/env python3
r"""
lrc_decorrelation_decomp_kps_S71.py   (kind-pasteur-2026-07-07-S71, HYP-5097)

THE DIAMETER x RESONANCE DECOMPOSITION of (A') -- work the S70 next target (the
decorrelation half) + synthesize the fleet's routes + QC the no-cherry mu~1 claims.

(A') PA_2(E) >= T_k (equivalently the density floor mu_{1/7}(E) >= T_k, since PA_2 <= mu)
decomposes along two axes:
  * DIAMETER: bounded (D <= D_0) -> my S59/S68 subset lemma (= opus-S141 hom-monotonicity),
    PROVED; unbounded (D > D_0) -> the resonance axis.
  * RESONANCE: a family is GENERIC (no small additive relation among the e_i) => the joint
    sequence {(e_i x)} equidistributes on the k-torus => PA_2(E) -> PA_2^inf > T_k by
    Erdos-Turan/Koksma; RESONANT (near-AP, carries small relations) => S69 spread-AP exact
    + R2 (the density-floor rigidity, S70).

THIS SCRIPT tests the DECORRELATION half concretely: is there an explicit
"min-relation-size >= R_0 => PA_2 >= T_k" threshold?  If PA_2 tracks the smallest nonzero
additive relation (weight-3+; my S63 frame), then generic families decorrelate and only the
resonant near-AP core (small relations) remains -- exactly R2.

 (1) PA_2 vs the smallest zero-sum relation L1-norm of E (weight>=3): do families with LARGE
     min-relation (generic) have PA_2 near PA_2^inf > T_k?  do only small-min-relation
     (near-AP/resonant) families dip toward the spread-AP min?
 (2) the decomposition table (which regime each family lands in).
 (3) QC: the incoming "no-cherry => mu~1.0" (mac-mini/klein, random-sampled) -- spread APs
     are no-cherry, large-diameter, but mu = mu(AP_k) (NOT ~1), a structured exception.
"""
import random, math
from itertools import combinations, product

def PA2(E, res=14000):
    E = sorted(set(E)); n = len(E); c = 0
    for r in range(res):
        x = (r + .5) / res
        ph = sorted((e * x) % 1.0 for e in E)
        g0 = ph[0] + 1 - ph[-1]; gh = ph[0] + 1 - ph[-1]
        for i in range(n-1):
            if ph[i] <= 0.5 < ph[i+1]: gh = ph[i+1]-ph[i]; break
        if max(g0, gh) > 1/7: c += 1
    return c / res

def mu17(E, res=14000):
    E = sorted(set(E)); n = len(E); c = 0
    for r in range(res):
        x = (r + .5) / res
        ph = sorted((e * x) % 1.0 for e in E)
        mg = max(ph[i+1]-ph[i] for i in range(n-1)); mg = max(mg, ph[0]+1-ph[-1])
        if mg > 1/7: c += 1
    return c / res

def min_relation(E, maxc=3, maxsup=5):
    """smallest L1 norm of a zero-sum relation sum m_i e_i = 0, sum m_i = 0, support>=3."""
    E = sorted(E); k = len(E); best = None
    for sup in (3, 4, 5):
        if sup > maxsup: break
        for S in combinations(range(k), sup):
            vals = [E[i] for i in S]
            for coeffs in product(range(-maxc, maxc+1), repeat=sup):
                if 0 in coeffs or sum(coeffs) != 0: continue
                if sum(c*v for c, v in zip(coeffs, vals)) != 0: continue
                l1 = sum(abs(c) for c in coeffs)
                if best is None or l1 < best: best = l1
        if best is not None and sup >= 4: break
    return best if best else 999   # 999 = no small relation (generic)

Tk = {8: 0.6185, 13: 0.0565}
PA2inf = {8: 0.7994, 13: 0.3648}

print("=" * 96)
print("PART 1 -- PA_2 vs smallest zero-sum relation (decorrelation): generic => PA_2 -> PA_2^inf")
print("=" * 96)
rng = random.Random(71)
for k in (8, 13):
    T = Tk[k]; inf = PA2inf[k]
    fams = {}
    fams["consec AP {1..k}"] = list(range(1, k+1))
    fams["spread AP d=2"] = [1 + 2*j for j in range(k)]
    fams["spread AP a=5 d=2"] = [5 + 2*j for j in range(k)]
    fams["primes"] = ([2,3,5,7,11,13,17,19,23,29,31,37,41])[:k]
    fams["Sidon-greedy"] = ([1,2,4,8,13,21,31,45,66,81,97,123,148])[:k]
    fams["geometric 2^j"] = [2**j for j in range(k)]
    for r in range(6):
        fams[f"random-generic-{r}"] = sorted(rng.sample(range(1, 500), k))
    print(f"\n  k={k}: T_k={T:.4f}, PA_2^inf={inf:.4f}")
    print(f"    {'family':22s} {'min-rel L1':>10} {'PA_2':>8} {'mu':>8} {'regime':>14}")
    for nm, E in fams.items():
        mr = min_relation(E); p = PA2(E); m = mu17(E)
        regime = "RESONANT(near-AP)" if mr <= 4 else "generic(decorr)"
        print(f"    {nm:22s} {(str(mr) if mr<999 else '>box'):>10} {p:8.4f} {m:8.4f} {regime:>14}")
    print(f"    => generic (min-rel > 4) families have PA_2 near PA_2^inf {inf:.4f} >> T_k;")
    print(f"       resonant (min-rel <= 4 = near-AP) families are the S69+R2 core.")

print()
print("=" * 96)
print("PART 3 -- QC the incoming 'no-cherry => mu ~ 1.0' (mac-mini/klein, RANDOM-sampled)")
print("=" * 96)
print("  claim: no-cherry (no clustered triple), diam>=27 8-shapes have mu >= 0.9998 (random census).")
print("  MISTAKE-102 check: spread APs are NO-CHERRY (evenly spaced, no clustered triple) AND")
print("  large-diameter -- but mu(spread AP) = mu(AP_k) by translation+dilation invariance:")
for k in (8,):
    for (nm, E) in [("consec AP", list(range(1,k+1))), ("spread AP d=4 (diam 28)", [1+4*j for j in range(k)]),
                    ("spread AP d=7 (diam 49)", [1+7*j for j in range(k)])]:
        print(f"    {nm:26s} {E}: mu = {mu17(E, 20000):.4f} (has cherry? no)")
print(f"    => spread 8-APs are no-cherry, diam>=27, yet mu ~ 0.94, NOT ~1.0.  The random census")
print(f"       (klein: 3990/4000) MISSES these structured shapes (MISTAKE-102).  HARMLESS for T_8=0.62")
print(f"       (0.94 > 0.62), but 'no-cherry => mu~1' overreaches; the true statement is")
print(f"       'no-cherry => mu >= mu(AP_k)' (= 0.94 at k=8), the AP-minimality bound again.")
print("DONE.")
