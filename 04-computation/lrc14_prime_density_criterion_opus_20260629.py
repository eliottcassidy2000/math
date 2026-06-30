"""
LRC(14): the PRIME-DENSITY criterion -- live successor to the refuted Vmax-rho* route.
=====================================================================================
opus session 2026-06-29 (follow-up to lrc14_rho_star_floor_refutation: inf rho*=0).

CRITERION.  For a prime p, the exact lonely count is
    N(S,p) = #{ a in Z/p : 14*min(r, p-r) >= p for all s in S,  r = (s a) mod p }.
SOUND:  N(S,p) >= 1  ==>  exists witness t=a/p with ||s t|| >= 1/14 for all s  ==>  M(S) >= 1/14.

WHY THIS FIXES THE rho* ROUTE:
  * No blind spot: it certifies every covering set on which rho* (Vmax ruler) vanished.
  * Sidesteps THM-566 (witness denominator unbounded): we do NOT bound D. There are
    infinitely many primes p coprime to the (finitely many) speeds; at such p the
    divisor-loaded structure is invisible (a loaded runner is a generic unit), and
    N(S,p)/p  ->  lonely-measure(S) := meas{x in [0,1): ||s x||>=1/14 for all s}.
  * EXACT, not a weakening:  lonely-measure(S) > 0  <=>  M(S) > 1/14  (the target).
    So the criterion recasts the goal as an ARITHMETIC density over Z/p:
        N(S,p) = sum_{a} prod_s 1_safe(sa)
               = p - sum_T (-1)^|T| #{a: s a in band  for all s in T}   (incl-excl)
    = main term  p * lonely-measure(S)  +  exponential-sum error (Gauss/Kloosterman).
    The open core is: lower-bound this uniformly over covering 13-sets
    (== bound the character-sum error below the main term). This is precisely where
    mac-mini's BV-Fourier uniform bound (HYP-3529 / THM-578 R-tail) applies.

WHY NO CHEAP FLOOR (the union-bound obstruction): each unsafe set has measure 1/7,
  and 13 * (1/7) = 13/7 > 1, so a first-moment/union bound is always vacuous.
  Positivity genuinely requires the second-order cancellation -- the analytic core.

Verifications below are exact (rationals) / exact-integer (counts mod p).
"""
from math import gcd, lcm
from fractions import Fraction as F
from functools import reduce
import random

def N_mod(S, D):
    c = 0
    for a in range(D):
        ok = True
        for s in S:
            r = (s * a) % D
            if 14 * min(r, D - r) < D: ok = False; break
        if ok: c += 1
    return c

def lonely_measure(S):           # = lim_{p->inf, gcd(p,prod S)=1} N(S,p)/p
    bps = {F(0), F(1)}
    for s in S:
        for a in range(s + 1):
            bps.add(F(14 * a + 1, 14 * s)); bps.add(F(14 * a - 1, 14 * s))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi <= lo: continue
        mid = (lo + hi) / 2; ok = True
        for s in S:
            f = mid * s; f = f - (f.numerator // f.denominator)
            if min(f, 1 - f) < F(1, 14): ok = False; break
        if ok: tot += hi - lo
    return tot

def is_cov(S): return all(any(s % q == 0 for s in S) for q in range(2, 15))
def prim(S):   return reduce(gcd, S) == 1

core = [1,2,3,4,5,6,7,8,9,10,11,13]
rho0 = [[1,2,3,12,18,20,21,22,23,24,25,26,28],
        [1,2,3,6,18,20,21,22,23,24,26,27,28],
        [1,2,3,18,20,21,22,23,24,25,26,27,28]]

print("[A] The refuted route's blind-spot sets ARE certified by a prime witness:")
for S in rho0:
    # smallest prime witness
    p = 17
    while True:
        if all(gcd(s,p)==1 for s in S) and N_mod([s%p for s in S], p) >= 1:
            break
        p += 1
        while not all(p % d for d in range(2,int(p**.5)+1)): p += 1
    print(f"    S={S[:4]}...  rho*=0 but N(S,{p})={N_mod([s%p for s in S],p)}>=1  -> M>=1/14")

print("\n[B] THM-566 divisor-loaded S_B={1..11,13,84*lcm(1..B)}: density STABLE at a fixed large prime")
for p in (1009, 4999):
    row = []
    for B in [11, 26, 41, 60, 100]:
        big = 84 * lcm(*range(1, B + 1))
        if big % p == 0: row.append(f"B{B}:skip"); continue
        n = N_mod([s % p for s in core] + [big % p], p)
        row.append(f"B{B}({len(str(big))}d):N/p={n/p:.4f}")
    print(f"    p={p}: " + "  ".join(row))

print("\n[C] lonely-measure(S) = large-prime density limit, > 0  <=>  M(S) > 1/14:")
print(f"    core {{1..11,13}}        : {float(lonely_measure(core)):.5f}")
for S in rho0:
    print(f"    rho*=0 {S[:4]}... : {float(lonely_measure(S)):.5f}")
random.seed(3); samp = []; t = 0
while len(samp) < 300 and t < 2_000_000:
    t += 1; S = sorted(random.sample(range(1, 80), 13))
    if prim(S) and is_cov(S): samp.append(S)
ms = sorted(float(lonely_measure(S)) for S in samp)
print(f"    random covering sample (n={len(samp)}): min lonely-measure = {ms[0]:.5f}; "
      f"5 smallest = {[round(x,4) for x in ms[:5]]}")
print("\n    => the live target is: lonely-measure(S) >= c0 > 0 uniformly over covering 13-sets,")
print("       i.e. N(S,p) > 0 at some prime p -- an exponential-sum lower bound (the analytic core).")
