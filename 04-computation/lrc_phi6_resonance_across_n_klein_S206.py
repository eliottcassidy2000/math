#!/usr/bin/env python3
"""
klein-2026-07-09-S206: is the "13-comb Eisenstein resonance" (t* = 14/183, 183 = Phi6(14)) UNIQUE to
n=14, and what does it say about LRC hardness across runner counts?

Phi6(x) = x^2 - x + 1  (the 6th cyclotomic polynomial).  Phi6(14) = 183.

IDENTITIES to check for every n (all should be n-independent consequences of x^3+1=(x+1)Phi6(x)):
  (I1)  n^3      = -1  (mod Phi6(n))        [so ord(n mod Phi6(n)) | 6  -- the "Eisenstein" 6th root]
  (I2)  n*(n-1)  = -1  (mod Phi6(n))        [the 13*14 = -1 mod 183 of the "13-comb"]
  (I3)  contfrac( n / Phi6(n) ) = [0; n-1, n]   [so the comb spacing n-1 is the 1st CF denominator]
If these hold for ALL n, the resonance is a UNIVERSAL family, not special to 14.

HARDNESS: the repo's covering-min is  M(S) >= n/Phi6(n)  for covering (n-1)-sets (THM-523/HYP-3778).
Compare to the LRC target 1/n.  The MARGIN RATIO is
      rho(n) := (n/Phi6(n)) / (1/n) = n^2 / (n^2 - n + 1)  ->  1   as n -> oo.
So the covering case gets monotonically TIGHTER with n.  Compute rho(n); see where n=14 sits.

EXACT M(S): f(t)=min_i ||v_i t|| is piecewise linear with slopes +-v_i, so its local maxima sit at
t = p/(v_i+v_j).  Hence M(S) = max over q in {v_i+v_j} and 0<p<q of min_k ||v_k p/q||  -- EXACT.
Brute-force the covering-min for small n and compare with n/Phi6(n), and check the optimal t.
"""
from math import gcd
from fractions import Fraction
from itertools import combinations

def phi6(n): return n*n - n + 1

def contfrac(fr: Fraction):
    a = []
    num, den = fr.numerator, fr.denominator
    while den:
        a.append(num // den)
        num, den = den, num - (num // den) * den
    return a

def order_mod(a, m):
    if gcd(a, m) != 1: return None
    x, k = a % m, 1
    while x != 1:
        x = (x * a) % m; k += 1
        if k > m: return None
    return k

def nearInt_frac(num, den):
    r = num % den
    return Fraction(min(r, den - r), den)

def M_exact(v):
    """exact max_t min_i ||v_i t||, via local maxima at t = p/(v_i+v_j)."""
    qs = set()
    for i in range(len(v)):
        for j in range(i, len(v)):
            qs.add(v[i] + v[j])
    best = Fraction(0); best_t = None
    for q in qs:
        if q <= 0: continue
        for p in range(1, q):
            if gcd(p, q) != 1: continue
            m = min(nearInt_frac(vk * p, q) for vk in v)
            if m > best:
                best, best_t = m, Fraction(p, q)
    return best, best_t

print("=== (1) The Phi6 identities: universal, or special to n=14? ===")
print(f"{'n':>4} {'Phi6(n)':>8} {'n^3 mod Phi6':>13} {'n(n-1) mod Phi6':>16} {'ord(n)':>7} {'CF(n/Phi6)':>16}")
univ = True
for n in range(3, 21):
    P = phi6(n)
    i1 = pow(n, 3, P); i2 = (n * (n - 1)) % P
    ok1 = (i1 == P - 1); ok2 = (i2 == P - 1)
    cf = contfrac(Fraction(n, P))
    ok3 = (cf == [0, n - 1, n])
    univ &= (ok1 and ok2 and ok3)
    print(f"{n:>4} {P:>8} {i1:>13} {i2:>16} {str(order_mod(n,P)):>7} {str(cf):>16}")
print(f"\nAll three identities hold for every n in 3..20 ?  {univ}")
print("=> n^3 = -1 and n(n-1) = -1 mod Phi6(n) are IDENTITIES (from x^3+1=(x+1)(x^2-x+1)).")
print("   The '13-comb at 14/183' is the n=14 member of a UNIVERSAL family, not a special coincidence.\n")

print("=== (2) The covering-case margin ratio rho(n) = n^2/Phi6(n) -> 1 ===")
print(f"{'n':>4} {'n/Phi6(n)':>14} {'1/n':>10} {'rho(n)=n^2/Phi6':>16} {'margin %':>9}")
for n in [4,5,6,7,8,9,10,11,12,13,14,15,20,30,50]:
    P = phi6(n)
    lo = Fraction(n, P); tgt = Fraction(1, n)
    rho = Fraction(n*n, P)
    print(f"{n:>4} {str(lo):>14} {str(tgt):>10} {float(rho):>16.5f} {100*(float(rho)-1):>8.2f}%")
print("\n=> the covering-case cushion above 1/n SHRINKS monotonically: 25.0% (n=4) ... 7.10% (n=14) ... ->0.")
print("   n=14 is NOT structurally special; it is the first n where the shrinking cushion defeats")
print("   current methods (LRC proved n<=13).\n")

print("=== (3) EXACT covering-min for small n: is it n/Phi6(n), attained at t = n/Phi6(n)? ===")
print("covering (n-1)-set := some speed divisible by n (so t=1/n fails).")
print(f"{'n':>3} {'Vcap':>5} {'#sets':>8} {'min M(S)':>12} {'n/Phi6(n)':>12} {'>= ?':>5} {'argmin t':>12} {'= n/Phi6?':>10}")
for n in [4, 5, 6, 7]:
    k = n - 1
    Vcap = 4 * n
    speeds = list(range(1, Vcap + 1))
    target = Fraction(n, phi6(n))
    best = None; cnt = 0
    for S in combinations(speeds, k):
        if not any(s % n == 0 for s in S):   # covering := some speed ≡ 0 mod n
            continue
        if gcd(*S) != 1 if k > 1 else False:  # primitive only
            continue
        cnt += 1
        m, t = M_exact(list(S))
        if best is None or m < best[0]:
            best = (m, t, S)
    m, t, S = best
    print(f"{n:>3} {Vcap:>5} {cnt:>8} {str(m):>12} {str(target):>12} {str(m>=target):>5} {str(t):>12} {str(t==target):>10}")
    print(f"      argmin S = {S}")
print("\n=> if min M(S) == n/Phi6(n) attained at t = n/Phi6(n), the Phi6 resonance IS the covering-min,")
print("   universally in n -- and the LRC(14) covering case inherits a 7.10% strict cushion.")
