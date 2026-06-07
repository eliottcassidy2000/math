#!/usr/bin/env python3
"""
S649 — Working the full LRC at n=19 (prime), and why it differs from n=14 (composite).

n=14 = 2*7 was attacked by fibering over the 7-clock (S640: CRT Z/14=Z/2 x Z/7, the
doubling <2> = two 3-cycles since ord_7(2)=3). n=19 is PRIME: no nontrivial divisor,
no CRT fiber, and 2 is a PRIMITIVE ROOT mod 19 (ord=18) so the doubling is a single
18-cycle -- maximal mixing, the OPPOSITE of 14. The fiber toolkit does NOT transfer.

This script:
  (A) contrasts the divisor/doubling structure of 19 vs 14;
  (B) confirms the provable corner: the consecutive {1,..,18} is lonely EXACTLY at t=1/19
      (the tight extremal; max gap = 1/19, not more) and checks a few other configs;
  (C) maps 19's real structure: Heegner sqrt(-19) (the chi=5 step), hex(2)=19, 2n-1=37,
      Paley-19 (19=3 mod 4) but doubling != QR (19 != 7 mod 8).
No external libs.
"""
from fractions import Fraction
from math import gcd

def isprime(m): return m>1 and all(m%d for d in range(2,int(m**0.5)+1))
def ordp(a,p):
    x=a%p; k=1
    while x!=1: x=(x*a)%p; k+=1
    return k
def divisors(m): return [d for d in range(1,m+1) if m%d==0]

print("="*68)
print("(A) divisor / doubling structure: 19 (prime) vs 14 (composite)")
print("="*68)
for n in [14,19]:
    print(f"  n={n}: prime? {isprime(n)}; divisors {divisors(n)}; 2n-1={2*n-1} (prime? {isprime(2*n-1)})")
    # doubling structure mod the odd part / mod n's odd prime
    if n==14:
        print(f"        14=2*7: CRT Z/14=Z/2 x Z/7; ord_7(2)={ordp(2,7)} -> doubling = TWO 3-cycles")
        print(f"        => a FIBER bundle over the 7-clock (S640); reducible.")
    else:
        print(f"        19 prime: ord_19(2)={ordp(2,19)} = 18 = 19-1 -> 2 is a PRIMITIVE ROOT")
        print(f"        => doubling = a SINGLE 18-cycle; NO fiber, NO sub-shell. IRREDUCIBLE by")
        print(f"           the n=14 method. The hard case for divisor/CRT reductions.")

print("\n" + "="*68)
print("(B) the provable corner: {1,..,18} is lonely EXACTLY at t=1/19 (tight)")
print("="*68)
def clock(x):
    f=x-int(x); f+=1 if f<0 else 0
    return min(f,1-f)
def min_clock_at(speeds, t):
    return min(clock(v*t) for v in speeds)
n=19; delta=Fraction(1,19); speeds=list(range(1,19))
# at t=1/19 exactly (rational):
t=Fraction(1,19)
gaps=[abs(Fraction((v*1)%19,19)) for v in speeds]
mins=min(min(Fraction((v)%19,19), 1-Fraction((v)%19,19)) for v in speeds)
print(f"  {{1,..,18}} at t=1/19: min clock distance over runners = {mins} = 1/19 ? {mins==delta}")
print(f"  -> all 18 runners are >= 1/19 away; the unit runners (k=1,18) are EXACTLY 1/19.")
print(f"     So loneliness holds with EQUALITY: the bound 1/19 is achieved, not beaten (tight).")
# can any config BEAT 1/19? (search a few; the LRC conjecture says no for 19 runners)
import random
random.seed(19)
best=Fraction(0)
for _ in range(3000):
    sp=random.sample(range(1,60),18)
    # best gap over a fine rational grid of times t=j/(q), q up to 2n-1=37 and a few more
    bestcfg=Fraction(0)
    for q in [19,37,38,57]:
        for j in range(1,q):
            m=min(min(Fraction((v*j)%q,q),1-Fraction((v*j)%q,q)) for v in sp)
            if m>bestcfg: bestcfg=m
    if bestcfg>best: best=bestcfg
print(f"  random 18-speed search (rational witnesses): best achievable gap found = {float(best):.4f}")
print(f"  (LRC(19): EVERY 18-speed set should reach gap >= 1/19={float(delta):.4f}; this is OPEN")
print(f"   in general -- proven only up to 7 runners. We prove the {{1,..,18}} corner.)")

print("\n" + "="*68)
print("(C) where 19's real structure lives: CM / Heegner, not fibers")
print("="*68)
heegner=[1,2,3,7,11,19,43,67,163]
print(f"  19 is a HEEGNER number ({19 in heegner}); 19 = 4*5-1; Q(sqrt(-19)) class number 1.")
print(f"  -> the rotation field for Eisenstein norm N=5 (HYP-2277) = Q(sqrt(-19)); the")
print(f"     conjectural chi=5 chromatic step (Moser tower: sqrt(-3)->sqrt(-11)->sqrt(-19)).")
print(f"  19 = 1+6+12 = hex(2) (centered hexagonal); 2n-1=37=hex(3). Eisenstein-lattice radius-2.")
QR=sorted({(x*x)%19 for x in range(1,19)})
print(f"  19 = 3 mod 4 -> Paley-19 tournament EXISTS (self-converse, S638); QR mod 19 = {QR}")
print(f"  but 2 is a NON-residue (2 in QR? {2 in QR}) and <2>=all of (Z/19)* (primitive root),")
print(f"  so <2> != QR: 19 is NOT in the p=7 mod 8 'doubling=Paley' family (S640). 19 = {19%8} mod 8.")
print("  => n=19's leverage is the CM/Heegner/cube-root side (sqrt(-19), the chi=5 frontier),")
print("     NOT the divisor-fiber side. The full LRC(19) wants the cyclotomic-depth q* attack")
print("     (opus S704) at the prime modulus, with witnesses t=j/19 and t=j/37.")
