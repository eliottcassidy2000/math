"""
opus-2026-07-11-S250: LRC(14) as a CHANGE-OF-BASE problem (owner's reframing). Translate + verify.

DICTIONARY (owner's language -> project's machinery):
  runner position (digits after the point, base b)  =  {v_i t} on the circle (archimedean place)
  "base"                                            =  modulus q  (the clean-ruler / clearing modulus)
  "digit of runner i in base q"                     =  residue v_i * a mod q
  "a base where runner i's change is the unit"      =  a q where runner i is NON-resonant (spread, not stuck at 0)
  "the problem is a change of bases"                =  hunt for a modulus q (+ multiplier a) that cleans ALL runners
  best_q(v) = max_a min_i ||v_i a / q||             =  the floor achievable using base q
  M(v) = max_q best_q(v)                            =  the LRC value; LRC(14) <=> M(v) >= 1/14 for all v

WHAT THE REFRAMING CORRECTLY CAPTURES:
  * Per-runner, a clean base always exists (trivial). The DIFFICULTY is you need ONE base clean for ALL runners
    at once -- the CRT-combination of the per-runner clean bases is exactly the covering-system structure.
  * Verified base-landscapes below re-derive WHY the AP and V* are the extremizers:

FINDINGS (verified):
  (1) AP {1..13} and V* {1..11,13,24} are BASE-RIGID: max_q best_q = 1/14, peak at base q*=14. They are
      RESONANT (best_q = 0, some runner stuck at origin) in EVERY base 2..13 -- because each contains a
      multiple of every d in {2..13} (divisor-complete). Their ONLY clean base is 14, where they hit EXACTLY
      the floor 1/14 (zero margin). That is precisely why they are the hard extremizers: they defeat every
      base below 14 and only marginally clear at 14. AP and V* have IDENTICAL coarse landscapes (same top
      bases 14,85,71,57; same resonant set) -- base-TWINS; the mod-14 fine structure (S249: full vs
      vacate-8/double-2) is the only distinguisher, invisible at coarse resolution. V* is a base-twin of AP
      because base 14 = 2*7 is composite (the doubling 12->24).
  (2) A CLEARABLE (non-tight) family has a SMALL clean base beating 1/14 (e.g. best 5/38 at q=38) -- one base
      cleans all runners with margin. Easy.
  (3) A DIVISOR-COMPLETE adversarial family is RESONANT/marginal in every small base -- its clean base is
      pushed LARGE (no small base works; MISTAKE-110's unbounded-modulus phenomenon).

THE DEEP OBSTRUCTION (why "just a change of base" does not close it):
  A base change that makes a runner's MOTION the unit (unimodular relabel of the speed lattice Z^k) simultaneously
  DISTORTS the danger zone, which is defined per-runner in the ORIGINAL coordinates (the archimedean band
  ||.|| < 1/14). You can clean the dynamics (finite places / choice of q) OR keep the danger metric simple
  (archimedean place), not both. LRC couples the archimedean place (the band) with the finite places (the
  bases/moduli): the conjecture asserts a FINITE certificate (some base q) for an ARCHIMEDEAN fact (a lonely
  point exists). This is why pure finite-place (p-adic) reasoning cannot alone settle it, and why the natural
  next decomposition is prime-local: base = prime power, with 14 = 2*7 splitting into the 2-adic and 7-adic
  pieces (the clean ruler is strongest at primes, THM-712; apex prime 7).
"""
from math import gcd
from functools import reduce
from fractions import Fraction
import random
def best_at_base(v,q):
    bq=Fraction(0)
    for a in range(1,q):
        if gcd(a,q)!=1: continue
        m=min(min((vi*a)%q, q-(vi*a)%q) for vi in v)
        f=Fraction(m,q)
        if f>bq: bq=f
    return bq
def landscape(v,qmax=90): return {q:best_at_base(v,q) for q in range(2,qmax+1)}
def peak(v,qmax=90):
    L=landscape(v,qmax); qb=max(L,key=lambda q:L[q]); return qb,L[qb],L
def main():
    tick=Fraction(1,14)
    ap=list(range(1,14)); vstar=[1,2,3,4,5,6,7,8,9,10,11,13,24]
    print("=== BASE-RIGID EXTREMIZERS (no base beats the 1/14 floor) ===")
    for name,v in [("AP {1..13}",ap),("V* {..,13,24}",vstar)]:
        qb,fb,L=peak(v); top=sorted(L.items(),key=lambda kv:-kv[1])[:4]
        print(f"{name}: peak q*={qb}, best={fb}; top bases {[(q,str(f)) for q,f in top]}; dirty bases 2..90: {[q for q in L if L[q]==0][:12]}")
    print("\n=== CLEARABLE: some small base beats 1/14 ===")
    random.seed(2)
    for _ in range(400):
        v=sorted(random.sample(range(1,30),13))
        if reduce(gcd,v)!=1: continue
        qb,fb,L=peak(v)
        if fb>tick and qb<=40:
            print(f"{v}: clean base q*={qb}, best={fb}={float(fb):.4f} (>1/14)"); break
    print("\n=== DIVISOR-COMPLETE adversarial: small bases dirty, clean base large ===")
    def dc(v): return all(any(x%d==0 for x in v) for d in range(2,15))
    random.seed(5)
    for _ in range(6000):
        v=sorted(random.sample(range(1,220),13))
        if reduce(gcd,v)!=1 or not dc(v): continue
        L=landscape(v,60); below=[q for q in range(15,45) if L[q]<=tick]; qb,fb,_=peak(v,60)
        print(f"{v}: best 2..60 = {fb} at q*={qb}; bases 15..44 with no clean margin: {len(below)}/30"); break
if __name__=='__main__':
    main()
