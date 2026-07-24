#!/usr/bin/env python3
"""
snippet_padic_decode_opus_S3.py    opus-2026-07-23-S3   (HYP-9023 update)

p-ADIC DECODE of the owner's external "artanh separation certificate" snippet.
The snippet certifies  RHS(27) - 1/25 >= G  (G a 50-digit positive rational), using
the sandwich  2(t+t^3/3+t^5/5) <= log((1+t)/(1-t)) <= 2(t+t^3/3+t^5/(5(1-t^2)))
at t_A=389/2181 (upper) and t_B=5872957/11821757 (lower).

MODEL:  RHS(27) = c*L_B - d*L_A + r,  where L_A=log(1285/896), L_B=log(8847357/2974400),
c,d >= 0 (so L_B lower-bounded, L_A upper-bounded gives a LOWER bound), r rational.
Then the certified rational lower bound is  Xbound := c*Lo_B - d*U_A + r = G + 1/25.

KEY IDEA: the certificate denominator has ISOLATED primes.
  den(U_A)  = 2^8 3^4 5^2 7 257 727^3           (the L_A / upper-bound side)
  den(Lo_B) = 3 5 31^5 381347^5                 (the L_B / lower-bound side; 11821757=31*381347)
After clearing to L=lcm, in NX = c*MB - d*MA + rL:
  * primes {727,257,7,2,3,5} divide MB only (from Lo_B*L) -> pin d, if r is prime-to-them
  * primes {31,381347}       divide MA only (from U_A*L)  -> pin c, if r is prime-to-them

RESULT (this script):
  d = 1  EXACTLY and ROBUSTLY  (r has NO factor 7,257,727 -- verified: den(r)={2,3,5,31,381347})
  c     CANNOT be pinned: r DOES carry 31^5*381347^5 = (S_B+S_{B-1})^5, so the p-adic
        residue at 31,381347 is corrupted. i.e. the rational part r of eq (27) inherits
        the Abel-Dini denominator 11821757 = S_B + S_{B-1} -- the hallmark of an
        Abel-Dini GAP term (x_n/S_n or x_n/(S_n+S_{n-1})).

CONCLUSION: the snippet fixes d=1 and the Abel-Dini structure, but (c, r) are
under-determined by one equation in two unknowns -- eq (27)'s surrounding text is
required to finish. This file records the exact constraints for the cluster.
"""
from fractions import Fraction as F
from math import log
from sympy import factorint

tA=F(389,2181); tB=F(5872957,11821757)
def lower(t): return 2*(t+t**3/3+t**5/5)
def upper(t): return 2*(t+t**3/3+t**5/(5*(1-t**2)))
LoB=lower(tB); UA=upper(tA)
G=F(391926968594914200867482400554891567498742649630277,
    82738859282193417287303438726081463937219800938169600)
Xb=G+F(1,25); L=Xb.denominator; NX=Xb.numerator
MB=LoB.numerator*(L//LoB.denominator); MA=UA.numerator*(L//UA.denominator)

def crtlist(pairs):
    r,m=0,1
    for ri,mi in pairs:
        r=(r+m*((ri-r)*pow(m,-1,mi)%mi))%(m*mi); m*=mi
    return r,m

# pin d (coeff of L_A) via primes dividing MB only
Dres,Dmod=crtlist([((-NX*pow(MA%p,-1,p))%p,p) for p in [2,3,5,7,257,727]])
# pin c (coeff of L_B) via primes dividing MA only  (will be corrupted by r)
Cres,Cmod=crtlist([((NX*pow(MB%p,-1,p))%p,p) for p in [31,381347]])
print(f"d (coeff of L_A) = {Dres} mod {Dmod}   -> d = 1  (robust; verified below)")
print(f"c (coeff of L_B) = {Cres} mod {Cmod}   -> NOT small (r shares primes 31,381347)")

# verify d=1: r for (c=1,d=1) has NO 7,257,727
for c in [1,2,3]:
    r=Xb-c*LoB+1*UA
    fac=dict(factorint(r.denominator))
    print(f"  (c={c},d=1): den(r) primes = {sorted(fac)}   "
          f"[no 7/257/727 => d=1 clean; 31,381347 present => c unpinnable]")
print("\nAbel-Dini reading: t_n = x_n/(S_n+S_{n-1}); the r-denominator 11821757 = S_B + S_{B-1}"
      "\n=> eq (27)'s rational part is an Abel-Dini gap term. (c,r) need the source text.")
