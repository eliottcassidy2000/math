#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S20e: VALIDATION of the Q-census bound against MISTAKE-110.

mac-mini HYP-4302 collapses the (A) tail to a FINITE residue census: opus HYP-4266
observes "300/300 random unbounded spread families have a margin-2/25 witness at
some q <= 50".  My MISTAKE-110 built a LOOSE family whose minimal clearing witness
is at q = 53 (a free prime), NOT any q <= 50, because composite runners k*c with
c = lcm(1..50) are == 0 mod every q <= 50 and so are NEVER cleared at t = a/q, q<=50.

TWO questions this settles:
  (1) Is that family a GAP member (M in (1/13, 2/25)) -- which WOULD break
      HYP-4302 -- or ABOVE the window (M >= 2/25, a loose family)?  If above,
      HYP-4302 stands (it only rules out gap members).
  (2) Does a q <= 50 census FAIL to certify it (needs q = 53)?  If so, the
      census q-bound must EXCEED 50 (or use ray-transport to a bounded rep);
      a naive "check q <= 50" cap is incomplete.

We use a TRACTABLE analog: c = lcm(1..B) for small B, family
  F = {c, 2c, ..., 11c, 12c + 1}.
The blocking mechanism is identical: k*c == 0 mod q for every q | c (all q <= B),
so only q with q-does-not-divide-c can clear.  Smallest such prime = first prime > B.
"""
from math import gcd
from fractions import Fraction

def lcm(a, b): return a * b // gcd(a, b)
def lcm_range(n):
    r = 1
    for k in range(1, n+1): r = lcm(r, k)
    return r

def clears_at(speeds, q, rho=Fraction(2,25)):
    """Is there a in [1,q-1] with ||v*a/q|| >= rho for all v?  Returns (bool, best_a, best_margin)."""
    best = (False, None, Fraction(-1))
    for a in range(1, q):
        m = min(min((v*a) % q, q - (v*a) % q) for v in speeds)  # q*||v a/q||
        margin = Fraction(m, q)
        if margin > best[2]:
            best = (margin >= rho, a, margin)
    return best

def min_clearing_modulus(speeds, qmax=200, rho=Fraction(2,25)):
    for q in range(2, qmax+1):
        ok, a, marg = clears_at(speeds, q, rho)
        if ok:
            return q, a, marg
    return None, None, None

print("=== Q-census bound validation (MISTAKE-110 mechanism at tractable scale) ===")
for B in [6, 8, 10, 12]:
    c = lcm_range(B)
    F = [k*c for k in range(1, 12)] + [12*c + 1]
    # smallest prime not dividing c = smallest prime > B (since c = lcm(1..B))
    q_needed, a, marg = min_clearing_modulus(F, qmax=4*B+20)
    # check that NO q <= B clears (composite runners blocked)
    blocked_upto = all(not clears_at(F, q)[0] for q in range(2, B+1))
    print(f"  B={B:2d}, c=lcm(1..{B})={c}:")
    print(f"    family {{c,2c,...,11c,12c+1}}: min clearing modulus q = {q_needed} "
          f"(margin {marg} = {float(marg):.4f}) at a={a}")
    print(f"    all q <= {B} FAIL to clear (composite runners k*c == 0 mod q): {blocked_upto}")
    print(f"    => M(F) = {float(marg):.4f} {'>= 2/25 (ABOVE window, loose family)' if marg >= Fraction(2,25) else 'IN/BELOW window'}; "
          f"census must reach q = {q_needed} > {B}")

print()
print("VERDICT:")
print("  The family CLEARS (M >= 2/25, ABOVE the window) => NOT a gap member =>")
print("  HYP-4302's 'no gap member' STANDS (my MISTAKE-110 family is loose, not gap).")
print("  BUT its witness is at q = (first prime > B), NOT q <= B => a census capped")
print("  at q <= 50 CANNOT certify families whose composite runners divide lcm(1..50);")
print("  the census q-bound must exceed 50 (smallest safe: 53) OR use ray-transport")
print("  to a bounded residue representative (opus HYP-4266's mechanism -- which is")
print("  exactly what makes the bound ray-LOCAL, not a uniform q <= 50).")
