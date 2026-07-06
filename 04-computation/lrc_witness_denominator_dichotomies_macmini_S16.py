#!/usr/bin/env python3
"""
[FIX S20]: removed the buggy `if gcd(a,q)!=1: continue` -- it skipped witnesses at
denominators that DIVIDE a pairwise sum/diff (e.g. q=11 via 22=4+18 as the non-coprime 2a/22),
underestimating M.  Now checks ALL a. See MISTAKE log.
mac-mini-2026-07-06-S16 -- LEVERAGING THE DICHOTOMIES on the (G) crux.

The tight witness is t=1/13 (on the 13-grid; lifts invisible there -- kps S22
character-sum bridge).  The open density floor lives at OFF-13-GRID witnesses.
M(S) reduced = c/q means the witness denominator is q.  This script extracts q
across the loneliness spectrum and stratifies it by the DICHOTOMIES the owner
named -- odd/even, add/mult (factorization, 13|q, small-prime structure),
pos/neg (reflection symmetry of the phase config) -- to see WHICH denominators
carry the near-tight values and what constrains them.

KEY DICHOTOMY LEVERAGE (developed):
 * ODD/EVEN (the '2' of 14=2.7): at an EVEN witness q=2p, an even runner 2w has
   phase (wa mod p)/p -- so the EVEN runners, HALVED, are a config at modulus p
   with clearance ceil(c/2).  Owner's seed E_p={0,+-2}, O_p={+-1}: for c=3 the
   halved even config avoids {0,+-1} mod p (clearance 2).  A self-similar DESCENT
   q -> q/2 preserving the value.  Depth is bounded by covering (only a multiple
   of 8=2^3 is guaranteed) => at most 3 halvings.
 * ADD/MULT: q's factorization.  For the 12-base the boundary is 2/25 (q=25=5^2)
   and the deep well 14/169 (q=169=13^2) -- the extremal denominators are prime
   SQUARES of covering-relevant primes.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import sys
sys.path.insert(0,'04-computation')
from lonely_profile import profile

def log(m=""): print(m,flush=True)
def M_exact(S):
    S=sorted(set(abs(x) for x in S))
    for cap in (14,11,8,6,4,3,2):
        p=profile(S,F(1,cap)); m=p.M()
        if m is not None: return m
def factor(n):
    f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
def fstr(f): return "*".join(f"{p}^{e}" if e>1 else f"{p}" for p,e in sorted(f.items())) or "1"

spectrum={
 "AP {1..12}":            list(range(1,13)),
 "doubled-apex {1..11,24}":list(range(1,12))+[24],
 "block {..17,19}":       [1,2,3,5,7,8,9,10,11,12,17,19],
 "{1..11,23}":            list(range(1,12))+[23],
 "{1..11,13}":            list(range(1,12))+[13],
 "deep well {1..11,168}": list(range(1,12))+[168],
 "{1..10,12,25}":         list(range(1,11))+[12,25],
 "{2..13}":               list(range(2,14)),
 "{1..11,14}":            list(range(1,12))+[14],
}
log("family                       M=c/q       q      q factored   even?  13|q  5|q")
log("-"*82)
for name,S in spectrum.items():
    M=M_exact(S); q=M.denominator; c=M.numerator; f=factor(q)
    log(f"{name:28s} {str(M):>7}={float(M):.4f} {q:>5}   {fstr(f):12s} "
        f"{'EVEN' if q%2==0 else 'odd ':4s}  {'yes' if q%13==0 else 'no ':3s}  {'yes' if q%5==0 else 'no'}")

log("\n=== ODD/EVEN DESCENT verification (even-q witness => even runners halve to clearance ceil(c/2) at q/2):")
# construct a family with an EVEN witness and verify the descent
for name,S in spectrum.items():
    M=M_exact(S); q=M.denominator; c=M.numerator
    if q%2==0:
        p=q//2
        # find witness a: some a with min_i ||v_i a/q|| = c/q
        awit=None
        for a in range(1,q):
            mn=min(min((v*a)%q, q-((v*a)%q)) for v in S)
            if mn==c: awit=a; break
        if awit is None:
            log(f"  {name}: q={q} even, witness a not found at clearance {c}"); continue
        evens=[v for v in S if v%2==0]
        halved=[v//2 for v in evens]
        # halved config phases at modulus p=q/2 with witness a
        clhalved=min(min((w*awit)%p, p-((w*awit)%p)) for w in halved) if halved else None
        log(f"  {name}: q={q}=2*{p}, c={c}, witness a={awit}; {len(evens)} even runners; "
            f"halved clearance at mod {p} = {clhalved} (predicted ceil({c}/2)={-(-c//2)})")

log("\nREADING: the near-tight spectrum lives on denominators that are prime-square /")
log("small-prime (5,13) structured; the parity descent turns an even witness into a")
log("HALF-modulus config (the '2' of 14=2.7 acting as renormalization).  The dichotomies")
log("STRATIFY the off-grid density floor: 13|q (character-sum/pinning), q even (halve),")
log("q coprime-to-26 (the generic residual).")
