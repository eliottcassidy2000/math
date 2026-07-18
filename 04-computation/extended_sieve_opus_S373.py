# opus-2026-07-17-S373 -- HYP-7560: THE EXTENDED RATIONAL-POINT SIEVE.
#
# THE OBSERVATION THE SEVEN-MODULI SIEVE UNDER-USES.  t = p/q is a lonely
# point iff  ||v p / q|| >= 1/14  for EVERY speed v, i.e. writing r = vp mod q,
#        min(r, q-r) * 14 >= q.
# This depends ONLY ON v mod q.  So "does a lonely rational with denominator
# <= Q exist" is a condition on the speeds mod lcm(1..Q) -- a FINITE condition,
# not an analytic one.  That is a different KIND of problem from the uniformity
# that blocked the ledger (THM-1095).
#
# The classical sieve is the special case q <= 14, where the band condition
# degenerates to "q divides no speed" (for q <= 14 every nonzero residue
# passes).  For q > 14 the band [q/14, 13q/14] is a genuine constraint but a
# WIDE one -- a random speed fails with probability only ~1/7 -- and there are
# phi(q) numerators p to try for each q.  So the sieve should get much stronger.
# THE TEST: how does the kill rate grow with Q?
from math import gcd
import random

def lonely_at(V, p, q):
    for v in V:
        r = (v*p) % q
        if min(r, q-r)*14 < q: return False
    return True

def min_denominator(V, Qmax):
    """smallest q admitting a lonely p/q, or None."""
    for q in range(1, Qmax+1):
        for p in range(1, q):
            if gcd(p,q)!=1: continue
            if lonely_at(V,p,q): return q,p
    return None

def blocks_seven(V):
    return all(any(v % q == 0 for v in V) for q in range(8,15))

print("(1) KILL RATE vs Q  (random comparable 13-families)")
random.seed(373)
fams=[sorted(random.sample(range(2,400),13)) for _ in range(400)]
for Q in [14,20,30,50,80,120,200]:
    killed=sum(1 for V in fams if min_denominator(V,Q))
    print(f"    Q <= {Q:4d}:  {killed:4d}/400 = {100*killed/400:5.1f}%")

print()
print("(2) THE HARD STRATUM -- families that BLOCK all seven classical moduli")
hard=[V for V in fams if blocks_seven(V)]
print(f"    {len(hard)} of 400 block all seven (the classical ~11% residue)")
for Q in [14,20,30,50,80,120,200]:
    killed=sum(1 for V in hard if min_denominator(V,Q))
    print(f"      Q <= {Q:4d}:  {killed:3d}/{len(hard)} = {100*killed/max(len(hard),1):5.1f}%")

print()
print("(3) THE NAMED ADVERSARIAL FAMILIES")
named=[("THM-1055 primitive failure",[27,36,46,70,101,114,117,121,140,160,194,277,293]),
       ("{1,...,13} (tight)",list(range(1,14))),
       ("2*{1,...,13} (dilate)",[2*i for i in range(1,14)]),
       ("AP d=8",[1+8*i for i in range(13)]),
       ("odd {1,3,...,25}",[2*i+1 for i in range(13)]),
       ("THM-1060 L=31 family",[248+8*i for i in range(13)])]
for nm,V in named:
    r=min_denominator(V,400)
    print(f"    {nm:30s} -> {('lonely at p/q = %d/%d'%(r[1],r[0])) if r else 'NONE with q <= 400'}")
