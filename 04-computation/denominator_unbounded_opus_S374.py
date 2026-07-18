# opus-2026-07-17-S374 -- HYP-7570: IS THE BOUNDED-DENOMINATOR CONJECTURE TRUE?
#
# THM-1100 conjectured an absolute Q0 with every primitive 13-family admitting a
# lonely p/q, q <= Q0.  I should have checked one construction before proposing it.
#
# THE CONSTRUCTION.  Blocking modulus q needs only ONE speed v with q | v: then
# at t = p/q that runner sits exactly at the origin, ||v p/q|| = 0 < 1/14, for
# EVERY p.  So a SINGLE speed divisible by lcm(1..Q) blocks EVERY q <= Q at once.
# Take V = {lcm(1..Q)} u {12 speeds coprime to it}: primitive, 13 speeds, and no
# lonely rational of denominator <= Q.  Q is arbitrary, so the minimal
# denominator is UNBOUNDED and the conjecture is FALSE.
# This is exactly why my adversarial searches kept climbing 25 -> 32 -> 39: the
# supremum is infinite, and the hill-climbs were crawling toward it.
from math import gcd
from functools import reduce
def lcm(a,b): return a*b//gcd(a,b)
def lcm_range(n): return reduce(lcm, range(1,n+1))
def lonely_at(V,p,q):
    for v in V:
        r=(v*p)%q
        if min(r,q-r)*14 < q: return False
    return True
def min_den(V,Qmax):
    for q in range(1,Qmax+1):
        for p in range(1,q):
            if gcd(p,q)!=1: continue
            if lonely_at(V,p,q): return q
    return None

print("(1) THE REFUTATION: one speed divisible by lcm(1..Q) blocks every q <= Q")
print("     Q    lcm(1..Q)          gcd(V)  min denominator (searched to 3Q)")
for Q in [10,15,20,25,30,40]:
    L=lcm_range(Q)
    V=sorted([L]+[p for p in [3,5,11,13,17,19,23,29,31,37,41,43]])
    g=reduce(gcd,V)
    md=min_den(V,3*Q)
    print(f"    {Q:4d}  {L:16d}   {g:5d}   {md if md else '> '+str(3*Q)}")

print()
print("(2) CHECK: is the blocking really by divisibility, as claimed?")
Q=30; L=lcm_range(Q)
V=sorted([L]+[3,5,11,13,17,19,23,29,31,37,41,43])
print(f"    V has max speed {max(V)} and gcd {reduce(gcd,V)} (primitive)")
bad=[q for q in range(1,Q+1) if L % q != 0]
print(f"    moduli q <= {Q} NOT dividing the big speed: {bad}  (empty => all blocked)")
print(f"    min denominator for this V, searched to 200: {min_den(V,200)}")

print()
print("(3) SO WHAT SURVIVES?  the minimal denominator vs the DIVISIBILITY of V")
print("    define D(V) = largest Q with every q <= Q dividing some speed.")
print("    the construction forces min-denominator > D(V).  Is it ~ D(V)?")
import random
random.seed(374)
for trial in range(8):
    V=sorted(random.sample(range(2,200),13))
    if reduce(gcd,V)!=1: continue
    D=0
    for q in range(1,60):
        if any(v%q==0 for v in V): D=q
        else: break
    md=min_den(V,300)
    print(f"    D(V) = {D:3d}   min denominator = {md}   V[:5] = {V[:5]}...")
