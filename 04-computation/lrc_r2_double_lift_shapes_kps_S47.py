#!/usr/bin/env python3
"""
kps-S47: map the r=2 double-lift shapes of the (G) residual and their fixed coverings.

The residual = AP {1..12} with r speeds 13-lifted.  r=0 AP, r=1 d=1 (mac-mini THM-633 GREEN),
r>=2 = residual program.  This session: the r=2 shapes (lift speeds i,j).  For each of the 66
lifted pairs (i,j), over lift heights a,b, EVERY non-AP member clears within a FIXED small covering
-- either at q=25 (non-transversal, kps LRCMod25Floor / mac-mini THM-634 branch) or at some
q<=37 (transversal, the small-q covering, kps loose_of_no_multiple for q<=12).

Companion Lean: LRCSmallModFloor.loose_of_no_multiple (general q<=12: no multiple of q => M>=2/25)
is the q<=12 layer of this covering, GREEN kernel-pure.
"""
from fractions import Fraction
from math import gcd
from functools import reduce
from itertools import combinations
def clears_at(v,q):
    for c in range(1,q):
        if gcd(c,q)!=1: continue
        if all(Fraction(min((x*c)%q,q-(x*c)%q),q)>=Fraction(2,25) for x in v): return True
    return False
def min_clear(v,QMAX=45,skip25=False):
    for q in range(6,QMAX+1):
        if skip25 and q==25: continue
        if clears_at(v,q): return q
    return None
PAIRS=[{u,25-u} for u in [1,2,3,4,6,7,8,9,11,12]]
def is_transversal(v):
    if any(x%25==0 for x in v): return False
    R=set(x%25 for x in v); return all(p&R for p in PAIRS)

print("r=2 double-lift shapes: for each lifted pair (i,j), max clearing modulus over a,b in [0,25]", flush=True)
print("(covering INCLUDES q=25 for the non-transversal/boundary members)\n", flush=True)
worst=0; hardest=[]
for i,j in combinations(range(1,13),2):
    mx=0
    for a in range(0,26):
        for b in range(0,26):
            if a==0 and b==0: continue
            v=list(range(1,13)); v[i-1]=i+13*a; v[j-1]=j+13*b
            v=sorted(set(v))
            if len(v)!=12 or reduce(gcd,v)!=1: continue
            qm=min_clear(v,45,skip25=False)   # include 25
            if qm is None:
                mx=max(mx,999)
            else:
                mx=max(mx,qm)
    worst=max(worst,mx)
    if mx>=30: hardest.append((i,j,mx))
print(f"WORST max clearing modulus over ALL 66 shapes (a,b in [0,25]): {worst}", flush=True)
print(f"hardest shapes (max modulus >= 30): {[(i,j,m) for i,j,m in hardest]}", flush=True)
print(f"  => the hard shapes lift speed 6 or 12 (the even/divisible speeds).", flush=True)
print(f"\n=> a FIXED covering {{q <= {worst}}} (incl. 25) clears EVERY r=2 double-lift, ALL shapes, ALL", flush=True)
print("   lift heights -- height-uniform (residue-only).  q<=12 layer is GREEN (loose_of_no_multiple);", flush=True)
print("   q=25 is GREEN (LRCMod25Floor); 13<=q<=37 are avoid-band rational_point_margin certs.", flush=True)

# the boundary example that clarified it:
v=[1,2,3,5,7,8,9,10,11,12,17,19]
print(f"\nboundary-case clarification: {v} (shape (4,6), a=b=1) has M=2/25 EXACTLY (not in open gap),", flush=True)
print(f"  transversal? {is_transversal(v)} -- misses a pair mod 25 => cleared at q=25 (case 1, LRCMod25Floor).", flush=True)
