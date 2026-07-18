from fractions import Fraction as F
from functools import reduce
from math import gcd
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def uncovered_exact(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    u=union(allv)
    return 1-sum(b-a for a,b in u), u
V=[1,2,3,4,5,6,7,8,9,10,11,13,24]
u,comp=uncovered_exact(V)
print("VERIFYING the candidate second tight family")
print(f"  V = {V}   (13 speeds, distinct: {len(set(V))==13}, gcd = {reduce(gcd,V)})")
print(f"  exact uncovered = {u}   float = {float(u)}")
print(f"  union has {len(comp)} component(s): {comp[:3]}{'...' if len(comp)>3 else ''}")
print()
print("  CROSS-CHECK by the rational-point test: is ANY p/q lonely, q <= 400?")
def lonely_at(V,p,q):
    for v in V:
        r=(v*p)%q
        if min(r,q-r)*14 < q: return False
    return True
hits=[(p,q) for q in range(1,401) for p in range(1,q) if gcd(p,q)==1 and lonely_at(V,p,q)]
print(f"    lonely rationals with q <= 400: {hits[:6]}{'...' if len(hits)>6 else ''}  (count {len(hits)})")
print()
print("  CONTRAST with the known tight family {1,...,13}:")
u2,c2=uncovered_exact(list(range(1,14)))
hits2=[(p,q) for q in range(1,60) for p in range(1,q) if gcd(p,q)==1 and lonely_at(list(range(1,14)),p,q)]
print(f"    uncovered = {u2}, lonely rationals q<=60: {hits2[:6]}  (count {len(hits2)})")
print()
print("  SUP OF THE GAP: for a tight family the sup of min_v ||v t|| should be exactly 1/14.")
print("  Sample min_v ||v t|| at the candidate lonely points of each family:")
def mingap(V,t):
    return min(min(F(v*t)%1, 1-F(v*t)%1) for v in V)
for nm,W,pts in [("{1..13}",list(range(1,14)),[F(1,14)]),
                 ("candidate",V,[F(p,q) for p,q in hits[:3]] or [F(1,14)])]:
    for t in pts:
        print(f"    {nm:12s} t={t}:  min_v ||v t|| = {mingap(W,t)} = {float(mingap(W,t)):.8f}  (1/14 = {float(F(1,14)):.8f})")
