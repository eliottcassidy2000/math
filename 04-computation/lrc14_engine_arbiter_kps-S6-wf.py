"""
Arbiter: which mu engine is correct for k=9 consec?  A=247/294=0.840136, B=4829/5880=0.821259.
Method: exact rational fine grid. Sample x=j/N for huge N, count fraction with maxgap>1/7.
Because mu is a union of finitely many rational intervals with denominators dividing
lcm of (e_i-e_j) and the theta-crossings, a grid at N = (big multiple) brackets the truth.
We use the EXACT maxgap_at and a deterministic dense grid plus midpoint-of-cell exact method
on a SMALL spread so the cell decomposition is unambiguous.

Decisive test: build the EXACT set of breakpoints as the union of BOTH engines' breakpoints
(superset), evaluate maxgap>theta at each cell midpoint exactly, sum. This refined partition
is a common refinement of both => it is the ground truth (any cell where maxgap-status is
constant; the union of breakpoints guarantees constancy within each cell for BOTH the order
and the theta-threshold). Compare to A and B.
"""
from fractions import Fraction as F
from math import floor

def maxgap_at(E, x):
    pts = sorted(set((F(e)*x) % 1 for e in E))
    if len(pts)==1:
        return F(1)
    g=F(0)
    for i in range(len(pts)):
        if i==len(pts)-1:
            d=(pts[0]+1)-pts[i]
        else:
            d=pts[i+1]-pts[i]
        if d>g: g=d
    return g

def breakpoints_union(E,theta):
    E=sorted(set(E)); n=len(E)
    bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(n):
            if i==j: continue
            d=E[j]-E[i]
            if d==0: continue
            # order changes: x=m/d
            for m in range(0,abs(d)+1):
                v=F(m,abs(d))
                if 0<=v<=1: bp.add(v)
            # theta crossings: (e_j-e_i)x = k+theta
            sa=abs(d)
            for k in range(-2,sa+3):
                val=(F(k)+theta)/d
                if 0<val<1: bp.add(val)
    return sorted(bp)

def mu_ground(E,theta):
    bp=breakpoints_union(E,theta)
    total=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2
        if maxgap_at(E,mid)>theta:
            total+=b-a
    return total

# Engine A
def mu_theta_A(E,theta):
    E=sorted(set(E)); n=len(E); bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1); total=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2; order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        ks=[floor(E[order[t]]*mid) for t in range(n)]; subs=[]
        for t in range(n):
            o1=order[t];o2=order[(t+1)%n];k1=ks[t];k2=ks[(t+1)%n];wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1];c=F(k1-k2+wrap)
            if s==0:
                if c>theta: subs.append((a,b))
            elif s>0:
                lo=max(a,(theta-c)/s)
                if lo<b: subs.append((lo,b))
            else:
                hi=min(b,(theta-c)/s)
                if a<hi: subs.append((a,hi))
        subs.sort(); cur=cb=None
        for lo,hi in subs:
            if cur is None: cur,cb=lo,hi
            elif lo<=cb: cb=max(cb,hi)
            else: total+=cb-cur; cur,cb=lo,hi
        if cur is not None: total+=cb-cur
    return total

theta=F(1,7)
for k in [8,9,10]:
    E=list(range(k))
    g=mu_ground(E,theta)
    a=mu_theta_A(E,theta)
    print(f"consec k={k}: GROUND={g}={float(g):.6f}  A={a}={float(a):.6f}  A==ground:{a==g}")

# also a couple random
import random
random.seed(11)
print("--- random arbiter ---")
for _ in range(8):
    k=random.randint(8,9)
    pts={0}
    while len(pts)<k: pts.add(random.randint(1,14))
    E=sorted(pts)
    g=mu_ground(E,theta); a=mu_theta_A(E,theta)
    print(E,"GROUND",float(g),"A",float(a),"eq",a==g)
