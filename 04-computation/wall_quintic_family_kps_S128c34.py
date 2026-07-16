#!/usr/bin/env python3
"""Dedekind-hygienic shapes (skip p | disc) + (2,q) family robustness for the level-5 wall."""
import sys, math
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
def polmul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
def wall_poly(k,s,q):
    coeffs=[F(0)]*(k+1)
    for j in range(k+1):
        fall=[F(1)]
        for i in range(j): fall=polmul(fall,[F(-i),F(1)])
        fact=math.factorial(j)
        c=F((-1)**j*s**j,fact)*q**(k-j)
        for d,co in enumerate(fall): coeffs[d]+=c*co
    den=reduce(lambda a,b:a*b//gcd(a,b),[c.denominator for c in coeffs],1)
    ic=[int(c*den) for c in coeffs]
    g=0
    for c in ic: g=gcd(g,abs(c))
    ic=[c//g for c in ic]
    if ic[-1]<0: ic=[-c for c in ic]
    return ic
def polydiv(a,b,p):
    a=[c%p for c in a]; b=[c%p for c in b]
    while b and b[-1]==0: b.pop()
    inv=pow(b[-1],p-2,p); q=[0]*(max(0,len(a)-len(b)+1)); r=a[:]
    while len(r)>=len(b):
        if r[-1]==0: r.pop(); continue
        fq=r[-1]*inv%p; sh=len(r)-len(b); q[sh]=fq
        for i,c in enumerate(b): r[sh+i]=(r[sh+i]-fq*c)%p
        while r and r[-1]==0: r.pop()
    return q,r
def deflate(a,r,p):
    hi=a[::-1]; out=[hi[0]%p]
    for c in hi[1:]: out.append((c+out[-1]*r)%p)
    out.pop(); return out[::-1]
def shape_mod_p(a,p):
    a=[c%p for c in a]
    if a[-1]==0: return None
    def gcdp(x,y):
        while y and any(y):
            _,rr=polydiv(x,y,p); x,y=y,rr
        while x and x[-1]==0: x.pop()
        return x
    der=[(i*c)%p for i,c in enumerate(a)][1:]
    if len(gcdp(a[:],der))-1>0: return None
    shape=[]; cur=a[:]
    changed=True
    while changed and len(cur)-1>0:
        changed=False
        for r in range(p):
            v=0
            for c in reversed(cur): v=(v*r+c)%p
            if v==0: cur=deflate(cur,r,p); shape.append(1); changed=True; break
    while len(cur)-1>=4:
        found=False
        for b in range(p):
            for c in range(p):
                q,r=polydiv(cur,[c,b,1],p)
                if not r or not any(r): cur=q; shape.append(2); found=True; break
            if found: break
        if not found: break
    if len(cur)-1>0: shape.append(len(cur)-1)
    return tuple(sorted(shape))
def sylv_disc(f):
    a=f; b=[(i*c) for i,c in enumerate(f)][1:]
    n=len(a)-1; m=len(b)-1; N=n+m
    M=[[F(0)]*N for _ in range(N)]
    for i in range(m):
        for j,c in enumerate(reversed(a)): M[i][i+j]=F(c)
    for i in range(n):
        for j,c in enumerate(reversed(b)): M[m+i][i+j]=F(c)
    det=F(1)
    for col in range(N):
        piv=next((r for r in range(col,N) if M[r][col]!=0),None)
        if piv is None: return 0
        if piv!=col: M[col],M[piv]=M[piv],M[col]; det=-det
        det*=M[col][col]; inv=1/M[col][col]
        for r in range(col+1,N):
            if M[r][col]!=0:
                fq=M[r][col]*inv
                for cc in range(col,N): M[r][cc]-=fq*M[col][cc]
    return int(det/F(f[-1]))*(-1)**(n*(n-1)//2)
primes=[3,7,11,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73]
print("q  | irred | S5 | unramified shapes")
for q in (5,7,9,11,13,15,17):
    W=wall_poly(5,2,q)
    disc=sylv_disc(W)
    shapes={}
    for p in primes:
        if disc%p==0 or W[-1]%p==0: continue
        s=shape_mod_p(W,p)
        if s is None: continue
        shapes[s]=shapes.get(s,0)+1
    irred=(5,) in shapes
    s5=irred and (((2,3) in shapes) or ((1,1,1,2) in shapes))
    mark=" <== LRC(14)" if q==13 else ""
    print("%2d |  %s  | %s | %s%s"%(q,"Y" if irred else "n","Y" if s5 else "?",dict(shapes),mark))
print("DONE")
