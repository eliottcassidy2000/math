from tandem import coeff
from collections import defaultdict
from math import comb
from fractions import Fraction as Q
def ck(v):
    if not v:raise ArithmeticError('check failed')
def xpoly(a,constant):
    # Q_r(m)=sum_{i=1}^r (-1)^(r-i) binom(r,i) i^m for positive m
    d=defaultdict(int);d[1,1]=constant
    for (r,s),v in a.items():
        for i in range(1,r+1):
            for j in range(1,s+1):
                d[i,j]+=(-1)**(r+s-i-j)*comb(r,i)*comb(s,j)*v
    return {z:v for z,v in d.items() if v}
def trim(a): return {z:v for z,v in a.items() if v}
def sides(a,axis,m):
    d=defaultdict(int)
    for (i,j),v in a.items():
        if axis:i,j=j,i
        d[j]+=v*i**m
    return trim(d)
def diag(a):
    d=defaultdict(int)
    for (i,j),v in a.items():d[i*j]+=v
    return trim(d)
def stop(a):
    # Constant or nonconstant polynomial in distinct positive bases, not zero.
    ck(a);b=max(a);i=1
    while abs(a[b])*b**i<=sum(abs(v)*j**i for j,v in a.items() if j<b):i+=1
    return i
def val(a,t):return sum(v*i**t for i,v in a.items())
def f(a,t,u):return sum(v*i**t*j**u for (i,j),v in a.items())
for p,q in [(12,31),(43,75)]:
    H,P=[xpoly(a,c) for a,c in zip(coeff(p,q),(-10,8))]
    if (p,q)==(12,31):
        ck(H[4,4]==18 and all(i<=4 and j<=4 for i,j in H))
        c=sum(abs(v) for z,v in H.items() if z!=(4,4));r=Q(3,4);bound=18
    else:
        ck(H[4,5]==48 and H[5,4]==-48 and all((i<=4 and j<=5) or (i<=5 and j<=4) for i,j in H))
        c=sum(abs(v) for z,v in H.items() if z not in ((4,5),(5,4)));r=Q(4,5);bound=Q(48,5)
    t0=1
    while c*r**t0>=bound:t0+=1
    roots=set();longest=0
    for z in range(1,t0):
        for axis in (0,1):
            h,p2=sides(H,axis,z),sides(P,axis,z);t1=stop(h or p2)
            longest=max(longest,t1)
            for u in range(1,t1):
                if val(h,u)==val(p2,u)==0:roots.add((u,z) if axis else (z,u))
    if q==75:
        h,p2=diag(H),diag(P)
        ck(max(h)==16 and h[16]==-84)
        t1=stop(h or p2)
        longest=max(longest,t1)
        for z in range(t0,t1):
            if val(h,z)==val(p2,z)==0:roots.add((z,z))
    print(p,q,'large quadrant threshold',t0,'maximum strip bound',longest,'solutions',roots)
    ck(roots=={(1,1)})
