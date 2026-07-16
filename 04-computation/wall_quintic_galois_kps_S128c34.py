#!/usr/bin/env python3
"""wall_quintic_galois_kps_S128c34.py -- kind-pasteur S128 cont.34.
The Bonferroni wall polynomial W_k(x) = sum_{j=0}^{k} (-1)^j C(x,j) 2^j 13^{k-j}:
exact coefficients, wall locations, irreducibility, Galois group of the level-5 wall
(Frobenius/Dedekind shapes), discriminant, and the (2,q)-family robustness."""
import sys
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def polmul(a,b,p=None):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            r[i+j]+=x*y
    return [c%p for c in r] if p else r

def wall_poly(k, s=2, q=13):
    # W_k(x) = sum_j (-1)^j C(x,j) s^j q^(k-j); C(x,j) = x(x-1)..(x-j+1)/j!
    coeffs=[F(0)]*(k+1)
    for j in range(k+1):
        # falling factorial poly
        fall=[F(1)]
        for i in range(j):
            fall=polmul(fall,[F(-i),F(1)])
        fact=1
        for i in range(2,j+1): fact*=i
        c=F((-1)**j * s**j, fact) * q**(k-j)
        for d,co in enumerate(fall):
            coeffs[d]+=c*co
    # clear denominators -> primitive integer polynomial
    from functools import reduce
    den=reduce(lambda a,b: a*b//gcd(a,b), [c.denominator for c in coeffs],1)
    ic=[int(c*den) for c in coeffs]
    g=0
    for c in ic: g=gcd(g,abs(c))
    return [c//g for c in ic]

def polmod(a,m,p):
    a=[c%p for c in a]
    dm=len(m)-1; inv=pow(m[-1],p-2,p)
    while len(a)-1>=dm and any(a):
        if a[-1]==0: a.pop(); continue
        f=a[-1]*inv%p; sh=len(a)-1-dm
        for i,c in enumerate(m): a[sh+i]=(a[sh+i]-f*c)%p
        while a and a[-1]==0: a.pop()
    return a if a else [0]

def polpowmod(base,e,m,p):
    r=[1]; b=polmod(base,m,p)
    while e:
        if e&1: r=polmod(polmul(r,b,p),m,p)
        b=polmod(polmul(b,b,p),m,p); e>>=1
    return r

def polgcd(a,b,p):
    a=[c%p for c in a]; b=[c%p for c in b]
    while b and any(b):
        a,b=b,polmod(a,b,p)
        while b and b[-1]==0: b.pop()
    if not a: return [1]
    inv=pow(a[-1],p-2,p)
    return [c*inv%p for c in a]

def ddf_shape(f,p):
    """distinct-degree factorization shape of squarefree f mod p: list of (d, total_degree_of_deg-d_part)."""
    f=[c%p for c in f]
    # check squarefree: gcd(f,f')
    fp=[(i*c)%p for i,c in enumerate(f)][1:]
    if len(polgcd(f,fp,p))-1>0: return None
    shape=[]; g=f[:]
    d=0; x=[0,1]
    while len(g)-1>0:
        d+=1
        h=polpowmod(x,p**d,g,p)
        hh=h[:]
        # h - x
        while len(hh)<2: hh.append(0)
        hh[1]=(hh[1]-1)%p
        while hh and hh[-1]==0: hh.pop()
        if not hh: hh=[0]
        c=polgcd(g,hh,p)
        deg=len(c)-1
        if deg>0:
            shape += [d]*(deg//d)
            # g //= c  (exact division mod p)
            q=[]; r=g[:]
            inv=pow(c[-1],p-2,p)
            while len(r)>=len(c) and any(r):
                if r[-1]==0: r.pop(); continue
                fq=r[-1]*inv%p; sh=len(r)-len(c)
                q.append((sh,fq))
                for i,cc in enumerate(c): r[sh+i]=(r[sh+i]-fq*cc)%p
                while r and r[-1]==0: r.pop()
            gg=[0]*(max(s for s,_ in q)+1)
            for s,fq in q: gg[s]=fq
            g=gg
    if len(g)-1>0: shape.append(len(g)-1)
    return tuple(sorted(shape))

def resultant(a,b):
    # exact integer resultant via Euclid with fractions
    a=[F(c) for c in a]; b=[F(c) for c in b]
    r=F(1); s=1
    while len(b)-1>0:
        da,db=len(a)-1,len(b)-1
        # a = q*b + rem
        rem=a[:]
        while len(rem)-1>=db and any(rem):
            if rem[-1]==0: rem.pop(); continue
            f=rem[-1]/b[-1]; sh=len(rem)-1-db
            for i,c in enumerate(b): rem[sh+i]-=f*c
            while rem and rem[-1]==0: rem.pop()
        if not rem or not any(rem): return 0
        dr=len(rem)-1
        r*= b[-1]**(da-dr) * (-1)**(da*db)
        a,b=b,rem
    return r * b[0]**(len(a)-1)

primes=[3,5,7,11,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113]
for k in range(1,6):
    W=wall_poly(k)
    print("W_%d primitive coeffs (low->high): %s"%(k,W))
    # wall location: largest real root (numeric)
    import cmath
    # numeric roots by companion-matrix-free Durand-Kerner
    n=len(W)-1
    ws=[complex(0.4+0.9j)**i*(1+abs(W[0]/W[-1])**(1/n)) for i in range(n)]
    for _ in range(500):
        new=[]
        for i,wi in enumerate(ws):
            num=sum(W[j]*wi**j for j in range(n+1))
            den=W[-1]
            for j,wj in enumerate(ws):
                if j!=i: den*= (wi-wj)
            new.append(wi-num/den)
        ws=new
    reals=sorted(w.real for w in ws if abs(w.imag)<1e-8)
    print("   roots: %s"%["%.6f"%w.real+("%+.6fi"%w.imag if abs(w.imag)>1e-8 else "") for w in sorted(ws,key=lambda z:(abs(z.imag)>1e-8,z.real))])
    if k==5:
        shapes={}
        irred=False
        for p in primes:
            if W[-1]%p==0: continue
            sh=ddf_shape(W,p)
            if sh is None: continue
            shapes[sh]=shapes.get(sh,0)+1
            if sh==(5,): irred=True
        print("   Frobenius shapes over %d primes: %s"%(sum(shapes.values()),shapes))
        print("   irreducible over Q:", "YES (irreducible mod some p)" if irred else "check needed")
        disc=resultant(W,[(i)*c for i,c in enumerate(W)][1:])
        lead=F(W[-1])
        disc=disc/lead
        disc=int(disc)*(-1)**(5*4//2)
        print("   discriminant:", disc)
        # square?
        import math
        a=abs(disc); r=math.isqrt(a)
        print("   |disc| square? ", r*r==a, "(sign %+d)"%(1 if disc>0 else -1))
        # verdict
        has23=any(s==(2,3) for s in shapes)
        has1112=any(s==(1,1,1,2) for s in shapes)
        has113=any(s==(1,1,3) for s in shapes)
        has14=any(s==(1,4) for s in shapes)
        if irred and (has23 or has1112):
            print("   GALOIS GROUP: S5 (transposition/(2,3)-class present) -> WALL LOCATION NOT SOLVABLE BY RADICALS")
        elif irred and has113 and r*r==a:
            print("   GALOIS GROUP: A5 -> unsolvable")
        elif irred:
            print("   GALOIS GROUP: solvable transitive (C5/D5/F20) -- shapes %s"%list(shapes))
print("== family robustness: level-5 wall for (s=2, q) ==")
for q in (5,7,9,11,13,15,17):
    W=wall_poly(5,2,q)
    shapes=set()
    irred=False
    for p in primes:
        if W[-1]%p==0 or q%p==0: continue
        sh=ddf_shape(W,p)
        if sh is None: continue
        shapes.add(sh)
        if sh==(5,): irred=True
    s5 = irred and (((2,3) in shapes) or ((1,1,1,2) in shapes))
    print("  q=%2d: irred=%s S5=%s shapes=%s"%(q,irred,s5,sorted(shapes)))
print("DONE")
