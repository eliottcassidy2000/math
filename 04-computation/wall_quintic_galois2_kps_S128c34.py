#!/usr/bin/env python3
"""wall_quintic_galois2_kps_S128c34.py -- Galois group of the level-5 wall quintic,
brute-force mod-p factor shapes (no clever DDF), Dedekind + discriminant."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
# W5 primitive (low->high), sign-normalized to positive leading coeff
f=[-5569395, 930376, -77680, 4300, -170, 4]
print("f(x) = 4x^5 -170x^4 +4300x^3 -77680x^2 +930376x -5569395")

def polmodp(a,p): return [c%p for c in a]
def deflate(a,r,p):
    # divide by (x-r) mod p, a monic-ish ok; synthetic division low->high: reverse to high->low
    hi=a[::-1]; out=[hi[0]%p]
    for c in hi[1:]:
        out.append((c+out[-1]*r)%p)
    rem=out.pop()
    return out[::-1],rem
def polydiv(a,b,p):
    # a,b low->high; returns (q,r)
    a=[c%p for c in a]; b=[c%p for c in b]
    while b and b[-1]==0: b.pop()
    inv=pow(b[-1],p-2,p); q=[0]*(max(0,len(a)-len(b)+1))
    r=a[:]
    while len(r)>=len(b):
        if r[-1]==0: r.pop(); continue
        fq=r[-1]*inv%p; sh=len(r)-len(b)
        q[sh]=fq
        for i,c in enumerate(b): r[sh+i]=(r[sh+i]-fq*c)%p
        while r and r[-1]==0: r.pop()
    return q,r
def shape_mod_p(a,p):
    """factor shape of squarefree quintic mod p by brute force."""
    a=polmodp(a,p)
    if a[-1]%p==0: return None
    # squarefree check via gcd(f,f') by brute: skip -- detect repeated roots during root extraction; for deg2 rep factors use resultant... simpler: compute gcd via Euclid
    def gcdp(x,y):
        x=[c%p for c in x]; y=[c%p for c in y]
        while y and any(y):
            _,rr=polydiv(x,y,p); x,y=y,rr
        while x and x[-1]==0: x.pop()
        return x
    der=[(i*c)%p for i,c in enumerate(a)][1:]
    g=gcdp(a,der)
    if len(g)-1>0: return None  # not squarefree, skip prime
    shape=[]; cur=a[:]
    # linear factors
    changed=True
    while changed and len(cur)-1>0:
        changed=False
        for r in range(p):
            # eval
            v=0
            for c in reversed(cur): v=(v*r+c)%p
            if v==0:
                cur,rem=deflate(cur,r,p); assert rem==0
                shape.append(1); changed=True; break
    d=len(cur)-1
    if d==0: return tuple(sorted(shape))
    # quadratic factors by trial division (monic x^2+bx+c)
    found=True
    while found and len(cur)-1>=4:
        found=False
        for b in range(p):
            for c in range(p):
                q,r=polydiv(cur,[c,b,1],p)
                if not r or not any(r):
                    cur=q; shape.append(2); found=True; break
            if found: break
    d=len(cur)-1
    if d>0:
        # remaining is irreducible of degree d (no linear roots, no quadratic divisors => deg 3,4,5 irreducible... deg 4 could be product of two irreducible quadratics -- handled above; deg 4 irreducible fine)
        shape.append(d)
    return tuple(sorted(shape))

primes=[3,7,11,17,19,23,29,31,37,41,43,47,53,59]
shapes={}
for p in primes:
    s=shape_mod_p(f,p)
    if s is None: continue
    shapes[s]=shapes.get(s,0)+1
print("Frobenius shapes:",shapes)
irred=(5,) in shapes
print("irreducible over Q:", "YES (irreducible mod p for some p)" if irred else "NOT CERTIFIED -- check factorizations")
hasS5=((2,3) in shapes) or ((1,1,1,2) in shapes) or ((1,4) in shapes and (1,1,3) in shapes)
# discriminant via Sylvester determinant with fractions
def sylvester_det(a,b):
    n=len(a)-1; m=len(b)-1
    N=n+m
    M=[[F(0)]*N for _ in range(N)]
    for i in range(m):
        for j,c in enumerate(reversed(a)): M[i][i+j]=F(c)
    for i in range(n):
        for j,c in enumerate(reversed(b)): M[m+i][i+j]=F(c)
    det=F(1)
    for col in range(N):
        piv=next((r for r in range(col,N) if M[r][col]!=0),None)
        if piv is None: return F(0)
        if piv!=col: M[col],M[piv]=M[piv],M[col]; det=-det
        det*=M[col][col]
        inv=1/M[col][col]
        for r in range(col+1,N):
            if M[r][col]!=0:
                fq=M[r][col]*inv
                for cc in range(col,N): M[r][cc]-=fq*M[col][cc]
    return det
der=[(i*c) for i,c in enumerate(f)][1:]
res=sylvester_det(f,der)
disc=res/F(f[-1])
disc=disc*(-1)**(5*4//2)
disc=int(disc)
print("disc(f) =",disc)
import math
a=abs(disc); r=math.isqrt(a)
print("|disc| perfect square?",r*r==a,"; sign","+" if disc>0 else "-")
# factor disc smallish
dd=a; fac={}
d=2
while d*d<=dd and d<200000:
    while dd%d==0: fac[d]=fac.get(d,0)+1; dd//=d
    d+=1
print("disc factorization (partial):",fac,"cofactor",dd)
if irred and hasS5:
    print("VERDICT: Galois group = S5 -- the level-5 wall location 15.6306... is NOT solvable by radicals")
elif irred and (1,1,3) in shapes and r*r==a:
    print("VERDICT: A5")
elif irred:
    print("VERDICT: solvable transitive group (C5/D5/F20); shapes:",sorted(shapes))
else:
    print("VERDICT: reducible -- wall root may still be radical-expressible; inspect")
print("DONE")
