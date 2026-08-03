#!/usr/bin/env python3
"""Exact reflected-channel arithmetic used by THM-3211's certificate.

This is a dependency-free floor-sum engine, a stabilized-ray scout, and an
exact evaluator for the periodic Bernoulli formula.  Its main block is a
finite hostile control, not the proof of the theorem's uniform quantifiers.
"""

from fractions import Fraction as F
from math import gcd

L=168
CELL=90
H=(1,2,3,4,6,12)
EDGES=((0,1),(0,2),(0,3),(0,4),(0,5),(1,2),(1,3),(1,4),(1,5))

def floor_moments(n,m,a,b):
    if n==0: return (0,0,0)
    s1=n*(n-1)//2
    s2=n*(n-1)*(2*n-1)//6
    qa,a0=divmod(a,m); qb,b0=divmod(b,m)
    base0=qa*s1+qb*n
    base1=qa*s2+qb*s1
    base2=qa*qa*s2+2*qa*qb*s1+qb*qb*n
    if a0==0: return base0,base1,base2
    height=(a0*(n-1)+b0)//m
    if height==0: return base0,base1,base2
    u0,u1,u2=floor_moments(height,a0,m,m-b0+a0-1)
    r0=n*height-u0
    r1=height*s1-(u2-u0)//2
    r2=n*height*height-2*u1-u0
    return base0+r0, base1+r1, base2+2*qa*r1+2*qb*r0+r2

def residue_prefix(n,m,a,b,threshold,base):
    shifted=floor_moments(n,m,a,b+m-threshold)
    d0=shifted[0]-base[0]; d1=shifted[1]-base[1]
    y0d=(shifted[2]-base[2]-d0)//2
    high_sum=a*d1+b*d0-m*y0d
    total=a*n*(n-1)//2+b*n-m*base[0]
    return n-d0,total-high_sum

def triangle_sum(n,m,a,b,peak,base,total):
    if peak<=0: return 0
    radius=(peak-1)//L
    if 2*radius>=m: raise RuntimeError((n,m,peak))
    lc,ls=residue_prefix(n,m,a,b,radius+1,base)
    bc,bs=residue_prefix(n,m,a,b,m-radius,base)
    hc,hs=n-bc,total-bs
    return peak*lc-L*ls+(peak-L*m)*hc+L*hs

def mass(e,p,f,q):
    if p>q: return mass(f,q,e,p)
    z,w=L*p-e,L*q-f
    r,s=e*CELL%L,f*CELL%L
    det=r*w-s*z
    if det%L: raise RuntimeError((e,p,f,q,det))
    b=(det//L)%z; a=w%z
    base=floor_moments(p,z,a,b)
    total=a*p*(p-1)//2+b*p-z*base[0]
    num=triangle_sum(p,z,a,b,12*(z+w),base,total)-triangle_sum(p,z,a,b,12*(w-z),base,total)
    return F(num,z*w)

def cleared_num(e,p,f,q):
    x=mass(e,p,f,q)
    return x.numerator*((L*p-e)*(L*q-f)//x.denominator)

def interpolate_g_poly(points):
    # Three exact (g,N) points -> coefficients in g.
    (x0,y0),(x1,y1),(x2,y2)=points
    den=(x0-x1)*(x0-x2)
    a=F(y0,den)+F(y1,(x1-x0)*(x1-x2))+F(y2,(x2-x0)*(x2-x1))
    b=-F(y0*(x1+x2),den)-F(y1*(x0+x2),(x1-x0)*(x1-x2))-F(y2*(x0+x1),(x2-x0)*(x2-x1))
    c=F(y0*x1*x2,den)+F(y1*x0*x2,(x1-x0)*(x1-x2))+F(y2*x0*x1,(x2-x0)*(x2-x1))
    return a,b,c

def exact_period_scout(e,P,f,Q,h0=20):
    M=abs(Q*e-P*f) or 1
    pol=[]
    for r in range(1,M+1):
        pts=[]
        for j in range(4):
            g=r+M*(h0+j)
            pts.append((g,cleared_num(e,g*P,f,g*Q)))
        q=interpolate_g_poly(pts[:3])
        if sum(q[i]*F(pts[3][0])**(2-i) for i in range(3))!=pts[3][1]:
            raise RuntimeError(('unstable',e,P,f,Q,r,M,pts,q))
        pol.append(q)
    ds=[d for d in range(1,M+1) if M%d==0]
    for d in ds:
        if all(pol[r]==pol[(r+d)%M] for r in range(M)):
            return M,d,pol
    raise RuntimeError('no period')

def ray_limit(e,P,f,Q):
    # Leading coefficient is recovered exactly from any three values in one
    # sufficiently stabilized residue class.  Large h is a scouting choice.
    C=abs(Q*e-P*f) or 1
    h=10000
    gs=[1+C*(h+i) for i in range(4)]
    nums=[]; dens=[]
    for g in gs:
        z=L*g*P-e; w=L*g*Q-f
        x=mass(e,g*P,f,g*Q)
        nums.append(x.numerator*(z*w//x.denominator))
        dens.append(z*w)
    # quadratic leading coefficient under unit h steps; g step C
    if nums[3]-3*nums[2]+3*nums[1]-nums[0] != 0:
        raise RuntimeError(('unstable cubic difference',e,P,f,Q,C,nums))
    if dens[3]-3*dens[2]+3*dens[1]-dens[0] != 0:
        raise RuntimeError(('denominator nonquadratic',e,P,f,Q,C,dens))
    an=F(nums[2]-2*nums[1]+nums[0],2)
    ad=F(dens[2]-2*dens[1]+dens[0],2)
    return an/ad

def fracpart(x):
    return x - (x.numerator//x.denominator)

def B3(x):
    x=fracpart(x)
    return x*x*x-F(3,2)*x*x+F(1,2)*x

def bernoulli_limit(e,P,f,Q):
    C=Q*e-P*f
    if C==0: raise RuntimeError(('rank-one case',e,P,f,Q))
    R=e*CELL%L; S=f*CELL%L
    D=Q*R-P*S
    a=F(P,14); b=F(Q,14)
    u=F(D+C,L); v=F(-D,L)
    combo=(B3(u+a-b)+B3(u-a+b)+B3(v+a-b)+B3(v-a+b)
           -B3(u+a+b)-B3(u-a-b)-B3(v+a+b)-B3(v-a-b))
    return F(1,49)+F(28,P*Q*C)*combo

if __name__=='__main__':
    for edge,(i,j) in enumerate(EDGES):
        for o,(e,f) in enumerate(((H[i],H[j]),(H[j],H[i]))):
            best=None
            for Q in range(2,101):
                for P in range((Q+1)//2,Q):
                    if gcd(P,Q)>1 or P+Q<8: continue
                    x=ray_limit(e,P,f,Q)
                    y=bernoulli_limit(e,P,f,Q)
                    if x!=y: raise RuntimeError(('formula mismatch',e,P,f,Q,x,y))
                    if best is None or x<best[0]: best=(x,P,Q)
            print(edge,o,e,f,best,float(best[0]),'minus1/105',best[0]-F(1,105))
