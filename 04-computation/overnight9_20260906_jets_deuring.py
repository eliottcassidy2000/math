#!/usr/bin/env python3
"""Exact one-digit Deuring precision loss, with independent controls.

No producer/repository imports. Integer polynomial ODE/Wronskian checks,
all lifts above small residue roots, and literal reciprocal/Smith checks
are controls for the analytic all-prime proof in the companion report.
"""
from fractions import Fraction as Q
from math import comb
from itertools import combinations
import sys
from sympy import Matrix,ZZ
from sympy.matrices.normalforms import smith_normal_form

sys.stdout.reconfigure(newline='\n')
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)


def trim(a):
    while len(a)>1 and not a[-1]:a.pop()
    return a


def add(a,b):
    c=[0]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return trim(c)


def scale(a,c):return trim([c*x for x in a])


def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)


def deriv(a):return trim([i*a[i] for i in range(1,len(a))] or [0])


def ev(a,x,mod=None):
    out=0
    for c in reversed(a):
        out=out*x+c
        if mod:out%=mod
    return out


def reflection(a):
    return trim([sum(a[j]*comb(j,i)*(-1)**i for j in range(i,len(a))) for i in range(len(a))])


def ode0(a,k):
    return add(add(mul([0,1,-1],deriv(deriv(a))),mul([1,-2],deriv(a))),scale(a,k*(k+1)))


def vp(x,p):
    x=Q(x)
    if not x:return float('inf')
    a,b=abs(x.numerator),x.denominator;v=0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v


def literal_jets(nodes,m):
    jets=[]
    for i,x in enumerate(nodes):
        b=[Q(1)]+[Q(0)]*(m-1)
        for j,y in enumerate(nodes):
            if i==j:continue
            f=[Q((-1)**r*comb(m+r-1,r),(x-y)**(m+r)) for r in range(m)]
            b=[sum(b[s]*f[r-s] for s in range(r+1)) for r in range(m)]
        jets.append(b)
    return jets


def smith(nodes,m,p):
    M=len(nodes)*m
    J=Matrix([[comb(j,r)*x**(j-r) if j>=r else 0 for j in range(M)] for x in nodes for r in range(m)])
    D=smith_normal_form(J,domain=ZZ)
    return tuple(vp(D[i,i],p) for i in range(M))


def polynomial_family(p):
    k=(p-1)//2
    P=[comb(k+j,j)*comb(2*k-j,k-j) for j in range(k+1)]
    R=scale(reflection(P),(-1)**k)
    difference=add(P,scale(R,-1))
    need(all(c%p==0 for c in difference),('integral divided companion',p))
    G=trim([c//p for c in difference])
    need(len(G)<=k,('strict companion degree',p))
    need(add(ode0(P,k),scale(deriv(P),-p))==[0],('exact P hypergeometric ODE',p))
    need(add(ode0(R,k),scale(deriv(R),p))==[0],('exact reflected ODE',p))
    need(add(ode0(G,k),scale(add(deriv(P),deriv(R)),-1))==[0],('exact divided inhomogeneous ODE',p))
    W=add(mul(P,deriv(G)),scale(mul(deriv(P),G),-1))
    F=add(add(mul([0,1,-1],W),scale(mul(P,P),-1)),[P[0]**2])
    need(len(F)<=p,('Wronskian degree strictly below characteristic',p,len(F)-1))
    need(all(c%p==0 for c in F),('Wronskian constant identity',p))
    need(P[0]%p==(-1)**k%p and ev(P,1,p)==1,'both endpoint constants are units')
    H=[comb(k,j)**2 for j in range(k+1)]
    need(all((P[j]-(-1)**k*H[j])%p==0 for j in range(k+1)),('Deuring reduction',p))
    third=[0]*(k+1)
    for j,c in enumerate(P):
        for r in range(j+1):third[k-j+r]+=c*comb(j,r)*(-1)**(j-r)
    need(all((third[j]-(-1)**k*P[j])%p==0 for j in range(k+1)),('third endpoint S3 reduction',p))
    # A=[X^2k]F_a(X)^k=(-1)^k H; B=[X^(2k-1)]F_a(X)^k.
    AA=scale(H,(-1)**k)
    B=[0]*(k+1)
    for r in range(k):
        s=k-1-r
        B[k-s]+=comb(k,r)*comb(k,s)*(-1)**(2*k-r-s)
    rhs=scale(mul([0,1],add(scale(AA,k),mul([1,-1],deriv(AA)))),-1)
    need(add(scale(B,k+1),scale(rhs,-1))==[0],('exact adjacent coefficient derivative identity',p))
    return P,R,G,H,third,B


def run():
    all_lifts=0;root_count=0;matrix_count=0;inverse_count=0
    for p in (3,5,7,11,13,19,23,31,43,47,59,67,71,79,83,103):
        k=(p-1)//2;m=k+1
        P,R,G,H,third,B=polynomial_family(p)
        roots=[a for a in range(2,p) if ev(H,a,p)==0]
        root_count+=len(roots)
        Pd=deriv(P)
        for a in roots:
            need(a*(1-a)*ev(Pd,a,p)*ev(G,a,p)%p==1,('root derivative/companion unit identity',p,a))
            need(ev(B,a,p)!=0,('adjacent coefficient unit at supersingular root',p,a))
            j=-(ev(P,a,p*p)//p)*pow(ev(Pd,a,p),-1,p)%p
            lift=a+p*j
            need(ev(P,lift,p*p)==0 and ev(R,lift,p*p)!=0,('unique first numerator Hensel lift misses companion',p,a,j))
            if p<=31:
                for digit in range(p):
                    b=a+p*digit
                    top=[ev(poly,b,p*p) for poly in (P,R,third)]
                    need(all(v%p==0 for v in top) and any(v for v in top),('complete p2 root-fiber cap',p,a,digit))
                    all_lifts+=1
        if p<=31:
            ordinary=next((a for a in range(2,p) if a not in roots),None)
            examples=sorted(set(([roots[0]] if roots else [])+([ordinary] if ordinary else [])
                                +([roots[0]+p*p] if roots else [])))
            for a in examples:
                normalized=literal_jets((0,1,a),m)
                sigma=int(ev(H,a,p)==0)
                need(min(vp(row[k],p) for row in normalized)==sigma,('actual rational top coefficients',p,a))
                if sigma:
                    need(all(vp(row[k-1],p)==0 for row in normalized),('all three actual adjacent jets are units',p,a))
                for e in (0,1,2,4):
                    nodes=(0,p**e,a*p**e)
                    actual=max(-vp(c,p) for row in literal_jets(nodes,m) for c in row if c)
                    predicted=0 if e==0 else (3*m-1)*e-sigma
                    need(actual==predicted,('actual all-depth largest-loss control',p,a,e))
                    inverse_count+=1
                    if p<=11 and e in (0,1,2) and a<p:
                        spectrum=smith(nodes,m,p);matrix_count+=1
                        need(spectrum[-1]==predicted,('full observer largest exponent',p,a,e))
                        need(sum(spectrum)==3*m*m*e,('confluent determinant normalization',p,a,e))
                        need(sum(spectrum[:-1])==3*m*m*e-predicted,('penultimate determinantal ideal',p,a,e))
        print('PRIME',p,'m',m,'Deuring_roots',len(roots),'all_companion_obstructions_unit',True)
    need(all_lifts==3+21+33+57+207+279,('explicit complete p2 lift universe',all_lifts))
    print('ROOT_CONTROLS',root_count,'rootclasses; completep2fibers',all_lifts)
    print('EXACT_OBSERVER_CONTROLS inverse',inverse_count,'fullSmith',matrix_count)
    print('EXACT_LOSS e>=1: (3m-1)e-[H_p(a)=0]; e0:0; all supersingular next-order jets unit')
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)


if __name__=='__main__':run()
