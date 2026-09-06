#!/usr/bin/env python3
"""Scratch exact audit for all four sorted signed (1,2,2) speed cones."""

from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations
from math import gcd

R = Q(3, 14)
H = 170
CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def fam_u(w):
    a, b, c = w
    out = []
    if 2*a + 2*b == c:
        out.append((1, (2, 2, -1)))
    if 2*a + c == 2*b:
        out.append((2, (2, -2, 1)))
    if a + 2*b == 2*c:
        out.append((3, (1, 2, -2)))
    if 2*a + b == 2*c:
        out.append((4, (2, 1, -2)))
    return out


def roofs(w):
    return tuple((3*(sum(w)-w[i])-1)//14 for i in range(3))


def line_row(w, u):
    a, b, c = w
    bounds = roofs(w)
    K = min(bounds[i]//abs(u[i]) for i in range(3))
    addresses = tuple(k for k in range(1, K+1) if k % 3)
    E = [Q(0)]*3
    mass = Q(0)
    for k in addresses:
        terms = []
        for i in range(3):
            j, ell = [z for z in range(3) if z != i]
            term = min(Q(3, 7*c),
                       Q(3*(w[j]+w[ell])-14*abs(u[i])*k,
                         14*w[j]*w[ell]))
            need(term > 0, ('positive line term', w, u, k, i, term))
            terms.append(term)
            E[i] += 2*term
        mass += 2*min(terms)
    carriers = {(s*k*u[0], s*k*u[1], s*k*u[2])
                for k in addresses for s in (-1, 1)}
    return tuple(E), mass, carriers, addresses


def raw(w):
    a, b, c = w
    bounds = roofs(w)
    live = set()
    d = gcd(b, c)
    modulus = c//d
    inv = pow(b//d, -1, modulus) if modulus > 1 else 0
    for x in range(-bounds[0], bounds[0]+1):
        if x % 3 == 0 or a*x % d:
            continue
        residue = (-a*x//d*inv) % modulus if modulus > 1 else 0
        low = -bounds[1] + (residue+bounds[1]) % modulus
        for y in range(low, bounds[1]+1, modulus):
            if y % 3 == 0:
                continue
            z = -(a*x+b*y)//c
            if z % 3 and abs(z) <= bounds[2]:
                live.add((x,y,z))
    E = [Q(0)]*3
    mass = Q(0)
    for C in live:
        terms=[]
        for i in range(3):
            j,ell=[z for z in range(3) if z!=i]
            term=min(Q(3,7*c),Q(3*(w[j]+w[ell])-14*abs(C[i]),14*w[j]*w[ell]))
            need(term>0, ('positive raw term', w, C, i, term))
            E[i]+=term
            terms.append(term)
        mass+=min(terms)
    return tuple(E),mass,live


def fvals(w, u, x):
    scale = w[2]
    W = tuple(Q(v, scale) for v in w)
    vals=[]
    for i in range(3):
        j,k=[z for z in range(3) if z!=i]
        vals.append(min(2*R,(R*(W[j]+W[k])-abs(u[i])*x)/(W[j]*W[k])))
    return tuple(vals)


def integral_piecewise(w, u, physical=False):
    W=tuple(Q(v,w[2]) for v in w)
    B=min(R*(W[j]+W[k])/abs(u[i])
          for i in range(3) for j,k in [tuple(z for z in range(3) if z!=i)])
    # Every profile is min(2R, alpha-beta*x). Build all exact breakpoints,
    # then select by an interior point and trapezoid each affine segment.
    lines=[]
    points={Q(0),B}
    for i in range(3):
        j,k=[z for z in range(3) if z!=i]
        alpha=R*(W[j]+W[k])/(W[j]*W[k])
        beta=Q(abs(u[i]),1)/(W[j]*W[k])
        lines.append((alpha,beta))
        cap=(alpha-2*R)/beta
        if 0<cap<B:
            points.add(cap)
    for i in range(3):
        for j in range(i):
            ai,bi=lines[i]; aj,bj=lines[j]
            if bi!=bj:
                x=(ai-aj)/(bi-bj)
                if 0<x<B:
                    points.add(x)
    points=sorted(points)
    values=[Q(0)]*3
    pvalue=Q(0)
    for lo,hi in zip(points,points[1:]):
        mid=(lo+hi)/2
        def v(i,x):
            a,b=lines[i]
            return min(2*R,a-b*x)
        if physical:
            owner=min(range(3),key=lambda i:v(i,mid))
            pvalue+=(hi-lo)*(v(owner,lo)+v(owner,hi))/2
        else:
            for i in range(3):
                values[i]+=(hi-lo)*(v(i,lo)+v(i,hi))/2
    return pvalue if physical else tuple(values)


def update(rec, key, value, payload):
    if key not in rec or value > rec[key][0]:
        rec[key]=(value,[payload])
    elif value == rec[key][0]:
        rec[key][1].append(payload)


rows=[]
for w in combinations([x for x in range(1,H+1) if x%3],3):
    if gcd(gcd(*w[:2]),w[2])!=1:
        continue
    ff=fam_u(w)
    if ff:
        need(len(ff)==1, ('unique signed cone', w, ff))
        rows.append((w,*ff[0]))

leaders={}
counts=defaultdict(int)
formula_forms=defaultdict(set)
for w,f,u in rows:
    counts[f]+=1
    E,mass,C,K=line_row(w,u)
    Er,mr,Cr=raw(w)
    need((E,mass,C)==(Er,mr,Cr),('line versus raw', w,(E,mass,C),(Er,mr,Cr)))
    for k in K:
        vals=fvals(w,u,Q(k,w[2]))
        terms=tuple(Q(w[2],1)*Q(3*(w[j]+w[ell])-14*abs(u[i])*k,
                                     14*w[j]*w[ell])
                    for i in range(3) for j,ell in [tuple(z for z in range(3) if z!=i)])
        terms=tuple(min(2*R,x) for x in terms)
        need(vals==terms,('normalized exact term', w,k,vals,terms))
    ints=integral_piecewise(w,u)
    pint=integral_piecewise(w,u,True)
    t=Q(w[0],w[2]) if f!=4 else Q(w[1],w[2])
    r2=R*R
    if f==1:
        claimed=(r2*(Q(7,8)+t/4), r2*(1-t/4), r2*(1-t+2*t*t))
    elif f==2:
        third = (r2*(1+t) if t <= Q(1,4) else
                 r2*(-16*t**3+16*t**2+11*t+4)/(4*(2*t+1)))
        claimed=(r2*Q(7+16*t,8*(1+2*t)), r2, third)
    elif f==3:
        claimed=(r2*(Q(7,8)+9*t/16), r2*(1-t/16), r2*(1-t/2+t*t/2))
    else:
        claimed=(r2*(1-t/16), r2*(Q(7,8)+9*t/16), r2*(1-t/2+t*t/2))
    for i,x in enumerate(claimed):
        if x is not None:
            need(ints[i]==x,('projection integral', w,i,ints[i],x))
    need(pint==r2*Q(7,8),('universal norm-five physical integral', f,w,pint))
    update(leaders,('N',f),min(E),(w,E,mass,K))
    update(leaders,('P',f),mass,(w,E,mass,K))
    parity=sum(v%2==0 for v in w)
    update(leaders,('N',f,parity),min(E),(w,E,mass,K))
    update(leaders,('P',f,parity),mass,(w,E,mass,K))

print('H',H,'rows',len(rows),'counts',dict(counts))
for key in sorted(leaders):
    print(key,'leader',leaders[key])
print('tail constants')
for f,L in [(1,Q(87,1568)),(2,Q(45,784)),(3,Q(363,6272)),(4,Q(363,6272))]:
    print('N',f,L)
print('all_family_physical_bulk',Q(3,56))
print('explicit_checks',CHECKS)
print('PASS')
