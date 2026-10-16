#!/usr/bin/env python3
"""Exact cubic-extension unit graphs; full induced edges and fixed old ports.

The integer power basis is zeta55^i w^j, i<40,j<3,
w^3+2w^2-4w-4=0, w=2theta. No numerical geometry is used.
"""
import os
os.environ['SYMPY_GROUND_TYPES'] = 'python'
import argparse
import hashlib
import json
from fractions import Fraction as F
from itertools import combinations
from math import lcm
import numpy as np
import sympy as s
from pysat.solvers import Solver
import charge_ports_eleven_20261011 as G

CHECKS = 0
def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)

O = (0,)*120
def lift(p): return tuple(p)+(0,)*80
I = lift(G.I)
W = (0,)*40+G.I+(0,)*40
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sc(a,c): return tuple(c*x for x in a)
def sub(a,b): return add(a,sc(b,-1))
def mul(a,b):
    c=[G.O for _ in range(5)]
    aa=[a[40*i:40*i+40] for i in range(3)]
    bb=[b[40*i:40*i+40] for i in range(3)]
    for i in range(3):
        if not any(aa[i]): continue
        for j in range(3):
            if any(bb[j]): c[i+j]=G.add(c[i+j],G.mul(aa[i],bb[j]))
    for d in [4,3]:
        for k,v in enumerate([4,4,-2]):
            c[d-3+k]=G.add(c[d-3+k],G.sc(c[d],v))
    return sum(c[:3],())
def conj(a): return sum((G.conj(a[40*i:40*i+40]) for i in range(3)),())
def power(a,k):
    if k<0:return power(conj(a),-k)
    out=I
    while k:
        if k&1:out=mul(out,a)
        a=mul(a,a);k//=2
    return out
def unit(a):return mul(a,conj(a))==I

def galois(a,h):
    e=next(e for e in range(55) if e%5==1 and e%11==h)
    blocks=[]
    for k in range(3):
        out=G.O
        for i,c in enumerate(a[40*k:40*k+40]):
            if c:out=G.add(out,G.sc(G.powers[e*i%55],c))
        blocks.extend(out)
    return tuple(blocks)

def generator():
    # Coefficients in increasing theta power: (rational part, j part).
    cc=[[(440,84),(-47,-51),(-42,20)],
        [(-263,272),(655,55),(-136,-82)],
        [(-487,-135),(791,137),(838,188)],
        [(927,219),(-838,-188),(-880,-168)]]
    v=O
    for h,c in enumerate(cc):
        for i,(r,j) in enumerate(c):
            b=G.mul(G.add(G.sc(G.I,r),G.sc(G.b,j)),G.powers[5*h])
            term=(0,)*(40*i)+G.sc(b,F(1,1541*2**i))+(0,)*(40*(2-i))
            v=add(v,term)
    theta=sc(W,F(1,2));xi=lift(G.powers[5])
    den=add(mul(theta,sub(I,conj(xi))),power(xi,-2))
    num=add(mul(theta,sub(I,xi)),power(xi,2))
    check(mul(v,den)==num,'independent tensor division identity')
    check(unit(v),'cubic direction norm one')
    check(any(v[40:]),'direction really leaves the original field')
    check(add(add(power(W,3),sc(power(W,2),2)),add(sc(W,-4),sc(I,-4)))==O,
          'integer cubic relation')
    return v

def graph(points):
    den=lcm(*(F(x).denominator for p in points for x in p))
    pts=[tuple(int(x*den) for x in p) for p in points]
    filters=[]
    for prime in [331,661,991]:
        check(den%prime!=0,'filter prime avoids global denominator')
        r=pow(int(s.primitive_root(prime)),(prime-1)//55,prime)
        roots=[t for t in range(prime) if (t**3+2*t*t-4*t-4)%prime==0]
        check(roots and int(G.poly.eval(r))%prime==0,'two defining relations vanish')
        t=roots[0]
        rp=[pow(t,j,prime)*pow(r,i,prime)%prime for j in range(3) for i in range(40)]
        rc=[pow(t,j,prime)*pow(r,-i,prime)%prime for j in range(3) for i in range(40)]
        ev=np.array([(sum(x*y for x,y in zip(p,rp))%prime,
                      sum(x*y for x,y in zip(p,rc))%prime) for p in pts],dtype=np.int64)
        filters.append((prime,t,ev))
    edges=[];surv=0
    for i in range(len(points)):
        ids=np.arange(i+1,len(points),dtype=np.int64)
        for prime,t,ev in filters:
            ok=((ev[ids,0]-ev[i,0])*(ev[ids,1]-ev[i,1])-den*den%prime)%prime==0
            ids=ids[ok]
        surv+=len(ids)
        for j in ids:
            d=sub(pts[i],pts[int(j)])
            if mul(d,conj(d))==sc(I,den*den):edges.append((i,int(j)))
    # Independent direct norm controls on every old port pair.
    for i,j in combinations(range(8),2):
        check(((i,j) in edges)==unit(sub(points[i],points[j])), 'all port pairs directly checked')
    return edges,dict(denominator=str(den),survivors=surv,
                      filters=[dict(prime=p,w_root=t) for p,t,ev in filters])

def colour(n,edges,q):
    clauses=[[q*i+c+1 for c in range(q)] for i in range(n)]
    clauses += [[-q*i-c-1,-q*j-c-1] for i,j in edges for c in range(q)]
    with Solver(name='g3',bootstrap_with=clauses) as solver:
        if not solver.solve():return None
        model=set(solver.get_model())
    c=[next(k for k in range(q) if q*i+k+1 in model) for i in range(n)]
    check(all(c[i]!=c[j] for i,j in edges),'all actual edges properly coloured')
    return c

def degeneracy(n,edges):
    adj=[set() for _ in range(n)]
    for i,j in edges:adj[i].add(j);adj[j].add(i)
    left=set(range(n));order=[];peak=0
    while left:
        v=min(left,key=lambda k:(len(adj[k]&left),k))
        peak=max(peak,len(adj[v]&left));order.append(v);left.remove(v)
    return peak,order

def circuit(v,xi):
    z,t=s.symbols('z t')
    P=1-z-t+t/z;Q=t/z**2-z**2
    polynomial=s.Poly(s.expand(z**6*(2*Q**3+2*Q**2*P-2*Q*P**2-P**3)),t,z)
    terms={}
    for (i,j),c in polynomial.terms():terms[i,j%11]=terms.get((i,j%11),0)+int(c)
    terms={a:c for a,c in terms.items() if c}
    value=O;positive=[];negative=[];directions=[]
    for (i,j),c in sorted(terms.items()):
        d=mul(power(v,i),power(xi,j));directions.append(d)
        check(unit(d),'every cubic circuit step has unit norm')
        value=add(value,sc(d,c))
        (positive if c>0 else negative).extend([(i,j,d)]*abs(c))
    check(value==O and len(positive)==len(negative)==32,'two literal 32-step walks meet')
    check(s.rem(polynomial.as_expr().subs(t,z*z),s.cyclotomic_poly(11,z),z)!=0,
          'theta=1 is norm-one but does not close this circuit')
    points=[];paths=[];endpoints=[]
    orderings=[lambda x:(x[0],x[1]),lambda x:(x[1],x[0]),
               lambda x:(-x[0],x[1]),lambda x:(x[1]%3,x[0],x[1])]
    for key in orderings:
        for side in [positive,negative]:
            p=O;walk=[p]
            for i,j,d in sorted(side,key=key):p=add(p,d);walk.append(p)
            check(len(walk)==33,'literal circuit length')
            paths.append([[i,j] for i,j,d in sorted(side,key=key)])
            points+=walk;endpoints.append(p)
    check(len(set(endpoints))==1,'all eight prescribed paths share their endpoint')
    return list(dict.fromkeys(points)),directions,dict(
        relation=[[i,j,c] for (i,j),c in sorted(terms.items())],
        positive_length=32,negative_length=32,
        common_endpoint=[str(x) for x in endpoints[0]],
        path_steps=paths)

def algebra(v):
    x=s.symbols('x');f=x**3+2*x*x-4*x-4
    check(s.discriminant(f,x)==592,'cubic discriminant')
    check(all(f.subs(x,a)*f.subs(x,b)<0 for a,b in [(-4,-2),(-1,0),(1,2)]),
          'three distinct real roots in disjoint intervals')
    check(f.subs(x,4)%19==0 and s.diff(f,x).subs(x,4)%19==3,
          'simple rational residue root at nineteen')
    check(pow(8,9,19)==18,'minus eleven is a nonsquare at nineteen')
    table=[4,3,2,4,3,2,1,0,4,1,0,4,3,2,1,0,2,1,0]
    directions=[(a,b) for a in range(19) for b in range(19) if (a*a+b*b)%19==1]
    for a,b in directions:
        for c in range(19):check(table[c]!=table[(c+a+b)%19],'recovered five-colour residue table')
    def trace(u):
        return G.add(G.sc(u[:40],3),G.add(G.sc(u[40:80],-2),G.sc(u[80:],12)))
    mean=G.sc(trace(v),F(1,3));dev=sub(v,lift(mean))
    defect=G.add(G.I,G.sc(G.mul(mean,G.conj(mean)),-1))
    check(trace(mul(dev,conj(dev)))==G.sc(defect,3),'exact trace variance identity')
    check(defect!=G.O,'trace does not preserve the new unit edge')
    return dict(cubic=str(f),discriminant=592,residue_prime=19,residue_root=4,
                residue_colouring=table,mean_direction=[str(x) for x in mean],
                squared_length_defect=[str(x) for x in defect])

def run(mode):
    v=generator();xi=lift(G.powers[5]);alpha=lift(G.powers[11])
    extra={};alg=algebra(v)
    if mode=='powers':
        positive=[mul(power(v,k),power(xi,h)) for k in [-1,0,1] for h in range(11)]
    elif mode=='mixed':
        eta=lift(G.eta);gamma=lift(G.sc(G.add(G.sc(G.I,5),G.b),F(1,6)))
        positive=[mul(power(v,k),power(xi,h)) for k in [0,1] for h in range(11)]
        positive += [mul(power(v,k),power(alpha,h)) for k in [0,1] for h in range(5)]
        positive += [eta,conj(eta),gamma,conj(gamma)]
    elif mode=='conjugates':
        positive=[mul(galois(v,h),power(xi,r)) for h in [1,3,4,5,9] for r in range(11)]
        positive += [power(xi,r) for r in range(11)]+[power(alpha,r) for r in range(5)]
    elif mode in ['circuit','circuit_conjugates']:
        seed,positive,extra=circuit(v,xi)
        if mode=='circuit_conjugates':
            positive=[galois(d,h) for h in [1,3,4,5,9] for d in positive]
    else:raise ValueError(mode)
    positive=list(dict.fromkeys(positive))
    check(all(unit(d) for d in positive),'all declared directions have unit norm')
    directions=positive+[sc(p,-1) for p in positive]
    addresses={}
    if mode=='circuit':
        points=list(dict.fromkeys([lift(p) for p in G.ports]+
                                 [add(lift(p),x) for p in G.ports for x in seed]))
        for i,p in enumerate(G.ports):
            for j,x in enumerate(seed):addresses.setdefault(add(lift(p),x),set()).add((i,j))
    elif mode=='circuit_conjugates':
        points=list(dict.fromkeys([lift(p) for p in G.ports]+
                                 [galois(x,h) for h in [1,3,4,5,9] for x in seed]))
        for i,h in enumerate([1,3,4,5,9]):
            for j,x in enumerate(seed):addresses.setdefault(galois(x,h),set()).add((i,j))
    else:
        points=list(dict.fromkeys([lift(p) for p in G.ports]+directions+
                                 [add(a,b) for i,a in enumerate(directions) for b in directions[i:]]))
    check(points[:8]==[lift(p) for p in G.ports],'all original ports literally fixed')
    edges,meta=graph(points)
    c=colour(len(points),edges,3)
    d,order=degeneracy(len(points),edges)
    dirs=set(directions)
    cross=[(i,j) for i,j in edges if sub(points[j],points[i]) not in dirs]
    if mode=='circuit':
        mixed=[(i,j) for i,j in edges if not any(a==c or b==d
               for a,b in addresses[points[i]] for c,d in addresses[points[j]])]
        extra.update(seed_vertices=len(seed),identifications=8*len(seed)-len(points),
                     mixed_edges=mixed)
    elif mode=='circuit_conjugates':
        mixed=[(i,j) for i,j in edges if points[i] in addresses and points[j] in addresses and
               not any(a==c for a,b in addresses[points[i]] for c,d in addresses[points[j]])]
        extra.update(seed_vertices=len(seed),cross_copy_edges=mixed,
                     cubic_cross_copy_edges=[(i,j) for i,j in mixed
                                             if any(points[i][40:]) or any(points[j][40:])])
    coordinates=[[str(x) for x in p] for p in points]
    coordinate_hash=hashlib.sha256(json.dumps(coordinates,separators=(',',':')).encode()).hexdigest()
    return dict(mode=mode,universe=('eight port translates of eight cubic-relation unit walks' if mode=='circuit'
                else 'five K-fixing Galois copies of the cubic-relation walk bank, plus fixed ports' if mode=='circuit_conjugates'
                else 'old eight ports and all sums of at most two symmetric directions'),
                vertices=len(points),edges=edges,three_colouring=c,degeneracy=d,peel_order=order,
                additional_edges=cross,directions=[[str(x) for x in p] for p in positive],
                coordinate_sha256=coordinate_hash,geometry=meta,
                source_port_coordinates=coordinates[:8],algebra=alg,
                source_ports=list(range(8)),circuit=extra,checks=CHECKS)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['powers','mixed','conjugates','circuit','circuit_conjugates'],default='powers')
    a=ap.parse_args();print(json.dumps(run(a.mode),indent=2,sort_keys=True))
