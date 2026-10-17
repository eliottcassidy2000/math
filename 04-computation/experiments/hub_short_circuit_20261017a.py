#!/usr/bin/env python3
"""The incoming 20-step relation joined to the new hub, retaining all 8 ports.
Full arithmetic in Q(zeta55,w,sqrt7,sqrt11), degree 480.
"""
import argparse,hashlib,json,random
from fractions import Fraction as F
from itertools import combinations
from math import lcm
import numpy as np
import sympy as s
from pysat.solvers import Solver
import cubic_crossedge_packets_20261016d as C
import cubic_hub_coupling_20261017a as H

checks=0
def need(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)
O=(0,)*480
def lift(x):return tuple(x)+(0,)*360
I=lift(C.I)
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def sc(x,c):return tuple(c*a for a in x)
def sub(x,y):return add(x,sc(y,-1))
def conj(x):return sum((C.conj(x[120*k:120*k+120]) for k in range(4)),())
def mul(x,y):
    out=[C.O for _ in range(4)]
    xx=[x[120*k:120*k+120] for k in range(4)];yy=[y[120*k:120*k+120] for k in range(4)]
    for r in range(4):
        if not any(xx[r]):continue
        for t in range(4):
            if any(yy[t]):out[r^t]=C.add(out[r^t],C.sc(C.mul(xx[r],yy[t]),(7 if r&t&1 else 1)*(11 if r&t&2 else 1)))
    return sum(out,())
def unit(x):return mul(x,conj(x))==I
def base_hub():
    out=[]
    for h in H.hub():
        p=O
        for k,v in enumerate(h):
            if not v:continue
            # sqrt7, sqrt11, i=j*sqrt11/11.
            block=k&3;coef=C.I
            if k&4:
                coef=C.lift(C.G.b)
                v*=F(1,11)
                if block&2:v*=11
                block^=2
            p=add(p,(0,)*(120*block)+C.sc(coef,v)+(0,)*(120*(3-block)))
        out.append(p)
    need(len(set(out))==9,'literal nine-vertex hub')
    need(sum(unit(sub(x,y)) for x,y in combinations(out,2))==14,'independent full hub edge count')
    return out

def circuit():
    v=C.generator();V={r:lift(C.galois(v,r)) for r in (1,3,4)}
    raw=[(0,4,-1),(1,4,1),(1,19,-1),(1,29,1),(1,39,1),
         (3,0,-1),(3,11,-1),(3,19,-2),(3,22,-1),(3,33,-1),(3,39,-1),
         (4,0,1),(4,4,1),(4,11,1),(4,14,1),(4,22,1),(4,24,-1),(4,33,1),(4,34,1)]
    directions=[]
    for r,k,c in raw:
        d=mul(I if r==0 else V[r],lift(C.lift(C.G.powers[k])))
        need(unit(d),'every imported circuit step unit')
        directions.extend([sc(d,1 if c>0 else -1)]*abs(c))
    total=O
    for d in directions:total=add(total,d)
    need(len(directions)==20 and total==O,'incoming short relation closes in our arithmetic')
    rng=random.Random(20261017);points=[O];orders=[]
    for _ in range(4):
        order=list(range(20));rng.shuffle(order);x=O;orders.append(order)
        for k in order:x=add(x,directions[k]);points.append(x)
        need(x==O,'each prescribed cyclic walk closes')
    return list(dict.fromkeys(points)),V,orders

def graph(points):
    den=lcm(*(F(v).denominator for p in points for v in p));pts=[tuple(int(v*den) for v in p) for p in points]
    filters=[]
    for prime in s.primerange(331,20000):
        if prime%55!=1 or den%prime==0 or pow(7,(prime-1)//2,prime)!=1 or pow(11,(prime-1)//2,prime)!=1:continue
        roots=[t for t in range(prime) if (t**3+2*t*t-4*t-4)%prime==0]
        if not roots:continue
        root=pow(int(s.primitive_root(prime)),(prime-1)//55,prime);w=roots[0]
        r7=next(t for t in range(prime) if t*t%prime==7);r11=next(t for t in range(prime) if t*t%prime==11)
        ev=[]
        for sign in (1,-1):
            powers=[(r7 if k&1 else 1)*(r11 if k&2 else 1)*pow(w,j,prime)*pow(root,sign*i,prime)%prime
                    for k in range(4) for j in range(3) for i in range(40)]
            ev.append([sum(v*c for v,c in zip(p,powers))%prime for p in pts])
        filters.append((prime,np.array(ev,dtype=np.int64).T))
        if len(filters)==3:break
    need(len(filters)==3,'three valid extended-field homomorphisms')
    edges=[];survivors=0;norm_cache={}
    for i in range(len(points)):
        ids=np.arange(i+1,len(points),dtype=np.int64)
        for p,ev in filters:ids=ids[((ev[ids,0]-ev[i,0])*(ev[ids,1]-ev[i,1])-den*den%p)%p==0]
        survivors+=len(ids)
        for k in ids:
            d=sub(pts[i],pts[int(k)])
            key=min(d,sc(d,-1))
            if key not in norm_cache:norm_cache[key]=mul(key,conj(key))==sc(I,den*den)
            if norm_cache[key]:edges.append((i,int(k)))
    for i,k in combinations(range(8),2):need(((i,k) in edges)==unit(sub(points[i],points[k])),'all old ports directly checked')
    return edges,dict(denominator=str(den),filter_primes=[p for p,e in filters],survivors=survivors,
                      exact_displacements=len(norm_cache))

def boundary(n,E):
    port_edges=[(i,j) for i,j in E if i<8 and j<8]
    patterns=[]
    def visit(t):
        k=len(t)
        if k==8:patterns.append(tuple(t));return
        for c in range(min(6,max(t,default=-1)+2)):
            if all(t[i]!=c for i,j in port_edges if j==k):visit(t+[c])
    visit([])
    clauses=[[6*i+c+1 for c in range(6)] for i in range(n)]+[[-6*i-c-1,-6*j-c-1] for i,j in E for c in range(6)]
    clauses += [[-6*i-c-1,-6*i-d-1] for i in range(n) for c,d in combinations(range(6),2)]
    removed=[]
    with Solver(name='g3',bootstrap_with=clauses) as solver:
        for p in patterns:
            if not solver.solve(assumptions=[6*i+c+1 for i,c in enumerate(p)]):removed.append(p);continue
            m=set(solver.get_model());colour=[next(c for c in range(6) if 6*i+c+1 in m) for i in range(n)]
            need(tuple(colour[:8])==p and all(colour[i]!=colour[j] for i,j in E),'each extension witness replayed')
    return dict(source_patterns=len(patterns),removed_patterns=removed)

def run(mode):
    hub=base_hub();seed,V,orders=circuit();addresses={};ports=[lift(C.lift(p)) for p in C.G.ports]
    dirs=[I] if mode=='translated' else [I,V[1],V[3],V[4]]
    points=ports[:];shapes=list(dict.fromkeys(mul(d,h) for d in dirs for h in hub))
    for hi,h in enumerate(shapes):
        for si,q in enumerate(seed):
            p=add(h,q);points.append(p);addresses.setdefault(p,set()).add((hi,si))
    points=list(dict.fromkeys(points));need(points[:8]==ports,'all original source ports fixed')
    E,meta=graph(points);c=C.colour(len(points),E,3);colour_bound=3
    if c is None:c=C.colour(len(points),E,4);colour_bound=4
    need(c is not None and all(c[i]!=c[j] for i,j in E),'positive colour upper bound verified')
    mixed=[];genuine=[];detail=[]
    for i,j in E:
        if points[i] not in addresses or points[j] not in addresses:continue
        if any(h==k or q==t for h,q in addresses[points[i]] for k,t in addresses[points[j]]):continue
        mixed.append([i,j]);d=sub(points[i],points[j])
        if any(any(d[120*k+40:120*k+120]) for k in range(4)):
            genuine.append([i,j])
            detail.append(dict(edge=[i,j],left_addresses=sorted(addresses[points[i]]),
                         right_addresses=sorted(addresses[points[j]]),
                         extension_endpoint=bool(any(points[i][120:]) or any(points[j][120:])),
                         extension_displacement=bool(any(d[120:])),
                         displacement=[str(v) for v in d]))
    b=boundary(len(points),E)
    return dict(mode=mode,vertices=len(points),edges=E,colouring=c,colour_bound=colour_bound,seed_vertices=len(seed),
          boundary=b,mixed_address_edges=mixed,cubic_mixed_address_edges=genuine,mixed_details=detail,geometry=meta,
          orders=orders,coordinate_hash=hashlib.sha256(json.dumps([[str(v) for v in p] for p in points]).encode()).hexdigest(),
          source_ports=[[str(v) for v in p] for p in ports],checks=checks,engine_checks=C.CHECKS)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['translated','rotated'],default='translated');args=ap.parse_args()
    print(json.dumps(run(args.mode),sort_keys=True,separators=(',',':')))
