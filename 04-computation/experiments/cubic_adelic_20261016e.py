#!/usr/bin/env python3
"""Exact mixed-prime cubic direction, joint circuit, and induced unit graphs."""
import argparse
import hashlib
import json
from fractions import Fraction as F
import cubic_crossedge_packets_20261016d as C
import sympy as s

checks=0
def need(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def new_direction():
    j=C.lift(C.G.b);w=C.W
    u=C.sc(C.add(C.add(C.sc(w,11),C.sc(C.I,-6)),
        C.mul(C.add(C.sub(C.power(w,2),w),C.sc(C.I,-12)),j)),F(1,38))
    a=C.sub(C.add(w,C.sc(C.I,3)),j)
    need(C.mul(u,C.conj(a))==a,'Hilbert 90 identity with literal numerator')
    need(C.unit(u),'nineteen direction is unit')
    need(any(u[40:]),'direction leaves L')
    D=C.add(C.add(C.power(w,2),C.sc(w,6)),C.sc(C.I,20))
    need(C.mul(D,C.sub(C.sc(C.I,4),w))==C.sc(C.I,76),'inverse denominator identity')
    x=s.symbols('x');f=x**3+2*x*x-4*x-4
    need(all(int(c)%19==0 for c in s.Poly(f-(x-4)*(x*x+6*x+1),x).all_coeffs()),
         'one plus two factorization modulo nineteen')
    need(pow(8,9,19)==18,'quadratic residue extension at nineteen')
    need(s.resultant(f,x*x+6*x+20,x)==76**2,'denominator norm')
    need(C.sc(C.add(C.sc(C.add(u,C.conj(u)),19),C.sc(C.I,6)),F(1,11))==w,
         'new real cubic coordinate is in the unit-path ring')
    return u

def connector(u):
    alpha=C.lift(C.G.powers[11]);short=[]
    for r in range(1,5):
        a=C.power(alpha,r)
        beta=C.sc(C.add(C.add(C.sc(C.power(a,3),-24),C.sc(C.power(a,2),-32)),
                  C.add(C.sc(a,-75),C.sc(C.I,72))),F(1,121))
        need(C.unit(beta),'inherited denominator-eleven connector direction')
        short.extend([beta]*11)
    short.extend([C.sc(C.I,-1)]*38)
    need(len(short)==82 and sum_vectors(short)==C.sc(C.I,F(1,11)),
         'literal eighty-two-edge connector to one eleventh')
    directions=[C.mul(a,b) for a in [u,C.conj(u),C.sc(C.I,-1),C.sc(C.I,-1)] for b in short]
    need(len(directions)==328 and all(C.unit(d) for d in directions),'literal connector unit steps')
    end=C.sc(C.sub(C.W,C.sc(C.I,4)),F(1,19))
    need(sum_vectors(directions)==end,'allowed cubic displacement with rational denominator nineteen')
    walks=[];points=[]
    # Two fully specified exact orders; these are paths, not assertions of induced paths.
    for order in [list(range(328)),[82*i+j for j in range(82) for i in range(4)]]:
        x=C.O;walk=[x]
        for k in order:x=C.add(x,directions[k]);walk.append(x)
        need(x==end,'both literal connector orders reach the marked endpoint')
        walks.append(order);points+=walk
    return list(dict.fromkeys(points)),dict(length=328,endpoint=[str(x) for x in end],
             direction_orderings=walks,directions=[[str(x) for x in d] for d in directions])

def sum_vectors(xs):
    a=C.O
    for x in xs:a=C.add(a,x)
    return a

def joint_circuit(u,v,xi):
    U,V,z,j=s.symbols('U V z j')
    A=-3*(1-U)+j*(1+U)
    P=s.expand(V*(A*(z*z-z)+2*(1-U))-(A*(z*z-z**3)+2*(1-U)*z**4))
    gauss=sum(int(s.legendre_symbol(h,11))*z**h for h in range(1,11))
    poly=s.Poly(s.expand(P.subs(j,gauss)),U,V,z)
    raw={}
    for (a,b,h),val in poly.terms():raw[a,b,h%11]=raw.get((a,b,h%11),0)+int(val)
    # Each Phi_11 block may be shifted by a constant without changing its value.
    terms={}
    shifts=[]
    for a in range(2):
        for b in range(2):
            row=[raw.get((a,b,h),0) for h in range(11)]
            med=sorted(row)[5];shifts.append([a,b,med])
            for h,c in enumerate(row):
                if c-med:terms[a,b,h]=c-med
    directions={key:C.mul(C.mul(C.power(u,key[0]),C.power(v,key[1])),C.power(xi,key[2])) for key in terms}
    total=C.O
    for key,c in terms.items():
        need(C.unit(directions[key]),'each joint circuit step has unit norm')
        total=C.add(total,C.sc(directions[key],c))
    need(total==C.O,'the joint bilinear circuit closes exactly')
    hostile=C.O
    for (a,b,h),c in terms.items():
        hostile=C.add(hostile,C.sc(C.mul(C.power(u,a),C.power(xi,2*b+h)),c))
    need(hostile!=C.O,'independent theta=1 substitute for v does not close')
    words=[[],[]]
    for key,c in sorted(terms.items()):words[c<0].extend([key]*abs(c))
    points=[];paths=[];ends=[]
    orders=[lambda k:k,lambda k:(k[2],k[0],k[1]),lambda k:(k[1],-k[2],k[0]),lambda k:(k[2]%3,k)]
    for order in orders:
        for word in words:
            word=sorted(word,key=order);point=C.O;points.append(point)
            for key in word:point=C.add(point,directions[key]);points.append(point)
            paths.append([list(k) for k in word]);ends.append(point)
    need(len(set(ends))==1,'all circuit paths have the same endpoint')
    return list(dict.fromkeys(points)),dict(terms=[[*k,c] for k,c in sorted(terms.items())],
        block_shifts=shifts,lengths=list(map(len,words)),paths=paths,
        endpoint=[str(x) for x in ends[0]])

def run(mode):
    u=new_direction();v=C.generator();xi=C.lift(C.G.powers[5]);alpha=C.lift(C.G.powers[11])
    seed,circuit=joint_circuit(u,v,xi);memberships={};addresses={};extra_metadata={}
    if mode=='bank':
        positive=[C.mul(C.mul(C.power(u,a),C.power(v,b)),C.power(xi,h))
                  for a,b in [(0,0),(1,0),(0,1),(1,1)] for h in range(11)]
        positive += [C.power(alpha,h) for h in range(5)]
        positive=list(dict.fromkeys(positive));dirs=positive+[C.sc(x,-1) for x in positive]
        points=list(dict.fromkeys([C.lift(p) for p in C.G.ports]+dirs+
                     [C.add(a,b) for i,a in enumerate(dirs) for b in dirs[i:]]))
    elif mode=='circuit':
        points=list(dict.fromkeys([C.lift(p) for p in C.G.ports]+
                  [C.add(C.lift(p),x) for p in C.G.ports for x in seed]))
        dirs=[]
    elif mode=='two_circuits':
        old,dirs,unused=C.circuit(v,xi)
        copies=[seed,old,[C.mul(u,p) for p in old]]
        for tag,copy in enumerate(copies):
            for p in copy:memberships.setdefault(p,set()).add(tag)
        points=list(dict.fromkeys([C.lift(p) for p in C.G.ports]+sum(copies,[])))
    elif mode=='connector':
        conn,extra_metadata=connector(u)
        copies=[[C.add(C.lift(p),x) for x in conn] for p in C.G.ports]
        for tag,copy in enumerate(copies):
            for i,p in enumerate(copy):
                memberships.setdefault(p,set()).add(tag)
                addresses.setdefault(p,set()).add((tag,i))
        points=list(dict.fromkeys([C.lift(p) for p in C.G.ports]+sum(copies,[])));dirs=[]
    else:raise ValueError(mode)
    need(points[:8]==[C.lift(p) for p in C.G.ports],'all eight old ports fixed')
    edges,meta=C.graph(points);colour=C.colour(len(points),edges,3)
    need(colour is not None,'positive three-colour witness exists in declared finite universe')
    need(all(colour[i]!=colour[j] for i,j in edges),'independent colour replay')
    degree,order=C.degeneracy(len(points),edges)
    bank=set(dirs)
    extra=[(i,j) for i,j in edges if C.sub(points[j],points[i]) not in bank and C.sub(points[i],points[j]) not in bank] if dirs else []
    cubic=[(i,j) for i,j in extra if any(points[i][40:]) or any(points[j][40:])]
    mixed=[(i,j) for i,j in edges if points[i] in memberships and points[j] in memberships
           and not memberships[points[i]]&memberships[points[j]]]
    new_mixed=[(i,j) for i,j in mixed if any(points[i][40:]) or any(points[j][40:])]
    off_axes=[(i,j) for i,j in edges if points[i] in addresses and points[j] in addresses
        and not any(a==c or b==d for a,b in addresses[points[i]] for c,d in addresses[points[j]])]
    new_off_axes=[(i,j) for i,j in off_axes if any(C.sub(points[i],points[j])[40:])]
    return dict(status='FINITE-EXACT; full induced finite graph, not the full field graph',mode=mode,
        vertices=len(points),edges=edges,three_colouring=colour,degeneracy=degree,
        degeneracy_order=order,filters=meta,nonbank_edges=extra,cubic_nonbank_edges=cubic,
        direction=[str(x) for x in u],joint_circuit=circuit,
        disjoint_copy_edges=mixed,cubic_disjoint_copy_edges=new_mixed,extra=extra_metadata,
        mixed_address_edges=off_axes,cubic_displacement_mixed_edges=new_off_axes,
        coordinate_sha256=hashlib.sha256(json.dumps([[str(x) for x in p] for p in points]).encode()).hexdigest(),
        checks=checks,engine_checks=C.CHECKS)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['bank','circuit','two_circuits','connector'],default='bank')
    args=ap.parse_args();print(json.dumps(run(args.mode),separators=(',',':'),sort_keys=True))
