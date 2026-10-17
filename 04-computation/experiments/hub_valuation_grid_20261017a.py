#!/usr/bin/env python3
"""Pay the prime-19 obstruction explicitly, with literal old source ports.

V=Gamma9 union {b*conj(a), -b*a}; grid V+uV+u^2V. The guaranteed
mixed edge 0--b*u^3 really retains the even hub extension in displacement.
Full induced edges use the independent w-basis field of degree 480.
"""
import argparse,hashlib,json
from itertools import product
from fractions import Fraction as F
import cubic_hub_coupling_20261017a as H
import hub_short_circuit_20261017a as S
C=S.C

def embed(h):
    p=S.O
    for k,v in enumerate(h):
        if not v:continue
        block=k&3;coef=C.I
        if k&4:
            coef=C.lift(C.G.b);v*=F(1,11)
            if block&2:v*=11
            block^=2
        p=S.add(p,(0,)*(120*block)+C.sc(coef,v)+(0,)*(120*(3-block)))
    return p

def run(mode):
    G=H.hub();b=G[2]
    directions=[b] if mode=='one' else G[1:5]
    V=list(dict.fromkeys(G+[p for h in directions for p in (H.mul(h,H.conj(H.a)),H.scale(H.mul(h,H.a),-1))]))
    need=S.need
    need(len(V)==9+2*len(directions),'two denominator-bearing coordinates per direction')
    jb=C.lift(C.G.b);w=C.W
    u=S.lift(C.sc(C.add(C.add(C.sc(w,11),C.sc(C.I,-6)),
                       C.mul(C.add(C.add(C.mul(w,w),C.sc(w,-1)),C.sc(C.I,-12)),jb)),F(1,38)))
    u2=S.mul(u,u);u3=S.mul(u2,u)
    need(S.unit(u),'independent w-basis unit')
    need(S.add(S.add(u3,S.mul(embed(H.a),u2)),S.sc(S.mul(embed(H.conj(H.a)),u),-1))==S.I,
         'cubic relation in independent w basis')
    base=[embed(x) for x in V];ports=[S.lift(C.lift(p)) for p in C.G.ports]
    points=ports[:];addresses={}
    for t in product(range(len(V)),repeat=3):
        p=S.add(S.add(base[t[0]],S.mul(u,base[t[1]])),S.mul(u2,base[t[2]]))
        need(p not in addresses,'cubic grid has unique base coordinates')
        addresses[p]=t;points.append(p)
    points=list(dict.fromkeys(points));need(points[:8]==ports,'all old ports unchanged')
    target=S.mul(embed(b),u3)
    need(target in addresses and S.unit(target) and any(target[120:]),'positive new-hub mixed edge')
    E,meta=S.graph(points);mixed=[];extension=[];details=[];directions=set()
    for i,k in E:
        if points[i] not in addresses or points[k] not in addresses:continue
        count=sum(x!=y for x,y in zip(addresses[points[i]],addresses[points[k]]))
        need(count in (1,3),'no two-coordinate unit displacement')
        if count==3:
            mixed.append([i,k])
            delta=tuple(H.sub(V[x],V[y]) for x,y in zip(addresses[points[i]],addresses[points[k]]))
            need(H.unit_direct(delta),'coupled edge checked in independent u quotient')
            negative=tuple(H.scale(x,-1) for x in delta);directions.add(min(delta,negative))
            details.append(dict(edge=[i,k],source=addresses[points[i]],target=addresses[points[k]],
                                coefficients=[[str(c) for c in x] for x in delta]))
            if any(S.sub(points[i],points[k])[120:]):extension.append([i,k])
    colouring=None;bound=None
    for q in range(3,7):
        colouring=C.colour(len(points),E,q)
        if colouring is not None:bound=q;break
    need(colouring is not None and all(colouring[i]!=colouring[k] for i,k in E),'proper colouring replay')
    base_colours=[0,1,1,1,1,0,0,2,2]+[0]*(len(V)-9)
    formula=[sum(base_colours[t] for t in addresses[p])%3 if p in addresses else None for p in points]
    need(all(formula[i]!=formula[k] for i,k in E if formula[i] is not None and formula[k] is not None),
         'explicit additive product colouring of the entire cubic grid')
    free=[i for i,c in enumerate(formula) if c is None]
    need(len(free)==6,'only six old ports remain outside the grid')
    port_edges=[(i,k) for i,k in E if i in free or k in free]
    for choices in product(range(3),repeat=len(free)):
        for i,c in zip(free,choices):formula[i]=c
        if all(formula[i]!=formula[k] for i,k in port_edges):break
    else:raise RuntimeError('additive grid colouring does not extend to original ports')
    need(all(formula[i]!=formula[k] for i,k in E),'independent formula-plus-port colouring replay')
    boundary=S.boundary(len(points),E)
    return dict(mode=mode,vertices=len(points),edges=E,colour_bound=bound,colouring=colouring,
                mixed_edges=mixed,extension_mixed_edges=extension,boundary=boundary,
                mixed_details=details,unoriented_mixed_directions=len(directions),
                product_base_colours=base_colours,product_port_colours=formula[:8],
                geometry=meta,base_vertices=[[str(v) for v in x] for x in V],
                guaranteed_edge=[points.index(S.O),points.index(target)],
                guaranteed_address=addresses[target],checks=S.checks,
                coordinate_hash=hashlib.sha256(json.dumps([[str(v) for v in p] for p in points]).encode()).hexdigest())

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['one','four'],default='one')
    print(json.dumps(run(parser.parse_args().mode),sort_keys=True,separators=(',',':')))
