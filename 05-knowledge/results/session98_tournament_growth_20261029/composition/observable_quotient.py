#!/usr/bin/env python3
"""Exact univariate surplus-jet observation ranks; no symbolic dependencies."""
from collections import defaultdict,Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json,argparse
from inherited_recreate import partitions

def require(ok,msg):
    if not ok:raise RuntimeError(msg)

@lru_cache(None)
def column(part,h):
    if not part:return {(0,()):1}
    out=defaultdict(int);d=part[0]
    for (r,bs),v in column(part[1:],h-d).items():
        for i in range(1,(h-r)//d+1):
            out[r+d*i,tuple(sorted(bs+(i,)))]+=v
    return dict(out)

@lru_cache(None)
def reduce_power(e,cap):
    if e<=cap:return {e:1}
    # relation product_(j=0)^cap (b-j), monic of degree cap+1.
    rel=[1]
    for j in range(cap+1):
        new=[0]*(len(rel)+1)
        for i,v in enumerate(rel):new[i]-=j*v;new[i+1]+=v
        rel=new
    result={0:1}
    for _ in range(e):
        nxt=defaultdict(int)
        for degree,v in result.items():
            if degree+1<=cap:nxt[degree+1]+=v
            else:
                for j,c in enumerate(rel[:-1]):nxt[j]-=v*c
        result={i:v for i,v in nxt.items() if v}
    return result

def grid_column(col):
    out=defaultdict(int)
    caps={1:1,2:2,3:4} # for h<=12, all other tournament-count caps exceed degree.
    for (r,bs),v in col.items():
        poly={():v}
        for i,e in Counter(bs).items():
            powers=reduce_power(e,caps[i]) if i in caps else {e:1}
            nxt=defaultdict(int)
            for mon,c in poly.items():
                for p,q in powers.items():nxt[tuple(sorted(mon+(i,)*p))]+=c*q
            poly=nxt
        for b,c in poly.items():out[r,b]+=c
    return {q:v for q,v in out.items() if v}

def echelon(columns):
    basis={};kernel=[]
    for j,col in enumerate(columns):
        v={i:F(x) for i,x in col.items() if x};rel={j:F(1)}
        while v:
            p=min(v)
            if p not in basis:
                scale=v[p];v={i:x/scale for i,x in v.items()};rel={i:x/scale for i,x in rel.items()}
                basis[p]=(v,rel);break
            bv,br=basis[p];factor=v[p]
            for dst,src in ((v,bv),(rel,br)):
                for i,x in src.items():
                    dst[i]=dst.get(i,F(0))-factor*x
                    if not dst[i]:del dst[i]
        if not v:
            require(bool(rel),'empty dependence')
            kernel.append(rel)
    # Every returned dependence is checked against original columns.
    for rel in kernel:
        summed=defaultdict(F)
        for j,v in rel.items():
            for row,w in columns[j].items():summed[row]+=v*w
        require(all(not v for v in summed.values()),'bad kernel')
    return len(basis),kernel

def modular_rank(columns,p):
    basis={}
    for col in columns:
        v={r:int(c)%p for r,c in col.items() if int(c)%p}
        while v:
            pivot=max(v)
            if pivot not in basis:
                inv=pow(v[pivot],-1,p)
                basis[pivot]={r:c*inv%p for r,c in v.items()}
                break
            factor=v[pivot]
            for r,c in basis[pivot].items():
                v[r]=(v.get(r,0)-factor*c)%p
                if not v[r]:del v[r]
    return len(basis)


def encoded(rel,parts):
    return {','.join(map(str,parts[j])) or 'constant':str(v) for j,v in rel.items()}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--h',type=int,default=12);args=parser.parse_args()
    require(0<=args.h<=12,'ordinary-grid caps audited only through12')
    rows=[]
    for h in range(args.h+1):
        parts=[p for r in range(h+1) for p in partitions(r)]
        cols=[column(p,h) for p in parts];gcols=[grid_column(c) for c in cols]
        rank,ker=echelon(cols);grank,gker=echelon(gcols)
        for prime in (1000003,1000033):
            require(modular_rank(cols,prime)==rank,'formal modular lower bound differs')
            require(modular_rank(gcols,prime)==grank,'grid modular lower bound differs')
        rows.append({'h':h,'dimension':len(parts),'formal_rank':rank,'ordinary_palette_rank':grank,
          'formal_kernel':[encoded(k,parts) for k in ker],
          'ordinary_palette_kernel':[encoded(k,parts) for k in gker]})
        # Compatibility with truncating the observation horizon.
        if h:
            for p,c,g in zip(parts,cols,gcols):
                require({q:v for q,v in c.items() if q[0]<h}==column(p,h-1),'formal projection')
                require({q:v for q,v in g.items() if q[0]<h}==grid_column(column(p,h-1)),'grid projection')
    # Direct evaluation versus reduced polynomial on the complete relevant low grid.
    for h in range(min(6,args.h)+1):
        for p in [p for r in range(h+1) for p in partitions(r)]:
            c=column(p,h);gc=grid_column(c)
            for b1,b2,b3 in product(range(2),range(3),range(5)):
                bs={1:b1,2:b2,3:b3,4:2,5:3,6:1}
                def val(poly):
                    out=defaultdict(int)
                    for (r,mon),v in poly.items():
                        for i in mon:v*=bs[i]
                        out[r]+=v
                    return {r:v for r,v in out.items() if v}
                require(val(c)==val(gc),'grid evaluation check')
    print(json.dumps({'status':'PASS','independent_modular_rank_primes':[1000003,1000033],'ranks':rows},indent=2))
if __name__=='__main__':main()
