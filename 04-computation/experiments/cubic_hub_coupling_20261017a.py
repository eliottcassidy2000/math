#!/usr/bin/env python3
"""Exact elimination of the middle cubic coordinate; finite hub alphabets.

Base B=Q(sqrt7,sqrt11,i), basis indexed by square-root bit masks 0..7.
The cubic basis is 1,u,u^2, with u^3+a u^2-conj(a)u-1=0.
All modular filters reject only. Every accepted displacement is replayed
using the independent degree-24 quotient algebra.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import product
import json
from math import isqrt
import numpy as np

O=(F(0),)*8
I=(F(1),)+(F(0),)*7
checks=0
def need(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def scale(x,c):return tuple(c*a for a in x)
def sub(x,y):return add(x,scale(y,-1))
def conj(x):return tuple(-a if k&4 else a for k,a in enumerate(x))
TABLE=[[(a^b,(7 if a&b&1 else 1)*(11 if a&b&2 else 1)*(-1 if a&b&4 else 1)) for b in range(8)] for a in range(8)]
def mul(x,y):
    out=[F(0)]*8
    for a,xa in enumerate(x):
        if not xa:continue
        for b,yb in enumerate(y):
            if yb:
                k,c=TABLE[a][b];out[k]+=c*xa*yb
    return tuple(out)
def norm(x):return mul(x,conj(x))
def basis(k):return tuple(F(t==k) for t in range(8))
j=basis(6)
a=scale(add(scale(I,20),scale(j,11)),F(1,19))

def inverse_real(x):
    need(x!=O and all(v==0 for v in x[4:]),'nonzero real denominator')
    p=I
    for flip in (1,2,3):
        p=mul(p,tuple(-c if (k&flip).bit_count()%2 else c for k,c in enumerate(x)))
    den=mul(p,x)
    need(den[0]!=0 and all(v==0 for v in den[1:]),'field norm of denominator')
    return scale(p,1/den[0])

def middle(A,C):
    delta=sub(norm(A),norm(C))
    if delta==O:return None
    K=sub(mul(a,mul(C,conj(A))),mul(conj(C),A))
    return mul(sub(mul(A,K),mul(C,conj(K))),inverse_real(delta))

def metric(A,B,C):
    D1=add(mul(B,conj(A)),mul(C,conj(B)));D2=mul(C,conj(A))
    return (add(sub(D1,mul(a,D2)),conj(D2))==O and
            add(add(add(norm(A),norm(B)),norm(C)),
                add(mul(conj(a),D2),mul(a,conj(D2))))==I)

def cubic_mul(x,y):
    v=[O for _ in range(5)]
    for r in range(3):
        for s in range(3):v[r+s]=add(v[r+s],mul(x[r],y[s]))
    for k in (4,3):
        for t,c in enumerate((I,conj(a),scale(a,-1))):v[k-3+t]=add(v[k-3+t],mul(v[k],c))
    return tuple(v[:3])

def cubic_conj(x):
    # u^-1 = u^2+a*u-conj(a), authenticated by the cubic relation.
    inv=(scale(conj(a),-1),a,I);inv2=cubic_mul(inv,inv)
    return tuple(add(add(mul(conj(x[0]),I if k==0 else O),mul(conj(x[1]),inv[k])),
                     mul(conj(x[2]),inv2[k])) for k in range(3))

def unit_direct(x):return cubic_mul(x,cubic_conj(x))==(I,O,O)

def hub():
    z=O;aa=I;bb=scale(sub(basis(1),scale(basis(4),3)),F(1,4));cc=scale(basis(4),-1)
    dd=scale(add(scale(I,3),basis(5)),F(1,4))
    return [z,aa,bb,cc,dd,add(aa,bb),add(aa,cc),add(bb,dd),add(cc,dd)]

def mask19(x):
    mask=0
    for k,(r,s) in enumerate(product((8,-8),(7,-7))):
        re=im=0
        for idx,c in enumerate(x):
            value=(c.numerator*pow(c.denominator,-1,19))%19
            value=value*(r if idx&1 else 1)*(s if idx&2 else 1)%19
            if idx&4:im+=value
            else:re+=value
        if re%19==0 and im%19==0:mask|=1<<k
    return mask

def field_filter(prime,D):
    r=next(t for t in range(prime) if t*t%prime==7)
    s=next(t for t in range(prime) if t*t%prime==11)
    ii=next(t for t in range(prime) if t*t%prime==prime-1)
    def ev(x,sign):
        return sum((v.numerator*pow(v.denominator,-1,prime))*(r if k&1 else 1)*
                   (s if k&2 else 1)*(sign*ii if k&4 else 1) for k,v in enumerate(x))%prime
    vals=np.array([(ev(x,1),ev(x,-1)) for x in D],dtype=np.int64)
    return prime,vals,(ev(a,1),ev(a,-1)),{tuple(row) for row in vals}

def controls():
    need(mul(j,j)==scale(I,-11),'literal j square')
    need(unit_direct((O,I,O)),'u is norm one')
    # u^3 has all three coordinates; it tests the solver independently of the search alphabet.
    v=(I,conj(a),scale(a,-1))
    need(unit_direct(v) and metric(*v),'nonempty coupled positive control')
    need(middle(v[0],v[2])==v[1],'middle coordinate recovered exactly')
    need(not metric(v[0],add(v[1],I),v[2]),'changed middle coordinate rejected')
    # Independent quotient and metric agree over the whole declared small box.
    for x in product((O,I,scale(I,-1),j),repeat=3):need(unit_direct(x)==metric(*x),'metric versus quotient')
    roots=[t for t in range(19) if (t**4-55*t*t+605)%19==0]
    need(roots==[1,4,15,18] and pow(8,9,19)==18,'incoming all-finite infinity claim hostile')
    need(all((4*t**3-110*t)%19 for t in roots),'four simple degree-one primes')
    D=10;values=[D*(1-D),D*(4-D),(1-D)*(D+3)]
    square_roots=[next(r for r in range(1,19) if (11*r*r-v)%19==0) for v in values]
    need(all((22*r)%19 for r in square_roots),'finite local point Hensel derivatives')
    return dict(F4_roots_mod19=roots,minus11_legendre_mod19=-1,
                finite_local_point_mod19=dict(D=D,z=square_roots))

def run(mode):
    G=hub()
    if mode=='hub':V=G
    elif mode=='jhub':V=list(dict.fromkeys(add(x,mul(j,y)) for x in G for y in G))
    elif mode=='jhub8':V=list(dict.fromkeys(add(x,scale(mul(j,y),8)) for x in G for y in G))
    elif mode=='jhub19':V=list(dict.fromkeys(add(x,scale(mul(j,y),19)) for x in G for y in G))
    else:raise ValueError(mode)
    D=sorted({sub(x,y) for x in V for y in V});masks=[mask19(x) for x in D]
    nonzero=[i for i,x in enumerate(D) if x!=O]
    pairs=[(i,k) for i in nonzero for k in nonzero if masks[i]|masks[k]==15]
    # These primes split all three quadratic generators and avoid all input denominators.
    primes=[p for p in range(503,2000) if all(p%d for d in range(2,isqrt(p)+1))
            and p%4==1 and pow(7,(p-1)//2,p)==pow(11,(p-1)//2,p)==1][:3]
    need(len(primes)==3,'three valid independent screens')
    filters=[field_filter(p,D) for p in primes]
    survivors=[];singular=0
    for ia,ic in pairs:
        keep=True
        for p,val,(av,ac),Dset in filters:
            A,Ab=map(int,val[ia]);C,Cb=map(int,val[ic]);delta=(A*Ab-C*Cb)%p
            K=(av*C*Ab-Cb*A)%p;Kb=(ac*Cb*A-C*Ab)%p
            if not delta:
                if (A*K-C*Kb)%p:keep=False;break
                continue
            inv=pow(delta,-1,p)
            B=(A*K-C*Kb)*inv%p;Bb=(Ab*Kb-Cb*K)*inv%p
            if (B,Bb) not in Dset or (A*Ab+B*Bb+C*Cb+ac*C*Ab+av*Cb*A-1)%p:
                keep=False;break
        if keep:survivors.append((ia,ic))
    exact=[];Dindex={x:i for i,x in enumerate(D)}
    for ia,ic in survivors:
        A,C=D[ia],D[ic];B=middle(A,C)
        if B is None:
            singular+=1
            # Do not silently discard a singular pair: compatibility must fail.
            K=sub(mul(a,mul(C,conj(A))),mul(conj(C),A))
            need(mul(A,K)!=mul(C,conj(K)),'singular compatibility forbidden by irreducibility')
            continue
        if B in Dindex and B!=O and metric(A,B,C):
            need(unit_direct((A,B,C)),'every accepted mixed edge independently replayed')
            exact.append([ia,Dindex[B],ic])
    base_edges=[(i,k) for i in range(len(V)) for k in range(i+1,len(V)) if norm(sub(V[i],V[k]))==I]
    hub_colours=[0,1,1,1,1,0,0,2,2]
    colours=hub_colours if mode=='hub' else [c for c in hub_colours for _ in G]
    need(len(colours)==len(V) and all(colours[i]!=colours[k] for i,k in base_edges),
         'explicit base graph colouring, independently checked on every pair')
    return dict(mode=mode,vertices_in_base=len(V),difference_count=len(D),zero_mask_histogram=dict(Counter(masks)),
                ordered_nonzero_pairs=len(nonzero)**2,after_fixed_prime_gate=len(pairs),filter_primes=primes,
                after_modular_filters=len(survivors),singular_survivors=singular,mixed_directions=exact,
                differences=[[str(v) for v in x] for x in D],base_vertices=[[str(v) for v in x] for x in V],
                base_edges=base_edges,base_colouring=colours)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['hub','jhub','jhub8','jhub19'],default='jhub');args=p.parse_args()
    hostile=controls();out=run(args.mode)
    print(json.dumps(dict(status='FINITE-EXACT complete difference-pair search in specified alphabet',
          run=out,controls=hostile,checks=checks),sort_keys=True,separators=(',',':')))
