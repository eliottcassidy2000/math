#!/usr/bin/env python3
"""Independent bounded interface checks; NOT an incoming-packet replay."""
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial
import hashlib, json
from pathlib import Path


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def norm(v):
    x,y=v
    return x*x-x*y+y*y


def plus(v,w):
    return (v[0]+w[0],v[1]+w[1])


def polarization(v,w):
    return norm(plus(v,w))-norm(v)-norm(w)


def native(s):
    require(s>0 and s%2 and s%8!=5,'primitive source')
    a=2 if s%8==1 else 1
    v=(3*s+1)//(2**a); j=0
    while v%8==5:
        v=(v-1)//4; j+=1
    return v,a+2*j,a,j


def odd_word(s):
    word=[]
    while s!=1:
        z=3*s+1;a=0
        while z%2==0:z//=2;a+=1
        s=z;word.append(a)
    return word


seeds=[1,3,7,9,11,15,17,19,23,25,27,31,33,35,39,41,43,47,49,51,55,1009,1011,1015]
for s in seeds:
    y,t,a,j=native(s)
    require((3*s+1)//2**a==4**j*y+(4**j-1)//3,'B expansion')
    require(2**t*(3*y+1)==3*(3*s+1)+2**a,'corrected carry')
    b=F(3*(3*s+1),2**t*(3*y+1)); d=1-b*b
    ell=F(2**a,18*s+6+2**a)
    require(ell==(1-b)/(1+b),'price formula')
    for e in ((1,0),(0,1),(-2,-1)):
        for sign in (-1,1):
            f=tuple((sign-b)*v for v in e)
            r=tuple(d*v-b*w for v,w in zip(e,f))
            require(norm(e)+norm(f)/d==norm(e)+norm(r)/d,'constant-charge/sign-gauge energy')
            require(norm(f)/d==(ell if sign==1 else 1/ell)*norm(e),'signed forcing price')
require(native(3)==(1,3,1,1),'s3 P clock')
require(odd_word(3)==[1,4],'s3 U clocks')
word=odd_word(3); L=len(word)-1; K=sum((a-1)//2 for a in word)
source_probability=F(2*factorial(K)*factorial(L+1),factorial(L+K+2))
require(source_probability==F(1,6) and F(1,31)!=source_probability,'local price differs from source probability at 3')
require(all(((10*64**k-1)//3)%8==5 for k in range(1,9)),'U hostile n primitive exclusion')
require(all(((64**(k+1)-1)//3)%8==5 for k in range(9)),'U hostile z primitive exclusion')

# A deliberately small synthetic source ledger, not a native N-census.
# (source label, terminal label, exact clock, original Eisenstein charge)
records=[(1,'A',0,(1,0)),(3,'A',1,(1,0)),(7,'A',2,(0,1)),
         (9,'A',0,(1,1)),(11,'B',0,(-2,-1)),(15,'B',1,(1,2)),
         (17,'A',0,(0,0)),(19,'B',2,(1,-1))]
signed=[(s,u,tuple((-1)**t*v for v in z)) for s,u,t,z in records]
energies={}
for q in (1,2,4,8,16):
    bins={}
    for s,u,z in signed:
        key=(u,((s-1)//2)%q)
        bins[key]=plus(bins.get(key,(0,0)),z)
    energies[q]=sum(map(norm,bins.values()))
owner_counts={q:0 for q in (1,2,4,8)}
for q in (1,2,4,8):
    pair_debt=0
    for (s,u,z),(r,v,w) in combinations(signed,2):
        if u!=v:continue
        delta=abs((s-r)//2);owner=delta & -delta
        if owner==q:pair_debt+=polarization(z,w);owner_counts[q]+=1
    require(energies[q]-energies[2*q]==pair_debt,'first-differing-bit debt')
source=sum(norm(z) for s,u,z in signed)
require(energies[16]==source,'source separation')
deltas=[energies[q]-energies[2*q] for q in (1,2,4,8)]
require(energies[1]==source+sum(deltas),'telescoping debt')
require(energies[1]<=source+sum(max(d,0) for d in deltas),'positive-level upper bound')


def defect_stats(word, internal_cycle=False):
    n=len(word); edges=[]; parents=list(range(n))
    def root(i):
        while parents[i]!=i:i=parents[i]
        return i
    for i,j in combinations(range(n),2):
        ci,li=word[i];cj,lj=word[j]
        forward=((cj-ci)%3==1) if ci!=cj else (((lj-li)%3==1) if internal_cycle else li<lj)
        if not forward:
            edges.append((i,j));ri,rj=root(i),root(j);parents[ri]=rj
    components=len({root(i) for i in range(n)})
    return len(edges),len(edges)-n+components


def interleavings(counts,prefix=()):
    if not sum(counts):yield prefix;return
    for color,count in enumerate(counts):
        if count:
            left=list(counts);left[color]-=1
            yield from interleavings(tuple(left),prefix+(color,))


def labelled(colors):
    seen=[0]*3;word=[]
    for c in colors:word.append((c,seen[c]));seen[c]+=1
    return word

clone={}
for sizes in ((1,2,3),(2,3,3),(3,3,3),(3,3,4)):
    rows=[defect_stats(labelled(w)) for w in interleavings(sizes)]
    edge=min(e for e,b in rows);cycle=min(b for e,b in rows)
    require(edge==sizes[0]*sizes[1],'clone edge minimum')
    require(cycle==max(sizes[0]-2,0)*(sizes[1]-1),'clone cycle minimum')
    clone[str(sizes)]={'interleavings':len(rows),'minimum_edges':edge,'minimum_cycles':cycle,
                      'cycles_at_min_edges':sorted({b for e,b in rows if e==edge}),
                      'edges_at_min_cycles':sorted({e for e,b in rows if b==cycle})}
require(clone['(3, 3, 4)']['cycles_at_min_edges']==[4],'edge/cycle optimizer separation')
require(10 in clone['(3, 3, 4)']['edges_at_min_cycles'],'ten-edge beta-two witness')
tower_order=[(0,0)]+[(1,j) for j in range(3)]+[(2,j) for j in range(3)]+[(0,1),(0,2)]
require(defect_stats(tower_order,True)==(12,4),'distinct order-nine tower witness')

for x,y in product(range(-2,3),repeat=2):
    xp,yp=3*x-y,x+2*y
    require(norm((xp,yp))==7*norm((x,y)),'norm-seven matrix')
    require((2*xp+yp)%7==0,'index-seven image')
    require(F(5*xp-yp,7)==2*x-y,'transformed readout')
    require((F(2*xp+yp,7),F(-xp+3*yp,7))==(x,y),'inverse matrix')
require(norm((3,1))==7 and norm((4,-1))==21,'norm shells')
# (1-w)(3+w) = 4-w, using w^2=-w-1.
require((3+1,1-3+1)==(4,-1),'factorization')

result={'scope':'bounded independent interface checks only; no full packet or large scale census replay',
        'native_seed_controls':len(seeds),'s3':{'P_braid_clock':3,'U_valuations':[1,4],'U_shortcut_clock':5,'U_literal_clock':7,'local_price':'1/31','U_source_probability':'1/6'},
        'synthetic_pair_ledger':{'sources':len(records),'energies':energies,'level_debts':deltas,'pair_owner_counts':owner_counts,'source_norm':source},
        'clone_controls':clone,'tower9_constructed_frame':{'edges':12,'cycles':4},
        'norm_matrix_vectors':25,'status':'PASS'}
print(json.dumps(result,sort_keys=True,indent=2))
