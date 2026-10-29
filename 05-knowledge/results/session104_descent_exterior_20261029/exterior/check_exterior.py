#!/usr/bin/env python3
"""Exact, stdlib-only checks; all gates survive python -O."""
from fractions import Fraction as F
from itertools import combinations, permutations
import json
from pathlib import Path

def gate(x,msg):
    if not x: raise RuntimeError(msg)
def T(a): return list(map(list,zip(*a)))
def mm(a,b): return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]
def eye(n): return [[int(i==j) for j in range(n)] for i in range(n)]
def rank(a):
    a=[list(map(F,r)) for r in a]; k=0
    for j in range(len(a[0])):
        p=next((i for i in range(k,len(a)) if a[i][j]),None)
        if p is None: continue
        a[k],a[p]=a[p],a[k]; z=a[k][j]; a[k]=[x/z for x in a[k]]
        for i in range(len(a)):
            if i!=k:
                z=a[i][j]; a[i]=[x-z*y for x,y in zip(a[i],a[k])]
        k+=1
        if k==len(a): break
    return k
pairs=list(combinations(range(4),2))
def wedge(a):
    return [[a[i][k]*a[j][l]-a[i][l]*a[j][k] for k,l in combinations(range(len(a[0])),2)] for i,j in combinations(range(len(a)),2)]
def det(a):
    n=len(a)
    return sum((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))*prod(a[i][p[i]] for i in range(n)) for p in permutations(range(n)))
def prod(xs):
    z=1
    for x in xs:z*=x
    return z
def transition(images): return [[int(images[j]==i) for j in range(len(images))] for i in range(len(images))]
# state order safe A,B0,B1; failed A,B; illegal X
trans={ 'a':[0,5,5,3,5,5], 'b':[2,5,5,4,5,5], 'c':[5,0,0,5,3,5], 'd':[5,1,1,5,4,5], 'e':[5,2,4,5,4,5] }
ms={a:transition(v) for a,v in trans.items()}
# Reachability and Moore distinguishability, not just number of named states.
seen={1}; change=True
while change:
    new=seen|{v[s] for v in trans.values() for s in seen}; change=new!=seen; seen=new
gate(len(seen)==6,'Moore reachability')
classes=[0,0,0,1,1,2]
while True:
    sig=[(classes[s],tuple(classes[v[s]] for v in trans.values())) for s in range(6)]
    ids={v:i for i,v in enumerate(dict.fromkeys(sig))}; new=[ids[v] for v in sig]
    if new==classes:break
    classes=new
gate(len(set(classes))==6,'Moore minimality')
quot=[0,1,2,3,3,3]; acc={}
for a,v in trans.items():
    out=[]
    for q in range(4):
        images={quot[v[s]] for s in range(6) if quot[s]==q}
        gate(len(images)==1,'acceptance quotient well-defined'); out.append(images.pop())
    acc[a]=transition(out)
    gate(rank(ms[a])==3 and rank(acc[a])==2 and rank(wedge(acc[a]))==1,'letter rank obstruction')
e2=mm(ms['e'],ms['e']); gate(rank(e2)==2,'universal exterior-square obstruction')
# K4 boundary: edge i<j has boundary e_j-e_i.
boundary=[[int(v==j)-int(v==i) for i,j in pairs] for v in range(4)]
gate(rank(boundary)==3,'boundary rank')
H=[[-1,-1,-1],[1,0,0],[0,1,0],[0,0,1]]
cycles=wedge(H); cuts=mm(T(boundary),H)
gate(mm(boundary,cycles)==[[0]*3 for _ in range(4)],'cycles exact')
gate(rank(cycles)==3 and rank(cuts)==3,'cut/cycle dimensions')
gate(mm(T(cuts),cycles)==[[0]*3 for _ in range(3)],'orthogonal cut/cycle')
gate(mm(boundary,cuts)==[[4*x for x in r] for r in H],'Laplacian eigenvalue 4')
basis=[cuts[i]+cycles[i] for i in range(6)]; gate(rank(basis)==6,'rational splitting')
# Hodge star for e0,e1,e2,e3, Euclidean standard metric.
star=[[0]*6 for _ in range(6)]
for j,p in enumerate(pairs):
    q=tuple(i for i in range(4) if i not in p); order=p+q
    star[pairs.index(q)][j]=(-1)**sum(order[i]>order[k] for i in range(4) for k in range(i+1,4))
gate(mm(star,star)==eye(6),'star square')
for s in [-1,1]:gate(rank([[eye(6)[i][j]+s*star[i][j] for j in range(6)] for i in range(6)])==3,'star eigenspace')
gate(mm(boundary,mm(star,cuts))==[[0]*3 for _ in range(4)],'star exchanges cuts/cycles')
# Honest synthetic lift: 1 direct-sum C, not actual grammar transition.
C=[[1,0,1],[1,1,1],[1,1,0]]
L=[[1,0,0,0]]+[[0]+r for r in C]; W=wedge(L); Wc=wedge(C)
expected=[[0]*6 for _ in range(6)]
for i in range(3):
    for j in range(3):expected[i][j]=C[i][j]; expected[i+3][j+3]=Wc[i][j]
gate(W==expected,'synthetic 3+3 lift')
gate(det(C)==-1,'det C')
for t in range(-4,5):
    char=lambda A:det([[t*int(i==j)-A[i][j] for j in range(3)] for i in range(3)])
    gate(char(C)==t**3-2*t*t-t+1,'C polynomial')
    gate(char(Wc)==t**3+t*t-2*t-1,'dual polynomial')
# Integral duality and the order-three cyclotomic root transformation.
D=[[C[i][j]-int(i==j) for j in range(3)] for i in range(3)]
S=[[0,-1,2],[-1,-2,4],[-1,-3,5]]
gate(det(S)==1 and mm(Wc,S)==mm(S,D),'integral exterior duality')
D2=mm(D,D); phi=[[D2[i][j]-2*int(i==j) for j in range(3)]for i in range(3)]
gate(mm(phi,C)==[[-int(i==j)for j in range(3)]for i in range(3)],'Mobius equals quadratic root transform')
iterate=D
for _ in range(3):
    square=mm(iterate,iterate)
    iterate=[[square[i][j]-2*int(i==j)for j in range(3)]for i in range(3)]
gate(iterate==D and phi!=D,'cyclotomic root transform order three')
# Determinantal Pluecker identity for all small integer potential quadruples.
from itertools import product
for x in product(range(-2,3),repeat=4):
    d=[x[j]-x[i] for i,j in pairs]
    gate(d[0]*d[5]-d[1]*d[4]+d[2]*d[3]==0,'difference Pluecker relation')
# Constructive universality: every primitive positive abc triple is carried.
from math import gcd
abc_cases=0
for A in range(1,101):
    for B in range(A,101):
        if gcd(A,B)!=1: continue
        C=A+B
        if A==B:
            gate((A,B,C)==(1,1,2),'primitive equal boundary')
            x=(0,2,3,6)
        else:
            d=B-A; x=(0,d,2*d,2*B)
        gate(all(x[i]<x[i+1] for i in range(3)),'strictly increasing abc potentials')
        ds=[x[j]-x[i] for i,j in pairs]
        products=(ds[0]*ds[5],ds[2]*ds[3],ds[1]*ds[4])
        g=gcd(gcd(products[0],products[1]),products[2])
        gate(tuple(z//g for z in products)==(A,B,C),'primitive abc carrier reconstruction')
        if A<B:gate(g==2*(B-A),'exact abc common factor')
        abc_cases+=1
result={'integral_duality_intertwiner':S,'cyclotomic_root_transform_order':3,'primitive_abc_carrier_cases':abc_cases,'status':'PASS','moore_states':6,'letter_ranks':{a:rank(m) for a,m in ms.items()},'e_squared_rank':rank(e2),'acceptance_letter_ranks':{a:rank(m) for a,m in acc.items()},'acceptance_wedge_letter_ranks':{a:rank(wedge(m)) for a,m in acc.items()},'boundary_rank':rank(boundary),'cycle_dimension':rank(cycles),'cut_dimension':rank(cuts),'hodge_eigenspace_dimensions':[3,3],'synthetic_lift':'C direct-sum exterior_square(C); not Moore representation'}
print(json.dumps(result,indent=2))
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
