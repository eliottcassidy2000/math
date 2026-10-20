#!/usr/bin/env python3
"""Independent matrix-action checks of marked-flip output and finite group identities.
No repository mathematical imports. Does not use the C++ affine-mask implementation.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import gcd, lcm
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent

def require(test, why):
    if not test: raise RuntimeError(why)

def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inverse(a):
    b=[0]*len(a)
    for i,x in enumerate(a):b[x]=i
    return tuple(b)
def order(a):
    seen=set();o=1
    for i in range(len(a)):
        if i in seen:continue
        j=i;k=0
        while j not in seen:seen.add(j);k+=1;j=a[j]
        o=lcm(o,k)
    return o

def matrix(n,s):
    a=[[0]*n for _ in range(n)]
    for j in range(1,n):a[0][j]=1
    k=0
    for i in range(1,n):
        for j in range(i+1,n):
            a[i][j]=1-((s>>k)&1); a[j][i]=1-a[i][j];k+=1
    return a

def normalized(a):
    n=len(a);b=[row[:] for row in a]
    switch=[not b[0][i] if i else False for i in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if switch[i]!=switch[j]:b[i][j]^=1;b[j][i]^=1
    require(all(b[0][j]==1 for j in range(1,n)), 'root normalization')
    s=k=0
    for i in range(1,n):
        for j in range(i+1,n):s|=(1-b[i][j])<<k;k+=1
    return s

def relabel(a,p):
    n=len(a);b=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):b[p[i]][p[j]]=a[i][j]
    return normalized(b)

def flip(a,g):
    b=[r[:] for r in a]
    for i,j in enumerate(g):
        if i<j:b[i][j]^=1;b[j][i]^=1
    return normalized(b)

def centralizer_fpf(g):
    pairs=[(i,j) for i,j in enumerate(g) if i<j]
    for targets in permutations(pairs):
        for bits in product((0,1),repeat=len(pairs)):
            h=[None]*len(g)
            for (a,b),(c,d),bit in zip(pairs,targets,bits):
                if bit:c,d=d,c
                h[a]=c;h[b]=d
            h=tuple(h)
            require(mul(h,g)==mul(g,h),'centralizer construction')
            yield h

def group_from_gens(gens):
    identity=tuple(range(len(gens[0])));seen={identity};todo=[identity]
    while todo:
        x=todo.pop()
        for g in gens:
            y=mul(g,x)
            if y not in seen:seen.add(y);todo.append(y)
    return sorted(seen)

def group_audit(G):
    G=set(G);remaining={x for x in G if order(x)>1 and order(x)&(order(x)-1)==0};weights=[];rows=[]
    while remaining:
        x=min(remaining);cl={mul(mul(h,x),inverse(h)) for h in G};remaining-=cl
        C={h for h in G if mul(h,x)==mul(x,h)}
        weight=Fraction(sum(order(h)%2 for h in C),len(C));weights.append(weight)
        rows.append({'element_order':order(x),'class_size':len(cl),'centralizer_size':len(C),'weight':str(weight)})
    lhs=Fraction(sum(order(x)%2==0 for x in G),len(G));require(sum(weights)==lhs,'general weighted identity')
    return {'order':len(G),'even_proportion':str(lhs),'classes':rows}

j=json.loads((ROOT/'independent_carrier.json').read_text());checks=[];pair_details=[];alltests=0
for c in j['cases']:
    n=c['n'];fixed=0;marks={}
    for r in c['rows']:
        s=r['state'];a=matrix(n,s); inv=[x for x in r['two_element_classes'] if x['order']==2]
        for x,flag,target in zip(inv,r['decorated_fixed'],r['flip_targets']):
            g=tuple(x['element']); require(relabel(a,g)==s,'marked involution must fix class')
            t=flip(a,g); candidates=list(centralizer_fpf(g));witnesses=[]
            for h in candidates:
                alltests+=1
                if relabel(a,h)==t:witnesses.append(h)
            require(bool(witnesses)==bool(flag),'marked fixed flag disagrees')
            fixed+=bool(witnesses)
            key=(s,g);tar=(target['state'],tuple(target['involution']));marks[key]=(tar,Fraction(x['weight']),bool(flag))
            pair_details.append({'n':n,'state':s,'involution':g,'flipped_state':t,'centralizer_size_in_Sn':len(candidates),'commuting_isomorphism_count':len(witnesses),'one_witness':witnesses[:1]})
    for key,(target,weight,flag) in marks.items():
        require(marks[target][0]==key,'marked flip squared')
        require(marks[target][1]==weight,'marked flip preserves centralizer weight')
        require((target==key)==flag,'marked permutation fixed-point flag')
    require(fixed==c['decorated_fixed'],'total fixed')
    checks.append({'n':n,'involution_pairs':len(marks),'fixed':fixed,'transpositions':(len(marks)-fixed)//2,'involution_mass':str(sum(x[1] for x in marks.values())),'fixed_mass':str(sum(x[1] for x in marks.values() if x[2]))})
# Entire S8 stabilizer and unmarked isomorphism witnesses for the decisive D8 hostile.
n=8;s=4722;a=matrix(n,s);G=[p for p in permutations(range(n)) if relabel(a,p)==s]
require(len(G)==8,'D8 hostile stabilizer');require(Counter(order(g) for g in G)=={1:1,2:5,4:2},'D8 order signature')
row=next(r for r in j['cases'][-1]['rows'] if r['state']==s)
reflection_rows=[x for x in row['two_element_classes'] if x['order']==2 and x['centralizer_size']==4]
hostiles=[]
for x in reflection_rows:
    g=tuple(x['element']);t=flip(a,g);h=next(p for p in permutations(range(n)) if relabel(a,p)==t)
    transported=mul(mul(inverse(h),g),h)
    same={mul(mul(k,g),inverse(k)) for k in G}
    require(transported not in same,'coarse witness changed marked conjugacy class')
    hostiles.append({'state':s,'g':g,'flipped_state':t,'coarse_isomorphism':h,'transported_g_back':transported,'commutes':mul(h,g)==mul(g,h),'weight':x['weight']})
D8=group_from_gens([(1,2,3,0),(0,3,2,1)])
S4=list(permutations(range(4)))
A5=[p for p in permutations(range(5)) if sum(p[i]>p[k] for i in range(5) for k in range(i+1,5))%2==0]
# Direct products on disjoint supports.
prod=lambda G,H:[tuple(g)+tuple(len(g)+h[i] for i in range(len(h))) for g in G for h in H]
C2=[(0,1),(1,0)];C3=group_from_gens([(1,2,0)])
finite={name:group_audit(G) for name,G in [('D8',D8),('C3_times_D8',prod(C3,D8)),('S4',S4),('A5',A5),('C2_times_A5',prod(C2,A5))]}
print(json.dumps({'method':'Direct adjacency matrices, explicit wreath-product centralizers, exhaustive S8 check of D8 hostile; independent finite permutation groups','marked_checks':checks,'centralizer_constrained_isomorphism_tests':alltests,'D8_hostiles':hostiles,'finite_group_controls':finite,'marked_pair_details':pair_details},indent=2,sort_keys=True))
