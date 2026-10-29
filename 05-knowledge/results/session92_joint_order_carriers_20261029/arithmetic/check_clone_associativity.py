#!/usr/bin/env python3
"""Exact clone weights, hostile regrouping, and all binary association trees.
No third-party packages. All checks remain active under python -O.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import gcd, prod
checks=0

def check(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(label)

def divs(n):
    return tuple(d for d in range(1,n+1) if n%d==0)

def factor(n):
    out={}
    p=2
    while p*p<=n:
        while n%p==0:
            out[p]=out.get(p,0)+1; n//=p
        p+=1
    if n>1: out[n]=out.get(n,0)+1
    return out

def trees(lo,hi):
    if hi-lo==1:
        return [lo]
    return [(a,b) for k in range(lo+1,hi) for a in trees(lo,k) for b in trees(k,hi)]

def messages(tree, factors):
    if isinstance(tree,int):
        return Counter({d:1 for d in divs(factors[tree])})
    l,r=(messages(t,factors) for t in tree)
    out=Counter()
    for a,ca in l.items():
        for b,cb in r.items():
            out[a*b]+=ca*cb
    return out

def law(tree, factors, total, corrected=True):
    if isinstance(tree,int):
        return {(total,):F(1)}
    l,r=tree; cl=messages(l,factors); cr=messages(r,factors)
    splits=[(a,b) for a in cl for b in cr if a*b==total]
    denom=sum(cl[a]*cr[b] for a,b in splits) if corrected else len(splits)
    out={}
    for a,b in splits:
        p=F(cl[a]*cr[b] if corrected else 1,denom)
        for x,px in law(l,factors,a,corrected).items():
            for y,py in law(r,factors,b,corrected).items():
                out[x+y]=p*px*py
    return out

for A in range(1,41):
    for B in range(1,41):
        ds=list(product(divs(A),divs(B)))
        c=Counter(a*b for a,b in ds)
        w={(a,b):F(a*b,c[a*b]) for a,b in ds}
        check(set(c)==set(divs(A*B)), 'surjectivity')
        for d in c:
            fa,fb,fd=factor(A),factor(B),factor(d)
            formula=prod(min(fa.get(p,0),fd.get(p,0))-max(0,fd.get(p,0)-fb.get(p,0))+1 for p in fa.keys()|fb.keys())
            check(c[d]==formula, 'valuation clone formula')
            check(sum(w[a,b] for a,b in ds if a*b==d)==d, 'native pushforward')
        check(sum(w.values())==sum(divs(A*B)), 'sigma repair')
        check(sum(w.values())-A*B==sum(divs(A*B))-A*B,'proper divisor repair')
        rank_one=all(w[a,b]*w[1,1]==w[a,1]*w[1,b] for a,b in ds)
        check(rank_one==(gcd(A,B)==1), 'independence iff coprime')

families=[(2,2,2),(2,2,2,2),(4,8,2),(6,6,6),(6,10,15),(2,3,5,7),(4,6,9,2),(2,2,2,2,2),(2,2,2,2,2,2)]
num_trees=0
for factors in families:
    raw=list(product(*(divs(n) for n in factors)))
    c=Counter(prod(t) for t in raw)
    for tree in trees(0,len(factors)):
        num_trees+=1
        check(messages(tree,factors)==c,'message counts independent of shape')
        for d,cd in c.items():
            target={t:F(1,cd) for t in raw if prod(t)==d}
            got=law(tree,factors,d)
            check(got==target,'all leaf probabilities preserve sources')
            for t,p in got.items():
                check(p*d==F(d,cd),'weighted native source')

left=law(((0,1),2),(2,2,2),2,False)
right=law((0,(1,2)),(2,2,2),2,False)
check(left=={(1,1,2):F(1,2),(1,2,1):F(1,4),(2,1,1):F(1,4)},'left hostile')
check(right=={(1,1,2):F(1,4),(1,2,1):F(1,4),(2,1,1):F(1,2)},'right hostile')
check(left!=right,'association failure exposed')
D=divs(6); c=Counter(a*b for a in D for b in D)
w={(a,b):F(a*b,c[a*b]) for a in D for b in D}
x=sum(v for (a,b),v in w.items() if a<6 and b<6)
y=sum(v for (a,b),v in w.items() if a==6 and b<6)
z=sum(v for (a,b),v in w.items() if a<6 and b==6)
check((x,y,z)==(22,F(33,2),F(33,2)),'repaired flags')
check(y*z/((x+y)*(x+z))==F(9,49),'repaired flag angle')
check(F(9,49)!=F(36,91),'abundancy angle no longer valid')
for k in range(2,8):
    tuples=[tuple(p if i==j else 1 for i in range(k)) for j in range(k) for p in [2]]
    g=[F(t[0]==2) for t in tuples]; h=[F(t[1]==2) for t in tuples]
    eg=sum(g)/k; eh=sum(h)/k
    check(eg==F(1,k) and eh==F(1,k), 'guard means')
    check(sum(a*b for a,b in zip(g,h))==0 and eg*eh==F(1,k*k), 'joint guard hostile')
    check(sum(a*a for a in g)/k-eg*eg==F(k-1,k*k), 'conditional variance')
print(f'PASS {checks} active exact checks; {num_trees} association trees; 1600 factor pairs')
print('6x6: naive total 144, native/repaired 91; proper 108 versus 55; repaired flag angle 9/49 != 36/91')
print('p,p,p total p: naive left (1/4,1/4,1/2), right (1/2,1/4,1/4), repaired (1/3,1/3,1/3)')
