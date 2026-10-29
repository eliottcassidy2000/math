#!/usr/bin/env python3
"""Actual conjugated marking triangles, with a linear edge budget."""
from itertools import product
import random
import json

CHECKS=0
def ck(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(msg)
def comp(a,b):return tuple(a[x] for x in b)
def inv(a):
    b=[0]*len(a)
    for i,j in enumerate(a):b[j]=i
    return tuple(b)
def identity(n):return tuple(range(n))
def neg(k):return tuple(sum((-(v//3**i)%3)*3**i for i in range(k)) for v in range(3**k))
def root(a):
    out=[0]*len(a);seen=set()
    for x in range(len(a)):
        if x in seen:continue
        cyc=[];y=x
        while y not in seen:seen.add(y);cyc.append(y);y=a[y]
        ck(len(cyc)%2==1,'odd root')
        for i,v in enumerate(cyc):out[v]=cyc[(i+(len(cyc)+1)//2)%len(cyc)]
    return tuple(out)
def conjugate(a,b):return comp(a,comp(b,inv(a)))
def local(k,path,r=1):
    n=3**k;d=len(path);prefix=sum(v*3**(d-1-i) for i,v in enumerate(path));tail=3**(k-d-1)
    return tuple(x+(((x//tail)%3+r)%3-(x//tail)%3)*tail if x//(3*tail)==prefix else x for x in range(n))
def coordinates(g,path=()):
    if len(g)==1:return []
    s=len(g)//3;r=g[0]//s;out=[(path,r)] if r else []
    for b in range(3):out+=coordinates(tuple(g[b*s+j]%s for j in range(s)),path+(b,))
    return out
def paired(k,a):
    s=3**(k-1);t=neg(k-1)
    return identity(s)+tuple(s+x for x in a)+tuple(2*s+x for x in conjugate(t,a))
def embed0(k,a):
    s=3**(k-1)
    return tuple(a)+tuple(range(s,3*s))
def canonical_triangle(k,depth):
    alpha=neg(k);x=local(k,());y=local(k,(0,)*(depth+1))
    t=comp(comp(x,x),alpha);u=comp(comp(y,y),alpha)
    return t,u
def holonomy(k,t,u):
    alpha=neg(k)
    return comp(inv(root(comp(u,alpha))),comp(root(comp(u,t)),root(comp(t,alpha))))
def transport_prefix(k,path):
    a=identity(3**k)
    for d,r in enumerate(path):
        if r:a=comp(a,local(k,(0,)*d,r))
    return a
def repair(k,c):
    if k<=1:return []
    s=3**(k-1);alpha=neg(k);child=c[:s];b=tuple(c[s+j]-s for j in range(s))
    triangles=[]
    for t,u in repair(k-1,child):
        # Both maps are anti on child0; outside it retain alpha, not the identity.
        triangles.append((tuple(t)+alpha[s:],tuple(u)+alpha[s:]))
    for path,r in coordinates(b):
        t,u=canonical_triangle(k,len(path));a=paired(k,transport_prefix(k-1,path))
        t,u=conjugate(a,t),conjugate(a,u)
        if r==2:t,u=u,t
        triangles.append((t,u))
    return triangles
def aut(k):
    if not k:return [(0,)]
    s=3**(k-1);A=aut(k-1)
    return [tuple(((b+r)%3)*s+abc[b][j] for b in range(3) for j in range(s))
            for r in range(3) for abc in product(A,repeat=3)]
def cent(k):
    if not k:return [(0,)]
    return [comp(embed0(k,c),paired(k,a)) for c in cent(k-1) for a in aut(k-1)]
def random_aut(k,rng):
    if not k:return (0,)
    s=3**(k-1);abc=[random_aut(k-1,rng) for _ in range(3)];r=rng.randrange(3)
    return tuple(((b+r)%3)*s+abc[b][j] for b in range(3) for j in range(s))
def random_cent(k,rng):
    if not k:return (0,)
    return comp(embed0(k,random_cent(k-1,rng)),paired(k,random_aut(k-1,rng)))
def forward(k,x,y):
    for d in range(k-1,-1,-1):
        a,b=x//3**d%3,y//3**d%3
        if a!=b:return (b-a)%3==1
    return False
def main():
    rows=[];rng=random.Random(270243)
    for k in range(2,6):
        n=3**k;alpha=neg(k)
        # Validate all canonical triangle types against the actual tournament arcs.
        for depth in range(k-1):
            t,u=canonical_triangle(k,depth)
            for mark in [t,u]:
                ck(comp(mark,mark)==identity(n),'canonical mark involution')
                for x in range(n):
                    for y in range(x):ck(forward(k,mark[x],mark[y])!=forward(k,x,y),'anti-isomorphism')
        cases=cent(k) if k<=3 else [random_cent(k,rng) for _ in range(24)]
        bound=(n-2*k-1)//4;maximum=0
        for c in cases:
            triangles=repair(k,c);maximum=max(maximum,len(triangles));actual=identity(n)
            ck(len(triangles)<=bound,'linear triangle count')
            for t,u in triangles:
                ck(comp(t,t)==identity(n) and comp(u,u)==identity(n),'transport endpoints involutory')
                h=holonomy(k,t,u)
                ck(comp(h,alpha)==comp(alpha,h),'based loop centralizes source mark')
                actual=comp(actual,h)
            ck(actual==c,'full port permutation restored')
        rows.append(dict(k=k,vertices=n,universe='complete' if k<=3 else '24 seeded controls',
                         tested=len(cases),maximum_triangles=maximum,triangle_bound=bound,edge_bound=3*bound))
    print(json.dumps(dict(rows=rows,active_checks=CHECKS),indent=2,sort_keys=True))
if __name__=='__main__':main()
