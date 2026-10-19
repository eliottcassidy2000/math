#!/usr/bin/env python3
"""Exact independent controls for simultaneous mod-seven two-subdivisions."""
from math import isqrt
import json, random
CHECKS=0

def require(value,msg):
    global CHECKS
    CHECKS+=1
    if not value: raise RuntimeError(msg)

def sq(x): return isqrt(x)**2==x

def construct(n, edges):
    original=[14*i+3 for i in range(n)]
    H=max(original)
    m=len(edges)
    threshold=100*(H+m+1)
    R=threshold+(1-threshold)%14
    labels=original[:]
    target=set()
    for i,(a,b) in enumerate(edges):
        A,B=original[a],original[b]
        r=R+14*i
        v=r*r-B
        s=(A+B+1-r*r)//2
        w=s*s-A
        require((A+B+1-r*r)%2==0,'parity')
        require(v%7==5 and w%7==6,'residue')
        require(H<v<w<5*R**4,'height and disjoint bands')
        vi,wi=len(labels),len(labels)+1
        labels.extend([v,w])
        target.update([tuple(sorted((a,wi))), (vi,wi), tuple(sorted((b,vi)))])
    require(len(labels)==len(set(labels)),'collision')
    for i in range(len(labels)):
        for j in range(i):
            require(sq(labels[i]+labels[j]) == ((j,i) in target),('edge mismatch',n,edges,i,j))
    K=n+2*m
    require(max(labels)<=5*(1514*K)**4,'uniform quartic')
    return max(labels)

counts={}
for n in range(1,6):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for mask in range(1<<len(pairs)):
        construct(n,[e for i,e in enumerate(pairs) if mask>>i&1])
    counts[n]=1<<len(pairs)
complete=[]
for n in range(6,25):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    h=construct(n,pairs)
    complete.append({'n':n,'m':len(pairs),'height':h})
rng=random.Random(20261019)
for z in range(100):
    n=rng.randrange(5,41)
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    rng.shuffle(pairs)
    pairs=pairs[:min(len(pairs),rng.randrange(1,101))]
    pairs=[e if rng.randrange(2) else e[::-1] for e in pairs]
    construct(n,pairs)
# Strong local bridge hostile: arbitrary full seed interval, endpoints can coincide.
bridge_counts=0
for H in list(range(1,21))+[50,100]:
    choices=range(1,H+1) if H<=20 else [1,2,H//2,H-1,H]
    for A in choices:
        for B in choices:
            r=2*H+5
            if r%2 != (A+B+1)%2: r+=1
            v=r*r-B; s=(A+B+1-r*r)//2; w=s*s-A
            require(H<v<w,'bridge bands')
            require(sq(v+w),'middle edge')
            for y in range(1,H+1):
                require(sq(v+y)==(y==B),'v forbidden adjacency')
                require(sq(w+y)==(y==A),'w forbidden adjacency')
            bridge_counts+=1
print(json.dumps({'status':'PASS','always_active_checks':CHECKS,'all_simple_graphs_through_five_vertices':counts,'complete_graphs':complete,'deterministic_random_edge_orders':100,'arbitrary_seed_bridges':bridge_counts},indent=2,sort_keys=True))
