#!/usr/bin/env python3
"""Exact boundary cube and prescribed unit geometry; no removable assertions."""
from collections import Counter,defaultdict
from fractions import Fraction as Q
from itertools import combinations,product
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
def need(test,message):
    if not test:
        raise RuntimeError(message)
def chi(s,x):
    return (-1)**((s&x).bit_count())
def fourier(w):
    return tuple(sum(w[x]*chi(s,x) for x in range(8)) for s in range(8))
def margins(w):
    return tuple(sum(w[x] for x in range(8)
                     if (x&mask)==value)
                 for mask in (3,5,6)
                 for value in range(8) if value&mask==value)
def support_margins(w):
    return tuple(int(v>0) for v in margins(w))

relations=[tuple((mask>>x)&1 for x in range(8)) for mask in range(256)]
tables=[fourier(w) for w in relations]
full_hist=Counter()
numerical_classes=defaultdict(list)
for mask,w in enumerate(relations):
    missing=[x for x in range(8) if not w[x]]
    independent=all((x^y).bit_count()!=1 for x,y in combinations(missing,2))
    full=all(support_margins(w))
    need(full==independent,'pair shadow/independence equivalence')
    if full:
        full_hist[sum(w)]+=1
    numerical_classes[margins(w)].append(mask)
need(dict(full_hist)=={4:2,5:8,6:16,7:8,8:1},'wrong35relation spectrum')
ambiguous=[v for v in numerical_classes.values() if len(v)>1]
parity_masks=[sum(1<<x for x in range(8) if x.bit_count()%2==p) for p in (0,1)]
need(len(ambiguous)==1 and set(ambiguous[0])==set(parity_masks),'Boolean marginal ambiguity')

gluing_checks=0
for a,b in product(range(256),repeat=2):
    w=relations[a&b]
    actual=tables[a&b]
    for s in range(8):
        convolution=sum(tables[a][t]*tables[b][s^t] for t in range(8))
        need(convolution==8*actual[s],'Fourier gluing convolution')
    need(sum(w)*8==sum(tables[a][s]*tables[b][s] for s in range(8)),
         'Parseval total count')
    gluing_checks+=1

weighted_classes=defaultdict(list)
for w in product(range(3),repeat=8):
    weighted_classes[margins(w)].append(w)
weighted_hist=Counter()
for group in weighted_classes.values():
    first=group[0]
    even=[first[x] for x in range(8) if chi(7,x)==1]
    odd=[first[x] for x in range(8) if chi(7,x)==-1]
    lo=max(-min(even),max(odd)-2)
    hi=min(min(odd),2-max(even))
    expected={tuple(first[x]+lam*chi(7,x) for x in range(8)) for lam in range(lo,hi+1)}
    need(set(group)==expected,'exact integer marginal fiber interval')
    weighted_hist[len(group)]+=1
need(sum(k*v for k,v in weighted_hist.items())==3**8,'weighted universe count')

# Exact rational coordinates in the basis 1,alpha, with norm a^2+ab+b^2.
def norm(v):
    a,b=v
    return a*a+a*b+b*b
def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))
def add(a,b):
    return tuple(x+y for x,y in zip(a,b))
def unit_direction(t):
    d=1+t+t*t
    v=(Q(1-t*t,d),Q(t*(t+2),d))
    need(norm(v)==1,'unit direction parameterization')
    return v
def all_edges(points):
    return {(i,j) for i,j in combinations(range(len(points)),2) if norm(sub(points[i],points[j]))==1}

base=[(Q(0),Q(0)),(Q(1),Q(0)),(Q(0),Q(1)),
      (Q(0),Q(-1)),(Q(2),Q(0)),(Q(0),Q(2))]
base_edges={(0,1),(0,2),(1,2),(0,3),(1,4),(2,5)}
need(all_edges(base)==base_edges,'inherited induced triangle/spoke geometry')

def attach_leaf(points,edges,fixed,parent,color):
    t=0
    while True:
        candidate=add(points[parent],unit_direction(t))
        if candidate not in points and all(j==parent or norm(sub(candidate,p))!=1
                                           for j,p in enumerate(points)):
            break
        t+=1
    idx=len(points)
    points.append(candidate)
    edges.add((min(parent,idx),max(parent,idx)))
    fixed[idx]=color
    return t

def make_graph(kind):
    if kind=='F':
        points=list(base); edges=set(base_edges); ports=(3,4,5)
        forbids=[(j,c) for j in range(3) for c in (3,4,5)]
    else:
        points=list(base[3:]);edges=set();ports=(0,1,2);forbids=[]
    fixed={};params=[]
    for port in ports:
        forbids.extend((port,c) for c in (2,3,4,5))
    for parent,color in forbids:
        params.append(attach_leaf(points,edges,fixed,parent,color))
    center=len(points)
    points.append((Q(100),Q(0)))
    for color in ((2,3,4,5) if kind=='F' else (3,4,5)):
        params.append(attach_leaf(points,edges,fixed,center,color))
    need(all_edges(points)==edges,'extra/missing unit edge in '+kind)
    need(len(set(points))==len(points),'coincident vertices in '+kind)
    return points,edges,fixed,ports,params

def extensions(n,edges,fixed):
    adjacency=[set() for _ in range(n)]
    for a,b in edges:
        adjacency[a].add(b);adjacency[b].add(a)
    colors=dict(fixed)
    if any(a in colors and b in colors and colors[a]==colors[b] for a,b in edges):
        return 0
    def search():
        if len(colors)==n:
            return 1
        menus=[(tuple(c for c in range(6) if c not in {colors[j] for j in adjacency[i] if j in colors}),i)
               for i in range(n) if i not in colors]
        menu,i=min(menus,key=lambda item:len(item[0]))
        total=0
        for c in menu:
            colors[i]=c
            total+=search()
            del colors[i]
        return total
    return search()

geometry={}
weights={}
for kind in ('F','G'):
    points,edges,fixed,ports,params=make_graph(kind)
    w=[]
    for x in range(8):
        pin=dict(fixed)
        pin.update({ports[j]:(x>>j)&1 for j in range(3)})
        w.append(extensions(len(points),edges,pin))
    weights[kind]=tuple(w)
    geometry[kind]=dict(vertices=len(points),edge_count=len(edges),
                         points=[[str(a),str(b)] for a,b in points],
                         edges=sorted(edges),fixed_colors=fixed,ports=ports,
                         direction_parameters=params,weights=w,walsh=fourier(w))
need((geometry['F']['vertices'],geometry['F']['edge_count'])==(32,31),'F sizes')
need((geometry['G']['vertices'],geometry['G']['edge_count'])==(19,15),'G sizes')
need(weights['F']==tuple(0 if x in (0,7) else 4 for x in range(8)),'exact geometric NAE count')
need(weights['G']==(3,)*8,'exact geometric free count')
for kind,w in weights.items():
    W=fourier(w)
    need(W[0]==24 and all(W[s]==0 for s in (1,2,4,7)),'coarse mode agreement')
    need(all(support_margins(w)),'full existential pair shadows')
need(fourier(weights['F'])==tuple(24 if s==0 else -8 if s in (3,5,6) else 0
                                for s in range(8)),'NAE pair-mode spectrum')
need(margins(weights['F'])!=margins(weights['G']),'counting marginals falsely equated')
need(weights['F'][0]==0 and weights['G'][0]==3,'pinned continuation hostile')

M=[[-2 if x==y else 2 if (x^y).bit_count()==1 else 0 for y in range(8)] for x in range(8)]
for s in range(8):
    v=[chi(s,x) for x in range(8)]
    need([sum(M[x][y]*v[y] for y in range(8)) for x in range(8)]
         ==[4*(1-s.bit_count())*z for z in v],'K3 spectrum')
for s in (1,2,4):
    plus=[1+chi(s,x) for x in range(8)]
    minus=[1-chi(s,x) for x in range(8)]
    need(plus!=minus,'distinct singleton tables')
    need([sum(M[x][y]*plus[y] for y in range(8)) for x in range(8)]
         ==[sum(M[x][y]*minus[y] for y in range(8)) for x in range(8)],'K3 nullspace hostile')

receipt=dict(status='FINITE-EXACT',all_cube_relations=256,
             full_pair_shadow_histogram=dict(sorted(full_hist.items())),
             numerical_marginal_collision=parity_masks,
             exact_gluing_pairs=gluing_checks,weighted_tables=3**8,
             weighted_fiber_size_histogram=dict(sorted(weighted_hist.items())),
             unit_graphs=geometry,
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
out=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
(HERE/'receipt.json').write_text(out,encoding='utf-8')
print(json.dumps(dict(status='PASS',relations=256,full_pair_shadows=35,
                     gluing_pairs=gluing_checks,weighted_tables=3**8,
                     F=(32,31),G=(19,15),receipt_sha256=hashlib.sha256(out.encode()).hexdigest()),sort_keys=True))
