#!/usr/bin/env python3
"""Actual fixed-source cubes, Gray-cycle port obstruction, and buffer lifts."""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations, combinations, product
import json

A = ((0,0,0,0,1,1,1,0),(1,0,0,0,0,1,1,0),
     (1,1,0,0,0,0,1,0),(1,1,1,0,0,0,0,0),
     (0,1,1,1,0,0,0,0),(0,0,1,1,1,0,0,0),
     (0,0,0,1,1,1,0,1),(1,1,1,1,1,1,0,0))
P = ((3,0,6,7,1,5,4,2),(3,1,6,7,0,5,4,2))
G = (4,0,2,3,1,5,6,7)
PINS = (0,2,3,5,6,7)
CHECKS = 0

def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)

def source(r):
    return tuple(tuple(0 if u==v else A[u%8][v%8] if u//8==v//8
                       else int(u//8<v//8) for v in range(8*r)) for u in range(8*r))

def path(r, word):
    return tuple(8*i+u for i in range(r) for u in P[(word>>i)&1])

def partner(B, p):
    # Independent unordered-edge implementation; transport tests use all ordered pairs.
    E = {frozenset((u,v)) for u,v in zip(p,p[1:])}
    return tuple(tuple(0 if u==v else B[u][v] if frozenset((u,v)) in E
                       else B[v][u] for v in range(len(B))) for u in range(len(B)))

def route(r, x, y):
    out=[x]
    for i in range(r):
        if (x^y)&(1<<i):
            x ^= 1<<i
            out.append(x)
    return out

def core():
    paths=[]; autos=[]
    for p in permutations(range(8)):
        if all(A[u][v] for u,v in zip(p,p[1:])):
            paths.append(p)
        if all(A[u][v]==A[p[u]][p[v]] for u,v in product(range(8),repeat=2)):
            autos.append(p)
    need(len(paths)==409 and autos==[tuple(range(8))], 'actual core census')
    need([p for p in paths if all(p[j]==P[0][j] for j in PINS)]==list(P), 'complete pinned fibre')
    for start in range(8):
        seen={start}
        while True:
            new=seen|{v for u in seen for v in range(8) if A[u][v]}
            if new==seen: break
            seen=new
        need(len(seen)==8,'core strongly connected: ordinal blocks are intrinsic SCCs')
    return paths

def cube(r):
    B=source(r); N=1<<r
    paths=[path(r,x) for x in range(N)]
    partners=[partner(B,p) for p in paths]
    free=[(i,j) for i in range(8*r) for j in range(i+2,8*r)]
    words={tuple(B[p[i]][p[j]] for i,j in free) for p in paths}
    converse={tuple(B[p[8*r-1-j]][p[8*r-1-i]] for i,j in free) for p in paths}
    inverse={tuple(1-b for b in w) for w in words}
    inverse_converse={tuple(1-b for b in w) for w in converse}
    need(len(words)==N and len(converse)==N,'actual ordinary and converse source ports inject')
    need(len(words|converse|inverse|inverse_converse)==4*N,'full source/converse/inverse endpoint reservoir disjoint')
    edges=0
    for x,p in enumerate(paths):
        need(len(set(p))==8*r and all(B[u][v] for u,v in zip(p,p[1:])), 'literal Hamiltonian path')
        need(all(p[8*i+j]==8*i+P[0][j] for i in range(r) for j in PINS), 'actual pins')
        need(partner(tuple(zip(*B)),p[::-1])==tuple(zip(*partners[x])), 'converse and reversed integer clocks')
        for i in range(r):
            y=x^(1<<i)
            if x>y: continue
            g=tuple(8*i+G[u%8] if u//8==i else u for u in range(8*r))
            need(all(partners[x][u][v]==partners[y][g[u]][g[v]] for u,v in product(range(8*r),repeat=2)), 'literal global target isomorphism')
            need([j for j in range(8*r) if p[j]!=paths[y][j]]==[8*i+1,8*i+4], 'two actual changed clocks')
            edges+=1
    # Exact Walsh eigenfunctions of random-coordinate fair resampling.
    for mask in range(N):
        f=lambda x: (-1)**((x&mask).bit_count())
        for x in range(N):
            value=F(f(x),2)+sum((F(f(x^(1<<i)),2*r) for i in range(r)),F(0))
            need(value==F(r-mask.bit_count(),r)*f(x),'exact Walsh spectrum')
    loads=Counter(); counts=Counter()
    for x,y in product(range(N), repeat=2):
        p=route(r,x,y); ell=len(p)-1
        need(ell==(x^y).bit_count() and p[-1]==y,'native shortest integer route')
        for u,v in zip(p,p[1:]):
            loads[u,v]+=ell; counts[u,v]+=1
    need(set(counts.values())=={N//2},'uniform oriented-edge pair count')
    congestion=max(F(2*r*c,N) for c in loads.values())
    need(congestion==F(r*(r+1),2),'exact length-weighted cube congestion')
    return dict(r=r,owners=8*r,states=N,edges=edges,gap=str(F(1,r)),
                all_source_paths=409**r,pinned_fraction=str(F(2,409)**r),
                exact_congestion=str(congestion),uniform_return_congestion_bound=4*r*r)

def assignments(eligible):
    out=[]; current=[]
    def visit(i, used):
        if i==len(eligible):
            out.append(tuple(current));return
        for p in sorted(eligible[i]-used):
            current.append(p);visit(i+1,used|{p});current.pop()
    visit(0,set())
    return out

def gray_case(r):
    m=1<<r; gray=[i^(i>>1) for i in range(m)]
    E=[]
    for i in range(m):
        a,b=gray[i],gray[(i+1)%m]; delta=a^b
        need(delta.bit_count()==1,'cyclic Gray edge')
        fixed=((1<<r)-1)^delta
        # These word pins are literal owner pins at clocks 8*j+1 in native paths.
        candidates={x for x in range(m) if (x&fixed)==(a&fixed)}
        need(candidates=={a,b},'complete eligibility from actual clock pins')
        for x in candidates:
            need(all(path(r,x)[8*j+1]==path(r,a)[8*j+1]
                     for j in range(r) if (fixed>>j)&1),'word pins are actual owner pins')
        E.append(candidates)
    states=assignments(E)
    first=tuple(gray); last=tuple(gray[1:]+gray[:1])
    need(set(states)=={first,last},'exactly two complete injective assignments')
    need(sum(a!=b for a,b in zip(first,last))==m,'sharp exponential changed-slot lower bound')
    # A genuine extra eight-owner block supplies a second layer of actual paths.
    E2=[e|{x+m for x in e} for e in E]
    current=list(first); history=[tuple(current)]
    def move(i,p):
        old=current[i]
        need(p in E2[i] and p not in current,'eligible unused actual port')
        need((old^p).bit_count()==1,'one authenticated native cube event')
        need(sum(u!=v for u,v in zip(path(r+1,old),path(r+1,p)))==2,'two actual owner-clock changes')
        current[i]=p; history.append(tuple(current))
        need(len(set(current))==m,'injectivity at every integer clock')
    move(0,gray[0]+m)
    for i in range(m-1,0,-1): move(i,gray[(i+1)%m])
    move(0,gray[1]+m); move(0,gray[1])
    need(tuple(current)==last and len(history)-1==m+2,'buffer return restores boundary and all slot endpoints')
    return dict(r=r,source_owners=8*r,slots=m,assignments=2,sharp_atomic_slot_budget=m,
                enlarged_source_owners=8*(r+1),buffer_route_events=m+2,changed_owner_clocks=2*(m+2))

def stabilized_route(r, E, start, target):
    """Every E_i is a cube face. Four layers retain all owned slots throughout."""
    m=1<<r; current=list(start); history=[tuple(current)]
    need(len(set(start))==len(start) and len(set(target))==len(target),'endpoint injections')
    def move(i, x):
        old=current[i]
        need((x%m) in E[i] and x not in current,'eligible unoccupied full resource')
        need((x^old).bit_count()==1,'one actual coordinate move')
        need(sum(u!=v for u,v in zip(path(r+2,x),path(r+2,old)))==2,'two changed owner clocks')
        current[i]=x;history.append(tuple(current))
    for i,(a,b) in enumerate(zip(start,target)):
        move(i,a+m)  # 00 -> 01; only this slot occupies the workspace 01.
        for c in route(r,a,b)[1:]:move(i,c+m)
        move(i,b+3*m)  # 01 -> 11, the distinct target-storage layer.
    for i,b in enumerate(target):
        move(i,b+2*m)  # 11 -> 10 -> 00, now the input layer is empty.
        move(i,b)
    expected=sum((a^b).bit_count()+4 for a,b in zip(start,target))
    need(tuple(current)==tuple(target) and len(history)-1==expected,'exact final assignment and integer event cost')
    return len(history)-1

def main():
    core()
    cubes=[cube(r) for r in range(1,7)]
    gray=[gray_case(r) for r in range(2,9)]
    return dict(status='FINITE-EXACT PASS; arbitrary r proved separately',cubes=cubes,gray=gray,checks=CHECKS)

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
