#!/usr/bin/env python3
"""Bounded exact tests of tight-suffix decomposition and buffer lifts.
The strict-rank subroutine is BFS, not an efficient general implementation.
"""
from collections import deque,defaultdict
from itertools import product
from functools import lru_cache
from random import Random
import json

CHECKS=0;MOVES=0

def require(x,msg):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(msg)

@lru_cache(maxsize=500000)
def rank0(values):
    basis={}
    for v in values:
        while v:
            bit=v.bit_length()-1
            if bit in basis:v^=basis[bit]
            else:basis[bit]=v;break
    return len(basis)

def rank(values):return rank0(tuple(sorted(set(values))))
def g(n):return n*n-n+1 if n else 0
def c(n):return 1 if n==1 else g(n-1)
def union(classes,ids):return [v for i in ids for v in classes[i]]
def vals(classes,state):return [classes[i][a] for i,a in enumerate(state)]

def enumerate_states(classes,base=()):
    br=rank(base);states=[]
    def rec(i,chosen,code):
        if i==len(classes):states.append(tuple(code));return
        for a,v in enumerate(classes[i]):
            if rank(tuple(base)+tuple(chosen)+(v,))==br+i+1:rec(i+1,chosen+[v],code+[a])
    rec(0,[],[]);return states

def strict_path(classes,indices,base,start,target):
    n=len(indices);br=rank(base)
    if start==target:return [start]
    # Check every nonempty subset rather than infer the inherited condition.
    for mask in range(1,1<<n):
        jj=[indices[z] for z in range(n) if mask>>z&1]
        require(rank(tuple(base)+tuple(union(classes,jj)))-br>len(jj),'strict Hall premise')
    def valid(st):return rank(tuple(base)+tuple(classes[i][a] for i,a in zip(indices,st)))==br+n
    require(valid(start) and valid(target),'strict endpoint independence')
    pred={start:None};queue=deque([start])
    while queue and target not in pred:
        st=queue.popleft()
        for j,i in enumerate(indices):
            for a in range(len(classes[i])):
                if a==st[j]:continue
                nxt=st[:j]+(a,)+st[j+1:]
                if nxt not in pred and valid(nxt):pred[nxt]=st;queue.append(nxt)
    require(target in pred,'strict rank connected')
    rev=[];st=target
    while st is not None:rev.append(st);st=pred[st]
    answer=list(reversed(rev));require(len(answer)-1<=g(n),'cited strict bound finite control')
    return answer

def route(classes,source,target,base=()):
    global MOVES
    s=len(classes);br=rank(base)
    require(all(rank(tuple(base)+tuple(A))-br>=s-i for i,A in enumerate(classes)),'rank staircase')
    require(rank(tuple(base)+tuple(vals(classes,source)))==br+s,'source independent')
    require(rank(tuple(base)+tuple(vals(classes,target)))==br+s,'target independent')
    tight=[i for i in range(s) if rank(tuple(base)+tuple(union(classes,range(i,s))))-br==s-i]
    # Exhaustively independently confirm the suffix-only assertion.
    for mask in range(1,1<<s):
        ids=[i for i in range(s) if mask>>i&1]
        if rank(tuple(base)+tuple(union(classes,ids)))-br==len(ids):
            require(ids==list(range(min(ids),s)),'all tight sets suffixes')
    a=tight[0] if tight else s
    starts=tight+[s];blocks=list(zip(starts,starts[1:]));bound=g(a)+sum(c(j-i) for i,j in blocks)
    current=list(source);path=[source]
    def add(nxt):
        old=tuple(current);new=tuple(nxt)
        if old==new:return
        changed=sum(x!=y for x,y in zip(old,new))
        require(changed in (1,2),'move width')
        ov=tuple(base)+tuple(vals(classes,old));nv=tuple(base)+tuple(vals(classes,new))
        require(rank(nv)==br+s,'lifted full independence')
        if changed==2:require(rank(ov+nv)==br+s,'double closure preserved')
        current[:]=new;path.append(new)
    for i,j in reversed(blocks):
        b=j-i;inner=tuple(base)+tuple(union(classes,range(j,s)));ir=rank(inner)
        require(rank(inner+tuple(union(classes,range(i,j))))-ir==b,'block rank')
        require(rank(inner+tuple(classes[i]))-ir==b,'full-rank first color')
        if b==1:
            nxt=list(current);nxt[i]=target[i];add(nxt);continue
        indices=list(range(i+1,j));start=tuple(current[z] for z in indices);end=tuple(target[z] for z in indices)
        pp=strict_path(classes,indices,inner,start,end)
        before=len(path)
        if len(pp)==1:
            nxt=list(current);nxt[i]=target[i];add(nxt)
        for k,partial in enumerate(pp[1:],1):
            nxt=list(current)
            for z,val in zip(indices,partial):nxt[z]=val
            candidates=[target[i]] if k==len(pp)-1 else [current[i]]+list(range(len(classes[i])))
            found=False
            for aa in candidates:
                blockvalues=[classes[i][aa]]+[classes[z][nxt[z]] for z in indices]
                if rank(inner+tuple(blockvalues))==ir+b:nxt[i]=aa;found=True;break
            require(found,'buffer exists');add(nxt)
        require(len(path)-before<=c(b),'block move bound')
        require(current[i:j]==list(target[i:j]),'block exactly installed')
    if a:
        contracted=tuple(base)+tuple(union(classes,range(a,s)))
        pp=strict_path(classes,list(range(a)),contracted,tuple(current[:a]),tuple(target[:a]))
        for partial in pp[1:]:add(list(partial)+current[a:])
    require(tuple(current)==target,'target reached')
    require(len(path)-1<=bound<=g(s),'quadratic bound')
    MOVES+=len(path)-1
    return len(path)-1,bound

def main():
    rng=Random(20261019);systems=routes=0;maxlen=0;graphic_fibers=graphic_routes=0
    # Binary linear matroids with classes of distinct actual colored elements;
    # equal values in different colors are parallel, and zero values are loops.
    for s in range(1,6):
        for trial in range(40):
            dim=s+(trial%2);classes=[]
            for i in range(s):
                need=s-i;pool=list(range(1,1<<dim));rng.shuffle(pool);A=[]
                for v in pool:
                    if rank(A)<need:A.append(v)
                if trial%3==0:A.append(0)
                if trial%5==0:A.append(A[0]) # distinct parallel realizations
                classes.append(tuple(A))
            states=enumerate_states(classes);require(bool(states),'nonempty staircase')
            candidates=[(states[0],states[-1])]+[(rng.choice(states),rng.choice(states)) for _ in range(4)]
            for x,y in candidates:
                got,bound=route(classes,x,y);maxlen=max(maxlen,got);routes+=1
            systems+=1
    # All arithmetic graceful states on labels0..m, all upper-tail groupings,
    # and deterministic endpoint pairs from every nonempty fiber through m=7.
    for m in range(1,8):
        classes=[tuple((1<<a)^(1<<(a+d)) for a in range(m+1-d)) for d in range(1,m+1)]
        allstates=enumerate_states(classes)
        for r in range(1,m+1):
            grouped=defaultdict(list)
            for st in allstates:grouped[st[r:]].append(st[:r])
            for tail,states in grouped.items():
                base=tuple(classes[r+i][a] for i,a in enumerate(tail));free=classes[:r]
                pairs=[(states[0],states[-1])]
                if len(states)>2:pairs.append((states[len(states)//2],states[0]))
                for x,y in pairs:
                    got,bound=route(free,x,y,base)
                    if r>=3:require(bound<=r*r-3*r+3,'graceful uniform bound')
                    graphic_routes+=1;maxlen=max(maxlen,got)
                graphic_fibers+=1
    # The published target-preservation hostile, represented by parallel spokes.
    A=((1,2,4,8),(1,2,8),(2,4),(4,1));states=enumerate_states(A)
    require(len(states)==5,'five-state hostile')
    source=(3,0,0,0);target=(3,1,1,1)
    got,bound=route(A,source,target)
    require(got>=4,'known hostile lower bound retained')
    # The cyclic three-color obstruction fails the staircase and is disconnected.
    hostile=((1,2),(2,4),(4,1))
    require(rank(hostile[0])<3,'canonical hostile outside hypothesis')
    hs=enumerate_states(hostile);require(len(hs)==2 and all(x!=y for x,y in zip(*hs)),'cyclic hostile three-coordinate gap')
    print(json.dumps({'status':'PASS','checks':CHECKS,'validated_moves':MOVES,'binary_systems':systems,'binary_routes':routes,'all_graceful_intervals_m_max':7,'graceful_upper_tail_fibers':graphic_fibers,'graceful_routes':graphic_routes,'maximum_produced_length':maxlen,'target_preservation_hostile_route':got,'scope':'Bounded exact verification; BFS implements the cited strict-rank path interface only for these finite controls and is not a polynomial-time general solver.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
