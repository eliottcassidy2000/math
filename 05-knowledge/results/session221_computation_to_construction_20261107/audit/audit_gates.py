#!/usr/bin/env python3
"""Independent literal gate and continuant audit; no producer imports."""
from itertools import product
from functools import lru_cache
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def data(out,root):
    n=len(out);full=(1<<n)-1
    reverse=tuple(sum(1<<v for v in range(n) if out[v]>>u&1) for u in range(n))
    @lru_cache(None)
    def count(mask,last,rev):
        rest=mask^(1<<last)
        if not rest:return 1
        incoming=out[last] if rev else reverse[last]
        bits=rest&incoming;total=0
        while bits:
            b=bits&-bits;bits-=b;total+=count(rest,b.bit_length()-1,rev)
        return total
    @lru_cache(None)
    def h(mask):
        if not mask:return 1
        return sum(count(mask,v,False) for v in range(n) if mask>>v&1)
    R=sum(count(mask,root,True)*h(full^mask) for mask in range(1,1<<n) if mask>>root&1)
    return h(full),R,count(full,root,False)

def core(n,code):
    out=[0]*n;t=0
    for i in range(n):
        for j in range(i+1,n):
            u,v=(i,j) if code>>t&1 else (j,i)
            out[u]|=1<<v;t+=1
    return tuple(out)

def gate(out,root,m):
    n=len(out);N=n+m+2;b=N-2;c=N-1
    # Specify every pair afresh, separate from the producer's sink/J helpers.
    new=[0]*N
    for u in range(N):
        for v in range(u+1,N):
            if v<n:forward=bool(out[u]>>v&1)
            else:forward=not(v==c and u==root)
            if forward:new[u]|=1<<v
            else:new[v]|=1<<u
    return tuple(new),b

def physical(ms):
    graph=(0,);root=0
    for m in ms:graph,root=gate(graph,root,m)
    return graph,root

def main():
    tested=0
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            old=core(n,code)
            for root in range(n):
                h,R,e=data(old,root)
                for m in range(4):
                    out,r=gate(old,root,m);hh,RR,_=data(out,r);q=2**(m+1)
                    expected=(h+q*R,2*h+q*R-(e if m==0 else 0))
                    if (hh,RR)!=expected:raise ValueError(('gate',old,root,m,hh,RR,expected))
                    tested+=1
    programs=0
    for k in range(4):
        for ms in product((1,2,3),repeat=k):
            if 1+sum(m+2 for m in ms)>13:continue
            h=R=1
            for m in ms:h,R=h+2**(m+1)*R,2*h+2**(m+1)*R
            out,r=physical(ms)
            if data(out,r)[:2]!=(h,R):raise ValueError(('program',ms,h,R))
            programs+=1
    seenpairs={};hist={};stack=[(1,1,())];cap=10**6;words=0;collision=None
    while stack:
        h,R,ms=stack.pop();words+=1
        if (h,R) in seenpairs:raise ValueError('pair collision')
        seenpairs[h,R]=ms
        if h in hist and (collision is None or h<collision['h']):collision={'h':h,'programs':[list(hist[h]),list(ms)]}
        hist.setdefault(h,ms)
        if h%4!=1 or h%8!=(1 if not ms or ms[0]>1 else 5):raise ValueError('congruence')
        q=4;m=1
        while h+q*R<=cap:
            stack.append((h+q*R,2*h+q*R,ms+(m,)));q*=2;m+=1
    if words!=4047 or len(hist)!=2737 or collision['h']!=89:raise ValueError(('census',words,len(hist),collision))
    literal=[]
    for ms in ((2,2),(3,1)):
        out,r=physical(ms);h,R,e=data(out,r)
        after,rr=gate(out,r,1);hh=data(after,rr)[0]
        literal.append({'program':list(ms),'h':h,'R':R,'next_h':hh})
    if [(d['h'],d['R'],d['next_h']) for d in literal]!=[(89,98,481),(89,106,513)]:raise ValueError('collision')
    result={'status':'PASS','independent_literal_macro_checks':tested,'independent_literal_programs':programs,'complete_program_count_below_1000000':words,'distinct_outputs':len(hist),'least_collision':89,'pair_collisions':0,'physical_collision_controls':literal,'scope':'Inherited O^m J^2 macro; audit supports exact reduced-language recurrence, pair uniqueness and scoped sparsity, not all-target acceptance.'}
    (ROOT/'gates_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
