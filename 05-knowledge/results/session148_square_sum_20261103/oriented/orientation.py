#!/usr/bin/env python3
"""Exact orientation distribution via compatible directed path event inclusion-exclusion."""
from math import isqrt,comb
from itertools import permutations

def graph(n):
    E=[(u,v) for u in range(1,n+1) for v in range(u+1,n+1) if isqrt(u+v)**2==u+v]
    A={i:set() for i in range(1,n+1)}
    for u,v in E:A[u].add(v);A[v].add(u)
    return E,A

def paths(n):
    E,A=graph(n); leaves=[v for v in A if len(A[v])==1]
    if len(leaves)>2:return E,[]
    starts=leaves[:1] if leaves else list(A); out=[]
    def walk(v,used,seq):
        if len(seq)==n:
            if leaves or seq[0]<seq[-1]:out.append(tuple(seq))
            return
        for w in sorted(A[v]-used):walk(w,used|{w},seq+[w])
    for s in starts:walk(s,{s},[s])
    return E,out

def moments(E,paths):
    idx={e:i for i,e in enumerate(E)}; events=[]
    for p in paths:
      for p2 in (p,p[::-1]):
        mask=bits=0
        for a,b in zip(p2,p2[1:]):
          i=idx[tuple(sorted((a,b)))];mask|=1<<i
          if a<b:bits|=1<<i
        events.append((mask,bits))
    B=[0]*(len(events)+1); B[0]=1<<len(E)
    def visit(i,mask,bits,k):
      for j in range(i,len(events)):
        em,eb=events[j]
        if (bits^eb)&mask&em:continue
        nm=mask|em;nb=bits|eb
        B[k+1]+=1<<(len(E)-nm.bit_count())
        visit(j+1,nm,nb,k+1)
    visit(0,0,0,0)
    c=[sum((-1)**(k-j)*comb(k,j)*B[k] for k in range(j,len(B))) for j in range(len(B))]
    if sum(c)!=1<<len(E) or any(x<0 for x in c):raise ValueError('invalid distribution')
    return events,B,c

def brute(n):
    E,p=paths(n); idx={e:i for i,e in enumerate(E)}; ev,_,c=moments(E,p); obs=[0]*len(c)
    for x in range(1<<len(E)):
      h=sum((x&m)==b for m,b in ev);obs[h]+=1
    if obs!=c:raise ValueError('distribution mismatch')
    # Independent permutation enumeration checks each unoriented path.
    perms=[s for s in permutations(range(1,n+1)) if s[0]<s[-1] and all(tuple(sorted(pair)) in E for pair in zip(s,s[1:]))]
    if set(p)!=set(perms):raise ValueError('path mismatch')

if __name__=='__main__':
    for n in range(2,9):brute(n)
    for n in (15,16,17,18,23,24,25,26,27):
        E,p=paths(n); ev,B,c=moments(E,p)
        inc={e:sum(e in {tuple(sorted(t)) for t in zip(path,path[1:])} for path in p) for e in E}
        print('n',n,'m',len(E),'U',len(p),'unused',([e for e in E if not inc[e]] if p else 'no Hamilton paths'), flush=True)
        print('orientation_counts', [(i,x) for i,x in enumerate(c) if x], 'compatible_subset_moments',[(i,x) for i,x in enumerate(B) if x],flush=True)
