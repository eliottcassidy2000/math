#!/usr/bin/env python3
"""Portable exact graph-birth checks. No external dependencies or assert gates."""
from math import isqrt
from itertools import combinations
from collections import Counter
import json
from pathlib import Path

def require(ok,msg):
    if not ok: raise ValueError(msg)
def graph(n):
    require(type(n) is int and n>=1,'positive integer n required')
    return {a:tuple(b for b in range(1,n+1) if a!=b and isqrt(a+b)**2==a+b) for a in range(1,n+1)}
def canon(p): return min(tuple(p),tuple(reversed(p)))
def paths(G):
    out=set();leaves=[v for v in G if len(G[v])==1]
    starts=[min(leaves)] if leaves else sorted(G)
    def rec(p,left):
        if not left:out.add(canon(p));return
        for v in G[p[-1]]:
            if v in left:rec(p+[v],left-{v})
    for s in starts:rec([s],set(G)-{s})
    return sorted(out)
def cover_key(pieces):return tuple(sorted(canon(p) for p in pieces if p))
def validate_cover(G,pieces):
    flat=[v for p in pieces for v in p]
    require(len(flat)==len(set(flat)) and set(flat)==set(G),'not a vertex partition')
    require(all(p for p in pieces),'empty path')
    require(all(b in G[a] for p in pieces for a,b in zip(p,p[1:])),'nonedge')
def endpoints(p):return set((p[0],p[-1]))
def edge_covers(G,k):
    """Independent spanning linear-forest edge search with rollback DSU."""
    V=sorted(G);E=sorted((a,b) for a in G for b in G[a] if a<b)
    target=len(V)-k;deg={v:0 for v in V};parent={v:v for v in V};size={v:1 for v in V};selected=[];out=[]
    def root(a):
        while parent[a]!=a:a=parent[a]
        return a
    def emit():
        A={v:[] for v in V}
        for a,b in selected:A[a].append(b);A[b].append(a)
        seen=set();pieces=[]
        for s in V:
            if s in seen or len(A[s])==2:continue
            p=[];last=None;v=s
            while True:
                p.append(v);seen.add(v);nxt=[w for w in A[v] if w!=last]
                if not nxt:break
                last,v=v,nxt[0]
            pieces.append(p)
        require(len(seen)==len(V),'cycle escaped DSU');out.append(cover_key(pieces))
    def rec(i):
        need=target-len(selected)
        if need==0:emit();return
        if len(E)-i<need:return
        if sum(2-deg[v] for v in V)<2*need:return
        a,b=E[i];ra,rb=root(a),root(b)
        if deg[a]<2 and deg[b]<2 and ra!=rb:
            if size[ra]<size[rb]:ra,rb=rb,ra
            parent[rb]=ra;size[ra]+=size[rb];deg[a]+=1;deg[b]+=1;selected.append((a,b))
            rec(i+1)
            selected.pop();deg[a]-=1;deg[b]-=1;size[ra]-=size[rb];parent[rb]=rb
        rec(i+1)
    rec(0)
    require(len(out)==len(set(out)),'duplicate cover');return sorted(out)
def birth(n,old,new):
    N=set(graph(n)[n]);extensions=[];internal=[];fibres=Counter(); endpoint_counts=Counter()
    for p in old:
        for e in endpoints(p):endpoint_counts[e]+=1
    for p in new:
        j=p.index(n)
        pieces=cover_key((p[:j],p[j+1:]));validate_cover(graph(n-1),pieces)
        if j in (0,len(p)-1):extensions.append({'path':p,'old_path':pieces[0],'port':p[1] if j==0 else p[-2]})
        else:internal.append({'path':p,'cover':pieces,'ports':sorted((p[j-1],p[j+1]))});fibres[pieces]+=1
    require(len(extensions)==sum(endpoint_counts[a] for a in N),'extension identity')
    return {'N':n,'neighbors':sorted(N),'old_count':len(old),'new_count':len(new),'old_endpoint_counts':dict(sorted(endpoint_counts.items())), 'extensions':extensions,'internal':internal,'internal_count':len(internal),'extension_count':len(extensions),'distinct_internal_deletion_covers':len(fibres),'deletion_multiplicity_histogram':dict(sorted(Counter(fibres.values()).items()))}
def main():
    data=[];P={n:paths(graph(n)) for n in range(24,29)}
    require([len(P[n]) for n in range(24,29)]==[0,10,12,35,52],'known counts')
    for n in range(25,29):
        row=birth(n,P[n-1],P[n])
        covers=edge_covers(graph(n-1),2)
        good=[c for c in covers if responses(c,row['neighbors'])[1]]
        deleted={cover_key(x['cover']) for x in row['internal']}
        require(set(good)==deleted,'independent cover census mismatch')
        row['old_two_cover_count']=len(covers)
        row['independent_compatible_cover_count']=len(good)
        data.append(row)
    print(json.dumps(data,sort_keys=True,indent=2))

# An endpoint state is an unordered collection of unordered pairs.
# Pair (a,a) records a singleton path, not two choices of the same port.
def endpoint_state(cover):return tuple(sorted(tuple(sorted((p[0],p[-1]))) for p in cover))
def state_key(pairs):return tuple(sorted(tuple(sorted(p)) for p in pairs))
def state_birth(state,v,neighbors):
    require(all(v not in p for p in state),'new label already exposed')
    N=set(neighbors);out=Counter();out[state_key(state+((v,v),))]+=1
    for i,p in enumerate(state):
        for a in sorted(set(p)&N):
            b=p[1] if p[0]==a else p[0]
            out[state_key(state[:i]+state[i+1:]+((b,v),))]+=1
    for i,j in combinations(range(len(state)),2):
        for a in sorted(set(state[i])&N):
            for b in sorted(set(state[j])&N):
                aa=state[i][1] if state[i][0]==a else state[i][0]
                bb=state[j][1] if state[j][0]==b else state[j][0]
                rest=tuple(p for z,p in enumerate(state) if z not in (i,j))
                out[state_key(rest+((aa,bb),))]+=1
    return out

def state_transfer(distribution,v,neighbors):
    out=Counter()
    for s,c in distribution.items():
        for t,m in state_birth(s,v,neighbors).items():out[t]+=c*m
    return out

def responses(cover,neighbors):
    E=[endpoints(p)&set(neighbors) for p in cover]
    return sum(map(len,E)),sum(len(E[i])*len(E[j]) for i,j in combinations(range(len(E)),2))

if __name__=='__main__':main()
