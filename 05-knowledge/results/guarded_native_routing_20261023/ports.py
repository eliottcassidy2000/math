#!/usr/bin/env python3
"""Exact distinct-terminal coupling and native permutation schedules."""
from collections import deque
from fractions import Fraction as F
from itertools import permutations
import json

CHECKS=0

def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)

def native_graph(r):
    words=[x for x in range(1<<r) if x.bit_count()%3==0]
    pos={w:i for i,w in enumerate(words)};adj=[]
    for w in words:
        neighbours=set()
        for j in range(r-1):
            if ((w>>j)&3) in (1,2):neighbours.add(pos[w^(3<<j)])
        for j in range(r-2):
            mask=7<<j
            if w&mask in (0,mask):neighbours.add(pos[w^mask])
        adj.append(sorted(neighbours))
    paths={}
    for x in range(len(words)):
        prev={x:None};todo=deque([x])
        while todo:
            u=todo.popleft()
            for v in adj[u]:
                if v not in prev:prev[v]=u;todo.append(v)
        need(len(prev)==len(words),'complete native move connectivity')
        for y in range(len(words)):
            p=[];v=y
            while v is not None:p.append(v);v=prev[v]
            paths[x,y]=p[::-1]
    return words,adj,paths

def matching(M):
    n=len(M);owner={}
    def insert(row,seen):
        for col in range(n):
            if M[row][col]>0 and col not in seen:
                seen.add(col)
                if col not in owner or insert(owner[col],seen):
                    owner[col]=row;return True
        return False
    for row in range(n):
        if not insert(row,set()):return None
    sigma=[None]*n
    for col,row in owner.items():sigma[row]=col
    return sigma

def decompose(R):
    M=[row[:] for row in R];n=len(M);terms=[]
    need(all(sum(row)==1 for row in R),'unit row demands')
    need(all(sum(R[i][j] for i in range(n))==1 for j in range(n)),'actual unit terminal capacities')
    while any(any(row) for row in M):
        sigma=matching(M);need(sigma is not None,'positive-support perfect matching')
        epsilon=min(M[i][sigma[i]] for i in range(n))
        need(epsilon>0,'positive coefficient')
        terms.append((epsilon,sigma))
        for i in range(n):M[i][sigma[i]]-=epsilon
        t=sum(M[0])
        need(all(sum(row)==t for row in M),'residual row balance')
        need(all(sum(M[i][j] for i in range(n))==t for j in range(n)),'residual terminal capacities')
    need(len(terms)<=n*n and sum(w for w,p in terms)==1,'finite exact convex decomposition')
    for i in range(n):
        for j in range(n):
            need(sum(w for w,p in terms if p[i]==j)==R[i][j],'reconstruct every actual return demand')
    return terms

def schedule(sigma,adj,paths,D):
    n=len(sigma);at=list(range(n));wanted=[None]*n
    for token,dest in enumerate(sigma):wanted[dest]=token
    slots=0;pairs=0
    for dest in range(n):
        current=at.index(wanted[dest])
        if current==dest:continue
        before=at[:];p=paths[dest,current]
        edges=list(zip(p,p[1:]));moves=edges+list(reversed(edges[:-1]))
        for u,v in moves:
            need(v in adj[u] and u in adj[v],'both opposite moves authenticated')
            at[u],at[v]=at[v],at[u];slots+=1
            need(sorted(at)==list(range(n)),'all actual terminal owners distinct at every slot')
        before[dest],before[current]=before[current],before[dest]
        need(at==before,'interior tokens restored after endpoint exchange')
        pairs+=1
    need(at==wanted,'entire terminal permutation retained')
    need(pairs<=n-1 and slots<=(n-1)*(2*D-1),'explicit integer-clock compilation bound')
    return slots

def main():
    reports=[]
    for r in (4,6):
        words,adj,paths=native_graph(r);n=len(words);D=4*r
        need(max(len(p)-1 for p in paths.values())<=D,'native diameter bound')
        pp=list(permutations(range(n))) if r==4 else [tuple(range(n)),tuple(reversed(range(n))),tuple((i+1)%n for i in range(n))]
        worst=max(schedule(p,adj,paths,D) for p in pp)
        R=[[F(i==j,3)+F(j==(i+1)%n,3)+F(1,3*n) for j in range(n)] for i in range(n)]
        need(any(R[i][j]!=R[j][i] for i in range(n) for j in range(n)),'nonreversible positive control')
        terms=decompose(R)
        compiled=max(schedule(p,adj,paths,D) for w,p in terms)
        reports.append(dict(r=r,states=n,permutations_checked=len(pp),maximum_slots=worst,
                            return_decomposition_terms=len(terms),maximum_return_schedule_slots=compiled,
                            certified_slot_bound=(n-1)*(2*D-1)))
    need(matching([[F(1),F(0)],[F(1),F(0)]]) is None,'scarce actual target defeats distinct coupling')
    eligible=({0,1},{1,2},{0,2})
    choices=[p for p in permutations(range(3)) if all(p[i] in eligible[i] for i in range(3))]
    need(len(choices)==2,'incoming heterogeneous-slot triangle has exactly two assignments')
    need(sum(a!=b for a,b in zip(*choices))==3,'no legal pair-only schedule for heterogeneous eligibility')
    return dict(status='FINITE-EXACT PASS; unit-capacity complete-fibre permutations only',
                checks=CHECKS,cases=reports)

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
