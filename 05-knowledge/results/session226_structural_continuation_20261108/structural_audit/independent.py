#!/usr/bin/env python3
"""Independent controls: no imports from either producer verifier."""
import itertools as it
import json

def need(p,why):
    if not p: raise RuntimeError(why)

def bank(n):
    pairs=list(it.combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        m=[[False]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            m[i][j]=bool(mask&(1<<k)); m[j][i]=not m[i][j]
        yield m

def from_rows(r): return [[bool(x>>j&1) for j in range(len(r))] for x in r]

def closure(m):
    q=[row[:] for row in m]
    for i in range(len(m)): q[i][i]=True
    for k in range(len(m)):
        for i in range(len(m)):
            for j in range(len(m)): q[i][j]=q[i][j] or (q[i][k] and q[k][j])
    return q

def graft(a,r,b):
    n=len(a); z=len(b); q=[[False]*(n+z) for _ in range(n+z)]
    for i in range(n):
        for j in range(n): q[i][j]=a[i][j]
    for i in range(z):
        for j in range(z): q[n+i][n+j]=b[i][j]
        q[n+i][r]=True
        for j in range(n):
            if j!=r: q[j][n+i]=True
    return q

def pathcount(m,verts=None,start=None,end=None):
    vs=range(len(m)) if verts is None else verts
    return sum((start is None or p[0]==start) and (end is None or p[-1]==end) and all(m[i][j] for i,j in zip(p,p[1:])) for p in it.permutations(vs))

def profile(m,r):
    n=len(m); nonroot=[i for i in range(n) if i!=r]
    D=pathcount(m); E=pathcount(m,end=r); d=pathcount(m,nonroot); U=0
    # Enumerate the two path vertex sets, then orders in each set independently.
    for mask in range(1<<len(nonroot)):
        P=[v for i,v in enumerate(nonroot) if mask>>i&1]
        R=[r]+[v for v in nonroot if v not in P]
        U+=pathcount(m,P)*pathcount(m,R,start=r)
    return D,E,U,d

def move(s,k):
    D,E,U,d=s
    if k=='O': return D,0,2*U,d
    if k=='Z': return D-E+U,d,2*U-d,d
    return D-E+U,D-E,D+U,D

def literal(m,r,k):
    if k!='O':
        out=graft(m,r,[[False]])
        return out, r if k=='Z' else len(m)
    n=len(m); out=[row+[True] for row in m]+[[False]*(n+1)]
    return out,r

counts={}; gadgets=[m for n in range(1,5) for m in bank(n)]
pairs=0
for n in range(1,4):
    for a in bank(n):
        c=closure(a)
        for r in range(n):
            for b in gadgets:
                got=all(all(row) for row in closure(graft(a,r,b)))
                need(got==(n>=2 and all(c[r])), 'independent strongness criterion')
                pairs+=1
counts['warshall_pairs_all_source_roots']=pairs
C=from_rows([8,1,3,6]); hist={}
for label,rows,expected in [('A',[14,28,24,16,1],[423,605,357,185,261]),('B',[10,28,25,16,1],[315,541,551,217,207])]:
    m=from_rows(rows)
    values=[pathcount(graft(m,r,C)) for r in range(5)]
    need(values==expected, 'permutation H9 count')
    hist[label]=values
counts['H9_permutation_counts']=hist
seeds=set(); literal_checks=0
for n in range(1,5):
    for a in bank(n):
        for r in range(n):
            s=profile(a,r); D,E,U,d=s
            need(0<=E<=min(D,d) and 0<=D-E<=U-d and 1<=D<=U,'subset-bin cone')
            seeds.add(s)
            for k in 'OZJ':
                b,t=literal(a,r,k)
                need(profile(b,t)==move(s,k),'literal operation')
                literal_checks+=1
counts['subset_bin_literal_operation_checks']=literal_checks
# Deduplicate STATES, but verify the exact envelope rather than odd-part only.
checked=0
for seed in seeds:
    states={seed}
    for depth in range(7):
        for s in states:
            D,E,U,d=s
            for cap in (1,2,3,5,9):
                if D>cap: continue
                allowed=any(U==u*(1<<k) for u in range(1,2*cap+1) for k in range(U.bit_length()))
                allowed|=D==seed[0] and any(U==seed[2]*(1<<k) for k in range(U.bit_length()))
                need(allowed,'exact dyadic sequence envelope'); checked+=1
        states={move(s,k) for s in states for k in 'OZJ'}
counts['exact_envelope_state_cap_checks']=checked
hostile=0
for N in range(1,1000,2):
    cs=[c for c in (0,-2,-4,-6,-8) if (2*N-c)%10==0]
    need(len(cs)==1,'hostile residue partition')
    c=cs[0]; g=(2*N-c)//5; h=N-g
    if N>=13:
        need(h>0 and g>0 and h>g and 3*g>=2*h,'hostile positive chamber')
    hostile+=1
counts['hostile_candidate_checks']=hostile
print(json.dumps({'status':'PASS',**counts},sort_keys=True,indent=2))
