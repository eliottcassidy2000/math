#!/usr/bin/env python3
"""Structural controls for integer gauge packets, not a coverage census."""
from itertools import product, permutations
from collections import Counter
from math import comb, factorial
import json
from pathlib import Path

def mask_tournament(n,m):
    rows=[0]*n; t=0
    for i in range(n):
        for j in range(i+1,n):
            if m>>t&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
            t+=1
    return rows

def ordered(a,b):
    n=len(a); m=len(b)
    return [r|(((1<<m)-1)<<n) for r in a]+[r<<n for r in b]

def H(a):
    n=len(a)
    if n==0:return 1
    dp=[{} for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,c in dp[mask].items():
            todo=a[v]&~mask
            while todo:
                bit=todo&-todo; todo-=bit; w=bit.bit_length()-1
                d=dp[mask|bit]; d[w]=d.get(w,0)+c
    return sum(dp[-1].values())

def tt(n):return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
tri=[2,4,1]
G={}
for i,j in product(range(2),repeat=2):
    a=ordered(mask_tournament(5,8+2*i),mask_tournament(5,8+2*j))
    a[0]^=1<<5; a[5]^=1
    G[i,j]=a
W={(0,0):215,(0,1):179,(1,0):245,(1,1):201}
for key,a in G.items():assert H(a)==W[key]

def filler(K,k,i,j):
    r=k if (i,j)==(0,1) else K-k if (i,j)==(1,0) else 0
    a=[]
    for _ in range(r):a=ordered(a,tri)
    return ordered(a,tt(3*(K-r)))

def signature(K,k,x):
    # Exact ordered SCC signature. G's unique H labels identify their classes.
    s=[]; N=len(x)
    for pos,i in enumerate(x):
        j=x[(pos+1)%N]
        r=k if (i,j)==(0,1) else K-k if (i,j)==(1,0) else 0
        s.append(('G',i,j)); s.extend([('C3',)]*r); s.extend([('pt',)]*(3*(K-r)))
    return tuple(s)

def hprod(K,k,x):
    h=1; N=len(x)
    for p,i in enumerate(x):
        j=x[(p+1)%N]
        r=k if (i,j)==(0,1) else K-k if (i,j)==(1,0) else 0
        h*=W[i,j]*3**r
    return h

# Independent subset DP on physical 13-vertex cells.
cell_tests=0
for k in range(2):
    for i,j in product(range(2),repeat=2):
        a=ordered(G[i,j],filler(1,k,i,j))
        r=k if (i,j)==(0,1) else 1-k if (i,j)==(1,0) else 0
        assert H(a)==W[i,j]*3**r
        cell_tests+=1

# All words K1..4, N1..9. Not giant graphs: SCC theorem supplies exact signatures.
cases=[]
for K in range(1,5):
  for N in range(1,10):
    seen={}
    for x in product(range(2),repeat=N):
        values={hprod(K,k,x) for k in range(K+1)}
        assert len(values)==1
        for k in range(K+1):seen[signature(K,k,x)]=hprod(K,k,x)
    got=Counter(seen.values())
    want=Counter({215**N:1,201**N:1})
    for n0 in range(1,N):
      for r in range(1,min(n0,N-n0)+1):
        h=215**(n0-r)*201**(N-n0-r)*(179*245*3**K)**r
        mult=(K+1)*N*comb(n0-1,r-1)*comb(N-n0-1,r-1)//r
        want[h]+=mult
    assert got==want
    assert len(got)==N*N//4+2
    assert len(seen)==2+(K+1)*(2**N-2)
    cases.append({'K':K,'N':N,'classes':len(seen),'values':len(got)})

# Boundary control: positive valuation drift fails for an unclosed walk.
def open_weight(K,k,x):
    h=1
    for i,j in zip(x,x[1:]):
        r=k if (i,j)==(0,1) else K-k if (i,j)==(1,0) else 0
        h*=W[i,j]*3**r
    return h
for K in range(1,5):
 for k in range(K+1):
  for n in range(1,7):
   for x in product(range(2),repeat=n):
    a=open_weight(K,k,x);b=open_weight(K,0,x)
    e=k*(x[-1]-x[0])
    assert (a==b*3**e) if e>=0 else (a*3**(-e)==b)

out={'status':'PASS','physical_subset_DP_cells':cell_tests,'whole_word_checks':'exact SCC signatures and products, not independent whole-graph DP','packet_censuses':cases,'boundary_formula':'PASS','scope':'K1..4,N1..9 controls; theorems all K>=1,N>=1; no Chan-Pak coverage claim'}
path=Path(__file__).with_name('RESULTS.json');path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

# Physical endpoint probe on smaller prefix controls, to verify local exposure.
def seam(a,v):
    n=len(a);b=ordered(a,[0]);b[v]^=1<<n;b[n]^=1<<v
    return b
probe=[]
for prefix in ([0],tri,G[0,1]):
    a=ordered(prefix,tri)
    b=ordered(prefix,tt(3))
    hp=H(prefix)
    assert H(seam(a,len(a)-1))==5*hp
    assert H(seam(b,len(b)-1))==hp
    probe.append({'prefix_order':len(prefix),'prefix_H':hp,'triangle_output':5*hp,'transitive_output':hp})
out['physical_endpoint_seam_controls']=probe
out['26_vertex_packet_seam_outputs']={'before':[131565,131565],'after':[219275,131565],'method':'proved SCC product and independently tested local seam, not 26-vertex subset DP'}
path.write_text(json.dumps(out,indent=2)+'\n')
print('endpoint exposure: PASS')
