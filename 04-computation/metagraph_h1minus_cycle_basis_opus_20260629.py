"""
Find a basis of the 7-dim R-odd homology H_1^-(metagraph, n=5) and test for FANO/heptagon structure.
Fano = cyclic heptagon, 7 lines = translates of difference set {1,2,4} mod 7 (each pair of lines meets in 1 pt).
Robbins(graph): cycles = strongly-connected closed walks. Duocylinder: two-scale (2*7) product.
"""
from itertools import permutations
import numpy as np
from fractions import Fraction as F
def canon(n, arcs):
    best=None
    for p in permutations(range(n)):
        key=tuple(sorted((p[i],p[j]) for (i,j) in arcs))
        if best is None or key<best: best=key
    return best
n=5
base=set((k,k-1) for k in range(n-1,0,-1))
tiles=sorted([(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in base and (j,i) not in base])
m=len(tiles)
def tour(bits):
    arcs=set(base)
    for idx,(i,j) in enumerate(tiles):
        arcs.add((i,j) if (bits>>idx)&1 else (j,i))
    return arcs
cl={}; tc=[]
for bits in range(1<<m):
    c=canon(n,tour(bits))
    if c not in cl: cl[c]=len(cl)
    tc.append(cl[c])
V=len(cl); inv={i:c for c,i in cl.items()}
# edges
edgeset=set()
for bits in range(1<<m):
    a=tc[bits]
    for idx in range(m):
        b=tc[bits^(1<<idx)]
        if a!=b: edgeset.add((min(a,b),max(a,b)))
edges=sorted(edgeset); E=len(edges); eidx={e:i for i,e in enumerate(edges)}
# R = complement
def comp(arcs): return set((j,i) for (i,j) in arcs)
R={i:cl[canon(n,comp(set(c)))] for c,i in cl.items()}
print(f"n=5 metagraph: V={V} E={E}")
# boundary matrix (over Q via numpy float then verify ranks)
B=np.zeros((V,E))
for k,(u,v) in enumerate(edges): B[u,k]=-1; B[v,k]=1
# cycle space = null space of B
from numpy.linalg import svd
U,s,Vt=svd(B)
rank=sum(1 for x in s if x>1e-9); b1=E-rank
# R on edges: signed permutation
Redge=np.zeros((E,E))
for k,(u,v) in enumerate(edges):
    ru,rv=R[u],R[v]
    e=(min(ru,rv),max(ru,rv)); kk=eidx[e]
    Redge[kk,k]= 1 if (ru<rv)==(u<v) else -1
# cycle space basis (null space of B)
ns=Vt[rank:].T   # E x b1, columns span cycle space
# R restricted to cycle space: project Redge into ns coords
# solve ns @ X = Redge @ ns  => X = pinv(ns) @ Redge @ ns
Rcyc=np.linalg.pinv(ns) @ Redge @ ns
# -1 eigenspace
w,vec=np.linalg.eig(Rcyc)
odd=[i for i in range(len(w)) if abs(w[i]+1)<1e-6]
print(f"b1={b1}; R-eigenvalues on cycle space: +1 x{sum(1 for x in w if abs(x-1)<1e-6)}, -1 x{len(odd)}  => b1^- = {len(odd)}")
# R-odd cycles in edge space
oddcycles=[ns @ vec[:,i].real for i in odd]
# find SHORT integer R-odd cycles by rounding/support
print(f"\nNS classes (R-pairs, non-self-converse):", [i for i in range(V) if R[i]!=i])
print(f"SC classes (R-fixed):", [i for i in range(V) if R[i]==i])
# support sizes of a natural integer basis: use the edges flipped by R (R-odd edges)
Rodd_edges=[k for k,(u,v) in enumerate(edges) if (min(R[u],R[v]),max(R[u],R[v]))==(u,v) and R[u]!=u]  # R-fixed edges reversed
print(f"R-reversed fixed edges (candidate 'points'): {[(edges[k]) for k in Rodd_edges]}")
# incidence: build the 7-dim space, look at which edges appear
supp=np.zeros(E)
for c in oddcycles: supp += np.abs(c)
active=[k for k in range(E) if supp[k]>1e-6]
print(f"edges active in R-odd homology: {len(active)} of {E}")
print(f"active edges: {[edges[k] for k in active]}")
