#!/usr/bin/env python3
"""Signed even-subgraph counts Phi_X(T) (X even, T tournament types) for n <= 5: orthogonality relations
 M^T D_T M = 2^m D_X  and  M D_X^{-1} M^T = 2^m D_T^{-1}  (D_T = diag(n!/|Aut T|), D_X = diag(n!/|Aut X|), m = C(n,2)),
 det M = +-2^k with k = mN/2 + (1/2) log2(prod |Aut T| / prod |Aut X|)."""
import itertools, subprocess
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher, DiGraphMatcher
from sympy import Matrix, Rational, log, factorial
from math import prod
def geng(n):
    return [nx.from_graph6_bytes(s.encode()) for s in subprocess.run(['geng','-q',str(n)],capture_output=True,text=True).stdout.split()]
def gentourng(n):
    Ts=[]
    for line in subprocess.run(['gentourng','-q',str(n)],capture_output=True,text=True).stdout.split():
        D=nx.DiGraph(); D.add_nodes_from(range(n)); k=0
        for i in range(n):
            for j in range(i+1,n):
                D.add_edge(i,j) if line[k]=='1' else D.add_edge(j,i); k+=1
        Ts.append(D)
    return Ts
def sign(E,g):
    s=0
    for u,v in E:
        if g[u]>g[v]: s^=1
    return -1 if s else 1
for n in range(3,6):
    perms=list(itertools.permutations(range(n))); m=n*(n-1)//2
    evens=[]
    for G in geng(n):
        E=[(min(u,v),max(u,v)) for u,v in G.edges()]
        auts=[tuple(mp[i] for i in range(n)) for mp in GraphMatcher(G,G).isomorphisms_iter()]
        if all(sign(E,g)==1 for g in auts): evens.append((E,len(auts)))
    tours=[]
    for T in gentourng(n):
        auts=sum(1 for _ in DiGraphMatcher(T,T).isomorphisms_iter())
        tours.append(({(u,v):(1 if T.has_edge(u,v) else -1) for u in range(n) for v in range(u+1,n)},auts))
    rows=[]
    for tv,_ in tours:
        row=[]
        for E,_ in evens:
            copies={}
            for g in perms:
                img=frozenset((min(g[u],g[v]),max(g[u],g[v])) for u,v in E); copies[img]=sign(E,g)
            row.append(sum(eps*prod(tv[e] for e in img) for img,eps in copies.items()))
        rows.append(row)
    M=Matrix(rows); N=len(rows); nf=factorial(n)
    DT=Matrix.diag(*[Rational(nf,a) for _,a in tours]); DX=Matrix.diag(*[Rational(nf,a) for _,a in evens])
    o1=(M.T*DT*M==2**m*DX); o2=(M*DX.inv()*M.T==2**m*DT.inv())
    d=M.det(); pa=prod(a for _,a in tours); px=prod(a for _,a in evens)
    k=Rational(m*N,2)+Rational(1,2)*log(Rational(pa,px),2)
    print(f"n={n}: N={N}, orthogonality 1 {o1}, orthogonality 2 {o2}, det M = {d} = +-2^{k} -> {abs(d)==2**int(k)}; prod|Aut T|={pa}, prod|Aut X|={px}, ratio = 2^{log(Rational(pa,px),2)}")
