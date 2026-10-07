#!/usr/bin/env python3
"""(c) book Ramsey vs the cycle-clique extremal family.
1. Turan-type colourings (red = disjoint cliques V_1..V_k, blue = complete multipartite) with no red B_{n-1}
   (red codegree <= n-2) and no blue B_n (blue codegree <= n-1): maximal order, exhaustive over partitions.
2. Stored n=12 two-level circulant witness (R(B_11,B_12) >= 47): verify; alpha/omega of red and blue;
   what the cycle-clique theorem R(C_m,K_s)=(m-1)(s-1)+1 (m>=s>=3) implies for it."""
import itertools, sys
import networkx as nx
import numpy as np

def partitions(N, maxpart=None):
    if maxpart is None: maxpart=N
    if N==0: yield []; return
    for p in range(min(N,maxpart),0,-1):
        for rest in partitions(N-p,p): yield [p]+rest

def turan_ok(parts,n):
    N=sum(parts)
    if max(parts)>n: return False                 # red K_s: edge codegree s-2 <= n-2
    k=len(parts)
    if k>=2:
        for i in range(k):
            for j in range(i+1,k):
                if N-parts[i]-parts[j]>n-1: return False   # blue edge between parts i,j
    return True
print("1. max N of a Turan-type (B_{n-1},B_n)-good colouring vs target 4n-2:")
for n in range(2,11):
    best=0; arg=None
    for N in range(1,4*n):
        found=None
        for P in partitions(N, n):
            if turan_ok(P,n): found=P; break
        if found: best=N; arg=found
    print(f"   n={n:2d}: max N={best:3d} (parts {arg}); max(2n,3n-3)={max(2*n,3*n-3):3d}; target 4n-2={4*n-2}")
n=100
print(f"   n=100 (formula): Turan-type max = {max(2*n,3*n-3)} vertices vs 398 needed")

# brute verification of the Turan colouring for small n via explicit codegrees
def build_turan(parts):
    N=sum(parts); lab=[]
    for i,p in enumerate(parts): lab+= [i]*p
    A=np.array([[1 if (a!=b and lab[a]==lab[b]) else 0 for b in range(N)] for a in range(N)])
    return A
def book_ok(A,n):
    N=len(A); M=A.astype(int); C=M@M; Bm=1-M; np.fill_diagonal(Bm,0); Cb=Bm@Bm
    for i in range(N):
        for j in range(i+1,N):
            if M[i,j] and C[i,j]>n-2: return False
            if not M[i,j] and Cb[i,j]>n-1: return False
    return True
for n in [4,5,6,7]:
    P=[n-1]*3
    print(f"   explicit check n={n}: three red K_{n-1} (N={3*n-3}) valid: {book_ok(build_turan(P),n)}; adding one vertex to a part valid: {book_ok(build_turan([n,n-1,n-1]),n)}")

# 2. stored witness
fn=sys.argv[1]
with open(fn) as f:
    head=f.readline().split(); nn=int(head[1]); m=int(head[3])
    U0=[int(c) for c in f.readline().split()[1]]
    U1=[int(c) for c in f.readline().split()[1]]
    D=[int(c) for c in f.readline().split()[1]]
N=4*nn-2
A=np.zeros((N,N),dtype=int)
for i in range(m):
    for j in range(m):
        dd=(j-i)%m
        if i!=j:
            if U0[dd]: A[i,j]=1
            if U1[dd]: A[m+i,m+j]=1
        if D[dd]: A[i,m+j]=1; A[m+j,i]=1
assert (A==A.T).all()
print(f"2. stored witness n={nn}, N={N}: |U0|={sum(U0)} |U1|={sum(U1)} |D|={sum(D)}; valid (no red B_{nn-1}, no blue B_{nn}): {book_ok(A,nn)}")
Gr=nx.from_numpy_array(A); Gb=nx.complement(Gr)
def omega(G): return max(len(c) for c in nx.find_cliques(G))
wr,wb=omega(Gr),omega(Gb)
ar,ab=wb,wr
degs=sorted(set(dict(Gr.degree()).values()))
print(f"   red degrees {degs}; omega(red)={wr}, alpha(red)={ar}; omega(blue)={wb}, alpha(blue)={ab}")
for name,alpha in [("red",ar),("blue",ab)]:
    s=alpha+1
    ms=[mm for mm in range(max(s,3),N+1) if (mm-1)*(s-1)+1<=N]
    print(f"   cycle-clique => {name} graph (alpha={alpha}) contains C_m for m in {ms} (needs m>=s={s}); all these pairs (m,s) have s<=7 -> pre-2026 literature: {all(s<=7 for _ in ms)}")
# direct cycle check for those lengths (DFS for a cycle of exact length)
def has_cycle_len(G,L):
    nodes=list(G.nodes())
    adj={v:set(G[v]) for v in nodes}
    for start in nodes:
        stack=[(start,[start])]
        while stack:
            v,path=stack.pop()
            if len(path)==L:
                if start in adj[v]: return True
                continue
            for w in adj[v]:
                if w>start and w not in path: stack.append((w,path+[w]))
            if len(stack)>200000: break
    return False
print("   direct red cycle lengths 3..12 present:", [L for L in range(3,13) if has_cycle_len(Gr,L)])
