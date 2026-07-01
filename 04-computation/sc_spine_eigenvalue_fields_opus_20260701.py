"""Thread: are the blue (SC-spine) eigenvalues QUADRATIC, and which fields Q(sqrt d)? Factor the char poly exactly.
Plus the flip-line CHAIN COMPLEX: Betti numbers of blue (T-join), black (even graph=cycle), full metagraph."""
import numpy as np, sympy as sp
from itertools import permutations
def build_graphs(n):
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; m=len(TILES); ti={t:i for i,t in enumerate(TILES)}
    TRANS=[ti[(n-y+1,n-x+1)] for (x,y) in TILES]; perms=list(permutations(range(n)))
    def adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=n-xL; yi=n-yL; A[xi][yi]=1 if bits[i]==0 else 0; A[yi][xi]=1-A[xi][yi]
        return A
    def canon(A):
        b=None
        for p in perms:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if b is None or s<b: b=s
        return b
    allt=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]; allt.append((mask,bits))
    cls={}; node={}
    def cid(bits):
        c=canon(adj(bits))
        if c not in cls: cls[c]=len(cls)
        return cls[c]
    for mask,bits in allt: node[mask]=cid(bits)
    N=len(cls); Blue=np.zeros((N,N),int); Black=np.zeros((N,N),int); seen=set()
    for mask,bits in allt:
        gs=all(bits[i]==bits[TRANS[i]] for i in range(m) if TRANS[i]!=i)
        fm=mask^((1<<m)-1); key=(min(mask,fm),max(mask,fm))
        if key in seen: continue
        seen.add(key); a=node[mask]; b=node[fm]; G=Blue if gs else Black
        if a==b: G[a,a]+=1
        else: G[a,b]+=1; G[b,a]+=1
    return Blue,Black,N,cls
def betti(G):  # simple-graph Betti over the support (nonzero off-diag); H0=comp, H1=E-V+comp
    N=G.shape[0]; adj=[[j for j in range(N) if j!=i and G[i,j]>0] for i in range(N)]
    nodes=[i for i in range(N) if adj[i] or G[i,i]>0]
    seen=set(); comp=0
    for s in nodes:
        if s in seen: continue
        comp+=1; st=[s]
        while st:
            u=st.pop()
            if u in seen: continue
            seen.add(u)
            for w in adj[u]:
                if w not in seen: st.append(w)
    E=sum(1 for i in range(N) for j in range(i+1,N) if G[i,j]>0)  # simple edges
    V=len(nodes)
    return comp, E-V+comp if V else 0
for n in [4,5,6]:
    Blue,Black,N,cls=build_graphs(n)
    M=sp.Matrix(Blue.tolist())
    cp=M.charpoly()
    facs=sp.factor_list(cp.as_expr())
    # collect irreducible factor degrees + quadratic fields
    degs=[]; fields=set()
    for f,mult in facs[1]:
        d=sp.degree(f, sp.Symbol('lambda'))
        degs.append((int(d),int(mult)))
        if d==2:
            # discriminant of the quadratic
            p=sp.Poly(f, sp.Symbol('lambda')); a,b,c=p.all_coeffs()
            disc=sp.nsimplify(b*b-4*a*c); fields.add(sp.sqrtdenest(sp.sqrt(disc)))
    b0B,b1B=betti(Blue); b0K,b1K=betti(Black)
    print(f"n={n}: {N} SC blue-nodes; char-poly irreducible-factor (degree x mult): {sorted(set(degs))}")
    print(f"   quadratic fields (sqrt disc): {sorted(str(x) for x in fields)}")
    print(f"   CHAIN COMPLEX Betti: blue(H0,H1)=({b0B},{b1B}); black even-graph(H0,H1)=({b0K},{b1K})")
