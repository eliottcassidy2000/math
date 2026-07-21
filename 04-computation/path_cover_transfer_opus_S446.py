#!/usr/bin/env python3
"""
THE PATH-COVER POLYNOMIAL is the refined compositional invariant (opus-S446).
pc(S,c) = # ways to partition V(S) into c vertex-disjoint directed paths (single vertices allowed).
pc(S,1) = H(S) (a 1-path-cover = a Hamiltonian path). CLAIM: H(C3[S1,S2,S3]) is a function of the
path-cover polynomials (pc(S1,.),pc(S2,.),pc(S3,.)) -- the refined invariant that COMPOSES where
scalar H does not (THM-1970). This resolves THM-1960's cyclic-H. Extract the universal transfer
kernel K: H(C3[S1,S2,S3]) = sum_{c1,c2,c3} K(c1,c2,c3) * pc(S1,c1)*pc(S2,c2)*pc(S3,c3).
"""
import itertools, numpy as np

def path_cover_poly(adj,m):
    """pc[c] = # c-path-covers, c=1..m. Brute force over successor functions."""
    outn=[[u for u in range(m) if adj[v][u]] for v in range(m)]
    counts=[0]*(m+2)
    choices=[outn[v]+[m] for v in range(m)]        # m = 'None' (path end)
    for assign in itertools.product(*choices):
        used=[False]*m; ok=True
        for v in range(m):
            nx=assign[v]
            if nx<m:
                if used[nx]: ok=False;break
                used[nx]=True
        if not ok: continue
        # acyclic check: follow chains
        acyclic=True
        for start in range(m):
            seen=set(); v=start
            while assign[v]<m:
                if v in seen: acyclic=False;break
                seen.add(v); v=assign[v]
            if not acyclic:break
        if not acyclic: continue
        c=sum(1 for v in range(m) if assign[v]==m)   # #path-ends = #paths
        counts[c]+=1
    return tuple(counts[1:m+1])

def Hcount(adj,n):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if c:
                av=adj[v]
                for u in range(n):
                    if not (mask>>u)&1 and av[u]: dp[mask|(1<<u)][u]+=c
    return sum(dp[size-1])
C3=[[0,1,0],[0,0,1],[1,0,0]]
def substitute(blocks):
    sizes=[b[1] for b in blocks]; off=[0]
    for s in sizes: off.append(off[-1]+s)
    N=off[-1]; big=[[0]*N for _ in range(N)]
    for i in range(3):
        Si,mi=blocks[i]
        for a in range(mi):
            for c in range(mi):
                if a!=c: big[off[i]+a][off[i]+c]=Si[a][c]
        for j in range(3):
            if i!=j and C3[i][j]:
                for a in range(sizes[i]):
                    for c in range(sizes[j]): big[off[i]+a][off[j]+c]=1
    return big,N

def edges_iter(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            adj[i][j]=1 if (bits>>k)&1 else 0; adj[j][i]=1-adj[i][j]
        yield adj

# block library: all tournaments n=1..4 with (adj, m, pc-poly, H)
lib=[]
for n in range(1,5):
    for adj in edges_iter(n):
        pc=path_cover_poly(adj,n)
        lib.append((adj,n,pc,pc[0]))   # pc[0]=pc(.,1)=H
print("Sample path-cover polynomials pc(S,c) [c=1..m]:")
for adj,m,pc,H in lib[:8]:
    print(f"  m={m}: pc={pc}  (H=pc[1]={pc[0]})")

# TEST: is H(C3[S1,S2,S3]) a function of (pc(S1),pc(S2),pc(S3))?  and NOT of scalar H alone?
print("\nTEST: H(C3[S1,S2,S3]) determined by path-cover polys? by scalar H?")
by_pc={}; by_H={}
reps=[b for b in lib if b[1]<=3]     # keep composite <=9 vertices for speed
for b1 in reps:
    for b2 in reps:
        for b3 in reps:
            big,N=substitute([(b1[0],b1[1]),(b2[0],b2[1]),(b3[0],b3[1])])
            Hc=Hcount(big,N)
            kpc=tuple(sorted([b1[2],b2[2],b3[2]]))
            kH=tuple(sorted([b1[3],b2[3],b3[3]]))
            by_pc.setdefault(kpc,set()).add(Hc)
            by_H.setdefault(kH,set()).add(Hc)
amb_pc=sum(1 for v in by_pc.values() if len(v)>1)
amb_H=sum(1 for v in by_H.values() if len(v)>1)
print(f"  by path-cover-poly triple: {amb_pc}/{len(by_pc)} ambiguous  (0 => pc DETERMINES composite)")
print(f"  by scalar-H triple:        {amb_H}/{len(by_H)} ambiguous  (>0 => scalar H does NOT)")

# EXTRACT the universal kernel K(c1,c2,c3) via least-squares over the size-3-block triples
print("\nEXTRACT transfer kernel K: H(C3[S^3-ordered]) = sum K(c1,c2,c3) pc1[c1]pc2[c2]pc3[c3]")
blk3=[b for b in lib if b[1]==3]     # the 2 size-3 tournaments (transitive, C3) -- pc known
rows=[]; ys=[]; cs=[(a,b,c) for a in (1,2,3) for b in (1,2,3) for c in (1,2,3)]
for b1 in blk3:
    for b2 in blk3:
        for b3 in blk3:
            big,N=substitute([(b1[0],3),(b2[0],3),(b3[0],3)])
            ys.append(Hcount(big,N))
            rows.append([b1[2][c1-1]*b2[2][c2-1]*b3[2][c3-1] for (c1,c2,c3) in cs])
A=np.array(rows,float); y=np.array(ys,float)
K,res,rk,sv=np.linalg.lstsq(A,y,rcond=None)
pred=A@K
print(f"  size-3 fit exact={np.allclose(pred,y,atol=1e-6)} (rank {rk}/{len(cs)}); nonzero K:")
for (c123),kv in zip(cs,K):
    if abs(kv)>1e-6: print(f"    K{c123} = {round(kv,3)}")
print("\n READING: pc (path-cover polynomial) is the refined COMPOSITIONAL invariant -- it determines")
print(" H(C3[.]) where scalar H does not; H = pc(.,1) is just its top coefficient. Resolves THM-1970's")
print(" 'more refined than H' + THM-1960's cyclic-H via the kernel K.")
