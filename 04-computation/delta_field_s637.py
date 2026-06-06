#!/usr/bin/env python3
"""
S637 — the POLARIZED DELTA FIELD: delta = dH, the discrete gradient of H on the tournament hypercube
(vertices=tournaments, edges=arc flips). H = #Hamiltonian paths = I(Omega,2) = OCF (Redei: H odd).
H is ANTIFERROMAGNETIC (THM-290): external field lambda_k=log(1+2^{s-1}), couplings<0 (frustrated),
Walsh-even & bandlimited (THM-163/259/260). So the delta field is the gradient of a frustrated Ising
energy; its POLARIZATION = sign field (toward the regular tournament = H-max); its DEFECTS = level
edges (delta=0) and reversed edges, which grow with n.
We compute, at n=5,6,7: the H-spectrum; the gradient field; polarization; level/reversed edges; the
critical tournaments (local extrema = metastable states); and check Walsh even/bandlimited.
"""
import itertools
from collections import Counter

def arcs(n): return [(i,j) for i in range(n) for j in range(n) if i<j]

def H_hampaths(adjbits, n, A):
    # adj[i][j]=True iff i->j ; H = #Hamiltonian paths via Held-Karp DP
    adj=[[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(A):
        if (adjbits>>k)&1: adj[i][j]=True
        else: adj[j][i]=True
    full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            c=dp[mask][v]
            if not c: continue
            for w in range(n):
                if not (mask>>w)&1 and adj[v][w]:
                    dp[mask|(1<<w)][w]+=c
    return sum(dp[full])

def field_stats(n, sample=None):
    A=arcs(n); m=len(A); M=1<<m
    import random
    rng=random.Random(0)
    if sample and M>sample:
        idxs=[rng.randrange(M) for _ in range(sample)]
    else:
        idxs=range(M)
    spec=Counter(); level_edges=0; reversed_edges=0; updown=0
    n_local_max=0; n_local_min=0; pol=Counter()  # polarization buckets
    seen={}
    for bits in idxs:
        H0=seen.get(bits) or H_hampaths(bits,n,A); seen[bits]=H0; spec[H0]+=1
        ds=[]
        for k in range(m):
            Hn=H_hampaths(bits^(1<<k),n,A)
            ds.append(Hn-H0)
        pos=sum(1 for d in ds if d>0); neg=sum(1 for d in ds if d<0); zer=sum(1 for d in ds if d==0)
        level_edges+=zer
        if pos==0 and neg+zer==m and neg>0: n_local_max+=1   # all neighbors <= (and some <)
        if neg==0 and pos+zer==m and pos>0: n_local_min+=1
        # polarization = (pos-neg)/m  bucket
        pol[round((pos-neg)/m,1)]+=1
    return spec, level_edges, n_local_max, n_local_min, pol, len(list(idxs))

if __name__=="__main__":
    for n,samp in [(5,None),(6,None),(7,4000)]:
        A=arcs(n); m=len(A)
        spec,level,lmax,lmin,pol,cnt=field_stats(n,samp)
        forb=[h for h in range(1,max(spec)+1,2) if h not in spec]
        tag="(full)" if samp is None or (1<<m)<=samp else f"(sample {cnt})"
        print(f"n={n} {tag}: m={m} arcs.  H-spectrum={sorted(spec)}  forbidden-odd={forb}")
        print(f"   level edges (delta=0) per config avg = {level/cnt:.2f}  (total {level})")
        print(f"   local MAXima of H (all-neighbors-<=): {lmax}   local MINima: {lmin}")
        print(f"   polarization (pos-neg)/m distribution: {dict(sorted(pol.items()))}")
        print()
    print("READING: delta=dH = gradient of the antiferromagnetic energy on the spin cube (THM-290).")
    print("  polarization = sign of the gradient (toward H-max=regular tournament); level edges (delta=0)")
    print("  = the FRUSTRATION DEFECTS where the field cancels; local extrema = metastable states.")
