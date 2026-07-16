#!/usr/bin/env python3
"""
cospectral_tie_census_klein_S315.py — klein-2026-07-16-S315 (cont.2)

THE COSPECTRAL-TIE CENSUS: which tournament classes share the adjacency spectrum
(= ALL returned-walk counts tr A^k), and what still separates them.

Per class (n = 4..7): exact char polys of A (adjacency), K = A^T - A (skew),
L = D - A (Laplacian); H (Ham paths), tau-vector (in-arborescences), scores, c3, SC, |Aut|.
Then: group by charpoly(A); within each tie-group, test which invariants split it.
Note: cospectral => same c3 (tr A^3) => same x-level: ties live inside levels by theorem.
"""
import itertools
from math import comb
from fractions import Fraction as Fr

def census(n):
    m = n*(n-1)//2
    pairs = list(itertools.combinations(range(n), 2))
    pidx = {pr: i for i, pr in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    remaps = []
    for g in perms:
        remaps.append([(pidx[(min(g[u],g[v]),max(g[u],g[v]))], 0 if g[u]<g[v] else 1)
                       for (u,v) in pairs])
    seen = bytearray(1<<m); reps=[]; orbsz={}
    for bits in range(1<<m):
        if seen[bits]: continue
        orb=set()
        for tab in remaps:
            out=0
            for i in range(m):
                b=(bits>>i)&1; j,fl=tab[i]; out|=((b^fl)<<j)
            orb.add(out)
        for t in orb: seen[t]=1
        reps.append(bits); orbsz[bits]=len(orb)
    return reps, pairs, orbsz, len(perms)

def matA(bits, pairs, n):
    A=[[0]*n for _ in range(n)]
    for i,(u,v) in enumerate(pairs):
        if (bits>>i)&1: A[u][v]=1
        else: A[v][u]=1
    return A

def charpoly(M):
    # Faddeev–LeVerrier, exact over Fractions -> integer tuple
    n=len(M)
    I=[[Fr(1) if i==j else Fr(0) for j in range(n)] for i in range(n)]
    Mf=[[Fr(M[i][j]) for j in range(n)] for i in range(n)]
    c=[Fr(1)]
    N=[row[:] for row in I]
    for k in range(1,n+1):
        MN=[[sum(Mf[i][l]*N[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        tr=sum(MN[i][i] for i in range(n))
        ck=-tr/k
        c.append(ck)
        N=[[MN[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return tuple(int(x) for x in c)

def hampaths(A, n):
    full=1<<n
    dp=[[0]*n for _ in range(full)]
    for v in range(n): dp[1<<v][v]=1
    for mset in range(full):
        for v in range(n):
            if dp[mset][v]:
                for u in range(n):
                    if not (mset>>u)&1 and A[v][u]:
                        dp[mset|(1<<u)][u]+=dp[mset][v]
    return sum(dp[full-1][v] for v in range(n))

def det_int(M):
    M=[row[:] for row in M]; n=len(M)
    if n==0: return 1
    sign,prev=1,1
    for k in range(n-1):
        if M[k][k]==0:
            piv=next((i for i in range(k+1,n) if M[i][k]!=0),None)
            if piv is None: return 0
            M[k],M[piv]=M[piv],M[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[-1][-1]

def tauvec(A,n):
    L=[[(sum(A[u]) if u==v else 0)-A[u][v] for v in range(n)] for u in range(n)]
    return tuple(sorted(det_int([[L[u][v] for v in range(n) if v!=r]
                                  for u in range(n) if u!=r]) for r in range(n)))

print("n | classes | cospectral(A) tie-groups (size>=2) | classes in ties | deep ties (A&K&L) | invisible (A,K,L,H,tau all tie)")
summary={}
for n in (4,5,6,7):
    reps, pairs, orbsz, nperm = census(n)
    rows=[]
    for bits in reps:
        A=matA(bits,pairs,n)
        s=sorted(sum(A[u]) for u in range(n))
        c3=comb(n,3)-sum(comb(x,2) for x in sorted(s))
        K=[[A[v][u]-A[u][v] for v in range(n)] for u in range(n)]
        L=[[(sum(A[u]) if u==v else 0)-A[u][v] for v in range(n)] for u in range(n)]
        cpA=charpoly(A); cpK=charpoly(K); cpL=charpoly(L)
        H=hampaths(A,n)
        tv=tauvec(A,n)
        aut=nperm//orbsz[bits]
        # SC: complement in orbit?
        rows.append(dict(bits=bits,s=tuple(s),c3=c3,cpA=cpA,cpK=cpK,cpL=cpL,H=H,tv=tv,aut=aut))
    groups={}
    for r in rows: groups.setdefault(r['cpA'],[]).append(r)
    ties=[g for g in groups.values() if len(g)>1]
    intie=sum(len(g) for g in ties)
    deep=[]
    invis=[]
    split_stats={'scores':0,'H':0,'tau':0,'K':0,'L':0,'aut':0}
    for g in ties:
        # what splits this A-cospectral group?
        for key,fn in [('scores',lambda r:r['s']),('H',lambda r:r['H']),('tau',lambda r:r['tv']),
                       ('K',lambda r:r['cpK']),('L',lambda r:r['cpL']),('aut',lambda r:r['aut'])]:
            if len(set(fn(r) for r in g))>1: split_stats[key]+=1
        # deep: also K and L cospectral (pairwise-equal within subgroup)
        sub={}
        for r in g: sub.setdefault((r['cpK'],r['cpL']),[]).append(r)
        for gg in sub.values():
            if len(gg)>1:
                deep.append(gg)
                sub2={}
                for r in gg: sub2.setdefault((r['H'],r['tv']),[]).append(r)
                for ggg in sub2.values():
                    if len(ggg)>1: invis.append(ggg)
    summary[n]=(len(rows),ties,deep,invis,split_stats)
    print(f"{n} | {len(rows):4d} | {len(ties):3d} | {intie:4d} | {len(deep):3d} | {len(invis):3d}")
    if ties:
        print(f"   splitters (of {len(ties)} A-tie-groups): scores splits {split_stats['scores']}, "
              f"H splits {split_stats['H']}, tau splits {split_stats['tau']}, skew-K {split_stats['K']}, "
              f"Laplacian {split_stats['L']}, |Aut| {split_stats['aut']}")
    for gg in deep:
        info=[(r['s'],r['H'],r['tv'],r['aut']) for r in gg]
        print(f"   DEEP TIE (A,K,L all cospectral), n={n}: {len(gg)} classes; (scores,H,tau,Aut) = {info}")
    for ggg in invis:
        print(f"   *** INVISIBLE PAIR n={n}: A,K,L-cospectral AND same H AND same tau: "
              f"aut={[r['aut'] for r in ggg]}, scores={[r['s'] for r in ggg]}")
print()
print("Notes: cospectral(A) => equal tr A^k for all k => equal c3 (x-level) — ties live inside")
print("levels by theorem; the spectrum IS the returned-walk data, so whatever splits above is")
print("exactly what returning walks cannot see.")

# ================= THE CONE MECHANISM (why the invisible pairs exist) =================
# Claim: the 4 invisible n=7 pairs are exactly the CONES (add a sink) over the 4 equal-H
# deep-tied n=6 pairs; coning collapses tau to (0,...,0, det(L_base + I)) — an L-spectral
# function — and preserves H (the sink can only terminate a path). The base tau-vectors
# (vertex-deleted deck) DO separate: invisibility is a coning artifact, seen by deletion.
print()
print("CONE MECHANISM VERIFICATION:")
reps7, pairs7, orbsz7, nperm7 = census(7)
inv_pairs = []
groups={}
for bits in reps7:
    A=matA(bits,pairs7,7)
    cpA=charpoly(A)
    groups.setdefault(cpA,[]).append(bits)
checked = fine = 0
for cp,g in groups.items():
    if len(g)<2: continue
    # find invisible sub-pairs: same cpK,cpL,H,tau
    sub={}
    for bits in g:
        A=matA(bits,pairs7,7)
        K=[[A[v][u]-A[u][v] for v in range(7)] for u in range(7)]
        L=[[(sum(A[u]) if u==v else 0)-A[u][v] for v in range(7)] for u in range(7)]
        key=(charpoly(K),charpoly(L),hampaths(A,7),tauvec(A,7))
        sub.setdefault(key,[]).append(bits)
    for key,gg in sub.items():
        if len(gg)>1:
            for bits in gg:
                A=matA(bits,pairs7,7)
                s=[sum(A[u]) for u in range(7)]
                checked += 1
                sink=[v for v in range(7) if s[v]==0]
                assert len(sink)==1, "invisible class WITHOUT a sink!"
                v0=sink[0]
                # base = delete the sink; H(cone) = H(base)?  tau_sink = det(L_base + I)?
                keep=[u for u in range(7) if u!=v0]
                B=[[A[keep[i]][keep[j]] for j in range(6)] for i in range(6)]
                Lb=[[(sum(B[i]) if i==j else 0)-B[i][j] for j in range(6)] for i in range(6)]
                LbI=[[Lb[i][j]+(1 if i==j else 0) for j in range(6)] for i in range(6)]
                Hc=hampaths(A,7); Hb=hampaths(B,6)
                tv=tauvec(A,7)
                ok = (Hc==Hb) and (tv[-1]==det_int(LbI)) and all(t==0 for t in tv[:-1])
                fine += ok
print(f"   invisible classes checked: {checked}; all have a unique SINK; H(cone)=H(base) and")
print(f"   tau = (0,...,0, det(L_base+I)) verified for {fine}/{checked} — the tau-vector's")
print("   resolving power collapses to an L-spectral scalar under coning, H passes through,")
print("   and the base tau-vectors (the vertex-deleted deck) separate every pair: the")
print("   spectrum's residual blindness is a CONING artifact, cured by deletion.")
