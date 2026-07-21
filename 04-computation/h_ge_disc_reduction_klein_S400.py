#!/usr/bin/env python3
"""
klein-2026-07-21-S400.  Toward proving HYP-8636: H(T) >= disc(T), equality iff transitive.
  H(T)   = # Hamiltonian paths (Redei, #P-hard, odd).
  disc(T)= |det(I+K)|/2^{n-1},  K=A-A^T skew,  = prod_j (1+mu_j^2)/2^{n-1}  (+-i*mu_j = eig K),
           with the FIXED energy constraint  sum_j mu_j^2 = C(n,2)  (=#arcs, same for all T of order n).

STEP 1 (the reduction, mirror of c3<=H THM-1860): is disc MULTIPLICATIVE over strong components?
  H is: H(T)=prod H(SCC_i) (Ham path traverses the linearly-ordered strong comps in order).
  If disc(T)=prod disc(SCC_i) too, then H>=disc REDUCES TO STRONGLY CONNECTED tournaments, and the
  equality-iff-transitive locus (transitive = all singleton SCCs, H=disc=1 each) drops out cleanly.

Also: (a) reconfirm H>=disc + min ratio cases; (b) tabulate H,disc on strong tournaments; (c) the
principal-minor / Pfaffian structure of det(I+K).
"""
import itertools, math
try:
    import numpy as np
except Exception:
    np = None

def all_tournaments(n):
    """iterate labelled tournaments as adjacency (A[i][j]=1 iff i->j), upper-triangle choices."""
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        A = [[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def ham_paths(A, n):
    full=(1<<n)-1
    out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        row=dp[mask]
        for v in range(n):
            c=row[v]
            if c:
                ov=out[v]
                for w in range(n):
                    if not (mask>>w&1) and (ov>>w&1):
                        dp[mask|1<<w][w]+=c
    return sum(dp[full][v] for v in range(n))

def disc(A, n):
    K=[[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            K[i][j]=A[i][j]-A[j][i]
    M=[[(1.0 if i==j else 0.0)+K[i][j] for j in range(n)] for i in range(n)]
    d=abs(np.linalg.det(np.array(M)))
    return d/(2**(n-1))

def strong_components(A, n):
    """Tarjan-ish via reachability; return list of components (as frozensets), condensation order."""
    # reachability closure
    reach=[set([i]) for i in range(n)]
    for _ in range(n):
        for i in range(n):
            for j in list(reach[i]):
                for k in range(n):
                    if A[j][k]: reach[i].add(k)
    comp={}
    seen=set()
    order=[]
    for i in range(n):
        if i in seen: continue
        c=frozenset(k for k in range(n) if k in reach[i] and i in reach[k])
        for v in c:
            comp[v]=c; seen.add(v)
        order.append(c)
    # unique comps
    uniq=[]
    us=set()
    for c in order:
        if c not in us: us.add(c); uniq.append(c)
    return uniq

def subA(A, S):
    S=sorted(S); m=len(S)
    return [[A[S[a]][S[b]] for b in range(m)] for a in range(m)], m

# ---------------- STEP 1: multiplicativity of disc over SCC ----------------
print("="*78)
print("STEP 1: is disc MULTIPLICATIVE over strong components?  (H is)")
print("="*78)
bad_disc=0; bad_H=0; tested=0
mult_examples=[]
for n in range(3,7):
    for A in all_tournaments(n):
        comps=strong_components(A,n)
        if len(comps)==1:  # already strong, trivial
            continue
        tested+=1
        Hf=ham_paths(A,n); Df=disc(A,n)
        Hp=1.0; Dp=1.0
        for c in comps:
            sa,m=subA(A,c)
            Hp*=ham_paths(sa,m); Dp*=disc(sa,m)
        if abs(Df-Dp)>1e-6*max(1,Df): bad_disc+=1
        if Hf!=round(Hp): bad_H+=1
    print(f"  n={n}: tested {tested} reducible labelled tournaments so far; "
          f"disc-mult violations={bad_disc}, H-mult violations={bad_H}")
print(f"RESULT: disc multiplicative over SCC? {'YES' if bad_disc==0 else 'NO ('+str(bad_disc)+' viol)'}"
      f"; H multiplicative? {'YES' if bad_H==0 else 'NO'}")
print("  => if YES, H>=disc REDUCES TO STRONGLY CONNECTED tournaments.\n")

# ---------------- STEP 2: H vs disc on STRONG tournaments, min ratio ----------------
print("="*78)
print("STEP 2: H vs disc on STRONGLY CONNECTED tournaments (the reduced problem)")
print("="*78)
seen_iso={}
def canon(A,n):
    best=None
    for p in itertools.permutations(range(n)):
        key=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or key<best: best=key
    return best

for n in range(3,7):
    recs=[]
    seen=set()
    for A in all_tournaments(n):
        comps=strong_components(A,n)
        if len(comps)!=1: continue
        k=canon(A,n)
        if k in seen: continue
        seen.add(k)
        Hf=ham_paths(A,n); Df=disc(A,n)
        recs.append((Hf,Df,k))
    recs.sort(key=lambda r:r[0]/max(r[1],1e-9))
    minr=recs[0]
    print(f"  n={n}: {len(recs)} strong iso classes. min H/disc = {minr[0]/minr[1]:.4f} "
          f"(H={minr[0]}, disc={minr[1]:.4f}); max disc among strong = {max(r[1] for r in recs):.4f}, "
          f"max H = {max(r[0] for r in recs)}")
    # is min ratio always > 1 strictly for strong?
    allgt=all(r[0] >= r[1]-1e-6 for r in recs)
    strict=all(r[0] > r[1]+1e-6 for r in recs)
    print(f"        H>=disc all strong: {allgt}; STRICT (H>disc) all strong: {strict}")
print("\nDONE.")
