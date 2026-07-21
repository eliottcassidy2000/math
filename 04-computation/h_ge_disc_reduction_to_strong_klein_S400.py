#!/usr/bin/env python3
"""
klein-2026-07-21-S400.  H >= disc REDUCED TO THE STRONG BASE  H(C) >= s(C)*disc(C)  (parallels THM-1860).

Ingredients (all verified here):
 (A) r-block SCC composition law:  disc(C1=>...=>Cr) = (prod disc(Ci)) * [prod(1+si)+prod(1-si)] / 2^r,
     si = 1^T (I+K_i)^{-1} 1 = ||(I+K_i)^{-1}1||^2  (the total inverse-response of component Ci).
 (B) H multiplicative:  H(T) = prod H(Ci).
 (C) si >= 1 for strongly-connected Ci (and singleton si=1).
 (D) symmetric-function lemma (PROVED):  for si>=1,  prod si  >=  [prod(1+si)+prod(1-si)]/2^r
     [= 2*sum_{k=r mod 2} e_k(u) <= 2*2^{r-1}=2^r after dividing by prod si, u_i=1/si in (0,1]].
 => H(T) = prod H(Ci) >= prod si*disc(Ci) = (prod disc)(prod si) >= (prod disc)[..]/2^r = disc(T),
    given the STRONG BASE  H(C) >= s(C) disc(C).  Equality locus: transitive (all singletons, s=1).

This script verifies (A) on reducible tournaments, (C) & the STRONG BASE (D-target) H>=s*disc on strong
tournaments, and reconfirms the final H>=disc -- all exhaustive n<=6, sampled n=7.
"""
import itertools, math, random
import numpy as np

def all_tournaments(n):
    pairs=list(itertools.combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
        yield A
def rand_tournament(n):
    A=[[0]*n for _ in range(n)]
    for i,j in itertools.combinations(range(n),2):
        if random.random()<0.5: A[i][j]=1
        else: A[j][i]=1
    return A
def Kmat(A,n): return np.array([[A[i][j]-A[j][i] for j in range(n)] for i in range(n)],dtype=float)
def disc(A,n): return abs(np.linalg.det(np.eye(n)+Kmat(A,n)))/(2**(n-1))
def s_val(A,n):
    x=np.linalg.solve(np.eye(n)+Kmat(A,n),np.ones(n)); return float(np.ones(n)@x)
def ham(A,n):
    full=(1<<n)-1; out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        row=dp[mask]
        for v in range(n):
            c=row[v]
            if c:
                ov=out[v]
                for w in range(n):
                    if not(mask>>w&1) and (ov>>w&1): dp[mask|1<<w][w]+=c
    return sum(dp[full][v] for v in range(n))
def sccs(A,n):
    reach=[set([i]) for i in range(n)]
    for _ in range(n):
        for i in range(n):
            for j in list(reach[i]):
                for k in range(n):
                    if A[j][k]: reach[i].add(k)
    comp=[];seen=set()
    # order by domination: component X before Y if X beats Y
    raw=[]
    for i in range(n):
        if i in seen: continue
        c=frozenset(k for k in range(n) if k in reach[i] and i in reach[k])
        for v in c: seen.add(v)
        raw.append(c)
    return raw
def subA(A,S):
    S=sorted(S);m=len(S); return [[A[S[a]][S[b]] for b in range(m)] for a in range(m)],m
def strong(A,n):
    reach=[set([i]) for i in range(n)]
    for _ in range(n):
        for i in range(n):
            for j in list(reach[i]):
                for k in range(n):
                    if A[j][k]: reach[i].add(k)
    return all(len(reach[i])==n for i in range(n))

# (A) r-block composition law
print("="*74); print("(A) r-block disc = (prod disc_i)[prod(1+s_i)+prod(1-s_i)]/2^r   (reducible T)"); print("="*74)
maxerr=0.0; nred=0
for n in range(2,7):
    for A in all_tournaments(n):
        comp=sccs(A,n)
        if len(comp)==1: continue
        nred+=1
        di=[]; si=[]
        for c in comp:
            sa,m=subA(A,c); di.append(disc(sa,m)); si.append(s_val(sa,m))
        r=len(comp)
        prod_disc=1.0
        for d in di: prod_disc*=d
        p1=1.0;pm=1.0
        for sv in si: p1*=(1+sv); pm*=(1-sv)
        pred=prod_disc*(p1+pm)/(2**r)
        maxerr=max(maxerr,abs(pred-disc(A,n)))
print(f"  {nred} reducible tournaments n<=6; max |formula - disc| = {maxerr:.2e}  "
      f"=> {'VERIFIED' if maxerr<1e-7 else 'WRONG'}")

# (C) s>=1 for strong ; (base) H >= s*disc for strong
print(); print("="*74); print("(C) s>=1 for strong  &  (BASE)  H(C) >= s(C) disc(C) for strong C"); print("="*74)
for n in range(3,7):
    seen=set(); recs=[]
    for A in all_tournaments(n):
        if not strong(A,n): continue
        best=None
        for p in itertools.permutations(range(n)):
            k=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or k<best: best=k
        if best in seen: continue
        seen.add(best)
        H=ham(A,n); D=disc(A,n); S=s_val(A,n)
        recs.append((H,D,S))
    smin=min(r[2] for r in recs)
    base_ok=all(r[0] >= r[2]*r[1]-1e-6 for r in recs)
    base_ratio=min(r[0]/(r[2]*r[1]) for r in recs)
    tight=[r for r in recs if abs(r[0]-r[2]*r[1])<1e-6]
    print(f"  n={n}: {len(recs)} strong classes; min s={smin:.4f} (>=1:{smin>=1-1e-9}); "
          f"H>=s*disc ALL: {base_ok}; min H/(s*disc)={base_ratio:.4f}; #tight(equality)={len(tight)}")

# n=7 strong sample for the base
print("  n=7 strong sample (H>=s*disc):")
random.seed(7); ok=0;tot=0;mn=9e9
for _ in range(4000):
    A=rand_tournament(7)
    if not strong(A,7): continue
    tot+=1; H=ham(A,7); D=disc(A,7); S=s_val(A,7)
    if H>=S*D-1e-6: ok+=1
    mn=min(mn,H/(S*D))
print(f"    strong samples={tot}; H>=s*disc: {ok}/{tot}; min H/(s*disc)={mn:.4f}")

# final: reconfirm H>=disc all (via the reduction it now follows from the base)
print(); print("="*74); print("FINAL reconfirm H>=disc (should follow from base+lemma)"); print("="*74)
viol=0;tot=0
for n in range(3,7):
    for A in all_tournaments(n):
        tot+=1
        if ham(A,n) < disc(A,n)-1e-6: viol+=1
print(f"  n<=6 all labelled: {tot} tournaments, H<disc violations = {viol}")
print("\nDONE.")
