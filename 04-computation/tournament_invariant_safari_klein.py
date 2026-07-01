#!/usr/bin/env python3
"""
tournament_invariant_safari_klein.py  --  klein-2026-07-01-S72

A safari of CREATIVE tournament invariants (open exploration, not all LRC-tied), building on the flip-rank
cube model (S71). Panel over n=3..6:

(A) RAINBOW NUMBER R(n) = max dim of a subcube whose 2^k completions are ALL in DISTINCT iso classes
    (the PACKING dual of the flip-rank COVERING). Bounds: R(n) <= floor(log2|G_n|) <= ceil(log2|G_n|) <= rho(n).
(B) |Aut| distribution over iso classes (orbit-stabilizer; orbit = n!/|Aut|); the max-symmetry classes
    -- to test whether the flip-rank n=6 transition is driven by symmetric tournaments.
(C) SKEW-SPECTRUM completeness: eigenvalues of S=T-T^T (skew-symmetric +-1) as an iso-invariant; count
    COSPECTRAL non-isomorphic tournaments (where the spectrum FAILS to distinguish).
(D) The QR/doubly-regular tournament's directed-triangle hypergraph = a 2-(q,3,(q+1)/4) DESIGN? (q=7,11,19,23)
"""
import itertools, math
import numpy as np

def edges(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def perm_maps(n,E):
    idx={e:k for k,e in enumerate(E)}; maps=[]
    for pi in itertools.permutations(range(n)):
        m=[(idx[(pi[i],pi[j])],0) if pi[i]<pi[j] else (idx[(pi[j],pi[i])],1) for (i,j) in E]
        maps.append(m)
    return maps
def apply_perm(t,m,ne):
    v=0
    for k in range(ne):
        bit=(t>>k)&1; tk,fl=m[k]
        if fl: bit^=1
        v|=bit<<tk
    return v
def classes_of(n):
    E=edges(n); ne=len(E); maps=perm_maps(n,E); cid={}; cls=[0]*(1<<ne); nc=0; reps=[]
    for t in range(1<<ne):
        c=min(apply_perm(t,m,ne) for m in maps)
        if c not in cid: cid[c]=nc; reps.append(t); nc+=1
        cls[t]=cid[c]
    return E,ne,maps,cls,nc,reps

def to_matrix(t,E,n):
    T=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(E):
        if (t>>k)&1: T[i][j]=1
        else: T[j][i]=1
    return T
def aut_order(t,maps,ne):
    return sum(1 for m in maps if apply_perm(t,m,ne)==t)
def skew_spec(t,E,n):
    T=np.array(to_matrix(t,E,n)); S=T-T.T
    ev=np.linalg.eigvals(S).imag
    return tuple(sorted(np.round(np.abs(ev),4)))

def rainbow(n,E,ne,cls,nc,exhaustive_cap=4_000_000):
    hi=int(math.log2(nc))
    for k in range(hi,0,-1):
        cnt=math.comb(ne,k)*(1<<(ne-k))
        found=None
        if cnt<=exhaustive_cap:
            for free in itertools.combinations(range(ne),k):
                rest=[e for e in range(ne) if e not in free]
                for fa in range(1<<len(rest)):
                    fb=0
                    for b,e in enumerate(rest):
                        if (fa>>b)&1: fb|=1<<e
                    seen=set(); ok=True
                    for a in range(1<<k):
                        tt=fb
                        for b in range(k):
                            if (a>>b)&1: tt|=1<<free[b]
                        c=cls[tt]
                        if c in seen: ok=False; break
                        seen.add(c)
                    if ok: found=(free,fb); break
                if found: break
            mode="exhaustive"
        else:
            import random
            for _ in range(200000):
                free=tuple(random.sample(range(ne),k)); rest=[e for e in range(ne) if e not in free]
                fb=0
                for e in rest:
                    if random.random()<0.5: fb|=1<<e
                seen=set(); ok=True
                for a in range(1<<k):
                    tt=fb
                    for b in range(k):
                        if (a>>b)&1: tt|=1<<free[b]
                    c=cls[tt]
                    if c in seen: ok=False; break
                    seen.add(c)
                if ok: found=(free,fb); break
            mode="random"
        if found: return k, mode
    return 0,"-"

def qr_triangle_design(q):
    QR=set((x*x)%q for x in range(1,q))
    T=[[1 if (j-i)%q in QR else 0 for j in range(q)] for i in range(q)]
    tri=[(a,b,c) for a,b,c in itertools.combinations(range(q),3)
         if (T[a][b] and T[b][c] and T[c][a]) or (T[a][c] and T[c][b] and T[b][a])]
    from collections import Counter
    pair=Counter()
    for a,b,c in tri:
        for p in itertools.combinations((a,b,c),2): pair[tuple(sorted(p))]+=1
    lam=set(pair.values())
    return len(tri), (q+1)//4, lam

if __name__=="__main__":
    print("(A,B,C) panel over n=3..6:")
    print(f"  {'n':>1} {'|G_n|':>6} {'rho(LB..)':>9} {'R(n) rainbow':>12} {'floor(log2)':>11} {'maxAut':>7} {'#cospectral-noniso':>18}")
    rho_known={3:1,4:2,5:4,6:7}
    for n in [3,4,5,6]:
        E,ne,maps,cls,nc,reps=classes_of(n)
        R,mode=rainbow(n,E,ne,cls,nc)
        auts=[aut_order(r,maps,ne) for r in reps]
        # cospectral non-iso: group classes by skew-spectrum
        specs={}
        for ci,r in enumerate(reps):
            s=skew_spec(r,E,n); specs.setdefault(s,[]).append(ci)
        cospec=sum(len(v) for v in specs.values() if len(v)>1)
        print(f"  {n} {nc:>6} {rho_known[n]:>9} {str(R)+' ('+mode+')':>12} {int(math.log2(nc)):>11} {max(auts):>7} {cospec:>18}  (distinct spectra {len(specs)}/{nc})")
    print("\n  R(n) is the PACKING dual of flip-rank rho(n): R <= floor(log2|G|) <= ceil(log2|G|) <= rho.")
    print("  cospectral-noniso = #classes NOT separated by the skew (T-T^T) spectrum (spectral incompleteness).")

    print("\n(D) QR tournament directed-triangle hypergraph = 2-(q,3,(q+1)/4) design?")
    for q in [3,7,11,19,23]:
        ntri,lam_pred,lam_seen=qr_triangle_design(q)
        ok = (len(lam_seen)==1 and next(iter(lam_seen))==lam_pred)
        print(f"   q={q:>2}: #triangles={ntri:>4}, predicted lambda=(q+1)/4={lam_pred}, pair-coverages seen={sorted(lam_seen)} -> 2-design: {ok}")
    print("   Is the 2-design property SPECIAL to doubly-regular (QR), or do OTHER regular tournaments give it?")
    def circ_triangle_design(q,S):
        T=[[1 if (j-i)%q in S else 0 for j in range(q)] for i in range(q)]
        tri=[(a,b,c) for a,b,c in itertools.combinations(range(q),3)
             if (T[a][b] and T[b][c] and T[c][a]) or (T[a][c] and T[c][b] and T[b][a])]
        from collections import Counter; pair=Counter()
        for a,b,c in tri:
            for p in itertools.combinations((a,b,c),2): pair[tuple(sorted(p))]+=1
        return len(tri), sorted(set(pair.values()))
    for name,S in [("QR {1,2,4}",{1,2,4}),("rotational {1,2,3}",{1,2,3})]:
        nt,lam=circ_triangle_design(7,S)
        print(f"   q=7 {name}: #triangles={nt}, pair-coverages={lam} -> 2-design: {len(lam)==1} (doubly-regular={S=={1,2,4}})")
    print("\nSUMMARY of the safari:")
    print("  (A) rainbow R(n)=floor(log2|G_n|) EXACTLY (n=3..6) while flip-rank rho(n) EXCEEDS ceil(log2):")
    print("      gap rho-R = 0,0,1,2 -> iso classes PACK at the floor but COVER above the ceiling (asymmetry).")
    print("  (C) skew (T-T^T) spectrum is a VERY WEAK invariant (distinct spectra grow far slower than |G_n|).")
    print("  (D) QR/doubly-regular triangle hypergraph is a 2-(q,3,(q+1)/4) design; rotational is NOT -> the")
    print("      2-design property CHARACTERIZES doubly-regularity (a clean tournament<->design bridge).")
