#!/usr/bin/env python3
"""
permanental_determines_H_hunt_monad.py  (monad-explorer-2026-06-15-S6)  [pure-python]

SHARP TEST of "does the pair (char poly, perm poly) DETERMINE H?".  The random (char,perm)
bucketing in permanental_companion_monad.py [5] is WEAK at n>=8 (collisions on the full
(char,perm) key are rare).  Here we bucket by CHAR POLY first (cospectral classes are well
populated), then within each cospectral class ask whether perm poly determines H.

Predicted break points (carriers that (char,perm) is forecast to MERGE):
   n=8:  D44 (disjoint 4-cycle pairs) merged with D35 inside E_8 = D44+D35;
         within a (char,perm) class  H = const - 4*D44, so (char,perm)->H iff D44 fixed.
   n=9:  T333 (disjoint triangle triples) merged with c9 inside O_9 = c9+T333;
         within a (char,perm) class  H = const + 6*T333, so (char,perm)->H iff T333 fixed.

FAST primitives: char poly (Faddeev int), perm poly (Ryser over Z[x]), H (bitmask DP).
Diagnostic carriers (D44,D35,c9,T333) computed by packing enumeration only on a witness.
"""
import sys, random
from fractions import Fraction

def random_tournament(n, rng):
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if rng.randint(0,1): A[i][j]=1
            else: A[j][i]=1
    return A

def matmul(A,B,n):
    C=[[0]*n for _ in range(n)]
    for i in range(n):
        Ai=A[i];Ci=C[i]
        for k in range(n):
            a=Ai[k]
            if a:
                Bk=B[k]
                for j in range(n): Ci[j]+=a*Bk[j]
    return C

def charpoly_int(A,n):
    M=[[1 if i==j else 0 for j in range(n)] for i in range(n)]; co=[1]
    for k in range(1,n+1):
        AM=matmul(A,M,n); tr=sum(AM[i][i] for i in range(n))
        ck=Fraction(-tr,k).numerator; co.append(ck)
        if k<n: M=[[AM[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return tuple(co)

def polymul(p,q):
    r=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        if a:
            for j,b in enumerate(q): r[i+j]+=a*b
    return r

def permpoly_int(A,n):
    """per(xI+A): returns tuple e_m^unsigned for m=0..n (coeff of x^{n-m}).
    Ryser:  per(M)=sum_S (-1)^{n-|S|} prod_i rowsum_i(S),  M=xI+A.
    rowsum_i(S) = [i in S]*x + a_i(S),  a_i(S)=popcount(out-neighbours of i in S).
    Term = sign * (prod_{i not in S} a_i(S)) * prod_{i in S}(x + a_i(S))  [0 if any out-row empty]."""
    rowmask=[0]*n
    for i in range(n):
        Ai=A[i]; m=0
        for j in range(n):
            if Ai[j]: m|=(1<<j)
        rowmask[i]=m
    full=[0]*(n+1)
    for Sbits in range(1<<n):
        scalar=1; poly=[1]; dead=False
        for i in range(n):
            ai=(rowmask[i]&Sbits).bit_count()
            if (Sbits>>i)&1:
                poly=polymul(poly,[ai,1])
            else:
                if ai==0: dead=True; break
                scalar*=ai
        if dead: continue
        sign=1 if (n-bin(Sbits).count('1'))%2==0 else -1
        sc=sign*scalar
        for t,c in enumerate(poly): full[t]+=sc*c
    return tuple(full[n-m] for m in range(n+1))

def count_ham_paths(A,n):
    full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        row=dp[mask]
        if not any(row): continue
        for v in range(n):
            cv=row[v]
            if not cv: continue
            Av=A[v]
            for w in range(n):
                if not (mask>>w)&1 and Av[w]: dp[mask|(1<<w)][w]+=cv
    return sum(dp[full][v] for v in range(n))

def carriers(A,n):
    """Named disjoint-packing carriers via cycle/packing enumeration (only on witnesses)."""
    cyc=[]
    for start in range(n):
        path=[start]; vis={start}
        def dfs(u):
            for w in range(n):
                if w==start:
                    if len(path)>=3 and A[u][start]: cyc.append((len(path),frozenset(path)))
                elif w>start and w not in vis and A[u][w]:
                    vis.add(w);path.append(w);dfs(w);path.pop();vis.discard(w)
        dfs(start)
    vs=[s for (_,s) in cyc]; ln=[L for (L,_) in cyc]; nc=len(cyc); Nlam={}
    def rec(start,used,lam):
        Nlam[lam]=Nlam.get(lam,0)+1
        for i in range(start,nc):
            if not (vs[i]&used): rec(i+1, used|vs[i], tuple(sorted(lam+(ln[i],))))
    rec(0,frozenset(),())
    return {'c6':Nlam.get((6,),0),'c7':Nlam.get((7,),0),'c8':Nlam.get((8,),0),'c9':Nlam.get((9,),0),
            'D33':Nlam.get((3,3),0),'D35':Nlam.get((3,5),0),'D44':Nlam.get((4,4),0),
            'T333':Nlam.get((3,3,3),0)}

def hunt(n,n_samples,seed=101):
    rng=random.Random(seed); by_char={}
    for _ in range(n_samples):
        A=random_tournament(n,rng)
        cp=charpoly_int(A,n)
        by_char.setdefault(cp,[]).append((permpoly_int(A,n), count_ham_paths(A,n), A))
    cospectral=[(cp,recs) for cp,recs in by_char.items() if len(recs)>=2]
    H_split=0; perm_fail=0; witness=None
    for cp,recs in cospectral:
        if len({r[1] for r in recs})>1: H_split+=1
        permmap={}
        for pp,H,A in recs:
            if pp in permmap and permmap[pp][0]!=H:
                perm_fail+=1
                if witness is None: witness=(permmap[pp][1], A, permmap[pp][0], H)
            else:
                permmap.setdefault(pp,(H,A))
    print(f" n={n}: {len(by_char)} cospectral classes ({len(cospectral)} with >=2 members) "
          f"from {n_samples} samples")
    print(f"     H non-spectral (class splits H):      {H_split}")
    print(f"     (char,perm) FAILS to determine H:     {perm_fail}   "
          f"-> {'DETERMINES H' if perm_fail==0 else 'does NOT determine H'}")
    if witness:
        A1,A2,H1,H2=witness
        c1=carriers(A1,n); c2=carriers(A2,n)
        print(f"     WITNESS (same char, same perm, different H): H={H1} vs H={H2}")
        print(f"        carriers T1: {{{', '.join(f'{k}={v}' for k,v in c1.items() if v)}}}")
        print(f"        carriers T2: {{{', '.join(f'{k}={v}' for k,v in c2.items() if v)}}}")
        diff={k:(c1[k],c2[k]) for k in c1 if c1[k]!=c2[k]}
        print(f"        DIFFERING carriers: {diff}")
    return perm_fail

if __name__=='__main__':
    print("="*82)
    print(" SHARP HUNT: does (char poly, perm poly) determine H?  (bucket by char first)")
    print("="*82)
    # optional argv: n  n_samples   (else default battery)
    if len(sys.argv)>=2:
        n=int(sys.argv[1]); ns=int(sys.argv[2]) if len(sys.argv)>=3 else 150000
        hunt(n,ns)
    else:
        for n,ns in [(7,40000),(8,150000),(9,120000)]:
            hunt(n,ns); sys.stdout.flush()
    print("DONE.")
