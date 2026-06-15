#!/usr/bin/env python3
"""
permanental_determines_H_hunt_monad.py  (monad-explorer-2026-06-15-S6)  [pure-python]

SHARP TEST of "does the pair (char poly, perm poly) DETERMINE H?".  The random (char,perm)
bucketing in permanental_companion_monad.py [5] is WEAK at n>=8 (collisions on the full
(char,perm) key are rare under uniform sampling -> few within-key comparisons).  Here we
bucket by CHAR POLY first (cospectral classes are well populated), then within each
cospectral class ask whether the perm poly determines H.  We ALSO track the exact carriers
that (char,perm) is predicted to merge:
   n=8:  D44 (disjoint 4-cycle pairs) -- merged with D35 inside E_8 = D44+D35;
         H = const - 4*D44 within a (char,perm) class, so  (char,perm)->H  iff D44 fixed.
   n=9:  T333 (disjoint triangle triples) -- merged with c9 inside O_9 = c9+T333;
         H = const + 6*T333 within a (char,perm) class, so  (char,perm)->H  iff T333 fixed.

For each n we report, over cospectral classes (>=2 members observed):
   - how many split H (H non-spectral),
   - how many have perm poly that FAILS to determine H (a real (char,perm) H-split),
   - and a concrete witness with the differing carrier.
"""
import sys, random
from fractions import Fraction

def random_tournament(n, rng):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
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
        ck=Fraction(-tr,k); ck=ck.numerator; co.append(ck)
        if k<n: M=[[AM[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return tuple(co)

def all_cycles(A,n):
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
    return cyc

def carriers_and_perm(A,n):
    """Return (perm_poly e_unsigned vector, H, dict of named carriers)."""
    cyc=all_cycles(A,n)
    vs=[s for (_,s) in cyc]; ln=[L for (L,_) in cyc]; nc=len(cyc)
    eu={}; H=0
    Nlam={}
    out=[]
    def rec(start,used,k,cov,lam):
        eu[cov]=eu.get(cov,0)+1
        Nlam[lam]=Nlam.get(lam,0)+1
        if all(L%2==1 for L in lam):  # odd packing (empty lam counts: all() True)
            H+=0  # placeholder, fixed below
        for i in range(start,nc):
            if not (vs[i]&used):
                rec(i+1, used|vs[i], k+1, cov+len(vs[i]), tuple(sorted(lam+(ln[i],))))
    # recompute H cleanly with 2^k weights over odd packings
    Hval=[0]
    def rec2(start,used,k,lam_all_odd):
        if lam_all_odd: Hval[0]+=2**k
        for i in range(start,nc):
            if not (vs[i]&used):
                rec2(i+1, used|vs[i], k+1, lam_all_odd and (ln[i]%2==1))
    rec(0,frozenset(),0,0,())
    rec2(0,frozenset(),0,True)
    H=Hval[0]
    eu_vec=tuple(eu.get(m,0) for m in range(n+1))
    car={
        'c6':Nlam.get((6,),0),'c7':Nlam.get((7,),0),'c8':Nlam.get((8,),0),'c9':Nlam.get((9,),0),
        'D33':Nlam.get((3,3),0),'D35':Nlam.get((3,5),0),'D44':Nlam.get((4,4),0),
        'T333':Nlam.get((3,3,3),0),
    }
    return eu_vec, H, car

def hunt(n, n_samples, seed=101):
    rng=random.Random(seed)
    by_char={}   # char -> list of (perm, H, car)
    for _ in range(n_samples):
        A=random_tournament(n,rng)
        cp=charpoly_int(A,n)
        eu,H,car=carriers_and_perm(A,n)
        by_char.setdefault(cp,[]).append((eu,H,car))
    cospectral=[(cp,recs) for cp,recs in by_char.items() if len(recs)>=2]
    H_split=0; perm_fail=0; witness=None
    carrier_var={'D44':0,'T333':0,'D35':0,'c9':0}
    for cp,recs in cospectral:
        Hs={r[1] for r in recs}
        if len(Hs)>1: H_split+=1
        # within this cospectral class, does perm determine H?
        permmap={}
        for eu,H,car in recs:
            if eu in permmap and permmap[eu]!=H:
                perm_fail+=1
                if witness is None:
                    # find the two records
                    a=[r for r in recs if r[0]==eu]
                    witness=(cp,a)
                break
            permmap[eu]=H
        # track whether D44/T333 vary within (char,perm) sub-buckets
        sub={}
        for eu,H,car in recs:
            sub.setdefault(eu,[]).append(car)
        for eu,cars in sub.items():
            for key in carrier_var:
                if len({c[key] for c in cars})>1: carrier_var[key]+=1
    print(f" n={n}: {len(by_char)} cospectral classes ({len(cospectral)} with >=2 members) "
          f"from {n_samples} samples")
    print(f"     H non-spectral (class splits H):            {H_split}")
    print(f"     (char,perm) FAILS to determine H:           {perm_fail}   "
          f"-> {'DETERMINES H' if perm_fail==0 else 'does NOT determine H'}")
    print(f"     within (char,perm): D44 varies {carrier_var['D44']}, "
          f"D35 varies {carrier_var['D35']}, c9 varies {carrier_var['c9']}, "
          f"T333 varies {carrier_var['T333']}")
    if witness:
        cp,a=witness
        print(f"     WITNESS (same char, same perm, different H):")
        for eu,H,car in a[:2]:
            rel={k:v for k,v in car.items() if v}
            print(f"        H={H}  carriers={rel}")
    return perm_fail

if __name__=='__main__':
    print("="*82)
    print(" SHARP HUNT: does (char poly, perm poly) determine H?  (bucket by char first)")
    print("="*82)
    for n,ns in [(7,30000),(8,120000),(9,120000)]:
        hunt(n,ns)
    print("DONE.")
