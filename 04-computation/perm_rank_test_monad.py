#!/usr/bin/env python3
"""
perm_rank_test_monad.py  (monad-explorer-2026-06-15-S6)  [pure-python, exact Q-rank]

Decisive, data-efficient test of "does (char poly, perm poly) determine H?", using the
project's WITHIN-COSPECTRAL-CLASS RANK method (robust to tiny classes; no need for exact
(char,perm) collisions).  Within a cospectral class the char poly (e_signed) is fixed, so
H and the perm-poly coordinates e_unsigned are both linear in the non-spectral CARRIERS.
(char,perm) determines H  <=>  H is an affine function of e_unsigned within classes
                            <=>  ker(carriers->e_unsigned) subset ker(carriers->H)
                            <=>  rank[ De_unsigned | DH ] == rank[ De_unsigned ]
where D* are within-class deltas vs a class reference.

If rank grows by adding DH, there is a within-class carrier direction that MOVES H while
keeping e_unsigned (hence the perm poly) FIXED  =>  (char,perm) does NOT determine H, and
we exhibit the offending carrier (predicted: D44 at n=8, T333 at n=9).
"""
import sys, random
from fractions import Fraction

def random_tournament(n,rng):
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
    rowmask=[0]*n
    for i in range(n):
        Ai=A[i];m=0
        for j in range(n):
            if Ai[j]: m|=(1<<j)
        rowmask[i]=m
    full=[0]*(n+1)
    for Sbits in range(1<<n):
        scalar=1;poly=[1];dead=False
        for i in range(n):
            ai=(rowmask[i]&Sbits).bit_count()
            if (Sbits>>i)&1: poly=polymul(poly,[ai,1])
            else:
                if ai==0: dead=True;break
                scalar*=ai
        if dead: continue
        sign=1 if (n-Sbits.bit_count())%2==0 else -1
        sc=sign*scalar
        for t,c in enumerate(poly): full[t]+=sc*c
    return tuple(full[n-m] for m in range(n+1))

def count_ham_paths(A,n):
    full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]
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
    cyc=[]
    for start in range(n):
        path=[start];vis={start}
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
            if not (vs[i]&used): rec(i+1,used|vs[i],tuple(sorted(lam+(ln[i],))))
    rec(0,frozenset(),())
    return {'c3':Nlam.get((3,),0),'c4':Nlam.get((4,),0),'c5':Nlam.get((5,),0),
            'c6':Nlam.get((6,),0),'c7':Nlam.get((7,),0),'c8':Nlam.get((8,),0),'c9':Nlam.get((9,),0),
            'D33':Nlam.get((3,3),0),'D34':Nlam.get((3,4),0),'D35':Nlam.get((3,5),0),
            'D44':Nlam.get((4,4),0),'D45':Nlam.get((4,5),0),'D36':Nlam.get((3,6),0),
            'T333':Nlam.get((3,3,3),0)}

def qrank(rows):
    """exact rank over Q of a list of lists (Fractions/ints)."""
    R=[[Fraction(x) for x in row] for row in rows]
    if not R: return 0
    ncol=len(R[0]); rank=0; r=0
    for c in range(ncol):
        piv=None
        for i in range(r,len(R)):
            if R[i][c]!=0: piv=i;break
        if piv is None: continue
        R[r],R[piv]=R[piv],R[r]
        pv=R[r][c]
        R[r]=[x/pv for x in R[r]]
        for i in range(len(R)):
            if i!=r and R[i][c]!=0:
                f=R[i][c]; R[i]=[a-f*b for a,b in zip(R[i],R[r])]
        r+=1; rank+=1
        if r==len(R): break
    return rank

def test(n,n_samples,seed=55,do_carriers=True):
    rng=random.Random(seed); by_char={}
    for _ in range(n_samples):
        A=random_tournament(n,rng)
        cp=charpoly_int(A,n)
        by_char.setdefault(cp,[]).append(A)
    # within-class deltas
    rows_perm=[]; rows_perm_H=[]; H_split=0; nclass2=0
    car_rows=[]; carH=[]
    carnames=None
    for cp,As in by_char.items():
        if len(As)<2: continue
        nclass2+=1
        ref=As[0]
        pp_ref=permpoly_int(ref,n); H_ref=count_ham_paths(ref,n)
        cr_ref=carriers(ref,n) if do_carriers else None
        if do_carriers and carnames is None: carnames=sorted(cr_ref)
        Hs={H_ref}
        for A in As[1:]:
            pp=permpoly_int(A,n); H=count_ham_paths(A,n)
            Hs.add(H)
            dperm=[pp[m]-pp_ref[m] for m in range(n+1)]
            dH=H-H_ref
            rows_perm.append(dperm); rows_perm_H.append(dperm+[dH])
            if do_carriers:
                cr=carriers(A,n)
                car_rows.append([cr[k]-cr_ref[k] for k in carnames]); carH.append(dH)
        if len(Hs)>1: H_split+=1
    rp=qrank(rows_perm); rpH=qrank(rows_perm_H)
    print(f"\n n={n}: {len(by_char)} cospectral classes, {nclass2} with >=2 members, "
          f"{len(rows_perm)} within-class deltas; {H_split} classes split H")
    print(f"   rank[ De_unsigned ]        = {rp}")
    print(f"   rank[ De_unsigned | DH ]   = {rpH}   "
          f"-> H {'IS' if rpH==rp else 'is NOT'} affine in perm poly  "
          f"=> (char,perm) {'DETERMINES' if rpH==rp else 'does NOT determine'} H")
    # carrier diagnosis: which single carrier, added to e_unsigned, captures H?
    if rpH>rp and car_rows:
        # find carriers k with ΔH NOT explained by De_unsigned but explained by De_unsigned+Δcarrier_k
        base=rows_perm_H  # already rank rpH
        for idx,k in enumerate(carnames):
            aug=[row+[cr[idx]] for row,cr in zip(rows_perm,car_rows)]
            augH=[row+[cr[idx],dh] for row,cr,dh in zip(rows_perm,car_rows,carH)]
            if qrank(augH)==qrank(aug):
                print(f"   ... adding carrier {k:5s} to perm poly RECOVERS H "
                      f"(H affine in e_unsigned + {k})")
    # also: is the merged-carrier prediction visible? report rank of carrier space vs perm space
    return rp,rpH

if __name__=='__main__':
    print("="*82)
    print(" WITHIN-CLASS RANK TEST: is H affine in the permanental-poly coordinates?")
    print("="*82)
    if len(sys.argv)>=2:
        do_c = not (len(sys.argv)>=4 and sys.argv[3]=='nocar')
        test(int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv)>=3 else 8000, do_carriers=do_c)
    else:
        for n,ns in [(6,4000),(7,6000),(8,9000),(9,12000)]:
            test(n,ns); sys.stdout.flush()
    print("\nDONE.")
