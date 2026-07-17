#!/usr/bin/env python3
"""
n9_wild_hunt_klein_S317.py — klein-2026-07-16-S317

THE n=9 INVISIBLE CENSUS + THE WILD-PAIR HUNT (owner directive).

Extends THM-925's method one rung: all 191,536 n=9 iso classes via 6880 x 256 =
1,761,280 extension tournaments; exact extended panels (cpA, cpL, H, tau_in, tau_out);
panel buckets; EXACT canonical forms (individualize-and-refine, parallelized) inside
buckets.  Classification of every invisible pair: CONE (sink present, base pair tied),
COCONE (source present, dual), TOWER (both), or WILD (neither sink nor source in some
member — a native tie the laundering principle cannot manufacture).

Predictions computed from the n=8 panels en route:
  one-eyed cones  = #equal-H (cpA,cpL)-tied n=8 pairs (tau_in collapse is free),
  one-eyed cocones = #one-eyed-invisible n=8 pairs (27),
  extended cones  = #(cpA,cpL,H,tau_out)-tied n=8 pairs,
  extended cocones = #(cpA,cpL,H,tau_in)-tied n=8 pairs (27; dual laundering).
Any found pair with a claim of wildness gets EXACT re-verification (Fraction
Faddeev-LeVerrier + Bareiss taus) before being reported.
"""
import itertools, time, sys, hashlib, os
import numpy as np
from fractions import Fraction as Fr

T0=time.time()
def log(m): print(f"[{time.time()-T0:8.1f}s] {m}", flush=True)

# ---------------- shared exact tools ----------------
def charpoly_frac(M):
    n=len(M)
    Mf=[[Fr(M[i][j]) for j in range(n)] for i in range(n)]
    c=[Fr(1)]; N=[[Fr(1) if i==j else Fr(0) for j in range(n)] for i in range(n)]
    for k in range(1,n+1):
        MN=[[sum(Mf[i][l]*N[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        tr=sum(MN[i][i] for i in range(n)); ck=-tr/k; c.append(ck)
        N=[[MN[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return tuple(int(x) for x in c)

def det_bareiss(M):
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

def exact_panel(A,n):
    cpA=charpoly_frac(A)
    s=[sum(A[u]) for u in range(n)]
    L=[[(s[u] if u==v else 0)-A[u][v] for v in range(n)] for u in range(n)]
    cpL=charpoly_frac(L)
    full=1<<n
    dp=[[0]*n for _ in range(full)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(full):
        for v in range(n):
            if dp[m][v]:
                for u in range(n):
                    if not (m>>u)&1 and A[v][u]:
                        dp[m|(1<<u)][u]+=dp[m][v]
    H=sum(dp[full-1][v] for v in range(n))
    def tauv(Lm):
        return tuple(sorted(det_bareiss([[Lm[u][v] for v in range(n) if v!=r]
                                         for u in range(n) if u!=r]) for r in range(n)))
    tin=tauv(L)
    Ar=[[A[v][u] for v in range(n)] for u in range(n)]
    sr=[sum(Ar[u]) for u in range(n)]
    Lr=[[(sr[u] if u==v else 0)-Ar[u][v] for v in range(n)] for u in range(n)]
    tout=tauv(Lr)
    return cpA,cpL,H,tin,tout

# ---------------- vectorized batch panel machinery ----------------
def charpoly_batch(M,n):
    N=M.shape[0]
    I=np.eye(n,dtype=np.int64)[None,:,:]
    coeffs=np.zeros((N,n+1),dtype=np.int64); coeffs[:,0]=1
    Nk=np.repeat(I,N,axis=0)
    for k in range(1,n+1):
        MN=np.matmul(M,Nk)
        tr=np.trace(MN,axis1=1,axis2=2)
        assert (tr%k==0).all()
        ck=-(tr//k); coeffs[:,k]=ck
        Nk=MN+ck[:,None,None]*I
    return coeffs

def H_batch(A,n):
    NT=A.shape[0]; out=np.zeros(NT,dtype=np.int64)
    full=1<<n
    trans=[(m,v,u) for m in range(full) for v in range(n) if (m>>v)&1
           for u in range(n) if not (m>>u)&1]
    SB=8192
    for b0 in range(0,NT,SB):
        b1=min(NT,b0+SB); B=b1-b0
        Ab=A[b0:b1]
        dp=np.zeros((B,full,n),dtype=np.int64)
        for v in range(n): dp[:,1<<v,v]=1
        for m,v,u in trans:
            dp[:,m|(1<<u),u]+=dp[:,m,v]*Ab[:,v,u]
        out[b0:b1]=dp[:,full-1,:].sum(axis=1)
    return out

def tau_batch(L,n):
    outs=[]; idx=np.arange(n)
    for r in range(n):
        keep=idx[idx!=r]
        minor=L[:,keep][:,:,keep].astype(np.float64)
        outs.append(np.rint(np.linalg.det(minor)).astype(np.int64))
    T=np.stack(outs,axis=1); T.sort(axis=1)
    return T

def panels_for_stack(A,n):
    cpA=charpoly_batch(A,n)
    s=A.sum(axis=2)
    L=-A.copy(); L[:,np.arange(n),np.arange(n)]=s
    cpL=charpoly_batch(L,n)
    H=H_batch(A,n)
    tin=tau_batch(L,n)
    Ar=np.transpose(A,(0,2,1)).copy()
    sr=Ar.sum(axis=2)
    Lr=-Ar; Lr[:,np.arange(n),np.arange(n)]=sr
    tout=tau_batch(Lr,n)
    return cpA,cpL,H,tin,tout,s

# ---------------- canonical labeling (individualize & refine) ----------------
def refine(rows,labels,n):
    for _ in range(3):
        new=[]
        for v in range(n):
            rv=rows[v]
            outs=sorted(labels[u] for u in range(n) if (rv>>u)&1)
            ins=sorted(labels[u] for u in range(n) if (rows[u]>>v)&1)
            new.append((labels[v],tuple(outs),tuple(ins)))
        comp={x:i for i,x in enumerate(sorted(set(new)))}
        nl=tuple(comp[x] for x in new)
        if nl==labels: return nl
        labels=nl
    return labels

def canon(rows,n):
    best=[None]
    def rec(labels,depth):
        cells={}
        for v,l in enumerate(labels): cells.setdefault(l,[]).append(v)
        nonsing=[l for l in sorted(cells) if len(cells[l])>1]
        if not nonsing:
            order=sorted(range(n),key=lambda v:labels[v])
            enc=0; bit=0
            for i in range(n):
                oi=rows[order[i]]
                for j in range(i+1,n):
                    if (oi>>order[j])&1: enc|=(1<<bit)
                    bit+=1
            if best[0] is None or enc<best[0]: best[0]=enc
            return
        for v in cells[nonsing[0]]:
            nl=list(labels); nl[v]=-1-depth
            rec(refine(rows,tuple(nl),n),depth+1)
    rec(refine(rows,(0,)*n,n),0)
    return best[0]

def canon_from_A(A,n):
    rows=tuple(int(sum((1<<v) if A[u][v] else 0 for v in range(n))) for u in range(n))
    return canon(rows,n)

# ---------------- n=7 -> n=8 class reps ----------------
def census7():
    n=7; m=21
    pairs=list(itertools.combinations(range(n),2))
    pidx={pr:i for i,pr in enumerate(pairs)}
    remaps=[[(pidx[(min(g[u],g[v]),max(g[u],g[v]))],0 if g[u]<g[v] else 1)
             for (u,v) in pairs] for g in itertools.permutations(range(n))]
    seen=bytearray(1<<m); reps=[]
    for bits in range(1<<m):
        if seen[bits]: continue
        orb=set()
        for tab in remaps:
            out=0
            for i in range(m):
                b=(bits>>i)&1; j,fl=tab[i]; out|=((b^fl)<<j)
            orb.add(out)
        for t in orb: seen[t]=1
        reps.append(bits)
    return reps,pairs

def build8():
    """returns list of 6880 n=8 class reps as row-bitmask tuples + their panel arrays"""
    reps7,pairs7=census7()
    log(f"n=7 classes: {len(reps7)}")
    A_all=np.zeros((len(reps7)*128,8,8),dtype=np.int64)
    for ci,bits in enumerate(reps7):
        A7=np.zeros((7,7),dtype=np.int64)
        for i,(u,v) in enumerate(pairs7):
            if (bits>>i)&1: A7[u,v]=1
            else: A7[v,u]=1
        base=np.zeros((8,8),dtype=np.int64); base[:7,:7]=A7
        for p in range(128):
            M=base.copy()
            for u in range(7):
                if (p>>u)&1: M[u,7]=1
                else: M[7,u]=1
            A_all[ci*128+p]=M
    cpA,cpL,H,tin,tout,s=panels_for_stack(A_all,8)
    log("n=8 panels done")
    buckets={}
    for t in range(A_all.shape[0]):
        key=hashlib.blake2b(cpA[t].tobytes()+cpL[t].tobytes()+H[t].tobytes()
                            +tin[t].tobytes()+tout[t].tobytes(),digest_size=16).digest()
        buckets.setdefault(key,[]).append(t)
    reps8=[]; rep_panel_idx=[]
    for key,mem in buckets.items():
        forms={}
        for t in mem:
            cf=canon_from_A(A_all[t],8)
            if cf not in forms: forms[cf]=t
        for cf,t in forms.items():
            reps8.append(cf); rep_panel_idx.append(t)
    log(f"n=8 classes: {len(reps8)}")
    assert len(reps8)==6880
    # panels per class rep (for predictions)
    ridx=np.array(rep_panel_idx)
    return reps8, (cpA[ridx],cpL[ridx],H[ridx],tin[ridx],tout[ridx],s[ridx])

def n8_predictions(P8):
    cpA,cpL,H,tin,tout,s=P8
    NC=cpA.shape[0]
    def group(keys):
        d={}
        for i in range(NC):
            d.setdefault(keys[i],[]).append(i)
        return d
    kA=[cpA[i].tobytes() for i in range(NC)]
    kL=[cpL[i].tobytes() for i in range(NC)]
    kH=[int(H[i]) for i in range(NC)]
    kti=[tin[i].tobytes() for i in range(NC)]
    kto=[tout[i].tobytes() for i in range(NC)]
    def pairs_tied(idx_keys):
        d={}
        for i in range(NC): d.setdefault(idx_keys(i),[]).append(i)
        return sum(len(g)*(len(g)-1)//2 for g in d.values() if len(g)>1)
    pred={}
    pred['cones_oneeyed']=pairs_tied(lambda i:(kA[i],kL[i],kH[i]))
    pred['cocones_oneeyed']=pairs_tied(lambda i:(kA[i],kL[i],kH[i],kti[i]))
    pred['cones_ext']=pairs_tied(lambda i:(kA[i],kL[i],kH[i],kto[i]))
    pred['cocones_ext']=pairs_tied(lambda i:(kA[i],kL[i],kH[i],kti[i]))
    pred['both_ext']=pairs_tied(lambda i:(kA[i],kL[i],kH[i],kti[i],kto[i]))
    return pred

# ---------------- main ----------------
def main():
    reps8,P8=build8()
    pred=n8_predictions(P8)
    log(f"n=8 tie-pair predictions: {pred}")
    log("  => predicted n=9 one-eyed invisible >= cones "
        f"{pred['cones_oneeyed']} + cocones {pred['cocones_oneeyed']} (overlap = both-tied)")
    log(f"  => predicted n=9 extended invisible >= cones(tau_out-tied) {pred['cones_ext']}"
        f" + cocones(tau_in-tied) {pred['cocones_ext']} - both {pred['both_ext']}")

    # rebuild reps8 as matrices
    reps8_A=np.zeros((6880,9,9),dtype=np.int64)
    for ci,enc in enumerate(reps8):
        bit=0
        for i in range(8):
            for j in range(i+1,8):
                if (enc>>bit)&1: reps8_A[ci,i,j]=1
                else: reps8_A[ci,j,i]=1
                bit+=1
    log("extension sweep over 6880 x 256 = 1,761,280 ...")
    ext_buckets={}; one_buckets={}
    sinkflags=np.zeros(6880*256,dtype=bool); srcflags=np.zeros(6880*256,dtype=bool)
    CH=256   # classes per chunk
    for c0 in range(0,6880,CH):
        c1=min(6880,c0+CH); NC=c1-c0
        A=np.zeros((NC*256,9,9),dtype=np.int64)
        for k,ci in enumerate(range(c0,c1)):
            base=reps8_A[ci]
            for p in range(256):
                M=base.copy()
                for u in range(8):
                    if (p>>u)&1: M[u,8]=1
                    else: M[8,u]=1
                A[k*256+p]=M
        cpA,cpL,H,tin,tout,s=panels_for_stack(A,9)
        gid0=c0*256
        sinkflags[gid0:gid0+NC*256]=(s==0).any(axis=1)
        srcflags[gid0:gid0+NC*256]=(s==8).any(axis=1)
        for t in range(NC*256):
            base_b=cpA[t].tobytes()+cpL[t].tobytes()+H[t].tobytes()+tin[t].tobytes()
            k1=hashlib.blake2b(base_b,digest_size=16).digest()
            k2=hashlib.blake2b(base_b+tout[t].tobytes(),digest_size=16).digest()
            one_buckets.setdefault(k1,[]).append(gid0+t)
            ext_buckets.setdefault(k2,[]).append(gid0+t)
        if (c0//CH)%4==0:
            log(f"  chunk {c0//CH+1}/{(6880+CH-1)//CH}: buckets ext={len(ext_buckets)} one={len(one_buckets)}")
    log(f"panel sweep done: ext buckets {len(ext_buckets)}, one-eyed {len(one_buckets)}")

    # exact validation of float taus on a random sample
    rng=np.random.default_rng(11)
    for gid in rng.integers(0,6880*256,25):
        ci,p=divmod(int(gid),256)
        M=reps8_A[ci].copy()
        for u in range(8):
            if (p>>u)&1: M[u,8]=1
            else: M[8,u]=1
        Alist=[[int(M[i,j]) for j in range(9)] for i in range(9)]
        ep=exact_panel(Alist,9)
        # recompute batch values for this single tournament
        cpA1,cpL1,H1,tin1,tout1,_=panels_for_stack(M[None,:,:],9)
        assert tuple(cpA1[0])==ep[0] and tuple(cpL1[0])==ep[1] and int(H1[0])==ep[2]
        assert tuple(tin1[0])==ep[3] and tuple(tout1[0])==ep[4]
    log("exact spot-validation of batch panels OK (25 samples)")

    # canonicalize within buckets (the ext buckets partition everything)
    canon_cache={}
    def canon_gid(gid):
        cf=canon_cache.get(gid)
        if cf is None:
            ci,p=divmod(gid,256)
            M=reps8_A[ci].copy()
            for u in range(8):
                if (p>>u)&1: M[u,8]=1
                else: M[8,u]=1
            rows=tuple(int(sum((1<<v) if M[uu,v] else 0 for v in range(9))) for uu in range(9))
            cf=canon(rows,9)
            canon_cache[gid]=cf
        return cf
    def dedupe(buckets,name):
        total=0; groups=[]
        done=0
        for key,mem in buckets.items():
            forms={}
            for gid in mem:
                cf=canon_gid(gid)
                if cf not in forms: forms[cf]=gid
            total+=len(forms)
            if len(forms)>1: groups.append((key,forms))
            done+=1
            if done%20000==0: log(f"  {name}: {done}/{len(buckets)} buckets, classes so far {total}")
        return total,groups
    log("canonicalizing (extended panel buckets = full partition)...")
    total_ext,groups_ext=dedupe(ext_buckets,"ext")
    log(f"TOTAL n=9 classes: {total_ext} (expect 191536)")
    total_one,groups_one=dedupe(one_buckets,"one")
    log(f"one-eyed grouping total (consistency): {total_one}")

    def npairs(groups): return sum(len(f)*(len(f)-1)//2 for _,f in groups)
    print()
    print("="*90)
    print(f"n=9 ONE-EYED panel (cpA,cpL,H,tau_in): {len(groups_one)} invisible groups, "
          f"{npairs(groups_one)} pairs")
    print(f"n=9 EXTENDED panel (+tau_out):        {len(groups_ext)} invisible groups, "
          f"{npairs(groups_ext)} pairs")
    print(f"predictions from n=8: {pred}")
    print("="*90)
    wild=[]
    for key,forms in sorted(groups_ext,key=lambda kf:-len(kf[1])):
        gids=list(forms.values())
        tags=[]
        for gid in gids:
            snk="SINK" if sinkflags[gid] else "";  src="SRC" if srcflags[gid] else ""
            tags.append((snk+"+"+src).strip("+") or "CORE")
        iswild=any(t=="CORE" for t in tags)
        if iswild: wild.append((key,forms))
        ci,p=divmod(gids[0],256)
        M=reps8_A[ci].copy()
        for u in range(8):
            if (p>>u)&1: M[u,8]=1
            else: M[8,u]=1
        sc=tuple(sorted(int(x) for x in M.sum(axis=1)))
        print(f"  group size {len(forms)}: tags={tags} scores[0]={sc} "
              f"{'*** WILD ***' if iswild else ''}")
    print()
    print(f"WILD candidate groups: {len(wild)}")
    # exact re-verification of wild candidates
    for key,forms in wild:
        gids=list(forms.values())
        panels=[]
        for gid in gids:
            ci,p=divmod(gid,256)
            M=reps8_A[ci].copy()
            for u in range(8):
                if (p>>u)&1: M[u,8]=1
                else: M[8,u]=1
            Alist=[[int(M[i,j]) for j in range(9)] for i in range(9)]
            panels.append(exact_panel(Alist,9))
        allsame=all(p==panels[0] for p in panels)
        print(f"  WILD group exact re-verify: panels identical = {allsame}; "
              f"members gid={gids}")
        if allsame:
            for gid in gids:
                ci,p=divmod(gid,256)
                print(f"    member: class8={ci} pattern={p:08b}")
    log("DONE")

if __name__=="__main__":
    main()
