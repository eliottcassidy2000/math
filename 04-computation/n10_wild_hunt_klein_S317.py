#!/usr/bin/env python3
"""
n10_wild_hunt_klein_S317.py — klein-2026-07-16-S317 (overnight job, harvest next session)

THE n=10 WILD HUNT: is every wild invisible pair a chirality pair (THM-931's question)?

Design (disk-backed, memory-light):
  1. rebuild the 191,536 n=9 reps (extension census, ~20 min);
  2. sweep all 191,536 x 512 = 98,066,432 n=10 extensions in chunks; per tournament
     compute the H-FREE panel (cpA, cpL, tau_in, tau_out) exactly (integer FL + float
     dets); append "digest gid" lines to dump files;
  3. external-sort the extended dump; scan adjacent equal digests -> buckets;
     canonicalize members (I-R); non-isomorphic panel-tied pairs get H computed
     individually — the extended panel = (cpA, cpL, H, tau_in, tau_out) as in THM-931;
  4. classify every surviving pair: SINK/SRC tags, chirality (member2 ~= member1^op),
     exact Fraction re-verification.
"""
import sys, os, time, hashlib, subprocess, itertools
sys.path.insert(0,'04-computation')
import numpy as np
from n9_wild_hunt_klein_S317 import (build8, charpoly_batch, canon, exact_panel)

T0=time.time()
def log(m): print(f"[{time.time()-T0:9.1f}s] {m}", flush=True)

SCRATCH="/private/tmp/claude-501/-Users-e-Documents-GitHub-ephrepos-math/165e91a8-4641-46c3-95dd-8a288f00d110/scratchpad"

def panels_noH(A,n):
    cpA=charpoly_batch(A,n)
    s=A.sum(axis=2)
    L=-A.copy(); L[:,np.arange(n),np.arange(n)]=s
    cpL=charpoly_batch(L,n)
    idx=np.arange(n)
    def taus(Lm):
        outs=[]
        for r in range(n):
            keep=idx[idx!=r]
            minor=Lm[:,keep][:,:,keep].astype(np.float64)
            outs.append(np.rint(np.linalg.det(minor)).astype(np.int64))
        Tm=np.stack(outs,axis=1); Tm.sort(axis=1); return Tm
    tin=taus(L)
    Ar=np.transpose(A,(0,2,1)).copy()
    sr=Ar.sum(axis=2)
    Lr=-Ar; Lr[:,np.arange(n),np.arange(n)]=sr
    tout=taus(Lr)
    return cpA,cpL,tin,tout,s

def build9():
    reps8,_=build8()
    reps8_A=np.zeros((6880,9,9),dtype=np.int64)
    for ci,enc in enumerate(reps8):
        bit=0
        for i in range(8):
            for j in range(i+1,8):
                if (enc>>bit)&1: reps8_A[ci,i,j]=1
                else: reps8_A[ci,j,i]=1
                bit+=1
    log("building n=9 reps (panel buckets + canon)...")
    buckets={}
    CH=256
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
        cpA,cpL,tin,tout,s=panels_noH(A,9)
        gid0=c0*256
        for t in range(NC*256):
            key=hashlib.blake2b(cpA[t].tobytes()+cpL[t].tobytes()+tin[t].tobytes()
                                +tout[t].tobytes(),digest_size=12).digest()
            buckets.setdefault(key,[]).append(gid0+t)
    log(f"n=9 panel buckets: {len(buckets)}")
    reps9=[]
    for key,mem in buckets.items():
        forms=set()
        for gid in mem:
            ci,p=divmod(gid,256)
            M=reps8_A[ci].copy()
            for u in range(8):
                if (p>>u)&1: M[u,8]=1
                else: M[8,u]=1
            rows=tuple(int(sum((1<<v) if M[uu,v] else 0 for v in range(9))) for uu in range(9))
            forms.add(canon(rows,9))
        reps9.extend(forms)
    log(f"n=9 classes: {len(reps9)}")
    assert len(reps9)==191536
    return sorted(reps9)

def main():
    reps9=build9()
    np.save(f"{SCRATCH}/reps9.npy",np.array(reps9,dtype=np.int64))
    def rep9_matrix(enc):
        M=np.zeros((10,10),dtype=np.int64)
        bit=0
        for i in range(9):
            for j in range(i+1,9):
                if (enc>>bit)&1: M[i,j]=1
                else: M[j,i]=1
                bit+=1
        return M
    NCLS=191536; PAT=512
    dumpE=open(f"{SCRATCH}/n10_ext.dump","w")
    CH=128
    nch=(NCLS+CH-1)//CH
    for c0 in range(0,NCLS,CH):
        c1=min(NCLS,c0+CH); NC=c1-c0
        A=np.zeros((NC*PAT,10,10),dtype=np.int64)
        for k,ci in enumerate(range(c0,c1)):
            base=rep9_matrix(reps9[ci])
            for p in range(PAT):
                M=base.copy()
                for u in range(9):
                    if (p>>u)&1: M[u,9]=1
                    else: M[9,u]=1
                A[k*PAT+p]=M
        cpA,cpL,tin,tout,s=panels_noH(A,10)
        gid0=c0*PAT
        for t in range(NC*PAT):
            base_b=cpA[t].tobytes()+cpL[t].tobytes()+tin[t].tobytes()
            k2=hashlib.blake2b(base_b+tout[t].tobytes(),digest_size=12).hexdigest()
            dumpE.write(f"{k2} {gid0+t}\n")
        if (c0//CH)%25==0:
            log(f"chunk {c0//CH+1}/{nch}")
    dumpE.close()
    log("sweep done; sorting dump...")
    subprocess.run(["sort","-o",f"{SCRATCH}/n10_ext.sorted",f"{SCRATCH}/n10_ext.dump"],check=True)
    log("ext sorted")

    def H_of(A,n=10):
        full=1<<n
        dp=[[0]*n for _ in range(full)]
        for v in range(n): dp[1<<v][v]=1
        for m in range(full):
            row=dp[m]
            for v in range(n):
                c=row[v]
                if c:
                    for u in range(n):
                        if not (m>>u)&1 and A[v][u]:
                            dp[m|(1<<u)][u]+=c
        return sum(dp[full-1][v] for v in range(n))

    def tourn(gid):
        ci,p=divmod(gid,PAT)
        M=rep9_matrix(reps9[ci])
        for u in range(9):
            if (p>>u)&1: M[u,9]=1
            else: M[9,u]=1
        return M

    def canon10(M):
        rows=tuple(int(sum((1<<v) if M[uu,v] else 0 for v in range(10))) for uu in range(10))
        return canon(rows,10)

    log("scanning ext buckets...")
    found=[]; nb=0
    with open(f"{SCRATCH}/n10_ext.sorted") as f:
        cur=None; mem=[]
        for line in f:
            k,g=line.split()
            if k!=cur:
                if len(mem)>1: nb+=1; found.append(list(mem))
                cur=k; mem=[int(g)]
            else:
                mem.append(int(g))
        if len(mem)>1: nb+=1; found.append(list(mem))
    log(f"ext multi-member buckets: {nb}")
    groups=[]
    for mem in found:
        forms={}
        for gid in mem:
            cf=canon10(tourn(gid))
            if cf not in forms: forms[cf]=gid
        if len(forms)>1: groups.append(forms)
    log(f"ext panel-tied groups with >=2 distinct classes (pre-H): {len(groups)}")
    print()
    print("="*90)
    wild=0
    for forms in groups:
        gids=list(forms.values())
        mats=[tourn(g) for g in gids]
        Hs=[H_of([[int(x) for x in r] for r in M]) for M in mats]
        byH={}
        for g,M,h in zip(gids,mats,Hs): byH.setdefault(h,[]).append((g,M))
        for h,members in byH.items():
            if len(members)<2: continue
            tags=[]
            for g,M in members:
                sc=M.sum(axis=1)
                t=("SINK" if (sc==0).any() else "")+("+SRC" if (sc==9).any() else "")
                tags.append(t.strip("+") or "CORE")
            iswild=any(t=="CORE" for t in tags)
            chir=""
            if len(members)==2:
                M1,M2=members[0][1],members[1][1]
                chir=" CHIRAL(T,T^op)" if canon10(M1.T.copy())==canon10(M2) else " NON-CHIRAL"
            print(f"INVISIBLE group (H={h}): tags={tags}"
                  f"{' *** WILD ***'+chir if iswild else ''} gids={[g for g,_ in members]}",
                  flush=True)
            if iswild:
                wild+=1
                ps=[exact_panel([[int(x) for x in row] for row in M],10) for _,M in members]
                print(f"  exact re-verify: panels identical = {all(p==ps[0] for p in ps)}",flush=True)
    print(f"WILD groups at n=10: {wild}")
    log("DONE")

if __name__=="__main__":
    main()
