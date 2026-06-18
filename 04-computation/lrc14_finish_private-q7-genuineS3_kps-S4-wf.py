# LRC(14) — GENUINE S3 (k>=2 big runners) test of the clearing-crossing / private-q angle.
# kps-S4-wf.  All prints flushed so output streams even under pipe buffering.
#
# S3 def: k=#{v>13}>=2 AND Vmax>=13*Vmin.  Covering: mult of every q in 2..14.
# We test:
#   (A) Does M(S)>=1/14 always? (LRC on S3)
#   (B) Is tau* a binding-pair crossing j/D with j>=D/14? (the clearing-crossing claim)
#   (C) Is the parked (unique mult of 14) runner ITSELF binding at tau*, or does a
#       small/non-parked pair govern? (the mechanism question)
#   (D) Private q-debt: for sets where the mult-of-14 runner uniquely covers q=14
#       (and possibly q=7), is M>=1/14 driven by that runner, or independent of it?
#
# Efficient generator: pick a small core P (>=2 elements, subset of 1..13) that with one
# big block can be made covering; force >=2 large runners.

import sys
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random

def pr(*a): print(*a, flush=True)

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def Mval(S):
    b=F(0); bt=None
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v; bt=t
    return b,bt
def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def is_prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def binding_set(S,t,M): return [v for v in S if nrm(v*t)==M]
def owners(S): return {q:[v for v in S if v%q==0] for q in range(2,15)}

def crossing_reps(S,t,M):
    B=binding_set(S,t,M); reps=[]
    for va,vb in combinations(sorted(set(B)),2):
        for D in (va+vb,abs(va-vb)):
            if D and (t*D).denominator==1:
                reps.append((va,vb,D,int(t*D)))
    return B,reps

def gen_S3(num, seed, maxbig):
    """k>=2 big runners, covering, primitive, Vmax>=13 Vmin."""
    random.seed(seed)
    out=[]; tries=0
    while len(out)<num and tries<num*400:
        tries+=1
        kbig=random.choice([2,2,2,3,3,4])
        ksmall=13-kbig
        if ksmall<1 or ksmall>13: continue
        if ksmall>13: continue
        P=sorted(random.sample(range(1,14),min(ksmall,13)))
        if len(P)!=ksmall:  # ksmall may exceed 13
            continue
        L=sorted(random.sample(range(14,maxbig+1),kbig))
        S=sorted(set(P+L))
        if len(S)!=13: continue
        if not is_prim(S): continue
        if not is_cov(S): continue
        Vmin=min(S);Vmax=max(S)
        if sum(1 for v in S if v>13)<2: continue
        if Vmax<13*Vmin: continue
        out.append(S)
    return out

if __name__=="__main__":
    pr("="*72)
    pr("GENUINE S3 (k>=2 big runners, Vmax>=13 Vmin), covering primitive sets")
    pr("="*72)
    S3=gen_S3(num=3000, seed=2024, maxbig=300)
    pr(f"generated {len(S3)} genuine-S3 covering primitive sets")
    n_break=0; minM=F(1); minMset=None
    n_cross=0; n_j_ok=0; n_nocross=0
    n_uniq14=0; n_parked_binding=0
    n_remove14_drops=0   # removing the unique-mult-of-14 runner pushes M<1/14
    worst=[]
    for S in S3:
        M,t=Mval(S)
        if M<F(1,14): n_break+=1
        if M<minM: minM=M; minMset=S
        B,reps=crossing_reps(S,t,M)
        if reps:
            n_cross+=1
            if all(14*j>=D for *_,D,j in [(a,b,D,j) for a,b,D,j in reps]):
                n_j_ok+=1
        else:
            n_nocross+=1
        ow=owners(S)
        if len(ow[14])==1:
            n_uniq14+=1
            w=ow[14][0]
            if w in B: n_parked_binding+=1
            Sm=[v for v in S if v!=w]
            Mm,_=Mval(Sm)
            if Mm<F(1,14): n_remove14_drops+=1
        worst.append((M,S))
    worst.sort(key=lambda x:x[0])
    pr(f"LRC breaks (M<1/14): {n_break}")
    pr(f"min M = {minM}  (={float(minM):.5f}; 1/14={float(F(1,14)):.5f}) on {minMset}")
    pr(f"tau* IS a binding-pair crossing j/D: {n_cross}/{len(S3)} (no-cross {n_nocross})")
    pr(f"  of those, ALL crossing reps satisfy j>=D/14: {n_j_ok}/{n_cross}")
    pr(f"sets with unique mult-of-14 (parked) runner: {n_uniq14}")
    pr(f"  parked runner is BINDING at tau*: {n_parked_binding}/{n_uniq14}")
    pr(f"  removing parked runner pushes M<1/14 (parked essential): {n_remove14_drops}/{n_uniq14}")
    pr("\n10 hardest (smallest M) genuine-S3 sets:")
    for M,S in worst[:10]:
        B,reps=crossing_reps(S,*Mval(S))  # recompute t
        Mv,t=Mval(S)
        ow=owners(S)
        w14 = ow[14][0] if len(ow[14])==1 else ow[14]
        pr(f"  M={M} ({float(M):.5f}) tau*={t} binding={B} mult14={w14} S={S}")
