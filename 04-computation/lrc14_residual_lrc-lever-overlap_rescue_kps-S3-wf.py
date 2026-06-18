#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_rescue_kps-S3-wf

HONEST EXAMINATION of the ruler-w0 FAILURES, and rescue attempts.

Two failures found:
  (A) AP-cluster S3 set: S=[4,7,8,9,11,12,28,32,36,40,44,48,52] -- 0 P-safe periods for ruler w0.
  (B) Synthetic extreme P=range(1,14), Delta={0,1}: 0 good periods at all w0.
      (But is (B) a REAL covering S3 set?  P=1..13 is 13 speeds already = full set, leaving NO
       room for a cluster!  So (B) is NOT a valid 13-set -- it's an artifact of the synthetic
       probe.  CHECK THIS.)

We test:
  (1) Is each failure a genuine covering primitive 13-set in case S3?
  (2) Does C(S) (truth: W(S\{v})>1/(7v) for SOME v) still hold?  (it must, by the verified S2/full)
  (3) RESCUE: does a DIFFERENT cluster member w' (not min) as ruler give a P-safe period?
      i.e. for some w' in L', do the periods of w' yield a >thresh window safe for all of A?
  (4) RESCUE 2: does removing a DIFFERENT v (not Vmax) make C fire trivially?
  (5) The exact M(S) and which v gives the widest W(S\{v}).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import sys

def flush(*a):
    print(*a); sys.stdout.flush()
C=F(1,14)

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def teeth_intervals(u,h=C):
    iv=[]
    for j in range(0,u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merged_intervals(iv):
    iv=sorted(iv); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    return merged
def danger_merged(A,h=C):
    iv=[]
    for u in A: iv+=teeth_intervals(u,h)
    return merged_intervals(iv)
def safe_from_danger(merged):
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def safe_components(A,h=C): return safe_from_danger(danger_merged(A,h))
def Wwidth(A):
    sc=safe_components(A)
    if not sc: return F(0)
    ws=[b-a for a,b in sc]
    if sc[0][0]==0 and sc[-1][1]==1 and len(sc)>1: ws.append(sc[0][1]+(1-sc[-1][0]))
    return max(ws)
def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def gcd_all(S): return reduce(gcd,S)

def cand(S):
    S=sorted(set(S)); Cs=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): Cs.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): Cs.add(F(k,d)); k+=1
    Cs.add(F(1,2)); return Cs
def Mval(S):
    b=F(0); at=None
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v; at=t
    return b,at

def ruler_good(P, Lp, Vmax, ruler):
    """Count periods of `ruler` giving a >thresh window safe for all of P u Lp."""
    thresh=F(1,7*Vmax)
    GA = safe_components(list(P)+list(Lp))
    cnt=0
    for j in range(0,ruler):
        lo=F(14*j+1,14*ruler); hi=F(14*j+13,14*ruler)
        if hi>1: hi=F(1)
        if lo>=hi: continue
        # any >thresh GA piece inside [lo,hi]?
        for a,b in GA:
            x=max(a,lo); y=min(b,hi)
            if y-x>thresh: cnt+=1; break
    return cnt

def examine(S):
    S=sorted(set(S))
    flush(f"\n  S={S}")
    flush(f"   |S|={len(S)} covering={is_covering(S)} gcd={gcd_all(S)}")
    Vmin=min(S); Vmax=max(S); k=sum(1 for v in S if v>13)
    flush(f"   Vmin={Vmin} Vmax={Vmax} k(>13)={k}  S3? {k>=2 and Vmax>=13*Vmin}")
    if len(S)!=13:
        flush("   *** NOT a 13-set -- artifact, skip ***"); return
    M,at=Mval(S)
    flush(f"   M(S)={M}={float(M):.5f}  (>=1/14={M>=C})  attained at tau={at}")
    # which v gives widest W(S\{v}) and does C fire?
    bestv=None; bestratio=F(0)
    for v in S:
        W=Wwidth([x for x in S if x!=v]); ratio=W*7*v
        if ratio>bestratio: bestratio=ratio; bestv=v
    flush(f"   best C-firing v={bestv}, W(S\\v)*7v = {float(bestratio):.4f}  (C fires if >1)")
    # ruler rescue: try every cluster member as ruler
    P=[v for v in S if v<=13]; L=sorted(v for v in S if v>13); Lp=[v for v in L if v!=Vmax]
    flush(f"   P={P}  Lp(cluster minus Vmax)={Lp}")
    for ruler in Lp:
        g=ruler_good(P,Lp,Vmax,ruler)
        flush(f"     ruler={ruler}: #good periods = {g}")

if __name__=="__main__":
    flush("="*70)
    flush("RESCUE / honest examination of ruler-w0 failures")
    flush("="*70)
    examine([4,7,8,9,11,12,28,32,36,40,44,48,52])   # AP-cluster failure
    examine([1,2,3,4,5,6,7,8,9,10,11,12,13])         # the synthetic 'failure' (NOT S3-valid)
    # a real S3 set with P=most of 1..13 and tiny cluster
    examine([1,2,3,4,5,6,7,8,9,10,11,420,421])
    examine([1,2,3,4,5,6,7,8,9,10,11,840,841])
