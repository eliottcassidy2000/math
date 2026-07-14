#!/usr/bin/env python3
"""
lrc14_band_certificate_klein_S312.py
====================================
klein-2026-07-14-S312 (owner: prove the crude M>=0.14 bound for the band residual).

RESULT (honest): the crude ANALYTIC bound does NOT exist. Two natural forms both fail:
  1. Bonferroni (B1..B7 at delta=1/14): all truncations NEGATIVE (oscillates) despite G_true~0.12.
     => no low-order inclusion-exclusion certificate.
  2. Absolute relation-lattice: rel_abs = int prod b(ct) dt - (6/7)^13 ~ 27000 (b has a log singularity,
     prod blows up at rationals). The relation series converges only CONDITIONALLY (signed), not absolutely.
     => main - rel_abs << 0. Same "signed not absolute" cancellation wall as the large-diameter route.

CORRECTION (MISTAKE-143 / THM-762): the run below found q in [15,25] on 120/120 RANDOMLY SAMPLED
rows.  It did not prove a uniform band theorem or ``loose <=> q<=25``.  Exact S312-predicate residuals
exist whose first rational witnesses have q=26 and q=27, and exact replay of the separate historical
S105 bank finds 91/8260 rows with q_min>25.  The per-family certificates printed below remain valid for
the sampled rows; only their universal interpretation was false.  See the exact refutation artifact
``lrc14_q25_uniformity_refutation_codex_S3.py``.
"""
import numpy as np, random
from math import gcd,pi,comb
NG=1<<19; t=np.arange(NG)/NG; delta=1.0/14
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
def gmeas(W):
    g=np.ones(NG,bool)
    for w in W:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); g&=(d>=delta)
    gi=g.astype(np.int8); r=int(np.sum(np.diff(np.concatenate([gi,gi[:1]]))==1)); return g.mean(),r
def vstar(P):
    m,r=gmeas(P); return (r/(pi*m)) if m>0 else 1e9
def is_band_residual(S):
    S=sorted(set(S))
    if len([x for x in S if x>14])<4: return False
    if max(S)>vstar([x for x in S if x!=max(S)]): return False
    if max(S)-min(S)>339: return False
    for w in [x for x in S if x>14]:
        if w>vstar([x for x in S if x!=w]): return False
    return True
def bonferroni(S,MAX=7):
    N=np.zeros(NG,dtype=np.int32)
    for c in S:
        fr=(c*t)%1.0; d=np.minimum(fr,1.0-fr); N+=(d<delta)
    bc=np.bincount(N)/NG; K=len(bc)
    Sk=[sum(bc[n]*comb(n,k) for n in range(K)) for k in range(MAX+1)]
    return [sum((-1)**k*Sk[k] for k in range(m+1)) for m in range(MAX+1)],bc[0]
def rational_witness(S,Qmax=200):
    for q in range(15,Qmax+1):
        lo=q/14.0; hi=q-lo
        for a in range(1,q):
            if all(lo<=(c*a)%q<=hi for c in S): return q,a
    return None,None
random.seed(31); found=0; qs=[]; ex_bonf=None
for _ in range(60000):
    C=sorted(random.sample(range(1,random.choice([80,250,450])+1),13))
    if not iscov(C) or not is_band_residual(C): continue
    if ex_bonf is None: tr,G=bonferroni(C); ex_bonf=(C,G,tr)
    q,a=rational_witness(C); 
    if q: qs.append(q)
    found+=1
    if found>=120: break
C,G,tr=ex_bonf
print("Bonferroni FAILS: Gtrue=%.4f but B1=%.2f B3=%.3f B5=%.3f B7=%.3f (oscillates negative)"%(G,tr[1],tr[3],tr[5],tr[7]))
print("Rational witness sample: %d/%d sampled band-residual families, q in [%d,%d] median~%d (not uniform; MISTAKE-143)"
      %(len(qs),found,min(qs),max(qs),sorted(qs)[len(qs)//2]))
