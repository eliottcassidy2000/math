#!/usr/bin/env python3
"""
lrc14_pairwise_overlap_klein_S293.py
====================================
klein-2026-07-13-S293 (owner: prove the pairwise coprime-overlap bound <= 1/49).

THM-739 (PROVED). For coprime speeds c,c', the two 1/14-bad sets overlap with the EXACT closed form
    |bad_c cap bad_c'| = 1/49 + (1/cc')*[ B2({(c'-c)/14}) - B2({(c'+c)/14}) ],   B2(x)=x^2-x+1/6.
Derivation: |bad_c cap bad_c'| = sum_{nc+mc'=0} 1B^(n) 1B^(m); coprime => n=c'k,m=-ck; the k-sum is a
cos/k^2 Bernoulli series = 2*pi^2*B2. Since B2 in [-1/12,1/6], the bracket in [-1/4,1/4], so
    1/49 - 1/(4cc') <= |bad_c cap bad_c'| <= 1/49 + 1/(4cc')   ->  1/49 as cc'->inf.
This verifies the closed form (direct grid overlap vs formula) and the +-1/(4cc') envelope.
"""
import numpy as np
from math import gcd
NG=1<<22; THR=1.0/14.0; t=np.arange(NG)/NG
def badmask(c):
    fr=(c*t)%1.0; d=np.minimum(fr,1.0-fr); return d<THR
def B2(x):
    x=x%1.0; return x*x-x+1.0/6.0
def direct(c,cp): return (badmask(c)&badmask(cp)).sum()/NG
def formula(c,cp): return 1.0/49.0 + (1.0/(c*cp))*(B2((cp-c)/14.0)-B2((cp+c)/14.0))
print("THM-739: |bad_c cap bad_c'| = 1/49 + (1/cc')[B2({(c'-c)/14})-B2({(c'+c)/14})]  (coprime)")
print("%10s %12s %12s %10s %s"%("(c,c')","direct","formula","|diff|","<=1/49+1/4cc'"))
maxerr=0
for c,cp in [(3,5),(5,8),(11,13),(15,22),(90,101),(23,45),(100,171),(2,3),(6,7),(13,15),(17,45),(50,99)]:
    if gcd(c,cp)!=1: continue
    dV=direct(c,cp); fV=formula(c,cp); e=abs(dV-fV); maxerr=max(maxerr,e)
    print("  (%2d,%3d)  %.8f   %.8f  %.2e   %s"%(c,cp,dV,fV,e,dV<=1.0/49.0+1.0/(4*c*cp)+1e-6))
print("max |direct-formula| = %.2e (NG=%d) ; 1/49=%.8f"%(maxerr,NG,1.0/49.0))
print("done.")
