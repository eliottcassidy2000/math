#!/usr/bin/env python3
"""initial_segment_core_kps_S128c58.py -- kind-pasteur S128 cont.58.
THE RESIDUAL STRATUM: core = {1,...,11} (the initial segment / AP extremizer), two killers
supplying 12-, 13-, 14-divisibility.  These are the only covering wide clusters that the
small-modulus avoidance criterion (q <= 28) misses.  Find their exact witnesses / M.
PRINT DATA ONLY."""
import sys
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def best_witness(V,qhi=3000):
    """largest min||v a/q|| over a/q with q<=qhi; early-exit once >=1/14"""
    best=(F(0),None)
    for q in range(2,qhi+1):
        for a in range(1,q//2+1):
            if gcd(a,q)!=1: continue
            t=F(a,q); mn=min(nd(v*t) for v in V)
            if mn>best[0]:
                best=(mn,(q,a))
                if mn>=F(1,14): return best
    return best
CORE=list(range(1,12))
print("### the 6 misses: core {1..11} ###")
for K in [[312,364],[312,392],[168,338],[252,338],[338,420],[338,504]]:
    V=CORE+K
    cov=all(any(v%q==0 for v in V) for q in range(2,15))
    mn,qa=best_witness(V)
    print("  K=%-12s covering(2..14):%-5s  best t=%-10s min||vt||=%-9s (%.5f)  >=1/14: %s"%(
        str(K),cov,("%d/%d"%(qa[1],qa[0])) if qa else "-",mn,float(mn),mn>=F(1,14)))
print()
print("### full sweep of the drop=12 stratum: core {1..11}, killers > 143 ###")
print("  enumerate killer pairs supplying 12|,13|,14| ; report any with NO witness q<=3000")
rows=0; bad=[]; worst=(F(1),None)
seen=set()
for k1 in range(144,900):
    for k2 in range(k1+1,900):
        V=CORE+[k1,k2]
        if len(set(V))!=13: continue
        if not all(any(v%q==0 for v in V) for q in range(2,15)): continue
        rows+=1
        # quick pass: small-modulus avoidance up to q=40 with threshold ceil(q/14)
        found=None
        for q in range(15,41):
            thr=-(-q//14)   # ceil(q/14)
            for a in range(1,q):
                if all(min((v*a)%q, q-((v*a)%q))>=thr for v in V): found=(q,a); break
            if found: break
        if not found:
            mn,qa=best_witness(V,qhi=1200)
            if mn<F(1,14): bad.append((k1,k2,mn,qa))
            if mn<worst[0]: worst=(mn,(k1,k2))
print("  covering families with core {1..11}: %d"%rows)
print("  with NO small-modulus (q<=40) certificate AND no witness q<=1200 at level 1/14: %d"%len(bad))
if bad:
    for k1,k2,mn,qa in bad[:10]:
        print("    K=[%d,%d] best=%s (%.5f) at %s"%(k1,k2,mn,float(mn),qa))
print("  smallest best-witness value seen among q<=40-uncertified: %s at %s"%(worst[0],worst[1]))
print("DONE")
