#!/usr/bin/env python3
"""mixed_killer_gap_kps_S128c59.py -- kind-pasteur S128 cont.59.
THE MIXED CASE, which neither horn covers as stated: k1 SMALL (< 874) and k2 LARGE.
The crude union bound over all of I fails for small k1 because its bad intervals are wide
relative to |I| = 1/1092.  FIX: do not bound k1 -- REMOVE its bad set exactly.  I \ Bad_k1
is a finite union of subintervals; take the LARGEST, call its length ell1, and then k2 only
needs  2/(7 k2) < ell1*(1 - 1/7)  i.e.  k2 > 2/(6*ell1) = 1/(3*ell1).
Measure ell1 over every small k1 and every core, and read off the k2 threshold.
PRINT DATA ONLY."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
rho=F(1,2184); LO=F(1,13)-rho; HI=F(1,13)+rho; ELL=HI-LO
def bad_free_subintervals(k):
    """I minus {t : ||k t|| < 1/14}, as a list of (a,b) with a<b"""
    # bad centres j/k with |k t - j| < 1/14  ->  t in (j/k - 1/(14k), j/k + 1/(14k))
    jlo=int(LO*k)-1; jhi=int(HI*k)+1
    bad=[]
    for j in range(jlo,jhi+1):
        a=F(j,k)-F(1,14*k); b=F(j,k)+F(1,14*k)
        if b<=LO or a>=HI: continue
        bad.append((max(a,LO),min(b,HI)))
    bad.sort()
    free=[]; cur=LO
    for a,b in bad:
        if a>cur: free.append((cur,a))
        cur=max(cur,b)
    if cur<HI: free.append((cur,HI))
    return free
print("### largest bad-free subinterval of I after removing ONE small killer ###")
print("  |I| = ell = 1/1092 = %.8f"%float(ELL))
print("   k1     #free pieces  largest piece      k2 threshold 1/(3*ell1)")
worst=None
for k1 in [144,157,169,182,200,250,300,400,437,500,700,873]:
    fr=bad_free_subintervals(k1)
    if not fr:
        print("  %-6d %-13d %-18s %s"%(k1,0,"NONE","-- k1 covers all of I"))
        continue
    ell1=max(b-a for a,b in fr)
    thr=F(1,3)/ell1
    if worst is None or ell1<worst[0]: worst=(ell1,k1)
    print("  %-6d %-13d %-18.8f %.1f"%(k1,len(fr),float(ell1),float(thr)))
print("  smallest largest-piece over the scanned k1: %.8f at k1=%d"%(float(worst[0]),worst[1]))
print()
print("### exhaustive over ALL small k1: is the largest free piece always positive? ###")
mn=None; zero=[]
for k1 in range(144,874):
    fr=bad_free_subintervals(k1)
    if not fr: zero.append(k1); continue
    e=max(b-a for a,b in fr)
    if mn is None or e<mn[0]: mn=(e,k1)
print("  k1 values whose bad set swallows ALL of I: %d %s"%(len(zero),zero[:10]))
print("  minimum largest-free-piece over k1 in [144,874): %.8f at k1=%d  -> k2 threshold %.1f"%(
    float(mn[0]),mn[1],float(F(1,3)/mn[0])))
print()
print("### direct check: mixed families, k1 small and k2 large, find t in I ###")
import random
random.seed(591)
tested=0; ok=0; fails=[]
for _ in range(400):
    drop=random.randint(1,12); P=[x for x in range(1,13) if x!=drop]; M=max(P)
    k1=random.randint(13*M+1,873); k2=random.randint(874,200000)
    V=P+[k1,k2]
    if len(set(V))!=13: continue
    if not all(any(v%q==0 for v in V) for q in range(2,15)): continue
    tested+=1
    fr=bad_free_subintervals(k1)
    got=None
    for a,b in fr:
        N=600
        for j in range(N+1):
            t=a+(b-a)*F(j,N)
            if min(nd(v*t) for v in V)>=F(1,14): got=t; break
        if got: break
    if got: ok+=1
    else: fails.append((drop,k1,k2))
print("  covering mixed families tested: %d ; certified by a t in I: %d"%(tested,ok))
if fails:
    print("  not found on the grid (first 5):",fails[:5])
print("DONE")
