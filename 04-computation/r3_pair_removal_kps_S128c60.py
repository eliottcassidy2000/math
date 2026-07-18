#!/usr/bin/env python3
"""r3_pair_removal_kps_S128c60.py -- kind-pasteur S128 cont.60 (background).
Worst largest-surviving-interval after removing TWO small killers exactly, over all 66
ten-speed cores.  The r=2 experience says the worst case sits at the TOP of the killer
range (dense bad intervals), so scan all pairs in the top window exactly, then sample the
rest to confirm.  PRINT DATA ONLY."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def safe_set(P):
    bps={F(0),F(1)}
    for p in P:
        for j in range(p+1):
            for s in (F(1,14*p),-F(1,14*p)):
                v=F(j,p)+s
                if 0<=v<=1: bps.add(v)
    B=sorted(bps); out=[]
    for i in range(len(B)-1):
        a,b=B[i],B[i+1]
        if b<=a: continue
        if all(nd(p*((a+b)/2))>=F(1,14) for p in P):
            if out and out[-1][1]==a: out[-1]=(out[-1][0],b)
            else: out.append((a,b))
    return out
def remove_bad_f(iv,k):
    """float version for speed; iv is list of (a,b) floats"""
    out=[]; w=1.0/(14.0*k)
    for a,b in iv:
        jlo=int(a*k)-1; jhi=int(b*k)+1; cur=a
        for j in range(jlo,jhi+1):
            x=j/k-w; y=j/k+w
            if y<=a or x>=b: continue
            x=max(x,a); y=min(y,b)
            if x>cur: out.append((cur,x))
            cur=max(cur,y)
        if cur<b: out.append((cur,b))
    return out
CORES=[sorted(c) for c in itertools.combinations(range(1,13),10)]
KB=900; TOP=110
print("### worst L after removing TWO small killers, all 66 ten-speed cores ###")
print("  killer range (13*maxP, %d) ; exact top-%d window all-pairs + random sample"%(KB,TOP))
random.seed(600)
worst=None; percore=[]
for P in CORES:
    ivF=[(float(a),float(b)) for a,b in safe_set(P)]
    M=max(P); lo=13*M+1
    ks=list(range(lo,KB))
    cand=ks[-TOP:]
    extra=random.sample(ks,min(90,len(ks)))
    pool=sorted(set(cand+extra))
    w=None
    for i in range(len(pool)):
        r1=remove_bad_f(ivF,pool[i])
        if not r1: w=(0.0,pool[i],None); break
        for j in range(i+1,len(pool)):
            rem=remove_bad_f(r1,pool[j])
            L=max((b-a for a,b in rem),default=0.0)
            if w is None or L<w[0]: w=(L,pool[i],pool[j])
    percore.append((w[0],tuple(P),w[1],w[2]))
    if worst is None or w[0]<worst[0]: worst=(w[0],tuple(P),w[1],w[2])
    print("  core %-40s worst L=%.8f at (%s,%s) -> thr %.1f"%(
        str(list(P)),w[0],w[1],w[2],(1.0/(3*w[0]) if w[0]>0 else -1)))
print()
L,Pw,a,b=worst
print("  GLOBAL WORST L = %.8f at core %s killers (%s,%s)"%(L,list(Pw),a,b))
print("  => r=3 measure-horn threshold for the third killer: 1/(3L) = %.1f"%(1.0/(3*L) if L>0 else -1))
percore.sort()
print("  five worst cores:", [(list(p),round(l,8)) for l,p,_,_ in percore[:5]])
print("DONE")
