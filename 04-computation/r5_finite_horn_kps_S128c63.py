#!/usr/bin/env python3
"""r5_finite_horn_kps_S128c63.py -- kind-pasteur S128 cont.63.
r=5 FINITE HORN.  The measure horn FAILS at r=5 (max R = 1.285 > 1), so unlike r<=3 the
finite horn is MANDATORY here.  It only has to cover the region where R >= 1, which the
crossing analysis puts at the very bottom of the killer range.
STEP 1: find, over all 495 eight-speed cores, the largest killer at which R >= 1.
STEP 2: run the finite horn over quintuples below that bound (plus margin), with the same
sound covering-necessary pruning (sum frac >= 1).  PRINT DATA ONLY."""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def la(r,q):
    r%=q; return min(r,q-r)
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
        if all(nd(p*((a+b)/2))>=F(1,14) for p in P): out.append((float(a),float(b)))
    mg=[]
    for a,b in out:
        if mg and abs(mg[-1][1]-a)<1e-15: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    return mg
def rm(iv,k):
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
def thresh(iv):
    if not iv: return None
    mu=sum(b-a for a,b in iv); N=len(iv); L=max(b-a for a,b in iv)
    if mu<=0 or L<=0: return None
    return min(N/(6*mu), 1.0/(3*L))
C8=[sorted(c) for c in itertools.combinations(range(1,13),8)]
print("### STEP 1: the region where R >= 1 (measure horn fails) ###")
W=22
maxk=0; cases=0
for P in C8:
    M=max(P); lo=13*M+1; iv0=safe_set(P)
    ks=list(range(lo,lo+W))
    for i in range(len(ks)):
        r1=rm(iv0,ks[i])
        if not r1: continue
        for j in range(i+1,len(ks)):
            r2=rm(r1,ks[j])
            if not r2: continue
            for h in range(j+1,len(ks)):
                r3=rm(r2,ks[h])
                if not r3: continue
                for m in range(h+1,len(ks)):
                    T=thresh(rm(r3,ks[m]))
                    if T is None: continue
                    if T/ks[m]>=1.0:
                        cases+=1
                        if ks[m]>maxk: maxk=ks[m]
print("  quadruples with R >= 1: %d ; largest k4 among them: %d"%(cases,maxk))
KB=maxk+20
print("  => finite horn bound KB = %d (largest R>=1 killer + 20 margin)"%KB)
print()
print("### STEP 2: the r=5 finite horn below KB ###")
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
tested=0; cert=0; bad=[]
for ci,P in enumerate(C8):
    M=max(P); lo=13*M+1
    if lo>=KB: continue
    bits=[i for i,(q,a) in enumerate(QS) if all(la(p*a,q)>=-(-q//14) for p in P)]
    n=len(bits)
    if n==0: continue
    FULL=(1<<n)-1
    ent=[]
    for k in range(lo,KB):
        km=0
        for jj,i in enumerate(bits):
            q,a=QS[i]
            if la(k*a,q)<-(-q//14): km|=(1<<jj)
        ent.append((bin(km).count("1")/n,k,km))
    ent.sort(key=lambda e:-e[0])
    fr=[e[0] for e in ent]; kv=[e[1] for e in ent]; km=[e[2] for e in ent]
    L=len(ent)
    for i in range(L):
        if fr[i]*5<1.0: break
        for j in range(i+1,L):
            if fr[i]+fr[j]*4<1.0: break
            u2=km[i]|km[j]
            for k in range(j+1,L):
                if fr[i]+fr[j]+fr[k]*3<1.0: break
                u3=u2|km[k]
                for l in range(k+1,L):
                    if fr[i]+fr[j]+fr[k]+fr[l]*2<1.0: break
                    u4=u3|km[l]
                    for m in range(l+1,L):
                        if fr[i]+fr[j]+fr[k]+fr[l]+fr[m]<1.0: break
                        tested+=1
                        if (u4|km[m])!=FULL: cert+=1; continue
                        K=tuple(sorted((kv[i],kv[j],kv[k],kv[l],kv[m])))
                        V=P+list(K)
                        if len(set(V))!=13: continue
                        if not all(any(v%q==0 for v in V) for q in range(2,15)): continue
                        bad.append((tuple(P),K))
    if ci%100==0: print("  ... core %d/%d tested=%d bad=%d"%(ci,len(C8),tested,len(bad)))
print()
print("  quintuples passing the necessary condition: %d"%tested)
print("  certified: %d ; UNCERTIFIED: %d"%(cert,len(bad)))
if bad:
    for P,K in bad[:10]: print("    core=%s K=%s"%(list(P),list(K)))
print("DONE")
