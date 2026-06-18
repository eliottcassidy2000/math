"""
Fast targeted adversarial probe for LRC(14) case S3.
Goal: hunt for a covering S3 13-set with exact M < 1/14 (would REFUTE LRC14),
and for any case where 'Lemma A fires' but the witness is not actually safe.

Keep cluster speeds modest (<=400) so exact Mval is fast; sweep MANY offset
patterns and small-part choices, including deliberately near-tight S3 sets
(Vmax just over 13*Vmin, tight clusters).
"""
import sys
from fractions import Fraction as F
from math import gcd
import random, itertools

def flush(*a): print(*a); sys.stdout.flush()

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r

def safe_components(A, h=F(1,14)):
    iv=[]
    for u in A:
        for j in range(0,u):
            c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
            if a<b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe

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
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b

def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def gcd_list(xs):
    g=0
    for x in xs: g=gcd(g,x)
    return g
H=F(1,14)

def lemmaA_window(K,Vmin,Vmax): return F(14*K+1,14*Vmin), F(14*K+13,14*Vmax)
def lemmaA_fires(P,L):
    Vmin=min(L); Vmax=max(L); Ps=safe_components(P,H); s=Vmax-Vmin
    if s==0: Kmax=2 if (13*Vmin-Vmax)>0 else -1
    else:
        if 13*Vmin-Vmax<=0: return None
        Kmax=(13*Vmin-Vmax-1)//(14*s)
    for K in range(0,Kmax+1):
        lo,hi=lemmaA_window(K,Vmin,Vmax)
        if not(lo<hi): continue
        for (a,b) in Ps:
            ov=(max(lo,a),min(hi,b))
            if ov[0]<ov[1]: return (ov[0]+ov[1])/2
    return None

def run(seed=0, trials=30000, maxV=400):
    rnd=random.Random(seed)
    minM=None; minS=None; below=[]; fire_fail=[]; ncov=0
    for _ in range(trials):
        psz=rnd.randint(1,9)
        P=sorted(rnd.sample(range(1,14), psz))
        csz=rnd.randint(2, max(2,13-len(P)))
        V0=rnd.randint(20,maxV)
        # tight cluster: spread small to put many in same residue band
        spread=rnd.choice([1,2,3,5,8,12,20,30,45,60])
        L=set()
        while len(L)<csz*4:
            L.add(V0+rnd.randint(0,spread))
            if rnd.random()<0.05: break
        L=sorted(L)[:csz]
        if len(L)<2: continue
        S=sorted(set(P)|set(L))
        if len(S)<13: continue
        S=S[:13]
        if len(set(S))!=13: continue
        if gcd_list(S)!=1: continue
        if not is_covering(S): continue
        Vmin=min(S); Vmax=max(S)
        if sum(1 for v in S if v>13)<2: continue
        if Vmax<13*Vmin: continue
        ncov+=1
        m=Mval(S)
        if minM is None or m<minM: minM=m; minS=S
        if m<H: below.append((S,m))
        Pp=[v for v in S if v<=13]; Ll=[v for v in S if v>13]
        w=lemmaA_fires(Pp,Ll)
        if w is not None:
            if not all(nrm(v*w)>=H for v in S) or m<H:
                fire_fail.append((S,float(w),m))
    return ncov,minM,minS,below,fire_fail

if __name__=="__main__":
    for seed in range(4):
        ncov,minM,minS,below,ff=run(seed=seed, trials=25000, maxV=400)
        flush(f"seed={seed}: covering S3 tested={ncov}  minM={minM} ({float(minM) if minM else None})  below1/14={len(below)}  fireFail={len(ff)}")
        if below: flush("   BELOW:", below[:3])
        if ff: flush("   FIREFAIL:", ff[:3])
        flush(f"   witness minS={minS}")
    flush("1/14 =", float(H))
