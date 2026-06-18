"""
Probe the NO-WITNESS cases found by the bulk hunt. For each:
 - confirm it's REALLY ALL-MULT7-LARGE (every mult of 7 > V*=max non-mult7), covering, primitive, S3
 - recompute V*, the mult-of-7 subsystem W={v//7}, and WHY safe_u returns None
 - compute Mval to see if M>=1/14 anyway (closure conclusion) vs the WITNESS mechanism failing
This distinguishes: (a) closure FALSE (M<1/14 found) vs (b) my safe_u/window scope wrong.
kind-pasteur S4-wf.
"""
from fractions import Fraction as F
from math import gcd
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def is_prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
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
    b=F(0); arg=None
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v; arg=t
    return b,arg
def safe_u(W,Vstar):
    R=F(1,2*Vstar); h=F(1,14); forb=[]
    for w in set(W):
        j=0
        while True:
            c=F(j,w); lo=c-h/F(w)
            if lo>R: break
            a=max(F(0),lo); b=min(R,c+h/F(w))
            if a<=b: forb.append((a,b))
            j+=1
    forb.sort(); merged=[]
    for a,b in forb:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    prev=F(0); gaps=[]
    for a,b in merged:
        if a>prev: gaps.append((prev,a))
        prev=max(prev,b)
    if prev<R: gaps.append((prev,R))
    for lo,hi in gaps:
        if hi<=0: continue
        u=(lo+hi)/2
        if u<=0: u=hi/2
        if u>0 and all(nrm(F(w)*u)>=h for w in W): return u,gaps,R
    for lo,hi in gaps:
        for u in (lo,hi):
            if u>0 and u<=R and all(nrm(F(w)*u)>=h for w in W): return u,gaps,R
    return None,gaps,R

nowit=[
[10, 24, 38, 44, 58, 65, 79, 136, 137, 147, 693, 1848, 1988],
[12, 19, 20, 23, 37, 52, 54, 94, 110, 112, 686, 1407, 1526],
[1, 5, 34, 39, 54, 74, 113, 120, 137, 140, 252, 1386, 1904],
[4, 8, 12, 37, 40, 46, 51, 55, 60, 63, 1547, 1736, 1841],
[19, 23, 24, 26, 29, 33, 36, 37, 39, 40, 42, 308, 1155],
[1, 6, 8, 10, 15, 16, 18, 22, 24, 26, 28, 378, 581],
[10, 17, 19, 25, 32, 45, 102, 109, 182, 588, 1344, 1456, 1463],  # near worst region
]
for S in nowit:
    S=sorted(S)
    m7=[v for v in S if v%7==0]; nm7=[v for v in S if v%7!=0]
    Vstar=max(nm7); allbig=all(x>Vstar for x in m7)
    W=[v//7 for v in m7]; wmin=min(W) if W else None
    R=F(1,2*Vstar)
    u,gaps,_=safe_u(W,Vstar)
    M,arg=Mval(S)
    print("S=",S)
    print("   cov=",is_cov(S)," prim=",is_prim(S)," m7=",m7," V*=",Vstar," ALLBIG=",allbig)
    print("   W=",sorted(W)," wmin=",wmin," 7*wmin=",7*wmin if wmin else None," R=1/(2V*)=",R,"~",float(R))
    print("   first-tooth half-width 1/(14 wmin)=",F(1,14*wmin),"~",float(F(1,14*wmin))," vs R=",float(R),
          "  (does first tooth of wmin swallow window? ", F(1,14*wmin)>=R,")")
    print("   safe_u in window:", u, "  #gaps=",len(gaps))
    print("   Mval=",M,"~",float(M)," >=1/14?",M>=F(1,14)," arg tau=",arg)
    print()
