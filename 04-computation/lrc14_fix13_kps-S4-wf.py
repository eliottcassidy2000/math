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
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b
def gaps(W,Vstar):
    R=F(1,2*Vstar); h=F(1,14); forb=[]
    for w in set(W):
        j=0
        while True:
            c=F(j,w); lo=c-h/F(w)
            if lo>R: break
            forb.append((max(F(0),lo),min(R,c+h/F(w)))); j+=1
    forb.sort(); merged=[]
    for a,b in forb:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    prev=F(0); g=[]
    for a,b in merged:
        if a>prev: g.append((prev,a))
        prev=max(prev,b)
    if prev<R: g.append((prev,R))
    return [(lo,hi) for lo,hi in g if hi>lo and hi>0]
# exactly 13 elements: 10 nonmult7 (incl 19=V*) + 3 mult7 (21,28,273)
# nonmult7 must cover q in {2,3,4,5,6,8,9,10,11,12,13} together with nothing from m7
# (m7 gives 7,14 and also 3(21),4(28),13(273),6(?). but require covering overall.)
nm7=[1,2,5,9,11,12,13,16,18,19]  # 10 elements, max 19
m7=[21,28,273]
S=sorted(nm7+m7)
print("S=",S,"|S|=",len(S))
print("cov=",is_cov(S),"prim=",is_prim(S))
Vmin=min(S);Vmax=max(S);k=sum(1 for v in S if v>13)
print("S3?",k>=2 and Vmax>=13*Vmin,"V*=",max(v for v in S if v%7!=0))
W=[v//7 for v in S if v%7==0]
print("W=",sorted(W),"window gaps=",gaps(W,max(v for v in S if v%7!=0)))
print("Mval=",Mval(S),"~",float(Mval(S)),">=1/14?",Mval(S)>=F(1,14))
