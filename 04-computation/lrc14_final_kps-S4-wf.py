"""
FINAL consolidation. kind-pasteur S4-wf.
(1) Canonical 13-set counterexample to the WITNESS-EXISTENCE sub-claim.
(2) Large targeted hunt: does ANY ALL-MULT7-LARGE covering primitive S3 set have
    M < 1/14 (which would REFUTE the closure CONCLUSION, not just the mechanism)?
    Focus on no-witness sets (where the claim's mechanism fails) and tight builds.
"""
from fractions import Fraction as F
from math import gcd
import random
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
def has_witness(S):
    m7=[v for v in S if v%7==0]; nm7=[v for v in S if v%7!=0]
    Vstar=max(nm7); W=[v//7 for v in m7]; R=F(1,2*Vstar); h=F(1,14); forb=[]
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
    return any(hi>lo and hi>0 for lo,hi in g)

print("(1) CANONICAL counterexample to witness-existence sub-claim:")
S0=[1, 6, 8, 10, 15, 16, 18, 22, 24, 26, 28, 378, 581]
print("    S=",S0," cov=",is_cov(S0)," prim=",is_prim(S0),
      " m7=",[v for v in S0 if v%7==0]," V*=",max(v for v in S0 if v%7!=0))
print("    ALL-MULT7-LARGE:",all(v>max(x for x in S0 if x%7!=0) for v in S0 if v%7==0),
      " has tau=k/7+u/7 witness:",has_witness(S0)," Mval=",Mval(S0),"~",float(Mval(S0)))

print("\n(2) Hunt for ANY ALL-MULT7-LARGE covering prim S3 set with M<1/14 ...")
random.seed(99)
worst=None; nbad=0; ntested=0; nnowit=0; worst_nowit=None
for _ in range(25000):
    Vstar=random.randint(13,160)
    if Vstar%7==0: continue
    pool=[v for v in range(1,Vstar+1) if v%7!=0]
    if len(pool)<8: continue
    core=set([Vstar])
    while len(core)<random.randint(8,11):
        core.add(random.choice(pool))
    if max(core)!=Vstar: continue
    nslots=13-len(core)
    if nslots<2: continue
    m7pool=[7*w for w in range(1,300) if 7*w>Vstar]
    m7_14=[x for x in m7pool if x%14==0]
    if not m7_14: continue
    m7=set([random.choice(m7_14)])
    while len(m7)<nslots: m7.add(random.choice(m7pool))
    if len(m7)!=nslots: continue
    S=sorted(core)+sorted(m7)
    if len(S)!=13 or not is_cov(S) or not is_prim(S): continue
    if not all(x>Vstar for x in m7): continue
    Vmin=min(S);Vmax=max(S);k=sum(1 for v in S if v>13)
    if not(k>=2 and Vmax>=13*Vmin): continue
    ntested+=1
    M=Mval(S)
    if M<F(1,14):
        nbad+=1
        if nbad<=5: print("    *** M<1/14 ***",S,"M=",M)
    if worst is None or M<worst[1]: worst=(S,M)
    if not has_witness(S):
        nnowit+=1
        if worst_nowit is None or M<worst_nowit[1]: worst_nowit=(S,M)
print("    tested=",ntested," sets with M<1/14:",nbad," no-witness sets:",nnowit)
if worst: print("    overall worst Mval:",worst[1],"~",float(worst[1])," set:",worst[0])
if worst_nowit: print("    worst Mval AMONG no-witness sets:",worst_nowit[1],"~",float(worst_nowit[1]),
                      " set:",worst_nowit[0])
print("    1/14 =",float(F(1,14)))
