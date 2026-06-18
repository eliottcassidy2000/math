# LRC(14) — focused diagnostic on KNOWN hard S3 sets for the private-q / clearing-crossing angle.
# kps-S4-wf.  Does the parked (mult-of-14) runner drive M, or do small pairs govern?

from fractions import Fraction as F
from math import gcd
from itertools import combinations
def pr(*a): print(*a, flush=True)
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
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
    b=F(0); bt=None
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v; bt=t
    return b,bt
def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def is_prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def binding_set(S,t,M): return [v for v in S if nrm(v*t)==M]
def owners(S): return {q:[v for v in S if v%q==0] for q in range(2,15)}
def private_qs(S):
    ow=owners(S); pq={}
    for q in range(2,15):
        if len(ow[q])==1: pq[q]=ow[q][0]
    return pq  # q -> the unique runner covering q
def crossing_reps(S,t,M):
    B=binding_set(S,t,M); reps=[]
    for va,vb in combinations(sorted(set(B)),2):
        for D in (va+vb,abs(va-vb)):
            if D and (t*D).denominator==1:
                reps.append((va,vb,D,int(t*D)))
    return B,reps

HARD = {
 "S* (MISTAKE-076)": [1,2,3,5,7,8,9,10,11,12,13,38,42],
 "G2 (HYP-2581)": [1,2,3,5,7,8,9,10,11,12,13,27,28],
 "champion 84": [1,2,3,4,5,6,7,8,9,10,11,13,84],
 "champion 182 (k=1!)": [1,2,3,4,5,6,7,8,9,10,11,12,182],
 "wide cluster": [2,6,8,12,20,78,104,130,156,182,208,234,260],
}
for name,S in HARD.items():
    S=sorted(S)
    if len(S)!=13:
        pr(f"--- {name}: NOT 13 elements ({len(S)}) skipping"); continue
    M,t=Mval(S)
    B,reps=crossing_reps(S,t,M)
    pq=private_qs(S)
    ow=owners(S)
    kbig=sum(1 for v in S if v>13)
    Vmin,Vmax=min(S),max(S)
    case = "S1" if kbig<=1 else ("S2" if Vmax<13*Vmin else "S3")
    mult14 = ow[14]
    mult7  = ow[7]
    # parked = unique mult of 14
    parked = ow[14][0] if len(ow[14])==1 else None
    parked_binding = (parked in B) if parked is not None else None
    pr(f"=== {name}: S={S}")
    pr(f"    cov={is_cov(S)} prim={is_prim(S)} k(big)={kbig} case={case} Vmin={Vmin} Vmax={Vmax}")
    pr(f"    M={M} (M*14={M*14}, >=1/14? {M>=F(1,14)}) tau*={t}")
    pr(f"    binding runners at tau*: {B}")
    pr(f"    crossing reps (va,vb,D,j): {reps}")
    pr(f"    j>=D/14 for all reps: {all(14*j>=D for *_,D,j in [(a,b,D,j) for a,b,D,j in reps])}")
    pr(f"    mult-of-14: {mult14}  mult-of-7: {mult7}")
    pr(f"    private q-owners (q->runner): {pq}")
    pr(f"    parked(unique mult14)={parked}  parked in binding pair? {parked_binding}")
    # does the parked runner's q-debt force the crossing?  Check: is parked in ANY crossing rep?
    parked_in_rep = parked is not None and any(parked in (a,b) for a,b,D,j in reps)
    pr(f"    parked appears in some crossing rep: {parked_in_rep}")
    # remove parked -> does M drop?
    if parked is not None:
        Sm=[v for v in S if v!=parked]
        Mm,_=Mval(Sm)
        pr(f"    M(S minus parked {parked}) = {Mm} (>=1/14? {Mm>=F(1,14)})  -> parked essential to LRC? {Mm<F(1,14)}")
    pr("")
