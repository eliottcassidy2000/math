"""Exact overlapping-admission tests, independent of a single-parent map."""
from check import Q,subsets,weight,coefficients,elementary,audit_kernel,need
import check as base
from itertools import combinations
from collections import defaultdict
from math import comb
import json

def choose(n,d):return comb(n,d) if 0<=d<=n else 0

def bound(m,k,S,a,b,A,Ba,Bb,K):
    ta=k-a;tb=k-b
    lo=sum(Q(A[a-d],Ba)*K**d*Q(choose(m-ta,d),comb(ta+d,d)) for d in range(1,a+1)) if 0<=ta<=m else Q(0)
    hi=sum(Q(A[b+d],Bb)*K**d*Q(choose(tb,d),comb(m-tb+d,d)) for d in range(1,S-b+1)) if 0<=tb<=m else Q(0)
    return 1+max(lo,hi)

def small(w,k,q):
    # Rebuild the reference quotient independently from literal auxiliary states.
    m=len(w);S=4;a=1;b=3;d=2;u=Q(2);rw=[Q(1),Q(2),Q(3),Q(1,2)]
    Z=coefficients([0,1,0,1],rw,2);labs=[i%2 for i in range(m)]
    charge=lambda J:sum(x for i,x in enumerate(labs) if J>>i&1)%2
    states=[J for J in subsets(m) if 0<=k-J.bit_count()<=S and Z[k-J.bit_count()].get((q-charge(J))%2,0)]
    mass={J:weight(J,w)*Z[k-J.bit_count()][(q-charge(J))%2] for J in states}
    C=[J for J in states if a<=k-J.bit_count()<=b]
    aux=subsets(m+d,k-a);W=w+[u]*d;den=comb(m+d,2)
    proj=lambda X:X&((1<<m)-1)
    grouped=defaultdict(list)
    for X in aux:grouped[proj(X)].append(X)
    need(set(grouped)==set(C),'literal auxiliary states have exactly the full interior support')
    F={J:mass[J]/sum(weight(X,W) for X in grouped[J]) for J in C}
    H0={J:defaultdict(Q) for J in C}
    for J in C:
        # Strong lumpability: all individual auxiliary representatives give same row.
        last=None
        for X in grouped[J]:
            row=defaultdict(Q)
            for i,j in combinations(range(m+d),2):
                if ((X>>i)&1)+((X>>j)&1)!=1:continue
                occupied=i if X>>i&1 else j;empty=j if occupied==i else i
                Y=X^(1<<i)^(1<<j);L=proj(Y)
                if L==J:continue
                row[L]+=Q(1,2*den)*W[empty]/(W[i]+W[j])*min(Q(1),F[L]/F[J])
            row[J]=1-sum(row.values())
            if last is not None:need(dict(row)==last,'strong lumpability checked on every auxiliary representative')
            last=dict(row)
        H0[J].update(last)
    kap=max(F.values())/min(F.values());gamma=Q(1,2*den)/kap
    admissions={J:{J:Q(1)} for J in C}
    for J in states:
        if J in admissions:continue
        r=k-J.bit_count();choices=[]
        if r<a:
            for I in combinations([i for i in range(m) if J>>i&1],a-r):choices.append(J^sum(1<<i for i in I))
        else:
            for I in combinations([i for i in range(m) if not J>>i&1],r-b):choices.append(J|sum(1<<i for i in I))
        need(choices and all(L in C for L in choices),'every distributed admission reaches a complete fibre')
        admissions[J]={L:Q(1,len(choices)) for L in choices}
    loads={J:Q(0) for J in C}
    for J in states:
        for L,v in admissions[J].items():loads[L]+=mass[J]*v/mass[L]
    A=[max(z.values()) for z in Z];Ba=min(Z[a].values());Bb=min(Z[b].values());K=max(w+[1/x for x in w])
    H=bound(m,k,S,a,b,A,Ba,Bb,K);actual=max(loads.values())
    need(actual<=H,'all source-weighted admission loads satisfy closed binomial bound')
    # Also audit the tighter source-specific restricted symmetric sums, parent by parent.
    for L in C:
        r=k-L.bit_count();v=Q(1)
        if r==a:
            e=elementary([w[i] for i in range(m) if not L>>i&1],a)
            v+=sum(Q(A[a-j],Ba)*e[j]/comb(L.bit_count()+j,j) for j in range(1,a+1))
        if r==b:
            e=elementary([1/w[i] for i in range(m) if L>>i&1],S-b)
            v+=sum(Q(A[b+j],Bb)*e[j]/comb(m-L.bit_count()+j,j) for j in range(1,S-b+1))
        need(loads[L]<=v,'restricted elementary-symmetric mass bound')
    P={J:defaultdict(Q) for J in states}
    for J in C:
        for L,v in H0[J].items():
            if L!=J:P[J][L]+=v/4
    for J in states:
        if J in C:continue
        for L,v in admissions[J].items():P[J][L]+=v/(4*H);P[L][J]+=mass[J]/mass[L]*v/(4*H)
    for J in states:P[J][J]=1-sum(P[J].values())
    report=audit_kernel(P,mass,gamma/(8*H))
    return {'m':m,'k':k,'q':q,'outside_odds':list(map(str,w)),'H_bound':str(H),'max_actual_load':str(actual),**report}

def main():
    tests=[]
    for w in ([Q(1)]*5,[Q(1,3),Q(3),Q(2),Q(1,2),Q(1)]):
        for k in (1,2,4,5,7,8):
            for q in (0,1):tests.append(small(w,k,q))
    A=[1,1,2,4,9,16,33,61,82,105,113,105,82,61,33,16,9,4,2,1,1]
    native=[]
    for m in (60,61):
        kap=Q(4032,67);gamma=Q(1,2*comb(m+10,2))/kap
        hd=1+9*m+4*comb(m,2)+2*comb(m,3)+comb(m,4)+comb(m,5)
        for k in (5,10,20,30,40,50,m+10,m+15):
            H=bound(m,k,20,5,15,A,1,1,Q(1));native.append({'m':m,'k':k,'H_balanced':str(H),
                'full_gap_bound':str(gamma/(8*H)),'improvement_over_deterministic_parent':str(Q(hd)/H)})
    # The fixed-delta envelope does not depend on outside inventory size.
    envelope=1+9*3+4*9+2*27+81+243
    bulk=[]
    for m in (20,40,80,160,320,640,1280):
        k=m//2+10;H=bound(m,k,20,5,15,A,1,1,Q(1))
        need(H<=envelope,'uniform delta=1/4 bulk envelope')
        bulk.append({'m':m,'k':k,'H':str(H)})
    return {'status':'PASS exact overlapping-parent weighted kernel and complete auxiliary quotient checks',
            'checks':base.CHECKS,'toy_cases':tests,'native_uniform_bounds':native,
            'delta_one_quarter_envelope':envelope,'bulk_controls':bulk,
            'scope':'24 finite exact PSD controls; all-size bulk theorem is proved in BALANCED.md'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
