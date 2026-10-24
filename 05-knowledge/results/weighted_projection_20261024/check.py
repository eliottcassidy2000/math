"""Exact small-carrier audits of the weighted proof. No assert statements.

These are finite controls, not native ROOT weights. Standard library only.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import comb
from collections import defaultdict
from pathlib import Path
import json

CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

def subsets(n,k=None):
    return [x for x in range(1<<n) if k is None or x.bit_count()==k]

def product(xs):
    z=Q(1)
    for x in xs:z*=x
    return z

def weight(x,w):return product(a for i,a in enumerate(w) if x>>i&1)

def elementary(w,d):
    a=[Q(1)]+[Q(0)]*d
    for x in w:
        for j in range(d,0,-1):a[j]+=x*a[j-1]
    return a

def psd_minor(A):
    """Exact LDL on a zero-row-sum symmetric matrix, deleting its last row.

    Zero pivots must have zero remaining row, so semidefinite cases are kept.
    """
    n=len(A)
    need(all(sum(row)==0 for row in A),'quadratic form annihilates constants')
    need(all(A[i][j]==A[j][i] for i in range(n) for j in range(n)),'symmetric form')
    B=[row[:-1].copy() for row in A[:-1]]
    for i in range(n-1):
        p=B[i][i];need(p>=0,'nonnegative exact LDL pivot')
        if p==0:
            need(all(B[i][j]==0 for j in range(i+1,n-1)),'zero pivot has zero row')
            continue
        for j in range(i+1,n-1):
            for k in range(j,n-1):
                B[j][k]-=B[i][j]*B[i][k]/p;B[k][j]=B[j][k]

def audit_kernel(P,mass,bound):
    states=list(P);total=sum(mass.values());pi={x:mass[x]/total for x in states}
    for x in states:
        need(sum(P[x].values())==1 and min(P[x].values())>=0,'stochastic row')
        for y,v in P[x].items():need(pi[x]*v==pi[y]*P[y].get(x,0),'exact detailed balance')
    A=[[pi[x]*((x==y)-P[x].get(y,0))-bound*((pi[x] if x==y else 0)-pi[x]*pi[y])
        for y in states] for x in states]
    psd_minor(A)
    return {'states':len(states),'certified_gap_lower_bound':str(bound),'exact_LDL':'PASS'}

def coupling():
    reports=[]
    for w in ([Q(1)]*6,[Q(1,1000),Q(1000),Q(2),Q(1,7),Q(13),Q(5)]):
        n=len(w);pairs=list(combinations(range(n),2));tested=0
        def outcomes(x,i,j,u):
            if ((x>>i)&1)+((x>>j)&1)!=1:return x
            return (x&~((1<<i)|(1<<j)))|(1<<(i if u<w[i]/(w[i]+w[j]) else j))
        for k in range(1,n):
            for x in subsets(n,k):
                for i in range(n):
                    if not x>>i&1:continue
                    for j in range(n):
                        if x>>j&1:continue
                        y=x^(1<<i)^(1<<j);ed=Q(0)
                        for a,b in pairs:
                            t=w[a]/(w[a]+w[b])
                            for lo,hi in ((Q(0),t),(t,Q(1))):
                                xx=outcomes(x,a,b,(lo+hi)/2);yy=outcomes(y,a,b,(lo+hi)/2)
                                distance=(xx^yy).bit_count()//2
                                need(distance<=1,'synchronous heat bath never increases adjacent distance')
                                ed+=(hi-lo)*distance/len(pairs)
                        need(ed==1-Q(1,len(pairs)),'exact adjacent contraction, arbitrary odds');tested+=1
        reports.append({'weights':list(map(str,w)),'adjacent_ordered_pairs':tested})
    return reports

def coefficients(labels,w,mod):
    Z=[defaultdict(Q) for _ in range(len(labels)+1)]
    for x in subsets(len(labels)):
        g=sum(a for i,a in enumerate(labels) if x>>i&1)%mod
        Z[x.bit_count()][g]+=weight(x,w)
    return Z

def toy(m,w,k,q,u=Q(1)):
    labels=[0,1,0,1];rw=[Q(1),Q(2),Q(3),Q(1,2)];S=4;a=1;b=3;d=b-a
    outside=[i%2 for i in range(m)];Z=coefficients(labels,rw,2)
    charge=lambda x:sum(v for i,v in enumerate(outside) if x>>i&1)%2
    allstates=[x for x in subsets(m) if 0<=k-x.bit_count()<=S and Z[k-x.bit_count()].get((q-charge(x))%2,0)>0]
    mass={x:weight(x,w)*Z[k-x.bit_count()][(q-charge(x))%2] for x in allstates}
    C=[x for x in allstates if a<=k-x.bit_count()<=b]
    F={x:Z[k-x.bit_count()][(q-charge(x))%2]/comb(d,k-x.bit_count()-a)/u**(k-x.bit_count()-a) for x in C}
    allF=[Z[r][g]/comb(d,r-a)/u**(r-a) for r in range(a,b+1) for g in range(2)]
    kap=max(allF)/min(allF);D=comb(m+d,2);gamma=Q(1,2*D)/kap
    H0={x:defaultdict(Q) for x in C}
    for x in C:
        r=k-x.bit_count();s=r-a
        proposals=defaultdict(Q)
        for i in range(m):
            if x>>i&1:
                if r<b:proposals[x^(1<<i)]+=Q(d-s,D)*u/(u+w[i])
                for j in range(m):
                    if not x>>j&1:proposals[x^(1<<i)^(1<<j)]+=Q(1,D)*w[j]/(w[i]+w[j])
            elif r>a:proposals[x^(1<<i)]+=Q(s,D)*w[i]/(u+w[i])
        for y,v in proposals.items():
            need(y in F,'reference proposals stay in complete band')
            H0[x][y]+=v*min(Q(1),F[y]/F[x])/2
        H0[x][x]=1-sum(H0[x].values())
    inner=audit_kernel(H0,{x:mass[x] for x in C},gamma)
    ew=elementary(w,a);ei=elementary([1/x for x in w],S-b)
    H=1+max(sum(max(Z[a-j].values())/min(Z[a].values())*ew[j] for j in range(1,a+1)),
            sum(max(Z[b+j].values())/min(Z[b].values())*ei[j] for j in range(1,S-b+1)))
    parent={x:x for x in C}
    for x in allstates:
        if x in parent:continue
        r=k-x.bit_count();y=x
        if r<a:
            for i in [i for i in range(m) if x>>i&1][:a-r]:y^=1<<i
        else:
            for i in [i for i in range(m) if not x>>i&1][:r-b]:y|=1<<i
        need(y in C,'boundary parent is complete');parent[x]=y
    loads={x:Q(0) for x in C}
    for x in allstates:loads[parent[x]]+=mass[x]/mass[parent[x]]
    need(max(loads.values())<=H,'actual entire attached mass bounded')
    P={x:defaultdict(Q) for x in allstates}
    for c in C:
        for z,p in H0[c].items():
            if z!=c:P[c][z]+=p/4
    for x in allstates:
        c=parent[x]
        if x==c:continue
        P[x][c]+=Q(1)/(4*H);P[c][x]+=mass[x]/mass[c]/(4*H)
    for x in allstates:P[x][x]=1-sum(P[x].values())
    full=audit_kernel(P,mass,gamma/(8*H))
    need(min(P[x][x] for x in allstates)>=Q(1,2),'full chain is lazy')
    # Independent row-by-row energy identity after nonuniform conditional refresh.
    blocks={x:[(x,0),(x,1)] for x in allstates};cond={x:[Q(1,x+2),Q(x+1,x+2)] for x in allstates}
    total=sum(mass.values());f={(x,j):Q((3*x+7*j)%11,3) for x in allstates for j in (0,1)}
    mean={x:sum(cond[x][j]*f[x,j] for j in (0,1)) for x in allstates}
    cvar=sum(mass[x]/total*sum(cond[x][j]*(f[x,j]-mean[x])**2 for j in (0,1)) for x in allstates)
    energy=Q(0);macro=Q(0)
    for x in allstates:
        for y,p in P[x].items():
            macro+=mass[x]/total*p*(mean[x]-mean[y])**2/2
            for i in (0,1):
                for j in (0,1):energy+=mass[x]/total*cond[x][i]*p*cond[y][j]*(f[x,i]-f[y,j])**2/2
    need(energy==cvar+macro,'exact full conditional-refresh energy identity')
    # Global activity gauge: scaling ALL V odds scales coefficients by t^r.
    scale=Q(7);ZZ=coefficients(labels,[scale*x for x in rw],2)
    need(all(ZZ[r][g]==scale**r*v for r,z in enumerate(Z) for g,v in z.items()),'reservoir activity covariance')
    scaledF=[ZZ[r][g]/comb(d,r-a)/(scale*u)**(r-a) for r in range(a,b+1) for g in range(2)]
    need(max(scaledF)/min(scaledF)==kap,'activity-matched comparison constant is gauge invariant')
    eww=elementary([scale*x for x in w],a)
    need(all(max(ZZ[a-j].values())/min(ZZ[a].values())*eww[j]==
             max(Z[a-j].values())/min(Z[a].values())*ew[j] for j in range(1,a+1)),'boundary debt is gauge invariant')
    return {'m':m,'k':k,'q':q,'odds':list(map(str,w)),'auxiliary_odds':str(u),'kappa':str(kap),'H':str(H),
            'max_actual_cluster_ratio':str(max(loads.values())),'interior':inner,'full':full}

def hostile():
    # Small analogue: order4 charge, two groups of4 outside slots; budget<8.
    # The charge-complete reservoir is a single chosen residue 0,1,2,3.
    labs=[0]*4+[1]*4+[0,1,2,3];polynomial=[0]*5;heavy=[]
    for x in subsets(12,4):
        if sum(v for i,v in enumerate(labs) if x>>i&1)%4:continue
        r=(x>>8).bit_count();polynomial[r]+=1
        if not r:heavy.append(x)
    need(heavy==[15,240],'precisely two zero-reservoir heavy states')
    need((heavy[0]^heavy[1]).bit_count()==8,'bounded-support chain cannot directly join heavy states')
    bounds=[]
    for T in (10,100,1000,10000):
        Z=sum(c*T**(4-r) for r,c in enumerate(polynomial));p=Q(T**4,Z)
        outside=1-2*p;upper=outside/(p*(1-p));bounds.append(str(upper))
    need(all(Q(x)>Q(y) for x,y in zip(bounds,bounds[1:])),'mass barrier worsens despite equal outside odds')
    return {'label_carrier_only':True,'mass_coefficients_by_reservoir_grade':polynomial,
            'T':[10,100,1000,10000],'gap_upper_bounds':bounds}

def main():
    cp=coupling();cases=[]
    for w in ([Q(1)]*4,[Q(1,10),Q(10),Q(2),Q(1,3)]):
        for k in (1,2,4,6,7):
            for q in (0,1):
                for u in (Q(1),Q(2)):cases.append(toy(4,w,k,q,u))
    h=hostile()
    return {'status':'PASS exact rational coupling, detailed balance, PSD gap certificates and weight hostile',
            'checks':CHECKS,'coupling':cp,'toy_cases':cases,'hostile':h,
            'scope':'finite algebraic controls; universal bounds follow PROOF.md; no toy law is claimed to be a ROOT law'}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
